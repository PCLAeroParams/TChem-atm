/*
NOTE: The Batched GMRES solver specifically interfaces with CVODE's memory, so it currently does
not work with other SUNDIALS solvers. 
*/
#pragma once
#include <cstdlib> // std::getenv, std::atol, std::atoi (For diagnostics)
#include <KokkosBatched_GMRES.hpp>
#include "TChem.hpp"
#include "TChem_Impl_Jacvec_FD.hpp"
#include "TChem_Impl_AerosolChemistryRHS.hpp"
#include <sundials/sundials_types.h> // sunrealtype?
#include <nvector/nvector_kokkos.hpp>
#include <sundials/sundials_core.hpp>
#include <cvode/cvode.h>

#define BGMR_MAXITER_DEFAULT 20

namespace TChem{

static_assert(std::is_same<real_type, sunrealtype>::value,
            "SUNDIALS and TChem real types differ; the N_Vector reinterpretations are invalid");
using device_type = typename Tines::UseThisDevice<TChem::exec_space>::type;

// Krylov typedefs
using MagType          = Kokkos::ArithTraits<sunrealtype>::mag_type;
using Norm2DView       = Tines::value_type_2d_view<MagType, device_type>;
using Scalar3DView     = Tines::value_type_3d_view<sunrealtype, device_type>;
using IntView1D        = Tines::value_type_1d_view<ordinal_type, device_type>;
using KrylovHandleType = KokkosBatched::KrylovHandle<Norm2DView, IntView1D, Scalar3DView>;

template <typename RHSType, typename PreconType>
struct BGMRTeamFunctor { 

  // typedefs
  using real_type_1d_view_type = Tines::value_type_1d_view<real_type, device_type>;
  using real_type_2d_view_type = Tines::value_type_2d_view<real_type, device_type>;

  // members 
  RHSType f;                          // System-level RHS function
  PreconType precon;                  // Preconditioner (Identity or Schur)
  real_type_2d_view_type y0;          // system linearization point (nBatch, m)
  real_type_2d_view_type delta;       // [GENERATED IN SCALED SPACE (initially 0), RETURNED UNSCALED] the unknown, solution to the linear system formed by (I-gamma*J)delta = b (input 0, output solution)
  real_type_2d_view_type f_y0;        // cached f(y0), passed from CVodeGetNonlinearSystemData
  real_type_2d_view_type b;           // [SCALED PRIOR TO PASSING] RHS of the linear system we are solving against
  real_type gamma;                    // a constant specific to the CVODE method (shared for each system)
  real_type_2d_view_type ewt;         // error weights for calculating the WRMS-scaled difference quotient Jv size(nBatch, m)
  KrylovHandleType handle;            // Kokkos Kernels Krylov attributes 
  ordinal_type team_scratch_work_size;// team-level workspace size (in count of reals) to be carved up in the functor operator()
  long nsolves;                       // diagnostic: solve index, for labelling the print
  int diag_system;                    // diagnostic: league rank whose apply() calls print a Jv breakdown; -1 = none

  // constructor
  BGMRTeamFunctor(RHSType f_,
                PreconType precon_,
                real_type_2d_view_type y0_,
                real_type_2d_view_type delta_,
                real_type_2d_view_type f_y0_,
                real_type_2d_view_type b_,
                real_type gamma_,
                real_type_2d_view_type ewt_,
                KrylovHandleType handle_,
                ordinal_type team_scratch_work_size_,
                long nsolves_ = 0,
                int diag_system_ = -1
            ) :
            f(f_),
            precon(precon_), 
            y0(y0_), 
            delta(delta_), 
            f_y0(f_y0_), 
            b(b_), 
            gamma(gamma_), 
            ewt(ewt_), 
            handle(handle_), 
            team_scratch_work_size(team_scratch_work_size_),
            nsolves(nsolves_), diag_system(diag_system_) {}

  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION void operator()(const MemberType &member) const{ 

    const int system_size   = y0.extent(1);
    const int i_sys = static_cast<int>(member.league_rank());

    // --- Partition scratch (level-1) memory for RHS and preconditioning ---
    // claim workspace region of memory
    ordinal_type level = 1;
    TChem::Scratch<real_type_1d_view_type> work(member.team_scratch(level), team_scratch_work_size);
    // pointer to the start of the workspace
    auto wptr = work.data();
    // get the size of the rhs workspace
    const ordinal_type rhs_ws = TChem::Impl::Aerosol_RHS<real_type, device_type>::getWorkSpaceSize(f.kmcd, f.amcd);
    // assign real_type_1d_view_type for the rhs workspace
    auto pw = real_type_1d_view_type(wptr, rhs_ws);  
    // advance work pointer to next region of scratch memory and assign y_shift to this chunk
    wptr += rhs_ws;
    auto y_shift = real_type_1d_view_type(wptr, system_size);       
    wptr += system_size;
    // advance work pointer to next region of scratch memory and assign f_shift to this chunk
    auto f_shift = real_type_1d_view_type(wptr, system_size);       
    wptr += system_size;
    auto u = real_type_1d_view_type(wptr, system_size); // preconditioning buffer
    wptr += system_size;

    // generate 2D slices of the linearization point y0, the system rhs b, and the ODE RHS function f()
    // size (1, n) views (expected by GMRES)
    auto delta_slice    = Kokkos::subview(delta, Kokkos::make_pair(i_sys, i_sys+1), Kokkos::ALL);
    auto b_slice    = Kokkos::subview(b, Kokkos::make_pair(i_sys, i_sys+1), Kokkos::ALL);
    auto f_slice = f.team_slice(i_sys, pw); 

    // 1D slices of x and f(x)
    auto delta_i = Kokkos::subview(delta, member.league_rank(), Kokkos::ALL); // size (n) views
    auto y0_i = Kokkos::subview(y0, member.league_rank(), Kokkos::ALL); 
    auto f_y0_i = Kokkos::subview(f_y0, member.league_rank(), Kokkos::ALL);
    auto ewt_i = Kokkos::subview(ewt, member.league_rank(), Kokkos::ALL); // get error weights for the current team/system

    // declare operator for Matrix-free operations on the linear system
    using Operator = Impl::JacvecFDTeam<RHSType, PreconType>;
    Operator A(f_slice, precon, y0_i, f_y0_i, gamma, ewt_i, y_shift, f_shift, u,
                i_sys == diag_system, nsolves);
        
    // GMRES expects the 2D array versions
    KokkosBatched::TeamVectorGMRES<MemberType>::template invoke<Operator, real_type_2d_view_type>(
        member, A, b_slice, delta_slice, handle);
    member.team_barrier();

    // Unscale delta before returning
    Kokkos::parallel_for(Kokkos::TeamVectorRange(member, system_size), [&](int j){
        delta_i(j) /= ewt_i(j);
    });
    member.team_barrier();

    // Solved ~ ((I-gamma*J)P^{-1})delta_u = -F 
    // Note that delta_u = P*delta_true, now solve for Pdelta_true = delta_u 
    // where delta_true = P^{-1}delta_u
    precon.team_solve(member, delta_i); // input is delta_u, output is delta_true
    member.team_barrier();

  }
};

struct SUNLinearSolverContent_BatchedGMRES 
{
    // typedefs
    using real_type_2d_view_type = Tines::value_type_2d_view<real_type, device_type>;
    using real_type_1d_view_type = Tines::value_type_1d_view<real_type, device_type>;
    using problem_type = TChem::Impl::AerosolChemistry_Problem<real_type, device_type>;
    using policy_type = typename UseThisTeamPolicy<exec_space>::type;

    // members
    KrylovHandleType handle;
    Impl::AerosolChemistryRHS<problem_type> rhs_object;
    ordinal_type numiters = 0;
    ordinal_type maxit = 0;
    ordinal_type system_size;
    ordinal_type n_systems; // nBatch
    ordinal_type gmres_max_iter;
    ordinal_type team_scratch_work_size = 0;
    ordinal_type last_flag;
    real_type resnorm;
    
    // Diagnostic parameters for dumping state at a given solve
    long nsolves = 0;        
    bool verbose = false;    // diagnostic printing; off by default, set from the driver
    /*
    std::string dump_file = "";
    long        dump_solve = -1;      // which solve index to capture, <0 disables
    real_type   dump_rtol = 0, dump_atol = 0;
    bool        dumped = false;
    */

    void* cvode_mem;
    SUNPSetupFn Psetup;
    SUNPSolveFn Psolve;
    void* PData;
    int pretype;
    N_Vector s1 = nullptr;

    policy_type policy; 
    policy_type norm_policy;

    // constructor
    SUNLinearSolverContent_BatchedGMRES(
                Impl::AerosolChemistryRHS<problem_type> rhs_object_, 
                ordinal_type system_size_, 
                ordinal_type n_systems_,
                ordinal_type gmres_max_iter_,
                int pretype_,
                //real_type gmres_tol_,
                policy_type policy_
            ) :
            handle(n_systems_, /*systems_per_team=*/1, gmres_max_iter_, /*monitor_residual=*/true),
            rhs_object(rhs_object_),
            system_size(system_size_),
            n_systems(n_systems_),
            gmres_max_iter(gmres_max_iter_),
            pretype(pretype_),
            //gmres_tol(gmres_tol_),
            policy(policy_) 
            {
                // use the init function for initialization instead of the constructor since 
                // it becomes solver's init operator and gets called by CVODE
            }

    int init(){

        // check cvode_mem has been set upstream
        // assert cvode_mem is not null
        if (cvode_mem == nullptr){
            return -1; // setatimes did not run before initializing the solver
        }

        if (pretype == SUN_PREC_RIGHT){
            if (Psetup == nullptr){
                // user defined preconditioner setup function schur_psetup must be attached to CVODE 
                return -1; 
            }
        }

        norm_policy = policy_type(n_systems, Kokkos::AUTO());

        ordinal_type scratch_level = 1;
        // Size scratch memory for RHS views, preconditioner buffer, and GMRES
        // | ---`work`, carved in functor----||--------------handled by GMRES--------------|
        //  *-------- rhs --------* *-precon-* *------------------GMRES--------------------*
        // [ pw | y_shift | f_shift |   u    | resid. vector | working vector | mask | tmp ]

        // ---- RHS scratch memory -----        
        ordinal_type rhs_scratch_size = 
        TChem::Impl::Aerosol_RHS<real_type, device_type>::getWorkSpaceSize(rhs_object.kmcd, 
            rhs_object.amcd) + 2*system_size;
    
        // ---- Preconditioner scratch memory ----
        ordinal_type precon_scratch_size = system_size;

        team_scratch_work_size = rhs_scratch_size + precon_scratch_size;
        const size_t team_scratch_work_bytes =
        TChem::Scratch<real_type_1d_view_type>::shmem_size(team_scratch_work_size);

        // ---- GMRES scratch memory ----
        // originally used level 0 for GMRES scratch, but with large systems (n_particles > 100), the scratch allocation exceeds shared memory (~48 KB for L40S)
        ordinal_type systems_per_team = 1; // each team works with 1 chemical system (key assumption of this solver)
        maxit = (gmres_max_iter < system_size) ? gmres_max_iter : system_size; // for small systems with size < gmres_max_iter, we only need at most system_size # of iterations
        const size_t team_scratch_gmres_bytes  = real_type_2d_view_type::shmem_size(systems_per_team, (maxit + 1) + system_size + 2); // resid vector + working vector + mask + tmp
        
        policy.set_scratch_size(scratch_level, Kokkos::PerTeam(
            team_scratch_work_bytes + team_scratch_gmres_bytes));

        // Configure the Krylov handle and Arnoldi view
        handle.set_scratch_pad_level(1); // use global scratch memory
        handle.set_ortho_strategy(1); // use modified gram-schmidt for consistency with SUNDIALS
        handle.set_compute_last_residual(false);
        handle.Arnoldi_view = KrylovHandleType::ArnoldiViewType(
            "Arnoldi_view", n_systems, gmres_max_iter, system_size + gmres_max_iter + 3);

        /*
        if (verbose) {
          const double arnoldi_mb = double(n_systems) * double(gmres_max_iter)
                                  * double(system_size + gmres_max_iter + 3)
                                  * double(sizeof(real_type)) / (1024.0*1024.0);
          printf("..BatchedGMRES attached: n_systems = %d, system_size = %d, max_iter = %d "
                 "(GMRES clamps to min(max_iter, m) = %d)\n",
                 n_systems, system_size, gmres_max_iter, maxit);
          printf("..BatchedGMRES scratch : RHS+precon = %d reals/team (%zu B), "
                 "GMRES _TMPView = %zu B/team, total level-1 = %zu B/team\n",
                 team_scratch_work_size, (size_t)team_scratch_work_bytes, (size_t)team_scratch_gmres_bytes,
                 (size_t)team_scratch_work_bytes + (size_t)team_scratch_gmres_bytes);
          printf("..BatchedGMRES handle scratch_pad_level = %d (must be 1), ortho_strategy = %d "
                 "(1 = modified GS, matching SPGMR)\n",
                 handle.get_scratch_pad_level(), handle.get_ortho_strategy());
          printf("..BatchedGMRES Arnoldi_view = (%d, %d, %d) = %.2f MB\n",
                 n_systems, gmres_max_iter, system_size + gmres_max_iter + 3, arnoldi_mb);
        }
        */

        return SUN_SUCCESS; // alias for 0
        // typical return values
        //      success: 0, unrecoverable: -1, recoverable: 1
    }

    int solve(N_Vector x, N_Vector b, sunrealtype tol){

        // the scaling vector s1 is set per solve by CVODE, so it should no longer be the nullptr 
        // value set when this struct is instantiated
        if (s1 == nullptr){
            return -1; // setscalingvectors did not run before calling solve
        }

        // local variables for accessing the current CVODE non-linear system data
        sunrealtype tcur = 0, gamma = 0, rl1 = 0;
        N_Vector ypred = nullptr, yn = nullptr, fn = nullptr, zn1 = nullptr; // nullptr is a placeholder value
        void* userdata_ptr = nullptr;

        CVodeGetNonlinearSystemData(cvode_mem, &tcur, &ypred, &yn, &fn,
                                    &gamma, &rl1, &zn1, &userdata_ptr);

        // reinterpret N_Vectors as 2D views
        real_type_2d_view_type y0_2d(N_VGetDeviceArrayPointer(yn), n_systems, system_size);
        real_type_2d_view_type fy0_2d(N_VGetDeviceArrayPointer(fn), n_systems, system_size);
        real_type_2d_view_type b_2d(N_VGetDeviceArrayPointer(b), n_systems, system_size);
        real_type_2d_view_type delta_2d(N_VGetDeviceArrayPointer(x), n_systems, system_size);
        real_type_2d_view_type ewt_2d(N_VGetDeviceArrayPointer(s1), n_systems, system_size);

        /*
        // ---- optional state dump, before b is scaled in place ----
        if (!dump_file.empty() && !dumped && nsolves == dump_solve) {
          sunrealtype h_last = 0; int q_cur = 0; long int nst_cur = 0, nli_cur = 0, njt_cur = 0;
          CVodeGetLastStep(cvode_mem, &h_last);
          CVodeGetCurrentOrder(cvode_mem, &q_cur);
          CVodeGetNumSteps(cvode_mem, &nst_cur);
          CVodeGetNumLinIters(cvode_mem, &nli_cur);
          CVodeGetNumJtimesEvals(cvode_mem, &njt_cur);

          auto y_h = Kokkos::create_mirror_view(y0_2d);   Kokkos::deep_copy(y_h,  y0_2d);
          auto f_h = Kokkos::create_mirror_view(fy0_2d);  Kokkos::deep_copy(f_h,  fy0_2d);
          auto e_h = Kokkos::create_mirror_view(ewt_2d);  Kokkos::deep_copy(e_h,  ewt_2d);
          auto b_h = Kokkos::create_mirror_view(b_2d);    Kokkos::deep_copy(b_h,  b_2d);

          FILE* fd = fopen(dump_file.c_str(), "w");
          if (fd == nullptr) {
            printf("WARNING: could not open state dump file %s\n", dump_file.c_str());
          } else {
            fprintf(fd, "# BatchedGMRES state dump at solve %ld\n", nsolves);
            fprintf(fd, "# columns per row: y  ewt  f(y)  b\n");
            fprintf(fd, "nbatch %d\n", n_systems);
            fprintf(fd, "neq %d\n", system_size);
            fprintf(fd, "gamma %24.17e\n", gamma);
            fprintf(fd, "rtol %24.17e\n", dump_rtol);
            fprintf(fd, "atol %24.17e\n", dump_atol);
            fprintf(fd, "t %24.17e\n", tcur);
            fprintf(fd, "h %24.17e\n", h_last);
            fprintf(fd, "q %d\n", q_cur);
            fprintf(fd, "nst %ld\n", nst_cur);
            fprintf(fd, "nli %ld\n", nli_cur);
            fprintf(fd, "njtimes %ld\n", njt_cur);
            for (ordinal_type sidx = 0; sidx < n_systems; sidx++) {
              fprintf(fd, "sample %d\n", sidx);
              for (ordinal_type j = 0; j < system_size; j++) {
                fprintf(fd, "%24.17e %24.17e %24.17e %24.17e\n",
                        y_h(sidx,j), e_h(sidx,j), f_h(sidx,j), b_h(sidx,j));
              }
            }
            fclose(fd);
            printf("..BatchedGMRES dumped state at solve %ld to %s "
                   "(t = %.6e, h = %.6e, q = %d, nst = %ld, gamma = %.6e)\n",
                   nsolves, dump_file.c_str(), tcur, h_last, q_cur, nst_cur, gamma);
          }
          dumped = true;
        }
        */

        // scale linear system RHS b and compute GMRES scalar relative tolerance 
        // (CVODE tol is absolute for each system, so we need to divide tol by the 
        // *max* of the 2 norm of b across the batch of systems)
        const ordinal_type m = system_size;
        real_type bnorm_max = 0;
        Kokkos::parallel_reduce(norm_policy, KOKKOS_LAMBDA (const policy_type::member_type& member, real_type& lmax) {
            const ordinal_type i = member.league_rank();
            real_type bsqr = 0;
            Kokkos::parallel_reduce(Kokkos::TeamVectorRange(member, m), [=] (const int& j, real_type& sum) {
                const real_type bs = ewt_2d(i,j)*b_2d(i,j);
                b_2d(i,j) = bs; // scale b in place
                sum += bs*bs;
            }, bsqr);

            // a single thread in each team participates in this parallel operation
            Kokkos::single(Kokkos::PerTeam(member), [&]() {
                const real_type bnorm_i = Kokkos::ArithTraits<real_type>::sqrt(bsqr);
                lmax = lmax > bnorm_i ? lmax : bnorm_i;
            });
            
        }, Kokkos::Max<real_type>(bnorm_max));

        // return early if b is zero for every system---the exact solution is x = 0 (handled by CVODE)
        if (!(bnorm_max > 0)){ 
            numiters = 0;
            resnorm = 0;
            last_flag = SUN_SUCCESS;
            return last_flag; 
        }

        // TODO possible to set tolerance as a view with each entry tol (scalar)/(||b_i||) for each system i?
        handle.set_tolerance(tol/bnorm_max);
        handle.reset(); // reset the iteration numbers and residual norms

        
        // Diagnostic 
        static const long diag_solve  = std::getenv("BGMR_DIAG_SOLVE")  ? std::atol(std::getenv("BGMR_DIAG_SOLVE"))  : -1;
        static const int  diag_system = std::getenv("BGMR_DIAG_SYSTEM") ? std::atoi(std::getenv("BGMR_DIAG_SYSTEM")) : 0;
        const int diag_sys = (nsolves == diag_solve) ? diag_system : -1;

        // Branch solve approach on the appropriate preconditioner 
        if (pretype == SUN_PREC_RIGHT){
            using PreconType = SchurPrecon;
            auto udata = static_cast<UserData*>(userdata_ptr);
            PreconType precon_object(fy0_2d, y0_2d, udata->aero_num_conc,
                udata->const_tracers, udata->temperature, udata->pressure,
                udata->kmcd, udata->amcd, udata->n_particles_track,
                udata->Jgg, udata->Jgp, udata->Jpg, udata->Jpp, udata->D,
                udata->S, udata->W, udata->gamma);

            // Create GMRES functor
            BGMRTeamFunctor<Impl::AerosolChemistryRHS<problem_type>, PreconType> functor(
                rhs_object, precon_object, y0_2d, delta_2d, fy0_2d, b_2d, gamma, ewt_2d, 
                handle, team_scratch_work_size, nsolves, diag_sys);
            
            // GMRES solve in parallel over teams
            Kokkos::parallel_for(policy, functor);
            Kokkos::fence();
        }
        else { // pretype == SUN_PREC_NONE
            using PreconType = IdentityPrecon;
            PreconType precon_object; // Identity has empty constructor
            // Create GMRES functor
            BGMRTeamFunctor<Impl::AerosolChemistryRHS<problem_type>, PreconType> functor(
                rhs_object, precon_object, y0_2d, delta_2d, fy0_2d, b_2d, gamma, ewt_2d, 
                handle, team_scratch_work_size, nsolves, diag_sys);
            
            // GMRES solve in parallel over teams
            Kokkos::parallel_for(policy, functor);
            Kokkos::fence();
        }
        
        /*
        // Diagnostic: print residual history
        if (diag_sys >= 0 && diag_sys < n_systems) {
            const int it   = handle.get_iteration_host(diag_sys);
            const int jmax = (it >= 0) ? it : maxit;
            printf("  [bgmr-diag] solve %ld sys %d: KK iteration = %d, KK tolerance = %.16e, stored ||r0|| = %.16e\n",
                   nsolves, diag_sys, it, (double)(tol/bnorm_max), (double)handle.get_norm_host(diag_sys, 0));
            for (int j = 1; j <= jmax; ++j) {
                printf("  [bgmr-diag]    stored res_norm[%d] = % .16e\n", j, (double)handle.get_norm_host(diag_sys, j));
            }
        }
        */
        
        // Check for any unconverged systems
        auto h = handle;
        int not_conv = 0;
        Kokkos::parallel_reduce(Kokkos::RangePolicy<exec_space>(0, n_systems), KOKKOS_LAMBDA (const ordinal_type i, ordinal_type& itermax, real_type& resnorm_max, int& n_unconv) {
            itermax = itermax > h.get_iteration(i) ? itermax : h.get_iteration(i);
            resnorm_max = resnorm_max > h.get_last_norm(i) ? resnorm_max : h.get_last_norm(i);
            // GMRES returns iteration count = -1 if it did not converge within the max # of iterations
            if (h.get_iteration(i) == -1){
                n_unconv += 1;
            }
        }, Kokkos::Max<ordinal_type>(numiters), Kokkos::Max<real_type>(resnorm), Kokkos::Sum<int>(not_conv));

        // failed solves meant max # iterations was reach by at least 1 system
        if (not_conv != 0){
            numiters = maxit;
        }

        if ((not_conv == 0) && (resnorm <= tol)) {
            last_flag = SUN_SUCCESS;
        } else {
            last_flag = SUNLS_RES_REDUCED;
        }

        /*
        // Solver diagnostics:
        // nsolves: Linear solver iteration 
        // numiters: MAX number of iterations required by GMRES across the batch of systems
        // resnorm: residual Ax* - b~ for the scaled system with candidate solution x*
        // bnorm_max: The maximum 2-norm of linear system RHSs across all systems in the batch (used to set GMRES tolerance)
        // reduction: resnorm/||b~||, <<1 for successful linear solves; if ~ 1, GMRES was not able to converge on an adequate solution
        // tol: CVODE absolute tolerance
        // tol/bnorm_max: relative tolerance used by GMRES
        // not_conv: number of systems GMRES failed to converge before it hit the max number of Arnoldi iterations
        // n_systems: number of independent systems in the batch 
        if (verbose) {
          const real_type reduction = (bnorm_max > 0) ? resnorm/bnorm_max : 0.0;
          printf("  [bgmr] solve %ld: iters(max) = %d, resnorm(max) = %.16e, ||b~|| = %.16e, "
                 "reduction = %.16e, tol = %.16e (rel = %.16e), not converged = %d/%d%s\n",
                 nsolves, numiters, resnorm, bnorm_max, reduction, tol, tol/bnorm_max,
                 not_conv, n_systems,
                 (last_flag == SUN_SUCCESS) ? "" : "   <-- FAIL");
        }
        */

        nsolves++;

        return last_flag;
    }

};

SUNLinearSolver SUNLinSol_BatchedGMRES(N_Vector y, 
    Impl::AerosolChemistryRHS<SUNLinearSolverContent_BatchedGMRES::problem_type> rhs_object,
    ordinal_type system_size, ordinal_type n_systems, int pretype,
    ordinal_type gmres_max_iter, SUNLinearSolverContent_BatchedGMRES::policy_type policy,
    SUNContext sunctx);
SUNErrCode SUNLinSol_BatchedGMRESSetVerbose(SUNLinearSolver S, bool onoff);
SUNErrCode SUNLinSol_BatchedGMRESSetStateDump(SUNLinearSolver S, const char* filename,
                                              long solve_index, sunrealtype rtol, sunrealtype atol);
SUNLinearSolver_Type SUNLinSolGetType_BGMR(SUNLinearSolver S);
SUNLinearSolver_ID SUNLinSolGetID_BGMR(SUNLinearSolver S);
SUNErrCode SUNLinSolSetATimes_BGMR(SUNLinearSolver S, void* ATData,
                                    SUNATimesFn ATimes);
SUNErrCode SUNLinSolSetScalingVectors_BGMR(SUNLinearSolver S, N_Vector s1,
        N_Vector s2);
SUNErrCode SUNLinSolSetZeroGuess_BGMR(SUNLinearSolver S, sunbooleantype onff);
SUNErrCode SUNLinSolInitialize_BGMR(SUNLinearSolver S);
SUNErrCode SUNLinSolSetPreconditioner_BGMR(SUNLinearSolver S, void* P_data,
                                            SUNPSetupFn Pset, SUNPSolveFn Psol);
int SUNLinSolSetup_BGMR(SUNLinearSolver S, SUNMatrix A);
int SUNLinSolSolve_BGMR(SUNLinearSolver S, SUNMatrix A,
                         N_Vector x, N_Vector b, sunrealtype delta);
int SUNLinSolNumIters_BGMR(SUNLinearSolver S);
sunrealtype SUNLinSolResNorm_BGMR(SUNLinearSolver S);
sunindextype SUNLinSolLastFlag_BGMR(SUNLinearSolver S);
SUNErrCode SUNLinSolFree_BGMR(SUNLinearSolver S);

} // namespace TChem