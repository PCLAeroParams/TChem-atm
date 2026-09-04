/*
NOTE: The Batched GMRES solver specifically interfaces with CVODE's memory, so it currently does
not work with other SUNDIALS solvers. 
*/
#pragma once
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

template <typename RHSFunctor>
struct MatrixFreeBatchedGMRESTeamFunctor { 

  // typedefs
  using real_type_1d_view_type = Tines::value_type_1d_view<real_type, device_type>;
  using real_type_2d_view_type = Tines::value_type_2d_view<real_type, device_type>;

  // members 
  RHSFunctor f;                       // device-callable for returning the RHS given a state x0
  real_type_2d_view_type y0;          // system linearization point (nBatch, m)
  real_type_2d_view_type delta;       // [GENERATED IN SCALED SPACE (initially 0), RETURNED UNSCALED] the unknown, solution to the linear system formed by (I-gamma*J)delta = b (input 0, output solution)
  real_type_2d_view_type f_y0;        // cached f(y0), passed from CVodeGetNonlinearSystemData
  real_type_2d_view_type b;           // [SCALED PRIOR TO PASSING] RHS of the linear system we are solving against
  real_type gamma;                    // a constant specific to the CVODE method (shared for each system)
  real_type_2d_view_type ewt;         // error weights for calculating the WRMS-scaled difference quotient Jv size(nBatch, m)
  KrylovHandleType handle;            // Kokkos Kernels Krylov attributes 
  ordinal_type per_team_extent;       // level-1 scratch memory size used by the jac-vec calculation (in bits/bytes?)
  bool verbose;                       // diagnostic: recompute f(y0) and compare against the passed f_y0
  long nsolves;                       // diagnostic: solve index, for labelling the print

  // constructor 
  MatrixFreeBatchedGMRESTeamFunctor(RHSFunctor f_, 
                                    real_type_2d_view_type y0_, 
                                    real_type_2d_view_type delta_, 
                                    real_type_2d_view_type f_y0_, 
                                    real_type_2d_view_type b_,  
                                    real_type gamma_, 
                                    real_type_2d_view_type ewt_,
                                    KrylovHandleType handle_, 
                                    ordinal_type per_team_extent_,
                                    bool verbose_ = false,
                                    long nsolves_ = 0
                                ) :
                                f(f_), 
                                y0(y0_), 
                                delta(delta_), 
                                f_y0(f_y0_), 
                                b(b_), 
                                gamma(gamma_), 
                                ewt(ewt_), 
                                handle(handle_), 
                                per_team_extent(per_team_extent_),
                                verbose(verbose_),
                                nsolves(nsolves_)  {}

  template <typename MemberType>
  KOKKOS_INLINE_FUNCTION void operator()(const MemberType &member) const{ 

    const int system_size   = y0.extent(1);
    const int first = static_cast<int>(member.league_rank());
    const int last  = first + 1;

    // --- Partition scratch (level-1) memory [ pw | x_shift | f_shift ] ---
    // claim workspace region of memory
    ordinal_type level = 1;
    TChem::Scratch<real_type_1d_view_type> work(member.team_scratch(level), per_team_extent);
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

    // generate 2D slices of the linearization point y0, the system rhs b, and the ODE RHS function f()
    // size (1, n) views (expected by GMRES)
    auto delta_slice    = Kokkos::subview(delta,   Kokkos::make_pair(first, last), Kokkos::ALL);
    auto b_slice    = Kokkos::subview(b,   Kokkos::make_pair(first, last), Kokkos::ALL);
    auto f_slice = f.team_slice(first, last, pw); 

    // 1D slices of x and f(x)
    auto delta_i = Kokkos::subview(delta, member.league_rank(), Kokkos::ALL); // size (n) views
    auto y0_i = Kokkos::subview(y0, member.league_rank(), Kokkos::ALL); 
    auto f_y0_i = Kokkos::subview(f_y0, member.league_rank(), Kokkos::ALL);
    auto ewt_i = Kokkos::subview(ewt, member.league_rank(), Kokkos::ALL); // get error weights for the current team/system

    // declare operator for Matrix-free operations on the linear system
    using Operator = Impl::JacvecFDTeam<RHSFunctor>;
    Operator A(f_slice, y0_i, f_y0_i, gamma, ewt_i, y_shift, f_shift);
    
    // LATER: could pass preconditioner at some point too
    
    // GMRES expects the 2D array versions
    KokkosBatched::TeamVectorGMRES<MemberType>::template invoke<Operator, real_type_2d_view_type>(
        member, A, b_slice, delta_slice, handle);
    member.team_barrier();

    // Unscale delta before returning
    Kokkos::parallel_for(Kokkos::TeamVectorRange(member, system_size), [&](int j){
        delta_i(j) /= ewt_i(j);
    });
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
    policy_type policy; 
    policy_type norm_policy;
    //real_type_2d_view_type f_y0;
    ordinal_type numiters = 0;
    ordinal_type gmres_max_iter;
    //real_type gmres_tol; 
    ordinal_type system_size;
    ordinal_type n_systems; // nBatch
    ordinal_type per_team_extent = 0;
    ordinal_type last_flag;
    real_type resnorm;

    // Diagnostic parameters for dumping state at a given solve
    long nsolves = 0;        
    std::string dump_file = "";
    long        dump_solve = -1;      // which solve index to capture, <0 disables
    real_type   dump_rtol = 0, dump_atol = 0;
    bool        dumped = false;
    bool verbose = false;    // diagnostic printing; off by default, set from the driver

    void* cvode_mem;

    N_Vector s1 = nullptr;

    // constructor
    SUNLinearSolverContent_BatchedGMRES(Impl::AerosolChemistryRHS<problem_type> rhs_object_, 
                 ordinal_type system_size_, 
                 ordinal_type n_systems_,
                 ordinal_type gmres_max_iter_,
                 //real_type gmres_tol_,
                 policy_type policy_
                ) :
                handle(n_systems_, /*systems_per_team=*/1, gmres_max_iter_, /*monitor_residual=*/true),
                rhs_object(rhs_object_),
                system_size(system_size_),
                n_systems(n_systems_),
                gmres_max_iter(gmres_max_iter_),
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

        norm_policy = policy_type(n_systems, Kokkos::AUTO());


        ordinal_type scratch_level = 1;
        // ---- RHS scratch memory -----        
        // [ pw | y_shift | f_shift ]
        per_team_extent = TChem::Impl::Aerosol_RHS<real_type, device_type>::getWorkSpaceSize(rhs_object.kmcd, rhs_object.amcd) + 2*system_size;
        const size_t per_team_rhs_scratch =
        TChem::Scratch<real_type_1d_view_type>::shmem_size(per_team_extent);

        // ---- GMRES scratch memory ----
        // originally used level 0 for GMRES scratch, but with large systems (n_particles > 100), the scratch allocation exceeds shared memory (~48 KB for L40S)
        ordinal_type systems_per_team = 1; // each team works with 1 chemical system (key assumption of this solver)
        const ordinal_type maxit = (gmres_max_iter < system_size) ? gmres_max_iter : system_size; // for small systems with size < gmres_max_iter, we only need at most system_size # of iterations
        const size_t per_team_gmres_scratch  = real_type_2d_view_type::shmem_size(systems_per_team, (maxit + 1) + system_size + 2); // Givens rotation history + working vector + mask + tmp
        
        policy.set_scratch_size(scratch_level, Kokkos::PerTeam(per_team_rhs_scratch + per_team_gmres_scratch));

        // Configure the Krylov handle and Arnoldi view
        handle.set_scratch_pad_level(1); // use global scratch memory
        //handle = KrylovHandleType(n_systems, systems_per_team, gmres_max_iter, /*monitor_residual=*/true);
        handle.set_ortho_strategy(1); // use modified gram-schmidt for consistency with SUNDIALS
        //handle.set_tolerance(gmres_tol);
        handle.set_compute_last_residual(true);
        handle.Arnoldi_view = KrylovHandleType::ArnoldiViewType(
            "Arnoldi_view", n_systems, gmres_max_iter, system_size + gmres_max_iter + 3);

        if (verbose) {
          const double arnoldi_mb = double(n_systems) * double(gmres_max_iter)
                                  * double(system_size + gmres_max_iter + 3)
                                  * double(sizeof(real_type)) / (1024.0*1024.0);
          printf("..BatchedGMRES attached: n_systems = %d, system_size = %d, max_iter = %d "
                 "(GMRES clamps to min(max_iter, m) = %d)\n",
                 n_systems, system_size, gmres_max_iter, maxit);
          printf("..BatchedGMRES scratch : RHS+operator = %d reals/team (%zu B), "
                 "GMRES _TMPView = %zu B/team, total level-1 = %zu B/team\n",
                 per_team_extent, (size_t)per_team_rhs_scratch, (size_t)per_team_gmres_scratch,
                 (size_t)per_team_rhs_scratch + (size_t)per_team_gmres_scratch);
          printf("..BatchedGMRES handle scratch_pad_level = %d (must be 1), ortho_strategy = %d "
                 "(1 = modified GS, matching SPGMR)\n",
                 handle.get_scratch_pad_level(), handle.get_ortho_strategy());
          printf("..BatchedGMRES Arnoldi_view = (%d, %d, %d) = %.2f MB\n",
                 n_systems, gmres_max_iter, system_size + gmres_max_iter + 3, arnoldi_mb);
        }

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
        real_type_2d_view_type x_2d(N_VGetDeviceArrayPointer(x), n_systems, system_size);
        real_type_2d_view_type ewt_2d(N_VGetDeviceArrayPointer(s1), n_systems, system_size);


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

        handle.set_tolerance(tol/bnorm_max);
        handle.reset(); // reset the iteration numbers and residual norms

        // Create GMRES functor
        MatrixFreeBatchedGMRESTeamFunctor<Impl::AerosolChemistryRHS<problem_type>> functor(rhs_object, y0_2d, x_2d, fy0_2d, b_2d, 
        gamma, ewt_2d, handle, per_team_extent, verbose, nsolves);
        
        // GMRES solve in parallel over teams
        Kokkos::parallel_for(policy, functor);
        Kokkos::fence();
        
        // Evaluate diagnostics for GMRES batch solve 
        numiters = 0;
        resnorm = 0;
        bool all_conv = true;
        for (ordinal_type i = 0; i < n_systems; i++){
            numiters = std::max(numiters, handle.get_iteration_host(i)); // TODO mixing int with whatever ordinal_type is.. could be an issue--exlicitly set numiters as type int?
            all_conv &= handle.is_converged_host(i); // is converged just tells us that we exit before exhausting # of iterations, not that we actually converge to the tolerance
            resnorm = std::max(resnorm, handle.get_last_norm_host(i));
        }

        if (all_conv && (resnorm <= tol)) {
            last_flag = SUN_SUCCESS;
        } else {
            last_flag = SUNLS_RES_REDUCED;
        }

        // Solver diagnostics:
        // nsolves: Linear solver iteration 
        // numiters: MAX number of iterations required by GMRES across the batch of systems
        // resnorm: residual Ax* - b~ for the scaled system with candidate solution x*
        // bnorm_max: The maximum 2-norm of linear system RHSs across all systems in the batch (used to set GMRES tolerance)
        // reduction: resnorm/||b~||, <<1 for successful linear solves; if ~ 1, GMRES was not able to converge on an adequate solution
        // tol: CVODE absolute tolerance
        // tol/bnorm_max: relative tolerance used by GMRES
        // n_maxed: number of systems GMRES failed to converge before it hit the max number of Arnoldi iterations
        // n_systems: number of independent systems in the batch 
        if (verbose) {
          ordinal_type n_maxed = 0;
          for (ordinal_type i = 0; i < n_systems; i++) {
            if (!handle.is_converged_host(i)) ++n_maxed;
          }
          const real_type reduction = (bnorm_max > 0) ? resnorm/bnorm_max : 0.0;
          printf("  [bgmr] solve %ld: iters(max) = %d, resnorm(max) = %.3e, ||b~|| = %.3e, "
                 "reduction = %.3e, tol = %.3e (rel = %.3e), not converged = %d/%d%s\n",
                 nsolves, numiters, resnorm, bnorm_max, reduction, tol, tol/bnorm_max,
                 n_maxed, n_systems,
                 (last_flag == SUN_SUCCESS) ? "" : "   <-- FAIL");
        }
        nsolves++;

        return last_flag;
    }

};

SUNLinearSolver SUNLinSol_BatchedGMRES(N_Vector y, Impl::AerosolChemistryRHS<SUNLinearSolverContent_BatchedGMRES::problem_type> rhs_object,
                                       ordinal_type system_size, ordinal_type n_systems, ordinal_type gmres_max_iter,
                                       SUNLinearSolverContent_BatchedGMRES::policy_type policy,
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
int SUNLinSolSetup_BGMR(SUNLinearSolver S, SUNMatrix A);
int SUNLinSolSolve_BGMR(SUNLinearSolver S, SUNMatrix A,
                         N_Vector x, N_Vector b, sunrealtype delta);
int SUNLinSolNumIters_BGMR(SUNLinearSolver S);
sunrealtype SUNLinSolResNorm_BGMR(SUNLinearSolver S);
sunindextype SUNLinSolLastFlag_BGMR(SUNLinearSolver S);
SUNErrCode SUNLinSolFree_BGMR(SUNLinearSolver S);

} // namespace TChem