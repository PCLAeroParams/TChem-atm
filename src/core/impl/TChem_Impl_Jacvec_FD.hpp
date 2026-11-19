#pragma once
#include <Kokkos_Core.hpp>
#include <Kokkos_ArithTraits.hpp>

namespace TChem {
namespace Impl {

// Struct for computing the linear system solved by Kokkos Kernels Batched GMRES
// via matrix-free difference quotient evaluation of the Jacobian-vector product.
// The struct apply() operates on individual system slices of the batch that KK BGMR 
// solves. To compute the difference quotient, the perturbation step size sigma is determined 
// in the same manner as SUNDIALS, whereby sigma = 1 / || u ||_WRMS where the WRMS norm 
// weights u by the error weights of x provided by CVODE. 
// The system is first unscaled and then preconditioned by the templated preconditioner 
// before computing the Jacobian-vector product and returning the re-scaled linear system
// LHS I-gamma*Ju in Y.  
template <typename RHSType, typename PreconType>
struct JacvecFDTeam{

    // typedefs
    using real_type = typename RHSType::real_type;
    using real_type_1d_view_type = typename RHSType::real_type_1d_view_type;
    using MagnitudeType = typename Kokkos::ArithTraits<real_type>::mag_type;

    // members
    RHSType f;                          // System-level RHS function
    PreconType precon;                  // Preconditioner (Identity or Schur)
    real_type_1d_view_type  x0;         // system linearization point
    real_type_1d_view_type f_x;         // cached f(x0), calculated by CVODE
    real_type gamma;                    // CVODE integrator constant
    real_type_1d_view_type ewt;         // system error weights for WRMS-based DQ step size
    real_type_1d_view_type x_shift;     // perturbation around the linearization point
    real_type_1d_view_type f_shift;     // RHS at perturbed state
    real_type_1d_view_type u;           // buffer for preconditioned solve input/output

    // Diagnostic parameters
    bool diag = false;
    long diag_solve = -1;               // solve index, used only to label the print

    // constructor
    KOKKOS_INLINE_FUNCTION
    JacvecFDTeam(RHSType f_, PreconType precon_, real_type_1d_view_type x0_, 
        real_type_1d_view_type f_x_, real_type gamma_, real_type_1d_view_type ewt_, 
        real_type_1d_view_type x_shift_, real_type_1d_view_type f_shift_,
        real_type_1d_view_type u_, bool diag_ = false, long diag_solve_ = -1)
        : f(f_), precon(precon_), x0(x0_), f_x(f_x_), gamma(gamma_), ewt(ewt_), 
          x_shift(x_shift_), f_shift(f_shift_), u(u_), diag(diag_), 
          diag_solve(diag_solve_) {}
    
    template <typename ArgTrans, typename ArgMode, typename MemberType, typename XViewType, typename YViewType>
    KOKKOS_INLINE_FUNCTION void apply(const MemberType &member, const XViewType &X, const YViewType &Y, 
        MagnitudeType alpha = 1, MagnitudeType beta = 0) const{

            const int system_size = x0.extent(0); // the size of each system

            // Allocate 1D subviews for X and Y (passed as 1 x system_size 2D views)
            auto x = Kokkos::subview(X, 0, Kokkos::ALL);
            auto y = Kokkos::subview(Y, 0, Kokkos::ALL);

            real_type sum_xsqr = 0;
            Kokkos::parallel_reduce(Kokkos::TeamVectorRange(member, 0, system_size), 
              [&](const int& i, real_type& xsum){
                xsum += x(i)*x(i);
            }, sum_xsqr);
            member.team_barrier();

            // initial condition, return r_0 = b - A*x_0 where x_0 = 0
            if (!(sum_xsqr > 0)){
                // Performance note: this is redundant with beta = 1, but a GMRES breakdown
                // could cause X on next call to apply to be zero
                Kokkos::parallel_for(Kokkos::TeamVectorRange(member, 0, system_size), [&](const int& i ){
                    y(i) = beta*y(i);
                });
                return;
            }

            // unscale x into u
            Kokkos::parallel_for(Kokkos::TeamVectorRange(member, 0, system_size), [&](const int& i ){
                u(i) = (1/ewt(i))*x(i); // diag(ewt_i)^-1 * x
            });
            member.team_barrier();

            // Solve the preconditioned system Pu* = u where u is returned in-place
            // (just symbolically referring to u* as the solution to clarify intent of 
            // this operation)
            precon.team_solve(member, u);
            member.team_barrier(); // NOTE could place barrier inside team solve or here..

            // compute the WRMS norm of u
            // Note ||u||_WRMS = sqrt(1/N sum(u_i*ewt_i)^2) = ||u \dot ewt ||_2 / sqrt(N)
            real_type sumsq = 0; 
            Kokkos::parallel_reduce(Kokkos::TeamVectorRange(member, 0, system_size), 
              [&](const int& i, real_type& lsum){
                lsum += (u(i)*ewt(i))*(u(i)*ewt(i)); 
            }, sumsq);
            member.team_barrier();
            real_type nrm_u = Kokkos::ArithTraits<real_type>::sqrt(sumsq/system_size);

            // if NaNs in u, exit early to prevent sigma from blowing up
            if (!(nrm_u > 0)){
                if (diag) {
                    Kokkos::single(Kokkos::PerTeam(member), [&]() {
                        Kokkos::printf("  [bgmr-diag] solve %ld: apply(alpha=%g, beta=%g) early return, nrm_u = %e\n",
                                       diag_solve, (double)alpha, (double)beta, (double)nrm_u);
                    });
                }
                Kokkos::parallel_for(Kokkos::TeamVectorRange(member, 0, system_size), [&](const int& i ){
                    y(i) = beta*y(i);
                });
                return;
            }

            real_type sigma = 1/nrm_u;

            // Fill in x_shift.
            Kokkos::parallel_for(Kokkos::TeamVectorRange(member, 0, system_size), [&](const int& i ){
                x_shift(i) = x0(i) + sigma*u(i);
            });
            member.team_barrier();
            
            // overload operator() so we can pass f() to the underlying TChem problem 
            // computeFunction(). Note that we cannot 'return' f_shift since that 
            // reassigns the member view handle (not allowed since apply is const), and 
            // instead must write to a buffer
            f(member, x_shift, f_shift);  // (member, in, out)
            member.team_barrier(); // TChem RHS does have barrier; including as safeguard against future changes

            /*
            // ---- Diagnostic breakdown of this Jv product ----
            // Scaled space: x is the Krylov vector passed in by GMRES, and for alpha = 1, beta = 0 this call returns
            //   Ax = ewt*(w - gamma*Jv) = x - gamma*ewt*Jv.
            // For the first Arnoldi product (x = v1, ||x|| = 1), sin(Ax, x) = h21/sqrt(h11^2 + h21^2), which is
            // exactly GMRES's one-iteration relative residual estimate. Read-only: nothing below changes the result.
            if (diag) {
                using ATS = Kokkos::ArithTraits<real_type>;
                real_type fx2 = 0, df2 = 0, g2 = 0, Ax2 = 0, Axx = 0, max_rel_pert = 0;
                int n_same = 0, n_lost = 0, n_neg = 0, n_bad = 0;

                Kokkos::parallel_reduce(Kokkos::TeamVectorRange(member, 0, system_size), [&](const int& i, real_type& s){
                    s += f_x(i)*f_x(i);
                }, fx2);
                Kokkos::parallel_reduce(Kokkos::TeamVectorRange(member, 0, system_size), [&](const int& i, real_type& s){
                    const real_type d = f_shift(i) - f_x(i);
                    s += d*d;
                }, df2);
                Kokkos::parallel_reduce(Kokkos::TeamVectorRange(member, 0, system_size), [&](const int& i, real_type& s){
                    const real_type g = gamma*ewt(i)*(f_shift(i) - f_x(i))/sigma;
                    s += g*g;
                }, g2);
                Kokkos::parallel_reduce(Kokkos::TeamVectorRange(member, 0, system_size), [&](const int& i, real_type& s){
                    const real_type Ax = x(i) - gamma*ewt(i)*(f_shift(i) - f_x(i))/sigma;
                    s += Ax*Ax;
                }, Ax2);
                Kokkos::parallel_reduce(Kokkos::TeamVectorRange(member, 0, system_size), [&](const int& i, real_type& s){
                    const real_type Ax = x(i) - gamma*ewt(i)*(f_shift(i) - f_x(i))/sigma;
                    s += Ax*x(i);
                }, Axx);
                // components where the RHS did not change at all
                Kokkos::parallel_reduce(Kokkos::TeamVectorRange(member, 0, system_size), [&](const int& i, int& c){
                    if (f_shift(i) == f_x(i)) c += 1;
                }, n_same);
                // components where a nonzero perturbation was lost to rounding (x0 + sigma*w == x0)
                Kokkos::parallel_reduce(Kokkos::TeamVectorRange(member, 0, system_size), [&](const int& i, int& c){
                    if (x(i) != 0 && x_shift(i) == x0(i)) c += 1;
                }, n_lost);
                // components pushed negative by the perturbation
                Kokkos::parallel_reduce(Kokkos::TeamVectorRange(member, 0, system_size), [&](const int& i, int& c){
                    if (x_shift(i) < 0 && !(x0(i) < 0)) c += 1;
                }, n_neg);
                Kokkos::parallel_reduce(Kokkos::TeamVectorRange(member, 0, system_size), [&](const int& i, int& c){
                    if (!Kokkos::isfinite(f_shift(i))) c += 1;
                }, n_bad);
                // largest perturbation actually applied, relative to the state
                Kokkos::parallel_reduce(Kokkos::TeamVectorRange(member, 0, system_size), [&](const int& i, real_type& m){
                    if (x0(i) != 0) {
                        const real_type r = ATS::abs(x_shift(i) - x0(i))/ATS::abs(x0(i));
                        if (r > m) m = r;
                    }
                }, Kokkos::Max<real_type>(max_rel_pert));
                member.team_barrier();

                Kokkos::single(Kokkos::PerTeam(member), [&]() {
                    const real_type x2  = sumsq;                  // ||x||^2
                    const real_type h11 = Axx/x2;                 // <Ax, x>/||x||^2
                    real_type h21sq = Ax2 - Axx*Axx/x2;           // ||Ax - h11*x||^2
                    if (h21sq < 0) h21sq = 0;
                    const real_type sin_th = (Ax2 > 0) ? ATS::sqrt(h21sq/Ax2) : real_type(0);
                    Kokkos::printf("  [bgmr-diag] solve %ld: apply(alpha=%g, beta=%g) system_size=%d gamma=%.6e sigma=%.6e ||x||=%.6e\n",
                                   diag_solve, (double)alpha, (double)beta, system_size, (double)gamma, (double)sigma,
                                   (double)ATS::sqrt(x2));
                    Kokkos::printf("  [bgmr-diag]    ||f_x||=%.6e  ||f_shift-f_x||=%.6e (rel %.3e)  ||gamma*ewt*Jv||=%.6e  ||Ax||=%.6e\n",
                                   (double)ATS::sqrt(fx2), (double)ATS::sqrt(df2),
                                   (double)(fx2 > 0 ? ATS::sqrt(df2/fx2) : real_type(0)),
                                   (double)ATS::sqrt(g2), (double)ATS::sqrt(Ax2));
                    Kokkos::printf("  [bgmr-diag]    h11=%.6e  h21/||x||=%.6e  sin(Ax,x)=%.6e\n",
                                   (double)h11, (double)ATS::sqrt(h21sq/x2), (double)sin_th);
                    Kokkos::printf("  [bgmr-diag]    f_shift==f_x in %d/%d, perturbation lost to rounding in %d, "
                                   "new negatives %d, non-finite f_shift %d, max |dx|/|x0| = %.3e\n",
                                   n_same, system_size, n_lost, n_neg, n_bad, (double)max_rel_pert);
                });
                member.team_barrier();
            }
            */

            Kokkos::parallel_for(Kokkos::TeamVectorRange(member, 0, system_size), [&](const int& i ){
                const real_type Ju = ( f_shift(i) - f_x(i) ) / (sigma);
                y(i) = beta*y(i) + alpha*ewt(i)*(u(i) - gamma*Ju);
            });

    } // apply()

}; // struct JacvecFDTeam

} // namespace Impl
} // namespace TChem
