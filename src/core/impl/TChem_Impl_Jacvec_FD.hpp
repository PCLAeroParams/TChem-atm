#pragma once
// Matrix-free Jacobian-vector operator for batched GMRES (one system / team)
#include <Kokkos_Core.hpp>
#include <Kokkos_ArithTraits.hpp>

namespace TChem {
namespace Impl {

template <typename RHSFunctor>
struct JacvecFDTeam{

    // typedefs
    using real_type = typename RHSFunctor::real_type;
    using real_type_1d_view_type = typename RHSFunctor::real_type_1d_view_type;
    using MagnitudeType = typename Kokkos::ArithTraits<real_type>::mag_type;

    // members (note that instead of 2D views we pass in 1D since systems_per_team=1 always)
    RHSFunctor f;                       // device-callable for returning the RHS given a state x
    real_type_1d_view_type  x0;         // system linearization point
    real_type_1d_view_type f_x;         // cached f(x0) (computed prior to passing to struct)
    real_type gamma;                    // a constant specific to the CVODE method 
    real_type_1d_view_type ewt;         // system error weights (used for computing the WRMS-scaled difference quotient Jv)
    real_type_1d_view_type x_shift;     // x_0 + sigma*x (note stored at level-1 scratch, allocated prior to calling JacvecFD's apply method)
    real_type_1d_view_type f_shift;     // similar to the approach for x_shift where we allocate prior to calling apply

    // constructor
    KOKKOS_INLINE_FUNCTION 
    JacvecFDTeam(RHSFunctor f_, real_type_1d_view_type x0_, real_type_1d_view_type f_x_, real_type gamma_, 
        real_type_1d_view_type ewt_, real_type_1d_view_type x_shift_, real_type_1d_view_type f_shift_)
        : f(f_), x0(x0_), f_x(f_x_), gamma(gamma_), ewt(ewt_), x_shift(x_shift_), f_shift(f_shift_) {}

    template <typename ArgTrans, typename ArgMode, typename MemberType, typename XViewType, typename YViewType>
    KOKKOS_INLINE_FUNCTION void apply(const MemberType &member, const XViewType &X, const YViewType &Y, 
        MagnitudeType alpha = 1, MagnitudeType beta = 0) const{

            const int n = x0.extent(0); // the size of each system

            // Allocate 1D subviews for X and Y (GMRES natively expects 2D, so GMRES calls apply with the 2D views (1, n) and internally we downgrade to 1D views)
            auto x = Kokkos::subview(X, 0, Kokkos::ALL);
            auto y = Kokkos::subview(Y, 0, Kokkos::ALL);

            // compute the WRMS norm of w 
            // Note ||W||_WRMS = sqrt(1/N sum(w_i*ewt_i)^2) = sqrt(1/N sum((x_i/ewt_i)*ewt_i)^2) = ||x||_2 / sqrt(N) for x a vector with unit 2-norm
            real_type sumsq = 0; 
            Kokkos::parallel_reduce(Kokkos::TeamVectorRange(member, 0, n), [&](const int& i, real_type& lsum){
                lsum += x(i)*x(i); // sum[(x_i)^2]
            }, sumsq);
            member.team_barrier();
            real_type nrm_w = Kokkos::ArithTraits<real_type>::sqrt(sumsq/n);

            // catch NaNs or if x = 0 (intitial condition for GMRES) to prevent sigma from blowing up
            // we need to include this catch here because KokkosKernels TeamVectorGMRES does not check for the zero vector before 
            // calling this apply() to get the matrix-vector product
            if (!(nrm_w > 0)){
                Kokkos::parallel_for(Kokkos::TeamVectorRange(member, 0, n), [&](const int& i ){
                    y(i) = beta*y(i); // Return b
                });
                return;
            }

            real_type sigma = 1/nrm_w;

            // Team vector range to fill in x_shift.
            Kokkos::parallel_for(Kokkos::TeamVectorRange(member, 0, n), [&](const int& i ){
                const real_type w_i =  (1/ewt(i))*x(i); // diag(ewt_i)^-1 * x
                //x_shift(i) = x0(i) + sigma*x(i);
                x_shift(i) = x0(i) + sigma*w_i;
            });
            member.team_barrier();
            
            // overload operator() so we can pass f() to the underlying TChem problem computeFunction(). Note that we cannot 'return' f_shift since that reassigns the member view handle (not allowed since apply is const), and instead must write to a buffer
            f(member, x_shift, f_shift);  // (member, in, out)
            member.team_barrier(); // TChem's RHS does include a barrier at its end, but this is included here for generality

            // TeamVectorRange to fill in Jv_prod (which can be allocated within the buffer of the parallel_for) and then immediately also fill in Y
            Kokkos::parallel_for(Kokkos::TeamVectorRange(member, 0, n), [&](const int& i ){
                const real_type Jv_prod = ( f_shift(i) - f_x(i) ) / (sigma);
                // Now we can compute y with whatever values of alpha and beta are passed in
                const real_type w_i =  (1/ewt(i))*x(i);
                y(i) = beta*y(i) + alpha*ewt(i)*(w_i - gamma*Jv_prod);
            });

    } // apply()

}; // struct JacvecFDTeam

} // namespace Impl
} // namespace TChem
