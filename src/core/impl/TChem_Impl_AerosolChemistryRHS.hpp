#pragma once
// AerosolChemistryRHS constructs per-system problems via AerosolChemistryRHS::team_slice()
// and the corresponding system RHS f(x) can be evaluated via its operator() method.
// This is used by the BatchedGMRES solver and its matrix-free method (JacvecFDTeam) to 
// construct linear systems for each team.
#include <Kokkos_Core.hpp>
#include "TChem.hpp"   // TChem::Impl::StateVector, Tines views, ordinal_type, etc.
#include "TChem_Impl_Aerosol_RHS.hpp"

namespace TChem {
namespace Impl {

template <typename ProblemType>
struct AerosolChemistryRHS {

    // typedefs
    using device_type               = typename ProblemType::device_type;
    using real_type                 = typename ProblemType::real_type;
    using ordinal_type              = TChem::ordinal_type;
    using real_type_1d_view_type    = typename ProblemType::real_type_1d_view_type;
    using real_type_2d_view_type    = typename ProblemType::real_type_2d_view_type;
    using kinetic_model_type        = typename ProblemType::kinetic_model_type;
    using aerosol_model_data_type   = typename ProblemType::aerosol_model_data_type;
    using ordinal_type_1d_view_type = Tines::value_type_1d_view<ordinal_type,device_type>;

    // members: the batched form holds the whole multi-system batch; the sliced form
    // (returned by team_slice) additionally carries a per-system-configured `problem`.
    /*
    real_type_2d_view_type  state;        // full state vectors, batched (nBatch, stateVecDim)
    real_type_2d_view_type  number_conc;  // aerosol number concentration, batched
    kinetic_model_type      kmcd;         // shared kinetic model const data
    aerosol_model_data_type amcd;         // shared aerosol model const data
    ProblemType             problem;      // configured per-system by team_slice()
    */
    // batch-level member views
    real_type_2d_view_type          num_concentration;
    real_type_2d_view_type          const_tracers;
    real_type_1d_view_type          temperature;
    real_type_1d_view_type          pressure;
    kinetic_model_type              kmcd;
    aerosol_model_data_type         amcd;
    ordinal_type_1d_view_type       n_particles_track;

    // team-level member subviews / scalar values
    real_type_1d_view_type          num_concentration_at_i;
    real_type_1d_view_type          const_tracers_at_i;
    real_type                       temperature_at_i;
    real_type                       pressure_at_i;
    ordinal_type                    n_particles_track_at_i;
    real_type_1d_view_type          work_at_i;
    bool slice_flag;

    // constructor (batched form; `problem` left default-constructed)
    /*
    KOKKOS_INLINE_FUNCTION
    AerosolChemistryRHS(real_type_2d_view_type state_, real_type_2d_view_type number_conc_,
                        kinetic_model_type kmcd_, aerosol_model_data_type amcd_)
        : state(state_), number_conc(number_conc_), kmcd(kmcd_), amcd(amcd_) {}
    */
    AerosolChemistryRHS(
           const real_type_2d_view_type& num_concentration_in,
           const real_type_2d_view_type& const_tracers_in,
           const real_type_1d_view_type& temperature_in,
           const real_type_1d_view_type& pressure_in,
           const kinetic_model_type& kmcd_in,
           const aerosol_model_data_type& amcd_in,
           const ordinal_type_1d_view_type& n_particles_track_in)
         : num_concentration(num_concentration_in),
           const_tracers(const_tracers_in),
           temperature(temperature_in),
           pressure(pressure_in),
           kmcd(kmcd_in),
           amcd(amcd_in),
           n_particles_track(n_particles_track_in) {
            slice_flag = false; // flag for catching operator() called on unsliced version of AerosolChemistryRHS
    }

    // Per-system slice
    KOKKOS_INLINE_FUNCTION
    AerosolChemistryRHS team_slice(int i_member, real_type_1d_view_type work) const { // note I am deprecating last here and will need to also update call sites
        
        auto num_conc_at_i = Kokkos::subview(num_concentration, i_member, Kokkos::ALL());
        auto constYs_i = Kokkos::subview(const_tracers, i_member, Kokkos::ALL());
        const ordinal_type n_part_i =
            (n_particles_track.span() == 0 || n_particles_track(i_member) < 0) 
            ? amcd.nParticles : n_particles_track(i_member);

        // configure the per-system problem on a shallow copy of *this
        // allows us to return a member slice of each view and call operator() to evaluate the member's RHS
        AerosolChemistryRHS sliced = *this; // already assigns amcd, kmcd
        sliced.num_concentration_at_i       = num_conc_at_i;
        sliced.const_tracers_at_i           = constYs_i;
        sliced.temperature_at_i             = temperature(i_member);
        sliced.pressure_at_i                = pressure(i_member);
        sliced.n_particles_track_at_i       = n_part_i;
        sliced.work_at_i                    = work; 
        sliced.slice_flag = true;
        return sliced;
    }

    // team-collective whole-vector RHS: fills f_out with f(x) using the whole team.
    template <typename MemberType, typename InView, typename OutView>
    KOKKOS_INLINE_FUNCTION
    void operator()(const MemberType& member, const InView& x, const OutView& f_out) const {
        //problem.computeFunction(member, x, f_out);

        // check if trying to call operator on non-sliced AerosolChemistryRHS
        if (!slice_flag){
            Kokkos::abort("Calling AerosolChemistryRHS::operator() with unsliced views");
        }

        TChem::Impl::Aerosol_RHS<real_type, device_type>::team_invoke(
            member, temperature_at_i, pressure_at_i, num_concentration_at_i, 
            x, const_tracers_at_i, f_out, work_at_i, kmcd, amcd, n_particles_track_at_i);
        member.team_barrier(); 
    }

};

} // namespace Impl
} // namespace TChem