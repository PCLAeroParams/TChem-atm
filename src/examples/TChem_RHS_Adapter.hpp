#pragma once
// NOTE: this file is TEMPORARILY located in src/examples for validating the struct
// appropriately calls the RHS evaluation method and properly creates per-team slices
// of the problem state.
//
// Adapter presenting TChem's aerosol chemistry RHS to the JacvecFDTeam RHSFunctor
// contract. Templated on the TChem problem type
// (TChem::Impl::AerosolChemistry_Problem<real_type, device_type>); concrete view and
// model-data types are pulled from it so we stay in sync with TChem.
#include <Kokkos_Core.hpp>
#include "TChem.hpp"   // TChem::Impl::StateVector, Tines views, ordinal_type, etc.

template <typename ProblemType>
struct AerosolChemistryRHS {

    // Types pulled from the TChem problem type (matches what the driver/CVODE provide).
    using real_type               = typename ProblemType::real_type;
    using ordinal_type            = TChem::ordinal_type;
    using real_type_1d_view_type  = typename ProblemType::real_type_1d_view_type;
    using real_type_2d_view_type  = typename ProblemType::real_type_2d_view_type;
    using kinetic_model_type      = typename ProblemType::kinetic_model_type;
    using aerosol_model_data_type = typename ProblemType::aerosol_model_data_type;

    // members: the batched form holds the whole multi-system batch; the sliced form
    // (returned by team_slice) additionally carries a per-system-configured `problem`.
    real_type_2d_view_type  state;        // full state vectors, batched (nBatch, stateVecDim)
    real_type_2d_view_type  number_conc;  // aerosol number concentration, batched
    kinetic_model_type      kmcd;         // shared kinetic model const data
    aerosol_model_data_type amcd;         // shared aerosol model const data
    ProblemType             problem;      // configured per-system by team_slice()

    // constructor (batched form; `problem` left default-constructed)
    KOKKOS_INLINE_FUNCTION
    AerosolChemistryRHS(real_type_2d_view_type state_, real_type_2d_view_type number_conc_,
                        kinetic_model_type kmcd_, aerosol_model_data_type amcd_)
        : state(state_), number_conc(number_conc_), kmcd(kmcd_), amcd(amcd_) {}

    // Per-system slice: decode system `first`'s parameters from the full state and
    // configure a Problem for it. `last` is unused (one system per team); it is kept for
    // signature symmetry with LaplacianRHS::team_slice. `work` is the per-team RHS scratch
    // (size problem_type::getWorkSpaceSize), allocated by the driver's kernel.
    KOKKOS_INLINE_FUNCTION
    AerosolChemistryRHS team_slice(int first, int last, real_type_1d_view_type work) const {
        (void)last;
        const ordinal_type total_n_species = kmcd.nSpec + amcd.nParticles*amcd.nSpec;
        auto state_at_i = Kokkos::subview(state, first, Kokkos::ALL());
        TChem::Impl::StateVector<real_type_1d_view_type> sv_at_i(total_n_species, state_at_i);
        const auto T_i = sv_at_i.Temperature();
        const auto P_i = sv_at_i.Pressure();
        auto Ys = sv_at_i.MassFractions();
        const ordinal_type n_active_gas_species = kmcd.nSpec - kmcd.nConstSpec;
        auto constYs_i = Kokkos::subview(Ys, Kokkos::make_pair(n_active_gas_species, kmcd.nSpec));
        auto number_conc_at_i = Kokkos::subview(number_conc, first, Kokkos::ALL());

        // configure the per-system problem on a shallow copy of *this
        AerosolChemistryRHS sliced = *this;
        sliced.problem._kmcd                = kmcd;
        sliced.problem._amcd                = amcd;
        sliced.problem._temperature         = T_i;
        sliced.problem._pressure            = P_i;
        sliced.problem._number_conc         = number_conc_at_i;
        sliced.problem._const_concentration = constYs_i;
        sliced.problem._work                = work;
        return sliced;
    }

    // team-collective whole-vector RHS: fills f_out with f(x) using the whole team.
    template <typename MemberType, typename InView, typename OutView>
    KOKKOS_INLINE_FUNCTION
    void operator()(const MemberType& member, const InView& x, const OutView& f_out) const {
        problem.computeFunction(member, x, f_out);
    }

};
