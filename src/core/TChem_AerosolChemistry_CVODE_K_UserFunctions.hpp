/* =====================================================================================
TChem-atm version 2.0.0
Copyright (2025) NTESS
https://github.com/sandialabs/TChem-atm

Copyright 2025 National Technology & Engineering Solutions of Sandia, LLC
(NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
Government retains certain rights in this software.

This file is part of TChem-atm. TChem-atm is open source software: you can redistribute
it and/or modify it under the terms of BSD 2-Clause License
(https://opensource.org/licenses/BSD-2-Clause). A copy of the licese is also
provided under the main directory

Questions? Contact Oscar Diaz-Ibarra at <odiazib@sandia.gov>, or
           Cosmin Safta at <csafta@sandia.gov> or,
           Nicole Riemer at <nriemer@illinois.edu> or,
           Matthew West at <mwest@illinois.edu>

Sandia National Laboratories, New Mexico/Livermore, NM/CA, USA
=====================================================================================
*/
#ifndef __TCHEM_AEROSOL_CHEMISTRY_CVODE_K_USER_FUNCTIONS_HPP__
#define __TCHEM_AEROSOL_CHEMISTRY_CVODE_K_USER_FUNCTIONS_HPP__

#include "TChem_Util.hpp"
#include "TChem_KineticModelData.hpp"
#include "TChem_Impl_AerosolChemistry.hpp"
#include "TChem_Impl_Aerosol_RHS.hpp"
#include <cstdio>
#include <cvode/cvode.h>
#include <memory>
#include <nvector/nvector_kokkos.hpp>
#include <sundials/sundials_core.hpp>
#include <sunlinsol/sunlinsol_kokkosdense.hpp>
#include <sunlinsol/sunlinsol_spgmr.h>
#include <sunmatrix/sunmatrix_kokkosdense.hpp>
#include <vector>
#include <KokkosBatched_LU_Decl.hpp>
#include <KokkosBatched_SolveLU_Decl.hpp>

namespace TChem {


struct TChemAerosolChemistryRHS {

  using device_type      = typename Tines::UseThisDevice<exec_space>::type;
  using problem_type = TChem::Impl::AerosolChemistry_Problem<real_type,device_type>;
  using policy_type = typename UseThisTeamPolicy<exec_space>::type;
  using real_type_1d_view_type = Tines::value_type_1d_view<real_type,device_type>;
  using real_type_2d_view_type = Tines::value_type_2d_view<real_type,device_type>;
  using ordinal_type_1d_view_type = Tines::value_type_1d_view<ordinal_type,device_type>;

  real_type_2d_view_type rhs;
  real_type_2d_view_type vals;
  real_type_2d_view_type aero_num_conc;
  real_type_2d_view_type const_tracers;
  real_type_1d_view_type temperature;
  real_type_1d_view_type pressure;
  KineticModelNCAR_ConstData<device_type> kmcd;
  AerosolModel_ConstData<device_type> amcd;
  ordinal_type level;
  ordinal_type per_team_extent;
  ordinal_type m;
  // Per-batch number of particles to evaluate in the RHS (length nBatch)
  ordinal_type_1d_view_type n_particles_track;


  TChemAerosolChemistryRHS(const real_type_2d_view_type& rhs_in,
           const real_type_2d_view_type& vals_in,
           const real_type_2d_view_type& aero_num_conc_in,
           const real_type_2d_view_type& const_tracers_in,
           const real_type_1d_view_type& temperature_in,
           const real_type_1d_view_type& pressure_in,
           const KineticModelNCAR_ConstData<device_type>& kmcd_in,
           const AerosolModel_ConstData<device_type>& amcd_in,
           const ordinal_type_1d_view_type& n_particles_track_in)
   : rhs(rhs_in),
     vals(vals_in),
     aero_num_conc(aero_num_conc_in),
     const_tracers(const_tracers_in),
     temperature(temperature_in),
     pressure(pressure_in),
     kmcd(kmcd_in),
     amcd(amcd_in),
     n_particles_track(n_particles_track_in)
     {
      level = 1;
      per_team_extent
         = Impl::Aerosol_RHS<real_type, device_type>::getWorkSpaceSize(kmcd, amcd);
      m = problem_type::getNumberOfTimeODEs(kmcd,amcd);
     }

  KOKKOS_INLINE_FUNCTION
  void operator()(const typename policy_type::member_type& member) const {
    const ordinal_type i = member.league_rank();

    const real_type_1d_view_type rhs_at_i =
    Kokkos::subview(rhs, i, Kokkos::ALL());
    const real_type_1d_view_type vals_at_i =
    Kokkos::subview(vals, i, Kokkos::ALL());
    const real_type_1d_view_type number_conc_at_i =
    Kokkos::subview(aero_num_conc, i, Kokkos::ALL());
    TChem::Scratch<real_type_1d_view_type> work(member.team_scratch(level),
                                   per_team_extent);
    const real_type_1d_view_type constYs  = Kokkos::subview(const_tracers, i, Kokkos::ALL());
    auto wptr = work.data();
    auto pw = real_type_1d_view_type(wptr, per_team_extent);
    wptr +=per_team_extent;
    // Resolve the per-batch sentinel: -1 means track all particles.
    // Check unallocated or -1 (default, all particles in a system)
    const ordinal_type n_part_i =
        (n_particles_track.span() == 0 || n_particles_track(i) < 0) 
        ? amcd.nParticles : n_particles_track(i);
    TChem::Impl::Aerosol_RHS<real_type, device_type>
    ::team_invoke(member,
    temperature(i), pressure(i), number_conc_at_i, vals_at_i, constYs,
    rhs_at_i, pw, kmcd, amcd, n_part_i);
  }
};

struct UserData
{
  using host_device_type = typename Tines::UseThisDevice<host_exec_space>::type;
  using device_type      = typename Tines::UseThisDevice<exec_space>::type;

  using real_type_1d_view_type = Tines::value_type_1d_view<real_type,device_type>;
  using real_type_2d_view_type = Tines::value_type_2d_view<real_type,device_type>;
  using real_type_3d_view_type = Tines::value_type_3d_view<real_type,device_type>;
  using real_type_4d_view_type = Tines::value_type_4d_view<real_type,device_type>;
  using ordinal_type_1d_view_type = Tines::value_type_1d_view<ordinal_type,device_type>;
  using policy_type = typename UseThisTeamPolicy<exec_space>::type;

  // TODO consider sizing these both to -1 to ensure attribs get set before integration?
  int n_systems  = 100; // number of independent chemical systems (e.g., grid cells, box models)
  int system_size = 3;  // size of each system

  policy_type policy;
  real_type_2d_view_type aero_num_conc;
  real_type_1d_view_type temperature;
  real_type_1d_view_type pressure;
  real_type_2d_view_type const_tracers;
  real_type_2d_view_type fac;
#if defined(TCHEM_ATM_ENABLE_GPU)
  real_type_3d_view_type JacRL;
#endif
  ordinal_type team_size_recommended_rhs{-1};
  ordinal_type team_size_recommended_jac{-1};
  ordinal_type team_size_max_rhs{-1};
  ordinal_type team_size_max_jac{-1};

  TChem::KineticModelNCAR_ConstData<device_type> kmcd;
  TChem::AerosolModel_ConstData<device_type> amcd;
  // Per-batch number of particles for RHS eval (length nBatch; entry -1 = all)
  ordinal_type_1d_view_type n_particles_track;

  // ---- Preconditioner data structures ----
  // TODO: will need to allocate these views in TChem_Driver.cpp when UserData is populated

  // Views for Jacobian blocks 
  real_type_3d_view_type Jgg;     // upper left    (gas-gas)
  real_type_4d_view_type Jgp;     // upper right   (gas-particle) [TODO allocate only first m rows!]
  real_type_4d_view_type Jpg;     // bottom left   (particle-gas) [TODO allocate only first m columns!]
  real_type_4d_view_type Jpp;     // bottom right  (within-particle)

  // LU factor views 
  real_type_4d_view_type D; // [batch, i_part, aero_sp_idx, aero_sp_idx]
  real_type_3d_view_type S;       // [batch, gas_sp_idx, gas_sp_idx]

  // W = D^{-1} * Jpg (block for each particle i); Programming note: Jpg = C / -gamma, so 
  // to solve D*W_true=C for W_true, need to carry -gamma*W
  real_type_4d_view_type W; // [batch, i_part, aero_sp_idx, gas_sp_idx] [TODO allocate only first m columns!]

  // state, RHS views with expanded particle contrib to gas tendency
  real_type_2d_view_type aero_f_2d; 
  real_type_3d_view_type gas_f_3d;
  real_type_2d_view_type y_perturb_2d;
  real_type_2d_view_type f_perturb_2d;
  real_type_2d_view_type aero_f_perturb_2d;
  real_type_3d_view_type gas_f_perturb_3d;
  real_type_2d_view_type f_gasonly_2d;

  // CVODE memory needed for accessing error weights
  void* cvode_mem;

  sunrealtype gamma;

};


struct SchurPrecon {

  // typedefs
  using device_type            = typename Tines::UseThisDevice<exec_space>::type;
  using problem_type = TChem::Impl::AerosolChemistry_Problem<real_type,device_type>;
  using real_type_1d_view_type = Tines::value_type_1d_view<real_type,device_type>;
  using real_type_2d_view_type = Tines::value_type_2d_view<real_type,device_type>;
  using real_type_3d_view_type = Tines::value_type_3d_view<real_type,device_type>;
  using real_type_4d_view_type = Tines::value_type_4d_view<real_type,device_type>;
  using ordinal_type_1d_view_type = Tines::value_type_1d_view<ordinal_type,device_type>;
  using AT = Kokkos::ArithTraits<real_type>;
  using range_type = Kokkos::pair<ordinal_type, ordinal_type>;

  // members
  // Batched state
  real_type_2d_view_type rhs;
  real_type_2d_view_type vals;
  real_type_2d_view_type aero_num_conc;
  real_type_2d_view_type const_tracers;
  real_type_1d_view_type temperature;
  real_type_1d_view_type pressure;
  KineticModelNCAR_ConstData<device_type> kmcd;
  AerosolModel_ConstData<device_type> amcd;
  ordinal_type_1d_view_type n_particles_track;

  // Views for Jacobian blocks 
  real_type_3d_view_type Jgg;     // upper left    (gas-gas)
  real_type_4d_view_type Jgp;     // upper right   (gas-particle)
  real_type_4d_view_type Jpg;     // bottom left   (particle-gas)
  real_type_4d_view_type Jpp;     // bottom right  (within-particle)

  // LU factor views (not including pivots!)
  real_type_4d_view_type D; // [batch, i_part, aero_sp_idx, aero_sp_idx]
  real_type_3d_view_type S;       // [batch, gas_sp_idx, gas_sp_idx]

  // W_ipart = D^{-1}_ipart * Jpg_ipart
  real_type_4d_view_type W; // [batch, i_part, aero_sp_idx, gas_sp_idx]

  const sunrealtype gamma;
  //ordinal_type n_systems;
  ordinal_type n_active_gas_species;
  ordinal_type n_total_species;
  ordinal_type n_equations;
  ordinal_type per_team_extent;
  ordinal_type n_condensible_gases;
  const ordinal_type level = 1;

  // construct/destructor
  SchurPrecon(const real_type_2d_view_type& rhs_in,
              const real_type_2d_view_type& vals_in,
              const real_type_2d_view_type& aero_num_conc_in,
              const real_type_2d_view_type& const_tracers_in,
              const real_type_1d_view_type& temperature_in,
              const real_type_1d_view_type& pressure_in,
              const KineticModelNCAR_ConstData<device_type>& kmcd_in,
              const AerosolModel_ConstData<device_type>& amcd_in,
              const ordinal_type_1d_view_type& n_particles_track_in,
              const real_type_3d_view_type& Jgg_in,
              const real_type_4d_view_type& Jgp_in,
              const real_type_4d_view_type& Jpg_in,
              const real_type_4d_view_type& Jpp_in,
              const real_type_4d_view_type& D_in,
              const real_type_3d_view_type& S_in,
              const real_type_4d_view_type& W_in,
              const sunrealtype gamma_in
            )
   : rhs(rhs_in),
     vals(vals_in),
     aero_num_conc(aero_num_conc_in),
     const_tracers(const_tracers_in),
     temperature(temperature_in),
     pressure(pressure_in),
     kmcd(kmcd_in),
     amcd(amcd_in),
     n_particles_track(n_particles_track_in),
     Jgg(Jgg_in),
     Jgp(Jgp_in),
     Jpg(Jpg_in),
     Jpp(Jpp_in),
     D(D_in),
     S(S_in),
     W(W_in),
     gamma(gamma_in) {

      //n_systems = rhs.extent(0); 
      n_active_gas_species = kmcd.nSpec - kmcd.nConstSpec;
      n_total_species = kmcd.nSpec + amcd.nParticles*amcd.nSpec;
      n_equations = problem_type::getNumberOfTimeODEs(kmcd, amcd);
      n_condensible_gases = amcd.nSimpol_tran; // NOTE this currently gives the # of condensible species for the SIMPOL mechanism, should replace with general condensible species count 

      // TODO this is just for the RHS.. need to figure out workspace sizing and where to put it
      per_team_extent
         = Impl::Aerosol_RHS<real_type, device_type>::getWorkSpaceSize(kmcd, amcd);

     }


  // TODO add workspace sizing
  static inline ordinal_type getWorkSpaceSize(){
    return 0;
  }


  template<typename MemberType>
  KOKKOS_INLINE_FUNCTION
  void team_setup(const MemberType &member, real_type_1d_view_type ewt_isys,
    // views for state and per-particle RHS contribution to gas tendencies
    real_type_1d_view_type aero_f, real_type_2d_view_type gas_f,
    real_type_1d_view_type y_perturb, real_type_1d_view_type f_perturb,
    real_type_1d_view_type f_gasonly, real_type_1d_view_type aero_f_perturb, 
    real_type_2d_view_type gas_f_perturb, sunbooleantype jok) const {

    const ordinal_type i_sys = member.league_rank();

    // system-level Jacobian blocks
    real_type_2d_view_type Jgg_isys = Kokkos::subview(Jgg, i_sys, Kokkos::ALL(), Kokkos::ALL());
    real_type_3d_view_type Jgp_isys = Kokkos::subview(Jgp, i_sys, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
    real_type_3d_view_type Jpg_isys = Kokkos::subview(Jpg, i_sys, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
    real_type_3d_view_type Jpp_isys = Kokkos::subview(Jpp, i_sys, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());

    // per-team subview of rhs, state, y, number conc., etc.
    const real_type_1d_view_type rhs_isys = Kokkos::subview(rhs, i_sys, Kokkos::ALL());
    const real_type_1d_view_type vals_isys = Kokkos::subview(vals, i_sys, Kokkos::ALL()); 
    const real_type_1d_view_type aero_num_conc_isys = Kokkos::subview(aero_num_conc, i_sys, Kokkos::ALL());
    const ordinal_type n_part_isys =
      (n_particles_track.span() == 0 || n_particles_track(i_sys) < 0) 
      ? amcd.nParticles : n_particles_track(i_sys);
    const real_type_1d_view_type constYs = Kokkos::subview(const_tracers, i_sys, Kokkos::ALL());
    const real_type temperature_isys = temperature(i_sys);
    const real_type pressure_isys = pressure(i_sys);

    // Construct Jacobian blocks if CVODE determines Jacobian stale
    if (!jok){

      // work for gas RHS
      TChem::Scratch<real_type_1d_view_type> work(member.team_scratch(level), per_team_extent);
      auto wptr = work.data();
      const ordinal_type problem_workspace_size = TChem::Impl::Aerosol_RHS<real_type, device_type>::getWorkSpaceSize(kmcd, amcd);
      auto pw = real_type_1d_view_type(wptr, problem_workspace_size);
      wptr +=problem_workspace_size;

      // fill unperturbed RHS for per particle contrib to gas tendencies
      TChem::Impl::Aerosol_RHS<real_type, device_type>::team_invoke_per_particle_gas_contrib(
        member, temperature_isys, pressure_isys, aero_num_conc_isys, vals_isys, aero_f, gas_f, 
        amcd, n_part_isys
      );
      member.team_barrier();

      // fill unperturbed gas-only RHS
      TChem::Impl::Aerosol_RHS<real_type, device_type>::team_invoke_gas(
        member, temperature_isys, pressure_isys, aero_num_conc_isys, vals_isys, constYs,
        f_gasonly, pw, kmcd, amcd, n_part_isys
      );
      member.team_barrier();

      // Fill the Jgg, Jpg columns for each active gas species
      for (int gas_sp_idx = 0; gas_sp_idx < n_active_gas_species; gas_sp_idx++){

          // compute sigma for current gas species
          //real_type weight = 1 / (reltol*AT::abs(vals_isys(gas_sp_idx)) + abstol_gas);
          real_type sigma = 1 / ewt_isys(gas_sp_idx);

          // create perturbed system state y_perturb with y_i + sigma_i
          Kokkos::single(Kokkos::PerTeam(member), [&]() {
            y_perturb(gas_sp_idx) = vals_isys(gas_sp_idx) + sigma;
          });
          member.team_barrier();

          if (gas_sp_idx < n_condensible_gases){
            // call regular rhs function with y_perturb
            TChem::Impl::Aerosol_RHS<real_type, device_type>::team_invoke( 
              member, temperature_isys, pressure_isys, aero_num_conc_isys, y_perturb, constYs, 
              f_perturb, pw, kmcd, amcd, n_part_isys);
            member.team_barrier();

            // fill Jpg_isys (per-particle) blocks for column gas_sp_idx
            Kokkos::parallel_for(Kokkos::TeamVectorMDRange(member, n_part_isys, amcd.nSpec), 
              [&] (const ordinal_type& part_idx, const ordinal_type& aero_sp_idx){
                ordinal_type i_row = n_active_gas_species + part_idx*amcd.nSpec + aero_sp_idx;
                Jpg_isys(part_idx, aero_sp_idx, gas_sp_idx) = (f_perturb(i_row) - rhs_isys(i_row)) / sigma;
            });
            member.team_barrier();

            // Difference quotient calculations
            // fill Jgg_isys for column gas_sp_idx
            Kokkos::parallel_for(Kokkos::TeamVectorRange(member, n_active_gas_species), 
              [&](const ordinal_type& i_gas){
                Jgg_isys(i_gas, gas_sp_idx) = (f_perturb(i_gas) - rhs_isys(i_gas)) / sigma;
            });
            member.team_barrier();
            
          } else {
            // call gas-only rhs function with y_perturb
            TChem::Impl::Aerosol_RHS<real_type, device_type>::team_invoke_gas( 
              member, temperature_isys, pressure_isys, aero_num_conc_isys, y_perturb, constYs, 
              f_perturb, pw, kmcd, amcd, n_part_isys);
            member.team_barrier();

            // Difference quotient calculations
            // fill Jgg_isys for column gas_sp_idx
            Kokkos::parallel_for(Kokkos::TeamVectorRange(member, n_active_gas_species), 
              [&](const ordinal_type& i_gas){
                Jgg_isys(i_gas, gas_sp_idx) = (f_perturb(i_gas) - f_gasonly(i_gas)) / sigma;
            });
            member.team_barrier();
            
          }

          // restore y_perturb
          Kokkos::single(Kokkos::PerTeam(member), [&]() {
            y_perturb(gas_sp_idx) = vals_isys(gas_sp_idx);
          });
          member.team_barrier();

      }

      // Fill in the Jpp, Jgp columns in parallel over number of aerosol species
      for (int aero_sp_idx = 0; aero_sp_idx < amcd.nSpec; aero_sp_idx++){

          Kokkos::parallel_for(Kokkos::TeamVectorRange(member, n_part_isys),
            [&](const ordinal_type& i_part) {
                // compute sigma for species aero_sp_idx in particle i_part (different val for each particle)
                ordinal_type part_spec_idx = n_active_gas_species + i_part*amcd.nSpec + aero_sp_idx;
                //real_type weight = 1 / (reltol*AT::abs(vals_isys(part_spec_idx)) + abstol_aero);
                real_type sigma = 1 / ewt_isys(part_spec_idx);
                
                // perturbed system state y_perturb where y_particleidx^j  + sigma_particleidx
                y_perturb(part_spec_idx) += sigma; 
          });
          member.team_barrier(); // need to fill perturbed state before calling RHS

          // pass perturbed state to modified RHS function, fmod(y_perturb)
          TChem::Impl::Aerosol_RHS<real_type, device_type>::team_invoke_per_particle_gas_contrib(
            member, temperature_isys, pressure_isys, aero_num_conc_isys, y_perturb, aero_f_perturb, gas_f_perturb, 
            amcd, n_part_isys
          );
          member.team_barrier();

          // difference quotient, store in Jpp_isys, Jgp_isys blocks
          Kokkos::parallel_for(Kokkos::TeamVectorRange(member, n_part_isys),
            [&](const ordinal_type& i_part) {
              ordinal_type part_spec_idx = n_active_gas_species + i_part*amcd.nSpec + aero_sp_idx;
              //real_type weight = 1 / (reltol*AT::abs(vals_isys(part_spec_idx)) + abstol_aero);
              real_type sigma = 1 / ewt_isys(part_spec_idx); 
              
              // Jgp_isys at particle i_part: loop for each condensible gas species in current aero_species column
              for (int gas_sp_idx = 0; gas_sp_idx < n_condensible_gases; gas_sp_idx++){
                Jgp_isys(i_part, gas_sp_idx, aero_sp_idx) = (gas_f_perturb(i_part, gas_sp_idx) - gas_f(i_part, gas_sp_idx)) / sigma;
              }
              
              // Jpp_isys at particle i_part: loop over jth aerosol species row for current aero_species column
              for (int aero_sp_j = 0; aero_sp_j < amcd.nSpec; aero_sp_j++){
                ordinal_type part_spec_j = n_active_gas_species + i_part*amcd.nSpec + aero_sp_j;
                Jpp_isys(i_part, aero_sp_j, aero_sp_idx) = (aero_f_perturb(part_spec_j) - aero_f(part_spec_j)) / sigma;
              }
                
          });
          member.team_barrier();

          // restore yhat (TeamVectorRange since need to restore strided per-particle species to non-perturbed values)
          Kokkos::parallel_for(Kokkos::TeamVectorRange(member, n_part_isys),
            [&](const ordinal_type& i_part) {
                ordinal_type part_spec_idx = n_active_gas_species + i_part*amcd.nSpec + aero_sp_idx; 
                y_perturb(part_spec_idx) = vals_isys(part_spec_idx); 
          });
          member.team_barrier();

      }

    } // if !jok

    // Subview D, W, S at the system level
    // D [batch, i_part, aero_sp_idx, aero_sp_idx]
    // W [batch, i_part, aero_sp_idx, gas_sp_idx]
    // S [batch, gas_sp_idx, gas_sp_idx]
    real_type_3d_view_type D_isys = Kokkos::subview(D, i_sys, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL()); 
    real_type_3d_view_type W_isys = Kokkos::subview(W, i_sys, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL()); 
    real_type_2d_view_type S_isys = Kokkos::subview(S, i_sys, Kokkos::ALL(), Kokkos::ALL()); 

    // construct D for each particle (store I-gamma*Jpp_ipart and LU factor in place)
    // Use linear index with aero species column index varying fastest for coalesced accesses
    // Tines uses layout right for views, so the rightmost index (column) has stride 1 in memory
    Kokkos::parallel_for(Kokkos::TeamVectorRange(member, n_part_isys*amcd.nSpec*amcd.nSpec),
      [&](const ordinal_type& k) {
        int i_part = k / (amcd.nSpec*amcd.nSpec); // rounds down to int
        int aero_sp_i = (k - i_part*amcd.nSpec*amcd.nSpec) / amcd.nSpec;
        int aero_sp_j = (k - i_part*amcd.nSpec*amcd.nSpec) % amcd.nSpec;
        if (aero_sp_i == aero_sp_j){
          D_isys(i_part, aero_sp_i, aero_sp_j) = 1 - gamma*Jpp_isys(i_part, aero_sp_i, aero_sp_j);
        } else{
          D_isys(i_part, aero_sp_i, aero_sp_j) = -1*gamma*Jpp_isys(i_part, aero_sp_i, aero_sp_j);
        }           
    });
    member.team_barrier();
    
    // Performance note: since the particle index is the slowest index in the D system subviews,
    // memory accesses for adjacent threads (which are handling adjacent D_i blocks) are 
    // *not* coalesced. In order to do so, the particle index would need to be the right-most
    // index (assuming LayoutRight views)
    Kokkos::parallel_for(Kokkos::TeamVectorRange(member, n_part_isys),
      [&](const ordinal_type& i_part) {

        // LU factor D_isys(i_part, :, :) in place
        real_type_2d_view_type D_ipart = Kokkos::subview(D_isys, i_part, 
          Kokkos::ALL(), Kokkos::ALL());  
        real_type tiny = 1e-8; // boost diagonal by small amount to move pivots off zero
        KokkosBatched::SerialLU<KokkosBatched::Algo::LU::Unblocked>::invoke(D_ipart, tiny);

        // place Jpg (linear system RHS) into W, pass to LU solve for all columns of W
        real_type_2d_view_type W_ipart = Kokkos::subview(W_isys, i_part, 
          Kokkos::ALL(), Kokkos::ALL());
        for (int a=0; a < amcd.nSpec; a++){
          for (int g=0; g < n_condensible_gases; g++){
            W_ipart(a, g) = Jpg_isys(i_part, a, g);
          }
        }
        // solve D_ipart*W_ipart = J_pg(i_part)
        KokkosBatched::SerialSolveLU
          <KokkosBatched::Trans::NoTranspose, KokkosBatched::Algo::Trsm::Unblocked>
            ::invoke(D_ipart, W_ipart);

    });
    member.team_barrier();

    // construct S = (I-gamma*Jgg) - gamma^2*sum(Jgp_ipart*W_ipart) (ipart=0, n_part_isys)
    // S <- I-gamma*Jgg
    Kokkos::parallel_for(Kokkos::TeamVectorRange(member, n_active_gas_species*n_active_gas_species),
      [&](const ordinal_type& k) {
        int gas_sp_i = k / n_active_gas_species;
        int gas_sp_j = k % n_active_gas_species;
        if (gas_sp_i == gas_sp_j){
          S_isys(gas_sp_i, gas_sp_j) = 1 - gamma*Jgg_isys(gas_sp_i, gas_sp_j);
        } else{
          S_isys(gas_sp_i, gas_sp_j) = -1*gamma*Jgg_isys(gas_sp_i, gas_sp_j);
        }           
    });
    member.team_barrier();

    // S += -gamma**2(Jgp_ipart*W_ipart) [ipart=0,..., n_part_isys]
    /*
    Kokkos::parallel_for(Kokkos::TeamVectorRange(member, n_part_isys*n_condensible_gases*n_condensible_gases),
      [&](const ordinal_type& k) {
        int i_part = k / (n_condensible_gases*n_condensible_gases); // rounds down to int
        int gas_sp_i = (k - i_part*n_condensible_gases*n_condensible_gases) / n_condensible_gases;
        int gas_sp_j = (k - i_part*n_condensible_gases*n_condensible_gases) % n_condensible_gases;    
        
        real_type sum = 0.0;
        for (int m = 0; m < amcd.nSpec; m++){
          sum += -1*gamma*gamma*Jgp_isys(i_part, gas_sp_i, m)*W_isys(i_part, m, gas_sp_j);
        }
        // atomic since each threads modifying same memory address
        Kokkos::atomic_add(&S_isys(gas_sp_i, gas_sp_j), sum);
    });
    member.team_barrier();
     */
    // Test heirarchical reduction version of S += -gamma**2(Jgp_ipart*W_ipart)
    Kokkos::parallel_for(Kokkos::TeamThreadRange(member, n_condensible_gases*n_condensible_gases),
      [&](const ordinal_type& k) {

        int gas_sp_i = k / n_condensible_gases;
        int gas_sp_j = k % n_condensible_gases;

        // sum matrix-matrix products at element (i,j) for Jgp*W from particle 0,...,n_part_isys
        real_type sum_JpgW_ij = 0.0;
        Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(member, n_part_isys*amcd.nSpec),
          [&](const ordinal_type& l, real_type& sum) {  
            
            int i_part = l / amcd.nSpec;
            int m = l % amcd.nSpec;

            sum += -1*gamma*gamma*Jgp_isys(i_part, gas_sp_i, m)*W_isys(i_part, m, gas_sp_j);
        }, sum_JpgW_ij);

        // use thread 0 in each team thread's vector lane to add to S
        Kokkos::single(Kokkos::PerThread(member), [&]() {
            S_isys(gas_sp_i, gas_sp_j) += sum_JpgW_ij;
        });        
    });
    member.team_barrier();

    // LU factor S in place
    Kokkos::single(Kokkos::PerTeam(member), [&]() {
      real_type tiny = 1e-8; // boost diagonal by small amount to move pivots off zero
      KokkosBatched::SerialLU<KokkosBatched::Algo::LU::Unblocked>::invoke(S_isys, tiny);
    });    
    
  }


  template<typename MemberType>
  KOKKOS_INLINE_FUNCTION
  void team_solve(const MemberType &member, real_type_1d_view_type x) const {

    // Solve the preconditioned system Px = b for x
    // Input linear system rhs b to x; return solution in place

    // system index
    const ordinal_type i_sys = member.league_rank();
    
    // number of active computational particles for the current system
    const ordinal_type n_part_isys =
      (n_particles_track.span() == 0 || n_particles_track(i_sys) < 0) 
      ? amcd.nParticles : n_particles_track(i_sys);

    // subview b -> b_gas, b_part (subviewing on x which contains the input, and since 
    // solves modify the view in-place, it also will carry the output)
    real_type_1d_view_type b_gas = Kokkos::subview(x, range_type(0, n_active_gas_species));
    real_type_1d_view_type b_part = Kokkos::subview(x, 
      range_type(n_active_gas_species, n_active_gas_species + n_part_isys*amcd.nSpec));
      
    // subview S, D, W, Jgp at isys
    real_type_2d_view_type S_isys = Kokkos::subview(S, i_sys, Kokkos::ALL(), Kokkos::ALL()); 
    real_type_3d_view_type D_isys = Kokkos::subview(D, i_sys, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
    real_type_3d_view_type W_isys = Kokkos::subview(W, i_sys, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL()); 
    real_type_3d_view_type Jgp_isys = Kokkos::subview(Jgp, i_sys, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());

    // Solve D_i*t_i = b_ipart (SolveLU places t_i in b_ipart, no explicit t required)
    Kokkos::parallel_for(Kokkos::TeamVectorRange(member, n_part_isys),
      [&](const ordinal_type& ipart) {
        real_type_2d_view_type D_ipart = Kokkos::subview(D_isys, ipart, 
          Kokkos::ALL(), Kokkos::ALL());
        real_type_1d_view_type b_ipart = Kokkos::subview(b_part, 
          range_type(ipart*amcd.nSpec, amcd.nSpec*(ipart+1)));
        
        KokkosBatched::SerialSolveLU
          <KokkosBatched::Trans::NoTranspose, KokkosBatched::Algo::Trsm::Unblocked>
            ::invoke(D_ipart, b_ipart);
    });
    member.team_barrier();

    // b_gas <- b_gas + gamma*sum(Jgp_isys*t_i)
    for (int m = 0; m < n_condensible_gases; m++){ // accumulate to first m entries of b_gas

      // sum mth index of matrix-vector product Jgp*t_i from particle 0,...,n_part_isys
      real_type sum_matvec = 0.0;
      Kokkos::parallel_reduce(Kokkos::TeamVectorRange(member, n_part_isys*amcd.nSpec),
        [&](const ordinal_type& l, real_type& accum) {  
          
          int i_part = l / amcd.nSpec;
          int a = l % amcd.nSpec;

          accum += gamma*Jgp_isys(i_part, m, a)*b_part(i_part*amcd.nSpec + a);
      }, sum_matvec);

      Kokkos::single(Kokkos::PerTeam(member), [&]() {
        b_gas(m) += sum_matvec;
      });
    }

    // Solve S*x_gas=b_gas
    Kokkos::single(Kokkos::PerTeam(member), [&]() {
      // performance note: would be nice to replace with a teamvector solve LU!!!
      KokkosBatched::SerialSolveLU
            <KokkosBatched::Trans::NoTranspose, KokkosBatched::Algo::Trsm::Unblocked>
              ::invoke(S_isys, b_gas); // result stored in b_gas
    });
    member.team_barrier();

    // solve x_ipart = t_ipart + gamma*W_ipart*x_gas[:n_condensible_gases]
    // note that W was solved as D*W=Jpg (note the RHS is the jacobian block instead of 
    // C, -gamma*Jpg, so we need to multiply resulting W by -gamma, hence the +gamma*W...)
    Kokkos::parallel_for(Kokkos::TeamVectorRange(member, n_part_isys*amcd.nSpec),
    [&](const ordinal_type& k) {
      int i_part = k / amcd.nSpec;
      int a = k % amcd.nSpec;

      for (int m = 0; m < n_condensible_gases; m++){
        b_part(i_part*amcd.nSpec+a) +=  gamma*W_isys(i_part, a, m)*b_gas(m);
      }
    });
    
  } // team_solve
};


struct AerosolChemistry_CVODE_K
{
  using device_type      = typename Tines::UseThisDevice<exec_space>::type;
  using real_type_1d_view_type = Tines::value_type_1d_view<real_type,device_type>;
  using real_type_2d_view_type = Tines::value_type_2d_view<real_type,device_type>;
  using policy_type = typename UseThisTeamPolicy<exec_space>::type;
  using problem_type = TChem::Impl::AerosolChemistry_Problem<real_type,device_type>;
  using MatType   = sundials::kokkos::DenseMatrix<exec_space>;

  // User-supplied functions called by CVODE
  static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);

  static int Jac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
                 void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
  
  static int psetup(sunrealtype t, N_Vector y, N_Vector fy, sunbooleantype jok,
                    sunbooleantype* jcurPtr, sunrealtype gamma, void* user_data);

  static int psolve(sunrealtype t, N_Vector y, N_Vector fy, N_Vector r, 
                    N_Vector z, sunrealtype gamma, sunrealtype delta, int lr, 
                    void* user_data);

};
} // namespace TChem

#endif
