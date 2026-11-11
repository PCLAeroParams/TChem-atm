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
#include "TChem_AerosolChemistry_CVODE_K_UserFunctions.hpp"

namespace TChem
{

// Right hand side function dy/dt = f(t,y)
int AerosolChemistry_CVODE_K::f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  auto udata = static_cast<UserData*>(user_data);

  const auto n_systems  = udata->n_systems;
  const auto policy = udata->policy;
  const auto kmcd = udata->kmcd;
  const auto amcd = udata->amcd;
  const auto system_size = udata->system_size;
  const auto aero_num_conc = udata->aero_num_conc;
  const auto temperature= udata->temperature;
  const auto pressure= udata->pressure;
  const auto const_tracers= udata->const_tracers;

  real_type_2d_view_type vals(N_VGetDeviceArrayPointer(y), n_systems, system_size);
  real_type_2d_view_type rhs(N_VGetDeviceArrayPointer(ydot), n_systems, system_size);

  const ordinal_type level = 1;
  // const ordinal_type number_of_equations = system_size;
  const ordinal_type per_team_extent
         = TChem::Impl::Aerosol_RHS<real_type, device_type>::getWorkSpaceSize(kmcd, amcd);
  const std::string profile_name = "TChem::AerosolChemistry::RHS_evaluation";

  const auto n_particles_track = udata->n_particles_track;
  TChemAerosolChemistryRHS rhs_tchem(rhs, vals, aero_num_conc,
     const_tracers, temperature, pressure, kmcd, amcd, n_particles_track);

  Kokkos::Profiling::pushRegion(profile_name);
  Kokkos::parallel_for(profile_name, policy, rhs_tchem);
  Kokkos::Profiling::popRegion();

  udata->team_size_recommended_rhs = policy.team_size_recommended(
    rhs_tchem, Kokkos::ParallelForTag());
  udata->team_size_max_rhs = policy.team_size_max(
    rhs_tchem, Kokkos::ParallelForTag());

 return 0;
}


// Jacobian of f(t,y)
int AerosolChemistry_CVODE_K::Jac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J, void* user_data
, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  auto udata = static_cast<UserData*>(user_data);

  const auto n_systems  = udata->n_systems;
  const auto policy = udata->policy;
  const auto kmcd = udata->kmcd;
  const auto amcd = udata->amcd;
  const auto aero_num_conc = udata->aero_num_conc;
  const auto temperature= udata->temperature;
  const auto pressure= udata->pressure;
  const auto const_tracers= udata->const_tracers;
  const auto system_size = udata->system_size;
  const auto fac= udata->fac;
#if defined(TCHEM_ATM_ENABLE_GPU)
  const auto JacRL= udata->JacRL;
#endif
  auto J_data = sundials::kokkos::GetDenseMat<MatType>(J)->View();

  real_type_2d_view_type y2d(N_VGetDeviceArrayPointer(y), n_systems, system_size);
  real_type_2d_view_type vals(N_VGetDeviceArrayPointer(y), n_systems, system_size);

  const ordinal_type level = 1;
  const ordinal_type number_of_equations = problem_type::getNumberOfTimeODEs(kmcd, amcd);
  const ordinal_type per_team_extent = problem_type::getWorkSpaceSize(kmcd,amcd)
    + number_of_equations;
  const std::string profile_name = "TChem::AerosolChemistry::Jacobian_evaluation";
  auto jac_lambda =        KOKKOS_LAMBDA(const typename policy_type::member_type& member) {

    const ordinal_type i = member.league_rank();

    const ordinal_type m = problem_type::getNumberOfTimeODEs(kmcd,amcd);
    const real_type_1d_view_type vals_at_i =
    Kokkos::subview(vals, i, Kokkos::ALL());

    const real_type_1d_view_type fac_at_i =
    Kokkos::subview(fac, i, Kokkos::ALL());

    const real_type_2d_view_type jacobian_at_i =
#if defined(TCHEM_ATM_ENABLE_GPU)
    Kokkos::subview(JacRL, i, Kokkos::ALL(), Kokkos::ALL());
#else
    Kokkos::subview(J_data, i, Kokkos::ALL(), Kokkos::ALL());
#endif

    const real_type_1d_view_type number_conc_at_i =
    Kokkos::subview(aero_num_conc, i, Kokkos::ALL());

    TChem::Scratch<real_type_1d_view_type> work(member.team_scratch(level),
                                   per_team_extent);

    auto wptr = work.data();

    const real_type_1d_view_type constYs  = Kokkos::subview(const_tracers, i, Kokkos::ALL());

    const ordinal_type problem_workspace_size = problem_type::getWorkSpaceSize(kmcd,amcd);
    auto pw = real_type_1d_view_type(wptr, problem_workspace_size);
    wptr +=problem_workspace_size;
    problem_type problem;
    problem._kmcd = kmcd;
    problem._amcd = amcd;
    /// initialize problem
    problem._fac = fac_at_i;
    problem._work = pw;
    problem._temperature= temperature(i);
    problem._pressure =pressure(i);
    problem._const_concentration= constYs;
    problem._number_conc =number_conc_at_i;
    problem.computeNumericalJacobian(member,vals_at_i,jacobian_at_i);
#if defined(TCHEM_ATM_ENABLE_GPU)
    //NOTE: Sundials uses a left layout on the device, while we are using a right layout.
    // This incompatible layouts is producing a runtime error.
    for (int k=0;k<system_size;++k)
      for (int j=0;j<system_size;++j)
       J_data(i, k, j) = jacobian_at_i(k,j);
#endif
    };
  Kokkos::Profiling::pushRegion(profile_name);
  Kokkos::parallel_for(profile_name, policy, jac_lambda);
  Kokkos::Profiling::popRegion();

  udata->team_size_recommended_jac = policy.team_size_recommended(
    jac_lambda, Kokkos::ParallelForTag());
  udata->team_size_max_jac = policy.team_size_max(
    jac_lambda, Kokkos::ParallelForTag());

  return 0;

}


int AerosolChemistry_CVODE_K::psetup(sunrealtype t, N_Vector y, N_Vector fy, 
  sunbooleantype jok, sunbooleantype* jcurPtr, sunrealtype gamma, void* user_data)
{

  auto udata = static_cast<UserData*>(user_data);

  const auto n_systems  = udata->n_systems;
  const auto policy = udata->policy;
  const auto kmcd = udata->kmcd;
  const auto amcd = udata->amcd;
  const auto system_size = udata->system_size;
  const auto aero_num_conc = udata->aero_num_conc;
  const auto temperature= udata->temperature;
  const auto pressure= udata->pressure;
  const auto const_tracers= udata->const_tracers;
  const auto n_particles_track = udata->n_particles_track;

  // preconditioner data structures
  const auto Jgg = udata->Jgg;
  const auto Jgp = udata->Jgp;
  const auto Jpg = udata->Jpg;
  const auto Jpp = udata->Jpp;
  const auto D = udata->D;
  const auto S = udata->S;
  const auto W = udata->W;
  const auto aero_f_2d = udata->aero_f_2d;
  const auto gas_f_3d = udata->gas_f_3d;
  const auto y_perturb_2d = udata->y_perturb_2d; 
  const auto f_perturb_2d = udata->f_perturb_2d;
  const auto aero_f_perturb_2d = udata->aero_f_perturb_2d;
  const auto gas_f_perturb_3d = udata->gas_f_perturb_3d;
  const auto f_gasonly_2d = udata->f_gasonly_2d;
  //const auto gamma = udata->gamma; // TODO remove from userdata?
  const auto cvode_mem = udata->cvode_mem;

  real_type_2d_view_type vals(N_VGetDeviceArrayPointer(y), n_systems, system_size);
  real_type_2d_view_type rhs(N_VGetDeviceArrayPointer(fy), n_systems, system_size);

  //const ordinal_type level = 1;
  //const ordinal_type per_team_extent = TChem::Impl::Aerosol_RHS<real_type, device_type>::getWorkSpaceSize(kmcd, amcd);

  N_Vector ewt_nvec = N_VClone(y);
  CVodeGetErrWeights(cvode_mem, ewt_nvec);
  real_type_2d_view_type ewt(N_VGetDeviceArrayPointer(ewt_nvec), n_systems, system_size);

  SchurPrecon preconditioner(rhs, vals, aero_num_conc, const_tracers, temperature, 
    pressure, kmcd, amcd, n_particles_track, Jgg, Jgp, Jpg, Jpp, D, S, W, gamma);

  // initialize y_perturb_2d with vals if stale jacobian
  if (!jok){
    Kokkos::parallel_for(policy, 
      KOKKOS_LAMBDA(const typename policy_type::member_type& member){
        const ordinal_type i_sys = member.league_rank();
        Kokkos::parallel_for(Kokkos::TeamVectorRange(member, system_size), 
          [&] (const ordinal_type& i_spec){
            y_perturb_2d(i_sys, i_spec) = vals(i_sys, i_spec);
        });
    });
    *jcurPtr = SUNTRUE;
  } else {
    *jcurPtr = SUNFALSE;
  }
  
  const std::string profile_name = "TChem::AerosolChemistry_CVODE_K::psetup";
  
  Kokkos::Profiling::pushRegion(profile_name);
  Kokkos::parallel_for(profile_name, policy, 
    KOKKOS_LAMBDA(const typename policy_type::member_type& member){

      const ordinal_type i_sys = member.league_rank();

      // team-level subviews
      real_type_1d_view_type ewt_isys = Kokkos::subview(ewt, i_sys, Kokkos::ALL());
      real_type_1d_view_type y_perturb = Kokkos::subview(y_perturb_2d, i_sys, Kokkos::ALL());
      real_type_1d_view_type f_perturb = Kokkos::subview(f_perturb_2d, i_sys, Kokkos::ALL());
      real_type_1d_view_type f_gasonly = Kokkos::subview(f_gasonly_2d, i_sys, Kokkos::ALL());
      real_type_1d_view_type aero_f = Kokkos::subview(aero_f_2d, i_sys, Kokkos::ALL());
      real_type_2d_view_type gas_f = Kokkos::subview(gas_f_3d, i_sys, Kokkos::ALL(), Kokkos::ALL());
      real_type_1d_view_type aero_f_perturb = Kokkos::subview(aero_f_perturb_2d, i_sys, Kokkos::ALL());
      real_type_2d_view_type gas_f_perturb = Kokkos::subview(gas_f_perturb_3d, i_sys, Kokkos::ALL(), Kokkos::ALL());

      preconditioner.team_setup(member, ewt_isys, aero_f, gas_f, y_perturb, 
          f_perturb, f_gasonly, aero_f_perturb, gas_f_perturb, jok);

  });
  Kokkos::fence();
  Kokkos::Profiling::popRegion();
  

  N_VDestroy(ewt_nvec);
  return 0;

}


int AerosolChemistry_CVODE_K::psolve(sunrealtype t, N_Vector y, N_Vector fy, 
  N_Vector r, N_Vector z, sunrealtype gamma, sunrealtype delta, int lr, void* user_data)
{

  auto udata = static_cast<UserData*>(user_data);

  return 0;

}

}// namespace TChem
