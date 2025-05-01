/* =====================================================================================
TChem-atm version 1.0
Copyright (2024) NTESS
https://github.com/sandialabs/TChem-atm

Copyright 2024 National Technology & Engineering Solutions of Sandia, LLC
(NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
Government retains certain rights in this software.

This file is part of TChem-atm. TChem-atm is open source software: you can redistribute
it and/or modify it under the terms of BSD 2-Clause License
(https://opensource.org/licenses/BSD-2-Clause). A copy of the licese is also
provided under the main directory

Questions? Contact Oscar Diaz-Ibarra at <odiazib@sandia.gov>, or
           Mike Schmidt at <mjschm@sandia.gov>, or
           Cosmin Safta at <csafta@sandia.gov>

Sandia National Laboratories, New Mexico/Livermore, NM/CA, USA
=====================================================================================
*/
#ifndef __TCHEM_AEROSOL_CHEMISTRY_CVODE_RHS_JACOBIAN_HPP__
#define __TCHEM_AEROSOL_CHEMISTRY_CVODE_RHS_JACOBIAN_HPP__

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

namespace TChem {

struct UserData
{
  using host_device_type = typename Tines::UseThisDevice<host_exec_space>::type;
  using device_type      = typename Tines::UseThisDevice<exec_space>::type;

  using real_type_1d_view_type = Tines::value_type_1d_view<real_type,device_type>;
  using real_type_2d_view_type = Tines::value_type_2d_view<real_type,device_type>;
  using real_type_3d_view_type = Tines::value_type_3d_view<real_type,device_type>;
  using policy_type = typename UseThisTeamPolicy<exec_space>::type;

  int nbatches  = 100; // number of chemical networks
  int batchSize = 3;   // size of each network
  policy_type policy;
  real_type_2d_view_type num_concentration;
  real_type_1d_view_type temperature;
  real_type_1d_view_type pressure;
  real_type_2d_view_type const_tracers;
  real_type_2d_view_type fac;
#if defined(TCHEM_ATM_ENABLE_GPU)
  real_type_3d_view_type JacRL;
#endif

  TChem::KineticModelNCAR_ConstData<device_type> kmcd;
  TChem::AerosolModel_ConstData<device_type> amcd;
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

};
} // namespace TChem

#endif
