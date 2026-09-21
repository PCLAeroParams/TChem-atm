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
#ifndef __TCHEM_AEROSOL_CHEMISTRY_IMPLICIT_EULER_HPP__
#define __TCHEM_AEROSOL_CHEMISTRY_IMPLICIT_EULER_HPP__

#include "TChem_Util.hpp"
#include "TChem_KineticModelData.hpp"
#include "TChem_Impl_AerosolChemistry_Problem.hpp"
#include "TChem_Impl_TimeIntegratorImplicitEuler.hpp"
#include "TChem_AerosolChemistry.hpp"  // for TeamConfOutput

namespace TChem {

struct AerosolChemistry_ImplicitEuler
{
  using host_device_type = typename Tines::UseThisDevice<host_exec_space>::type;
  using device_type      = typename Tines::UseThisDevice<exec_space>::type;

  using real_type_0d_view_type = Tines::value_type_0d_view<real_type,device_type>;
  using real_type_1d_view_type = Tines::value_type_1d_view<real_type,device_type>;
  using real_type_2d_view_type = Tines::value_type_2d_view<real_type,device_type>;

  using real_type_0d_view_host_type = Tines::value_type_0d_view<real_type,host_device_type>;
  using real_type_1d_view_host_type = Tines::value_type_1d_view<real_type,host_device_type>;
  using real_type_2d_view_host_type = Tines::value_type_2d_view<real_type,host_device_type>;

  template<typename DeviceType>
  static inline ordinal_type getWorkSpaceSize(
    const KineticModelNCAR_ConstData<DeviceType>& kmcd,
    const AerosolModel_ConstData<DeviceType>& amcd)
  {
    using device_type = DeviceType;
    using problem_type = Impl::AerosolChemistry_Problem<real_type, device_type>;
    using time_integrator_type = Impl::TimeIntegratorImplicitEuler<real_type, device_type>;

    const ordinal_type m = problem_type::getNumberOfEquations(kmcd, amcd) + 1;

    ordinal_type work_size_problem(0);
#if defined(TCHEM_ATM_ENABLE_SACADO_JACOBIAN_AEROSOL_CHEMISTRY)
    if (m < 32) {
      using value_type = Sacado::Fad::SLFad<real_type,32>;
      using problem_value_type = Impl::AerosolChemistry_Problem<value_type, device_type>;
      work_size_problem = problem_value_type::getWorkSpaceSize(kmcd, amcd);
    } else if (m < 64) {
      using value_type = Sacado::Fad::SLFad<real_type,64>;
      using problem_value_type = Impl::AerosolChemistry_Problem<value_type, device_type>;
      work_size_problem = problem_value_type::getWorkSpaceSize(kmcd, amcd);
    } else if (m < 128) {
      using value_type = Sacado::Fad::SLFad<real_type,128>;
      using problem_value_type = Impl::AerosolChemistry_Problem<value_type, device_type>;
      work_size_problem = problem_value_type::getWorkSpaceSize(kmcd, amcd);
    } else if (m < 256) {
      using value_type = Sacado::Fad::SLFad<real_type,256>;
      using problem_value_type = Impl::AerosolChemistry_Problem<value_type, device_type>;
      work_size_problem = problem_value_type::getWorkSpaceSize(kmcd, amcd);
    } else if (m < 512) {
      using value_type = Sacado::Fad::SLFad<real_type,512>;
      using problem_value_type = Impl::AerosolChemistry_Problem<value_type, device_type>;
      work_size_problem = problem_value_type::getWorkSpaceSize(kmcd, amcd);
    } else if (m < 1024) {
      using value_type = Sacado::Fad::SLFad<real_type,1024>;
      using problem_value_type = Impl::AerosolChemistry_Problem<value_type, device_type>;
      work_size_problem = problem_value_type::getWorkSpaceSize(kmcd, amcd);
    } else {
      TCHEM_CHECK_ERROR(0,
                        "Error: Number of equations is bigger than size of sacado fad type");
    }
#else
    {
      work_size_problem = problem_type::getWorkSpaceSize(kmcd, amcd);
    }
#endif

    ordinal_type wlen(0);
    time_integrator_type::workspace(m - 1, wlen);

    return (m - 1) + wlen + work_size_problem;
  }

  static void
  runHostBatch( /// thread block size
           typename UseThisTeamPolicy<host_exec_space>::type& policy,
           /// input
           const real_type_1d_view_host& tol_newton,
           const real_type_2d_view_host& tol_time,
           const real_type_2d_view_host& fac,
           const time_advance_type_1d_view_host& tadv,
           const real_type_2d_view_host& state,
           const real_type_2d_view_host& number_conc,
           /// output
           const real_type_1d_view_host& t_out,
           const real_type_1d_view_host& dt_out,
           const real_type_2d_view_host& state_out,
           TeamConfOutput& team_conf_output,
           /// const data from kinetic model
           const KineticModelNCAR_ConstData<interf_host_device_type>& kmcd,
           const AerosolModel_ConstData<interf_host_device_type>& amcd
           );

  static void
  runDeviceBatch( /// thread block size
           typename UseThisTeamPolicy<exec_space>::type& policy,
           /// input
           const real_type_1d_view& tol_newton,
           const real_type_2d_view& tol_time,
           const real_type_2d_view& fac,
           const time_advance_type_1d_view& tadv,
           const real_type_2d_view& state,
           const real_type_2d_view& number_conc,
           /// output
           const real_type_1d_view& t_out,
           const real_type_1d_view& dt_out,
           const real_type_2d_view& state_out,
           TeamConfOutput& team_conf_output,
           /// const data from kinetic model
           const KineticModelNCAR_ConstData<device_type>& kmcd,
           const AerosolModel_ConstData<device_type>& amcd
           );

};

} // namespace TChem

#endif
