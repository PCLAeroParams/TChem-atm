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
#ifndef __TCHEM_AEROSOL_CHEMISTRY_HPP__
#define __TCHEM_AEROSOL_CHEMISTRY_HPP__


#include "TChem_Util.hpp"
#include "TChem_KineticModelData.hpp"
#include "TChem_Impl_AerosolChemistry.hpp"

namespace TChem {

struct TeamConfOutput
{
   ordinal_type team_size_recommended{-1};
   ordinal_type team_size_max{-1};
};

struct AerosolChemistry
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
    const ordinal_type m = problem_type::getNumberOfEquations(kmcd,amcd);
    ordinal_type work_size(0);
    work_size = Impl::AerosolChemistry<real_type, device_type>::getWorkSpaceSize(kmcd,amcd);
    return work_size + m;
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
