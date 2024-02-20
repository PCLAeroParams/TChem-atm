#ifndef __TCHEM_LINV3_STRAT_CHEM_SOLVER_HPP__
#define __TCHEM_LINV3_STRAT_CHEM_SOLVER_HPP__

#include "TChem_KineticModelData.hpp"
#include "TChem_Util.hpp"
#include "TChem_Impl_linv3_strat_chem_solve.hpp"

namespace TChem {
#if 1
struct Linv3StratosphereSolver
{
  using host_device_type = typename Tines::UseThisDevice<host_exec_space>::type;
  using device_type      = typename Tines::UseThisDevice<exec_space>::type;

  using real_type_1d_view_type = Tines::value_type_1d_view<real_type,device_type>;
  using real_type_2d_view_type = Tines::value_type_2d_view<real_type,device_type>;
  using ordinal_type_1d_view_type = Tines::value_type_1d_view<ordinal_type,device_type>;

  using real_type_1d_view_host_type = Tines::value_type_1d_view<real_type,host_device_type>;
  using real_type_2d_view_host_type = Tines::value_type_2d_view<real_type,host_device_type>;
  using ordinal_type_1d_view_host_type = Tines::value_type_1d_view<ordinal_type,host_device_type>;

  static void runHostBatch( /// input
    const real_type_1d_view_host_type& temperature,
  const real_type_1d_view_host_type& pressure,

  const real_type_2d_view_host_type& volume_mixing_ratio,
  const real_type dt, const real_type rlats, const real_type psc_T, const real_type  sza,
  const real_type chlorine_loading, 
  const real_type_1d_view_host_type& o3col,
  const ordinal_type_1d_view_host_type& tropFlag,
  const real_type_1d_view_host_type& water_vapor_volume_mixing_ratio,
  const linoz_input_parameters_1d_view_host& linoz_inputs, 
  const linoz_vmr_idx_type& linoz_vmr_idx,
  const real_type_2d_view_host_type& volume_mixing_ratio_out);

  static void runDeviceBatch( /// input
    typename UseThisTeamPolicy<exec_space>::type& policy,
   const real_type_1d_view_type& temperature,
  const real_type_1d_view_type& pressure,
  const real_type_2d_view_type& volume_mixing_ratio,
  const real_type dt, const real_type rlats, const real_type psc_T, const real_type  sza,
  const real_type chlorine_loading, 
  const real_type_1d_view_type& o3col,
  const ordinal_type_1d_view_type& tropFlag,
  const real_type_1d_view_type& water_vapor_volume_mixing_ratio,
  const linoz_input_parameters_1d_view& linoz_inputs, 
  const linoz_vmr_idx_type& linoz_vmr_idx,
  const real_type_2d_view_type& volume_mixing_ratio_out);

  static void runDeviceBatch( /// input
  const real_type_1d_view_type& temperature,
  const real_type_1d_view_type& pressure,
  const real_type_2d_view_type& volume_mixing_ratio,
  const real_type dt, const real_type rlats, const real_type psc_T, const real_type  sza,
  const real_type chlorine_loading, 
  const real_type_1d_view_type& o3col,
  const ordinal_type_1d_view_type& tropFlag,
  const real_type_1d_view_type& water_vapor_volume_mixing_ratio,
  const linoz_input_parameters_1d_view& linoz_inputs, 
  const linoz_vmr_idx_type& linoz_vmr_idx,
  const real_type_2d_view_type& volume_mixing_ratio_out);


};
#endif
} // namespace TChem

#endif
