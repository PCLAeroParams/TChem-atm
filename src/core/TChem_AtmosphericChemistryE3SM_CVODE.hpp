#ifndef __TCHEM_ATMOSPHERIC_CHEMISTRY_E3SM_CVODE_HPP__
#define __TCHEM_ATMOSPHERIC_CHEMISTRY_E3SM_CVODE_HPP__

#include "TChem_KineticModelData.hpp"
#include "TChem_Util.hpp"
#include "TChem_Impl_AtmosphericChemistryE3SM_CVODE.hpp"
#include "TChem_Impl_AtmosphericChemistryE3SM.hpp"

namespace TChem {

struct AtmosphericChemistryE3SM_CVODE
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
    const KineticModelNCAR_ConstData<DeviceType>& kmcd)
  {
    using device_type = DeviceType;
    using problem_type = Impl::AtmosphericChemistryE3SM_Problem<real_type, device_type>;
    const ordinal_type m = problem_type::getNumberOfEquations(kmcd);

    ordinal_type work_size(0);
    work_size = problem_type::getWorkSpaceSize(kmcd);

    return work_size + m*m;

  }

	   static void 
  runHostBatch( /// thread block size
           typename UseThisTeamPolicy<host_exec_space>::type& policy,
           /// input
           const real_type_2d_view_host& tol,
           const real_type_2d_view_host& fac,
           const time_advance_type_1d_view_host& tadv,
           const real_type_2d_view_host& state,
           const real_type_2d_view_host& photo_rates,
           const real_type_2d_view_host& external_sources,
           /// output
           const real_type_1d_view_host& t_out,
           const real_type_1d_view_host& dt_out,
           const real_type_2d_view_host& state_out,
           /// const data from kinetic model
           const Tines::value_type_1d_view<KineticModelNCAR_ConstData<interf_host_device_type>, interf_host_device_type>& kmcds,
           const Tines::value_type_1d_view<Tines::TimeIntegratorCVODE<real_type, host_device_type>, host_device_type>& cvodes);

};

 
} // namespace TChem

#endif