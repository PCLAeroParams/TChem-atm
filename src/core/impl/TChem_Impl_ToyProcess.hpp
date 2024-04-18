#ifndef __TCHEM_IMPL_TOY_PROCESS_HPP__
#define __TCHEM_IMPL_TOY_PROCESS_HPP__

#include "TChem_Util.hpp"
#include "TChem_Impl_SingleParticleUtils.hpp"

namespace TChem {
namespace Impl {
  template<typename ValueType, typename DeviceType>
struct ToyProcess
{
  using value_type = ValueType;
  using device_type = DeviceType;
  using scalar_type = typename ats<value_type>::scalar_type;

  using real_type = scalar_type;
  /// sacado is value type
  using value_type_1d_view_type = Tines::value_type_1d_view<value_type,device_type>;
  using real_type_1d_view_type = Tines::value_type_1d_view<real_type,device_type>;

  using aerosol_model_data_type= AerosolModel_ConstData<device_type>;

  KOKKOS_INLINE_FUNCTION static ordinal_type getWorkSpaceSize()
  {
    ordinal_type workspace_size=0;
    return workspace_size;
  }
    template<typename MemberType>
  KOKKOS_INLINE_FUNCTION static void team_invoke(
    const MemberType& member,
    const real_type& t,
    const real_type& p,
    const ordinal_type i_part,
    const ordinal_type i_toy,
    const real_type_1d_view_type& number_conc,
    const value_type_1d_view_type& state,
    const value_type_1d_view_type& omega,
    const aerosol_model_data_type& amcd

    )
    {
      // 1. Let's compute the reaction/mass transfer constant
      // use real type if all variables in the operation are real_type
      const auto& toy_params = amcd.toy_params(i_toy);
      const real_type k = toy_params.A * ats<real_type>::exp(toy_params.B/t)  ;
      ordinal_type GAS_SPEC_=toy_params.gas_sp_index;
      ordinal_type AERO_SPEC_i_phase=toy_params.aerosol_sp_index;

      value_type cond_rate = k*state(GAS_SPEC_);
      value_type evap_rate = k*state(AERO_SPEC_i_phase+i_part*amcd.nSpec);

      // // per-particle mass concentration rates
      value_type diff = - evap_rate + cond_rate;

      // update rhs
      Kokkos::atomic_add(&omega(GAS_SPEC_), -number_conc(i_part) * diff);
      Kokkos::atomic_add(&omega(AERO_SPEC_i_phase+i_part*amcd.nSpec), diff);
    }

};

} // namespace Impl
} // namespace TChem

#endif
