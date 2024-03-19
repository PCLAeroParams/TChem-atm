#ifndef __TCHEM_IMPL_AEROSOL_RHS_HPP__
#define __TCHEM_IMPL_AEROSOL_RHS_HPP__

#include "TChem_Impl_SIMPOL_phase_transfer.hpp"

namespace TChem {
namespace Impl {
  template<typename ValueType, typename DeviceType>
struct Aerosol_RHS
{
  using value_type = ValueType;
  using device_type = DeviceType;
  using scalar_type = typename ats<value_type>::scalar_type;

  using real_type = scalar_type;
  /// sacado is value type
  using value_type_1d_view_type = Tines::value_type_1d_view<value_type,device_type>;
  using real_type_1d_view_type = Tines::value_type_1d_view<real_type,device_type>;
  using ordinal_type_1d_view_type = Tines::value_type_1d_view<ordinal_type,device_type>;
  using aerosol_model_data_type= AerosolModel_ConstData<device_type>;

  using const_real_type_1d_view_type = Tines::value_type_1d_view<const real_type,device_type>;

  KOKKOS_INLINE_FUNCTION static ordinal_type getWorkSpaceSize()
  {
    ordinal_type workspace_size=0;
    return workspace_size;
  }
    // update RHS
  template<typename MemberType>
  KOKKOS_INLINE_FUNCTION static
  void team_invoke(const MemberType& member,
    const real_type& t,
    const real_type& p,
    const real_type_1d_view_type& number_conc,
    const value_type_1d_view_type& state,
    const value_type_1d_view_type& omega,
    const aerosol_model_data_type& amcd
    )
  {

    using SIMPOL_single_particle_type = TChem::Impl::SIMPOL_single_particle<value_type, device_type >;

    for (int i_part = 0; i_part < amcd.nParticles; i_part++)
    {
    // FIXME: number of simpol is hard-coded to one.
    ordinal_type i_simpol=0;
    SIMPOL_single_particle_type
    ::invoke_team(member, i_part,i_simpol,
                  t, p, number_conc,
                  state, omega,
                  amcd);
    }

  }

};
} // namespace Impl
} // namespace TChem

#endif
