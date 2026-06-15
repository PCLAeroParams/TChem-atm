#ifndef __TCHEM_IMPL_MOSAIC_HPP__
#define __TCHEM_IMPL_MOSAIC_HPP__

/*
Code based on the Model for Simulating Aerosol
Interactions and Chemistry (MOSAIC)

Zaveri, R. A., Easter, R. C., Fast, J. D., & Peters, L. K. (2008).
Model for simulating aerosol interactions and chemistry (MOSAIC).
Journal of Geophysical Research: Atmospheres.
*/

#include "TChem.hpp"
#include "TChem_Util.hpp"
#include "TChem_Impl_MOSAIC_ModelData.hpp"

namespace TChem {
namespace Impl {

using Kokkos::min;
using Kokkos::max;

template<typename ValueType, typename DeviceType>
struct MOSAIC{

  using value_type = ValueType;
  using device_type = DeviceType;
  using scalar_type = typename ats<value_type>::scalar_type;

  using real_type = scalar_type;

  using value_type_1d_view_type = Tines::value_type_1d_view<value_type,device_type>;
  using value_type_2d_view_type = Tines::value_type_2d_view<value_type,device_type>;
  using real_type_1d_view_type = Tines::value_type_1d_view<real_type,device_type>;
  using real_type_2d_view_type = Tines::value_type_2d_view<real_type,device_type>;

  KOKKOS_INLINE_FUNCTION static ordinal_type getWorkSpaceSize() {
    ordinal_type workspace_size=0;
    return workspace_size;
  }

  #include "TChem_Impl_MOSAIC_Util.hpp"

  #include "TChem_Impl_MOSAIC_MESA.hpp"

  #include "TChem_Impl_MOSAIC_ASTEM.hpp"

};

} // namespace Impl
} // namespace TChem

#endif
