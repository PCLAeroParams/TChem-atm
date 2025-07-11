#ifndef __TCHEM_TEST_MOSAIC_HPP__
#define __TCHEM_TEST_MOSAIC_HPP__

#include "TChem.hpp"
#include "TChem_Impl_MOSAIC.hpp"

namespace TChem {
namespace Test {
void static MOSAIC_Test() {

    using device_type = typename Tines::UseThisDevice<TChem::exec_space>::type;
    const auto host_exec_space = TChem::host_exec_space();
    const auto mmd = TChem::Impl::MosaicModelData<device_type>();
    auto b_mtem = mmd.b_mtem.template view<device_type>();
    const auto b_mtem_h = Kokkos::create_mirror_view_and_copy(host_exec_space, b_mtem);
    EXPECT_TRUE(b_mtem_h(0,mmd.jnh4so4,mmd.jnh4hso4) == -4.13219);
}

}// namespace Test
}// namespace TChem

TEST(MOSAIC, device) {
  TChem::Test::MOSAIC_Test();
}

#endif
