#ifndef __TCHEM_TEST_MOSAIC_HPP__
#define __TCHEM_TEST_MOSAIC_HPP__

#include "TChem.hpp"
#include "TChem_Impl_MOSAIC.hpp"

namespace TChem {
namespace Test {
void static MOSAIC_Test() {

    using device_type = typename Tines::UseThisDevice<TChem::exec_space>::type;
    using ordinal_type = TChem::ordinal_type;
    using real_type = TChem::real_type;
    using real_type_1d_view = TChem::real_type_1d_view;
    using real_type_2d_view = TChem::real_type_2d_view;
    using real_type_2d_view_host = TChem::real_type_2d_view_host;

    std::string profile_name ="Verification_test_MTEM_compute_log_gamZ";
    using policy_type =
          typename TChem::UseThisTeamPolicy<TChem::exec_space>::type;
    const auto exec_space_instance = TChem::exec_space();
    const auto host_exec_space = TChem::host_exec_space();
    policy_type policy(exec_space_instance, 1, Kokkos::AUTO());

    const auto mmd = TChem::Impl::MosaicModelData<device_type>();

    // Check this routines on GPUs.
     Kokkos::parallel_for(
    profile_name,
    policy,
    KOKKOS_LAMBDA(const typename policy_type::member_type& member) {

    auto b_mtem = mmd.b_mtem.template view<device_type>();

    EXPECT_TRUE(b_mtem(0,mmd.jnh4so4,mmd.jnh4hso4) == -4.13219);
    });
}

}// namespace Test
}// namespace TChem

TEST(MOSAIC, device) {
  TChem::Test::MOSAIC_Test();
}

#endif
