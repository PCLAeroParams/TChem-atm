#ifndef __TCHEM_TEST_MOSAIC_HPP__
#define __TCHEM_TEST_MOSAIC_HPP__

#include "TChem.hpp"
#include "TChem_Impl_MOSAIC.hpp"
using ordinal_type = TChem::ordinal_type;
using real_type = TChem::real_type;
using real_type_1d_view = TChem::real_type_1d_view;
using real_type_2d_view = TChem::real_type_2d_view;
using real_type_2d_view_host = TChem::real_type_2d_view_host;

namespace TChem {
namespace Test {
void static MOSAIC_Test() {

    using device_type = typename Tines::UseThisDevice<TChem::exec_space>::type;
    using ordinal_type = TChem::ordinal_type;
    using real_type = TChem::real_type;
    using real_type_1d_view = TChem::real_type_1d_view;
    using real_type_2d_view = TChem::real_type_2d_view;
    using real_type_2d_view_host = TChem::real_type_2d_view_host;

    const auto exec_space_instance = TChem::exec_space();

    const auto mmd = TChem::Impl::MosaicModelData<device_type>();

    EXPECT_TRUE(mmd.b_mtem.h_view(0,mmd.jnh4so4,mmd.jnh4hso4) == -4.13219);

  //   using policy_type =
  //         typename TChem::UseThisTeamPolicy<TChem::exec_space>::type;
  //   // team policy
  //   ordinal_type nBatch =1;
  //   policy_type policy(exec_space_instance, nBatch, Kokkos::AUTO());
  //   Kokkos::parallel_for
  //     ("MOSAIC",
  //      policy,
  //      KOKKOS_LAMBDA(const typename policy_type::member_type& member) {


  // });

}
}// namespace Test

}// namespace TChem

TEST(MOSAIC, device) {
  TChem::Test::MOSAIC_Test();
}

// TEST(MOSAIC, verification_device)
// {

// }

#endif
