#ifndef __TCHEM_TEST_MESA_HPP__
#define __TCHEM_TEST_MESA_HPP__

#include "TChem.hpp"
#include "TChem_Impl_MESA.hpp"
using ordinal_type = TChem::ordinal_type;
using real_type = TChem::real_type;
using real_type_1d_view = TChem::real_type_1d_view;
using real_type_2d_view = TChem::real_type_2d_view;
using real_type_2d_view_host = TChem::real_type_2d_view_host;

namespace TChem {
namespace Test {
void static MESA_Test()
{
    using value_type = real_type;
    using device_type = typename Tines::UseThisDevice<TChem::exec_space>::type;
    using MESA_type = TChem::Impl::MESA<value_type, device_type >;
    using value_type_1d_view_type = typename MESA_type::value_type_1d_view_type;
    using real_type_1d_view_type = typename MESA_type::real_type_1d_view_type;

    const auto exec_space_instance = TChem::exec_space();

  //   using policy_type =
  //         typename TChem::UseThisTeamPolicy<TChem::exec_space>::type;
  //   // team policy
  //   ordinal_type nBatch =1;
  //   policy_type policy(exec_space_instance, nBatch, Kokkos::AUTO());
  //   Kokkos::parallel_for
  //     ("MESA",
  //      policy,
  //      KOKKOS_LAMBDA(const typename policy_type::member_type& member) {


  // });

}
}// namespace Test

}// namespace TChem

TEST(MESA, Device)
{
  TChem::Test::MESA_Test();
}

// TEST(MESA, verification_device)
// {

// }

#endif
