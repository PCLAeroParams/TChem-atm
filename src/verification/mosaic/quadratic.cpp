#include "TChem.hpp"
#include "TChem_Impl_MOSAIC.hpp"
#include <verification.hpp>
#include "skywalker.hpp"

using device_type = typename Tines::UseThisDevice<TChem::exec_space>::type;
using real_type_1d_view = TChem::real_type_1d_view;
using ordinal_type = TChem::ordinal_type;
using namespace skywalker;
using namespace TChem;

void quadratic(Ensemble *ensemble) {
  ensemble->process([=](const Input &input, Output &output) {

    const auto a_arr = input.get_array("a");
    const auto b_arr = input.get_array("b");
    const auto c_arr = input.get_array("c");

    real_type_1d_view a("a", 1);
    verification::convert_1d_vector_to_1d_view_device(a_arr, a);

    real_type_1d_view b("b", 1);
    verification::convert_1d_vector_to_1d_view_device(b_arr, b);

    real_type_1d_view c("c", 1);
    verification::convert_1d_vector_to_1d_view_device(c_arr, c);

    real_type_1d_view outputs_quadratic("outputs_quadratic", 1);

    std::string profile_name = "Verification_test_quadratic";
    using policy_type =
          typename TChem::UseThisTeamPolicy<TChem::exec_space>::type;
    const auto exec_space_instance = TChem::exec_space();
    const auto host_exec_space = TChem::host_exec_space();
    policy_type policy(exec_space_instance, 1, Kokkos::AUTO());

    Kokkos::parallel_for(
    profile_name,
    policy,
    KOKKOS_LAMBDA(const typename policy_type::member_type& member) {
      Real& result = outputs_quadratic(0);
      TChem::Impl::MOSAIC<real_type, device_type>::quadratic(
        a(0),
        b(0),
        c(0),
        result);
    });

    const auto outputs_quadratic_h = Kokkos::create_mirror_view_and_copy(host_exec_space, outputs_quadratic);

    Real result = outputs_quadratic_h(0);

    output.set("result", result);
  });
}
