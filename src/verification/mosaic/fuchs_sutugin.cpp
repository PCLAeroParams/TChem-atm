#include "TChem.hpp"
#include "TChem_Impl_MOSAIC.hpp"
#include <verification.hpp>
#include "skywalker.hpp"

using device_type = typename Tines::UseThisDevice<TChem::exec_space>::type;
using real_type_1d_view = TChem::real_type_1d_view;
using oridinal_type_1d_view = TChem::ordinal_type_1d_view;
using ordinal_type = TChem::ordinal_type;
using namespace skywalker;
using namespace TChem;

void fuchs_sutugin(Ensemble *ensemble) {
  ensemble->process([=](const Input &input, Output &output) {

    const auto rkn_arr = input.get_array("rkn");
    const auto a_arr = input.get_array("a");

    real_type_1d_view rkn("rkn", 1);
    verification::convert_1d_vector_to_1d_view_device(rkn_arr, rkn);

    real_type_1d_view a("a", 1);
    verification::convert_1d_vector_to_1d_view_device(a_arr, a);

    // Reals or int that are defined outside of the parallel_for region are passed as const.
    real_type_1d_view outputs_fuchs_sutugin("outputs_fuchs_sutugin", 1);

    std::string profile_name ="Verification_test_fuchs_sutugin";
    using policy_type =
          typename TChem::UseThisTeamPolicy<TChem::exec_space>::type;
    const auto exec_space_instance = TChem::exec_space();
    const auto host_exec_space = TChem::host_exec_space();
    policy_type policy(exec_space_instance, 1, Kokkos::AUTO());

    Kokkos::parallel_for(
    profile_name,
    policy,
    KOKKOS_LAMBDA(const typename policy_type::member_type& member) {

      Real& fuchs_sutugin = outputs_fuchs_sutugin(0);

    // Perform the adjustment calculation
    TChem::Impl::MOSAIC<real_type, device_type>::fuchs_sutugin(
      rkn(0),
      a(0),
      fuchs_sutugin);
    });

    const auto outputs_fuchs_sutugin_h = Kokkos::create_mirror_view_and_copy(host_exec_space, outputs_fuchs_sutugin);

    Real fuchs_sutugin = outputs_fuchs_sutugin_h(0);

    output.set("fuchs_sutugin", fuchs_sutugin);

  });
}
