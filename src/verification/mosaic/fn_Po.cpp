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

void fn_Po(Ensemble *ensemble) {
  ensemble->process([=](const Input &input, Output &output) {

    const auto Po_298_arr = input.get_array("Po_298");
    const auto DH_arr = input.get_array("DH");
    const auto T_arr = input.get_array("T");

    real_type_1d_view Po_298("Po_298", 1);
    verification::convert_1d_vector_to_1d_view_device(Po_298_arr, Po_298);

    real_type_1d_view DH("DH", 1);
    verification::convert_1d_vector_to_1d_view_device(DH_arr, DH);

    real_type_1d_view T("T", 1);
    verification::convert_1d_vector_to_1d_view_device(T_arr, T);

    // Reals or int that are defined outside of the parallel_for region are passed as const.
    real_type_1d_view outputs_fn_Po("outputs_fn_Po", 1);

    std::string profile_name ="Verification_test_fn_Po";
    using policy_type =
          typename TChem::UseThisTeamPolicy<TChem::exec_space>::type;
    const auto exec_space_instance = TChem::exec_space();
    const auto host_exec_space = TChem::host_exec_space();
    policy_type policy(exec_space_instance, 1, Kokkos::AUTO());

     Kokkos::parallel_for(
    profile_name,
    policy,
    KOKKOS_LAMBDA(const typename policy_type::member_type& member) {

      Real& Po = outputs_fn_Po(0);

    // Perform the adjustment calculation
    TChem::Impl::MOSAIC<real_type, device_type>::fn_Po(
      Po_298(0),
      DH(0),
      T(0),
      Po);
    });

    const auto outputs_fn_Po_h = Kokkos::create_mirror_view_and_copy(host_exec_space, outputs_fn_Po);

    Real Po = outputs_fn_Po_h(0);

    output.set("Po", Po);

  });
}
