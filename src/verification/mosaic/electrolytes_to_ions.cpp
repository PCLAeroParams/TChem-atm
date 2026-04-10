#include "TChem.hpp"
#include "TChem_Impl_MOSAIC.hpp"
#include <verification.hpp>
#include "skywalker.hpp"

using device_type = typename Tines::UseThisDevice<TChem::exec_space>::type;
using real_type_1d_view = TChem::real_type_1d_view;
using ordinal_type = TChem::ordinal_type;
using namespace skywalker;
using namespace TChem;

void electrolytes_to_ions(Ensemble *ensemble) {
  ensemble->process([=](const Input &input, Output &output) {

    const auto aer_arr = input.get_array("aer");
    const auto electrolyte_arr = input.get_array("electrolyte");
    const auto aer_sum_arr = input.get_array("aer_sum");
    const auto aer_percent_arr = input.get_array("aer_percent");

    const auto mmd = TChem::Impl::MosaicModelData<device_type>();

    real_type_1d_view aer("aer", mmd.naer);
    verification::convert_1d_vector_to_1d_view_device(aer_arr, aer);

    real_type_1d_view electrolyte("electrolyte", mmd.nelectrolyte);
    verification::convert_1d_vector_to_1d_view_device(electrolyte_arr, electrolyte);

    real_type_1d_view aer_sum("aer_sum", 1);
    verification::convert_1d_vector_to_1d_view_device(aer_sum_arr, aer_sum);

    real_type_1d_view aer_percent("aer_percent", mmd.naer);
    verification::convert_1d_vector_to_1d_view_device(aer_percent_arr, aer_percent);

    std::string profile_name ="Verification_test_electrolytes_to_ions";
    using policy_type =
          typename TChem::UseThisTeamPolicy<TChem::exec_space>::type;
    const auto exec_space_instance = TChem::exec_space();
    policy_type policy(exec_space_instance, 1, Kokkos::AUTO());

    Kokkos::parallel_for(
    profile_name,
    policy,
    KOKKOS_LAMBDA(const typename policy_type::member_type& member) {
    TChem::Impl::MOSAIC<real_type, device_type>::electrolytes_to_ions(
      mmd,
      aer,
      electrolyte,
      aer_sum(0),
      aer_percent);
    });

    std::vector<Real> aer_out(mmd.naer), electrolyte_out(mmd.nelectrolyte), aer_sum_out(1), aer_percent_out(mmd.naer);
    verification::convert_1d_view_device_to_1d_vector(aer, aer_out);
    verification::convert_1d_view_device_to_1d_vector(electrolyte, electrolyte_out);
    verification::convert_1d_view_device_to_1d_vector(aer_sum, aer_sum_out);
    verification::convert_1d_view_device_to_1d_vector(aer_percent, aer_percent_out);

    output.set("aer", aer_out);
    output.set("electrolyte", electrolyte_out);
    output.set("aer_sum", aer_sum_out[0]);
    output.set("aer_percent", aer_percent_out);
  });
}
