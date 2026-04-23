#include "TChem.hpp"
#include "TChem_Impl_MOSAIC.hpp"
#include <verification.hpp>
#include "skywalker.hpp"

using device_type = typename Tines::UseThisDevice<TChem::exec_space>::type;
using real_type_1d_view = TChem::real_type_1d_view;
using ordinal_type = TChem::ordinal_type;
using namespace skywalker;
using namespace TChem;

void absorb_tiny_nh4no3(Ensemble *ensemble) {
  ensemble->process([=](const Input &input, Output &output) {

    const auto delta_nh3_max_arr   = input.get_array("delta_nh3_max");
    const auto delta_hno3_max_arr  = input.get_array("delta_hno3_max");
    const auto electrolyte_sum_arr = input.get_array("electrolyte_sum");
    const auto aer_solid_arr       = input.get_array("aer_solid");
    const auto aer_liquid_arr      = input.get_array("aer_liquid");
    const auto aer_total_arr       = input.get_array("aer_total");
    const auto gas_arr             = input.get_array("gas");

    const auto mmd = TChem::Impl::MosaicModelData<device_type>();

    real_type delta_nh3_max_val   = delta_nh3_max_arr[0];
    real_type delta_hno3_max_val  = delta_hno3_max_arr[0];
    real_type electrolyte_sum_val = electrolyte_sum_arr[0];

    real_type_1d_view aer_solid("aer_solid", mmd.naer);
    verification::convert_1d_vector_to_1d_view_device(aer_solid_arr, aer_solid);

    real_type_1d_view aer_liquid("aer_liquid", mmd.naer);
    verification::convert_1d_vector_to_1d_view_device(aer_liquid_arr, aer_liquid);

    real_type_1d_view aer_total("aer_total", mmd.naer);
    verification::convert_1d_vector_to_1d_view_device(aer_total_arr, aer_total);

    real_type_1d_view gas("gas", mmd.ngas_volatile);
    verification::convert_1d_vector_to_1d_view_device(gas_arr, gas);

    real_type_1d_view delta_nh3_max_view("delta_nh3_max", 1);
    real_type_1d_view delta_hno3_max_view("delta_hno3_max", 1);
    real_type_1d_view electrolyte_sum_view("electrolyte_sum", 1);
    Kokkos::deep_copy(delta_nh3_max_view, delta_nh3_max_val);
    Kokkos::deep_copy(delta_hno3_max_view, delta_hno3_max_val);
    Kokkos::deep_copy(electrolyte_sum_view, electrolyte_sum_val);

    std::string profile_name = "Verification_test_absorb_tiny_nh4no3";
    using policy_type =
          typename TChem::UseThisTeamPolicy<TChem::exec_space>::type;
    const auto exec_space_instance = TChem::exec_space();
    policy_type policy(exec_space_instance, 1, Kokkos::AUTO());

    Kokkos::parallel_for(
    profile_name,
    policy,
    KOKKOS_LAMBDA(const typename policy_type::member_type& member) {
      TChem::Impl::MOSAIC<real_type, device_type>::absorb_tiny_nh4no3(
        mmd,
        aer_solid,
        aer_liquid,
        aer_total,
        gas,
        delta_nh3_max_view(0),
        delta_hno3_max_view(0),
        electrolyte_sum_view(0));
    });

    std::vector<Real> aer_solid_out(mmd.naer), aer_liquid_out(mmd.naer),
                      aer_total_out(mmd.naer), gas_out(mmd.ngas_volatile),
                      delta_nh3_max_out(1), delta_hno3_max_out(1), electrolyte_sum_out(1);
    verification::convert_1d_view_device_to_1d_vector(aer_solid, aer_solid_out);
    verification::convert_1d_view_device_to_1d_vector(aer_liquid, aer_liquid_out);
    verification::convert_1d_view_device_to_1d_vector(aer_total, aer_total_out);
    verification::convert_1d_view_device_to_1d_vector(gas, gas_out);
    verification::convert_1d_view_device_to_1d_vector(delta_nh3_max_view, delta_nh3_max_out);
    verification::convert_1d_view_device_to_1d_vector(delta_hno3_max_view, delta_hno3_max_out);
    verification::convert_1d_view_device_to_1d_vector(electrolyte_sum_view, electrolyte_sum_out);

    output.set("aer_solid", aer_solid_out);
    output.set("aer_liquid", aer_liquid_out);
    output.set("aer_total", aer_total_out);
    output.set("gas", gas_out);
    output.set("delta_nh3_max", delta_nh3_max_out);
    output.set("delta_hno3_max", delta_hno3_max_out);
    output.set("electrolyte_sum", electrolyte_sum_out);
  });
}
