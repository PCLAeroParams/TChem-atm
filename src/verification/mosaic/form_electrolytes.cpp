#include "TChem.hpp"
#include "TChem_Impl_MOSAIC.hpp"
#include <verification.hpp>
#include "skywalker.hpp"

using device_type = typename Tines::UseThisDevice<TChem::exec_space>::type;
using real_type_1d_view = TChem::real_type_1d_view;
using ordinal_type = TChem::ordinal_type;
using namespace skywalker;
using namespace TChem;

void form_electrolytes(Ensemble *ensemble) {
  ensemble->process([=](const Input &input, Output &output) {

    const auto jp_arr = input.get_array("jp");
    const auto aer_curr_arr = input.get_array("aer_curr");
    const auto aer_solid_arr = input.get_array("aer_solid");
    const auto aer_liquid_arr = input.get_array("aer_liquid");
    const auto aer_total_arr = input.get_array("aer_total");
    const auto gas_arr = input.get_array("gas");
    const auto store_arr = input.get_array("store");
    const auto electrolyte_arr = input.get_array("electrolyte");
    const auto total_species_arr = input.get_array("total_species");
    const auto tot_cl_in_arr = input.get_array("tot_cl_in");
    const auto electrolyte_sum_arr = input.get_array("electrolyte_sum");
    const auto epercent_arr = input.get_array("epercent");
    const auto aer_sum_arr = input.get_array("aer_sum");
    const auto aer_percent_arr = input.get_array("aer_percent");

    const auto mmd = TChem::Impl::MosaicModelData<device_type>();

    real_type_1d_view jp("jp", 1);
    verification::convert_1d_vector_to_1d_view_device(jp_arr, jp);

    real_type_1d_view aer_curr("aer_curr", mmd.naer);
    verification::convert_1d_vector_to_1d_view_device(aer_curr_arr, aer_curr);

    real_type_1d_view aer_solid("aer_solid", mmd.naer);
    verification::convert_1d_vector_to_1d_view_device(aer_solid_arr, aer_solid);

    real_type_1d_view aer_liquid("aer_liquid", mmd.naer);
    verification::convert_1d_vector_to_1d_view_device(aer_liquid_arr, aer_liquid);

    real_type_1d_view aer_total("aer_total", mmd.naer);
    verification::convert_1d_vector_to_1d_view_device(aer_total_arr, aer_total);

    real_type_1d_view gas("gas", mmd.ngas_volatile);
    verification::convert_1d_vector_to_1d_view_device(gas_arr, gas);

    real_type_1d_view store("store", mmd.naer);
    verification::convert_1d_vector_to_1d_view_device(store_arr, store);

    real_type_1d_view electrolyte("electrolyte", mmd.nelectrolyte);
    verification::convert_1d_vector_to_1d_view_device(electrolyte_arr, electrolyte);

    real_type_1d_view total_species("total_species", mmd.ngas_volatile);
    verification::convert_1d_vector_to_1d_view_device(total_species_arr, total_species);

    real_type_1d_view tot_cl_in("tot_cl_in", 1);
    verification::convert_1d_vector_to_1d_view_device(tot_cl_in_arr, tot_cl_in);

    real_type_1d_view electrolyte_sum("electrolyte_sum", 1);
    verification::convert_1d_vector_to_1d_view_device(electrolyte_sum_arr, electrolyte_sum);

    real_type_1d_view epercent("epercent", mmd.nelectrolyte);
    verification::convert_1d_vector_to_1d_view_device(epercent_arr, epercent);

    real_type_1d_view aer_sum("aer_sum", 1);
    verification::convert_1d_vector_to_1d_view_device(aer_sum_arr, aer_sum);

    real_type_1d_view aer_percent("aer_percent", mmd.naer);
    verification::convert_1d_vector_to_1d_view_device(aer_percent_arr, aer_percent);

    std::string profile_name ="Verification_test_form_electrolytes";
    using policy_type =
          typename TChem::UseThisTeamPolicy<TChem::exec_space>::type;
    const auto exec_space_instance = TChem::exec_space();
    policy_type policy(exec_space_instance, 1, Kokkos::AUTO());

    Kokkos::parallel_for(
    profile_name,
    policy,
    KOKKOS_LAMBDA(const typename policy_type::member_type& member) {
    TChem::Impl::MOSAIC<real_type, device_type>::form_electrolytes(
      mmd,
      jp(0),
      aer_curr,
      aer_solid,
      aer_liquid,
      aer_total,
      gas,
      store,
      electrolyte,
      total_species,
      tot_cl_in(0),
      electrolyte_sum(0),
      epercent,
      aer_sum(0),
      aer_percent);
    });

    std::vector<Real> aer_curr_out(mmd.naer), aer_solid_out(mmd.naer), aer_liquid_out(mmd.naer), aer_total_out(mmd.naer), gas_out(mmd.ngas_volatile), store_out(mmd.naer), electrolyte_out(mmd.nelectrolyte), total_species_out(mmd.ngas_volatile), tot_cl_in_out(1), electrolyte_sum_out(1), epercent_out(mmd.nelectrolyte), aer_sum_out(1), aer_percent_out(mmd.naer);
    verification::convert_1d_view_device_to_1d_vector(aer_curr, aer_curr_out);
    verification::convert_1d_view_device_to_1d_vector(aer_solid, aer_solid_out);
    verification::convert_1d_view_device_to_1d_vector(aer_liquid, aer_liquid_out);
    verification::convert_1d_view_device_to_1d_vector(aer_total, aer_total_out);
    verification::convert_1d_view_device_to_1d_vector(gas, gas_out);
    verification::convert_1d_view_device_to_1d_vector(store, store_out);
    verification::convert_1d_view_device_to_1d_vector(electrolyte, electrolyte_out);
    verification::convert_1d_view_device_to_1d_vector(total_species, total_species_out);
    verification::convert_1d_view_device_to_1d_vector(tot_cl_in, tot_cl_in_out);
    verification::convert_1d_view_device_to_1d_vector(electrolyte_sum, electrolyte_sum_out);
    verification::convert_1d_view_device_to_1d_vector(epercent, epercent_out);
    verification::convert_1d_view_device_to_1d_vector(aer_sum, aer_sum_out);
    verification::convert_1d_view_device_to_1d_vector(aer_percent, aer_percent_out);

    output.set("aer_curr", aer_curr_out);
    output.set("aer_solid", aer_solid_out);
    output.set("aer_liquid", aer_liquid_out);
    output.set("aer_total", aer_total_out);  
    output.set("gas", gas_out);
    output.set("store", store_out);
    output.set("electrolyte", electrolyte_out);
    output.set("total_species", total_species_out);
    output.set("tot_cl_in", tot_cl_in_out[0]);
    output.set("electrolyte_sum", electrolyte_sum_out[0]);
    output.set("epercent", epercent_out);
    output.set("aer_sum", aer_sum_out[0]);
    output.set("aer_percent", aer_percent_out);
  });
}
