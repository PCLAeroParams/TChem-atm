#include "TChem.hpp"
#include "TChem_Impl_MOSAIC.hpp"
#include <verification.hpp>
#include "skywalker.hpp"

using device_type = typename Tines::UseThisDevice<TChem::exec_space>::type;
using real_type_1d_view = TChem::real_type_1d_view;
using ordinal_type = TChem::ordinal_type;
using namespace skywalker;
using namespace TChem;

void form_nacl(Ensemble *ensemble) {
  ensemble->process([=](const Input &input, Output &output) {

    const auto jp_arr = input.get_array("jp");
    const auto electrolyte_arr = input.get_array("electrolyte");
    const auto store_arr = input.get_array("store");
    const auto aer_curr_arr = input.get_array("aer_curr");
    const auto aer_arr = input.get_array("aer");
    const auto gas_arr = input.get_array("gas");
    const auto total_species_arr = input.get_array("total_species");
    const auto tot_cl_in_arr = input.get_array("tot_cl_in");

    const auto mmd = TChem::Impl::MosaicModelData<device_type>();

    real_type_1d_view jp("jp", 1);
    verification::convert_1d_vector_to_1d_view_device(jp_arr, jp);

    real_type_1d_view store("store", mmd.naer);
    verification::convert_1d_vector_to_1d_view_device(store_arr, store);

    real_type_1d_view electrolyte("electrolyte", mmd.nelectrolyte);
    verification::convert_1d_vector_to_1d_view_device(electrolyte_arr, electrolyte);

    real_type_1d_view aer_curr("aer_curr", mmd.naer);
    verification::convert_1d_vector_to_1d_view_device(aer_curr_arr, aer_curr);

    real_type_1d_view aer_solid("aer_solid", mmd.naer);
    real_type_1d_view aer_liquid("aer_liquid", mmd.naer);
    real_type_1d_view aer_total("aer_total", mmd.naer);

    std::vector<std::vector<real_type>> aer_arr_2d;
    for (size_t i = 0; i < 3*mmd.naer; i += mmd.naer) {
        std::vector<real_type> aer_temp(aer_arr.begin() + i, aer_arr.begin() + i + mmd.naer);
        aer_arr_2d.push_back(aer_temp);
    }

    verification::convert_1d_vector_to_1d_view_device(aer_arr_2d[0], aer_solid);
    verification::convert_1d_vector_to_1d_view_device(aer_arr_2d[1], aer_liquid);
    verification::convert_1d_vector_to_1d_view_device(aer_arr_2d[2], aer_total);

    real_type_1d_view gas("gas", mmd.ngas_volatile);
    verification::convert_1d_vector_to_1d_view_device(gas_arr, gas);

    real_type_1d_view total_species("total_species", mmd.ngas_volatile);
    verification::convert_1d_vector_to_1d_view_device(total_species_arr, total_species);

    real_type_1d_view tot_cl_in("tot_cl_in", 1);
    verification::convert_1d_vector_to_1d_view_device(tot_cl_in_arr, tot_cl_in);

    std::string profile_name ="Verification_test_form_nacl";
    using policy_type =
          typename TChem::UseThisTeamPolicy<TChem::exec_space>::type;
    const auto exec_space_instance = TChem::exec_space();
    policy_type policy(exec_space_instance, 1, Kokkos::AUTO());

    Kokkos::parallel_for(
    profile_name,
    policy,
    KOKKOS_LAMBDA(const typename policy_type::member_type& member) {
    TChem::Impl::MOSAIC<real_type, device_type>::form_nacl(
      mmd,
      jp(0),
      electrolyte,
      store,
      aer_curr,
      aer_solid,
      aer_liquid,
      aer_total,
      gas,
      total_species,
      tot_cl_in(0));
    });

    std::vector<Real> electrolyte_out(mmd.nelectrolyte), store_out(mmd.naer), gas_out(mmd.ngas_volatile), total_species_out(mmd.ngas_volatile), tot_cl_in_out(1);
    verification::convert_1d_view_device_to_1d_vector(electrolyte, electrolyte_out);
    verification::convert_1d_view_device_to_1d_vector(store, store_out);
    verification::convert_1d_view_device_to_1d_vector(gas, gas_out);
    verification::convert_1d_view_device_to_1d_vector(total_species, total_species_out);
    verification::convert_1d_view_device_to_1d_vector(tot_cl_in, tot_cl_in_out);

    verification::convert_1d_view_device_to_1d_vector(aer_solid, aer_arr_2d[0]);
    verification::convert_1d_view_device_to_1d_vector(aer_liquid, aer_arr_2d[1]);
    verification::convert_1d_view_device_to_1d_vector(aer_total, aer_arr_2d[2]);

    std::vector<real_type> aer_flattened;
    for (const auto& row : aer_arr_2d) {
        aer_flattened.insert(aer_flattened.end(), row.begin(), row.end());
    }

    output.set("electrolyte", electrolyte_out);
    output.set("store", store_out);
    output.set("aer", aer_flattened);
    output.set("gas", gas_out);
    output.set("total_species", total_species_out);
    output.set("tot_cl_in", tot_cl_in_out);

  });
}