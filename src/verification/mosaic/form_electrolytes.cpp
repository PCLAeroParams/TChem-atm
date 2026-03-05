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

    const auto jp_arr              = input.get_array("jp");
    const auto XT_arr              = input.get_array("XT");
    const auto aer_db              = input.get_array("aer");
    const auto electrolyte_arr     = input.get_array("electrolyte");
    const auto epercent_arr        = input.get_array("epercent");
    const auto aer_sum_arr         = input.get_array("aer_sum");
    const auto aer_percent_arr     = input.get_array("aer_percent");
    const auto electrolyte_sum_arr = input.get_array("electrolyte_sum");
    const auto gas_arr             = input.get_array("gas");

    const auto mmd = TChem::Impl::MosaicModelData<device_type>();

    // jp is a scalar stored as a length-1 array
    real_type_1d_view jp("jp", 1);
    verification::convert_1d_vector_to_1d_view_device(jp_arr, jp);

    // XT is a scalar stored as a length-1 array
    real_type_1d_view XT("XT", 1);
    verification::convert_1d_vector_to_1d_view_device(XT_arr, XT);

    // aer(:,:,ibin) arrives as a flat array of 3*naer values (solid, liquid, total)
    const ordinal_type n_phases = 3;
    const auto nsize_aero = static_cast<ordinal_type>(aer_db.size()) / n_phases;

    real_type_1d_view aer_solid ("aer_solid",  nsize_aero);
    real_type_1d_view aer_liquid("aer_liquid", nsize_aero);
    real_type_1d_view aer_total ("aer_total",  nsize_aero);

    std::vector<std::vector<real_type>> aer_db_2d;
    for (size_t i = 0; i < n_phases * nsize_aero; i += nsize_aero) {
      std::vector<real_type> tmp(aer_db.begin() + i, aer_db.begin() + i + nsize_aero);
      aer_db_2d.push_back(tmp);
    }
    verification::convert_1d_vector_to_1d_view_device(aer_db_2d[0], aer_solid);
    verification::convert_1d_vector_to_1d_view_device(aer_db_2d[1], aer_liquid);
    verification::convert_1d_vector_to_1d_view_device(aer_db_2d[2], aer_total);

    // electrolyte, epercent, aer_percent are for the jp phase (always liquid)
    real_type_1d_view electrolyte("electrolyte", mmd.nelectrolyte);
    verification::convert_1d_vector_to_1d_view_device(electrolyte_arr, electrolyte);

    real_type_1d_view epercent("epercent", mmd.nelectrolyte);
    verification::convert_1d_vector_to_1d_view_device(epercent_arr, epercent);

    real_type_1d_view aer_sum("aer_sum", 1);
    verification::convert_1d_vector_to_1d_view_device(aer_sum_arr, aer_sum);

    real_type_1d_view aer_percent("aer_percent", mmd.naer);
    verification::convert_1d_vector_to_1d_view_device(aer_percent_arr, aer_percent);

    real_type_1d_view electrolyte_sum("electrolyte_sum", 1);
    verification::convert_1d_vector_to_1d_view_device(electrolyte_sum_arr, electrolyte_sum);

    real_type_1d_view gas("gas", mmd.ngas_volatile);
    verification::convert_1d_vector_to_1d_view_device(gas_arr, gas);

    // total_species and tot_cl_in are needed only when form_nacl hits Cl-deficit;
    // allocate and zero-initialise so the function is safe even if unused.
    real_type_1d_view total_species("total_species", mmd.ngas_volatile);
    real_type_1d_view tot_cl_in_view("tot_cl_in", 1);

    std::string profile_name = "Verification_test_form_electrolytes";
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
          aer_liquid,     // aer(:,jp=jliquid,ibin)
          aer_solid,
          aer_liquid,
          aer_total,
          electrolyte,
          epercent,
          aer_sum(0),
          aer_percent,
          electrolyte_sum(0),
          XT(0),
          gas,
          total_species,
          tot_cl_in_view(0));
      });

    // Collect outputs
    std::vector<Real> XT_out(1);
    verification::convert_1d_view_device_to_1d_vector(XT, XT_out);

    verification::convert_1d_view_device_to_1d_vector(aer_solid,  aer_db_2d[0]);
    verification::convert_1d_view_device_to_1d_vector(aer_liquid, aer_db_2d[1]);
    verification::convert_1d_view_device_to_1d_vector(aer_total,  aer_db_2d[2]);
    std::vector<real_type> aer_flattened;
    for (const auto& row : aer_db_2d)
      aer_flattened.insert(aer_flattened.end(), row.begin(), row.end());

    std::vector<Real> electrolyte_out(mmd.nelectrolyte);
    verification::convert_1d_view_device_to_1d_vector(electrolyte, electrolyte_out);

    std::vector<Real> epercent_out(mmd.nelectrolyte);
    verification::convert_1d_view_device_to_1d_vector(epercent, epercent_out);

    std::vector<Real> aer_sum_out(1);
    verification::convert_1d_view_device_to_1d_vector(aer_sum, aer_sum_out);

    std::vector<Real> aer_percent_out(mmd.naer);
    verification::convert_1d_view_device_to_1d_vector(aer_percent, aer_percent_out);

    std::vector<Real> electrolyte_sum_out(1);
    verification::convert_1d_view_device_to_1d_vector(electrolyte_sum, electrolyte_sum_out);

    std::vector<Real> gas_out(mmd.ngas_volatile);
    verification::convert_1d_view_device_to_1d_vector(gas, gas_out);

    output.set("XT",              XT_out[0]);
    output.set("aer",             aer_flattened);
    output.set("electrolyte",     electrolyte_out);
    output.set("epercent",        epercent_out);
    output.set("aer_sum",         aer_sum_out[0]);
    output.set("aer_percent",     aer_percent_out);
    output.set("electrolyte_sum", electrolyte_sum_out[0]);
    output.set("gas",             gas_out);
  });
}
