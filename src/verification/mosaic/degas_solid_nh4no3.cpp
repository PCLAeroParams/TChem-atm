#include "TChem.hpp"
#include "TChem_Impl_MOSAIC.hpp"
#include <verification.hpp>
#include "skywalker.hpp"

using device_type = typename Tines::UseThisDevice<TChem::exec_space>::type;
using real_type_1d_view = TChem::real_type_1d_view;
using ordinal_type = TChem::ordinal_type;
using namespace skywalker;
using namespace TChem;

void degas_solid_nh4no3(Ensemble *ensemble) {
  ensemble->process([=](const Input &input, Output &output) {

    const auto gas_arr               = input.get_array("gas");
    const auto Keq_sg_arr            = input.get_array("Keq_sg");
    const auto electrolyte_solid_arr = input.get_array("electrolyte_solid");
    const auto aer_solid_arr         = input.get_array("aer_solid");
    const auto aer_liquid_arr        = input.get_array("aer_liquid");
    const auto aer_total_arr         = input.get_array("aer_total");

    const auto mmd = TChem::Impl::MosaicModelData<device_type>();

    real_type_1d_view gas("gas", mmd.ngas_volatile);
    verification::convert_1d_vector_to_1d_view_device(gas_arr, gas);

    real_type_1d_view Keq_sg("Keq_sg", mmd.nrxn_aer_sg);
    verification::convert_1d_vector_to_1d_view_device(Keq_sg_arr, Keq_sg);

    real_type_1d_view electrolyte_solid("electrolyte_solid", mmd.nelectrolyte);
    verification::convert_1d_vector_to_1d_view_device(electrolyte_solid_arr, electrolyte_solid);

    real_type_1d_view aer_solid("aer_solid", mmd.naer);
    verification::convert_1d_vector_to_1d_view_device(aer_solid_arr, aer_solid);

    real_type_1d_view aer_liquid("aer_liquid", mmd.naer);
    verification::convert_1d_vector_to_1d_view_device(aer_liquid_arr, aer_liquid);

    real_type_1d_view aer_total("aer_total", mmd.naer);
    verification::convert_1d_vector_to_1d_view_device(aer_total_arr, aer_total);

    std::string profile_name = "Verification_test_degas_solid_nh4no3";
    using policy_type =
          typename TChem::UseThisTeamPolicy<TChem::exec_space>::type;
    const auto exec_space_instance = TChem::exec_space();
    policy_type policy(exec_space_instance, 1, Kokkos::AUTO());

    Kokkos::parallel_for(
    profile_name,
    policy,
    KOKKOS_LAMBDA(const typename policy_type::member_type& member) {
      TChem::Impl::MOSAIC<real_type, device_type>::degas_solid_nh4no3(
        mmd,
        gas,
        Keq_sg,
        electrolyte_solid,
        aer_solid,
        aer_liquid,
        aer_total);
    });

    std::vector<real_type> gas_out(mmd.ngas_volatile);
    std::vector<real_type> Keq_sg_out(mmd.nrxn_aer_sg);
    std::vector<real_type> electrolyte_solid_out(mmd.nelectrolyte);
    std::vector<real_type> aer_solid_out(mmd.naer);
    std::vector<real_type> aer_liquid_out(mmd.naer);
    std::vector<real_type> aer_total_out(mmd.naer);

    verification::convert_1d_view_device_to_1d_vector(gas, gas_out);
    verification::convert_1d_view_device_to_1d_vector(Keq_sg, Keq_sg_out);
    verification::convert_1d_view_device_to_1d_vector(electrolyte_solid, electrolyte_solid_out);
    verification::convert_1d_view_device_to_1d_vector(aer_solid, aer_solid_out);
    verification::convert_1d_view_device_to_1d_vector(aer_liquid, aer_liquid_out);
    verification::convert_1d_view_device_to_1d_vector(aer_total, aer_total_out);

    output.set("gas", gas_out);
    output.set("Keq_sg", Keq_sg_out);
    output.set("electrolyte_solid", electrolyte_solid_out);
    output.set("aer_solid", aer_solid_out);
    output.set("aer_liquid", aer_liquid_out);
    output.set("aer_total", aer_total_out);

  });
}
