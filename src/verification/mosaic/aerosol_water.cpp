#include "TChem.hpp"
#include "TChem_Impl_MOSAIC.hpp"
#include <verification.hpp>
#include "skywalker.hpp"

using device_type = typename Tines::UseThisDevice<TChem::exec_space>::type;
using real_type_1d_view = TChem::real_type_1d_view;
using ordinal_type = TChem::ordinal_type;
using namespace skywalker;
using namespace TChem;

void aerosol_water(Ensemble *ensemble) {
  ensemble->process([=](const Input &input, Output &output) {

    const auto electrolyte_arr = input.get_array("electrolyte");
    const auto aH2O_a_arr = input.get_array("aH2O_a");
    auto jaerosolstate_arr = input.get_array("jaerosolstate");
    auto jphase_arr = input.get_array("jphase");
    auto jhyst_leg_arr = input.get_array("jhyst_leg");

    const auto mmd = TChem::Impl::MosaicModelData<device_type>();

    real_type_1d_view electrolyte("electrolyte", mmd.nelectrolyte);
    verification::convert_1d_vector_to_1d_view_device(electrolyte_arr, electrolyte);

    real_type_1d_view aH2O_a("aH2O_a", 1);
    verification::convert_1d_vector_to_1d_view_device(aH2O_a_arr, aH2O_a);

    real_type_1d_view jaerosolstate("jaerosolstate", 1);
    verification::convert_1d_vector_to_1d_view_device(jaerosolstate_arr, jaerosolstate);

    real_type_1d_view jphase("jphase", 1);
    verification::convert_1d_vector_to_1d_view_device(jphase_arr, jphase);

    real_type_1d_view jhyst_leg("jhyst_leg", 1);
    verification::convert_1d_vector_to_1d_view_device(jhyst_leg_arr, jhyst_leg);

    real_type_1d_view outputs_aerosol_water("outputs_aerosol_water", 1);

    real_type_1d_view molalities("molalities", mmd.nelectrolyte);

    std::string profile_name ="Verification_test_aerosol_water";
    using policy_type =
          typename TChem::UseThisTeamPolicy<TChem::exec_space>::type;
    const auto exec_space_instance = TChem::exec_space();
    const auto host_exec_space = TChem::host_exec_space();
    policy_type policy(exec_space_instance, 1, Kokkos::AUTO());

    // Check this routines on GPUs.
    Kokkos::parallel_for(
    profile_name,
    policy,
    KOKKOS_LAMBDA(const typename policy_type::member_type& member) {

      Real& aerosol_water = outputs_aerosol_water(0);

    // Perform the adjustment calculation
    TChem::Impl::MOSAIC<real_type, device_type>::aerosol_water(
      mmd,
      electrolyte,
      aH2O_a(0),
      molalities,
      jaerosolstate(0),
      jphase(0),
      jhyst_leg(0),
      aerosol_water);
    });

    const auto outputs_aerosol_water_h = Kokkos::create_mirror_view_and_copy(host_exec_space, outputs_aerosol_water);

    Real aerosol_water = outputs_aerosol_water_h(0);

    output.set("aerosol_water", aerosol_water);

    verification::convert_1d_view_device_to_1d_vector(jaerosolstate, jaerosolstate_arr);
    output.set("jaerosolstate", jaerosolstate_arr);

    verification::convert_1d_view_device_to_1d_vector(jhyst_leg, jhyst_leg_arr);
    output.set("jhyst_leg", jhyst_leg_arr);

    verification::convert_1d_view_device_to_1d_vector(jphase, jphase_arr);
    output.set("jphase", jphase_arr);

    std::vector<real_type> molalities_arr(mmd.nelectrolyte);
    verification::convert_1d_view_device_to_1d_vector(molalities, molalities_arr);
    output.set("molalities", molalities_arr);
  });
}