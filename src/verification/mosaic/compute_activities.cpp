#include "TChem.hpp"
#include "TChem_Impl_MOSAIC.hpp"
#include <verification.hpp>
#include "skywalker.hpp"

using device_type = typename Tines::UseThisDevice<TChem::exec_space>::type;
using real_type_1d_view = TChem::real_type_1d_view;
using real_type_2d_view = TChem::real_type_2d_view;
using ordinal_type = TChem::ordinal_type;
using namespace skywalker;
using namespace TChem;

void compute_activities(Ensemble *ensemble) {
  ensemble->process([=](const Input &input, Output &output) {

    const auto aer_arr = input.get_array("aer");
    const auto ma_arr = input.get_array("ma");
    const auto mc_arr = input.get_array("mc");
    const auto Keq_ll_arr = input.get_array("Keq_ll");
    const auto electrolyte_arr = input.get_array("electrolyte");
    const auto log_gam_arr = input.get_array("log_gam");
    const auto log_gamZ_arr = input.get_array("log_gamZ");
    auto jaerosolstate_arr = input.get_array("jaerosolstate");
    auto jphase_arr = input.get_array("jphase");
    auto jhyst_leg_arr = input.get_array("jhyst_leg");
    const auto aH2O_a_arr = input.get_array("aH2O_a");

    const auto mmd = TChem::Impl::MosaicModelData<device_type>();

    real_type_1d_view aer("aer", mmd.naer);
    verification::convert_1d_vector_to_1d_view_device(aer_arr, aer);

    real_type_1d_view ma("ma", mmd.nelectrolyte);
    verification::convert_1d_vector_to_1d_view_device(ma_arr, ma);

    real_type_1d_view mc("mc", mmd.nelectrolyte);
    verification::convert_1d_vector_to_1d_view_device(mc_arr, mc);

    real_type_1d_view Keq_ll("Keq_ll", mmd.nrxn_aer_ll);
    verification::convert_1d_vector_to_1d_view_device(Keq_ll_arr, Keq_ll);

    real_type_1d_view electrolyte_solid("electrolyte_solid", mmd.nelectrolyte);
    real_type_1d_view electrolyte_liquid("electrolyte_liquid", mmd.nelectrolyte);
    real_type_1d_view electrolyte_total("electrolyte_total", mmd.nelectrolyte);

    const ordinal_type n = 3;
    const auto nsize_electrolyte = static_cast<ordinal_type>(electrolyte_arr.size())/n;
    std::vector<std::vector<real_type>> electrolyte_arr_2d;
    for (size_t i = 0; i < n*nsize_electrolyte; i += nsize_electrolyte) {
        // Create a small vector from a segment of long_vector
        std::vector<real_type> electrolyte_temp(electrolyte_arr.begin() + i, electrolyte_arr.begin() + i + nsize_electrolyte);
        electrolyte_arr_2d.push_back(electrolyte_temp);
    }

    verification::convert_1d_vector_to_1d_view_device(electrolyte_arr_2d[0], electrolyte_solid);
    verification::convert_1d_vector_to_1d_view_device(electrolyte_arr_2d[1], electrolyte_liquid);
    verification::convert_1d_vector_to_1d_view_device(electrolyte_arr_2d[2], electrolyte_total);

    real_type_1d_view log_gam("log_gam", mmd.nelectrolyte);
    verification::convert_1d_vector_to_1d_view_device(log_gam_arr, log_gam);

    real_type_2d_view log_gamZ("log_gamZ", mmd.nelectrolyte, mmd.nelectrolyte);
    verification::convert_1d_vector_to_2d_view_device(log_gamZ_arr, log_gamZ);

    real_type_1d_view jaerosolstate("jaerosolstate", 1);
    verification::convert_1d_vector_to_1d_view_device(jaerosolstate_arr, jaerosolstate);

    real_type_1d_view jphase("jphase", 1);
    verification::convert_1d_vector_to_1d_view_device(jphase_arr, jphase);

    real_type_1d_view jhyst_leg("jhyst_leg", 1);
    verification::convert_1d_vector_to_1d_view_device(jhyst_leg_arr, jhyst_leg);

    real_type_1d_view aH2O_a("aH2O_a", 1);
    verification::convert_1d_vector_to_1d_view_device(aH2O_a_arr, aH2O_a);

    real_type_1d_view molalities("molalities", mmd.nelectrolyte);
    real_type_1d_view xmol("xmol", mmd.nelectrolyte);
    real_type_1d_view gam("gam", mmd.nelectrolyte);
    real_type_1d_view activity("activity", mmd.nelectrolyte);

    std::string profile_name ="Verification_test_compute_activites";
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

    // Perform the adjustment calculation
    TChem::Impl::MOSAIC<real_type, device_type>::compute_activities(
      mmd,
      molalities,
      xmol,
      aer,
      ma,
      mc,
      Keq_ll,
      electrolyte_solid, 
      electrolyte_liquid,
      electrolyte_total,
      log_gam,
      log_gamZ,
      gam,
      activity,
      jaerosolstate(0),
      jphase(0),
      jhyst_leg(0),
      aH2O_a(0));
    });

    verification::convert_1d_view_device_to_1d_vector(jaerosolstate, jaerosolstate_arr);
    output.set("jaerosolstate", jaerosolstate_arr);

    verification::convert_1d_view_device_to_1d_vector(jhyst_leg, jhyst_leg_arr);
    output.set("jhyst_leg", jhyst_leg_arr);

    verification::convert_1d_view_device_to_1d_vector(jphase, jphase_arr);
    output.set("jphase", jphase_arr);

    std::vector<real_type> molalities_arr(mmd.nelectrolyte);
    verification::convert_1d_view_device_to_1d_vector(molalities, molalities_arr);
    output.set("molalities", molalities_arr);

    std::vector<real_type> activity_arr(mmd.nelectrolyte);
    verification::convert_1d_view_device_to_1d_vector(activity, activity_arr);
    output.set("activity", activity_arr);
  });
}