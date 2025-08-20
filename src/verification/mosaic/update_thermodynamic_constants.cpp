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

void update_thermodynamic_constants(Ensemble *ensemble) {
  ensemble->process([=](const Input &input, Output &output) {

    const auto mmd = TChem::Impl::MosaicModelData<device_type>();

    const auto T_K_arr = input.get_array("T_K");
    const auto aH2O_arr = input.get_array("aH2O");
    const auto log_gamZ_arr = input.get_array("log_gamZ");
    auto Keq_gl_arr = input.get_array("Keq_gl");
    auto Keq_ll_arr = input.get_array("Keq_ll");
    auto Kp_nh3_arr = input.get_array("Kp_nh3");
    auto Kp_nh4no3_arr = input.get_array("Kp_nh4no3");
    auto Kp_nh4cl_arr = input.get_array("Kp_nh4cl");
    auto Keq_sg_arr = input.get_array("Keq_sg");
    auto Keq_sl_arr = input.get_array("Keq_sl");
    auto Po_soa_arr = input.get_array("Po_soa");
    auto sat_soa_arr = input.get_array("sat_soa");
    auto sigma_water_arr = input.get_array("sigma_water");
    auto MDRH_T_arr = input.get_array("MDRH_T");
    auto Kp_nh4no3_0_arr = input.get_array("Kp_nh4no3_0");
    auto Kp_nh4cl_0_arr = input.get_array("Kp_nh4cl_0");

    real_type_1d_view T_K("T_K", 1);
    verification::convert_1d_vector_to_1d_view_device(T_K_arr, T_K);

    real_type_1d_view aH2O("aH2O", 1);
    verification::convert_1d_vector_to_1d_view_device(aH2O_arr, aH2O);

    real_type_2d_view log_gamZ("log_gamZ", mmd.nelectrolyte, mmd.nelectrolyte);
    verification::convert_1d_vector_to_2d_view_device(log_gamZ_arr, log_gamZ);

    real_type_1d_view Keq_gl("Keq_gl", mmd.nrxn_aer_gl);
    verification::convert_1d_vector_to_1d_view_device(Keq_gl_arr, Keq_gl);

    real_type_1d_view Keq_ll("Keq_ll", mmd.nrxn_aer_ll);
    verification::convert_1d_vector_to_1d_view_device(Keq_ll_arr, Keq_ll);

    real_type_1d_view Kp_nh3("Kp_nh3", 1);
    verification::convert_1d_vector_to_1d_view_device(Kp_nh3_arr, Kp_nh3);

    real_type_1d_view Kp_nh4no3("Kp_nh4no3", 1);
    verification::convert_1d_vector_to_1d_view_device(Kp_nh4no3_arr, Kp_nh4no3);

    real_type_1d_view Kp_nh4cl("Kp_nh4cl", 1);
    verification::convert_1d_vector_to_1d_view_device(Kp_nh4cl_arr, Kp_nh4cl);

    real_type_1d_view Keq_sg("Keq_sg", mmd.nrxn_aer_sg);
    verification::convert_1d_vector_to_1d_view_device(Keq_sg_arr, Keq_sg);

    real_type_1d_view Keq_sl("Keq_sl", mmd.nrxn_aer_sl);
    verification::convert_1d_vector_to_1d_view_device(Keq_sl_arr, Keq_sl);
    
    real_type_1d_view Po_soa("Po_soa", mmd.ngas_volatile);
    verification::convert_1d_vector_to_1d_view_device(Po_soa_arr, Po_soa);

    real_type_1d_view sat_soa("sat_soa", mmd.ngas_volatile);
    verification::convert_1d_vector_to_1d_view_device(sat_soa_arr, sat_soa);

    real_type_1d_view sigma_water("sigma_water", 1);
    verification::convert_1d_vector_to_1d_view_device(sigma_water_arr, sigma_water);

    real_type_1d_view MDRH_T("MDRH_T", 63);
    verification::convert_1d_vector_to_1d_view_device(MDRH_T_arr, MDRH_T);

    real_type_1d_view Kp_nh4no3_0("Kp_nh4no3_0", 1);
    verification::convert_1d_vector_to_1d_view_device(Kp_nh4no3_0_arr, Kp_nh4no3_0);

    real_type_1d_view Kp_nh4cl_0("Kp_nh4cl_0", 1);
    verification::convert_1d_vector_to_1d_view_device(Kp_nh4cl_0_arr, Kp_nh4cl_0);

    std::string profile_name ="Verification_test_update_thermodynamic_constants";
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
    TChem::Impl::MOSAIC<real_type, device_type>::update_thermodynamic_constants(
      mmd, T_K(0), aH2O(0), log_gamZ,
      Keq_gl, Keq_ll, Kp_nh3(0), Kp_nh4no3(0),
      Kp_nh4cl(0), Keq_sg, Keq_sl, Po_soa,
      sat_soa, sigma_water(0), MDRH_T,
      Kp_nh4no3_0(0), Kp_nh4cl_0(0));
    });

    verification::convert_1d_view_device_to_1d_vector(Keq_gl, Keq_gl_arr);
    verification::convert_1d_view_device_to_1d_vector(Keq_ll, Keq_ll_arr);
    verification::convert_1d_view_device_to_1d_vector(Kp_nh3, Kp_nh3_arr);
    verification::convert_1d_view_device_to_1d_vector(Kp_nh4no3, Kp_nh4no3_arr);
    verification::convert_1d_view_device_to_1d_vector(Kp_nh4cl, Kp_nh4cl_arr);
    verification::convert_1d_view_device_to_1d_vector(Keq_sg, Keq_sg_arr);
    verification::convert_1d_view_device_to_1d_vector(Keq_sl, Keq_sl_arr);
    verification::convert_1d_view_device_to_1d_vector(Po_soa, Po_soa_arr);
    verification::convert_1d_view_device_to_1d_vector(sat_soa, sat_soa_arr);
    verification::convert_1d_view_device_to_1d_vector(sigma_water, sigma_water_arr);
    verification::convert_1d_view_device_to_1d_vector(MDRH_T, MDRH_T_arr);
    verification::convert_1d_view_device_to_1d_vector(Kp_nh4no3_0, Kp_nh4no3_0_arr);
    verification::convert_1d_view_device_to_1d_vector(Kp_nh4cl_0, Kp_nh4cl_0_arr);

    output.set("Keq_gl", Keq_gl_arr);
    output.set("Keq_ll", Keq_ll_arr);
    output.set("Kp_nh3", Kp_nh3_arr);
    output.set("Kp_nh4no3", Kp_nh4no3_arr);
    output.set("Kp_nh4cl", Kp_nh4cl_arr);
    output.set("Keq_sg", Keq_sg_arr);
    output.set("Keq_sl", Keq_sl_arr);
    output.set("Po_soa", Po_soa_arr);
    output.set("sat_soa", sat_soa_arr);
    output.set("sigma_water", sigma_water_arr);
    output.set("MDRH_T", MDRH_T_arr);
    output.set("Kp_nh4no3_0", Kp_nh4no3_0_arr);
    output.set("Kp_nh4cl_0", Kp_nh4cl_0_arr);

  });
}
