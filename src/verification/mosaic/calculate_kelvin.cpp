#include "TChem.hpp"
#include "TChem_Impl_MOSAIC.hpp"
#include <verification.hpp>
#include "skywalker.hpp"

using device_type = typename Tines::UseThisDevice<TChem::exec_space>::type;
using real_type_1d_view = TChem::real_type_1d_view;
using namespace skywalker;
using namespace TChem;

void calculate_kelvin(Ensemble *ensemble) {
  ensemble->process([=](const Input &input, Output &output) {

    const auto vol_wet_a_arr  = input.get_array("vol_wet_a");
    const auto num_a_arr      = input.get_array("num_a");
    const auto aH2O_a_arr     = input.get_array("aH2O_a");
    const auto sigma_water_arr= input.get_array("sigma_water");
    const auto T_K_arr        = input.get_array("T_K");

    const auto mmd = TChem::Impl::MosaicModelData<device_type>();

    real_type_1d_view vol_wet_a_view("vol_wet_a",   1);
    real_type_1d_view num_a_view("num_a",           1);
    real_type_1d_view aH2O_a_view("aH2O_a",         1);
    real_type_1d_view sigma_water_view("sigma_water",1);
    real_type_1d_view T_K_view("T_K",               1);

    verification::convert_1d_vector_to_1d_view_device(vol_wet_a_arr,   vol_wet_a_view);
    verification::convert_1d_vector_to_1d_view_device(num_a_arr,       num_a_view);
    verification::convert_1d_vector_to_1d_view_device(aH2O_a_arr,      aH2O_a_view);
    verification::convert_1d_vector_to_1d_view_device(sigma_water_arr, sigma_water_view);
    verification::convert_1d_vector_to_1d_view_device(T_K_arr,         T_K_view);

    real_type_1d_view volume_a_view("volume_a",   1);
    real_type_1d_view DpmV_view("DpmV",           1);
    real_type_1d_view sigma_soln_view("sigma_soln",1);
    real_type_1d_view kelvin_view("kelvin",        1);

    std::string profile_name = "Verification_test_calculate_kelvin";
    using policy_type =
          typename TChem::UseThisTeamPolicy<TChem::exec_space>::type;
    const auto exec_space_instance = TChem::exec_space();
    policy_type policy(exec_space_instance, 1, Kokkos::AUTO());

    Kokkos::parallel_for(
    profile_name,
    policy,
    KOKKOS_LAMBDA(const typename policy_type::member_type& member) {
      TChem::Impl::MOSAIC<real_type, device_type>::calculate_kelvin(
        vol_wet_a_view(0),
        num_a_view(0),
        aH2O_a_view(0),
        sigma_water_view(0),
        T_K_view(0),
        volume_a_view(0),
        DpmV_view(0),
        sigma_soln_view(0),
        kelvin_view(0));
    });

    std::vector<real_type> vol_wet_a_out(1);
    std::vector<real_type> num_a_out(1);
    std::vector<real_type> aH2O_a_out(1);
    std::vector<real_type> sigma_water_out(1);
    std::vector<real_type> T_K_out(1);
    std::vector<real_type> volume_a_out(1);
    std::vector<real_type> DpmV_out(1);
    std::vector<real_type> sigma_soln_out(1);
    std::vector<real_type> kelvin_out(1);

    verification::convert_1d_view_device_to_1d_vector(vol_wet_a_view,   vol_wet_a_out);
    verification::convert_1d_view_device_to_1d_vector(num_a_view,       num_a_out);
    verification::convert_1d_view_device_to_1d_vector(aH2O_a_view,      aH2O_a_out);
    verification::convert_1d_view_device_to_1d_vector(sigma_water_view, sigma_water_out);
    verification::convert_1d_view_device_to_1d_vector(T_K_view,         T_K_out);
    verification::convert_1d_view_device_to_1d_vector(volume_a_view,    volume_a_out);
    verification::convert_1d_view_device_to_1d_vector(DpmV_view,        DpmV_out);
    verification::convert_1d_view_device_to_1d_vector(sigma_soln_view,  sigma_soln_out);
    verification::convert_1d_view_device_to_1d_vector(kelvin_view,      kelvin_out);

    output.set("vol_wet_a",   vol_wet_a_out);
    output.set("num_a",       num_a_out);
    output.set("aH2O_a",      aH2O_a_out);
    output.set("sigma_water", sigma_water_out);
    output.set("T_K",         T_K_out);
    output.set("volume_a",    volume_a_out);
    output.set("DpmV",        DpmV_out);
    output.set("sigma_soln",  sigma_soln_out);
    output.set("kelvin",      kelvin_out);
  });
}
