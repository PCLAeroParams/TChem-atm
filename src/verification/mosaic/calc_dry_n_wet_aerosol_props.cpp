#include "TChem.hpp"
#include "TChem_Impl_MOSAIC.hpp"
#include <verification.hpp>
#include "skywalker.hpp"

using device_type = typename Tines::UseThisDevice<TChem::exec_space>::type;
using real_type_1d_view = TChem::real_type_1d_view;
using namespace skywalker;
using namespace TChem;

void calc_dry_n_wet_aerosol_props(Ensemble *ensemble) {
  ensemble->process([=](const Input &input, Output &output) {

    const auto aer_total_arr     = input.get_array("aer_total");
    const auto water_a_arr       = input.get_array("water_a");
    const auto num_a_arr         = input.get_array("num_a");
    const auto jaerosolstate_arr = input.get_array("jaerosolstate");

    const auto mmd = TChem::Impl::MosaicModelData<device_type>();

    real_type_1d_view aer_total("aer_total", mmd.naer);
    verification::convert_1d_vector_to_1d_view_device(aer_total_arr, aer_total);

    real_type_1d_view water_a_view("water_a", 1);
    verification::convert_1d_vector_to_1d_view_device(water_a_arr, water_a_view);

    real_type_1d_view num_a_view("num_a", 1);
    verification::convert_1d_vector_to_1d_view_device(num_a_arr, num_a_view);

    real_type_1d_view jaerosolstate_view("jaerosolstate", 1);
    verification::convert_1d_vector_to_1d_view_device(jaerosolstate_arr, jaerosolstate_view);

    real_type_1d_view mass_dry_a_view("mass_dry_a", 1);
    real_type_1d_view vol_dry_a_view ("vol_dry_a",  1);
    real_type_1d_view mass_wet_a_view("mass_wet_a", 1);
    real_type_1d_view vol_wet_a_view ("vol_wet_a",  1);
    real_type_1d_view dens_dry_a_view("dens_dry_a", 1);
    real_type_1d_view dens_wet_a_view("dens_wet_a", 1);
    real_type_1d_view Dp_dry_a_view  ("Dp_dry_a",   1);
    real_type_1d_view Dp_wet_a_view  ("Dp_wet_a",   1);
    real_type_1d_view area_dry_a_view("area_dry_a", 1);
    real_type_1d_view area_wet_a_view("area_wet_a", 1);

    std::string profile_name = "Verification_test_calc_dry_n_wet_aerosol_props";
    using policy_type =
          typename TChem::UseThisTeamPolicy<TChem::exec_space>::type;
    const auto exec_space_instance = TChem::exec_space();
    policy_type policy(exec_space_instance, 1, Kokkos::AUTO());

    Kokkos::parallel_for(
    profile_name,
    policy,
    KOKKOS_LAMBDA(const typename policy_type::member_type& member) {
      TChem::Impl::MOSAIC<real_type, device_type>::calc_dry_n_wet_aerosol_props(
        mmd,
        aer_total,
        water_a_view(0),
        num_a_view(0),
        jaerosolstate_view(0),
        mass_dry_a_view(0),
        vol_dry_a_view(0),
        mass_wet_a_view(0),
        vol_wet_a_view(0),
        dens_dry_a_view(0),
        dens_wet_a_view(0),
        Dp_dry_a_view(0),
        Dp_wet_a_view(0),
        area_dry_a_view(0),
        area_wet_a_view(0));
    });

    std::vector<real_type> aer_total_out(mmd.naer);
    std::vector<real_type> water_a_out(1), num_a_out(1), jaerosolstate_out(1);
    std::vector<real_type> mass_dry_a_out(1), vol_dry_a_out(1);
    std::vector<real_type> mass_wet_a_out(1), vol_wet_a_out(1);
    std::vector<real_type> dens_dry_a_out(1), dens_wet_a_out(1);
    std::vector<real_type> Dp_dry_a_out(1), Dp_wet_a_out(1);
    std::vector<real_type> area_dry_a_out(1), area_wet_a_out(1);

    verification::convert_1d_view_device_to_1d_vector(aer_total,         aer_total_out);
    verification::convert_1d_view_device_to_1d_vector(water_a_view,       water_a_out);
    verification::convert_1d_view_device_to_1d_vector(num_a_view,         num_a_out);
    verification::convert_1d_view_device_to_1d_vector(jaerosolstate_view, jaerosolstate_out);
    verification::convert_1d_view_device_to_1d_vector(mass_dry_a_view,    mass_dry_a_out);
    verification::convert_1d_view_device_to_1d_vector(vol_dry_a_view,     vol_dry_a_out);
    verification::convert_1d_view_device_to_1d_vector(mass_wet_a_view,    mass_wet_a_out);
    verification::convert_1d_view_device_to_1d_vector(vol_wet_a_view,     vol_wet_a_out);
    verification::convert_1d_view_device_to_1d_vector(dens_dry_a_view,    dens_dry_a_out);
    verification::convert_1d_view_device_to_1d_vector(dens_wet_a_view,    dens_wet_a_out);
    verification::convert_1d_view_device_to_1d_vector(Dp_dry_a_view,      Dp_dry_a_out);
    verification::convert_1d_view_device_to_1d_vector(Dp_wet_a_view,      Dp_wet_a_out);
    verification::convert_1d_view_device_to_1d_vector(area_dry_a_view,    area_dry_a_out);
    verification::convert_1d_view_device_to_1d_vector(area_wet_a_view,    area_wet_a_out);

    output.set("aer_total",     aer_total_out);
    output.set("water_a",       water_a_out);
    output.set("num_a",         num_a_out);
    output.set("jaerosolstate", jaerosolstate_out);
    output.set("mass_dry_a",    mass_dry_a_out);
    output.set("vol_dry_a",     vol_dry_a_out);
    output.set("mass_wet_a",    mass_wet_a_out);
    output.set("vol_wet_a",     vol_wet_a_out);
    output.set("dens_dry_a",    dens_dry_a_out);
    output.set("dens_wet_a",    dens_wet_a_out);
    output.set("Dp_dry_a",      Dp_dry_a_out);
    output.set("Dp_wet_a",      Dp_wet_a_out);
    output.set("area_dry_a",    area_dry_a_out);
    output.set("area_wet_a",    area_wet_a_out);
  });
}
