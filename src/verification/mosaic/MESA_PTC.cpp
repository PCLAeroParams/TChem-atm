#include "TChem.hpp"
#include "TChem_Impl_MOSAIC.hpp"
#include <verification.hpp>
#include "skywalker.hpp"

using device_type = typename Tines::UseThisDevice<TChem::exec_space>::type;
using real_type_1d_view = TChem::real_type_1d_view;
using real_type_2d_view = TChem::real_type_2d_view;
using namespace skywalker;
using namespace TChem;

void MESA_PTC(Ensemble *ensemble) {
  ensemble->process([=](const Input &input, Output &output) {

    const auto aer_liquid_arr         = input.get_array("aer_liquid");
    const auto aer_solid_arr          = input.get_array("aer_solid");
    const auto aer_total_arr          = input.get_array("aer_total");
    const auto electrolyte_solid_arr  = input.get_array("electrolyte_solid");
    const auto electrolyte_liquid_arr = input.get_array("electrolyte_liquid");
    const auto electrolyte_total_arr  = input.get_array("electrolyte_total");
    const auto epercent_liquid_arr    = input.get_array("epercent_liquid");
    const auto epercent_solid_arr     = input.get_array("epercent_solid");
    const auto epercent_total_arr     = input.get_array("epercent_total");
    const auto Keq_ll_arr             = input.get_array("Keq_ll");
    const auto Keq_sl_arr             = input.get_array("Keq_sl");
    const auto log_gamZ_arr           = input.get_array("log_gamZ");
    const auto jaerosolstate_arr      = input.get_array("jaerosolstate");
    const auto jphase_arr             = input.get_array("jphase");
    const auto jhyst_leg_arr          = input.get_array("jhyst_leg");
    const auto aH2O_a_arr             = input.get_array("aH2O_a");
    const auto water_a_arr            = input.get_array("water_a");
    const auto phi_salt_old_arr       = input.get_array("phi_salt_old");
    const auto alpha_salt_arr         = input.get_array("alpha_salt");
    const auto mass_wet_a_arr         = input.get_array("mass_wet_a");

    const auto mmd = TChem::Impl::MosaicModelData<device_type>();

    real_type_1d_view aer_liquid("aer_liquid", mmd.naer);
    verification::convert_1d_vector_to_1d_view_device(aer_liquid_arr, aer_liquid);

    real_type_1d_view aer_solid("aer_solid", mmd.naer);
    verification::convert_1d_vector_to_1d_view_device(aer_solid_arr, aer_solid);

    real_type_1d_view aer_total("aer_total", mmd.naer);
    verification::convert_1d_vector_to_1d_view_device(aer_total_arr, aer_total);

    real_type_1d_view electrolyte_solid("electrolyte_solid", mmd.nelectrolyte);
    verification::convert_1d_vector_to_1d_view_device(electrolyte_solid_arr, electrolyte_solid);

    real_type_1d_view electrolyte_liquid("electrolyte_liquid", mmd.nelectrolyte);
    verification::convert_1d_vector_to_1d_view_device(electrolyte_liquid_arr, electrolyte_liquid);

    real_type_1d_view electrolyte_total("electrolyte_total", mmd.nelectrolyte);
    verification::convert_1d_vector_to_1d_view_device(electrolyte_total_arr, electrolyte_total);

    real_type_1d_view epercent_liquid("epercent_liquid", mmd.nelectrolyte);
    verification::convert_1d_vector_to_1d_view_device(epercent_liquid_arr, epercent_liquid);

    real_type_1d_view epercent_solid("epercent_solid", mmd.nelectrolyte);
    verification::convert_1d_vector_to_1d_view_device(epercent_solid_arr, epercent_solid);

    real_type_1d_view epercent_total("epercent_total", mmd.nelectrolyte);
    verification::convert_1d_vector_to_1d_view_device(epercent_total_arr, epercent_total);

    real_type_1d_view Keq_ll("Keq_ll", mmd.nrxn_aer_ll);
    verification::convert_1d_vector_to_1d_view_device(Keq_ll_arr, Keq_ll);

    real_type_1d_view Keq_sl("Keq_sl", mmd.nsalt);
    verification::convert_1d_vector_to_1d_view_device(Keq_sl_arr, Keq_sl);

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

    real_type_1d_view water_a("water_a", 1);
    verification::convert_1d_vector_to_1d_view_device(water_a_arr, water_a);

    real_type_1d_view phi_salt_old("phi_salt_old", mmd.nsalt);
    verification::convert_1d_vector_to_1d_view_device(phi_salt_old_arr, phi_salt_old);

    real_type_1d_view alpha_salt("alpha_salt", mmd.nsalt);
    verification::convert_1d_vector_to_1d_view_device(alpha_salt_arr, alpha_salt);

    real_type_1d_view mass_wet_a("mass_wet_a", 1);
    verification::convert_1d_vector_to_1d_view_device(mass_wet_a_arr, mass_wet_a);

    real_type_1d_view phi_salt("phi_salt", mmd.nsalt);
    real_type_1d_view phi_bar("phi_bar", mmd.nsalt);
    real_type_1d_view hsalt("hsalt", mmd.nsalt);
    real_type_1d_view sat_ratio("sat_ratio", mmd.nsalt);
    real_type_1d_view flux_sl("flux_sl", mmd.nsalt);
    real_type_1d_view frac_salt_solid("frac_salt_solid", mmd.nsalt);
    real_type_1d_view frac_salt_liq("frac_salt_liq", mmd.nsalt);
    real_type_1d_view eleliquid("eleliquid", mmd.nelectrolyte);
    real_type_1d_view jsalt_present("jsalt_present", mmd.nsalt);

    real_type_1d_view store("store", mmd.naer);
    real_type_1d_view na("na", mmd.nanion);
    real_type_1d_view nc("nc", mmd.ncation);
    real_type_1d_view xeq_a("xeq_a", mmd.nanion);
    real_type_1d_view xeq_c("xeq_c", mmd.ncation);
    real_type_1d_view na_Ma("na_Ma", mmd.nanion);
    real_type_1d_view nc_Mc("nc_Mc", mmd.ncation);
    real_type_1d_view aer_percent("aer_percent", mmd.naer);
    real_type_1d_view molalities("molalities", mmd.nelectrolyte);
    real_type_1d_view xmol("xmol", mmd.nelectrolyte);
    real_type_1d_view ma("ma", mmd.nanion);
    real_type_1d_view mc_view("mc_view", mmd.ncation);
    real_type_1d_view log_gam("log_gam", mmd.nelectrolyte);
    real_type_1d_view gam("gam", mmd.nelectrolyte);
    real_type_1d_view activity("activity", mmd.nelectrolyte);

    real_type_1d_view tau_p("tau_p", mmd.nsalt);
    real_type_1d_view tau_d("tau_d", mmd.nsalt);

    real_type_1d_view mass_dry_a_view("mass_dry_a", 1);
    real_type_1d_view vol_dry_a_view("vol_dry_a", 1);
    real_type_1d_view vol_wet_a_view("vol_wet_a", 1);
    real_type_1d_view growth_factor_view("growth_factor", 1);
    real_type_1d_view electrolyte_sum_solid_view("electrolyte_sum_solid", 1);

    std::string profile_name = "Verification_test_MESA_PTC";
    using policy_type =
          typename TChem::UseThisTeamPolicy<TChem::exec_space>::type;
    const auto exec_space_instance = TChem::exec_space();
    policy_type policy(exec_space_instance, 1, Kokkos::AUTO());

    Kokkos::parallel_for(
    profile_name,
    policy,
    KOKKOS_LAMBDA(const typename policy_type::member_type& member) {
      TChem::Impl::MOSAIC<real_type, device_type>::MESA_PTC(
        mmd,
        aer_liquid, aer_solid, aer_total,
        electrolyte_solid, electrolyte_liquid, electrolyte_total,
        epercent_liquid, epercent_solid, epercent_total,
        Keq_sl, Keq_ll, log_gamZ,
        jaerosolstate(0), jphase(0), jhyst_leg(0), aH2O_a(0), water_a(0),
        mass_dry_a_view(0), vol_dry_a_view(0),
        mass_wet_a(0), vol_wet_a_view(0),
        growth_factor_view(0),
        electrolyte_sum_solid_view(0),
        jsalt_present, hsalt, phi_salt, phi_salt_old, phi_bar, alpha_salt,
        sat_ratio, flux_sl, frac_salt_solid, frac_salt_liq, eleliquid,
        store, na, nc, xeq_a, xeq_c, na_Ma, nc_Mc,
        aer_percent, molalities, xmol, ma, mc_view, log_gam, gam, activity,
        tau_p, tau_d);
    });

    std::vector<real_type> aer_liquid_out(mmd.naer);
    std::vector<real_type> aer_solid_out(mmd.naer);
    std::vector<real_type> aer_total_out(mmd.naer);
    std::vector<real_type> electrolyte_solid_out(mmd.nelectrolyte);
    std::vector<real_type> electrolyte_liquid_out(mmd.nelectrolyte);
    std::vector<real_type> electrolyte_total_out(mmd.nelectrolyte);
    std::vector<real_type> epercent_liquid_out(mmd.nelectrolyte);
    std::vector<real_type> epercent_solid_out(mmd.nelectrolyte);
    std::vector<real_type> epercent_total_out(mmd.nelectrolyte);
    std::vector<real_type> Keq_ll_out(mmd.nrxn_aer_ll);
    std::vector<real_type> Keq_sl_out(mmd.nsalt);
    std::vector<real_type> log_gamZ_out(mmd.nelectrolyte * mmd.nelectrolyte);
    std::vector<real_type> jaerosolstate_out(1);
    std::vector<real_type> jphase_out(1);
    std::vector<real_type> jhyst_leg_out(1);
    std::vector<real_type> aH2O_a_out(1);
    std::vector<real_type> water_a_out(1);
    std::vector<real_type> phi_salt_old_out(mmd.nsalt);
    std::vector<real_type> alpha_salt_out(mmd.nsalt);
    std::vector<real_type> phi_salt_out(mmd.nsalt);
    std::vector<real_type> flux_sl_out(mmd.nsalt);
    std::vector<real_type> sat_ratio_out(mmd.nsalt);
    std::vector<real_type> mass_dry_a_out(1);
    std::vector<real_type> vol_dry_a_out(1);
    std::vector<real_type> mass_wet_a_out(1);
    std::vector<real_type> vol_wet_a_out(1);
    std::vector<real_type> growth_factor_out(1);
    std::vector<real_type> electrolyte_sum_solid_out(1);

    verification::convert_1d_view_device_to_1d_vector(aer_liquid, aer_liquid_out);
    verification::convert_1d_view_device_to_1d_vector(aer_solid, aer_solid_out);
    verification::convert_1d_view_device_to_1d_vector(aer_total, aer_total_out);
    verification::convert_1d_view_device_to_1d_vector(electrolyte_solid, electrolyte_solid_out);
    verification::convert_1d_view_device_to_1d_vector(electrolyte_liquid, electrolyte_liquid_out);
    verification::convert_1d_view_device_to_1d_vector(electrolyte_total, electrolyte_total_out);
    verification::convert_1d_view_device_to_1d_vector(epercent_liquid, epercent_liquid_out);
    verification::convert_1d_view_device_to_1d_vector(epercent_solid, epercent_solid_out);
    verification::convert_1d_view_device_to_1d_vector(epercent_total, epercent_total_out);
    verification::convert_1d_view_device_to_1d_vector(Keq_ll, Keq_ll_out);
    verification::convert_1d_view_device_to_1d_vector(Keq_sl, Keq_sl_out);
    verification::convert_2d_view_device_to_1d_vector(log_gamZ, log_gamZ_out);
    verification::convert_1d_view_device_to_1d_vector(jaerosolstate, jaerosolstate_out);
    verification::convert_1d_view_device_to_1d_vector(jphase, jphase_out);
    verification::convert_1d_view_device_to_1d_vector(jhyst_leg, jhyst_leg_out);
    verification::convert_1d_view_device_to_1d_vector(aH2O_a, aH2O_a_out);
    verification::convert_1d_view_device_to_1d_vector(water_a, water_a_out);
    verification::convert_1d_view_device_to_1d_vector(phi_salt_old, phi_salt_old_out);
    verification::convert_1d_view_device_to_1d_vector(alpha_salt, alpha_salt_out);
    verification::convert_1d_view_device_to_1d_vector(phi_salt, phi_salt_out);
    verification::convert_1d_view_device_to_1d_vector(flux_sl, flux_sl_out);
    verification::convert_1d_view_device_to_1d_vector(sat_ratio, sat_ratio_out);
    verification::convert_1d_view_device_to_1d_vector(mass_dry_a_view, mass_dry_a_out);
    verification::convert_1d_view_device_to_1d_vector(vol_dry_a_view, vol_dry_a_out);
    verification::convert_1d_view_device_to_1d_vector(mass_wet_a, mass_wet_a_out);
    verification::convert_1d_view_device_to_1d_vector(vol_wet_a_view, vol_wet_a_out);
    verification::convert_1d_view_device_to_1d_vector(growth_factor_view, growth_factor_out);
    verification::convert_1d_view_device_to_1d_vector(electrolyte_sum_solid_view, electrolyte_sum_solid_out);

    output.set("aer_liquid",            aer_liquid_out);
    output.set("aer_solid",             aer_solid_out);
    output.set("aer_total",             aer_total_out);
    output.set("electrolyte_solid",     electrolyte_solid_out);
    output.set("electrolyte_liquid",    electrolyte_liquid_out);
    output.set("electrolyte_total",     electrolyte_total_out);
    output.set("epercent_liquid",       epercent_liquid_out);
    output.set("epercent_solid",        epercent_solid_out);
    output.set("epercent_total",        epercent_total_out);
    output.set("Keq_ll",                Keq_ll_out);
    output.set("Keq_sl",                Keq_sl_out);
    output.set("log_gamZ",              log_gamZ_out);
    output.set("jaerosolstate",         jaerosolstate_out);
    output.set("jphase",                jphase_out);
    output.set("jhyst_leg",             jhyst_leg_out);
    output.set("aH2O_a",                aH2O_a_out);
    output.set("water_a",               water_a_out);
    output.set("phi_salt_old",          phi_salt_old_out);
    output.set("alpha_salt",            alpha_salt_out);
    output.set("phi_salt",              phi_salt_out);
    output.set("flux_sl",               flux_sl_out);
    output.set("sat_ratio",             sat_ratio_out);
    output.set("mass_dry_a",            mass_dry_a_out);
    output.set("vol_dry_a",             vol_dry_a_out);
    output.set("mass_wet_a",            mass_wet_a_out);
    output.set("vol_wet_a",             vol_wet_a_out);
    output.set("growth_factor",         growth_factor_out);
    output.set("electrolyte_sum_solid", electrolyte_sum_solid_out);
  });
}
