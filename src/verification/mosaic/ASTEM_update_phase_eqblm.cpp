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

void ASTEM_update_phase_eqblm(Ensemble *ensemble) {
  ensemble->process([=](const Input &input, Output &output) {

    const auto mmd = TChem::Impl::MosaicModelData<device_type>();

    const auto read_view = [&](const std::string& name, const ordinal_type n) {
      real_type_1d_view v(name, n);
      verification::convert_1d_vector_to_1d_view_device(input.get_array(name), v);
      return v;
    };
    const auto read_view_2d = [&](const std::string& name, const ordinal_type n0, const ordinal_type n1) {
      real_type_2d_view v(name, n0, n1);
      verification::convert_1d_vector_to_2d_view_device(input.get_array(name), v);
      return v;
    };

    const auto aer_liquid = read_view("aer_liquid", mmd.naer);
    const auto aer_solid = read_view("aer_solid", mmd.naer);
    const auto aer_total = read_view("aer_total", mmd.naer);
    const auto electrolyte_solid = read_view("electrolyte_solid", mmd.nelectrolyte);
    const auto electrolyte_liquid = read_view("electrolyte_liquid", mmd.nelectrolyte);
    const auto electrolyte_total = read_view("electrolyte_total", mmd.nelectrolyte);
    const auto epercent_liquid = read_view("epercent_liquid", mmd.nelectrolyte);
    const auto epercent_solid = read_view("epercent_solid", mmd.nelectrolyte);
    const auto epercent_total = read_view("epercent_total", mmd.nelectrolyte);
    const auto Keq_ll = read_view("Keq_ll", mmd.nrxn_aer_ll);
    const auto Keq_sl = read_view("Keq_sl", mmd.nrxn_aer_sl);
    const auto log_gamZ = read_view_2d("log_gamZ", mmd.nelectrolyte, mmd.nelectrolyte);
    const auto MDRH_T = read_view("MDRH_T", 63);
    const auto aH2O = read_view("aH2O", 1);
    const auto jaerosolstate = read_view("jaerosolstate", 1);
    const auto jphase = read_view("jphase", 1);
    const auto jhyst_leg = read_view("jhyst_leg", 1);
    const auto aH2O_a = read_view("aH2O_a", 1);
    const auto water_a = read_view("water_a", 1);
    const auto mass_dry_a = read_view("mass_dry_a", 1);
    const auto vol_dry_a = read_view("vol_dry_a", 1);
    const auto mass_wet_a = read_view("mass_wet_a", 1);
    const auto vol_wet_a = read_view("vol_wet_a", 1);
    const auto growth_factor = read_view("growth_factor", 1);
    const auto electrolyte_sum_solid = read_view("electrolyte_sum_solid", 1);
    const auto MDRH = read_view("MDRH", 1);
    const auto jsalt_present = read_view("jsalt_present", mmd.nsalt);
    const auto phi_salt_old = read_view("phi_salt_old", mmd.nsalt);
    const auto alpha_salt = read_view("alpha_salt", mmd.nsalt);
    const auto phi_salt = read_view("phi_salt", mmd.nsalt);
    const auto flux_sl = read_view("flux_sl", mmd.nsalt);
    const auto sat_ratio = read_view("sat_ratio", mmd.nsalt);
    const auto ma = read_view("ma", mmd.nanion);
    const auto mc_view = read_view("mc", mmd.ncation);
    const auto gam = read_view("gam", mmd.nelectrolyte);
    const auto activity = read_view("activity", mmd.nelectrolyte);

    real_type_1d_view hsalt("hsalt", mmd.nsalt);
    real_type_1d_view phi_bar("phi_bar", mmd.nsalt);
    real_type_1d_view frac_salt_solid("frac_salt_solid", mmd.nsalt);
    real_type_1d_view frac_salt_liq("frac_salt_liq", mmd.nsalt);
    real_type_1d_view eleliquid("eleliquid", mmd.nelectrolyte);
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
    real_type_1d_view log_gam("log_gam", mmd.nelectrolyte);
    real_type_1d_view tau_p("tau_p", mmd.nsalt);
    real_type_1d_view tau_d("tau_d", mmd.nsalt);

    std::string profile_name = "Verification_test_ASTEM_update_phase_eqblm";
    using policy_type =
          typename TChem::UseThisTeamPolicy<TChem::exec_space>::type;
    const auto exec_space_instance = TChem::exec_space();
    policy_type policy(exec_space_instance, 1, Kokkos::AUTO());

    Kokkos::parallel_for(
    profile_name,
    policy,
    KOKKOS_LAMBDA(const typename policy_type::member_type& member) {
      using MOSAIC = TChem::Impl::MOSAIC<real_type, device_type>;
      Kokkos::single(Kokkos::PerTeam(member), [&]() {
        MOSAIC::ASTEM_update_phase_eqblm(mmd, aer_liquid, aer_solid, aer_total,
          electrolyte_solid, electrolyte_liquid, electrolyte_total,
          epercent_liquid, epercent_solid, epercent_total,
          Keq_sl, Keq_ll, log_gamZ, MDRH_T, aH2O(0),
          jaerosolstate(0), jphase(0), jhyst_leg(0), aH2O_a(0), water_a(0),
          mass_dry_a(0), vol_dry_a(0), mass_wet_a(0), vol_wet_a(0), growth_factor(0),
          electrolyte_sum_solid(0), MDRH(0),
          jsalt_present, hsalt, phi_salt, phi_salt_old, phi_bar,
          alpha_salt, sat_ratio, flux_sl, frac_salt_solid, frac_salt_liq,
          eleliquid, store, na, nc, xeq_a, xeq_c, na_Ma, nc_Mc, aer_percent,
          molalities, xmol, ma, mc_view, log_gam, gam, activity, tau_p, tau_d);
      });
    });

    const auto write_view = [&](const std::string& name, const real_type_1d_view& v) {
      std::vector<real_type> out(v.extent(0));
      verification::convert_1d_view_device_to_1d_vector(v, out);
      output.set(name, out);
    };
    const auto write_view_2d = [&](const std::string& name, const real_type_2d_view& v) {
      std::vector<real_type> out(v.extent(0)*v.extent(1));
      verification::convert_2d_view_device_to_1d_vector(v, out);
      output.set(name, out);
    };

    write_view("aer_liquid", aer_liquid);
    write_view("aer_solid", aer_solid);
    write_view("aer_total", aer_total);
    write_view("electrolyte_solid", electrolyte_solid);
    write_view("electrolyte_liquid", electrolyte_liquid);
    write_view("electrolyte_total", electrolyte_total);
    write_view("epercent_liquid", epercent_liquid);
    write_view("epercent_solid", epercent_solid);
    write_view("epercent_total", epercent_total);
    write_view("Keq_ll", Keq_ll);
    write_view("Keq_sl", Keq_sl);
    write_view_2d("log_gamZ", log_gamZ);
    write_view("MDRH_T", MDRH_T);
    write_view("aH2O", aH2O);
    write_view("jaerosolstate", jaerosolstate);
    write_view("jphase", jphase);
    write_view("jhyst_leg", jhyst_leg);
    write_view("aH2O_a", aH2O_a);
    write_view("water_a", water_a);
    write_view("mass_dry_a", mass_dry_a);
    write_view("vol_dry_a", vol_dry_a);
    write_view("mass_wet_a", mass_wet_a);
    write_view("vol_wet_a", vol_wet_a);
    write_view("growth_factor", growth_factor);
    write_view("electrolyte_sum_solid", electrolyte_sum_solid);
    write_view("MDRH", MDRH);
    write_view("jsalt_present", jsalt_present);
    write_view("phi_salt_old", phi_salt_old);
    write_view("alpha_salt", alpha_salt);
    write_view("phi_salt", phi_salt);
    write_view("flux_sl", flux_sl);
    write_view("sat_ratio", sat_ratio);
    write_view("ma", ma);
    write_view("mc", mc_view);
    write_view("gam", gam);
    write_view("activity", activity);
  });
}
