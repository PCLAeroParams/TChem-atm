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

void ASTEM_semi_volatiles(Ensemble *ensemble) {
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

    // number of particles (MOSAIC bins); the per-particle arrays are (nbin_a, n)
    const ordinal_type nbin_a = static_cast<ordinal_type>(input.get_array("nbin_a")[0]);

    const auto dtchem = read_view("dtchem", 1);
    const auto aH2O = read_view("aH2O", 1);
    const auto Kp_nh4no3 = read_view("Kp_nh4no3", 1);
    const auto Kp_nh4cl = read_view("Kp_nh4cl", 1);
    const auto Kp_nh4no3_0 = read_view("Kp_nh4no3_0", 1);
    const auto Kp_nh4cl_0 = read_view("Kp_nh4cl_0", 1);
    const auto Keq_ll = read_view("Keq_ll", mmd.nrxn_aer_ll);
    const auto Keq_gl = read_view("Keq_gl", mmd.nrxn_aer_gl);
    const auto Keq_sg = read_view("Keq_sg", mmd.nrxn_aer_sg);
    const auto Keq_sl = read_view("Keq_sl", mmd.nrxn_aer_sl);
    const auto log_gamZ = read_view_2d("log_gamZ", mmd.nelectrolyte, mmd.nelectrolyte);
    const auto MDRH_T = read_view("MDRH_T", 63);
    const auto gas = read_view("gas", mmd.ngas_volatile);
    const auto total_species = read_view("total_species", mmd.ngas_volatile);
    const auto tot_cl_in = read_view("tot_cl_in", 1);
    const auto sum_bin_s = read_view("sum_bin_s", mmd.ngas_volatile);
    const auto sum_vdf_s = read_view("sum_vdf_s", mmd.ngas_volatile);
    const auto sum_vol_s = read_view("sum_vol_s", mmd.ngas_volatile);
    const auto avg_df_gas_s = read_view("avg_df_gas_s", mmd.ngas_volatile);
    const auto kg = read_view_2d("kg", nbin_a, mmd.ngas_volatile);
    const auto kel = read_view_2d("kel", nbin_a, mmd.ngas_volatile);
    const auto delta_nh3_max = read_view("delta_nh3_max", nbin_a);
    const auto delta_hno3_max = read_view("delta_hno3_max", nbin_a);
    const auto delta_hcl_max = read_view("delta_hcl_max", nbin_a);
    const auto aer_solid = read_view_2d("aer_solid", nbin_a, mmd.naer);
    const auto aer_liquid = read_view_2d("aer_liquid", nbin_a, mmd.naer);
    const auto aer_total = read_view_2d("aer_total", nbin_a, mmd.naer);
    const auto electrolyte_solid = read_view_2d("electrolyte_solid", nbin_a, mmd.nelectrolyte);
    const auto electrolyte_liquid = read_view_2d("electrolyte_liquid", nbin_a, mmd.nelectrolyte);
    const auto electrolyte_total = read_view_2d("electrolyte_total", nbin_a, mmd.nelectrolyte);
    const auto electrolyte_sum_solid = read_view("electrolyte_sum_solid", nbin_a);
    const auto electrolyte_sum_liquid = read_view("electrolyte_sum_liquid", nbin_a);
    const auto electrolyte_sum_total = read_view("electrolyte_sum_total", nbin_a);
    const auto epercent_solid = read_view_2d("epercent_solid", nbin_a, mmd.nelectrolyte);
    const auto epercent_liquid = read_view_2d("epercent_liquid", nbin_a, mmd.nelectrolyte);
    const auto epercent_total = read_view_2d("epercent_total", nbin_a, mmd.nelectrolyte);
    const auto ma = read_view_2d("ma", nbin_a, mmd.nanion);
    const auto mc = read_view_2d("mc", nbin_a, mmd.ncation);
    const auto gam = read_view_2d("gam", nbin_a, mmd.nelectrolyte);
    const auto activity = read_view_2d("activity", nbin_a, mmd.nelectrolyte);
    const auto gam_ratio = read_view("gam_ratio", nbin_a);
    const auto jaerosolstate = read_view("jaerosolstate", nbin_a);
    const auto jphase = read_view("jphase", nbin_a);
    const auto jhyst_leg = read_view("jhyst_leg", nbin_a);
    const auto aH2O_a = read_view("aH2O_a", nbin_a);
    const auto water_a = read_view("water_a", nbin_a);
    const auto mass_dry_a = read_view("mass_dry_a", nbin_a);
    const auto vol_dry_a = read_view("vol_dry_a", nbin_a);
    const auto mass_wet_a = read_view("mass_wet_a", nbin_a);
    const auto vol_wet_a = read_view("vol_wet_a", nbin_a);
    const auto growth_factor = read_view("growth_factor", nbin_a);
    const auto MDRH = read_view("MDRH", nbin_a);
    const auto df_gas_s = read_view_2d("df_gas_s", nbin_a, mmd.ngas_volatile);
    const auto df_gas_l = read_view_2d("df_gas_l", nbin_a, mmd.ngas_volatile);
    const auto phi_volatile_s = read_view_2d("phi_volatile_s", nbin_a, mmd.ngas_volatile);
    const auto phi_volatile_l = read_view_2d("phi_volatile_l", nbin_a, mmd.ngas_volatile);
    const auto flux_s = read_view_2d("flux_s", nbin_a, mmd.ngas_volatile);
    const auto flux_l = read_view_2d("flux_l", nbin_a, mmd.ngas_volatile);
    const auto Heff = read_view_2d("Heff", nbin_a, mmd.ngas_volatile);
    const auto integrate_solid = read_view_2d("integrate_solid", nbin_a, mmd.ngas_volatile);
    const auto integrate_liquid = read_view_2d("integrate_liquid", nbin_a, mmd.ngas_volatile);
    const auto volatile_s = read_view_2d("volatile_s", nbin_a, mmd.ngas_volatile);
    const auto h_s_i_m = read_view_2d("h_s_i_m", nbin_a, mmd.ngas_volatile);
    const auto idry_case3a = read_view("idry_case3a", nbin_a);
    const auto cumul_steps_ASTEM = read_view("cumul_steps_ASTEM", 1);
    const auto isteps_ASTEM_max = read_view("isteps_ASTEM_max", 1);
    const auto jASTEM_fail = read_view("jASTEM_fail", 1);

    // outputs only: the routine resets them before using them
    real_type_1d_view ieqblm_ASTEM("ieqblm_ASTEM", 1);
    real_type_1d_view isteps_ASTEM("isteps_ASTEM", 1);

    // work arrays (local or recomputed before use in MOSAIC), one row per particle
    real_type_2d_view sfc_a("sfc_a", nbin_a, mmd.ngas_volatile);
    real_type_2d_view store("store", nbin_a, mmd.naer);
    real_type_2d_view aer_percent("aer_percent", nbin_a, mmd.naer);
    real_type_2d_view na("na", nbin_a, mmd.nanion);
    real_type_2d_view nc("nc", nbin_a, mmd.ncation);
    real_type_2d_view xeq_a("xeq_a", nbin_a, mmd.nanion);
    real_type_2d_view xeq_c("xeq_c", nbin_a, mmd.ncation);
    real_type_2d_view na_Ma("na_Ma", nbin_a, mmd.nanion);
    real_type_2d_view nc_Mc("nc_Mc", nbin_a, mmd.ncation);
    real_type_2d_view molalities("molalities", nbin_a, mmd.nelectrolyte);
    real_type_2d_view xmol("xmol", nbin_a, mmd.nelectrolyte);
    real_type_2d_view log_gam("log_gam", nbin_a, mmd.nelectrolyte);
    real_type_2d_view jsalt_present("jsalt_present", nbin_a, mmd.nsalt);
    real_type_2d_view hsalt("hsalt", nbin_a, mmd.nsalt);
    real_type_2d_view phi_salt("phi_salt", nbin_a, mmd.nsalt);
    real_type_2d_view phi_salt_old("phi_salt_old", nbin_a, mmd.nsalt);
    real_type_2d_view phi_bar("phi_bar", nbin_a, mmd.nsalt);
    real_type_2d_view alpha_salt("alpha_salt", nbin_a, mmd.nsalt);
    real_type_2d_view sat_ratio("sat_ratio", nbin_a, mmd.nsalt);
    real_type_2d_view flux_sl("flux_sl", nbin_a, mmd.nsalt);
    real_type_2d_view frac_salt_solid("frac_salt_solid", nbin_a, mmd.nsalt);
    real_type_2d_view frac_salt_liq("frac_salt_liq", nbin_a, mmd.nsalt);
    real_type_2d_view eleliquid("eleliquid", nbin_a, mmd.nelectrolyte);
    real_type_2d_view tau_p("tau_p", nbin_a, mmd.nsalt);
    real_type_2d_view tau_d("tau_d", nbin_a, mmd.nsalt);
    real_type_2d_view gas_bin("gas_bin", nbin_a, mmd.ngas_volatile);
    real_type_2d_view total_species_bin("total_species_bin", nbin_a, mmd.ngas_volatile);
    real_type_2d_view hcl_net("hcl_net", nbin_a, 3);
    real_type_1d_view aer_sum("aer_sum", nbin_a);

    std::string profile_name = "Verification_test_ASTEM_semi_volatiles";
    using policy_type =
          typename TChem::UseThisTeamPolicy<TChem::exec_space>::type;
    const auto exec_space_instance = TChem::exec_space();
    policy_type policy(exec_space_instance, 1, Kokkos::AUTO());

    Kokkos::parallel_for(
    profile_name,
    policy,
    KOKKOS_LAMBDA(const typename policy_type::member_type& member) {
      using MOSAIC = TChem::Impl::MOSAIC<real_type, device_type>;
      // all threads of the team call the routine; each gets the same values of these
      real_type cumul_steps_ASTEM_team = cumul_steps_ASTEM(0);
      real_type isteps_ASTEM_max_team = isteps_ASTEM_max(0);
      real_type jASTEM_fail_team = jASTEM_fail(0);
      real_type ieqblm_ASTEM_team = 0.0;
      real_type isteps_ASTEM_team = 0.0;
      MOSAIC::ASTEM_semi_volatiles(member, mmd, nbin_a, dtchem(0), aH2O(0), Kp_nh4no3(0),
        Kp_nh4cl(0), Kp_nh4no3_0(0), Kp_nh4cl_0(0), Keq_ll, Keq_gl, Keq_sg, Keq_sl, log_gamZ,
        MDRH_T, gas, total_species, tot_cl_in(0), sum_bin_s, sum_vdf_s, sum_vol_s, avg_df_gas_s,
        kg, kel, delta_nh3_max, delta_hno3_max, delta_hcl_max, aer_solid, aer_liquid, aer_total,
        electrolyte_solid, electrolyte_liquid, electrolyte_total, electrolyte_sum_solid,
        electrolyte_sum_liquid, electrolyte_sum_total, epercent_solid, epercent_liquid,
        epercent_total, ma, mc, gam, activity, gam_ratio, jaerosolstate, jphase, jhyst_leg, aH2O_a,
        water_a, mass_dry_a, vol_dry_a, mass_wet_a, vol_wet_a, growth_factor, MDRH, df_gas_s,
        df_gas_l, phi_volatile_s, phi_volatile_l, flux_s, flux_l, Heff, integrate_solid,
        integrate_liquid, volatile_s, h_s_i_m, idry_case3a, sfc_a, store, aer_percent, na, nc,
        xeq_a, xeq_c, na_Ma, nc_Mc, molalities, xmol, log_gam, jsalt_present, hsalt, phi_salt,
        phi_salt_old, phi_bar, alpha_salt, sat_ratio, flux_sl, frac_salt_solid, frac_salt_liq,
        eleliquid, tau_p, tau_d, gas_bin, total_species_bin, hcl_net, aer_sum,
        cumul_steps_ASTEM_team, isteps_ASTEM_max_team, jASTEM_fail_team, ieqblm_ASTEM_team,
        isteps_ASTEM_team);
      Kokkos::single(Kokkos::PerTeam(member), [&]() {
        cumul_steps_ASTEM(0) = cumul_steps_ASTEM_team;
        isteps_ASTEM_max(0) = isteps_ASTEM_max_team;
        jASTEM_fail(0) = jASTEM_fail_team;
        ieqblm_ASTEM(0) = ieqblm_ASTEM_team;
        isteps_ASTEM(0) = isteps_ASTEM_team;
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

    output.set("nbin_a", std::vector<real_type>{static_cast<real_type>(nbin_a)});
    // gam_ratio and electrolyte_sum_liquid are not compared. MESA_PTC updates them in MOSAIC but not here,
    // and MOSAIC does not use those values: the flux routines recompute gam_ratio before reading it,
    // and electrolyte_sum(jliquid) is never read.
    write_view("dtchem", dtchem);
    write_view("aH2O", aH2O);
    write_view("Kp_nh4no3", Kp_nh4no3);
    write_view("Kp_nh4cl", Kp_nh4cl);
    write_view("Kp_nh4no3_0", Kp_nh4no3_0);
    write_view("Kp_nh4cl_0", Kp_nh4cl_0);
    write_view("Keq_ll", Keq_ll);
    write_view("Keq_gl", Keq_gl);
    write_view("Keq_sg", Keq_sg);
    write_view("Keq_sl", Keq_sl);
    write_view_2d("log_gamZ", log_gamZ);
    write_view("MDRH_T", MDRH_T);
    write_view("gas", gas);
    write_view("total_species", total_species);
    write_view("tot_cl_in", tot_cl_in);
    write_view("sum_bin_s", sum_bin_s);
    write_view("sum_vdf_s", sum_vdf_s);
    write_view("sum_vol_s", sum_vol_s);
    write_view("avg_df_gas_s", avg_df_gas_s);
    write_view_2d("kg", kg);
    write_view_2d("kel", kel);
    write_view("delta_nh3_max", delta_nh3_max);
    write_view("delta_hno3_max", delta_hno3_max);
    write_view("delta_hcl_max", delta_hcl_max);
    write_view_2d("aer_solid", aer_solid);
    write_view_2d("aer_liquid", aer_liquid);
    write_view_2d("aer_total", aer_total);
    write_view_2d("electrolyte_solid", electrolyte_solid);
    write_view_2d("electrolyte_liquid", electrolyte_liquid);
    write_view_2d("electrolyte_total", electrolyte_total);
    write_view("electrolyte_sum_solid", electrolyte_sum_solid);
    write_view("electrolyte_sum_total", electrolyte_sum_total);
    write_view_2d("epercent_solid", epercent_solid);
    write_view_2d("epercent_liquid", epercent_liquid);
    write_view_2d("epercent_total", epercent_total);
    write_view_2d("ma", ma);
    write_view_2d("mc", mc);
    write_view_2d("gam", gam);
    write_view_2d("activity", activity);
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
    write_view("MDRH", MDRH);
    write_view_2d("df_gas_s", df_gas_s);
    write_view_2d("df_gas_l", df_gas_l);
    write_view_2d("phi_volatile_s", phi_volatile_s);
    write_view_2d("phi_volatile_l", phi_volatile_l);
    write_view_2d("flux_s", flux_s);
    write_view_2d("flux_l", flux_l);
    write_view_2d("Heff", Heff);
    write_view_2d("integrate_solid", integrate_solid);
    write_view_2d("integrate_liquid", integrate_liquid);
    write_view_2d("volatile_s", volatile_s);
    write_view_2d("h_s_i_m", h_s_i_m);
    write_view("idry_case3a", idry_case3a);
    write_view("cumul_steps_ASTEM", cumul_steps_ASTEM);
    write_view("isteps_ASTEM_max", isteps_ASTEM_max);
    write_view("jASTEM_fail", jASTEM_fail);
    write_view("ieqblm_ASTEM", ieqblm_ASTEM);
    write_view("isteps_ASTEM", isteps_ASTEM);
  });
}
