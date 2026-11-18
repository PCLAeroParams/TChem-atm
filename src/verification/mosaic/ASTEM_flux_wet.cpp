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

void ASTEM_flux_wet(Ensemble *ensemble) {
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

    const auto gas = read_view("gas", mmd.ngas_volatile);
    const auto kg = read_view("kg", mmd.ngas_volatile);
    const auto kel = read_view("kel", mmd.ngas_volatile);
    const auto Keq_ll = read_view("Keq_ll", mmd.nrxn_aer_ll);
    const auto Keq_gl = read_view("Keq_gl", mmd.nrxn_aer_gl);
    const auto log_gamZ = read_view_2d("log_gamZ", mmd.nelectrolyte, mmd.nelectrolyte);
    const auto Kp_nh4no3 = read_view("Kp_nh4no3", 1);
    const auto Kp_nh4cl = read_view("Kp_nh4cl", 1);
    const auto total_species = read_view("total_species", mmd.ngas_volatile);
    const auto tot_cl_in = read_view("tot_cl_in", 1);
    const auto delta_nh3_max = read_view("delta_nh3_max", 1);
    const auto delta_hno3_max = read_view("delta_hno3_max", 1);
    const auto delta_hcl_max = read_view("delta_hcl_max", 1);
    const auto aer_solid = read_view("aer_solid", mmd.naer);
    const auto aer_liquid = read_view("aer_liquid", mmd.naer);
    const auto aer_total = read_view("aer_total", mmd.naer);
    const auto electrolyte_solid = read_view("electrolyte_solid", mmd.nelectrolyte);
    const auto electrolyte_liquid = read_view("electrolyte_liquid", mmd.nelectrolyte);
    const auto electrolyte_total = read_view("electrolyte_total", mmd.nelectrolyte);
    const auto electrolyte_sum_liquid = read_view("electrolyte_sum_liquid", 1);
    const auto electrolyte_sum_total = read_view("electrolyte_sum_total", 1);
    const auto epercent_liquid = read_view("epercent_liquid", mmd.nelectrolyte);
    const auto aer_sum_liquid = read_view("aer_sum_liquid", 1);
    const auto aer_percent_liquid = read_view("aer_percent_liquid", mmd.naer);
    const auto ma = read_view("ma", mmd.nanion);
    const auto mc = read_view("mc", mmd.ncation);
    const auto gam = read_view("gam", mmd.nelectrolyte);
    const auto activity = read_view("activity", mmd.nelectrolyte);
    const auto gam_ratio = read_view("gam_ratio", 1);
    const auto jaerosolstate = read_view("jaerosolstate", 1);
    const auto jphase = read_view("jphase", 1);
    const auto jhyst_leg = read_view("jhyst_leg", 1);
    const auto aH2O_a = read_view("aH2O_a", 1);
    const auto water_a = read_view("water_a", 1);
    const auto Keq_nh4no3 = read_view("Keq_nh4no3", 1);
    const auto Keq_nh4cl = read_view("Keq_nh4cl", 1);
    const auto kelvin_nh4no3 = read_view("kelvin_nh4no3", 1);
    const auto kelvin_nh4cl = read_view("kelvin_nh4cl", 1);
    const auto sfc_a = read_view("sfc_a", mmd.ngas_volatile);
    const auto df_gas_s = read_view("df_gas_s", mmd.ngas_volatile);
    const auto df_gas_l = read_view("df_gas_l", mmd.ngas_volatile);
    const auto phi_volatile_s = read_view("phi_volatile_s", mmd.ngas_volatile);
    const auto phi_volatile_l = read_view("phi_volatile_l", mmd.ngas_volatile);
    const auto flux_s = read_view("flux_s", mmd.ngas_volatile);
    const auto Heff = read_view("Heff", mmd.ngas_volatile);
    const auto integrate_solid = read_view("integrate_solid", mmd.ngas_volatile);
    const auto integrate_liquid = read_view("integrate_liquid", mmd.ngas_volatile);
    const auto ieqblm_ASTEM = read_view("ieqblm_ASTEM", 1);

    // work arrays (local arrays in MOSAIC)
    real_type_1d_view store("store", mmd.naer);
    real_type_1d_view na("na", mmd.nanion);
    real_type_1d_view nc("nc", mmd.ncation);
    real_type_1d_view xeq_a("xeq_a", mmd.nanion);
    real_type_1d_view xeq_c("xeq_c", mmd.ncation);
    real_type_1d_view na_Ma("na_Ma", mmd.nanion);
    real_type_1d_view nc_Mc("nc_Mc", mmd.ncation);
    real_type_1d_view molalities("molalities", mmd.nelectrolyte);
    real_type_1d_view xmol("xmol", mmd.nelectrolyte);
    real_type_1d_view log_gam("log_gam", mmd.nelectrolyte);

    std::string profile_name = "Verification_test_ASTEM_flux_wet";
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
        MOSAIC::ASTEM_flux_wet(mmd, gas, kg, kel, Keq_ll, Keq_gl, log_gamZ, Kp_nh4no3(0),
          Kp_nh4cl(0), total_species, tot_cl_in(0), delta_nh3_max(0), delta_hno3_max(0),
          delta_hcl_max(0), aer_solid, aer_liquid, aer_total, electrolyte_solid,
          electrolyte_liquid, electrolyte_total, electrolyte_sum_liquid(0),
          electrolyte_sum_total(0), epercent_liquid, aer_sum_liquid(0), aer_percent_liquid,
          ma, mc, gam, activity, gam_ratio(0), jaerosolstate(0), jphase(0), jhyst_leg(0),
          aH2O_a(0), water_a(0), Keq_nh4no3(0), Keq_nh4cl(0), kelvin_nh4no3(0),
          kelvin_nh4cl(0), sfc_a, df_gas_s, df_gas_l, phi_volatile_s, phi_volatile_l,
          flux_s, Heff, integrate_solid, integrate_liquid, ieqblm_ASTEM(0), store, na, nc,
          xeq_a, xeq_c, na_Ma, nc_Mc, molalities, xmol, log_gam);
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

    write_view("gas", gas);
    write_view("kg", kg);
    write_view("kel", kel);
    write_view("Keq_ll", Keq_ll);
    write_view("Keq_gl", Keq_gl);
    write_view_2d("log_gamZ", log_gamZ);
    write_view("Kp_nh4no3", Kp_nh4no3);
    write_view("Kp_nh4cl", Kp_nh4cl);
    write_view("total_species", total_species);
    write_view("tot_cl_in", tot_cl_in);
    write_view("delta_nh3_max", delta_nh3_max);
    write_view("delta_hno3_max", delta_hno3_max);
    write_view("delta_hcl_max", delta_hcl_max);
    write_view("aer_solid", aer_solid);
    write_view("aer_liquid", aer_liquid);
    write_view("aer_total", aer_total);
    write_view("electrolyte_solid", electrolyte_solid);
    write_view("electrolyte_liquid", electrolyte_liquid);
    write_view("electrolyte_total", electrolyte_total);
    write_view("electrolyte_sum_liquid", electrolyte_sum_liquid);
    write_view("electrolyte_sum_total", electrolyte_sum_total);
    write_view("epercent_liquid", epercent_liquid);
    write_view("aer_sum_liquid", aer_sum_liquid);
    write_view("aer_percent_liquid", aer_percent_liquid);
    write_view("ma", ma);
    write_view("mc", mc);
    write_view("gam", gam);
    write_view("activity", activity);
    write_view("gam_ratio", gam_ratio);
    write_view("jaerosolstate", jaerosolstate);
    write_view("jphase", jphase);
    write_view("jhyst_leg", jhyst_leg);
    write_view("aH2O_a", aH2O_a);
    write_view("water_a", water_a);
    write_view("Keq_nh4no3", Keq_nh4no3);
    write_view("Keq_nh4cl", Keq_nh4cl);
    write_view("kelvin_nh4no3", kelvin_nh4no3);
    write_view("kelvin_nh4cl", kelvin_nh4cl);
    write_view("sfc_a", sfc_a);
    write_view("df_gas_s", df_gas_s);
    write_view("df_gas_l", df_gas_l);
    write_view("phi_volatile_s", phi_volatile_s);
    write_view("phi_volatile_l", phi_volatile_l);
    write_view("flux_s", flux_s);
    write_view("Heff", Heff);
    write_view("integrate_solid", integrate_solid);
    write_view("integrate_liquid", integrate_liquid);
    write_view("ieqblm_ASTEM", ieqblm_ASTEM);
  });
}
