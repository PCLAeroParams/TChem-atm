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

void ASTEM_secondary_organics(Ensemble *ensemble) {
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
    const auto gas = read_view("gas", mmd.ngas_volatile);
    const auto total_species = read_view("total_species", mmd.ngas_volatile);
    const auto sat_soa = read_view("sat_soa", mmd.ngas_volatile);
    const auto aer_total = read_view_2d("aer_total", nbin_a, mmd.naer);
    const auto kg = read_view_2d("kg", nbin_a, mmd.ngas_volatile);
    const auto jaerosolstate = read_view("jaerosolstate", nbin_a);
    const auto df_gas_o = read_view_2d("df_gas_o", nbin_a, mmd.ngas_volatile);
    const auto flux_o = read_view_2d("flux_o", nbin_a, mmd.ngas_volatile);
    const auto phi_volatile_o = read_view_2d("phi_volatile_o", nbin_a, mmd.ngas_volatile);
    const auto phi_volatile_l = read_view_2d("phi_volatile_l", nbin_a, mmd.ngas_volatile);
    const auto Heff = read_view_2d("Heff", nbin_a, mmd.ngas_volatile);
    const auto integrate_total = read_view_2d("integrate_total", nbin_a, mmd.ngas_volatile);
    // outputs only: the routine resets them before using them
    real_type_1d_view ieqblm_soa("ieqblm_soa", 1);
    real_type_1d_view isteps_SOA("isteps_SOA", 1);

    // work array (recomputed before use in MOSAIC), one row per particle
    real_type_2d_view sfc_a("sfc_a", nbin_a, mmd.ngas_volatile);

    std::string profile_name = "Verification_test_ASTEM_secondary_organics";
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
      real_type ieqblm_soa_team = 0.0;
      real_type isteps_SOA_team = 0.0;
      MOSAIC::ASTEM_secondary_organics(member, mmd, nbin_a, dtchem(0), gas, total_species,
        sat_soa, aer_total, kg, jaerosolstate, df_gas_o, flux_o, phi_volatile_o,
        phi_volatile_l, Heff, integrate_total, sfc_a, ieqblm_soa_team, isteps_SOA_team);
      Kokkos::single(Kokkos::PerTeam(member), [&]() {
        ieqblm_soa(0) = ieqblm_soa_team;
        isteps_SOA(0) = isteps_SOA_team;
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
    write_view("dtchem", dtchem);
    write_view("gas", gas);
    write_view("total_species", total_species);
    write_view("sat_soa", sat_soa);
    write_view_2d("aer_total", aer_total);
    write_view_2d("kg", kg);
    write_view("jaerosolstate", jaerosolstate);
    write_view_2d("df_gas_o", df_gas_o);
    write_view_2d("flux_o", flux_o);
    write_view_2d("phi_volatile_o", phi_volatile_o);
    write_view_2d("phi_volatile_l", phi_volatile_l);
    write_view_2d("Heff", Heff);
    write_view_2d("integrate_total", integrate_total);
    write_view("ieqblm_soa", ieqblm_soa);
    write_view("isteps_SOA", isteps_SOA);
  });
}
