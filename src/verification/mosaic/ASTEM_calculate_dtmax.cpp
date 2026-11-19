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

void ASTEM_calculate_dtmax(Ensemble *ensemble) {
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
    const auto kg = read_view_2d("kg", nbin_a, mmd.ngas_volatile);
    const auto flux_s = read_view_2d("flux_s", nbin_a, mmd.ngas_volatile);
    const auto df_gas_s = read_view_2d("df_gas_s", nbin_a, mmd.ngas_volatile);
    const auto phi_volatile_l = read_view_2d("phi_volatile_l", nbin_a, mmd.ngas_volatile);
    const auto integrate_liquid = read_view_2d("integrate_liquid", nbin_a, mmd.ngas_volatile);
    const auto electrolyte_solid = read_view_2d("electrolyte_solid", nbin_a, mmd.nelectrolyte);
    const auto aer_solid = read_view_2d("aer_solid", nbin_a, mmd.naer);
    const auto idry_case3a = read_view("idry_case3a", nbin_a);
    const auto jaerosolstate = read_view("jaerosolstate", nbin_a);
    const auto volatile_s = read_view_2d("volatile_s", nbin_a, mmd.ngas_volatile);
    const auto sum_bin_s = read_view("sum_bin_s", mmd.ngas_volatile);
    const auto sum_vdf_s = read_view("sum_vdf_s", mmd.ngas_volatile);
    const auto sum_vol_s = read_view("sum_vol_s", mmd.ngas_volatile);
    const auto avg_df_gas_s = read_view("avg_df_gas_s", mmd.ngas_volatile);
    const auto h_s_i_m = read_view_2d("h_s_i_m", nbin_a, mmd.ngas_volatile);
    // output only: MOSAIC sets dtmax before using it
    real_type_1d_view dtmax("dtmax", 1);

    std::string profile_name = "Verification_test_ASTEM_calculate_dtmax";
    using policy_type =
          typename TChem::UseThisTeamPolicy<TChem::exec_space>::type;
    const auto exec_space_instance = TChem::exec_space();
    policy_type policy(exec_space_instance, 1, Kokkos::AUTO());

    Kokkos::parallel_for(
    profile_name,
    policy,
    KOKKOS_LAMBDA(const typename policy_type::member_type& member) {
      using MOSAIC = TChem::Impl::MOSAIC<real_type, device_type>;
      // all threads of the team call the routine; each gets the same dtmax
      real_type dtmax_team = 0.0;
      MOSAIC::ASTEM_calculate_dtmax(member, mmd, nbin_a, dtchem(0), gas,
        kg, flux_s, df_gas_s, phi_volatile_l, integrate_liquid,
        electrolyte_solid, aer_solid, idry_case3a, jaerosolstate,
        volatile_s, sum_bin_s, sum_vdf_s, sum_vol_s, avg_df_gas_s,
        h_s_i_m, dtmax_team);
      Kokkos::single(Kokkos::PerTeam(member), [&]() {
        dtmax(0) = dtmax_team;
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
    write_view_2d("kg", kg);
    write_view_2d("flux_s", flux_s);
    write_view_2d("df_gas_s", df_gas_s);
    write_view_2d("phi_volatile_l", phi_volatile_l);
    write_view_2d("integrate_liquid", integrate_liquid);
    write_view_2d("electrolyte_solid", electrolyte_solid);
    write_view_2d("aer_solid", aer_solid);
    write_view("idry_case3a", idry_case3a);
    write_view("jaerosolstate", jaerosolstate);
    write_view_2d("volatile_s", volatile_s);
    write_view("sum_bin_s", sum_bin_s);
    write_view("sum_vdf_s", sum_vdf_s);
    write_view("sum_vol_s", sum_vol_s);
    write_view("avg_df_gas_s", avg_df_gas_s);
    write_view_2d("h_s_i_m", h_s_i_m);
    write_view("dtmax", dtmax);
  });
}
