#include "TChem.hpp"
#include "TChem_Impl_MOSAIC.hpp"
#include <verification.hpp>
#include "skywalker.hpp"

using device_type = typename Tines::UseThisDevice<TChem::exec_space>::type;
using real_type_1d_view = TChem::real_type_1d_view;
using ordinal_type = TChem::ordinal_type;
using namespace skywalker;
using namespace TChem;

void ASTEM_flux_wet_case3b(Ensemble *ensemble) {
  ensemble->process([=](const Input &input, Output &output) {

    const auto mmd = TChem::Impl::MosaicModelData<device_type>();

    const auto read_view = [&](const std::string& name, const ordinal_type n) {
      real_type_1d_view v(name, n);
      verification::convert_1d_vector_to_1d_view_device(input.get_array(name), v);
      return v;
    };

    const auto gas = read_view("gas", mmd.ngas_volatile);
    const auto kg = read_view("kg", mmd.ngas_volatile);
    const auto kel = read_view("kel", mmd.ngas_volatile);
    const auto Keq_ll = read_view("Keq_ll", mmd.nrxn_aer_ll);
    const auto Keq_gl = read_view("Keq_gl", mmd.nrxn_aer_gl);
    const auto Keq_nh4cl = read_view("Keq_nh4cl", 1);
    const auto gam = read_view("gam", mmd.nelectrolyte);
    const auto gam_ratio = read_view("gam_ratio", 1);
    const auto water_a = read_view("water_a", 1);
    const auto ma = read_view("ma", mmd.nanion);
    const auto mc = read_view("mc", mmd.ncation);
    const auto sfc_a = read_view("sfc_a", mmd.ngas_volatile);
    const auto df_gas_l = read_view("df_gas_l", mmd.ngas_volatile);
    const auto phi_volatile_l = read_view("phi_volatile_l", mmd.ngas_volatile);
    const auto Heff = read_view("Heff", mmd.ngas_volatile);
    const auto integrate_liquid = read_view("integrate_liquid", mmd.ngas_volatile);
    const auto ieqblm_ASTEM = read_view("ieqblm_ASTEM", 1);

    std::string profile_name = "Verification_test_ASTEM_flux_wet_case3b";
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
        MOSAIC::ASTEM_flux_wet_case3b(mmd, gas, kg, kel, Keq_ll, Keq_gl, Keq_nh4cl(0), gam,
          gam_ratio(0), water_a(0), ma, mc, sfc_a, df_gas_l, phi_volatile_l, Heff,
          integrate_liquid, ieqblm_ASTEM(0));
      });
    });

    const auto write_view = [&](const std::string& name, const real_type_1d_view& v) {
      std::vector<real_type> out(v.extent(0));
      verification::convert_1d_view_device_to_1d_vector(v, out);
      output.set(name, out);
    };

    write_view("gas", gas);
    write_view("kg", kg);
    write_view("kel", kel);
    write_view("Keq_ll", Keq_ll);
    write_view("Keq_gl", Keq_gl);
    write_view("Keq_nh4cl", Keq_nh4cl);
    write_view("gam", gam);
    write_view("gam_ratio", gam_ratio);
    write_view("water_a", water_a);
    write_view("ma", ma);
    write_view("mc", mc);
    write_view("sfc_a", sfc_a);
    write_view("df_gas_l", df_gas_l);
    write_view("phi_volatile_l", phi_volatile_l);
    write_view("Heff", Heff);
    write_view("integrate_liquid", integrate_liquid);
    write_view("ieqblm_ASTEM", ieqblm_ASTEM);
  });
}
