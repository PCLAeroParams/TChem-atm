#include "TChem.hpp"
#include "TChem_Impl_MOSAIC.hpp"
#include <verification.hpp>
#include "skywalker.hpp"

using device_type = typename Tines::UseThisDevice<TChem::exec_space>::type;
using real_type_1d_view = TChem::real_type_1d_view;
using ordinal_type = TChem::ordinal_type;
using namespace skywalker;
using namespace TChem;

void ASTEM_flux_dry(Ensemble *ensemble) {
  ensemble->process([=](const Input &input, Output &output) {

    const auto mmd = TChem::Impl::MosaicModelData<device_type>();

    const auto read_view = [&](const std::string& name, const ordinal_type n) {
      real_type_1d_view v(name, n);
      verification::convert_1d_vector_to_1d_view_device(input.get_array(name), v);
      return v;
    };

    const auto gas = read_view("gas", mmd.ngas_volatile);
    const auto kg = read_view("kg", mmd.ngas_volatile);
    const auto Keq_sg = read_view("Keq_sg", mmd.nrxn_aer_sg);
    const auto Kp_nh4no3_0 = read_view("Kp_nh4no3_0", 1);
    const auto Kp_nh4cl_0 = read_view("Kp_nh4cl_0", 1);
    const auto electrolyte_solid = read_view("electrolyte_solid", mmd.nelectrolyte);
    const auto epercent_solid = read_view("epercent_solid", mmd.nelectrolyte);
    const auto aer_solid = read_view("aer_solid", mmd.naer);
    const auto aer_liquid = read_view("aer_liquid", mmd.naer);
    const auto aer_total = read_view("aer_total", mmd.naer);
    const auto Keq_nh4no3_0 = read_view("Keq_nh4no3_0", 1);
    const auto Keq_nh4cl_0 = read_view("Keq_nh4cl_0", 1);
    const auto phi_nh4no3_s = read_view("phi_nh4no3_s", 1);
    const auto phi_nh4cl_s = read_view("phi_nh4cl_s", 1);
    const auto sfc_a = read_view("sfc_a", mmd.ngas_volatile);
    const auto df_gas_s = read_view("df_gas_s", mmd.ngas_volatile);
    const auto phi_volatile_s = read_view("phi_volatile_s", mmd.ngas_volatile);
    const auto flux_s = read_view("flux_s", mmd.ngas_volatile);
    const auto integrate_solid = read_view("integrate_solid", mmd.ngas_volatile);
    const auto idry_case3a = read_view("idry_case3a", 1);
    const auto ieqblm_ASTEM = read_view("ieqblm_ASTEM", 1);

    std::string profile_name = "Verification_test_ASTEM_flux_dry";
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
        MOSAIC::ASTEM_flux_dry(mmd, gas, kg, Keq_sg, Kp_nh4no3_0(0), Kp_nh4cl_0(0),
          electrolyte_solid, epercent_solid, aer_solid, aer_liquid, aer_total,
          Keq_nh4no3_0(0), Keq_nh4cl_0(0), phi_nh4no3_s(0), phi_nh4cl_s(0), sfc_a, df_gas_s,
          phi_volatile_s, flux_s, integrate_solid, idry_case3a(0), ieqblm_ASTEM(0));
      });
    });

    const auto write_view = [&](const std::string& name, const real_type_1d_view& v) {
      std::vector<real_type> out(v.extent(0));
      verification::convert_1d_view_device_to_1d_vector(v, out);
      output.set(name, out);
    };

    write_view("gas", gas);
    write_view("kg", kg);
    write_view("Keq_sg", Keq_sg);
    write_view("Kp_nh4no3_0", Kp_nh4no3_0);
    write_view("Kp_nh4cl_0", Kp_nh4cl_0);
    write_view("electrolyte_solid", electrolyte_solid);
    write_view("epercent_solid", epercent_solid);
    write_view("aer_solid", aer_solid);
    write_view("aer_liquid", aer_liquid);
    write_view("aer_total", aer_total);
    write_view("Keq_nh4no3_0", Keq_nh4no3_0);
    write_view("Keq_nh4cl_0", Keq_nh4cl_0);
    write_view("phi_nh4no3_s", phi_nh4no3_s);
    write_view("phi_nh4cl_s", phi_nh4cl_s);
    write_view("sfc_a", sfc_a);
    write_view("df_gas_s", df_gas_s);
    write_view("phi_volatile_s", phi_volatile_s);
    write_view("flux_s", flux_s);
    write_view("integrate_solid", integrate_solid);
    write_view("idry_case3a", idry_case3a);
    write_view("ieqblm_ASTEM", ieqblm_ASTEM);
  });
}
