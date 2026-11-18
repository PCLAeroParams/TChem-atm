#include "TChem.hpp"
#include "TChem_Impl_MOSAIC.hpp"
#include <verification.hpp>
#include "skywalker.hpp"

using device_type = typename Tines::UseThisDevice<TChem::exec_space>::type;
using real_type_1d_view = TChem::real_type_1d_view;
using ordinal_type = TChem::ordinal_type;
using namespace skywalker;
using namespace TChem;

void ASTEM_flux_soa(Ensemble *ensemble) {
  ensemble->process([=](const Input &input, Output &output) {

    const auto mmd = TChem::Impl::MosaicModelData<device_type>();

    const auto read_view = [&](const std::string& name, const ordinal_type n) {
      real_type_1d_view v(name, n);
      verification::convert_1d_vector_to_1d_view_device(input.get_array(name), v);
      return v;
    };

    const auto gas = read_view("gas", mmd.ngas_volatile);
    const auto aer_total = read_view("aer_total", mmd.naer);
    const auto sat_soa = read_view("sat_soa", mmd.ngas_volatile);
    const auto sfc_a = read_view("sfc_a", mmd.ngas_volatile);
    const auto df_gas_o = read_view("df_gas_o", mmd.ngas_volatile);
    const auto flux_o = read_view("flux_o", mmd.ngas_volatile);
    const auto phi_volatile_o = read_view("phi_volatile_o", mmd.ngas_volatile);
    const auto Heff = read_view("Heff", mmd.ngas_volatile);
    const auto integrate_total = read_view("integrate_total", mmd.ngas_volatile);
    const auto ieqblm_soa = read_view("ieqblm_soa", 1);

    std::string profile_name = "Verification_test_ASTEM_flux_soa";
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
        MOSAIC::ASTEM_flux_soa(mmd, gas, aer_total, sat_soa, sfc_a, df_gas_o, flux_o,
          phi_volatile_o, Heff, integrate_total, ieqblm_soa(0));
      });
    });

    const auto write_view = [&](const std::string& name, const real_type_1d_view& v) {
      std::vector<real_type> out(v.extent(0));
      verification::convert_1d_view_device_to_1d_vector(v, out);
      output.set(name, out);
    };

    write_view("gas", gas);
    write_view("aer_total", aer_total);
    write_view("sat_soa", sat_soa);
    write_view("sfc_a", sfc_a);
    write_view("df_gas_o", df_gas_o);
    write_view("flux_o", flux_o);
    write_view("phi_volatile_o", phi_volatile_o);
    write_view("Heff", Heff);
    write_view("integrate_total", integrate_total);
    write_view("ieqblm_soa", ieqblm_soa);
  });
}
