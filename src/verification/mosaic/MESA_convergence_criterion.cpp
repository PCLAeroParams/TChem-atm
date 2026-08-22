#include "TChem.hpp"
#include "TChem_Impl_MOSAIC.hpp"
#include <verification.hpp>
#include "skywalker.hpp"

using device_type = typename Tines::UseThisDevice<TChem::exec_space>::type;
using real_type_1d_view = TChem::real_type_1d_view;
using namespace skywalker;
using namespace TChem;

void MESA_convergence_criterion(Ensemble *ensemble) {
  ensemble->process([=](const Input &input, Output &output) {

    const auto electrolyte_solid_arr = input.get_array("electrolyte_solid");
    const auto mass_dry_salt_arr     = input.get_array("mass_dry_salt");
    const auto phi_salt_arr          = input.get_array("phi_salt");
    const auto flux_sl_arr           = input.get_array("flux_sl");
    const auto aer_solid_arr         = input.get_array("aer_solid");

    const auto mmd = TChem::Impl::MosaicModelData<device_type>();

    real_type_1d_view electrolyte_solid("electrolyte_solid", mmd.nelectrolyte);
    verification::convert_1d_vector_to_1d_view_device(electrolyte_solid_arr, electrolyte_solid);

    real_type_1d_view mass_dry_salt_view("mass_dry_salt", 1);
    verification::convert_1d_vector_to_1d_view_device(mass_dry_salt_arr, mass_dry_salt_view);

    real_type_1d_view phi_salt("phi_salt", mmd.nsalt);
    verification::convert_1d_vector_to_1d_view_device(phi_salt_arr, phi_salt);

    real_type_1d_view flux_sl("flux_sl", mmd.nsalt);
    verification::convert_1d_vector_to_1d_view_device(flux_sl_arr, flux_sl);

    real_type_1d_view aer_solid("aer_solid", mmd.naer);
    verification::convert_1d_vector_to_1d_view_device(aer_solid_arr, aer_solid);

    real_type_1d_view iconverge_mass_view("iconverge_mass", 1);
    real_type_1d_view iconverge_flux_view("iconverge_flux", 1);
    real_type_1d_view idissolved_view("idissolved", 1);

    std::string profile_name = "Verification_test_MESA_convergence_criterion";
    using policy_type =
          typename TChem::UseThisTeamPolicy<TChem::exec_space>::type;
    const auto exec_space_instance = TChem::exec_space();
    policy_type policy(exec_space_instance, 1, Kokkos::AUTO());

    Kokkos::parallel_for(
    profile_name,
    policy,
    KOKKOS_LAMBDA(const typename policy_type::member_type& member) {
      real_type iconverge_mass = 0;
      real_type iconverge_flux = 0;
      real_type idissolved     = 0;

      TChem::Impl::MOSAIC<real_type, device_type>::MESA_convergence_criterion(
        mmd,
        electrolyte_solid,
        mass_dry_salt_view(0),
        phi_salt,
        flux_sl,
        aer_solid,
        iconverge_mass,
        iconverge_flux,
        idissolved);

      iconverge_mass_view(0) = static_cast<real_type>(iconverge_mass);
      iconverge_flux_view(0) = static_cast<real_type>(iconverge_flux);
      idissolved_view(0)     = static_cast<real_type>(idissolved);
    });

    std::vector<real_type> electrolyte_solid_out(mmd.nelectrolyte);
    std::vector<real_type> mass_dry_salt_out(1);
    std::vector<real_type> phi_salt_out(mmd.nsalt);
    std::vector<real_type> flux_sl_out(mmd.nsalt);
    std::vector<real_type> aer_solid_out(mmd.naer);
    std::vector<real_type> iconverge_mass_out(1);
    std::vector<real_type> iconverge_flux_out(1);
    std::vector<real_type> idissolved_out(1);

    verification::convert_1d_view_device_to_1d_vector(electrolyte_solid, electrolyte_solid_out);
    verification::convert_1d_view_device_to_1d_vector(mass_dry_salt_view, mass_dry_salt_out);
    verification::convert_1d_view_device_to_1d_vector(phi_salt, phi_salt_out);
    verification::convert_1d_view_device_to_1d_vector(flux_sl, flux_sl_out);
    verification::convert_1d_view_device_to_1d_vector(aer_solid, aer_solid_out);
    verification::convert_1d_view_device_to_1d_vector(iconverge_mass_view, iconverge_mass_out);
    verification::convert_1d_view_device_to_1d_vector(iconverge_flux_view, iconverge_flux_out);
    verification::convert_1d_view_device_to_1d_vector(idissolved_view, idissolved_out);

    output.set("electrolyte_solid", electrolyte_solid_out);
    output.set("mass_dry_salt",     mass_dry_salt_out);
    output.set("phi_salt",          phi_salt_out);
    output.set("flux_sl",           flux_sl_out);
    output.set("aer_solid",         aer_solid_out);
    output.set("iconverge_mass",    iconverge_mass_out);
    output.set("iconverge_flux",    iconverge_flux_out);
    output.set("idissolved",        idissolved_out);
  });
}
