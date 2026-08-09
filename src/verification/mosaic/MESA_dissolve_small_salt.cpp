#include "TChem.hpp"
#include "TChem_Impl_MOSAIC.hpp"
#include <verification.hpp>
#include "skywalker.hpp"

using device_type = typename Tines::UseThisDevice<TChem::exec_space>::type;
using real_type_1d_view = TChem::real_type_1d_view;
using ordinal_type = TChem::ordinal_type;
using namespace skywalker;
using namespace TChem;

void MESA_dissolve_small_salt(Ensemble *ensemble) {
  ensemble->process([=](const Input &input, Output &output) {

    // js is Fortran 1-based; convert to 0-based for C++ views/mosaic constants
    // The Fortran yaml writes integers as 1-element arrays: "js: [1]"
    const auto js_arr = input.get_array("js");
    const ordinal_type js = static_cast<ordinal_type>(js_arr[0]) - 1;

    const auto aer_liquid_arr        = input.get_array("aer_liquid");
    const auto aer_solid_arr         = input.get_array("aer_solid");
    const auto electrolyte_solid_arr = input.get_array("electrolyte_solid");

    const auto mmd = TChem::Impl::MosaicModelData<device_type>();

    real_type_1d_view aer_liquid("aer_liquid", mmd.naer);
    real_type_1d_view aer_solid("aer_solid", mmd.naer);
    real_type_1d_view electrolyte_solid("electrolyte_solid", mmd.nelectrolyte);

    verification::convert_1d_vector_to_1d_view_device(aer_liquid_arr,        aer_liquid);
    verification::convert_1d_vector_to_1d_view_device(aer_solid_arr,         aer_solid);
    verification::convert_1d_vector_to_1d_view_device(electrolyte_solid_arr, electrolyte_solid);

    std::string profile_name = "Verification_test_MESA_dissolve_small_salt";
    using policy_type =
          typename TChem::UseThisTeamPolicy<TChem::exec_space>::type;
    const auto exec_space_instance = TChem::exec_space();
    const auto host_exec_space = TChem::host_exec_space();
    policy_type policy(exec_space_instance, 1, Kokkos::AUTO());

    Kokkos::parallel_for(
    profile_name,
    policy,
    KOKKOS_LAMBDA(const typename policy_type::member_type& member) {
      TChem::Impl::MOSAIC<real_type, device_type>::MESA_dissolve_small_salt(
        mmd,
        js,
        aer_liquid,
        aer_solid,
        electrolyte_solid);
    });

    std::vector<real_type> aer_liquid_out(mmd.naer);
    std::vector<real_type> aer_solid_out(mmd.naer);
    std::vector<real_type> electrolyte_solid_out(mmd.nelectrolyte);

    verification::convert_1d_view_device_to_1d_vector(aer_liquid,        aer_liquid_out);
    verification::convert_1d_view_device_to_1d_vector(aer_solid,         aer_solid_out);
    verification::convert_1d_view_device_to_1d_vector(electrolyte_solid, electrolyte_solid_out);

    // js is read-only but output per validation convention; use 1-element vector to match
    // Fortran write_output_var_int format: output.js=[[1],]
    output.set("js",                std::vector<real_type>{static_cast<real_type>(js + 1)});
    output.set("aer_liquid",        aer_liquid_out);
    output.set("aer_solid",         aer_solid_out);
    output.set("electrolyte_solid", electrolyte_solid_out);

  });
}
