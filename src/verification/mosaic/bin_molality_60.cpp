#include "TChem.hpp"
#include "TChem_Impl_MOSAIC.hpp"
#include <verification.hpp>
#include "skywalker.hpp"

using device_type = typename Tines::UseThisDevice<TChem::exec_space>::type;
using real_type_1d_view = TChem::real_type_1d_view;
using oridinal_type_1d_view = TChem::ordinal_type_1d_view;
using ordinal_type = TChem::ordinal_type;
using namespace skywalker;
using namespace TChem;

void bin_molality_60(Ensemble *ensemble) {
  ensemble->process([=](const Input &input, Output &output) {

    const auto je_arr = input.get_array("je");

    real_type_1d_view je("je", 1);
    verification::convert_1d_vector_to_1d_view_device(je_arr, je);

    const auto mmd = TChem::Impl::MosaicModelData<device_type>();

    // Prepare variables for output

    // Reals or int that are defined outside of the parallel_for region are passed as const.
    real_type_1d_view outputs_bin_molality_60("outputs_bin_molality_60", 1);

    std::string profile_name ="Verification_test_bin_molality_60";
    using policy_type =
          typename TChem::UseThisTeamPolicy<TChem::exec_space>::type;
    const auto exec_space_instance = TChem::exec_space();
    const auto host_exec_space = TChem::host_exec_space();
    policy_type policy(exec_space_instance, 1, Kokkos::AUTO());

    // Check this routines on GPUs.
     Kokkos::parallel_for(
    profile_name,
    policy,
    KOKKOS_LAMBDA(const typename policy_type::member_type& member) {

      Real& molality = outputs_bin_molality_60(0);

    // Perform the adjustment calculation
    TChem::Impl::MOSAIC<real_type, device_type>::bin_molality_60(
      mmd,
      je(0)-1,
      molality);
    });

    const auto outputs_bin_molality_60_h = Kokkos::create_mirror_view_and_copy(host_exec_space, outputs_bin_molality_60);

    Real molality = outputs_bin_molality_60_h(0);

    // Assuming the outputs are scalar and can be directly set in the ensemble
    output.set("molality", molality);

  });
}