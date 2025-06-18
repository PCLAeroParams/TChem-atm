#include "TChem.hpp"
#include "TChem_Impl_MOSAIC.hpp"
#include <verification.hpp>
#include "skywalker.hpp"

using real_type = TChem::real_type;
using value_type = real_type;
using device_type = typename Tines::UseThisDevice<TChem::exec_space>::type;
using real_type_1d_view = TChem::real_type_1d_view;
using real_type_2d_view = TChem::real_type_2d_view;
using ordinal_type = TChem::ordinal_type;
using namespace skywalker;
using namespace TChem;

void MTEM_compute_log_gamZ(Ensemble *ensemble) {
  ensemble->process([=](const Input &input, Output &output) {
    // Assuming the input arrays are 1D and of the same length

    const auto aH2O_arr = input.get_array("aH2O");

    const auto mmd = TChem::Impl::MosaicModelData<device_type>();

    real_type_2d_view log_gamZ("log_gamZ", mmd.nelectrolyte, mmd.nelectrolyte);

    real_type_1d_view aH2O("aH20", 1);
    verification::convert_1d_vector_to_1d_view_device(aH2O_arr, aH2O);

    std::string profile_name ="Verification_test_MTEM_compute_log_gamZ";
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

    // Perform the adjustment calculation
    TChem::Impl::MOSAIC<real_type, device_type>::MTEM_compute_log_gamZ(
      mmd,
      aH2O(0),
      log_gamZ);
    });

    std::vector<Real> log_gamZ_out(mmd.nelectrolyte*mmd.nelectrolyte, 0.0);
    verification::convert_2d_view_device_to_1d_vector(log_gamZ, log_gamZ_out);

    output.set("log_gamZ", log_gamZ_out);
  });
}
