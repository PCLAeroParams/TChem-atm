#include "TChem.hpp"
#include "TChem_Impl_MOSAIC.hpp"
#include <verification.hpp>
#include "skywalker.hpp"

using device_type = typename Tines::UseThisDevice<TChem::exec_space>::type;
using real_type_1d_view = TChem::real_type_1d_view;
using ordinal_type = TChem::ordinal_type;
using namespace skywalker;
using namespace TChem;

void form_nh4so4(Ensemble *ensemble) {
  ensemble->process([=](const Input &input, Output &output) {

    const auto store_arr = input.get_array("store");
    const auto electrolyte_arr = input.get_array("electrolyte");

    const auto mmd = TChem::Impl::MosaicModelData<device_type>();

    real_type_1d_view store("store", mmd.naer);
    verification::convert_1d_vector_to_1d_view_device(store_arr, store);

    real_type_1d_view electrolyte("electrolyte", mmd.nelectrolyte);
    verification::convert_1d_vector_to_1d_view_device(electrolyte_arr, electrolyte);

    std::string profile_name ="Verification_test_form_nh4so4";
    using policy_type =
          typename TChem::UseThisTeamPolicy<TChem::exec_space>::type;
    const auto exec_space_instance = TChem::exec_space();
    policy_type policy(exec_space_instance, 1, Kokkos::AUTO());

    Kokkos::parallel_for(
    profile_name,
    policy,
    KOKKOS_LAMBDA(const typename policy_type::member_type& member) {
    TChem::Impl::MOSAIC<real_type, device_type>::form_nh4so4(
      mmd,
      electrolyte,
      store);
    });

    std::vector<Real> store_out(mmd.naer), electrolyte_out(mmd.nelectrolyte);
    verification::convert_1d_view_device_to_1d_vector(store, store_out);
    verification::convert_1d_view_device_to_1d_vector(electrolyte, electrolyte_out);

    output.set("store", store_out);
    output.set("electrolyte", electrolyte_out);
  });
}
