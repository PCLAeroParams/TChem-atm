#include "TChem.hpp"
#include "TChem_Impl_MOSAIC.hpp"
#include <verification.hpp>
#include "skywalker.hpp"

using device_type = typename Tines::UseThisDevice<TChem::exec_space>::type;
using real_type_1d_view = TChem::real_type_1d_view;
using ordinal_type = TChem::ordinal_type;
using namespace skywalker;
using namespace TChem;

void MESA_estimate_eleliquid(Ensemble *ensemble) {
  ensemble->process([=](const Input &input, Output &output) {

    const auto aer_liquid_arr             = input.get_array("aer_liquid");
    const auto eleliquid_arr              = input.get_array("eleliquid");
    const auto na_arr                     = input.get_array("na");
    const auto nc_arr                     = input.get_array("nc");
    const auto xeq_a_arr                  = input.get_array("xeq_a");
    const auto xeq_c_arr                  = input.get_array("xeq_c");
    const auto na_Ma_arr                  = input.get_array("na_Ma");
    const auto nc_Mc_arr                  = input.get_array("nc_Mc");
    const auto store_arr                  = input.get_array("store");
    const auto electrolyte_sum_liquid_arr = input.get_array("electrolyte_sum_liquid");
    const auto epercent_liquid_arr        = input.get_array("epercent_liquid");

    const auto mmd = TChem::Impl::MosaicModelData<device_type>();

    real_type_1d_view aer_liquid("aer_liquid", mmd.naer);
    verification::convert_1d_vector_to_1d_view_device(aer_liquid_arr, aer_liquid);

    real_type_1d_view eleliquid("eleliquid", mmd.nelectrolyte);
    verification::convert_1d_vector_to_1d_view_device(eleliquid_arr, eleliquid);

    real_type_1d_view na("na", mmd.nanion);
    verification::convert_1d_vector_to_1d_view_device(na_arr, na);

    real_type_1d_view nc("nc", mmd.ncation);
    verification::convert_1d_vector_to_1d_view_device(nc_arr, nc);

    real_type_1d_view xeq_a("xeq_a", mmd.nanion);
    verification::convert_1d_vector_to_1d_view_device(xeq_a_arr, xeq_a);

    real_type_1d_view xeq_c("xeq_c", mmd.ncation);
    verification::convert_1d_vector_to_1d_view_device(xeq_c_arr, xeq_c);

    real_type_1d_view na_Ma("na_Ma", mmd.nanion);
    verification::convert_1d_vector_to_1d_view_device(na_Ma_arr, na_Ma);

    real_type_1d_view nc_Mc("nc_Mc", mmd.ncation);
    verification::convert_1d_vector_to_1d_view_device(nc_Mc_arr, nc_Mc);

    real_type_1d_view store("store", mmd.naer);
    verification::convert_1d_vector_to_1d_view_device(store_arr, store);

    real_type_1d_view electrolyte_sum_liquid("electrolyte_sum_liquid", 1);
    verification::convert_1d_vector_to_1d_view_device(electrolyte_sum_liquid_arr, electrolyte_sum_liquid);

    real_type_1d_view epercent_liquid("epercent_liquid", mmd.nelectrolyte);
    verification::convert_1d_vector_to_1d_view_device(epercent_liquid_arr, epercent_liquid);

    std::string profile_name = "Verification_test_MESA_estimate_eleliquid";
    using policy_type =
          typename TChem::UseThisTeamPolicy<TChem::exec_space>::type;
    const auto exec_space_instance = TChem::exec_space();
    const auto host_exec_space = TChem::host_exec_space();
    policy_type policy(exec_space_instance, 1, Kokkos::AUTO());

    Kokkos::parallel_for(
    profile_name,
    policy,
    KOKKOS_LAMBDA(const typename policy_type::member_type& member) {
      TChem::Impl::MOSAIC<real_type, device_type>::MESA_estimate_eleliquid(
        mmd,
        aer_liquid,
        eleliquid,
        electrolyte_sum_liquid(0),
        epercent_liquid,
        na,
        nc,
        xeq_a,
        xeq_c,
        na_Ma,
        nc_Mc,
        store);
    });

    std::vector<real_type> aer_liquid_out(mmd.naer);
    std::vector<real_type> eleliquid_out(mmd.nelectrolyte);
    std::vector<real_type> electrolyte_sum_liquid_out(1);
    std::vector<real_type> epercent_liquid_out(mmd.nelectrolyte);
    std::vector<real_type> na_out(mmd.nanion);
    std::vector<real_type> nc_out(mmd.ncation);
    std::vector<real_type> xeq_a_out(mmd.nanion);
    std::vector<real_type> xeq_c_out(mmd.ncation);
    std::vector<real_type> na_Ma_out(mmd.nanion);
    std::vector<real_type> nc_Mc_out(mmd.ncation);
    std::vector<real_type> store_out(mmd.naer);

    verification::convert_1d_view_device_to_1d_vector(aer_liquid, aer_liquid_out);
    verification::convert_1d_view_device_to_1d_vector(eleliquid, eleliquid_out);
    verification::convert_1d_view_device_to_1d_vector(electrolyte_sum_liquid, electrolyte_sum_liquid_out);
    verification::convert_1d_view_device_to_1d_vector(epercent_liquid, epercent_liquid_out);
    verification::convert_1d_view_device_to_1d_vector(na, na_out);
    verification::convert_1d_view_device_to_1d_vector(nc, nc_out);
    verification::convert_1d_view_device_to_1d_vector(xeq_a, xeq_a_out);
    verification::convert_1d_view_device_to_1d_vector(xeq_c, xeq_c_out);
    verification::convert_1d_view_device_to_1d_vector(na_Ma, na_Ma_out);
    verification::convert_1d_view_device_to_1d_vector(nc_Mc, nc_Mc_out);
    verification::convert_1d_view_device_to_1d_vector(store, store_out);

    output.set("aer_liquid",             aer_liquid_out);
    output.set("eleliquid",              eleliquid_out);
    output.set("electrolyte_sum_liquid", electrolyte_sum_liquid_out);
    output.set("epercent_liquid",        epercent_liquid_out);
    output.set("na",                     na_out);
    output.set("nc",                     nc_out);
    output.set("xeq_a",                  xeq_a_out);
    output.set("xeq_c",                  xeq_c_out);
    output.set("na_Ma",                  na_Ma_out);
    output.set("nc_Mc",                  nc_Mc_out);
    output.set("store",                  store_out);
  });
}
