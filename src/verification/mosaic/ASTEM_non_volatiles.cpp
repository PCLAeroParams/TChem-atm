#include "TChem.hpp"
#include "TChem_Impl_MOSAIC.hpp"
#include <verification.hpp>
#include "skywalker.hpp"

using device_type = typename Tines::UseThisDevice<TChem::exec_space>::type;
using real_type_1d_view = TChem::real_type_1d_view;
using ordinal_type = TChem::ordinal_type;
using namespace skywalker;
using namespace TChem;

void ASTEM_non_volatiles(Ensemble *ensemble) {
  ensemble->process([=](const Input &input, Output &output) {

    const auto mmd = TChem::Impl::MosaicModelData<device_type>();

    const auto read_view = [&](const std::string& name, const ordinal_type n) {
      real_type_1d_view v(name, n);
      verification::convert_1d_vector_to_1d_view_device(input.get_array(name), v);
      return v;
    };

    const auto dtchem         = read_view("dtchem",         1);
    const auto gas            = read_view("gas",            mmd.ngas_volatile);
    const auto jaerosolstate  = read_view("jaerosolstate",  1);
    const auto kg             = read_view("kg",             mmd.ngas_volatile);
    const auto epercent_total = read_view("epercent_total", mmd.nelectrolyte);
    const auto aer_total      = read_view("aer_total",      mmd.naer);
    const auto sumkg_h2so4    = read_view("sumkg_h2so4",    1);
    const auto sumkg_msa      = read_view("sumkg_msa",      1);
    const auto sumkg_nh3      = read_view("sumkg_nh3",      1);
    const auto sumkg_hno3     = read_view("sumkg_hno3",     1);
    const auto sumkg_hcl      = read_view("sumkg_hcl",      1);

    real_type_1d_view delta_h2so4("delta_h2so4", 1);
    real_type_1d_view delta_tmsa("delta_tmsa",   1);
    real_type_1d_view delta_nh3("delta_nh3",     1);
    real_type_1d_view delta_hno3("delta_hno3",   1);
    real_type_1d_view delta_hcl("delta_hcl",     1);
    real_type_1d_view delta_nh4("delta_nh4",     1);
    real_type_1d_view delta_nh3_max("delta_nh3_max",   1);
    real_type_1d_view delta_hno3_max("delta_hno3_max", 1);
    real_type_1d_view delta_hcl_max("delta_hcl_max",   1);

    std::string profile_name = "Verification_test_ASTEM_non_volatiles";
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
        MOSAIC::ASTEM_non_volatiles_gas(
          mmd, dtchem(0),
          sumkg_h2so4(0), sumkg_msa(0), sumkg_nh3(0), sumkg_hno3(0), sumkg_hcl(0),
          gas,
          delta_h2so4(0), delta_tmsa(0), delta_nh3(0), delta_hno3(0), delta_hcl(0));

        MOSAIC::ASTEM_non_volatiles(
          mmd, jaerosolstate(0), kg, epercent_total, aer_total,
          sumkg_h2so4(0), sumkg_msa(0), sumkg_nh3(0), sumkg_hno3(0), sumkg_hcl(0),
          delta_h2so4(0), delta_tmsa(0), delta_nh3(0), delta_hno3(0), delta_hcl(0),
          delta_nh3_max(0), delta_hno3_max(0), delta_hcl_max(0), delta_nh4(0));
      });
    });

    const auto write_view = [&](const std::string& name, const real_type_1d_view& v) {
      std::vector<real_type> out(v.extent(0));
      verification::convert_1d_view_device_to_1d_vector(v, out);
      output.set(name, out);
    };

    write_view("dtchem",         dtchem);
    write_view("gas",            gas);
    write_view("jaerosolstate",  jaerosolstate);
    write_view("kg",             kg);
    write_view("epercent_total", epercent_total);
    write_view("aer_total",      aer_total);
    write_view("sumkg_h2so4",    sumkg_h2so4);
    write_view("sumkg_msa",      sumkg_msa);
    write_view("sumkg_nh3",      sumkg_nh3);
    write_view("sumkg_hno3",     sumkg_hno3);
    write_view("sumkg_hcl",      sumkg_hcl);
    write_view("delta_h2so4",    delta_h2so4);
    write_view("delta_tmsa",     delta_tmsa);
    write_view("delta_nh3",      delta_nh3);
    write_view("delta_hno3",     delta_hno3);
    write_view("delta_hcl",      delta_hcl);
    write_view("delta_nh3_max",  delta_nh3_max);
    write_view("delta_hno3_max", delta_hno3_max);
    write_view("delta_hcl_max",  delta_hcl_max);
    write_view("delta_nh4",      delta_nh4);
  });
}
