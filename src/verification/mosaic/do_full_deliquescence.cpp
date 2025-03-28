#include "TChem.hpp"
#include "TChem_Impl_MOSAIC.hpp"
#include <verification.hpp>
#include "skywalker.hpp"

using device_type = typename Tines::UseThisDevice<TChem::exec_space>::type;
using real_type_1d_view = TChem::real_type_1d_view;
using ordinal_type = TChem::ordinal_type;
using namespace skywalker;
using namespace TChem;

void do_full_deliquescence(Ensemble *ensemble) {
    ensemble->process([=](const Input &input, Output &output) {

        const auto electrolyte_db = input.get_array("electrolyte");
        const auto aer_db = input.get_array("aer");

        const ordinal_type n = 3;
        const auto nsize_electrolyte = static_cast<ordinal_type>(electrolyte_db.size())/n;
        const auto nsize_aero = static_cast<ordinal_type>(aer_db.size())/n;

        const auto mmd = TChem::Impl::MosaicModelData<device_type>();

        real_type_1d_view electrolyte_solid("electrolyte_solid", nsize_electrolyte);
        real_type_1d_view electrolyte_liquid("electrolyte_liquid", nsize_electrolyte);
        real_type_1d_view electrolyte_total("electrolyte_total", nsize_electrolyte);
        real_type_1d_view aer_solid("aer_solid", nsize_aero);
        real_type_1d_view aer_liquid("aer_liquid", nsize_aero);
        real_type_1d_view aer_total("aer_total", nsize_aero);

        std::vector<std::vector<real_type>> electrolyte_db_2d;
        for (size_t i = 0; i < n*nsize_electrolyte; i += nsize_electrolyte) {
            std::vector<real_type> electrolyte_temp(electrolyte_db.begin() + i, electrolyte_db.begin() + i + nsize_electrolyte);
            electrolyte_db_2d.push_back(electrolyte_temp);
        }

        verification::convert_1d_vector_to_1d_view_device(electrolyte_db_2d[0], electrolyte_solid);
        verification::convert_1d_vector_to_1d_view_device(electrolyte_db_2d[1], electrolyte_liquid);
        verification::convert_1d_vector_to_1d_view_device(electrolyte_db_2d[2], electrolyte_total);

        std::vector<std::vector<real_type>> aer_db_2d;
        for (size_t i = 0; i < n*nsize_aero; i += nsize_aero) {
            std::vector<real_type> aer_temp(aer_db.begin() + i, aer_db.begin() + i + nsize_aero);
            aer_db_2d.push_back(aer_temp);
        }

        verification::convert_1d_vector_to_1d_view_device(aer_db_2d[0], aer_solid);
        verification::convert_1d_vector_to_1d_view_device(aer_db_2d[1], aer_liquid);
        verification::convert_1d_vector_to_1d_view_device(aer_db_2d[2], aer_total);

        std::string profile_name ="Verification_test_do_full_deliquescence";
        using policy_type =
            typename TChem::UseThisTeamPolicy<TChem::exec_space>::type;
        const auto exec_space_instance = TChem::exec_space();
        const auto host_exec_space = TChem::host_exec_space();
        policy_type policy(exec_space_instance, 1, Kokkos::AUTO());

        Kokkos::parallel_for(
        profile_name,
        policy,
        KOKKOS_LAMBDA(const typename policy_type::member_type& member) {

            TChem::Impl::MOSAIC<real_type, device_type>::do_full_deliquescence(
            mmd,
            electrolyte_solid, electrolyte_liquid, electrolyte_total,
            aer_solid, aer_liquid, aer_total);
        });

            verification::convert_1d_view_device_to_1d_vector(electrolyte_solid, electrolyte_db_2d[0]);
            verification::convert_1d_view_device_to_1d_vector(electrolyte_liquid, electrolyte_db_2d[1]);
            verification::convert_1d_view_device_to_1d_vector(electrolyte_total, electrolyte_db_2d[2]);

            std::vector<real_type> electrolyte_flattened;
            for (const auto& row : electrolyte_db_2d) {
                electrolyte_flattened.insert(electrolyte_flattened.end(), row.begin(), row.end());
            }
            output.set("electrolyte", electrolyte_flattened);

            verification::convert_1d_view_device_to_1d_vector(aer_solid, aer_db_2d[0]);
            verification::convert_1d_view_device_to_1d_vector(aer_liquid, aer_db_2d[1]);
            verification::convert_1d_view_device_to_1d_vector(aer_total, aer_db_2d[2]);

            std::vector<real_type> aer_flattened;
            for (const auto& row : aer_db_2d) {
                aer_flattened.insert(aer_flattened.end(), row.begin(), row.end());
            }
            output.set("aer", aer_flattened);

    });
}