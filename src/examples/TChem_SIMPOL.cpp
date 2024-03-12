#include "TChem.hpp"
#include "TChem_Impl_SIMPOL_constant.hpp"
#include "TChem_Impl_SIMPOL_phase_transfer.hpp"
using ordinal_type = TChem::ordinal_type;
using real_type = TChem::real_type;
using real_type_1d_view = TChem::real_type_1d_view;
using real_type_2d_view = TChem::real_type_2d_view;
using real_type_2d_view_host = TChem::real_type_2d_view_host;

int
main(int argc, char* argv[])
{
  Kokkos::initialize(argc, argv);
  {
    const bool detail = false;

    TChem::exec_space().print_configuration(std::cout, detail);
    TChem::host_exec_space().print_configuration(std::cout, detail);

    using host_device_type      = typename Tines::UseThisDevice<TChem::host_exec_space>::type;

    std::string chemFile="test.yaml";

    /// construct kmd and use the view for testing
    TChem::AerosolModelData amd = TChem::AerosolModelData(chemFile);
    const auto amcd = TChem::create_AerosolModelConstData<host_device_type>(amd);


    using SIMPOL_single_particle_type = TChem::Impl::SIMPOL_single_particle<real_type, host_device_type >;
    using SIMPOL_constant_type = TChem::Impl::SIMPOL_constant<real_type, host_device_type >;

    const auto member = Tines::HostSerialTeamMember();
    // inputs
    real_type t = 272.5; // K
    real_type p = 101253.3; // pa
    // outputs
    real_type alpha=0;
    real_type mfp_m=0;
    real_type KGM3_TO_PPM_=0;
    real_type equil_constant=0;

    SIMPOL_constant_type::team_invoke(member, t, p, alpha, mfp_m, KGM3_TO_PPM_, equil_constant);

    printf("alpha %e \n", alpha);
    printf("mfp_m %e \n", mfp_m);
    printf("KGM3_TO_PPM_ %e \n", KGM3_TO_PPM_);
    printf("equil_constant %e \n", equil_constant);

    using value_type_1d_view_type = typename SIMPOL_single_particle_type::value_type_1d_view_type;
    ordinal_type nSpec=3;
    real_type number_conc = 1.3e6; // particle number concentration (#/cc)
    // initial concentration
    real_type ethanol=0.1;
    real_type ethanol_aq = 1.0e-8 ;
    real_type H2O_aq = 1.4e-2;

    value_type_1d_view_type state("state", nSpec);
    state(0) = ethanol;
    state(1) = ethanol_aq/number_conc;
    state(2) = H2O_aq/number_conc;

    using real_type_1d_view_type = typename SIMPOL_single_particle_type::real_type_1d_view_type;
    real_type_1d_view_type molecular_weigths("molecular_weigths", nSpec);
    real_type_1d_view_type density("density", nSpec);
    Kokkos::deep_copy(density, 1e3);
    value_type_1d_view_type omega("omega", nSpec);

    molecular_weigths(0)=0.04607;
    molecular_weigths(1)=0.04607;
    molecular_weigths(2)=0.01801;


    using ordinal_type_1d_view_type = typename SIMPOL_single_particle_type::ordinal_type_1d_view_type;
    ordinal_type_1d_view_type species_type("species_type", nSpec);
    Kokkos::deep_copy(species_type, 1);
    species_type(0)=-1;


    SIMPOL_single_particle_type
    ::invoke_team(member, t, p, number_conc, state, molecular_weigths,
                  density, species_type, omega);





  }
  Kokkos::finalize();
  return 0;
}
