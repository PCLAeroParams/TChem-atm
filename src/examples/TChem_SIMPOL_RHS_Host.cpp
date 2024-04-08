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

    using host_device_type = typename Tines::UseThisDevice<TChem::host_exec_space>::type;

    std::string chemFile ="config_gas.yaml";
    std::string aerochemFile="test_SIMPOL_phase_transfer.yaml";

    // construct kmd and use the view for testing
    TChem::KineticModelData kmd(chemFile);
    const auto kmcd = TChem::createNCAR_KineticModelConstData<host_device_type>(kmd);

    TChem::AerosolModelData amd = TChem::AerosolModelData(aerochemFile, kmd);
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

    SIMPOL_constant_type::team_invoke(member, t,
     p, alpha, mfp_m, KGM3_TO_PPM_, equil_constant, amcd.simpol_params(0));
#if defined(TCHEM_ENABLE_SERIAL_TEST_OUTPUT)
    printf("amcd.simpol_params(0).B1 %e \n", amcd.simpol_params(0).B1);
    printf("alpha %e \n", alpha);
    printf("mfp_m %e \n", mfp_m);
    printf("KGM3_TO_PPM_ %e \n", KGM3_TO_PPM_);
    printf("equil_constant %e \n", equil_constant);
#endif
    using value_type_1d_view_type = typename SIMPOL_single_particle_type::value_type_1d_view_type;
    using real_type_1d_view_type = typename SIMPOL_single_particle_type::real_type_1d_view_type;
    real_type_1d_view_type number_conc("number_conc", amcd.nParticles);
    Kokkos::deep_copy(number_conc, 1.3e6); // particle number concentration (#/cc)
    // initial concentration
    real_type ethanol=0.1;
    real_type ethanol_aq = 1.0e-8 ;
    real_type H2O_aq = 1.4e-2;
    ordinal_type ntotal_species = amcd.nSpec_gas + amcd.nSpec*amcd.nParticles;

    value_type_1d_view_type state("state", ntotal_species);
    // gas species
    state(0) = ethanol;
    // aerosol_species
    for (int i_part = 0; i_part < amcd.nParticles; i_part++)
    {
      state(1+amcd.nSpec*i_part) = ethanol_aq/number_conc(i_part);
      state(2+amcd.nSpec*i_part) = H2O_aq/number_conc(i_part);
    }

    value_type_1d_view_type omega("omega", ntotal_species);

    ordinal_type i_simpol=0;
    for (int i_part = 0; i_part < amcd.nParticles; i_part++)
    {
#if defined(TCHEM_ENABLE_SERIAL_TEST_OUTPUT)
    printf("----Working on particle No %d ---\n", i_part);
#endif
    SIMPOL_single_particle_type
    ::team_invoke(member, i_part,i_simpol, t, p, number_conc, state, omega, amcd);
    }

  //save results to a file.
  std::string outputFile ="Simpol_RHS_HOST.txt";
  FILE *fout = fopen(outputFile.c_str(), "w");
  TChem::Test::writeReactionRates(outputFile, omega.extent(0), omega);
  fclose(fout);

  }
  Kokkos::finalize();
  return 0;
}
