#include "TChem.hpp"
#include "TChem_Impl_SingleParticleAerosolWater.hpp"
//#include "TChem_Impl_SIMPOL_phase_transfer.hpp"
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
    const bool detail = false; //  this doesnt really seem to do anything

    // Currently the execution space and host execution space are the same thing (I think)
    printf("[TChem_ZSR::main] Execution space details:\n");
    TChem::exec_space().print_configuration(std::cout, detail);
    printf("\n");

    printf("[TChem_ZSR::main] Host execution space details:\n");
    TChem::host_exec_space().print_configuration(std::cout, detail);
    printf("\n");

    using host_device_type = typename Tines::UseThisDevice<TChem::host_exec_space>::type;

    std::string chemFile ="config_gas.yaml";
    std::string aerochemFile="test_ZSR.yaml";

    // construct kmd and use the view for testing
    printf("[main] kmd parsing ...\n");
    TChem::KineticModelData kmd(chemFile);
    const auto kmcd = TChem::createNCAR_KineticModelConstData<host_device_type>(kmd);

    printf("[main] amd parsing ...\n");
    // TODO: Add parsing for ZSR model data
    TChem::AerosolModelData amd = TChem::AerosolModelData(aerochemFile, kmd);
    const auto amcd = TChem::create_AerosolModelConstData<host_device_type>(amd);

    using aerosol_water_single_particle_type = TChem::Impl::AerosolWater_SingleParticle<real_type, host_device_type>;
    
    // create a 1d view number_conc with length equal to the number of computational particles
    // (set via YAML AERO_REP_SINGLE_PARTICLE entry "maximum computational particles"), currently set = 1
    real_type_1d_view number_conc("number_conc", amcd.nParticles);

    // assuming constant number concentration
    // set each computational particle to have multiplicity 1.3e6
    Kokkos::deep_copy(number_conc, 1.3e6); // particle number concentration

    // Just run serial right now for simplicity
    const auto member = Tines::HostSerialTeamMember();

    ordinal_type ntotal_species = amcd.nSpec_gas + amcd.nSpec*amcd.nParticles;

    //printf("[TChem_ZSR::main] Number of gas species: %d\n", amcd.nSpec_gas);
    //printf("[TChem_ZSR::main] Number of aerosol species: %d\n", amcd.nSpec);
    // TChem's state data structure
    real_type_1d_view state("state", ntotal_species);

    //real_type state[amd.aerosol_sp_name_idx_.size()+1]; // add 1 for gas phase water
    //printf("[TChem_ZSR::main] Number of aerosol species according to aerosol_sp_name_idx_: %d\n", amd.aerosol_sp_name_idx_.size());
    //int rh_index = amd.aerosol_sp_name_idx_.size();
    int i_part = 0;
    //int number_conc = 1.0;

    
    //ordinal_type rh_idx = amd.aerosol_sp_name_idx_.size();

    // relative index of H2O_aq amongst aerosol species (not amongst both gas and aerosol + multiple particles)
    ordinal_type aqueous_water_idx = amd.aerosol_sp_name_idx_.at("H2O_aq");


    //printf("[TChem_ZSR::main] index of H2O: %d\n", rh_idx);
    //printf("[TChem_ZSR::main] index of H2O_aq: %d\n", aqueous_water_idx);
    //printf("[TChem_ZSR::main] index of Cl_m: %d\n", amd.aerosol_sp_name_idx_.at("Cl_m"));
    //printf("[TChem_ZSR::main] index of Ca_pp: %d\n", amd.aerosol_sp_name_idx_.at("Ca_pp"));
    //printf("[TChem_ZSR::main] index of Na_p: %d\n", amd.aerosol_sp_name_idx_.at("Na_p"));
            
    // Loop over RH
    for (int i; i<101; i++){

        // gas species
        ordinal_type rh_idx = 0;
        state(rh_idx) = i/100.0; //  H2O (gas phase)

        // aerosol species
        for (int i_part = 0; i_part < amcd.nParticles; i_part++){
            state(1+amcd.nSpec*i_part) = 2.5; // Na_p
            state(2+amcd.nSpec*i_part) = 5.3; // Cl_m
            state(3+amcd.nSpec*i_part) = 1.3; // Ca_pp
            state(4+amcd.nSpec*i_part) = 0.0; // H2O_aq
        }// i_part
       member.team_barrier();

        /*
        state[amd.aerosol_sp_name_idx_.at("Cl_m")] = 5.3;    // Cl
        state[amd.aerosol_sp_name_idx_.at("Ca_pp")] = 1.3;   // Ca
        state[rh_idx] = i/100.0; // H2O (gas phase) (note gas phase H2O is not part of the species mapping)
        state[amd.aerosol_sp_name_idx_.at("H2O_aq")] = 0;    // H2O (aerosol phase)
        */ 

        // TODO use parallel for over particles
        int i_part = 0;
        aerosol_water_single_particle_type::team_invoke(member, i_part, number_conc, state, amcd, rh_idx, aqueous_water_idx);
        //printf("[TChem_ZSR::main] RH %f\n", state[2]);
        //printf("[TChem_ZSR::main] Total aerosol water content %f\n\n", state[3]);
        printf("%f,%f\n", state(rh_idx), state(amcd.nSpec_gas + amcd.nSpec*i_part + aqueous_water_idx));
      
    }
    

    
        
    
#if defined(TCHEM_ENABLE_SERIAL_TEST_OUTPUT)
    printf("Some stuff in the serial test output")

#endif

    /*
    
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
    // aerosol species
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
  */ 

  }
  Kokkos::finalize();
  return 0;
}