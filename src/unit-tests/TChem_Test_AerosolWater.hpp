/* =====================================================================================
TChem-atm version 1.0
Copyright (2024) NTESS
https://github.com/sandialabs/TChem-atm

Copyright 2024 National Technology & Engineering Solutions of Sandia, LLC
(NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
Government retains certain rights in this software.

This file is part of TChem-atm. TChem-atm is open source software: you can redistribute
it and/or modify it under the terms of BSD 2-Clause License
(https://opensource.org/licenses/BSD-2-Clause). A copy of the licese is also
provided under the main directory

Questions? Contact Oscar Diaz-Ibarra at <odiazib@sandia.gov>, or
           Mike Schmidt at <mjschm@sandia.gov>, or
           Cosmin Safta at <csafta@sandia.gov>

Sandia National Laboratories, New Mexico/Livermore, NM/CA, USA
=====================================================================================
*/
#ifndef __TCHEM_TEST_AEROSOLWATER_HPP__
#define __TCHEM_TEST_AEROSOLWATER_HPP__

#include "TChem.hpp"
#include "TChem_Impl_SingleParticleAerosolWater.hpp"

using ordinal_type = TChem::ordinal_type;
using real_type = TChem::real_type;
using real_type_1d_view = TChem::real_type_1d_view;
using real_type_2d_view = TChem::real_type_2d_view;
using real_type_2d_view_host = TChem::real_type_2d_view_host;

namespace TChem {
namespace Test {
void static AerosolWater_test()
{
    using device_type = typename Tines::UseThisDevice<TChem::exec_space>::type;
    // input files need to located in same directory as the executable.
    std::string chemFile ="../examples/runs/atmospheric_chemistry/ZSR/config_gas.yaml";
    std::string aerochemFile="../examples/runs/atmospheric_chemistry/ZSR/test_ZSR.yaml";

    // construct kmd; gas phase
    printf("[main] kmd parsing ...\n");
    TChem::KineticModelData kmd(chemFile);
    const auto kmcd = TChem::createNCAR_KineticModelConstData<device_type>(kmd);
    // construct amd; aerosols
    printf("[main] amd parsing ...\n");
    TChem::AerosolModelData amd = TChem::AerosolModelData(aerochemFile, kmd);
    const auto amcd = TChem::create_AerosolModelConstData<device_type>(amd);

    using value_type = real_type;

    using aerosol_water_single_particle_type = TChem::Impl::AerosolWater_SingleParticle<value_type, device_type>;

    using value_type_1d_view_type = typename aerosol_water_single_particle_type::value_type_1d_view_type;
    using real_type_1d_view_type = typename aerosol_water_single_particle_type::real_type_1d_view_type;

    // create a 1d view number_conc with length equal to the number of computational particles
    // (set via YAML AERO_REP_SINGLE_PARTICLE entry "maximum computational particles"), currently set = 1
    real_type_1d_view_type number_conc("number_conc", amcd.nParticles);
    
    // assuming constant number concentration
    // set each computational particle to have multiplicity 1.3e6
    Kokkos::deep_copy(number_conc, 1.3e6); // particle number concentration (#/cc)

    ordinal_type ntotal_species = amcd.nSpec_gas + amcd.nSpec*amcd.nParticles;

    value_type_1d_view_type state("state", ntotal_species);
    // by default view are initialized with zeros.
    //value_type_1d_view_type omega("omega", ntotal_species);

    //auto omega_out = omega;

    // relative index of H2O_aq amongst aerosol species (not amongst both gas and aerosol + multiple particles)
    ordinal_type aqueous_water_idx = amd.aerosol_sp_name_idx_.at("H2O_aq");

    const auto exec_space_instance = TChem::exec_space();

    using policy_type =
          typename TChem::UseThisTeamPolicy<TChem::exec_space>::type;
    
    
    // FIXME:
    // is this a bug in nvcc?
    //https://github.com/google/googletest/issues/4104
 #if  1
    // team policy
    ordinal_type nBatch =1;
    policy_type policy(exec_space_instance, nBatch, Kokkos::AUTO());
    Kokkos::parallel_for
      ("AEROSOL WATER",
       policy,
       KOKKOS_LAMBDA(const typename policy_type::member_type& member) {

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
        Kokkos::parallel_for(
          Kokkos::TeamThreadRange(member, amcd.nParticles),
          [&](const ordinal_type& i_part) {
              aerosol_water_single_particle_type
              ::team_invoke(member, i_part, number_conc, state, amcd, rh_idx, aqueous_water_idx);
            });
        member.team_barrier();
        */

        
        ordinal_type i_part = 0;
       
        aerosol_water_single_particle_type::team_invoke(member, i_part, number_conc, state, amcd, rh_idx, aqueous_water_idx);
        //printf("[TChem_ZSR::main] RH %f\n", state[2]);
        //printf("[TChem_ZSR::main] Total aerosol water content %f\n\n", state[3]);
     
        printf("%f,%f\n", state(rh_idx), state(amcd.nSpec_gas + amcd.nSpec*i_part + aqueous_water_idx));
      
      }
        
  });
  

#endif

  /*
  // we need to copy data from device to host.
  auto omega_host = Kokkos::create_mirror_view_and_copy(TChem::host_exec_space(), omega_out);
  // save results to a file.
  std::string outputFile ="Simpol_RHS_DEVICE.txt";
  FILE *fout = fopen(outputFile.c_str(), "w");
  TChem::Test::writeReactionRates(outputFile, omega_host.extent(0), omega_host);
  fclose(fout);

  
    */

}
}// namespace Test


}// namespace TChem


TEST(AerosolWater, Device)
{
  TChem::Test::AerosolWater_test();
}


/*
TEST(AerosolWater, verification_device)
{
 /// pass test with a relative error of
  const real_type threshold =1e-12;
  EXPECT_TRUE(TChem::Test::compareFilesValues("Simpol_RHS_DEVICE.txt",
          "references/partmc_integration/simpol/Simpol_RHS_HOST.txt",threshold)
        );
}
*/

#endif
