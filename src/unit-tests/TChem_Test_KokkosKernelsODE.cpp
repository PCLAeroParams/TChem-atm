/* =====================================================================================
TChem-atm version 1.0
Copyright (2024) NTESS
https://github.com/sandialabs/TChem-atm

Copyright 2024 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains
certain rights in this software.

This file is part of TChem-atm. TChem-atm is open source software: you can redistribute it
and/or modify it under the terms of BSD 2-Clause License
(https://opensource.org/licenses/BSD-2-Clause). A copy of the licese is also
provided under the main directory

Questions? Contact Oscar Diaz-Ibarra at <odiazib@sandia.gov>, or
           Mike Schmidt at <mjschm@sandia.gov>, or
           Cosmin Safta at <csafta@sandia.gov>

Sandia National Laboratories, New Mexico/Livermore, NM/CA, USA
===================================================================================== */
#include <gtest/gtest.h>
#include "TChem.hpp"
int TestExamplesInternal(std::string path_inputs, std::string exec) {
  {
    int r_val(0);
    std::string logfile(exec + ".test-log");
    std::string rm=("rm -f " + logfile);
    const auto rm_c_str = rm.c_str();
    r_val = std::system(rm_c_str); 
    TCHEM_CHECK_ERROR(r_val, "system call rm -f returns non-zero return value");

    std::string chemFile=" --chemfile="+path_inputs+"config_gas.yaml ";
    std::string aerofile="--aerofile="+path_inputs+"test_SIMPOL_phase_transfer.yaml ";
    std::string inputfile_particles="--inputfile_particles="+path_inputs+"scenario_conditions_particle.yaml ";
    std::string outputfile="--outputfile=full_gas_KK.dat ";
    std::string dtmin="--dtmin=10 ";
    std::string tend="--tend=50 "; 
    
    std::string invoke=(exec + chemFile + aerofile +
                        inputfile_particles + outputfile + dtmin + tend + " > " + logfile);
    const auto invoke_c_str = invoke.c_str();
    printf("Tines testing : %s\n", invoke_c_str);
    r_val = std::system(invoke_c_str);
    TCHEM_CHECK_ERROR(r_val, "system call example executable returns non-zero return value");
    std::ifstream file(logfile);
    for (std::string line; getline(file, line); ) {
      printf("%s\n", line.c_str());
      EXPECT_TRUE(line.find("FAIL") == std::string::npos);
    }
  }
  return 0;
}

TEST(TimeIntegration,KokkosKernelsBDF) {
  TestExamplesInternal(
    "../examples/runs/atmospheric_chemistry/Simpol/",
    "../examples/TChem_AerosolChemistry_KokkosKernels.x");
}

int
main(int argc, char* argv[])
{
  int r_val(0);
  Kokkos::initialize(argc, argv);
  {
    const bool detail = false;

    TChem::exec_space().print_configuration(std::cout, detail);
    TChem::host_exec_space().print_configuration(std::cout, detail);

    ::testing::InitGoogleTest(&argc, argv);
    r_val = RUN_ALL_TESTS();
  }
  Kokkos::finalize();

  return r_val;
}
