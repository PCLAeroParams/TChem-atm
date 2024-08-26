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

#include "TChem_Util.hpp"
using ordinal_type = TChem::ordinal_type;
using real_type = TChem::real_type;

#include "TChem_Test_ReactionRates.hpp"
#include "TChem_Test_RateofProgress.hpp"
#include "TChem_Test_NetProductionRates.hpp"
#include "TChem_Test_AtmosphericChemistryE3SM.hpp"
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
