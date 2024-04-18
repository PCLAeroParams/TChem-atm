#ifndef __TCHEM_TEST_TOYPROCESS_HPP__
#define __TCHEM_TEST_TOYPROCESS_HPP__

#include "TChem.hpp"
#include "TChem_Impl_ToyProcess.hpp"
using ordinal_type = TChem::ordinal_type;
using real_type = TChem::real_type;
using real_type_1d_view = TChem::real_type_1d_view;
using real_type_2d_view = TChem::real_type_2d_view;
using real_type_2d_view_host = TChem::real_type_2d_view_host;

#define TCHEM_TEST_IMPL_SACADO

TEST(ToyProcess, Device)
{
    constexpr real_type zero=0.0;
    using device_type = typename Tines::UseThisDevice<TChem::exec_space>::type;
    // input files need to located in same directory as the executable.
    std::string chemFile ="../examples/runs/atmospheric_chemistry/ToyProcess/config_gas.yaml";
    std::string aerochemFile="../examples/runs/atmospheric_chemistry/ToyProcess/toy_process.yaml";

    // construct kmd; gas phase
    TChem::KineticModelData kmd(chemFile);
    const auto kmcd = TChem::createNCAR_KineticModelConstData<device_type>(kmd);
    // construct amd; aerosols
    TChem::AerosolModelData amd = TChem::AerosolModelData(aerochemFile, kmd);
    const auto amcd = TChem::create_AerosolModelConstData<device_type>(amd);

#if defined(TCHEM_TEST_IMPL_SACADO)
    // length of value type is needed at compilation time.
    using value_type = Sacado::Fad::SLFad<real_type,128>;
#else
   using value_type = real_type;
#endif
    using toy_process_type = TChem::Impl::ToyProcess<value_type, device_type >;

    using value_type_1d_view_type = typename toy_process_type::value_type_1d_view_type;
    using real_type_1d_view_type = typename toy_process_type::real_type_1d_view_type;

    real_type_1d_view_type number_conc("number_conc", amcd.nParticles);
    // assuming constant number concentration
    Kokkos::deep_copy(number_conc, 1.3e6); // particle number concentration (#/cc)

    ordinal_type ntotal_species = amcd.nSpec_gas + amcd.nSpec*amcd.nParticles;

    value_type_1d_view_type state("state", ntotal_species);
    // by default view are initialized with zeros.
    value_type_1d_view_type omega("omega", ntotal_species);

#if defined(TCHEM_TEST_IMPL_SACADO)
    // this view is used to copy value of omega; omegas has value and its derivatives.
    real_type_1d_view_type omega_out("omega", ntotal_species);
#else
    auto omega_out = omega;
#endif
    const auto exec_space_instance = TChem::exec_space();

    using policy_type =
          typename TChem::UseThisTeamPolicy<TChem::exec_space>::type;

    // team policy
    ordinal_type nBatch =1;
    policy_type policy(exec_space_instance, nBatch, Kokkos::AUTO());
    Kokkos::parallel_for
      ("Toy process",
       policy,
       KOKKOS_LAMBDA(const typename policy_type::member_type& member) {
       real_type t = 272.5; // K
       real_type p = 101253.3; // pa

       // initial concentration
       real_type A = 0.1;
       real_type B = 1.0e-8 ;

       // gas species
       state(0) = A;
       // aerosol species
       for (int i_part = 0; i_part < amcd.nParticles; i_part++)
       {
        state(1+amcd.nSpec*i_part) = B/number_conc(i_part);
       }// i_part
       member.team_barrier();

       ordinal_type itoy=0;
       Kokkos::parallel_for(
      Kokkos::TeamThreadRange(member, amcd.nParticles),
       [&](const ordinal_type& i_part) {
           toy_process_type
          ::team_invoke(member, i_part, itoy, t, p, number_conc, state, omega, amcd);
        });
        member.team_barrier();
#if defined(TCHEM_TEST_IMPL_SACADO)
        //copy values of omega
        for (ordinal_type i = 0; i < ntotal_species; i++)
          omega_out(i) = omega(i).val();
#endif

  });
  // we need to copy data from device to host.
  auto omega_host = Kokkos::create_mirror_view_and_copy(TChem::host_exec_space(), omega_out);
  // save results to a file.
  std::string outputFile ="ToyProcess_RHS_DEVICE.txt";
  FILE *fout = fopen(outputFile.c_str(), "w");
  TChem::Test::writeReactionRates(outputFile, omega_host.extent(0), omega_host);
  fclose(fout);
}

TEST(SimpolRHS, verification_device)
{
 /// pass test with a relative error of
  const real_type threshold =1e-12;
  EXPECT_TRUE(TChem::Test::compareFilesValues("ToyProcess_RHS_DEVICE.txt",
          "references/partmc_integration/toy_process/ToyProcess_RHS_HOST.txt",threshold)
        );
}

#endif
