#include "TChem.hpp"
#include "TChem_Impl_SIMPOL_constant.hpp"
#include "TChem_Impl_SIMPOL_phase_transfer.hpp"
using ordinal_type = TChem::ordinal_type;
using real_type = TChem::real_type;
using real_type_1d_view = TChem::real_type_1d_view;
using real_type_2d_view = TChem::real_type_2d_view;
using real_type_2d_view_host = TChem::real_type_2d_view_host;

#define TCHEM_TEST_IMPL_SACADO
int
main(int argc, char* argv[])
{
  Kokkos::initialize(argc, argv);
  {
    const bool detail = false;
    constexpr real_type zero=0.0;

    TChem::exec_space().print_configuration(std::cout, detail);
    TChem::host_exec_space().print_configuration(std::cout, detail);

    using device_type = typename Tines::UseThisDevice<TChem::exec_space>::type;
    // input files need to located in same directory as the executable.
    std::string chemFile ="config_gas.yaml";
    std::string aerochemFile="test_SIMPOL_phase_transfer.yaml";

    /// construct kmd; gas phase
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
    using SIMPOL_single_particle_type = TChem::Impl::SIMPOL_single_particle<value_type, device_type >;


    using value_type_1d_view_type = typename SIMPOL_single_particle_type::value_type_1d_view_type;
    using real_type_1d_view_type = typename SIMPOL_single_particle_type::real_type_1d_view_type;

    real_type_1d_view_type number_conc("number_conc", amcd.nParticles);
    // assuming constant number concentration
    Kokkos::deep_copy(number_conc, 1.3e6); // particle number concentration (#/cc)

    ordinal_type ntotal_species = amcd.nSpec_gas + amcd.nSpec*amcd.nParticles;

    value_type_1d_view_type state("state", ntotal_species);
    // by defaul view are initialized with zeros.
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
      ("SIMPOL RHS",
       policy,
       KOKKOS_LAMBDA(const typename policy_type::member_type& member) {
       real_type t = 272.5; // K
       real_type p = 101253.3; // pa

       // initial concentration
       real_type ethanol=0.1;
       real_type ethanol_aq = 1.0e-8 ;
       real_type H2O_aq = 1.4e-2;

       // gas species
       state(0) = ethanol;
       // aerosol_species
       for (int i = 0; i < amcd.nParticles; i++)
       {
        state(1+amcd.nSpec*i) = ethanol_aq/number_conc(i);
        state(2+amcd.nSpec*i) = H2O_aq/number_conc(i);
       }
       member.team_barrier();

       ordinal_type i_simpol=0;
       Kokkos::parallel_for(
      Kokkos::TeamThreadRange(member, amcd.nParticles),
       [&](const ordinal_type& i_part) {
           SIMPOL_single_particle_type
          ::team_invoke(member, i_part, i_simpol, t, p, number_conc, state, omega, amcd);
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
  //save results to a file.
  std::string outputFile ="Simpol_RHS.txt";
  FILE *fout = fopen(outputFile.c_str(), "w");
  TChem::Test::writeReactionRates(outputFile, omega_host.extent(0), omega_host);
  fclose(fout);

  }
  Kokkos::finalize();
  return 0;
}
