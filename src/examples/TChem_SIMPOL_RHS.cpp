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
    constexpr real_type zero=0.0;

    TChem::exec_space().print_configuration(std::cout, detail);
    TChem::host_exec_space().print_configuration(std::cout, detail);

    using device_type = typename Tines::UseThisDevice<TChem::exec_space>::type;

    std::string chemFile ="config_gas.yaml";
    std::string aerochemFile="test_SIMPOL_phase_transfer.yaml";

    /// construct kmd and use the view for testing
    TChem::KineticModelData kmd(chemFile);
    const auto kmcd = TChem::createNCAR_KineticModelConstData<device_type>(kmd);

    TChem::AerosolModelData amd = TChem::AerosolModelData(aerochemFile, kmd);
    const auto amcd = TChem::create_AerosolModelConstData<device_type>(amd);

    // using value_type = typename real_type;// Sacado::Fad::SLFad<real_type,12>;
    using SIMPOL_single_particle_type = TChem::Impl::SIMPOL_single_particle<real_type, device_type >;

    using value_type_1d_view_type = typename SIMPOL_single_particle_type::value_type_1d_view_type;
    using real_type_1d_view_type = typename SIMPOL_single_particle_type::real_type_1d_view_type;
    real_type_1d_view_type number_conc("number_conc", amcd.nParticles);
    Kokkos::deep_copy(number_conc, 1.3e6); // particle number concentration (#/cc)

    ordinal_type ntotal_species = amcd.nSpec_gas + amcd.nSpec*amcd.nParticles;

    value_type_1d_view_type state("state", ntotal_species);
    value_type_1d_view_type omega("omega", ntotal_species);
    Kokkos::deep_copy(omega, zero);
    // this view is used to copy value of omega; omegas has value and its derivices.
    real_type_1d_view_type omega_out("omega", ntotal_species);


    const auto exec_space_instance = TChem::exec_space();

    using policy_type =
          typename TChem::UseThisTeamPolicy<TChem::exec_space>::type;

    /// team policy
    ordinal_type nBatch =1;
    policy_type policy(exec_space_instance, nBatch, 1,1);
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
        for (int i_part = 0; i_part < amcd.nParticles; i_part++)
        {
          SIMPOL_single_particle_type
          ::team_invoke(member, i_part,i_simpol, t, p, number_conc, state, omega, amcd);
        }
        member.team_barrier();
        //copy values of omega
        for (ordinal_type i = 0; i < ntotal_species; i++)
          omega_out(i) =omega(i);//.val();


  });
  // we need to copy data from device to host.
  auto omega_host = Kokkos::create_mirror_view(omega_out);
  // deep_copy(out, in)
  Kokkos::deep_copy(omega_host,omega_out);

  printf("---RHSs--\n");
  std::cout << "omega(0) : "<< omega_host(0) << "\n";
  for (ordinal_type i_part = 0; i_part < amcd.nParticles; i_part++)
  {
    ordinal_type is = amcd.nSpec_gas + i_part*amcd.nSpec;
    for (ordinal_type i = 0; i < amcd.nSpec; i++)
    {
      // printf("omega(%d) %e \n",is+i,omega(is+i));
      std::cout << "omega("<<is+i<<") : "<< omega_host(is+i) << "\n";
    }
  }

  }
  Kokkos::finalize();
  return 0;
}
