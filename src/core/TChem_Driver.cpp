#include "TChem_Driver.hpp"
#include "TChem.hpp"
#include "TChem_KineticModelNCAR_ConstData.hpp"
#include "TChem_CommandLineParser.hpp" 

using real_type = TChem::real_type;

void normal_vec(int n, real_type x[]) {
  for (int i = 0; i<n; ++i) x[i] = 1.0; //rnorm();
}

void call_something(int n){
   std::string chemFile("chem.yaml");
   printf("Printing in TChem c++ code\n");
   printf("%f\n", TChem::CONV_PPM);
   printf("%f\n", TChem::PI());
}

void do_something(int n){
   printf("Printing in TChem c++ code doing something\n");
}

void initialize_kokkos(){

  Kokkos::InitializationSettings settings;
  settings.set_num_threads(8);
  settings.set_device_id(1);
  printf("Starting Kokkos\n");
  Kokkos::initialize(settings);

  using exec_space = Kokkos::DefaultExecutionSpace;
  using host_exec_space = Kokkos::DefaultHostExecutionSpace;

  using host_device_type = typename Tines::UseThisDevice<TChem::host_exec_space>::type;
  using device_type      = typename Tines::UseThisDevice<exec_space>::type;

  using time_integrator_cvode_type = Tines::TimeIntegratorCVODE<real_type,host_device_type>;
  Tines::value_type_1d_view<time_integrator_cvode_type,host_device_type> cvodes;
  using real_type_2d_view_host_type = Tines::value_type_2d_view<real_type,host_device_type>;
  using real_type_1d_view_type = Tines::value_type_1d_view<real_type,device_type>;
  using real_type_2d_view_type = Tines::value_type_2d_view<real_type,device_type>;

  using real_type_1d_view_host_type = Tines::value_type_1d_view<real_type,host_device_type>;
  using real_type_2d_view_host_type = Tines::value_type_2d_view<real_type,host_device_type>;

  using kinetic_model_type = TChem::KineticModelNCAR_ConstData<device_type>;
  using kinetic_model_host_type = TChem::KineticModelNCAR_ConstData<host_device_type>;

  using ordinal_type = TChem::ordinal_type;

  std::string chemFile("chem.yaml");

  TChem::KineticModelData kmd(chemFile);
    const auto kmcd = TChem::createNCAR_KineticModelConstData<device_type>(kmd);

    const ordinal_type stateVecDim = TChem::Impl::getStateVectorSize(kmcd.nSpec);

    printf("Number of Species %d \n", kmcd.nSpec);
    printf("Number of Reactions %d \n", kmcd.nReac);

    printf("done making the kmd\n");

    const auto speciesNamesHost = Kokkos::create_mirror_view(kmcd.speciesNames);
    Kokkos::deep_copy(speciesNamesHost, kmcd.speciesNames);
    for (ordinal_type k = 0; k < kmcd.nSpec; k++)
       printf("%s \n", &speciesNamesHost(k, 0));
  int nSpec = kmcd.nSpec;
  
}

void finalize_kokkos(){

  Kokkos::finalize();
  printf("Cleaning up Kokkos\n");
}
