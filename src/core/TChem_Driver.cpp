#include "TChem_Driver.hpp"
#include "TChem.hpp"
#include "TChem_KineticModelNCAR_ConstData.hpp"
#include "TChem_CommandLineParser.hpp" 

using real_type = TChem::real_type;
using ordinal_type = TChem::ordinal_type;

  using exec_space = Kokkos::DefaultExecutionSpace;
  using host_exec_space = Kokkos::DefaultHostExecutionSpace;

  using host_device_type = typename Tines::UseThisDevice<TChem::host_exec_space>::type;
  using device_type      = typename Tines::UseThisDevice<exec_space>::type;

static TChem::Driver *g_tchem = nullptr;

void normal_vec(int n, real_type x[]) {
  for (int i = 0; i<n; ++i) x[i] = 1.0; //rnorm();
}

void call_something(int n){
   printf("Printing in TChem c++ code\n");
   printf("%f\n", TChem::CONV_PPM);
   printf("%f\n", TChem::PI());
}

void do_something(int n){
//   printf("Kokkos initialized: %s\n", Kokkos::is_initialized ? "true" : "false");
}

void initialize_kokkos(const char * chemFile){

  g_tchem = new TChem::Driver();

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

  g_tchem->createGasKineticModel(chemFile);
  g_tchem->createGasKineticModelConstData();

/*  TChem::KineticModelData kmd(chemFile);
  const auto kmcd = TChem::createNCAR_KineticModelConstData<device_type>(kmd);
  const ordinal_type stateVecDim = TChem::Impl::getStateVectorSize(kmcd.nSpec);
  printf("Testing input %s \n", chemFile);
  printf("Number of Species %d \n", kmcd.nSpec);
  printf("Number of Reactions %d \n", kmcd.nReac);

  printf("done making the kmd\n"); */

/*  std::string species_names[kmcd.nSpec];

  printf("Get the gas species\n");

  const auto speciesNamesHost = Kokkos::create_mirror_view(kmcd.speciesNames);
  Kokkos::deep_copy(speciesNamesHost, kmcd.speciesNames);
  for (ordinal_type k = 0; k < kmcd.nSpec; k++)
  {
     species_names[k] = &speciesNamesHost(k,0);
     printf("%s \n", species_names[k].c_str());
  }

  printf("more testing\n"); */

/*  TChem::Driver tchem;
  tchem.createGasKineticModel(chemFile);
  const auto kmcd_new = TChem::createNCAR_KineticModelConstData<device_type>(tchem._kmd);
  printf("Number of Species %d \n", kmcd_new.nSpec);
  printf("Number of Reactions %d \n", kmcd_new.nReac);
  printf("%s \n", tchem._chem_file.c_str()); */
  
}

void finalize_kokkos(){

  g_tchem->freeAll();
  delete g_tchem;
  printf("Cleaning up TChem");
  Kokkos::finalize();
  printf("Cleaning up Kokkos\n");
}

void TChem::Driver::createGasKineticModel(const std::string &chem_file) {
  _chem_file = chem_file;
  _kmd = KineticModelData(_chem_file);
}

// FIXME: Add some checks
void TChem::Driver::createGasKineticModelConstData() {
  printf("Creating kmcd \n");
  using interf_host_device_type = typename Tines::UseThisDevice<host_exec_space>::type;
  _kmcd = TChem::createNCAR_KineticModelConstData<interf_device_type>(_kmd);

  printf("Number of Species %d \n", _kmcd.nSpec);
  printf("Number of Reactions %d \n", _kmcd.nReac);

  std::string species_names[_kmcd.nSpec];

  printf("Get the gas species\n");

  const auto speciesNamesHost = Kokkos::create_mirror_view(_kmcd.speciesNames);
  Kokkos::deep_copy(speciesNamesHost, _kmcd.speciesNames);
  for (ordinal_type k = 0; k < _kmcd.nSpec; k++)
  {
     species_names[k] = &speciesNamesHost(k,0);
     printf("%s \n", species_names[k].c_str());
  }

  printf("End creating kmcd \n");

}

// FIXME: Free things
void TChem::Driver::freeAll() {
  g_tchem->freeGasKineticModel();
}

void TChem::Driver::freeGasKineticModel() {
  _chem_file = std::string();
  _therm_file = std::string();
}
