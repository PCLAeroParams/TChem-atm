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

void initialize(const char * chemFile){

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
  g_tchem->createStateVector();
  g_tchem->getSpeciesNames();
  g_tchem->getLengthOfStateVector();

//  real_type_2d_view_host_type state;
//  g_tchem->getStateVectorHost();
//  g_tchem->getAllStateVectorHost(state);
  
}

void finalize(){

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

void TChem::Driver::createStateVector() {
  // FIXME: add error checking
  const ordinal_type len = TChem::Impl::getStateVectorSize(_kmcd.nSpec);
  const ordinal_type nBatch = 1;
  _state = real_type_2d_view_host("state dev", nBatch, len);
  for (ordinal_type k = 0; k < _kmcd.nSpec; k++){
    _state(0,k) = k*k;
  }
}

auto TChem::Driver::getStateVector() { //TChem::real_type *view){
  printf("in TChem::Driver::getStateVector\n");
  auto state_at_0 = Kokkos::subview(_state, 0, Kokkos::ALL); 
  for (ordinal_type k = 0; k < _kmcd.nSpec; k++){
     printf("%e\n", state_at_0(k));
  }
  return state_at_0;
}

void TChem_getStateVector(TChem::real_type *state){
  state[0] = 100.0;
  state[50] = 50.0;
  auto q = g_tchem->getStateVector(); //&state);
  for (ordinal_type k = 0; k < 67; k++){ 
  printf("in TChem_getStateVector %e\n", q[k]);
  state[k] = q[k];
  }
}

void TChem_setStateVector(double *array){
  g_tchem->setStateVector();
  for (ordinal_type k = 0; k < 67; k++){
     printf("in setStateVector %e\n", array[k]);
  }
}

void TChem::Driver::setStateVector() { //TChem::real_type *view){
  printf("in TChem::Driver::setStateVector\n");
  for (ordinal_type k = 0; k < _kmcd.nSpec; k++){
     _state(0,k) = 1.0;
  }
}

void TChem_getSpeciesNames(){
  g_tchem->getSpeciesNames();
}

void TChem::Driver::getSpeciesNames(){
  std::string species_names[_kmcd.nSpec];
  const auto speciesNamesHost = Kokkos::create_mirror_view(_kmcd.speciesNames);
  Kokkos::deep_copy(speciesNamesHost, _kmcd.speciesNames);
  for (ordinal_type k = 0; k < _kmcd.nSpec; k++)
  {
     species_names[k] = &speciesNamesHost(k,0);
     printf("%s \n", species_names[k].c_str());
  }
}

ordinal_type TChem::Driver::getNumberOfSpecies() { return _kmcd.nSpec;}

// FIXME: add error checking
ordinal_type TChem_getNumberOfSpecies(){
   ordinal_type nSpec = g_tchem->getNumberOfSpecies();
   return nSpec;
}

int TChem_getLengthOfStateVector() { return g_tchem == nullptr ? -1 : g_tchem->getLengthOfStateVector(); }

ordinal_type TChem::Driver::getLengthOfStateVector() const {
//  TCHEM_CHECK_ERROR(!_is_gasphase_kmcd_created, "const Kinetic model first needs to be created");
  return Impl::getStateVectorSize(_kmcd.nSpec);
}
