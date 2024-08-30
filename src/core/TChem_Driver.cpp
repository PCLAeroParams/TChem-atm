#include "TChem_Driver.hpp"
#include "TChem.hpp"
#include "TChem_KineticModelNCAR_ConstData.hpp"
#include "TChem_CommandLineParser.hpp" 

using real_type = TChem::real_type;
using ordinal_type = TChem::ordinal_type;

using exec_space = Kokkos::DefaultExecutionSpace;
using host_exec_space = Kokkos::DefaultHostExecutionSpace;

using host_device_type = typename Tines::UseThisDevice<TChem::host_exec_space>::type;
using device_type = typename Tines::UseThisDevice<exec_space>::type;

static TChem::Driver *g_tchem = nullptr;

/* Initialize the model given input YAML files */
void initialize(const char* chemFile, const char* aeroFile, const char* numericsFile){

  g_tchem = new TChem::Driver();

  Kokkos::InitializationSettings settings;
  settings.set_num_threads(8);
  settings.set_device_id(0);
  Kokkos::initialize(settings);

  using exec_space = Kokkos::DefaultExecutionSpace;
  using host_exec_space = Kokkos::DefaultHostExecutionSpace;

  const bool detail = false;
  TChem::exec_space().print_configuration(std::cout, detail);
  TChem::host_exec_space().print_configuration(std::cout, detail);

  using host_device_type = typename Tines::UseThisDevice<TChem::host_exec_space>::type;
  using device_type = typename Tines::UseThisDevice<exec_space>::type;

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
  g_tchem->getLengthOfStateVector();
  g_tchem->createNumerics(numericsFile);

}

/* Finalize the model by freeing memory and finalizing Kokkos */
void finalize(){
  g_tchem->freeAll();
  delete g_tchem;
  Kokkos::finalize();
}

/* Read in the solver settings */
void TChem::Driver::createNumerics(const std::string &numerics_file) {

   YAML::Node root = YAML::LoadFile(numerics_file);

   auto atol_newton = root["solver_info"]["atol_newton"];
   auto rtol_newton =  root["solver_info"]["rtol_newton"];
   auto dtmin = root["solver_info"]["dtmin"];
   auto atol_time = root["solver_info"]["atol_time"];
   auto rtol_time = root["solver_info"]["rtol_time"];
   auto max_num_newton_iterations = root["solver_info"]["max_newton_iterations"];
   auto num_time_iterations_per_interval = root["solver_info"]["max_time_iterations"];
   auto jacobian_interval = root["solver_info"]["jacobian_interval"];

   _atol_newton = atol_newton.as<real_type>(1e-10);
   _rtol_newton = rtol_newton.as<real_type>(1e-6);
   _dtmin = dtmin.as<real_type>(1e-8);
   _atol_time = atol_time.as<real_type>(1e-12);
   _rtol_time = rtol_time.as<real_type>(1e-4);
   _max_num_newton_iterations = max_num_newton_iterations.as<ordinal_type>(100);
   _num_time_iterations_per_interval = num_time_iterations_per_interval.as<ordinal_type>(1e3);
   _jacobian_interval = jacobian_interval.as<ordinal_type>(1);

}

/* Create the gas kinetic model from a YAML file */
void TChem::Driver::createGasKineticModel(const std::string &chem_file) {
  _chem_file = chem_file;
  _kmd = KineticModelData(_chem_file);
}

/* Create the gas kinetic constant model data */
void TChem::Driver::createGasKineticModelConstData() {
  printf("Creating kmcd \n");
  using interf_host_device_type = typename Tines::UseThisDevice<TChem::host_exec_space>::type;
  _kmcd_host = TChem::createNCAR_KineticModelConstData<interf_host_device_type>(_kmd);
  _kmcd_device = TChem::createNCAR_KineticModelConstData<device_type>(_kmd);
  printf("Number of Species %d \n", _kmcd_host.nSpec);
  printf("Number of Reactions %d \n", _kmcd_host.nReac);

  printf("End creating kmcd \n");

}

/* Free gas kinetic model */
void TChem::Driver::freeAll() {
  g_tchem->freeGasKineticModel();
}

void TChem::Driver::freeGasKineticModel() {
  _chem_file = std::string();
  _therm_file = std::string();
}

void TChem::Driver::createStateVector() {
  const ordinal_type len = TChem::Impl::getStateVectorSize(_kmcd_host.nSpec);
  const ordinal_type nBatch = 1;
  _state = real_type_2d_view_host("state dev", nBatch, len);
  for (ordinal_type k = 0; k < len; k++) {
    _state(0,k) = 0.0;
  }
}

/* Get the state vector */
void TChem_getStateVector(TChem::real_type *state){
  auto q = g_tchem->getStateVector();
  auto len = TChem_getLengthOfStateVector();
  for (ordinal_type k = 0; k < len; k++) {
     state[k] = q[k];
  }
}

auto TChem::Driver::getStateVector() {
  auto state_at_0 = Kokkos::subview(_state, 0, Kokkos::ALL);
  return state_at_0;
}

/* Set the values of the state vector */
void TChem_setStateVector(double *array){
  g_tchem->setStateVector(array);
}

void TChem::Driver::setStateVector(double *array) {
  auto len = TChem_getLengthOfStateVector();
  for (ordinal_type k = 0; k < len; k++){
     _state(0,k) = array[k];
  }
}

/* Return species name */
int TChem_getSpeciesName(int * index, char* result, const std::string::size_type buffer_size){
  std::string specName = g_tchem->getSpeciesName(index);  
  specName.copy(result, buffer_size);
  result[specName.length()] = '\0';
  return specName.length();
}

std::string TChem::Driver::getSpeciesName(int *index){
  const auto speciesNamesHost = Kokkos::create_mirror_view(_kmcd_host.speciesNames);
  Kokkos::deep_copy(speciesNamesHost, _kmcd_host.speciesNames);
  int k = *index;
  std::string species_name = &speciesNamesHost(k,0);
  return species_name;
}

/* Return number of species */
ordinal_type TChem_getNumberOfSpecies(){
   ordinal_type nSpec = g_tchem->getNumberOfSpecies();
   return nSpec;
}

ordinal_type TChem::Driver::getNumberOfSpecies() { return _kmcd_host.nSpec;}

/* Return length of the state vector */
int TChem_getLengthOfStateVector() { 
  return g_tchem == nullptr ? -1 : g_tchem->getLengthOfStateVector();
}

ordinal_type TChem::Driver::getLengthOfStateVector() const {
  return Impl::getStateVectorSize(_kmcd_host.nSpec);
}

/* Integrate a time step */
void TChem_doTimestep(const double &del_t){
  g_tchem->doTimestep(del_t);
}

/* Integrate a time step */
void TChem::Driver::doTimestep(const double del_t){

  const auto exec_space_instance = TChem::exec_space();
//  using device_type = typename Tines::UseThisDevice<exec_space>::type;
  using device_type = typename Tines::UseThisDevice<TChem::exec_space>::type;
  using interf_host_device_type = typename Tines::UseThisDevice<TChem::host_exec_space>::type;
  using problem_type = TChem::Impl::AtmosphericChemistry_Problem<real_type, interf_host_device_type>;
  using policy_type = typename TChem::UseThisTeamPolicy<TChem::exec_space>::type;
  policy_type policy(exec_space_instance, 1, Kokkos::AUTO());

  const ordinal_type level = 1;
  ordinal_type per_team_extent(0);

  per_team_extent = TChem::AtmosphericChemistry::getWorkSpaceSize(_kmcd_host);

  const ordinal_type per_team_scratch =
      TChem::Scratch<real_type_1d_view>::shmem_size(per_team_extent);
  policy.set_scratch_size(level, Kokkos::PerTeam(per_team_scratch));

  auto number_of_equations = problem_type::getNumberOfTimeODEs(_kmcd_host);
  real_type_2d_view tol_time("tol time", number_of_equations, 2);
  real_type_1d_view tol_newton("tol newton", 2);
  real_type_2d_view fac("fac", 1, number_of_equations);

  {
      auto tol_time_host = Kokkos::create_mirror_view(tol_time);
      auto tol_newton_host = Kokkos::create_mirror_view(tol_newton);
      for (ordinal_type i = 0, iend = tol_time.extent(0); i < iend; ++i) {
          tol_time_host(i, 0) = _atol_time;
          tol_time_host(i, 1) = _rtol_time;
      }
      tol_newton_host(0) = _atol_newton;
      tol_newton_host(1) = _rtol_newton;
      
      Kokkos::deep_copy(tol_time, tol_time_host);
      Kokkos::deep_copy(tol_newton, tol_newton_host);
  } 

  using time_advance_type = TChem::time_advance_type;
  time_advance_type tadv_default;
  tadv_default._tbeg = 0.0;
  tadv_default._tend = del_t;
  tadv_default._dt = _dtmin;
  tadv_default._dtmin = _dtmin;
  tadv_default._dtmax = del_t;
  tadv_default._max_num_newton_iterations = _max_num_newton_iterations;
  tadv_default._num_time_iterations_per_interval = _num_time_iterations_per_interval;
  tadv_default._jacobian_interval = _jacobian_interval;

  time_advance_type_1d_view tadv("tadv", 1);
  Kokkos::deep_copy(tadv, tadv_default);
  real_type_1d_view t("time", 1);
  Kokkos::deep_copy(t, tadv_default._tbeg);
  real_type_1d_view dt("delta time", 1);
  Kokkos::deep_copy(dt, tadvx_default._dt);

  int max_num_time_iterations = 1e5;
  real_type tend = del_t;
  ordinal_type iter = 0;
  real_type tsum(0);
  auto stateVecDim = TChem_getLengthOfStateVector();
  real_type_2d_view state("StateVector Devices", 1, stateVecDim);
  Kokkos::deep_copy(state, _state);
  for (; iter < max_num_time_iterations && tsum <= tend * 0.9999; ++iter) {
    TChem::AtmosphericChemistry::runDeviceBatch(policy, tol_newton, tol_time, fac, tadv,
              state, t, dt, state, _kmcd_device);

    tsum = 0.0;
    Kokkos::parallel_reduce(
        Kokkos::RangePolicy<TChem::exec_space>(0, 1),
        KOKKOS_LAMBDA(const ordinal_type &i, real_type &update) {
          tadv(i)._tbeg = t(i);
          tadv(i)._dt = dt(i);
          update += t(i);
     }, tsum);
  }
  Kokkos::deep_copy(_state, state);
}

/* Return the state vector from host */
void TChem::Driver::getStateVectorHost(real_type_2d_const_view_host &view) {
  TCHEM_CHECK_ERROR(_state.span() == 0, "State vector should be constructed");
  auto hv = Kokkos::create_mirror_view(_state);
  view = real_type_2d_const_view_host(&hv(0, 0), hv.extent(0), hv.extent(1));
} 
