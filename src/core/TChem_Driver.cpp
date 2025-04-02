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
void initialize(const char* chemFile, const char* aeroFile, const char* numericsFile,
  const ordinal_type nBatch){

  g_tchem = new TChem::Driver();

  Kokkos::InitializationSettings settings;
//  settings.set_num_threads(8);
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

  g_tchem->setBatchSize(nBatch);
  g_tchem->createGasKineticModel(chemFile);
  g_tchem->createGasKineticModelConstData();
  g_tchem->createAerosolModel(aeroFile);
  g_tchem->createAerosolModelConstData();
  g_tchem->createStateVector(nBatch);
  g_tchem->getLengthOfStateVector();
  g_tchem->createNumerics(numericsFile);
  g_tchem->createNumberConcentrationVector(nBatch);

}

/* Finalize the model by freeing memory and finalizing Kokkos */
void finalize(){
  g_tchem->freeAll();
  delete g_tchem;
  Kokkos::finalize();
}

void TChem::Driver::setBatchSize(const ordinal_type nBatch) {
   _nBatch = nBatch;
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
   auto max_num_time_iterations = root["solver_info"]["max_num_time_iterations"];
   auto num_time_iterations_per_interval = root["solver_info"]["num_time_iterations_per_interval"];
   auto jacobian_interval = root["solver_info"]["jacobian_interval"];

   auto team_size = root["solver_info"]["team_size"];
   auto vector_size = root["solver_info"]["vector_size"];


   _atol_newton = atol_newton.as<real_type>(1e-10);
   _rtol_newton = rtol_newton.as<real_type>(1e-6);
   _dtmin = dtmin.as<real_type>(1e-8);
   _atol_time = atol_time.as<real_type>(1e-12);
   _rtol_time = rtol_time.as<real_type>(1e-4);
   _max_num_newton_iterations = max_num_newton_iterations.as<ordinal_type>(100);
   _max_num_time_iterations = max_num_time_iterations.as<ordinal_type>(1e3);
   _num_time_iterations_per_interval = num_time_iterations_per_interval.as<ordinal_type>(1e1);
   _jacobian_interval = jacobian_interval.as<ordinal_type>(1);

   // If team_size and vector_size are not specified, default to -1
   _team_size = team_size.as<ordinal_type>(-1);
   _vector_size = vector_size.as<ordinal_type>(-1);

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
  _aero_file = std::string();
}

/* Create the aerosol model from a YAML file */
void TChem::Driver::createAerosolModel(const std::string &aero_file) {
  _aero_file = aero_file;
  _amd = AerosolModelData(_aero_file, _kmd);
}

/* Create the aerosol constant model data */
void TChem::Driver::createAerosolModelConstData() {
  printf("Creating amcd \n");
  using interf_host_device_type = typename Tines::UseThisDevice<TChem::host_exec_space>::type;
  _amcd_host = TChem::create_AerosolModelConstData<interf_host_device_type>(_amd);
  _amcd_device = TChem::create_AerosolModelConstData<device_type>(_amd);
  printf("Number of aerosol species %d \n", _amcd_host.nSpec);
  printf("Maximum number of particles %d \n", _amcd_host.nParticles);
  printf("End creating amcd \n");
}

void TChem::Driver::createStateVector(ordinal_type nBatch) {
  const ordinal_type len = TChem::Impl::getStateVectorSize(_kmcd_host.nSpec + _amcd_host.nSpec * _amcd_host.nParticles);
  printf("Length %d %d %d \n", len, _amcd_host.nSpec, _amcd_host.nParticles);
  _state = real_type_2d_view_host("state dev", nBatch, len);
}

/* Create number concentration vector */
void TChem::Driver::createNumberConcentrationVector(ordinal_type nBatch) {
  const ordinal_type len = _amcd_host.nParticles;
  _number_concentration = real_type_2d_view_host("number_concentration", nBatch, len);
}

/* Set the values of the state vector */
void TChem_setNumberConcentrationVector(double *array, const ordinal_type iBatch){
  g_tchem->setNumberConcentrationVector(array, iBatch);
}

void TChem::Driver::setNumberConcentrationVector(double *array, const ordinal_type iBatch) {
  auto len = _amcd_host.nParticles;
  for (ordinal_type k = 0; k < len; k++){
     _number_concentration(iBatch,k) = array[k];
  }
}

/* Get the state vector */
auto TChem::Driver::getStateVector(const ordinal_type iBatch) {
  auto state_at_i_batch = Kokkos::subview(_state, iBatch, Kokkos::ALL);
  return state_at_i_batch;
}

void TChem_getStateVector(TChem::real_type *state, const ordinal_type iBatch){
  auto q = g_tchem->getStateVector(iBatch);
  auto len = TChem_getLengthOfStateVector();
  for (ordinal_type k = 0; k < len; k++) {
     state[k] = q[k];
  }
}

/* Set the values of the state vector */
void TChem_setStateVector(double *array, const ordinal_type iBatch){
  g_tchem->setStateVector(array, iBatch);
}

void TChem::Driver::setStateVector(double *array, const ordinal_type iBatch) {
  auto len = TChem_getLengthOfStateVector();
  for (ordinal_type k = 0; k < len; k++){
     _state(iBatch,k) = array[k];
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
  int k = *index;
  std::string species_name = &_kmcd_host.speciesNames(k,0);
  return species_name;
}

/* Return number of species */
ordinal_type TChem_getNumberOfSpecies(){
   ordinal_type nSpec = g_tchem->getNumberOfSpecies();
   return nSpec;
}

ordinal_type TChem::Driver::getNumberOfSpecies() { return _kmcd_host.nSpec;}

ordinal_type TChem_getNumberOfAeroSpecies(){
   ordinal_type nSpec = g_tchem->getNumberOfAeroSpecies();
   return nSpec;
}

ordinal_type TChem::Driver::getNumberOfAeroSpecies() { return _amcd_host.nSpec;}

/* Get the density of the aerosol species at a given index */
real_type TChem_getAerosolSpeciesDensity(int *index) {
   auto density = g_tchem->getAerosolSpeciesDensity(index);
   return density;
}

real_type TChem::Driver::getAerosolSpeciesDensity(int *index) {
   printf("Species index %d \n", *index);
   return (_amcd_host.aerosol_density(*index));
}

/* Get the molecular weight of the aerosol species at a given index */
real_type TChem_getAerosolSpeciesMW(int *index) {
   auto mw = g_tchem->getAerosolSpeciesMW(index);
   return mw;
}

real_type TChem::Driver::getAerosolSpeciesMW(int *index) {
   printf("Species index %d \n", *index);
   return (_amcd_host.molecular_weights(*index));
}

/* Get the hygroscopicity parameter kappa of the aerosol species at a given index */
real_type TChem_getAerosolSpeciesKappa(int *index) {
   auto kappa = g_tchem->getAerosolSpeciesKappa(index);
   return kappa;
}

real_type TChem::Driver::getAerosolSpeciesKappa(int *index) {
   printf("Species index %d \n", *index);
   return (_amcd_host.aerosol_kappa(*index));
}

/* Return species name */
int TChem_getAerosolSpeciesName(int * index, char* result, const std::string::size_type buffer_size){
  std::string specName = g_tchem->getAerosolSpeciesName(index);
  specName.copy(result, buffer_size);
  result[specName.length()] = '\0';
  return specName.length();
}

std::string TChem::Driver::getAerosolSpeciesName(int *index){


  std::map<int, std::string> aero_idx_sp_name;
  for (std::map<std::string, int>::iterator
        i = _amd.aerosol_sp_name_idx_.begin();
        i != _amd.aerosol_sp_name_idx_.end(); ++i)
        aero_idx_sp_name[i->second] = i->first;

  std::string species_name = aero_idx_sp_name[*index];

//  const auto speciesNamesHost = Kokkos::create_mirror_view(_kmcd_host.speciesNames);
//  int k = *index;
//  std::string species_name = &_kmcd_host.speciesNames(k,0);
  return species_name;
}

/* Return length of the state vector */
int TChem_getLengthOfStateVector() { 
  return g_tchem == nullptr ? -1 : g_tchem->getLengthOfStateVector();
}

ordinal_type TChem::Driver::getLengthOfStateVector() const {
  return Impl::getStateVectorSize(_kmcd_host.nSpec + _amcd_host.nSpec * _amcd_host.nParticles);
}

int TChem_getNumberConcentrationVectorSize() {
  return g_tchem == nullptr ? -1 : g_tchem->getNumberConcentrationVectorSize();
}

ordinal_type TChem::Driver::getNumberConcentrationVectorSize() const{
  return _amcd_host.nParticles;
}

/* Integrate a time step */
void TChem_doTimestep(const double &del_t){
//  g_tchem->doTimestep(del_t);
  g_tchem->doTimestep_sparse(del_t);
}

/* Integrate a time step */
void TChem::Driver::doTimestep(const double del_t){
  const auto exec_space_instance = TChem::host_exec_space();
  using device_type = typename Tines::UseThisDevice<TChem::exec_space>::type;
  using interf_host_device_type = typename Tines::UseThisDevice<TChem::host_exec_space>::type;
  using problem_type = TChem::Impl::AerosolChemistry_Problem<real_type, interf_host_device_type>;
  using policy_type = typename TChem::UseThisTeamPolicy<TChem::exec_space>::type;

  policy_type policy(exec_space_instance, _nBatch, Kokkos::AUTO());

  if (_team_size > 0 && _vector_size > 0) {
      policy = policy_type(exec_space_instance,  _nBatch, _team_size, _vector_size);
  } else if (_team_size > 0 && _vector_size < 0) {
      policy = policy_type(exec_space_instance, _nBatch,  _team_size);
  }

  using time_integrator_cvode_type = Tines::TimeIntegratorCVODE<real_type,host_device_type>;
  Tines::value_type_1d_view<time_integrator_cvode_type,interf_host_device_type> cvodes;
  cvodes = Tines::value_type_1d_view<time_integrator_cvode_type,interf_host_device_type>("cvodes", _nBatch);

  const ordinal_type total_n_species = _kmcd_host.nSpec + _amcd_host.nSpec * _amcd_host.nParticles;
  const ordinal_type n_active_vars = total_n_species - _kmcd_host.nConstSpec;
  for (ordinal_type i=0;i<_nBatch;++i){
       cvodes(i).create(n_active_vars);
  }

  const ordinal_type level = 1;
  ordinal_type per_team_extent(0);

 // per_team_extent = TChem::AtmosphericChemistry::getWorkSpaceSize(_kmcd_device);
  per_team_extent = TChem::AerosolChemistry_CVODE::getWorkSpaceSize(_kmcd_host, _amcd_host);

  const ordinal_type per_team_scratch =
      TChem::Scratch<real_type_1d_view>::shmem_size(per_team_extent);
  policy.set_scratch_size(level, Kokkos::PerTeam(per_team_scratch));

 // auto number_of_equations = problem_type::getNumberOfTimeODEs(_kmcd_host);
  auto number_of_equations = problem_type::getNumberOfTimeODEs(_kmcd_host, _amcd_host);
  real_type_2d_view tol_time("tol time", number_of_equations, 2);
  real_type_1d_view tol_newton("tol newton", 2);
  real_type_2d_view fac("fac", _nBatch, number_of_equations);

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

  time_advance_type_1d_view tadv("tadv", _nBatch);
  Kokkos::deep_copy(tadv, tadv_default);
  real_type_1d_view t("time", _nBatch);
  Kokkos::deep_copy(t, tadv_default._tbeg);
  real_type_1d_view dt("delta time", _nBatch);
  Kokkos::deep_copy(dt, tadv_default._dt);

  real_type tend = del_t;
  ordinal_type iter = 0;
  real_type tsum(0);
  auto stateVecDim = TChem_getLengthOfStateVector();
  real_type_2d_view state("StateVector Devices", _nBatch, stateVecDim);
  Kokkos::deep_copy(state, _state);
  real_type_2d_view number_conc("NumberConcentration", _nBatch, _amcd_host.nParticles);
  Kokkos::deep_copy(number_conc, _number_concentration);

  const real_type zero(0);

  for (; iter < _max_num_time_iterations && tsum <= tend * 0.9999; ++iter) {
       TChem::AerosolChemistry_CVODE::runHostBatch(
              policy, tol_time, fac, tadv, state, number_conc, t, dt, state,
              _kmcd_host, _amcd_host, cvodes);

    tsum = zero;

    Kokkos::parallel_reduce(
        Kokkos::RangePolicy<TChem::exec_space>(0, _nBatch),
        KOKKOS_LAMBDA(const ordinal_type &i, real_type &update) {
          tadv(i)._tbeg = t(i);
          tadv(i)._dt = dt(i);
          update += t(i);
     }, tsum);
     Kokkos::fence();
  }
  Kokkos::deep_copy(_state, state);
}

/* Return the state vector from host */
void TChem::Driver::getStateVectorHost(real_type_2d_const_view_host &view) {
  TCHEM_CHECK_ERROR(_state.span() == 0, "State vector should be constructed");
  view = real_type_2d_const_view_host(&_state(0,0), _state.extent(0), _state.extent(1));
}

/* Integrate a time step */
void TChem::Driver::doTimestep_sparse(const double del_t){

  using VecType   = sundials::kokkos::Vector<TChem::exec_space>;
  using MatType   = sundials::kokkos::DenseMatrix<TChem::exec_space>;
  using LSType    = sundials::kokkos::DenseLinearSolver<TChem::exec_space>;
  using SizeType  = VecType::size_type;

  const auto exec_space_instance = TChem::exec_space();
  using device_type = typename Tines::UseThisDevice<TChem::exec_space>::type;
  using interf_host_device_type = typename Tines::UseThisDevice<TChem::host_exec_space>::type;
  using problem_type = TChem::Impl::AerosolChemistry_Problem<real_type, interf_host_device_type>;
  using policy_type = typename TChem::UseThisTeamPolicy<TChem::exec_space>::type;

  policy_type policy(exec_space_instance, _nBatch, Kokkos::AUTO());

  if (_team_size > 0 && _vector_size > 0) {
      policy = policy_type(exec_space_instance,  _nBatch, _team_size, _vector_size);
  } else if (_team_size > 0 && _vector_size < 0) {
      policy = policy_type(exec_space_instance, _nBatch,  _team_size);
  }

  auto number_of_equations = problem_type::getNumberOfTimeODEs(_kmcd_host, _amcd_host);

  // Create UserData
  TChem::UserData udata;

  udata.nbatches = _nBatch;
  udata.num_concentration = _number_concentration;
  udata.batchSize=number_of_equations;
  udata.kmcd=_kmcd_host;
  udata.amcd=_amcd_host;

  // Create the SUNDIALS context
  sundials::Context sunctx;
  using real_type_1d_view_type = Tines::value_type_1d_view<real_type, device_type>;
  using real_type_2d_view_type = Tines::value_type_2d_view<real_type, device_type>;
  using real_type_2d_view_host_type = Tines::value_type_2d_view<real_type, host_device_type>;

  // Create vector with the initial condition
  const sunrealtype T0 = SUN_RCONST(0.0);

  SizeType length{static_cast<SizeType>(_nBatch * number_of_equations)};
  VecType y{length, sunctx};
  real_type_2d_view_type y2d((y.View()).data(), _nBatch, number_of_equations);

  const ordinal_type n_active_gas_species = _kmcd_host.nSpec - _kmcd_host.nConstSpec;
  real_type_2d_view_type const_tracers("const_tracers", _nBatch, _kmcd_host.nConstSpec);

  real_type_1d_view_type temperature("temperature", _nBatch);
  real_type_1d_view_type pressure("pressure", _nBatch);
  using range_type = Kokkos::pair<ordinal_type, ordinal_type>;
  const ordinal_type level = 1;

   Kokkos::parallel_for(
      "fill_y", Kokkos::RangePolicy<TChem::exec_space>(0, _nBatch),
      KOKKOS_LAMBDA(const SizeType i) {
        const real_type_1d_view_type state_at_i =
               Kokkos::subview(_state, i, Kokkos::ALL());
        const ordinal_type total_n_species = _kmcd_host.nSpec + _amcd_host.nParticles*_amcd_host.nSpec;
        TChem::Impl::StateVector<real_type_1d_view_type> sv_at_i(total_n_species, state_at_i);
        temperature(i) = sv_at_i.Temperature();
        pressure(i) = sv_at_i.Pressure();
        const real_type_1d_view_type Ys = sv_at_i.MassFractions();

        const auto activeYs = Kokkos::subview(Ys, range_type(0, n_active_gas_species));
        const auto constYs  = Kokkos::subview(Ys, range_type(n_active_gas_species, _kmcd_host.nSpec));
        const real_type_1d_view_type partYs = Kokkos::subview(Ys, range_type(_kmcd_host.nSpec, total_n_species));

        for (ordinal_type j=0;j<n_active_gas_species;++j){
          y2d(i, j) = activeYs(j);
        }

        for (ordinal_type j=n_active_gas_species;j<total_n_species- _kmcd_host.nConstSpec;++j)
        {
          y2d(i, j) = partYs(j-n_active_gas_species);
        }

        for (ordinal_type j=0;j<_kmcd_host.nConstSpec;++j){
          const_tracers(i, j) = constYs(j);
        }

     });

/*  int i = 0;
  for (ordinal_type j=0;j<n_active_gas_species;++j){
    printf("%d %e :\n", j, y2d(i, j));
  }
*/
    udata.temperature = temperature;
    udata.pressure = pressure;
    udata.const_tracers = const_tracers;

    // Create vector of absolute tolerances
    VecType abstol{length, sunctx};
    N_VConst(SUN_RCONST(_atol_time), abstol);

    // Create CVODE using Backward Differentiation Formula methods
    void* cvode_mem = CVodeCreate(CV_BDF, sunctx);
//    if (check_ptr(cvode_mem, "CVodeCreate")) { return 1; }


    // Initialize the integrator and set the ODE right-hand side function
    int retval = CVodeInit(cvode_mem, TChem::AerosolChemistry_CVODE_K::f, T0, y);
//    if (check_flag(retval, "CVodeInit")) { return 1; }


    // Attach the user data structure
    retval = CVodeSetUserData(cvode_mem, &udata);
//    if (check_flag(retval, "CVodeSetUserData")) { return 1; }

    // Specify the scalar relative tolerance and vector absolute tolerances
    retval = CVodeSVtolerances(cvode_mem, SUN_RCONST(_rtol_time), abstol);
//    if (check_flag(retval, "CVodeSVtolerances")) { return 1; }

     // Create the matrix and linear solver objects
    std::unique_ptr<sundials::ConvertibleTo<SUNMatrix>> A;
    std::unique_ptr<sundials::ConvertibleTo<SUNLinearSolver>> LS;

    ordinal_type per_team_extent=0;

    // Create matrix-free GMRES linear solver
    LS = std::make_unique<sundials::experimental::SUNLinearSolverView>(
        SUNLinSol_SPGMR(y, SUN_PREC_NONE, 0, sunctx));

    // Attach the linear solver to CVODE
    retval = CVodeSetLinearSolver(cvode_mem, LS->Convert(), nullptr);
//    if (check_flag(retval, "CVodeSetLinearSolver")) { return 1; }
    per_team_extent
         = TChem::Impl::Aerosol_RHS<real_type, device_type>::getWorkSpaceSize(_kmcd_host, _amcd_host);

    const ordinal_type per_team_scratch =
      TChem::Scratch<real_type_1d_view_type>::shmem_size(per_team_extent);
      policy.set_scratch_size(level, Kokkos::PerTeam(per_team_scratch));
    udata.policy = policy;

    // Final time and time between outputs

    const sunrealtype Tf    = SUN_RCONST(del_t);
    const sunrealtype dTout = SUN_RCONST(del_t); //_dtmin);

    // Number of output times
    const int Nt_p = static_cast<int>(ceil(Tf / dTout));

    const int Nt =  _max_num_time_iterations > 0 ? _max_num_time_iterations: Nt_p;

    // Current time and first output time
    sunrealtype t    = T0;
    sunrealtype tout = T0 + dTout;

    // Initial output
    real_type_2d_view_host_type y2d_h((y.HostView()).data(), udata.nbatches, udata.batchSize);
    sundials::kokkos::CopyFromDevice(y);
    Kokkos::fence();

//    const auto density_host =
//    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace() , density);
    const auto temperature_host =
    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace() , temperature);
    const auto pressure_host =
    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), pressure);
    const auto const_tracers_host =
    Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),const_tracers);

  // Time stepping
    int iout = 0;
    for (iout = 0; iout < Nt; iout++)
    {
      // Advance in time
      retval = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
      exec_space_instance.fence();
//      if (check_flag(retval, "CVode")) { break; }

      // // Output solution from some batches
      sundials::kokkos::CopyFromDevice(y);
      Kokkos::fence();

      tout += dTout;
      tout = (tout > Tf) ? Tf : tout;
    }

  // Copy results back
  int i = 0;
  auto stateVecDim = TChem_getLengthOfStateVector();
  real_type_2d_view state("StateVector Devices", _nBatch, stateVecDim);
  for (ordinal_type j=0;j<n_active_gas_species;++j){
    _state(i,j+3) = y2d_h(i, j);
  }
  const ordinal_type total_n_species = _kmcd_host.nSpec + _amcd_host.nParticles*_amcd_host.nSpec;
  for (ordinal_type j=n_active_gas_species;j<total_n_species - _kmcd_host.nConstSpec;++j)
  {
    _state(i,j+3+_kmcd_host.nConstSpec) = y2d_h(i,j);
  }  
}
