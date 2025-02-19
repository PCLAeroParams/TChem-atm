/* =====================================================================================
TChem-atm version 1.0
Copyright (2024) NTESS
https://github.com/sandialabs/TChem-atm

Copyright 2024 National Technology & Engineering Solutions of Sandia, LLC
(NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
Government retains certain rights in this software.

This file is part of TChem-atm. TChem-atm is open source software: you can redistribute
it and/or modify it under the terms of BSD 2-Clause License
(https://opensource.org/licenses/BSD-2-Clause). A copy of the licese is also
provided under the main directory

Questions? Contact Oscar Diaz-Ibarra at <odiazib@sandia.gov>, or
           Mike Schmidt at <mjschm@sandia.gov>, or
           Cosmin Safta at <csafta@sandia.gov>

Sandia National Laboratories, New Mexico/Livermore, NM/CA, USA
=====================================================================================
*/
#include "TChem.hpp"
#include "TChem_CommandLineParser.hpp"
#include "TChem_Impl_AerosolChemistry.hpp"

#include <cstdio>
#include <cvode/cvode.h>
#include <memory>
#include <nvector/nvector_kokkos.hpp>
#include <sundials/sundials_core.hpp>
#include <sunlinsol/sunlinsol_kokkosdense.hpp>
#include <sunlinsol/sunlinsol_spgmr.h>
#include <sunmatrix/sunmatrix_kokkosdense.hpp>
#include <vector>

using ordinal_type = TChem::ordinal_type;
using real_type = TChem::real_type;
// User data structure available in user-supplied callback functions
using policy_type =
          typename TChem::UseThisTeamPolicy<TChem::exec_space>::type;
using device_type = typename Tines::UseThisDevice<TChem::exec_space>::type;
using problem_type = TChem::Impl::AerosolChemistry_Problem<real_type,device_type>;
using host_device_type = typename Tines::UseThisDevice<TChem::host_exec_space>::type;
using real_type_0d_view_type = Tines::value_type_0d_view<real_type, device_type>;
using real_type_1d_view_type = Tines::value_type_1d_view<real_type, device_type>;
using real_type_2d_view_type = Tines::value_type_2d_view<real_type, device_type>;
using real_type_3d_view_type = Tines::value_type_3d_view<real_type, device_type>;
using real_type_2d_view_host_type = Tines::value_type_2d_view<real_type, host_device_type>;

using VecType   = sundials::kokkos::Vector<TChem::exec_space>;
using MatType   = sundials::kokkos::DenseMatrix<TChem::exec_space>;
using LSType    = sundials::kokkos::DenseLinearSolver<TChem::exec_space>;
using SizeType  = VecType::size_type;

// Check for an unrecoverable (negative) return value from a SUNDIALS function
static int check_flag(const int flag, const std::string funcname)
{
  if (flag < 0)
  {
    std::cerr << "ERROR: " << funcname << " returned " << flag << std::endl;
    return 1;
  }
  return 0;
}

// Check if a function returned a NULL pointer
static int check_ptr(const void* ptr, const std::string funcname)
{
  if (ptr) { return 0; }
  std::cerr << "ERROR: " << funcname << " returned NULL" << std::endl;
  return 1;
}

struct UserData
{
  int nbatches  = 100; // number of chemical networks
  int batchSize = 3;   // size of each network
  policy_type policy;
  real_type_2d_view_type num_concentration;
  real_type_1d_view_type temperature;
  real_type_1d_view_type pressure;
  real_type_2d_view_type const_tracers;
  real_type_2d_view_type fac;
  TChem::KineticModelNCAR_ConstData<device_type> kmcd;
  TChem::AerosolModel_ConstData<device_type> amcd;
};

// User-supplied functions called by CVODE
static int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data);

static int Jac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J,
               void* user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);


int main(int argc, char *argv[]) {
  std::string chemFile( "chem.yaml");
  std::string aeroFile("aero.yaml");
  std::string outputFile("AerosolChemistry.dat");
  std::string inputFile("None");
  std::string input_file_particles("None");
  std::string outputFileTimes("wall_times.json");
  int nBatch(1), team_size(-1), vector_size(-1);
  bool verbose(true);
  bool use_cloned_samples(false);
  int number_of_particles(-1);
  real_type rtol_time(1e-4), atol_time(1e-12);
  real_type tbeg(0), tend(1);
  real_type dtmin(1e-8), dtmax(1e-1);

  /// parse command line arguments
  TChem::CommandLineParser opts("This example computes the net production rates with a given state vector");
  opts.set_option<std::string>("chemfile", "Chem file name e.g., chem.inp",
                               &chemFile);
  opts.set_option<std::string>("aerofile", "Aerosol Chem file name e.g., aero.yaml",
                               &aeroFile);
  opts.set_option<std::string>("inputfile", "Chem file name e.g., chem.inp",
                               &inputFile);
  opts.set_option<std::string>("inputfile_particles", "input file name e.g., chem.inp",
                               &input_file_particles);
  opts.set_option<std::string>("outputfile", "Output omega file name e.g., omega.dat", &outputFile);
  opts.set_option<std::string>("outputfile_times", "Wal times file name e.g., times.json", &outputFileTimes);
  opts.set_option<int>("team_thread_size", "time thread size ", &team_size);
  opts.set_option<int>("vector_thread_size", "vector thread size ", &vector_size);
  opts.set_option<int>("number_of_particles", "Set the number of particles; this will overwrite the value from the input file. ", &number_of_particles);
  opts.set_option<bool>("verbose", "If true, printout the first omega values", &verbose);
  opts.set_option<int>("batch_size", " number of batches or samples, e.g. 10  ", &nBatch);
  opts.set_option<bool>(
      "use_cloned_samples", "If true, one state vector will be cloned.", &use_cloned_samples);
  opts.set_option<real_type>(
      "rtol-time", "Tolerence used for adaptive time stepping", &rtol_time);
  opts.set_option<real_type>(
      "atol-time", "Tolerence used for adaptive time stepping", &atol_time);
  opts.set_option<real_type>("tbeg", "Time begin", &tbeg);
  opts.set_option<real_type>("tend", "Time end", &tend);
  opts.set_option<real_type>("dtmin", "Minimum time step size", &dtmin);

  const bool r_parse = opts.parse(argc, argv);
  if (r_parse)
    return 0; // print help return
  // Create the SUNDIALS context
  sundials::Context sunctx;
  Kokkos::initialize(argc, argv);
  {

    // by default use input condition from chemFile
    if (inputFile == "None") {
       inputFile=chemFile;
    }
    // check first in chem file
    if (input_file_particles == "None") {
       input_file_particles=chemFile;
    }

  	const bool detail = false;
    constexpr real_type zero =0;

    TChem::exec_space().print_configuration(std::cout, detail);
    TChem::host_exec_space().print_configuration(std::cout, detail);
    const auto exec_space_instance = TChem::exec_space();

    /// construct kmd and use the view for testing
    printf("kmd parsing ...\n");
    TChem::KineticModelData kmd(chemFile);
    const auto kmcd = TChem::createNCAR_KineticModelConstData<device_type>(kmd);

    printf("amd parsing ...\n");
    TChem::AerosolModelData amd(aeroFile, kmd);
    // if(number_of_particles > 0) {
    //   amd.setNumberofParticles(number_of_particles);
    // }
    const auto amcd = TChem::create_AerosolModelConstData<device_type>(amd);

    const ordinal_type total_n_species =kmcd.nSpec + amcd.nSpec * amcd.nParticles;
    const ordinal_type number_of_equations = problem_type::getNumberOfTimeODEs(kmcd, amcd);
    const ordinal_type stateVecDim =
        TChem::Impl::getStateVectorSize(total_n_species);
    const auto speciesNamesHost = Kokkos::create_mirror_view(kmcd.speciesNames);
    Kokkos::deep_copy(speciesNamesHost, kmcd.speciesNames);

    printf("Number of Gas Species %d \n", kmcd.nSpec);
    printf("Number of Gas Reactions %d \n", kmcd.nReac);

    printf("Number of Aerosol Species %d \n", amcd.nSpec);
    printf("Number of Aerosol Particles %d \n", amcd.nParticles);

    printf("stateVecDim %d \n", stateVecDim);

 // read scenario condition from yaml file
    real_type_2d_view_host_type state_scenario_host;
    ordinal_type nbatch_files=1;
    printf("conditions parsing ...\n");
    TChem::AtmChemistry::setScenarioConditions(inputFile,
     speciesNamesHost, kmcd.nSpec, stateVecDim, state_scenario_host, nbatch_files);
    printf("Number of const species %d \n", kmcd.nConstSpec);

    real_type_2d_view_host_type num_concentration_host;
    real_type_2d_view_host_type state_host;

    // scenario particles
    amd.scenarioConditionParticles(input_file_particles,
                                   nbatch_files,
                                   num_concentration_host,
                                   state_scenario_host);

    real_type_2d_view_type num_concentration;
    real_type_2d_view_type state;

    if (nbatch_files == 1 && use_cloned_samples && nBatch > 1) {
      // only clone samples if nbatch_files is 1
      printf("-------------------------------------------------------\n");
      printf("--------------------Warning----------------------------\n");
      printf("Using cloned samples ... only for numerical experiments\n");
      state = real_type_2d_view_type("StateVector Devices", nBatch, stateVecDim);
      state_host = Kokkos::create_mirror_view(state);
      auto state_scenario_host_at_0 = Kokkos::subview(state_scenario_host, 0, Kokkos::ALL);
      auto state_at_0 = Kokkos::subview(state, 0, Kokkos::ALL);
      Kokkos::deep_copy(state_at_0, state_scenario_host_at_0);
      TChem::Test::cloneView(state);
      //make sure to copy state to host view.
      Kokkos::deep_copy(state_host, state);
      //scenario particles
      auto num_concentration_host_at_0 = Kokkos::subview(num_concentration_host, 0, Kokkos::ALL);
      num_concentration = real_type_2d_view_type("num_concentration", nBatch, amd.nParticles_);
      auto num_concentration_at_0 = Kokkos::subview(num_concentration, 0, Kokkos::ALL);
      Kokkos::deep_copy(num_concentration_at_0, num_concentration_host_at_0);
      TChem::Test::cloneView(num_concentration);

    } else {
      nBatch = nbatch_files;
      state = real_type_2d_view_type("StateVector Devices", nBatch, stateVecDim);
      Kokkos::deep_copy(state, state_scenario_host);
      state_host = state_scenario_host;

      // scenario particles
      num_concentration = real_type_2d_view_type("num_concentration", nBatch, amd.nParticles_);
      Kokkos::deep_copy(num_concentration, num_concentration_host);
    }
    // output
    real_type_2d_view_type fac;

    using policy_type = typename TChem::UseThisTeamPolicy<TChem::exec_space>::type;

    // fprintf(fout_times, " Computation type, total time [sec], time per sample [sec/sample]\n ");

    policy_type policy(exec_space_instance, nBatch, Kokkos::AUTO());

    if (team_size > 0 && vector_size > 0) {
        policy = policy_type(exec_space_instance,  nBatch, team_size, vector_size);
    } else if (team_size > 0 && vector_size < 0) {
      // only set team size
       policy = policy_type(exec_space_instance, nBatch,  team_size);
    }

    using range_type = Kokkos::pair<ordinal_type, ordinal_type>;
    const ordinal_type level = 1;
    const ordinal_type per_team_extent = problem_type::getWorkSpaceSize(kmcd,amcd)
    + number_of_equations;
    const ordinal_type per_team_scratch =
    TChem::Scratch<real_type_1d_view_type>::shmem_size(per_team_extent);
    policy.set_scratch_size(level, Kokkos::PerTeam(per_team_scratch));

    // Create UserData
    UserData udata;

    udata.nbatches = nBatch;
    udata.policy = policy;
    udata.num_concentration = num_concentration;
    udata.batchSize=number_of_equations;
    udata.kmcd=kmcd;
    udata.amcd=amcd;
    // Create vector with the initial condition
    const sunrealtype T0 = SUN_RCONST(tbeg);

    SizeType length{static_cast<SizeType>(nBatch * number_of_equations)};
    VecType y{length, sunctx};
    real_type_2d_view_type y2d((y.View()).data(), nBatch, number_of_equations);

    const ordinal_type n_active_gas_species = kmcd.nSpec - kmcd.nConstSpec;
    real_type_2d_view_type const_tracers("const_tracers",nBatch, n_active_gas_species);


    real_type_1d_view_type temperature("temperature",nBatch);
    real_type_1d_view_type pressure("pressure",nBatch);

    Kokkos::parallel_for(
      "fill_y", Kokkos::RangePolicy<TChem::exec_space>(0, nBatch),
      KOKKOS_LAMBDA(const SizeType i) {
        const real_type_1d_view_type state_at_i =
        Kokkos::subview(state, i, Kokkos::ALL());
        const ordinal_type total_n_species = kmcd.nSpec + amcd.nParticles*amcd.nSpec;
        TChem::Impl::StateVector<real_type_1d_view_type> sv_at_i(total_n_species, state_at_i);
        temperature(i) = sv_at_i.Temperature();
        pressure(i) = sv_at_i.Pressure();
        const real_type_1d_view_type Ys = sv_at_i.MassFractions();

        const auto activeYs = Kokkos::subview(Ys, range_type(0, n_active_gas_species));
        const auto constYs  = Kokkos::subview(Ys, range_type(n_active_gas_species, kmcd.nSpec));
        const real_type_1d_view_type partYs = Kokkos::subview(Ys, range_type(kmcd.nSpec, total_n_species));

        for (ordinal_type j=0;j<n_active_gas_species;++j){
          y2d(i, j) = activeYs(j);
        }

        for (ordinal_type j=n_active_gas_species;j<total_n_species- kmcd.nConstSpec;++j)
        {
          y2d(i, j) = partYs(j-n_active_gas_species);
        }

        for (ordinal_type j=0;j<kmcd.nConstSpec;++j){
          const_tracers(i, j) = constYs(j);
        }

     });

    udata.temperature = temperature;
    udata.pressure = pressure;
    udata.const_tracers =const_tracers;

    // Create vector of absolute tolerances
    VecType abstol{length, sunctx};
    N_VConst(SUN_RCONST(atol_time), abstol);

    // Create CVODE using Backward Differentiation Formula methods
    void* cvode_mem = CVodeCreate(CV_BDF, sunctx);
    if (check_ptr(cvode_mem, "CVodeCreate")) { return 1; }


    // Initialize the integrator and set the ODE right-hand side function
    int retval = CVodeInit(cvode_mem, f, T0, y);
    if (check_flag(retval, "CVodeInit")) { return 1; }


    // Attach the user data structure
    retval = CVodeSetUserData(cvode_mem, &udata);
    if (check_flag(retval, "CVodeSetUserData")) { return 1; }

    // Specify the scalar relative tolerance and vector absolute tolerances
    retval = CVodeSVtolerances(cvode_mem, SUN_RCONST(rtol_time), abstol);
    if (check_flag(retval, "CVodeSVtolerances")) { return 1; }

     // Create the matrix and linear solver objects
    std::unique_ptr<sundials::ConvertibleTo<SUNMatrix>> A;
    std::unique_ptr<sundials::ConvertibleTo<SUNLinearSolver>> LS;


    // Linear solver type
    int solver_type = 0;


    if (solver_type == 0)
    {
      // Create Kokkos dense block diagonal matrix
      A = std::make_unique<MatType>(udata.nbatches, udata.batchSize, udata.batchSize, sunctx);

      // Create Kokkos batched dense linear solver
      LS = std::make_unique<LSType>(sunctx);

      // Attach the matrix and linear solver to CVODE
      retval = CVodeSetLinearSolver(cvode_mem, LS->Convert(), A->Convert());
      if (check_flag(retval, "CVodeSetLinearSolver")) { return 1; }

      // Set the user-supplied Jacobian function
      retval = CVodeSetJacFn(cvode_mem, Jac);
      if (check_flag(retval, "CVodeSetJacFn")) { return 1; }
      fac = real_type_2d_view_type("fac", nBatch, number_of_equations);
      udata.fac=fac;
    }
    else
    {
      // Create matrix-free GMRES linear solver
      LS = std::make_unique<sundials::experimental::SUNLinearSolverView>(
        SUNLinSol_SPGMR(y, SUN_PREC_NONE, 0, sunctx));

      // Attach the linear solver to CVODE
      retval = CVodeSetLinearSolver(cvode_mem, LS->Convert(), nullptr);
      if (check_flag(retval, "CVodeSetLinearSolver")) { return 1; }
    }

    // Final time and time between outputs

    const sunrealtype Tf    = SUN_RCONST(tend);
    const sunrealtype dTout = SUN_RCONST(dtmax);

    // Number of output times
    const int Nt = static_cast<int>(ceil(Tf / dTout));

    // Current time and first output time
    sunrealtype t    = T0;
    sunrealtype tout = T0 + dTout;

    // Initial output
    real_type_2d_view_host_type y2d_h((y.HostView()).data(), udata.nbatches, udata.batchSize);
    sundials::kokkos::CopyFromDevice(y);
    Kokkos::fence();

    std::cout << "At t = " << t << std::endl;
    for (int j = 0; j < udata.nbatches; j ++)
    {
      std::cout << "  batch " << j << ": y = " << y2d_h(j, 0) << " "
                << y2d_h(j, 1) << " " << y2d_h(j, 2) << std::endl;
    }


    // Loop over output times
    for (int iout = 0; iout < Nt; iout++)
    {
      // Advance in time
      retval = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
      if (check_flag(retval, "CVode")) { break; }

      // Output solution from some batches
      sundials::kokkos::CopyFromDevice(y);
      Kokkos::fence();
      std::cout << "At t = " << t << std::endl;
      for (int j = 0; j < udata.nbatches; j += 10)
      {
        std::cout << "  batch " << j << ": y = " << y2d_h(j, 0) << " "
                  << y2d_h(j, 1) << " " << y2d_h(j, 2) << std::endl;
      }

      tout += dTout;
      tout = (tout > Tf) ? Tf : tout;
    }

        // Print some final statistics
    long int nst, nfe, nsetups, nje, nni, ncfn, netf;

    retval = CVodeGetNumSteps(cvode_mem, &nst);
    check_flag(retval, "CVodeGetNumSteps");
    retval = CVodeGetNumRhsEvals(cvode_mem, &nfe);
    check_flag(retval, "CVodeGetNumRhsEvals");
    retval = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
    check_flag(retval, "CVodeGetNumLinSolvSetups");
    retval = CVodeGetNumErrTestFails(cvode_mem, &netf);
    check_flag(retval, "CVodeGetNumErrTestFails");
    retval = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
    check_flag(retval, "CVodeGetNumNonlinSolvIters");
    retval = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
    check_flag(retval, "CVodeGetNumNonlinSolvConvFails");
    retval = CVodeGetNumJacEvals(cvode_mem, &nje);
    check_flag(retval, "CVodeGetNumJacEvals");

    std::cout << "\nFinal Statistics:\n"
              << "  Steps            = " << nst << "\n"
              << "  RHS evals        = " << nfe << "\n"
              << "  LS setups        = " << nsetups << "\n"
              << "  Jac evals        = " << nje << "\n"
              << "  NLS iters        = " << nni << "\n"
              << "  NLS fails        = " << ncfn << "\n"
              << "  Error test fails = " << netf << "\n";

    // Free objects
    CVodeFree(&cvode_mem);

#if 0

    if (verbose) {


     if (use_cloned_samples) {
       auto rhs_at_0 = Kokkos::subview(rhs, 0, Kokkos::ALL());
       auto rhs_host_at_0 = Kokkos::create_mirror_view(rhs_at_0);
       Kokkos::deep_copy(rhs_host_at_0, rhs_at_0);
       FILE *fout = fopen(outputFile.c_str(), "w");
       TChem::Test::writeReactionRates(outputFile, number_of_equations, rhs_host_at_0);
       fclose(fout);
     } else {
      auto rhs_host= Kokkos::create_mirror_view(rhs);
      Kokkos::deep_copy(rhs_host, rhs);
      FILE *fout = fopen(outputFile.c_str(), "w");
      TChem::Test::write2DMatrix(outputFile, nBatch, number_of_equations,rhs_host);
      fclose(fout);
     }

    }
#endif

    // fprintf(fout_times, "}\n ");// end index time
    // fclose(fout_times);




  }
  Kokkos::finalize();

  return 0;
}
// Right hand side function dy/dt = f(t,y)
int f(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
{
  auto udata = static_cast<UserData*>(user_data);

  const auto nbatches  = udata->nbatches;


  const auto policy = udata->policy;
  const auto kmcd = udata->kmcd;
  const auto amcd = udata->amcd;
  const auto batchSize = udata->batchSize;
  const auto num_concentration = udata->num_concentration;
  const auto temperature= udata->temperature;
  const auto pressure= udata->pressure;
  const auto const_tracers= udata->const_tracers;

  real_type_2d_view_type y2d(N_VGetDeviceArrayPointer(y), nbatches, batchSize);
  real_type_2d_view_type vals(N_VGetDeviceArrayPointer(y), nbatches, batchSize);
  real_type_2d_view_type rhs(N_VGetDeviceArrayPointer(ydot), nbatches, batchSize);

  const ordinal_type level = 1;
  const ordinal_type number_of_equations = batchSize;
  const ordinal_type per_team_extent = problem_type::getWorkSpaceSize(kmcd,amcd)
    + number_of_equations;
  const std::string profile_name = "TChem::AerosolChemistry::RHS_evaluation";
  Kokkos::Profiling::pushRegion(profile_name);
  Kokkos::parallel_for
  (profile_name,
       policy,
       KOKKOS_LAMBDA(const typename policy_type::member_type& member) {

        const ordinal_type i = member.league_rank();

        const ordinal_type m = problem_type::getNumberOfTimeODEs(kmcd,amcd);
        const real_type_1d_view_type rhs_at_i =
        Kokkos::subview(rhs, i, Kokkos::ALL());
        const real_type_1d_view_type vals_at_i =
        Kokkos::subview(vals, i, Kokkos::ALL());

                const real_type_1d_view_type number_conc_at_i =
        Kokkos::subview(num_concentration, i, Kokkos::ALL());

        TChem::Scratch<real_type_1d_view_type> work(member.team_scratch(level),
                                       per_team_extent);

        auto wptr = work.data();

        const real_type_1d_view_type constYs  = Kokkos::subview(const_tracers, i, Kokkos::ALL());

        const ordinal_type problem_workspace_size = problem_type::getWorkSpaceSize(kmcd,amcd);
        auto pw = real_type_1d_view_type(wptr, problem_workspace_size);
        wptr +=problem_workspace_size;
        problem_type problem;
        problem._kmcd = kmcd;
        problem._amcd = amcd;
        /// initialize problem
        // problem._fac = fac_at_i;
        problem._work = pw;
        problem._temperature= temperature(i);
        problem._pressure =pressure(i);
        problem._const_concentration= constYs;
        problem._number_conc =number_conc_at_i;
        problem.computeFunction(member,vals_at_i,rhs_at_i);
      });


 return 0;
}

// Jacobian of f(t,y)
int Jac(sunrealtype t, N_Vector y, N_Vector fy, SUNMatrix J, void* user_data
, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  auto udata = static_cast<UserData*>(user_data);

  const auto nbatches  = udata->nbatches;
  const auto policy = udata->policy;
  const auto kmcd = udata->kmcd;
  const auto amcd = udata->amcd;
  const auto num_concentration = udata->num_concentration;
  const auto temperature= udata->temperature;
  const auto pressure= udata->pressure;
  const auto const_tracers= udata->const_tracers;
  const auto batchSize = udata->batchSize;
  const auto fac= udata->fac;
  auto J_data = sundials::kokkos::GetDenseMat<MatType>(J)->View();

  real_type_2d_view_type y2d(N_VGetDeviceArrayPointer(y), nbatches, batchSize);
  real_type_2d_view_type vals(N_VGetDeviceArrayPointer(y), nbatches, batchSize);
  // real_type_3d_view_type jacobian(SUNDenseMatrix_Data(J), nbatches, batchSize, batchSize);

  const ordinal_type level = 1;
  const ordinal_type number_of_equations = problem_type::getNumberOfTimeODEs(kmcd, amcd);
  const ordinal_type per_team_extent = problem_type::getWorkSpaceSize(kmcd,amcd)
    + number_of_equations;
  const std::string profile_name = "TChem::AerosolChemistry::Jacobian_evaluation";
  Kokkos::Profiling::pushRegion(profile_name);
  Kokkos::parallel_for
  (profile_name,
       policy,
       KOKKOS_LAMBDA(const typename policy_type::member_type& member) {

        const ordinal_type i = member.league_rank();

        const ordinal_type m = problem_type::getNumberOfTimeODEs(kmcd,amcd);
        const real_type_1d_view_type vals_at_i =
        Kokkos::subview(vals, i, Kokkos::ALL());

        const real_type_1d_view_type fac_at_i =
        Kokkos::subview(fac, i, Kokkos::ALL());

        const real_type_2d_view_type jacobian_at_i =
        Kokkos::subview(J_data, i, Kokkos::ALL(), Kokkos::ALL());

                const real_type_1d_view_type number_conc_at_i =
        Kokkos::subview(num_concentration, i, Kokkos::ALL());

        TChem::Scratch<real_type_1d_view_type> work(member.team_scratch(level),
                                       per_team_extent);

        auto wptr = work.data();

        const real_type_1d_view_type constYs  = Kokkos::subview(const_tracers, i, Kokkos::ALL());

        const ordinal_type problem_workspace_size = problem_type::getWorkSpaceSize(kmcd,amcd);
        auto pw = real_type_1d_view_type(wptr, problem_workspace_size);
        wptr +=problem_workspace_size;
        problem_type problem;
        problem._kmcd = kmcd;
        problem._amcd = amcd;
        /// initialize problem
        problem._fac = fac_at_i;
        problem._work = pw;
        problem._temperature= temperature(i);
        problem._pressure =pressure(i);
        problem._const_concentration= constYs;
        problem._number_conc =number_conc_at_i;
        problem.computeNumericalJacobian(member,vals_at_i,jacobian_at_i);
        // for (int k=0;k<m;++k)
        //   for (int j=0;j<m;++j)
        //    J_data(i, k, j) = jacobian_at_i(k,j);
      });


  return 0;

}
