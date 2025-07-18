/* =====================================================================================
TChem_atm version 2.0
Copyright (2024) NTESS
https://github.com/sandialabs/TChem-atm

Copyright 2024 National Technology & Engineering Solutions of Sandia, LLC
(NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
Government retains certain rights in this software.

This file is part of TChem-atm_atm. TChem_atm is open source software: you can
redistribute it and/or modify it under the terms of BSD 2-Clause License
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

using ordinal_type = TChem::ordinal_type;
using real_type = TChem::real_type;
using time_advance_type = TChem::time_advance_type;

using real_type_0d_view = TChem::real_type_0d_view;
using real_type_1d_view = TChem::real_type_1d_view;
using real_type_2d_view = TChem::real_type_2d_view;

using time_advance_type_0d_view = TChem::time_advance_type_0d_view;
using time_advance_type_1d_view = TChem::time_advance_type_1d_view;

using real_type_0d_view_host = TChem::real_type_0d_view_host;
using real_type_1d_view_host = TChem::real_type_1d_view_host;
using real_type_2d_view_host = TChem::real_type_2d_view_host;

using time_advance_type_0d_view_host = TChem::time_advance_type_0d_view_host;
using time_advance_type_1d_view_host = TChem::time_advance_type_1d_view_host;

int main(int argc, char *argv[]) {

  /// default inputs
  std::string prefixPath("");

  const real_type zero(0);
  real_type tbeg(0), tend(1);
  real_type dtmin(1e-8), dtmax(1e-1);
  real_type rtol_time(1e-4), atol_newton(1e-10), rtol_newton(1e-6);
  int num_time_iterations_per_interval(1e1), max_num_time_iterations(1e3),
      max_num_newton_iterations(100), jacobian_interval(1);
  real_type atol_time(1e-12);
  int nBatch(1), team_size(-1), vector_size(-1);
  std::string outputFileTimes("wall_times.json");
  bool verbose(true);
  std::string chemFile("chem.yaml");
  std::string aeroFile("aero.yaml");
  std::string outputFile("AerosolChemistry.dat");
  std::string inputFile("None");
  std::string input_file_particles("None");
  bool use_cloned_samples(false);
  int number_of_particles(-1);

  // For scaling studies, we must execute this example many times.
  // Thus, we do not want to write the solution to a file.
  // In those cases, we pass write_time_profiles false.
  bool write_time_profiles(true);
  /// parse command line arguments
  TChem::CommandLineParser opts(
      "This example computes the solution of gas chemistry problem");
  opts.set_option<std::string>(
      "inputsPath", "path to input files e.g., data/inputs", &prefixPath);
  opts.set_option<std::string>("chemfile", "Chem file name e.g., chem.inp",
                               &chemFile);
  opts.set_option<std::string>("aerofile", "Aerosol Chem file name e.g., aero.yaml",
                               &aeroFile);
  opts.set_option<std::string>("inputfile", "Chem file name e.g., chem.inp",
                               &inputFile);
  opts.set_option<std::string>("inputfile_particles", "input file name e.g., chem.inp",
                               &input_file_particles);
  opts.set_option<real_type>("tbeg", "Time begin", &tbeg);
  opts.set_option<real_type>("tend", "Time end", &tend);
  opts.set_option<real_type>("dtmin", "Minimum time step size", &dtmin);
  opts.set_option<real_type>("dtmax", "Maximum time step size", &dtmax);
  opts.set_option<real_type>(
      "atol-newton", "Absolute tolerance used in newton solver", &atol_newton);
  opts.set_option<real_type>(
      "rtol-newton", "Relative tolerance used in newton solver", &rtol_newton);
  opts.set_option<std::string>(
      "outputfile", "Output file name e.g., IgnSolution.dat", &outputFile);
  opts.set_option<real_type>(
      "tol-time", "Tolerence used for adaptive time stepping", &rtol_time);
  opts.set_option<int>("time-iterations-per-interval",
                       "Number of time iterations per interval to store qoi",
                       &num_time_iterations_per_interval);
  opts.set_option<std::string>("outputfile_times",
                               "Wal times file name e.g., times.json",
                               &outputFileTimes);
  opts.set_option<int>("max-time-iterations",
                       "Maximum number of time iterations",
                       &max_num_time_iterations);
  opts.set_option<int>("jacobian-interval",
                       "Jacobians are evaluated once in this interval",
                       &jacobian_interval);
  opts.set_option<int>("max-newton-iterations",
                       "Maximum number of newton iterations",
                       &max_num_newton_iterations);
  opts.set_option<bool>(
      "verbose", "If true, printout the first Jacobian values", &verbose);
  opts.set_option<bool>(
      "write-time-profiles", "If true, this example will write the time profile output to a file.", &write_time_profiles);
  opts.set_option<int>("team_thread_size", "time thread size ", &team_size);
  opts.set_option<int>("batch_size", " number of batches or samples, e.g. 10  ", &nBatch);
  opts.set_option<int>("vector_thread_size", "vector thread size ",
                       &vector_size);
  opts.set_option<int>("number_of_particles", "Set the number of particles; this will overwrite the value from the input file. ", &number_of_particles);
  opts.set_option<real_type>("atol-time", "Absolute tolerance used for adaptive time stepping", &atol_time);
  opts.set_option<bool>(
      "use_cloned_samples", "If true, one state vector will be cloned.", &use_cloned_samples);


  const bool r_parse = opts.parse(argc, argv);
  if (r_parse)
    return 0; // print help return

  Kokkos::initialize(argc, argv);
  {
    const bool detail = false;

    // by default use input condition from chemFile
    if (inputFile == "None") {
       inputFile=chemFile;
    }
    // check first in chem file
    if (input_file_particles == "None") {
       input_file_particles=chemFile;
    }

#if defined(TCHEM_EXAMPLE_KOKKOS_KERNELS_ODE)
    printf("   TChem is running ATM model box with KokkosKernels\n");
#else
    printf("   TChem is running ATM model box with Tines\n");
#endif
    TChem::exec_space().print_configuration(std::cout, detail);
    TChem::host_exec_space().print_configuration(std::cout, detail);

    /// construct kmd and use the view for testing
    printf("kmd parsing ...\n");
    TChem::KineticModelData kmd(chemFile);
    const auto kmcd = TChem::createNCAR_KineticModelConstData<TChem::interf_device_type>(kmd);

    printf("amd parsing ...\n");
    TChem::AerosolModelData amd(aeroFile, kmd);
    if(number_of_particles > 0) {
      amd.setNumberofParticles(number_of_particles);
    }
    const auto amcd = TChem::create_AerosolModelConstData<TChem::interf_device_type>(amd);

    const ordinal_type total_n_species =kmcd.nSpec + amcd.nSpec * amcd.nParticles;
    const ordinal_type stateVecDim =
        TChem::Impl::getStateVectorSize(total_n_species);

    printf("Number of Gas Species %d \n", kmcd.nSpec);
    printf("Number of Gas Reactions %d \n", kmcd.nReac);

    printf("Number of Aerosol Species %d \n", amcd.nSpec);
    printf("Number of Aerosol Particles %d \n", amcd.nParticles);

    printf("stateVecDim %d \n", stateVecDim);

    const auto speciesNamesHost = Kokkos::create_mirror_view(kmcd.speciesNames);
    Kokkos::deep_copy(speciesNamesHost, kmcd.speciesNames);

    FILE *fout;
    if (write_time_profiles) {
     fout = fopen(outputFile.c_str(), "w");
    }
    // read scenario condition from yaml file
   // read scenario condition from yaml file
    real_type_2d_view_host state_scenario_host;
    ordinal_type nbatch_files=1;
    printf("conditions parsing ...\n");
    TChem::AtmChemistry::setScenarioConditions(inputFile,
     speciesNamesHost, kmcd.nSpec, stateVecDim, state_scenario_host, nbatch_files);
    printf("Number of const species %d \n", kmcd.nConstSpec);

    real_type_2d_view_host num_concentration_host;
    real_type_2d_view_host state_host;


    // scenario particles
    amd.scenarioConditionParticles(input_file_particles, nbatch_files, num_concentration_host, state_scenario_host);

    real_type_2d_view num_concentration;
    real_type_2d_view state;

    if (nbatch_files == 1 && use_cloned_samples && nBatch > 1) {
      // only clone samples if nbatch_files is 1
      printf("-------------------------------------------------------\n");
      printf("--------------------Warning----------------------------\n");
      printf("Using cloned samples ... only for numerical experiments\n");
      state = real_type_2d_view("StateVector Devices", nBatch, stateVecDim);
      state_host = Kokkos::create_mirror_view(state);
      auto state_scenario_host_at_0 = Kokkos::subview(state_scenario_host, 0, Kokkos::ALL);
      auto state_at_0 = Kokkos::subview(state, 0, Kokkos::ALL);
      Kokkos::deep_copy(state_at_0, state_scenario_host_at_0);
      TChem::Test::cloneView(state);
      //make sure to copy state to host view.
      Kokkos::deep_copy(state_host, state);
      //scenario particles
      auto num_concentration_host_at_0 = Kokkos::subview(num_concentration_host, 0, Kokkos::ALL);
      num_concentration = real_type_2d_view("num_concentration", nBatch, amd.nParticles_);
      auto num_concentration_at_0 = Kokkos::subview(num_concentration, 0, Kokkos::ALL);
      Kokkos::deep_copy(num_concentration_at_0, num_concentration_host_at_0);
      TChem::Test::cloneView(num_concentration);

    } else {
      nBatch = nbatch_files;
      state = real_type_2d_view("StateVector Devices", nBatch, stateVecDim);
      Kokkos::deep_copy(state, state_scenario_host);
      state_host = state_scenario_host;

      // scenario particles
      num_concentration = real_type_2d_view("num_concentration", nBatch, amd.nParticles_);
      Kokkos::deep_copy(num_concentration, num_concentration_host);
    }
    printf("Number of nbatch %d \n",nBatch);
    auto writeState = [](const ordinal_type iter,
                         const real_type_1d_view_host _t,
                         const real_type_1d_view_host _dt,
                         const real_type_2d_view_host _state_at_i,
                         FILE *fout) { // sample, time, density, pressure,
                                       // temperature, mass fraction
      for (size_t sp = 0; sp < _state_at_i.extent(0); sp++) {
        fprintf(fout, "%d \t %15.10e \t  %15.10e \t ", iter, _t(sp), _dt(sp));
        //
        for (ordinal_type k = 0, kend = _state_at_i.extent(1); k < kend; ++k)
          fprintf(fout, "%15.10e \t", _state_at_i(sp, k));

        fprintf(fout, "\n");
      }

    };

    FILE *fout_times = fopen(outputFileTimes.c_str(), "w");
    fprintf(fout_times, "{\n");
    fprintf(fout_times, " \"Aerosol Chemistry\": \n {\n");


    TChem::TeamConfOutput team_conf_output;
    const auto exec_space_instance = TChem::exec_space();

    using policy_type =
          typename TChem::UseThisTeamPolicy<TChem::exec_space>::type;

#if defined(TCHEM_EXAMPLE_KOKKOS_KERNELS_ODE)

      /// team policy
      // The Kokkos-kernels BDF solver is designed for range policy. Therefore, I will use a team size of 1.
      policy_type policy(exec_space_instance, nBatch, 1);
#else
    policy_type policy(exec_space_instance, nBatch, Kokkos::AUTO());

    if (team_size > 0 && vector_size > 0) {
      // if team_size and vector size are positives
      // here we modified it.
      policy = policy_type(exec_space_instance,  nBatch, team_size, vector_size);
    } else if (team_size > 0 && vector_size < 0) {
      // here we let kokkos pick vector size
      // only set team size
      policy = policy_type(exec_space_instance, nBatch,  team_size);
    }
#endif
{
        ordinal_type number_of_equations(0);

        using problem_type =
            TChem::Impl::AerosolChemistry_Problem<real_type,
                                                          TChem::interf_device_type>;
        number_of_equations = problem_type::getNumberOfTimeODEs(kmcd, amcd);


      const ordinal_type level = 1;
      ordinal_type per_team_extent(0);

#if defined(TCHEM_EXAMPLE_KOKKOS_KERNELS_ODE)
      per_team_extent = TChem::AerosolChemistry_KokkosKernels::getWorkSpaceSize(kmcd, amcd);
#else
      per_team_extent = TChem::AerosolChemistry::getWorkSpaceSize(kmcd, amcd);
#endif
      const ordinal_type per_team_scratch =
          TChem::Scratch<real_type_1d_view>::shmem_size(per_team_extent);
      policy.set_scratch_size(level, Kokkos::PerTeam(per_team_scratch));

      { /// time integration
        real_type_1d_view t("time", nBatch);
        Kokkos::deep_copy(t, tbeg);
        real_type_1d_view dt("delta time", nBatch);
        Kokkos::deep_copy(dt, dtmax);

        real_type_2d_view tol_time("tol time", number_of_equations, 2);
        real_type_1d_view tol_newton("tol newton", 2);

        real_type_2d_view fac("fac", nBatch, number_of_equations);

        auto tol_time_host = Kokkos::create_mirror_view(tol_time);
        auto tol_newton_host = Kokkos::create_mirror_view(tol_newton);

        /// tune tolerence
        {
          for (ordinal_type i = 0, iend = tol_time.extent(0); i < iend; ++i) {
            tol_time_host(i, 0) = atol_time;
            tol_time_host(i, 1) = rtol_time;
          }
          tol_newton_host(0) = atol_newton;
          tol_newton_host(1) = rtol_newton;
        }

        Kokkos::deep_copy(tol_time, tol_time_host);
        Kokkos::deep_copy(tol_newton, tol_newton_host);

        time_advance_type tadv_default;
        tadv_default._tbeg = tbeg;
        tadv_default._tend = tend;
        tadv_default._dt = dtmin;
        tadv_default._dtmin = dtmin;
        tadv_default._dtmax = dtmax;
        tadv_default._max_num_newton_iterations = max_num_newton_iterations;
        tadv_default._num_time_iterations_per_interval =
            num_time_iterations_per_interval;
        tadv_default._jacobian_interval = jacobian_interval;

        time_advance_type_1d_view tadv("tadv", nBatch);
        Kokkos::deep_copy(tadv, tadv_default);

        real_type_1d_view_host t_host;
        real_type_1d_view_host dt_host;
        t_host = real_type_1d_view_host("time host", nBatch);
        dt_host = real_type_1d_view_host("dt host", nBatch);

        ordinal_type iter = 0;
        if (write_time_profiles) {
          fprintf(fout, "%s \t %s \t %s \t ", "iter", "t", "dt");
          fprintf(fout, "%s \t %s \t %s \t", "Density[kg/m3]", "Pressure[Pascal]",
                "Temperature[K]");

          for (ordinal_type k = 0; k < kmcd.nSpec; k++)
            fprintf(fout, "%s \t", &speciesNamesHost(k, 0));
          //given aero idx return species name
          std::map<int, std::string> aero_idx_sp_name;
          for (std::map<std::string, int>::iterator
             i = amd.aerosol_sp_name_idx_.begin();
            i != amd.aerosol_sp_name_idx_.end(); ++i)
            aero_idx_sp_name[i->second] = i->first;

          for (ordinal_type ipart = 0; ipart < amd.nParticles_; ipart++)
          {
            for (ordinal_type isp = 0; isp < amd.nSpec_; isp++)
            {
              // std::cout << "species Name : "<< aero_idx_sp_name[i] << "\n";
              auto aero_sp_name = aero_idx_sp_name[isp]+"_p"+std::to_string(ipart);
              fprintf(fout, "%s \t", aero_sp_name.c_str());
            }// isp
          }// ipar

          fprintf(fout, "\n");
          writeState(-1, t_host, dt_host, state_host, fout);
        }
        real_type tsum(0);
        Kokkos::Timer timer;


        for (; iter < max_num_time_iterations && tsum <= tend * 0.9999;
             ++iter) {

          timer.reset();

#if defined(TCHEM_EXAMPLE_KOKKOS_KERNELS_ODE)
          TChem::AerosolChemistry_KokkosKernels::runDeviceBatch(
              policy, tol_time, fac, tadv, state, num_concentration, t, dt, state,
              kmcd, amcd);
#else
          TChem::AerosolChemistry::runDeviceBatch(
              policy, tol_newton, tol_time, fac, tadv, state, num_concentration, t, dt, state,
              team_conf_output,
              kmcd, amcd);
#endif

          exec_space_instance.fence();

          const real_type t_device_batch = timer.seconds();
          fprintf(fout_times, "\"%s%d\": %20.14e, \n", "wall_time_iter_", iter,
                  t_device_batch);

          if (write_time_profiles) {
            Kokkos::deep_copy(dt_host, dt);
            Kokkos::deep_copy(t_host, t);
            Kokkos::deep_copy(state_host, state);
            writeState(iter, t_host, dt_host, state_host, fout);
          }
          /// carry over time and dt computed in this step
          tsum = zero;
          Kokkos::parallel_reduce(
              Kokkos::RangePolicy<TChem::exec_space>(0, nBatch),
              KOKKOS_LAMBDA(const ordinal_type &i, real_type &update) {
                tadv(i)._tbeg = t(i);
                tadv(i)._dt = dt(i);
                // printf("t %e, dt %e\n", t(i), dt(i));
                // printf("tadv t %e, tadv dt %e\n", tadv(i)._tbeg, tadv(i)._dt
                // );
                update += t(i);
              },
              tsum);
          Kokkos::fence();
          tsum /= nBatch;
        }// end iter
        fprintf(fout_times, "%s: %d, \n", "\"number_of_time_iters\"", iter);

      }
    }

    fprintf(fout_times, "%s: %d \n", "\"number_of_samples\"", nBatch);
    fprintf(fout_times, "} \n "); // reaction rates

    fprintf(fout_times, ", \n ");
    fprintf(fout_times, " \"Kokkos Settings\": \n  {\n");
    fprintf(fout_times, "%s: %d, \n", "\"team_size\"", team_size);
    fprintf(fout_times, "%s: %d, \n", "\"vector_size\"", vector_size);
    fprintf(fout_times, "%s: %d, \n", "\"vector_length_max\"", policy.vector_length_max());
    fprintf(fout_times, "%s: %d, \n", "\"impl_vector_length\"", policy.impl_vector_length());
    fprintf(fout_times, "%s: %d, \n", "\"team_size_recommended\"", team_conf_output.team_size_recommended);
    fprintf(fout_times, "%s: %d \n", "\"team_size_max\"", team_conf_output.team_size_max);
    fprintf(fout_times, "} \n "); // reaction rates


    if (write_time_profiles) {
      fclose(fout);
    }
    fprintf(fout_times, "}\n "); // end index time
    fclose(fout_times);
  }
  Kokkos::finalize();

  return 0;
}
