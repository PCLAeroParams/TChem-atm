/* =====================================================================================
TChem_atm version 2.0
Copyright (2024) NTESS
https://github.com/sandialabs/TChem-atm

Copyright 2024 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains
certain rights in this software.

This file is part of TChem-atm_atm. TChem_atm is open source software: you can redistribute it
and/or modify it under the terms of BSD 2-Clause License
(https://opensource.org/licenses/BSD-2-Clause). A copy of the licese is also
provided under the main directory

Questions? Contact Oscar Diaz-Ibarra at <odiazib@sandia.gov>, or
           Mike Schmidt at <mjschm@sandia.gov>, or
           Cosmin Safta at <csafta@sandia.gov>

Sandia National Laboratories, New Mexico/Livermore, NM/CA, USA
===================================================================================== */
#include "TChem_CommandLineParser.hpp" // done
#include "TChem_KineticModelData.hpp" // done
#include "TChem_Util.hpp"

#include "TChem_AtmosphericChemistry.hpp" // done

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
// FIXME: remove TPL YAML
 // FIXME: use consistent name format
 // FIXME: I will keep tchem name space for now; but I want to change tchem to tchem_atm  (or other name that we decided to keep)
#if defined(TCHEM_ATM_ENABLE_TPL_YAML_CPP)

  /// default inputs
  std::string prefixPath("");

  bool write_time_profiles(true);
  const real_type zero(0);
  real_type tbeg(0), tend(1);
  real_type dtmin(1e-8), dtmax(1e-1);
  real_type rtol_time(1e-4), atol_newton(1e-10), rtol_newton(1e-6);
  int num_time_iterations_per_interval(1e1), max_num_time_iterations(1e3), max_num_newton_iterations(100),
      jacobian_interval(1);

  int nBatch(1), team_size(-1), vector_size(-1);
  bool verbose(false);
  std::string chemFile("chem.yaml");
  std::string outputFile("AtmosphericChemistry.dat");
  std::string outputFileTimes("wall_times.json");
  bool use_cloned_samples(false);

  /// parse command line arguments
  TChem::CommandLineParser opts("This example computes the solution of an ignition problem");
  opts.set_option<std::string>("inputsPath", "path to input files e.g., data/inputs", &prefixPath);
  opts.set_option<std::string>("chemfile", "Chem file name e.g., chem.inp", &chemFile);
  opts.set_option<real_type>("tbeg", "Time begin", &tbeg);
  opts.set_option<real_type>("tend", "Time end", &tend);
  opts.set_option<real_type>("dtmin", "Minimum time step size", &dtmin);
  opts.set_option<real_type>("dtmax", "Maximum time step size", &dtmax);
  opts.set_option<real_type>("atol-newton", "Absolute tolerance used in newton solver", &atol_newton);
  opts.set_option<real_type>("rtol-newton", "Relative tolerance used in newton solver", &rtol_newton);
  opts.set_option<std::string>("outputfile", "Output file name e.g., IgnSolution.dat", &outputFile);
  opts.set_option<real_type>("tol-time", "Tolerence used for adaptive time stepping", &rtol_time);
  opts.set_option<int>("time-iterations-per-interval", "Number of time iterations per interval to store qoi",
                       &num_time_iterations_per_interval);
  opts.set_option<int>("max-time-iterations", "Maximum number of time iterations", &max_num_time_iterations);
  opts.set_option<int>("jacobian-interval", "Jacobians are evaluated once in this interval", &jacobian_interval);
  opts.set_option<int>("max-newton-iterations", "Maximum number of newton iterations", &max_num_newton_iterations);
  opts.set_option<bool>("verbose", "If true, printout the first Jacobian values", &verbose);
  opts.set_option<int>("team_thread_size", "User defined team size", &team_size);
  opts.set_option<int>("vector_thread_size", "User defined vector size", &vector_size);
  opts.set_option<int>("batch_size", " number of batches or samples, e.g. 10  ", &nBatch);
  opts.set_option<std::string>("outputfile_times",
                               "Wal times file name e.g., times.json",
                               &outputFileTimes);
  opts.set_option<bool>(
      "write-time-profiles", "If true, this example will write the time profile output to a file.", &write_time_profiles);
  opts.set_option<bool>(
      "use_cloned_samples", "If true, one state vector will be cloned.", &use_cloned_samples);
  const bool r_parse = opts.parse(argc, argv);
  if (r_parse)
    return 0; // print help return

  Kokkos::initialize(argc, argv);
  {
    const bool detail = false;

    TChem::exec_space().print_configuration(std::cout, detail);
    TChem::host_exec_space().print_configuration(std::cout, detail);

    using device_type = typename Tines::UseThisDevice<TChem::exec_space>::type;

    /// construct kmd and use the view for testing
    TChem::KineticModelData kmd(chemFile);
    const auto kmcd = TChem::createNCAR_KineticModelConstData<device_type>(kmd);

    const ordinal_type stateVecDim = TChem::Impl::getStateVectorSize(kmcd.nSpec);

    printf("Number of Species %d \n", kmcd.nSpec);
    printf("Number of Reactions %d \n", kmcd.nReac);
    const auto speciesNamesHost = Kokkos::create_mirror_view(kmcd.speciesNames);
    Kokkos::deep_copy(speciesNamesHost, kmcd.speciesNames);

    FILE *fout;
    if (write_time_profiles) {
     fout = fopen(outputFile.c_str(), "w");
    }

    FILE *fout_times = fopen(outputFileTimes.c_str(), "w");
    fprintf(fout_times, "{\n");
    fprintf(fout_times, " \"Aerosol Chemistry\": \n {\n");
    // read scenario condition from yaml file
    real_type_2d_view_host state_scenario_host;
    ordinal_type nbatch_files=1;
    TChem::AtmChemistry ::setScenarioConditions(chemFile, speciesNamesHost, kmcd.nSpec, stateVecDim, state_scenario_host, nbatch_files);

    real_type_2d_view state;
    real_type_2d_view_host state_host;
    if (nbatch_files == 1 && use_cloned_samples && nBatch > 1) {
      printf("-------------------------------------------------------\n");
      printf("--------------------Warning----------------------------\n");
      printf("Using cloned samples ... only for numerical experiments\n");
      state = real_type_2d_view("StateVector Devices", nBatch, stateVecDim);
      state_host = Kokkos::create_mirror_view(state);
      auto state_scenario_host_at_0 = Kokkos::subview(state_scenario_host, 0, Kokkos::ALL);
      auto state_at_0 = Kokkos::subview(state, 0, Kokkos::ALL);
      Kokkos::deep_copy(state_at_0, state_scenario_host_at_0);
      TChem::Test::cloneView(state);
    } else {
      nBatch = nbatch_files;
      state = real_type_2d_view("StateVector Devices", nBatch, stateVecDim);
      Kokkos::deep_copy(state, state_scenario_host);
      state_host = state_scenario_host;

    }
    if (verbose)
    {
      for (ordinal_type i = 0; i < nBatch; i++) {
        printf("Host::Initial condition sample No %d\n", i);
        const auto state_at_i_host = Kokkos::subview(state_host, i, Kokkos::ALL());
        for (ordinal_type k = 0, kend = state_at_i_host.extent(0); k < kend; ++k)
          printf(" %e", state_at_i_host(k));
        printf("\n");
      }
    }
    if (verbose)
    {
    Kokkos::parallel_for(
        Kokkos::RangePolicy<TChem::exec_space>(0, nBatch), KOKKOS_LAMBDA(const ordinal_type &i) {
          printf("Devices::Initial condition sample No %d\n", i);
          auto state_at_i = Kokkos::subview(state, i, Kokkos::ALL());
          for (ordinal_type k = 0, kend = state_at_i.extent(0); k < kend; ++k)
            printf(" %e", state_at_i(k));
          printf("\n");
        });
    }
    //

    auto writeState = [](const ordinal_type iter, const real_type_1d_view_host _t, const real_type_1d_view_host _dt,
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

    auto printState = [](const time_advance_type _tadv, const real_type _t, const real_type_1d_view_host _state_at_i)
    {
      printf("%e %e %e %e %e", _t, _t - _tadv._tbeg, _state_at_i(0), _state_at_i(1), _state_at_i(2));
      for (ordinal_type k = 3, kend = _state_at_i.extent(0); k < kend; ++k)
        printf(" %e", _state_at_i(k));
      printf("\n");
    /// iter, t, dt, rho, pres, temp, Ys ...
    };
    {
      const auto exec_space_instance = TChem::exec_space();

      using policy_type = typename TChem::UseThisTeamPolicy<TChem::exec_space>::type;

      /// team policy
      policy_type policy(exec_space_instance, nBatch, Kokkos::AUTO());

    if (team_size > 0 && vector_size > 0) {
        policy = policy_type(exec_space_instance,  nBatch, team_size, vector_size);
    } else if (team_size > 0 && vector_size < 0) {
      // only set team size
       policy = policy_type(exec_space_instance, nBatch,  team_size);
    }

      const ordinal_type level = 1;
      ordinal_type per_team_extent(0);

      per_team_extent = TChem::AtmosphericChemistry::getWorkSpaceSize(kmcd);

      const ordinal_type per_team_scratch = TChem::Scratch<real_type_1d_view>::shmem_size(per_team_extent);
      policy.set_scratch_size(level, Kokkos::PerTeam(per_team_scratch));

      { /// time integration
        real_type_1d_view t("time", nBatch);
        Kokkos::deep_copy(t, tbeg);
        real_type_1d_view dt("delta time", nBatch);
        Kokkos::deep_copy(dt, dtmax);

        real_type_1d_view_host t_host;
        real_type_1d_view_host dt_host;

        t_host = real_type_1d_view_host("time host", nBatch);
        dt_host = real_type_1d_view_host("dt host", nBatch);

        ordinal_type number_of_equations(0);

        using problem_type = TChem::Impl::AtmosphericChemistry_Problem<real_type, device_type>;
        number_of_equations = problem_type::getNumberOfTimeODEs(kmcd);

        real_type_2d_view tol_time("tol time", number_of_equations, 2);
        real_type_1d_view tol_newton("tol newton", 2);

        real_type_2d_view fac("fac", nBatch, number_of_equations);

        /// tune tolerence
        {
          auto tol_time_host = Kokkos::create_mirror_view(tol_time);
          auto tol_newton_host = Kokkos::create_mirror_view(tol_newton);

          const real_type atol_time = 1e-12;
          for (ordinal_type i = 0, iend = tol_time.extent(0); i < iend; ++i) {
            tol_time_host(i, 0) = atol_time;
            tol_time_host(i, 1) = rtol_time;
          }
          tol_newton_host(0) = atol_newton;
          tol_newton_host(1) = rtol_newton;

          Kokkos::deep_copy(tol_time, tol_time_host);
          Kokkos::deep_copy(tol_newton, tol_newton_host);
        }

        time_advance_type tadv_default;
        tadv_default._tbeg = tbeg;
        tadv_default._tend = tend;
        tadv_default._dt = dtmin;
        tadv_default._dtmin = dtmin;
        tadv_default._dtmax = dtmax;
        tadv_default._max_num_newton_iterations = max_num_newton_iterations;
        tadv_default._num_time_iterations_per_interval = num_time_iterations_per_interval;
        tadv_default._jacobian_interval = jacobian_interval;

        time_advance_type_1d_view tadv("tadv", nBatch);
        Kokkos::deep_copy(tadv, tadv_default);

        /// host views to print QOI
        const auto tadv_at_i = Kokkos::subview(tadv, 0);
        const auto t_at_i = Kokkos::subview(t, 0);
        const auto state_at_i = Kokkos::subview(state, 0, Kokkos::ALL());

        auto tadv_at_i_host = Kokkos::create_mirror_view(tadv_at_i);
        auto t_at_i_host = Kokkos::create_mirror_view(t_at_i);
        auto state_at_i_host = Kokkos::create_mirror_view(state_at_i);

        ordinal_type iter = 0;
        /// print of store QOI for the first sample
        if (verbose)
        {
          /// could use cuda streams
          Kokkos::deep_copy(tadv_at_i_host, tadv_at_i);
          Kokkos::deep_copy(t_at_i_host, t_at_i);
          Kokkos::deep_copy(state_at_i_host, state_at_i);
          printState(tadv_at_i_host(), t_at_i_host(), state_at_i_host);
        }


        Kokkos::deep_copy(dt_host, dt);
        Kokkos::deep_copy(t_host, t);

        if (write_time_profiles) {

        fprintf(fout, "%s \t %s \t %s \t ", "iter", "t", "dt");
        fprintf(fout, "%s \t %s \t %s \t", "Density[kg/m3]", "Pressure[Pascal]", "Temperature[K]");

        for (ordinal_type k = 0; k < kmcd.nSpec; k++)
          fprintf(fout, "%s \t", &speciesNamesHost(k, 0));
        fprintf(fout, "\n");
        writeState(-1, t_host, dt_host, state_host, fout);
        }

        Kokkos::Timer timer;
        real_type tsum(0);
        for (; iter < max_num_time_iterations && tsum <= tend * 0.9999; ++iter) {

          timer.reset();
          TChem::AtmosphericChemistry::runDeviceBatch(policy, tol_newton, tol_time, fac, tadv, state, t, dt, state,
                                                      kmcd);
          exec_space_instance.fence();
          const real_type t_device_batch = timer.seconds();
          fprintf(fout_times, "\"%s%d\": %20.14e, \n", "wall_time_iter_", iter,
                  t_device_batch);

          if (verbose) {
            /// could use cuda streams
            Kokkos::deep_copy(tadv_at_i_host, tadv_at_i);
            Kokkos::deep_copy(t_at_i_host, t_at_i);
            Kokkos::deep_copy(state_at_i_host, state_at_i);
            printState(tadv_at_i_host(), t_at_i_host(), state_at_i_host);
          }

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
                // printf("tadv t %e, tadv dt %e\n", tadv(i)._tbeg, tadv(i)._dt );
                update += t(i);
              },
              tsum);
          Kokkos::fence();
          tsum /= nBatch;
        }
        fprintf(fout_times, "%s: %d, \n", "\"number_of_time_iters\"", iter);
      }
    }
    Kokkos::fence(); /// timing purpose
    if (verbose)
    {
    Kokkos::parallel_for(
        Kokkos::RangePolicy<TChem::exec_space>(0, nBatch), KOKKOS_LAMBDA(const ordinal_type &i) {
          printf("Devices:: Solution sample No %d\n", i);
          auto state_at_i = Kokkos::subview(state, i, Kokkos::ALL());
          for (ordinal_type k = 0, kend = state_at_i.extent(0); k < kend; ++k)
            printf(" %e", state_at_i(k));
          printf("\n");
        });
    }
    fprintf(fout_times, "%s: %d \n", "\"number_of_samples\"", nBatch);
    fprintf(fout_times, "} \n "); // reaction rates

    if (write_time_profiles) {
      fclose(fout);
    }
    fprintf(fout_times, "}\n "); // end index time
    fclose(fout_times);

  }
  Kokkos::finalize();

#else
  printf("This example requires Yaml ...\n");
#endif

  return 0;
}
