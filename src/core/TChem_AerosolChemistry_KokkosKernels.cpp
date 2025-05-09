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
#include "TChem_AerosolChemistry_KokkosKernels.hpp"

namespace TChem
{
	  template<typename PolicyType,
           typename ValueType,
           typename DeviceType>
  void
  AerosolChemistry_KokkosKernels_TemplateRun( /// required template arguments
          const std::string& profile_name,
          const ValueType& dummyValueType,
          /// team size setting
          const PolicyType& policy,
          /// input
          const Tines::value_type_2d_view<real_type, DeviceType>& tol,
          const Tines::value_type_2d_view<real_type, DeviceType>& fac,
          const Tines::value_type_1d_view<time_advance_type, DeviceType>& tadv,
          const Tines::value_type_2d_view<real_type, DeviceType>& state,
          const Tines::value_type_2d_view<real_type, DeviceType>& number_conc,
          /// output
          const Tines::value_type_1d_view<real_type, DeviceType>& t_out,
          const Tines::value_type_1d_view<real_type, DeviceType>& dt_out,
          const Tines::value_type_2d_view<real_type, DeviceType>& state_out,
          /// const data from kinetic model
          const KineticModelNCAR_ConstData<DeviceType>& kmcd,
          const AerosolModel_ConstData<DeviceType>& amcd
          ) {
    Kokkos::Profiling::pushRegion(profile_name);
    using policy_type = PolicyType;

    using real_type_0d_view_type = Tines::value_type_0d_view<real_type, DeviceType>;
    using real_type_1d_view_type = Tines::value_type_1d_view<real_type, DeviceType>;
    using real_type_2d_view_type = Tines::value_type_2d_view<real_type, DeviceType>;
    using range_type = Kokkos::pair<ordinal_type, ordinal_type>;
    using problem_type = TChem::Impl::AerosolChemistry_Problem<real_type,DeviceType>;
    const ordinal_type level = 1;
    const ordinal_type per_team_extent = AerosolChemistry_KokkosKernels::getWorkSpaceSize(kmcd,amcd);
    Kokkos::parallel_for
      (profile_name,
       policy,
       KOKKOS_LAMBDA(const typename policy_type::member_type& member) {
  const ordinal_type i = member.league_rank();

  const auto tadv_at_i = tadv(i);
  const real_type t_end = tadv_at_i._tend;
  const real_type_0d_view_type t_out_at_i = Kokkos::subview(t_out, i);
  if (t_out_at_i() < t_end) {
    const real_type_1d_view_type state_at_i =
        Kokkos::subview(state, i, Kokkos::ALL());

      const real_type_1d_view_type number_conc_at_i =
        Kokkos::subview(number_conc, i, Kokkos::ALL());

      const real_type_1d_view_type state_out_at_i =
        Kokkos::subview(state_out, i, Kokkos::ALL());
      const real_type_1d_view_type fac_at_i =
        Kokkos::subview(fac, i, Kokkos::ALL());
      const real_type_0d_view_type dt_out_at_i = Kokkos::subview(dt_out, i);
      Scratch<real_type_1d_view_type> work(member.team_scratch(level),
                                       per_team_extent);
      auto wptr = work.data();
      // const real_type_1d_view_type ww(wptr, work.extent(0));
      const ordinal_type total_n_species = kmcd.nSpec + amcd.nParticles*amcd.nSpec;
      Impl::StateVector<real_type_1d_view_type> sv_at_i(total_n_species, state_at_i);
      Impl::StateVector<real_type_1d_view_type> sv_out_at_i(total_n_species,
                                                        state_out_at_i);
      TCHEM_CHECK_ERROR(!sv_at_i.isValid(),
                        "Error: input state vector is not valid");
      TCHEM_CHECK_ERROR(!sv_out_at_i.isValid(),
                        "Error: input state vector is not valid");
    {
      const ordinal_type max_num_time_iterations = tadv_at_i._num_time_iterations_per_interval;
      const real_type
        dt_in = tadv_at_i._dt,
        dt_min = tadv_at_i._dtmin,
        dt_max = tadv_at_i._dtmax;
      const real_type t_beg = tadv_at_i._tbeg;

      const real_type temperature = sv_at_i.Temperature();
      const real_type pressure = sv_at_i.Pressure();
      const real_type density = sv_at_i.Density();
      const real_type_1d_view_type Ys = sv_at_i.MassFractions();
      const ordinal_type n_active_gas_species = kmcd.nSpec - kmcd.nConstSpec;
      const real_type_1d_view_type activeYs = Kokkos::subview(Ys,
          range_type(0, n_active_gas_species));
      const real_type_1d_view_type constYs = Kokkos::subview(Ys,
          range_type(n_active_gas_species, kmcd.nSpec));
      const real_type_1d_view_type partYs = Kokkos::subview(Ys,
          range_type(kmcd.nSpec, total_n_species));

      const real_type_0d_view_type temperature_out(sv_out_at_i.TemperaturePtr());
      const real_type_0d_view_type pressure_out(sv_out_at_i.PressurePtr());
      const real_type_1d_view_type Ys_out = sv_out_at_i.MassFractions();
      const real_type_0d_view_type density_out(sv_out_at_i.DensityPtr());

      const ordinal_type m = problem_type::getNumberOfTimeODEs(kmcd,amcd);     

      /// problem setup
      const ordinal_type problem_workspace_size = problem_type::getWorkSpaceSize(kmcd,amcd);

      problem_type problem;
      problem._kmcd = kmcd;
      problem._amcd = amcd;

      /// problem workspace
      auto wptr = work.data();
      auto pw = real_type_1d_view_type(wptr, problem_workspace_size);
      wptr += problem_workspace_size;

      /// temporal var
      real_type_1d_view_type vals(wptr, m);
      wptr += m;

      real_type_1d_view_type rhs(wptr, m);
      wptr += m;

      real_type_1d_view_type update(wptr, m);
      wptr += m;

      ordinal_type subTemp_dims[4];
      AerosolChemistry_KokkosKernels::get_subTemp_dims(m, subTemp_dims);
      real_type_2d_view_type subTemp(wptr, subTemp_dims[0], subTemp_dims[1]);
      wptr += subTemp_dims[0]*subTemp_dims[1];
      //
      real_type_2d_view_type subTemp2(wptr,subTemp_dims[2], subTemp_dims[3]);
      wptr += subTemp_dims[2]*subTemp_dims[3];

      /// error check
      const ordinal_type workspace_used(wptr - work.data()), workspace_extent(work.extent(0));
      if (workspace_used > workspace_extent) {
        Kokkos::abort("Error Ignition ZeroD Sacado : workspace used is larger than it is provided\n");
      }

      /// initialize problem
      problem._fac = fac_at_i;
      problem._work = pw;
      problem._temperature = temperature;
      problem._pressure = pressure;
      problem._const_concentration = constYs;
      problem._number_conc = number_conc_at_i;
      // active gas species
      for (ordinal_type i=0;i<n_active_gas_species;++i){
        vals(i) = activeYs(i);
      }

      for (ordinal_type i=n_active_gas_species;i<total_n_species- kmcd.nConstSpec;++i)
      {
        vals(i) = partYs(i-n_active_gas_species);
      }

      real_type t = t_out_at_i();
      Impl::KokkosKernelsODE my_ode(m, problem, member);

      auto t_start = 0.0; //tadv_at_i._tbeg;
      auto dt = tadv_at_i._dt; 
      auto t_end = dt;//adv_at_i._tend;


      //NOTE: hard-coded value.
      // I got this from the example code. 
      real_type max_step = (t_end - t_start) / 10;
     
      KokkosODE::Experimental::BDFSolve(my_ode, t_start, t_end, dt, max_step,
                                      vals, vals, subTemp, subTemp2, rhs, update);

      t_out_at_i() += dt;
      dt_out_at_i() = dt;

      // active gas species
      for (ordinal_type i=0;i<n_active_gas_species;++i)
        Ys_out(i) = vals(i);
      // particle species
      // Ys_out also contains invariant species
      // total_n_species includes invariant species
      for (ordinal_type i=n_active_gas_species;i<total_n_species-kmcd.nConstSpec;++i)
        Ys_out(i+kmcd.nConstSpec)=vals(i);

      temperature_out() =temperature;

    }
  }
      });
    Kokkos::Profiling::popRegion();
  }



#define TCHEM_RUN_AEROSOL_CHEMISTRY()     \
  AerosolChemistry_KokkosKernels_TemplateRun(          \
          profile_name,       \
          value_type(),       \
          policy,       \
          tol,          \
          fac,          \
          tadv,         \
          state,        \
          number_conc,  \
          t_out,        \
          dt_out,       \
          state_out,    \
          kmcd,        \
          amcd )

  void
  AerosolChemistry_KokkosKernels::runHostBatch( /// thread block size
           typename UseThisTeamPolicy<host_exec_space>::type& policy,
           /// input
           const real_type_2d_view_host& tol,
           const real_type_2d_view_host& fac,
           const time_advance_type_1d_view_host& tadv,
           const real_type_2d_view_host& state,
           const real_type_2d_view_host& number_conc,
           /// output
           const real_type_1d_view_host& t_out,
           const real_type_1d_view_host& dt_out,
           const real_type_2d_view_host& state_out,
           /// const data from kinetic model
           const KineticModelNCAR_ConstData<interf_host_device_type>& kmcd,
           const AerosolModel_ConstData<interf_host_device_type>& amcd
           ) {
    const std::string profile_name = "TChem::AerosolChemistry::runHostBatch::kmcd array";
    // Note: we do not support SACADO and Kokkos-kernels. Thus, BDF::KokkosKernel solver uses a numerical Jacobian.
    using value_type = real_type;
    TCHEM_RUN_AEROSOL_CHEMISTRY();
  }// namespace TChem

    void
  AerosolChemistry_KokkosKernels::runDeviceBatch( /// thread block size
           typename UseThisTeamPolicy<exec_space>::type& policy,
           /// input
           const real_type_2d_view& tol,
           const real_type_2d_view& fac,
           const time_advance_type_1d_view& tadv,
           const real_type_2d_view& state,
           const real_type_2d_view& number_conc,
           /// output
           const real_type_1d_view& t_out,
           const real_type_1d_view& dt_out,
           const real_type_2d_view& state_out,
           /// const data from kinetic model
           const KineticModelNCAR_ConstData<device_type>& kmcd,
           const AerosolModel_ConstData<device_type>& amcd
           ) {
    const std::string profile_name = "TChem::AerosolChemistry::runDeviceBatch::kmcd array";
    using value_type = real_type;
    TCHEM_RUN_AEROSOL_CHEMISTRY();
  }// namespace TChem

}
