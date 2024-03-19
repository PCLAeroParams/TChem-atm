#include "TChem_Util.hpp"

#include "TChem_AerosolChemistry_CVODE.hpp"

namespace TChem
{
	  template<typename PolicyType,
           typename ValueType,
           typename DeviceType>
  void
  AerosolChemistry_CVODE_TemplateRun( /// required template arguments
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
          const Tines::value_type_1d_view<KineticModelNCAR_ConstData<DeviceType>,DeviceType>& kmcds,
          const Tines::value_type_1d_view<AerosolModel_ConstData<DeviceType>,DeviceType>& amcds,
          const Tines::value_type_1d_view<Tines::TimeIntegratorCVODE<real_type,DeviceType>,DeviceType>& cvodes) {
    Kokkos::Profiling::pushRegion(profile_name);
#if defined(TINES_ENABLE_TPL_SUNDIALS)
    using policy_type = PolicyType;

    using real_type_0d_view_type = Tines::value_type_0d_view<real_type, DeviceType>;
    using real_type_1d_view_type = Tines::value_type_1d_view<real_type, DeviceType>;

    using problem_type = TChem::Impl::AerosolChemistry_Problem<real_type,DeviceType>;
    const ordinal_type level = 1;
    const ordinal_type per_team_extent = AerosolChemistry_CVODE::getWorkSpaceSize(kmcds(0),amcds(0));

    /// this assumes host space
    Kokkos::parallel_for
      (profile_name,
       policy,
       KOKKOS_LAMBDA(const typename policy_type::member_type& member) {
  const ordinal_type i = member.league_rank();
  const auto kmcd = (kmcds.extent(0) == 1 ? kmcds(0) : kmcds(i));
  const auto amcd = (amcds.extent(0) == 1 ? amcds(0) : amcds(i));
        auto cvode = cvodes(i);

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
      const real_type_1d_view_type ww(wptr, work.extent(0));

      Impl::StateVector<real_type_1d_view_type> sv_at_i(kmcd.nSpec, state_at_i);
      Impl::StateVector<real_type_1d_view_type> sv_out_at_i(kmcd.nSpec,
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
      const auto activeYs = real_type_1d_view_type(Ys.data(),
                              kmcd.nSpec - kmcd.nConstSpec );
      const auto constYs  = real_type_1d_view_type(Ys.data()
                            + kmcd.nSpec - kmcd.nConstSpec,  kmcd.nSpec );

      const real_type_0d_view_type temperature_out(sv_out_at_i.TemperaturePtr());
      const real_type_0d_view_type pressure_out(sv_out_at_i.PressurePtr());
      const real_type_1d_view_type Ys_out = sv_out_at_i.MassFractions();
      const real_type_0d_view_type density_out(sv_out_at_i.DensityPtr());

      const ordinal_type m = problem_type::getNumberOfTimeODEs(kmcd,amcd);
      auto vals = cvode.getStateVector();

      /// problem setup
      const ordinal_type problem_workspace_size = problem_type::getWorkSpaceSize(kmcd,amcd);

      problem_type problem;
      problem._kmcd = kmcd;
      problem._amcd = amcd;

      /// problem workspace
      auto wptr = work.data();
      auto pw = real_type_1d_view_type(wptr, problem_workspace_size);
      wptr += problem_workspace_size;

      /// error check
      const ordinal_type workspace_used(wptr - work.data()), workspace_extent(work.extent(0));
      if (workspace_used > workspace_extent) {
        Kokkos::abort("Error Ignition ZeroD Sacado : workspace used is larger than it is provided\n");
      }

      /// time integrator workspace
      auto tw = real_type_1d_view_type(wptr, workspace_extent - workspace_used);

      /// initialize problem
      problem._fac = fac_at_i;
      problem._work_cvode = tw; // time integrator workspace
      problem._work = pw;
      problem._temperature= temperature;
      problem._pressure =pressure;
      problem._const_concentration= constYs;
      problem._number_conc =number_conc_at_i;
      for (ordinal_type i=0;i<m;++i)
        vals(i) = activeYs(i);

      real_type t = t_out_at_i(), dt = 0;
      cvode.initialize(t,
           dt_in, dt_min, dt_max,
           tol(0,0), tol(0,1),
           TChem::Impl::ProblemAerosolChemistry_ComputeFunctionCVODE,
           TChem::Impl::ProblemAerosolChemistry_ComputeJacobianCVODE);

      cvode.setProblem(problem);
      for (ordinal_type iter=0;iter<max_num_time_iterations && t<=t_end;++iter) {
        const real_type t_prev = t;
        cvode.advance(t_end, t, -1);
        dt = t - t_prev;
      }
      t_out_at_i() = t;
      dt_out_at_i() = dt;

      for (ordinal_type i=0;i<m;++i)
        Ys_out(i) = vals(i);
      temperature_out() =temperature;

    }
  }
      });
#endif
    Kokkos::Profiling::popRegion();
  }



#define TCHEM_RUN_AEROSOL_CHEMISTRY()     \
  AerosolChemistry_CVODE_TemplateRun(          \
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
          state_out,        \
          kmcds,        \
          amcds,        \
          cvodes)

  void
  AerosolChemistry_CVODE::runHostBatch( /// thread block size
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
           const Tines::value_type_1d_view<KineticModelNCAR_ConstData<interf_host_device_type>, interf_host_device_type>& kmcds,
           const Tines::value_type_1d_view<AerosolModel_ConstData<interf_host_device_type>, interf_host_device_type>& amcds,
           const Tines::value_type_1d_view<Tines::TimeIntegratorCVODE<real_type, host_device_type>, host_device_type>& cvodes) {
#if defined(TINES_ENABLE_TPL_SUNDIALS)
    const std::string profile_name = "TChem::AerosolChemistry::runHostBatch::kmcd array";
    // Note: we do not support SACADO and CVODE. Thus, CVODE uses a numerical Jacobian.
    using value_type = real_type;
    TCHEM_RUN_AEROSOL_CHEMISTRY();
#endif
  }// namespace TChem
}
