#ifndef __TCHEM_IMPL_SIMPOL_CONSTANT_HPP__
#define __TCHEM_IMPL_SIMPOL_CONSTANT_HPP__

#include "TChem_Util.hpp"
namespace TChem {
namespace Impl {
  template<typename ValueType, typename DeviceType>
struct SIMPOL_constant
{
  using value_type = ValueType;
  using device_type = DeviceType;
  using scalar_type = typename ats<value_type>::scalar_type;

  using real_type = scalar_type;
  /// sacado is value type
  using value_type_1d_view_type = Tines::value_type_1d_view<value_type,device_type>;
  using real_type_1d_view_type = Tines::value_type_1d_view<real_type,device_type>;

  using aerosol_model_data_type= AerosolModel_ConstData<device_type>;

  KOKKOS_INLINE_FUNCTION static ordinal_type getWorkSpaceSize()
  {
    ordinal_type workspace_size=0;
    return workspace_size;
  }

//FIXME: move to util?
/** Calculate the mean speed of a gas-phase species
 * [ \f$\mbox{m}\,\mbox{s}^{-1}\f$ ] :
 * \f[
 *   v = \sqrt{\frac{8RT}{\pi MW}}
 * \f]
 * where R is the ideal gas constant [\f$\mbox{J}\, \mbox{K}^{-1}\,
 * \mbox{mol}^{-1}\f$], T is temperature [K], and MW is the molecular weight of
 * the gas-phase species
 * [\f$\mbox{kg}\, \mbox{mol}^{-1}\f$]
 * @param temperature__K Temperature [K]
 * @param mw__kg_mol Molecular weight of the gas-phase species [\f$\mbox{kg}\,
 * \mbox{mol}^{-1}\f$]
 */
 KOKKOS_INLINE_FUNCTION static
 real_type mean_speed__m_s(const real_type t,
                          const real_type mw) {
  return sqrt(8.0 * RUNIV * t / (PI() * mw));
}

  /** Calculate the mean free path of a gas-phase species [m]
 * \f[
 *   \lambda = 3.0 D_g / v
 * \f]
 * where \f$D_g\f$ is the gas-phase diffusion coefficient
 * [\f$\mbox{m}^2\,\mbox{s}^{-1}\f$] and
 * \f$v\f$ is the mean speed of the gas-phase molecules
 * [ \f$\mbox{m}\,\mbox{s}^{-1}\f$ ].
 *
 * @param diffusion_coeff__m2_s Diffusion coefficient of the gas species
 * [\f$\mbox{m}^2\, \mbox{s}^{-1}\f$]
 * @param temperature__K Temperature [K]
 * @param mw__kg_mol Molecular weight of the gas-phase species [\f$\mbox{kg}\,
 * \mbox{mol}^{-1}\f$]
 */
 KOKKOS_INLINE_FUNCTION static real_type mean_free_path_m(const real_type diffusion_coeff__m2_s,
                                       const real_type t,
                                       const real_type mw) {
  return 3.0 * diffusion_coeff__m2_s /
         mean_speed__m_s(t, mw);
}

    template<typename MemberType>
  KOKKOS_INLINE_FUNCTION static void team_invoke(
    const MemberType& member,
    const real_type& t,
    const real_type& p,
    real_type& alpha,
    real_type& mfp_m,
    real_type& KGM3_TO_PPM_,
    real_type& equil_constant,
    const SIMPOL_PhaseTransferType& simpol_params
    )
    {

  real_type DIFF_COEFF_=simpol_params.diffusion_coeff;
  real_type MW_=simpol_params.molecular_weight;

  real_type N_star=simpol_params.N_star;
  real_type B1_=simpol_params.B1;
  real_type B2_=simpol_params.B2;
  real_type B3_=simpol_params.B3;
  real_type B4_=simpol_params.B4;
  bool compute_alpha=simpol_params.compute_alpha;

  alpha = 0.1;

  if (compute_alpha) {
    /*
    ! Mass accomodation equation is based on equations in:
    ! Ervens, B., et al., 2003. "CAPRAM 2.4 (MODAC mechanism): An extended
    ! and condensed tropospheric aqueous mechanism and its application."
    ! J. Geophys. Res. 108, 4426. doi:10.1029/2002JD002202
    */
    // enthalpy change (kcal mol-1)
    real_type N_star_2_3 = ats<real_type>::pow(N_star,2.0/3.0);
    real_type delta_h = - 10.0*(N_star-1.0) + 7.53*(N_star_2_3-1.0) - 1.0;
    // entropy change (cal mol-1)
    real_type delta_s = - 32.0*(N_star-1.0) +
              9.21*(N_star_2_3-1.0) - 1.3;
    // Convert dH and dS to (J mol-1)
    delta_h = delta_h * 4184.0;
    delta_s = delta_s * 4.184;

    real_type del_G = delta_h - t * delta_s;
    alpha =  ats<real_type>::exp(-del_G / (RUNIV * t));
    alpha = alpha / (1.0 + alpha);
  }

  /// save the mean free path [m] for calculating condensation rates
  mfp_m = mean_free_path_m(DIFF_COEFF_, t, MW_);

  // SIMPOL.1 vapor pressure [Pa]
  real_type vp = B1_ / t + B2_ + B3_ * t +
              B4_ * ats<real_type>::log(t);
  vp = ATMPA() * ats<real_type>::pow(10, vp);

  /* Set the kg/m3 -> ppm conversion prefactor (multiply by T/P to get
    ! conversion)
    ! (ppm_x*Pa_air*m^3/K/kg_x) = Pa_air*m^3/mol_air/K * mol_x/kg_x *
    !                   1.0e6ppm_x*mol_air/mol_x */
  const real_type CONV = RUNIV / MW_ * 1.0e6;
  // Calculate the conversion from kg_x/m^3 -> ppm_x
  KGM3_TO_PPM_ = CONV * t / p;

    // Calculate the partitioning coefficient K_eq (ppm_x/kg_x*kg_tot)
  // such that for partitioning species X at equilibrium:
  //   [X]_gas = [X]_aero * activity_coeff_X * K_eq * MW_tot_aero / [tot]_aero
  // where 'tot' indicates all species within an aerosol phase combined
  // with []_gas in (ppm) and []_aero in (kg/m^3)
  equil_constant  = vp              // Pa_x / (mol_x_aero/mol_tot_aero)
                 / p  // 1/Pa_air
                 / MW_           // mol_x_aero/kg_x_aero
                 * 1.0e6;        // ppm_x / (Pa_x/P_air)
                                 // = ppm_x * mol_tot_aero / kg_x_aero
  }


};

} // namespace Impl
} // namespace TChem

#endif
