#ifndef __TCHEM_IMPL_SIMPOL_PHASE_TRANSFER_HPP__
#define __TCHEM_IMPL_SIMPOL_PHASE_TRANSFER_HPP__


/// Minimum aerosol-phase mass concentration [kg m-3]
#define MINIMUM_MASS_ 1.0e-25L
/// Minimum mass assumed molecular weight [kg mol-1]
#define MINIMUM_MW_ 0.1L
/// Minimum mass assumed density [kg m-3]
#define MINIMUM_DENSITY_ 1800.0L

#define CHEM_SPEC_UNKNOWN_TYPE 0
#define CHEM_SPEC_VARIABLE 1
#define CHEM_SPEC_CONSTANT 2
#define CHEM_SPEC_PSSA 3
#define CHEM_SPEC_ACTIVITY_COEFF 4

// FIXME: NUM_STATE_VAR_
#define NUM_STATE_VAR_ 3

#include "TChem_Util.hpp"
#include "TChem_Impl_SIMPOL_constant.hpp"

namespace TChem {
namespace Impl {
  template<typename ValueType, typename DeviceType>
struct SIMPOL_single_particle
{
  using value_type = ValueType;
  using device_type = DeviceType;
  using scalar_type = typename ats<value_type>::scalar_type;

  using real_type = scalar_type;
  /// sacado is value type
  using value_type_1d_view_type = Tines::value_type_1d_view<value_type,device_type>;
  using real_type_1d_view_type = Tines::value_type_1d_view<real_type,device_type>;
  using ordinal_type_1d_view_type = Tines::value_type_1d_view<ordinal_type,device_type>;

  KOKKOS_INLINE_FUNCTION static ordinal_type getWorkSpaceSize()
  {
    ordinal_type workspace_size=0;
    return workspace_size;
  }

  template<typename MemberType>
  KOKKOS_INLINE_FUNCTION static void
  aero_phase_get_volume__m3_m3(const MemberType& member,
                        const value_type_1d_view_type& state,
                        const real_type_1d_view_type& density,
                        const ordinal_type_1d_view_type& species_type,
                        value_type& volume)
   {
     // Sum the mass and MW
   volume = MINIMUM_MASS_ / MINIMUM_DENSITY_;
  for (int i_spec = 0; i_spec < NUM_STATE_VAR_; i_spec++) {
    if (species_type(i_spec) == CHEM_SPEC_VARIABLE ||
        species_type(i_spec) == CHEM_SPEC_CONSTANT ||
        species_type(i_spec) == CHEM_SPEC_PSSA) {
      volume += state(i_spec) / density(i_spec);
    }
  }

  }
  template<typename MemberType>
  KOKKOS_INLINE_FUNCTION static void effective_radius(
    const MemberType& member,
    const real_type& t,
    const real_type& p,
    const value_type_1d_view_type& state,
    const ordinal_type& num_of_phases,
    const real_type_1d_view_type& density,
    const ordinal_type_1d_view_type& species_type,
    value_type& radius)
    {
    radius=0.0;
    for (int i_phase = 0; i_phase < num_of_phases; ++i_phase) {
    value_type volume=0.0;
    aero_phase_get_volume__m3_m3(member, state, density, species_type, volume );
    radius += volume;
    }//
    radius = ats<value_type>::pow((radius * 3.0 / 4.0 / PI()), 1.0 / 3.0);
    }
  template<typename MemberType>
  KOKKOS_INLINE_FUNCTION static void
  aero_phase_get_mass__kg_m3(const MemberType& member,
                             const value_type_1d_view_type& state,
                             const real_type_1d_view_type& molecular_weigths,
                             const ordinal_type_1d_view_type& species_type,
                             value_type&MW,
                             value_type&mass
                             )
  {
     // Sum the mass and MW
  mass = MINIMUM_MASS_;
  value_type moles = MINIMUM_MASS_ / MINIMUM_MW_;
  for (int i_spec = 0; i_spec < NUM_STATE_VAR_; i_spec++) {
    if (species_type(i_spec) == CHEM_SPEC_VARIABLE ||
        species_type(i_spec) == CHEM_SPEC_CONSTANT ||
        species_type(i_spec) == CHEM_SPEC_PSSA) {
      mass += state(i_spec);
      moles += state(i_spec) / molecular_weigths(i_spec);
    }
  }
  MW = mass / moles;
  }


/** Calculate the transition regime correction factor [unitless] \cite Fuchs1971
 * : \f[ f(K_n,\alpha) = \frac{0.75 \alpha ( 1 + K_n )}{K_n(1+K_n) + 0.283\alpha
 * K_n + 0.75 \alpha} \f] where the Knudsen Number \f$K_n = \lambda / r\f$
 * (where \f$\lambda\f$ is the mean free path [m] of the gas-phase species and
 * \f$r\f$ is the effective radius of the particles [m]), and \f$ \alpha \f$ is
 * the mass accomodation coefficient, which is typically assumed to equal 0.1
 * \cite Zaveri2008.
 *
 *  @param mean_free_path__m mean free path of the gas-phase species [m]
 *  @param radius__m Particle effective radius [m]
 *  @param alpha Mass accomodation coefficient [unitless]
 */

KOKKOS_INLINE_FUNCTION static
value_type
transition_regime_correction_factor(real_type mean_free_path__m,
                                    value_type radius__m,
                                    real_type alpha) {
  value_type K_n = mean_free_path__m / radius__m;
  return (0.75 * alpha * (1.0 + K_n)) /
         (K_n * K_n + (1.0 + 0.283 * alpha) * K_n + 0.75 * alpha);
}

/** Calculate the gas-aerosol reaction rate constant for the transition regime
 * [\f$\mbox{m}^3\, \mbox{particle}^{-1}\, \mbox{s}^{-1}\f$]
 *
 *  The rate constant \f$k_c\f$ is calculated according to \cite Zaveri2008 as:
 *  \f[
 *    k_c = 4 \pi r D_g f_{fs}( K_n, \alpha )
 *  \f]
 *  where \f$r\f$ is the radius of the particle(s) [m], \f$D_g\f$ is the
 * diffusion coefficient of the gas-phase species
 * [\f$\mbox{m}^2\,\mbox{s}^{-1}\f$], \f$f_{fs}( K_n, \alpha )\f$ is the Fuchs
 * Sutugin transition regime correction factor [unitless], \f$K_n\f$ is the
 * Knudsen Number [unitess], and \f$\alpha\f$ is the mass accomodation
 * coefficient.
 *
 *  Rates can be calculated as:
 *  \f[
 *    r_c = [G] N_a k_c
 *  \f]
 *  where \f$[G]\f$ is the gas-phase species concentration [ppm], \f$N_a\f$ is
 * the number concentration of particles [\f$\mbox{particle}\,\mbox{m}^{-3}\f$]
 * and the rate \f$r_c\f$ is in [\f$\mbox{ppm}\,\mbox{s}^{-1}\f$].
 *
 *  @param diffusion_coeff__m2_s Diffusion coefficent of the gas species
 *  [\f$\mbox{m}^2\, \mbox{s}^{-1}\f$]
 *  @param mean_free_path__m Mean free path of gas molecules [m]
 *  @param radius__m Particle radius [m]
 *  @param alpha Mass accomodation coefficient [unitless]
 */
KOKKOS_INLINE_FUNCTION
static real_type gas_aerosol_transition_rxn_rate_constant(
    real_type diffusion_coeff__m2_s, real_type mean_free_path__m,
    real_type radius__m, real_type alpha) {
  return 4.0 * PI() * radius__m * diffusion_coeff__m2_s *
         transition_regime_correction_factor(mean_free_path__m, radius__m,
                                             alpha);
}

  // FIXME:
KOKKOS_INLINE_FUNCTION
static void get_inputParameters(ordinal_type& num_of_phases,
                      ordinal_type& GAS_SPEC_,
                      ordinal_type& AERO_SPEC_i_phase)
   {
    num_of_phases =1;
    GAS_SPEC_=0;
    AERO_SPEC_i_phase=1;
   }

  // update RHS
  template<typename MemberType>
  KOKKOS_INLINE_FUNCTION static
  void invoke_team(const MemberType& member,
    const real_type& t,
    const real_type& p,
    const real_type& number_conc,
    const value_type_1d_view_type& state,
    const real_type_1d_view_type& molecular_weigths,
    const real_type_1d_view_type& density,
    const ordinal_type_1d_view_type& species_type,
    const value_type_1d_view_type& omega
    )
  {
    // FIXME: input parameters
  ordinal_type GAS_SPEC_=0;
  ordinal_type AERO_SPEC_i_phase=1;
  ordinal_type num_of_phases =1;
  get_inputParameters(num_of_phases, GAS_SPEC_, AERO_SPEC_i_phase);

  using SIMPOL_constant_type
   = TChem::Impl::SIMPOL_constant<value_type, device_type >;

  //FIXME: inputs
  real_type DIFF_COEFF_=1;
  real_type MW_=1;
  SIMPOL_constant_type::getInputParameters(DIFF_COEFF_, MW_ );

     // aerosol phase (m3/#/s)
  real_type mfp_m =0.0;
  value_type alpha= 0.0;
  real_type EQUIL_CONST_=0.0;
  real_type KGM3_TO_PPM_=0.0;
  SIMPOL_constant_type::team_invoke( member, t, p, alpha,
    mfp_m, KGM3_TO_PPM_, EQUIL_CONST_);

  // Compute radious
  value_type radius =1;
  effective_radius(member, t, p, state,
                  num_of_phases, density, species_type, radius);
  printf("radius %e \n", radius);

   // Get the particle number concentration (#/m3)
   // Get the total mass of the aerosol phase (kg/mol)
   value_type aero_phase_avg_MW=0;
   // Get the total mass of the aerosol phase (kg/m3)
   value_type aero_phase_mass=0;
   aero_phase_get_mass__kg_m3(member, state, molecular_weigths, species_type,
                              aero_phase_avg_MW, aero_phase_mass );

  printf("aero_phase_mass %e \n", aero_phase_mass);
  printf("aero_phase_avg_MW %e \n", aero_phase_avg_MW);

   // Calculate the rate constant for diffusion limited mass transfer to the

   value_type cond_rate =
   gas_aerosol_transition_rxn_rate_constant(DIFF_COEFF_, mfp_m, radius, alpha);

   printf("cond_rate %e \n", cond_rate);


  //FIXME: to be done
  // Get the activity coefficient (if one exists)
  value_type act_coeff = 1.0;
  // if (AERO_ACT_ID_(i_phase) > -1) {
  //   act_coeff = state[AERO_ACT_ID_(i_phase)];
  // }
  // Calculate the evaporation rate constant (ppm_x*m^3/kg_x/s)
  value_type evap_rate =
        cond_rate *  act_coeff * (EQUIL_CONST_ * aero_phase_avg_MW / aero_phase_mass);

  printf("act_coeff %e \n", act_coeff);
  printf("EQUIL_CONST_ %e \n", EQUIL_CONST_);
  printf("evap_rate %e \n", evap_rate);

  // Calculate the evaporation and condensation rates
  cond_rate *= state(GAS_SPEC_);
  evap_rate *= state(AERO_SPEC_i_phase);
  // // per-particle mass concentration rates
  value_type diff = - evap_rate + cond_rate;

  printf("state(GAS_SPEC_) %e \n", state(GAS_SPEC_));
  printf("state(%d) %e \n",AERO_SPEC_i_phase, state(AERO_SPEC_i_phase));
  printf("diff %e \n",diff);
  printf("A evap_rate %e \n", evap_rate);
  printf("A cond_rate %e \n", cond_rate);

  // Change in the gas-phase is evaporation - condensation (ppm/s)
  omega(GAS_SPEC_) = -number_conc * diff;
  omega(AERO_SPEC_i_phase) = diff / KGM3_TO_PPM_;

  for (int i = 0; i < 3; i++)
  {
    printf("omega(%d) %e \n",i,omega(i));
  }

#if 0
#endif
  }

};

} // namespace Impl
} // namespace TChem

#endif
