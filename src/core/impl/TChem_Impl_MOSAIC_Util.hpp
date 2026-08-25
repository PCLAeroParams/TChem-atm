#ifndef __TCHEM_IMPL_MOSAIC_UTIL_HPP__
#define __TCHEM_IMPL_MOSAIC_UTIL_HPP__

  KOKKOS_INLINE_FUNCTION static
  void update_thermodynamic_constants(const MosaicModelData<DeviceType>& mosaic,
                                      const real_type& T_K,
                                      real_type& aH2O,
                                      const real_type_2d_view_type& log_gamZ,
                                      const real_type_1d_view_type& Keq_gl,
                                      const real_type_1d_view_type& Keq_ll,
                                      real_type& Kp_nh3,
                                      real_type& Kp_nh4no3,
                                      real_type& Kp_nh4cl,
                                      const real_type_1d_view_type& Keq_sg,
                                      const real_type_1d_view_type& Keq_sl,
                                      const real_type_1d_view_type& Po_soa,
                                      const real_type_1d_view_type& sat_soa,
                                      real_type& sigma_water,
                                      const real_type_1d_view_type& MDRH_T,
                                      real_type& Kp_nh4no3_0,
                                      real_type& Kp_nh4cl_0) {

    const real_type tr = 298.15;        // reference temperature
    const real_type rt = 82.056*T_K/(1.0e9*1.0e6); // [m^3 atm/nmol]
    real_type K_eq = 0.0;

    // gas-liquid
    Keq_gl(0)= 1.0;                // Kelvin Effect (default)
    fn_Keq(57.64, 13.79, -5.39,T_K, K_eq);     // NH3(g)  <=> NH3(l)
    Keq_gl(1) = K_eq*rt;
    fn_Keq(2.63e6, 29.17, 16.83,T_K, K_eq);     // HNO3(g) <=> NO3- + H+
    Keq_gl(2) = K_eq*rt;
    fn_Keq(2.00e6, 30.20, 19.91,T_K, K_eq);     // HCl(g)  <=> Cl- + H+
    Keq_gl(3) = K_eq*rt;

    // liquid-liquid
    fn_Keq(1.0502e-2, 8.85, 25.14, T_K, K_eq);      // HSO4- <=> SO4= + H+
    Keq_ll(0) = K_eq;
    fn_Keq(1.805e-5, -1.50, 26.92, T_K, K_eq);      // NH3(l) + H2O = NH4+ + OH-
    Keq_ll(1) = K_eq;
    fn_Keq(1.01e-14,-22.52, 26.92, T_K, K_eq);      // H2O(l) <=> H+ + OH-
    Keq_ll(2) = K_eq;

    Kp_nh3    = Keq_ll(2)/(Keq_ll(1)*Keq_gl(1));
    Kp_nh4no3 = Kp_nh3/Keq_gl(2);
    Kp_nh4cl  = Kp_nh3/Keq_gl(3);

    // solid-gas
    fn_Keq(4.72e-17, -74.38, 6.12, T_K, K_eq); // NH4NO3<=>NH3(g)+HNO3(g)
    Keq_sg(0) = K_eq/ats<real_type>::pow(rt, 2.0);
    fn_Keq(8.43e-17, -71.0, 2.40, T_K, K_eq);  // NH4Cl <=>NH3(g)+HCl(g)
    Keq_sg(1) = K_eq/ats<real_type>::pow(rt, 2.0);


    // solid-liquid
    fn_Keq(1.04,-2.65, 38.57, T_K, K_eq);  // amSO4(s) = 2NH4+ + SO4=
    Keq_sl(mosaic.jnh4so4) = K_eq;
    fn_Keq(11.8, -5.19, 54.4, T_K, K_eq);   // lvcite(s)= 3NH4+ + HSO4- + SO4=
    Keq_sl(mosaic.jlvcite) = K_eq;
    fn_Keq(117.0,-2.87, 15.83, T_K, K_eq);  // amHSO4(s)= NH4+ + HSO4-
    Keq_sl(mosaic.jnh4hso4)= K_eq;
    Keq_sl(mosaic.jnh4msa) = 1.0e15;        // NH4MSA(s)= NH4+ + MSA-
    fn_Keq(12.21,-10.4, 17.56, T_K, K_eq);  // NH4NO3(s)= NH4+ + NO3-
    Keq_sl(mosaic.jnh4no3) = K_eq;
    fn_Keq(17.37,-6.03, 16.92, T_K, K_eq);  // NH4Cl(s) = NH4+ + Cl-
    Keq_sl(mosaic.jnh4cl)  = K_eq;
    fn_Keq(0.491, 0.98, 39.75, T_K, K_eq);  // Na2SO4(s)= 2Na+ + SO4=
    Keq_sl(mosaic.jna2so4) = K_eq;
    fn_Keq(313.0, 0.8,  14.79, T_K, K_eq);  // NaHSO4(s)= Na+ + HSO4-
    Keq_sl(mosaic.jnahso4) = K_eq;
    Keq_sl(mosaic.jna3hso4)= 1.0e15;        // Na3H(SO4)2(s) = 2Na+ + HSO4- + SO4=
    Keq_sl(mosaic.jnamsa)  = 1.0e15;        // NaMSA(s) = Na+ + MSA-
    fn_Keq(11.95,-8.22, 16.01, T_K, K_eq);  // NaNO3(s) = Na+ + NO3-
    Keq_sl(mosaic.jnano3)  = K_eq;
    fn_Keq(38.28,-1.52, 16.89, T_K, K_eq);  // NaCl(s)  = Na+ + Cl-
    Keq_sl(mosaic.jnacl)   = K_eq;
    fn_Keq(8.0e11, 32.84,44.79, T_K, K_eq); // CaCl2(s) = Ca++ + 2Cl-
    Keq_sl(mosaic.jcacl2)  = K_eq;
    fn_Keq(4.31e5, 7.83,42.01, T_K, K_eq);  // Ca(NO3)2(s) = Ca++ + 2NO3-
    Keq_sl(mosaic.jcano3)  = K_eq;
    Keq_sl(mosaic.jcamsa2) = 1.0e15;       // CaMSA2(s)= Ca+ + 2MSA-

    // vapor pressures of soa species
    real_type Po = 0.0;
    fn_Po(5.7e-5, 156.0, T_K, Po); // [Pascal]
    Po_soa(mosaic.iaro1_g) = Po;
    fn_Po(1.6e-3, 156.0, T_K, Po); // [Pascal]
    Po_soa(mosaic.iaro2_g) = Po;
    fn_Po(5.0e-6, 156.0, T_K, Po); // [Pascal]
    Po_soa(mosaic.ialk1_g) = Po;
    fn_Po(5.0e-6, 156.0, T_K, Po); // [Pascal]
    Po_soa(mosaic.iole1_g) = Po;
    fn_Po(4.0e-6, 156.0, T_K, Po); // [Pascal]
    Po_soa(mosaic.iapi1_g) = Po;
    fn_Po(1.7e-4, 156.0, T_K, Po); // [Pascal]
    Po_soa(mosaic.iapi2_g) = Po;
    fn_Po(2.5e-5, 156.0, T_K, Po); // [Pascal]
    Po_soa(mosaic.ilim1_g) = Po;
    fn_Po(1.2e-4, 156.0, T_K, Po); // [Pascal]
    Po_soa(mosaic.ilim2_g) = Po;

    real_type sat_factor = 0.5; // = 1.0 for original SORGAM parameters
    for (ordinal_type iv = 0; iv < mosaic.ngas_volatile; iv++) {
      sat_soa(iv) = sat_factor * 1.0e9*Po_soa(iv)/(RUNIV*T_K); // [nmol/m^3(air)]
    }

    // water surface tension
    real_type term = (647.15 - T_K)/647.15;
    sigma_water = 0.2358 * ats<real_type>::pow(term, 1.256) * (1.0 - 0.625*term); // surface tension of pure water in N/m

    // MDRH(T)
    real_type drh = 0.0;
    for (ordinal_type j_index = 0; j_index < 63; j_index++) {
      drh_mutual(mosaic, j_index, T_K, drh);
      MDRH_T(j_index) = drh;
    }

    MTEM_compute_log_gamZ(mosaic, aH2O, log_gamZ);

  // 6/25/2008 - start
    const real_type gam_nh4no3_0 = ats<real_type>::pow(10.0, log_gamZ(mosaic.jnh4no3, mosaic.jnh4no3));
    const real_type gam_nh4cl_0  = ats<real_type>::pow(10.0, log_gamZ(mosaic.jnh4cl, mosaic.jnh4cl));

    real_type molality = 0.0;
    molality_0(mosaic, mosaic.jnh4no3, aH2O, molality);
    const real_type m_nh4no3_0   = molality;
    molality_0(mosaic, mosaic.jnh4cl, aH2O, molality);
    const real_type m_nh4cl_0    = molality;

    Kp_nh4no3_0  = Kp_nh4no3*ats<real_type>::pow((m_nh4no3_0*gam_nh4no3_0), 2.0);
    Kp_nh4cl_0   = Kp_nh4cl*ats<real_type>::pow((m_nh4cl_0*gam_nh4cl_0), 2.0);
  // 6/25/2008 - end
  } // update_thermodynamic_constants

  KOKKOS_INLINE_FUNCTION static
  void check_aerosol_mass(const MosaicModelData<DeviceType>& mosaic,
                          const real_type_1d_view_type& aer_total,
                          real_type& mass_dry_a,
                          real_type& jphase,
                          real_type& jaerosolstate,
                          real_type& num_a) {

    mass_dry_a = 0.0;

    real_type aer_H = (2.0*aer_total(mosaic.iso4_a)  +
                           aer_total(mosaic.ino3_a)  +
                           aer_total(mosaic.icl_a)   +
                           aer_total(mosaic.imsa_a)  +
                       2.0*aer_total(mosaic.ico3_a)) -
                      (2.0*aer_total(mosaic.ica_a)   +
                           aer_total(mosaic.ina_a)   +
                           aer_total(mosaic.inh4_a));
    aer_H = max(aer_H, 0.0);

    auto mw_aer_mac = mosaic.mw_aer_mac.template view<DeviceType>();
    for (ordinal_type iaer = 0; iaer < mosaic.naer; iaer++) {
      mass_dry_a = mass_dry_a + aer_total(iaer)*mw_aer_mac(iaer); // ng/m^3(air)
    }
    mass_dry_a = mass_dry_a + aer_H;

    const real_type drymass = mass_dry_a; // ng/m^3(air)
    mass_dry_a = mass_dry_a*1.0e-15; // g/cc(air)

    if (drymass < mosaic.mass_cutoff) {
      jaerosolstate = mosaic.no_aerosol;
      jphase = 0;
      if (drymass == 0.0) {
        num_a = 0.0;
      }
    }
  } // check_aerosol_mass

  KOKKOS_INLINE_FUNCTION static
  void calculate_XT(const MosaicModelData<DeviceType>& mosaic,
                    const real_type_1d_view_type& aer,
                    real_type &XT) {

    //aer nmol/m^3
    if ((aer(mosaic.iso4_a) + aer(mosaic.imsa_a)) > 0.0) {
        XT = (aer(mosaic.inh4_a) + aer(mosaic.ina_a) +
             2.0 * aer(mosaic.ica_a)) /
            (aer(mosaic.iso4_a) + 0.5 * aer(mosaic.imsa_a));
    } else {
        XT = -1.0;
    }
  } // calculate_XT

  KOKKOS_INLINE_FUNCTION static
  void aerosol_water(const MosaicModelData<DeviceType>& mosaic,
                     const real_type_1d_view_type& electrolyte,
                     const real_type& aH2O_a,
                     const real_type_1d_view_type& molalities,
                     real_type& jaerosolstate,
                     real_type& jphase,
                     real_type& jhyst_leg,
                     real_type& aerosol_water) {

    for (ordinal_type je = 0; je < mosaic.nelectrolyte; je++) {
      real_type molality = 0.0;
      bin_molality(mosaic, je, aH2O_a, molality); // compute aH2O dependent binary molalities  EFFI
      molalities(je) = molality;
    }

    real_type dum = 0.0;
    for (ordinal_type je = 0; je < (mosaic.nsalt + 4); je++) { // include hno3 and hcl in water calculation
      dum += electrolyte(je) / molalities(je);
    }

    aerosol_water = dum * 1.0e-9;
    if (aerosol_water <= 0.0) {
      jaerosolstate = mosaic.all_solid;
      jphase = mosaic.jsolid;
      jhyst_leg = mosaic.jhyst_lo;
    }
  } // aerosol_water

  KOKKOS_INLINE_FUNCTION static
  void fnlog_gamZ(const MosaicModelData<DeviceType>& mosaic,
                  const ordinal_type& jA,
                  const ordinal_type& jE,
                  const real_type& aH2O,
                  real_type& log_gamZ_) {

    // FIXME: aH2O should not be local; make sure updated with RH

    auto b_mtem = mosaic.b_mtem.template view<DeviceType>();
    auto aw_min = mosaic.aw_min.template view<DeviceType>();

    const real_type aw = max(aH2O, aw_min(jE));

    log_gamZ_ = b_mtem(0,jA,jE) + aw *
               (b_mtem(1,jA,jE) + aw *
               (b_mtem(2,jA,jE) + aw *
               (b_mtem(3,jA,jE) + aw *
               (b_mtem(4,jA,jE) + aw *
                b_mtem(5,jA,jE) ))));
  } // fnlog_gamZ

  KOKKOS_INLINE_FUNCTION static
  void fn_Keq(const real_type& Keq_298,
              const real_type& a,
              const real_type& b,
              const real_type& T,
              real_type& Keq) {

    real_type tt = 298.15/T;

    Keq = Keq_298*ats<real_type>::exp(a*(tt - 1.0) + b*(1.0 + ats<real_type>::log(tt) - tt));
  } // fn_Keq

  KOKKOS_INLINE_FUNCTION static
  void fn_Po(const real_type& Po_298,
             const real_type& DH,
             const real_type& T_K,
             real_type& Po) {

    // Van't Hoff Equation
    Po = Po_298*ats<real_type>::exp(-(DH/(RUNIV/1000))*(1.0/T_K - (1/298.15)));
  } // fn_Po

  KOKKOS_INLINE_FUNCTION static
  void drh_mutual(const MosaicModelData<DeviceType>& mosaic,
                  const ordinal_type& j_index,
                  const real_type& T_K,
                  real_type& drh_mutual) {

    if (j_index == 6 || j_index == 7 || (j_index >= 33 && j_index <= 50)) {
      drh_mutual = 10.0; // CaNO3 or CaCl2 containing mixtures
    } else {
      auto d_mdrh = mosaic.d_mdrh.template view<DeviceType>();
      drh_mutual = d_mdrh(j_index,0) + T_K *
                  (d_mdrh(j_index,1) + T_K *
                  (d_mdrh(j_index,2) + T_K *
                   d_mdrh(j_index,3) )) + 1.0;
    }
  } // drh_mutual

  KOKKOS_INLINE_FUNCTION static
  void molality_0(const MosaicModelData<DeviceType>& mosaic,
                  const ordinal_type& je,
                  real_type& aw,
                  real_type& molality) {

    auto aw_min = mosaic.aw_min.template view<DeviceType>();
    auto a_zsr = mosaic.a_zsr.template view<DeviceType>();

    aw = max(aw, aw_min(je));
    aw = min(aw, 0.999999);

    if (aw < 0.97) {

      real_type xm = a_zsr(0,je) +
                 aw*(a_zsr(1,je) +
                 aw*(a_zsr(2,je) +
                 aw*(a_zsr(3,je) +
                 aw*(a_zsr(4,je) +
                 aw* a_zsr(5,je) ))));

        molality = 55.509*xm/(1.0 - xm);
    } else {
      auto b_zsr = mosaic.b_zsr.template view<DeviceType>();
      molality = -b_zsr(je)*ats<real_type>::log(aw);
    }
  } // molality_0

    KOKKOS_INLINE_FUNCTION static
  void bin_molality(const MosaicModelData<DeviceType>& mosaic,
                    const ordinal_type& je,
                    const real_type& aH2O_a,
                    real_type& molality) {

    auto aw_min = mosaic.aw_min.template view<DeviceType>();
    auto a_zsr = mosaic.a_zsr.template view<DeviceType>();

    real_type aw = max(aH2O_a, aw_min(je));
    aw = min(aw, 0.999999);

    if (aw < 0.97) {

      real_type xm = a_zsr(0,je) +
                 aw*(a_zsr(1,je) +
                 aw*(a_zsr(2,je) +
                 aw*(a_zsr(3,je) +
                 aw*(a_zsr(4,je) +
                 aw* a_zsr(5,je) ))));

        molality = 55.509*xm/(1.0 - xm);
    } else {
      auto b_zsr = mosaic.b_zsr.template view<DeviceType>();
      molality = -b_zsr(je)*ats<real_type>::log(aw);
    }
  } // bin_molality

  KOKKOS_INLINE_FUNCTION static
  void fuchs_sutugin(const real_type& rkn,
                     const real_type& a,
                     real_type& fuchs_sut) {

    const real_type rnum  = 0.75*a*(1. + rkn);
    const real_type denom = ats<real_type>::pow(rkn, 2.0) + rkn + 0.283*rkn*a + 0.75*a;
    fuchs_sut = rnum/denom;
  } // fuchs_sutugin

  KOKKOS_INLINE_FUNCTION static
  void gas_diffusivity(const real_type& T_K,
                       const real_type& P,
                       const real_type& MW,
                       const real_type& Vm,
                       real_type& gas_diff) {
    gas_diff = (1.0e-3 * ats<real_type>::pow(T_K, 1.75) * ats<real_type>::sqrt(1.0/MW + 0.035))/
                    (P * ats<real_type>::pow(ats<real_type>::pow(Vm, 0.333333) + 2.7189, 2.0));
  } // gas_diffusivity

  KOKKOS_INLINE_FUNCTION static
  void mean_molecular_speed(const real_type& T_K,
                            const real_type& MW,
                            real_type& mean_molec_speed) {

    mean_molec_speed = 1.455e4 * ats<real_type>::sqrt(T_K/MW);
  } // mean_molecular_speed

  KOKKOS_INLINE_FUNCTION static
  void bin_molality_60(const MosaicModelData<DeviceType>& mosaic,
                    const ordinal_type& je,
                    real_type& molality) {

    auto a_zsr = mosaic.a_zsr.template view<DeviceType>();

    const real_type aw = 0.6;

    real_type xm = a_zsr(0,je) +
               aw*(a_zsr(1,je) +
               aw*(a_zsr(2,je) +
               aw*(a_zsr(3,je) +
               aw*(a_zsr(4,je) +
               aw* a_zsr(5,je) ))));

    molality = 55.509*xm/(1.0 - xm);
  } // bin_molality_60

  KOKKOS_INLINE_FUNCTION static
  void quadratic(const real_type& a,
                 const real_type& b,
                 const real_type& c,
                 real_type& result) {
  
    using ats = Kokkos::ArithTraits<real_type>;
    real_type x;
    if (b != 0.0) {
      x = 4.0 * (a/b) * (c/b);
    } else {
      x = 1.0e+6;
    }

    if (ats::abs(x) < 1.0e-6) {
      real_type dum = (0.5*x) + (0.125*x*x) + (0.0625*x*x*x);
      result = (-0.5*b/a) * dum;
      if (result < 0.0) {
        result = -b/a - result;
      }
    } else {
      real_type disc = b*b - 4.0*a*c;
      real_type sqrtdisc = ats::sqrt(disc);
      real_type quad1 = ((-b) + sqrtdisc) / (2.0*a);
      real_type quad2 = ((-b) - sqrtdisc) / (2.0*a);
      result = max(quad1, quad2);
    }
  } // quadratic

  KOKKOS_INLINE_FUNCTION static
  void calc_dry_n_wet_aerosol_props(const MosaicModelData<DeviceType>& mosaic,
                                    const real_type_1d_view_type& aer_total,
                                    const real_type  water_a,
                                    const real_type  num_a,
                                    const real_type  jaerosolstate,
                                    real_type& mass_dry_a,
                                    real_type& vol_dry_a,
                                    real_type& mass_wet_a,
                                    real_type& vol_wet_a,
                                    real_type& dens_dry_a,
                                    real_type& dens_wet_a,
                                    real_type& Dp_dry_a,
                                    real_type& Dp_wet_a,
                                    real_type& area_dry_a,
                                    real_type& area_wet_a) {

    auto mw_aer_mac   = mosaic.mw_aer_mac.template view<DeviceType>();
    auto dens_aer_mac = mosaic.dens_aer_mac.template view<DeviceType>();

    mass_dry_a = 0.0;
    vol_dry_a  = 0.0;
    area_dry_a = 0.0;

    if (static_cast<ordinal_type>(jaerosolstate) != mosaic.no_aerosol) {

      real_type aer_H = (2.0 * aer_total(mosaic.iso4_a) +
                               aer_total(mosaic.ino3_a)  +
                               aer_total(mosaic.icl_a)   +
                               aer_total(mosaic.imsa_a)  +
                         2.0 * aer_total(mosaic.ico3_a)) -
                        (2.0 * aer_total(mosaic.ica_a)  +
                               aer_total(mosaic.ina_a)   +
                               aer_total(mosaic.inh4_a));
      aer_H = max(aer_H, (real_type)0.0);

      for (int iaer = 0; iaer < mosaic.naer; ++iaer) {
        mass_dry_a += aer_total(iaer) * mw_aer_mac(iaer);
        vol_dry_a  += aer_total(iaer) * mw_aer_mac(iaer) / dens_aer_mac(iaer);
      }
      mass_dry_a = (mass_dry_a + aer_H) * 1.e-15;
      vol_dry_a  = (vol_dry_a  + aer_H) * 1.e-15;

      mass_wet_a = mass_dry_a + water_a * 1.e-3;
      vol_wet_a  = vol_dry_a  + water_a * 1.e-3;

      dens_dry_a = mass_dry_a / vol_dry_a;
      dens_wet_a = mass_wet_a / vol_wet_a;

      Dp_dry_a = ats<real_type>::pow(1.90985 * vol_dry_a / num_a, 0.3333333);
      Dp_wet_a = ats<real_type>::pow(1.90985 * vol_wet_a / num_a, 0.3333333);

      area_dry_a = 0.785398 * num_a * Dp_dry_a * Dp_dry_a;
      area_wet_a = 0.785398 * num_a * Dp_wet_a * Dp_wet_a;

    } else {
      dens_dry_a = 1.0;
      dens_wet_a = 1.0;
    }
  } // calc_dry_n_wet_aerosol_props

#endif
