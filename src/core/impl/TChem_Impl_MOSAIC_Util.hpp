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
  void MTEM_compute_log_gamZ(const MosaicModelData<DeviceType>& mosaic,
                             const real_type& aH2O,
                             real_type_2d_view_type log_gamZ) {

    real_type log_gamZ_ = 0.0;

  // sulfate-poor species
    ordinal_type jA = mosaic.jhno3;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4no3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4no3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4cl, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4cl) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jna2so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jna2so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnano3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnano3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnacl, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnacl) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jcano3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jcano3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jcacl2, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jcacl2) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jhno3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jhno3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jhcl, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jhcl) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jh2so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jh2so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4hso4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4hso4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jlvcite, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jlvcite) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnahso4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnahso4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jna3hso4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jna3hso4) = log_gamZ_;

    jA = mosaic.jhcl;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4no3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4no3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4cl, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4cl) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jna2so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jna2so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnano3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnano3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnacl, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnacl) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jcano3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jcano3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jcacl2, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jcacl2) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jhno3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jhno3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jhcl, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jhcl) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jh2so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jh2so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4hso4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4hso4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jlvcite, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jlvcite) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnahso4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnahso4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jna3hso4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jna3hso4) = log_gamZ_;

    jA = mosaic.jnh4so4;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4no3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4no3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4cl, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4cl) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jna2so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jna2so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnano3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnano3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnacl, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnacl) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jcano3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jcano3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jcacl2, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jcacl2) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jhno3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jhno3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jhcl, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jhcl) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jh2so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jh2so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4hso4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4hso4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jlvcite, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jlvcite) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnahso4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnahso4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jna3hso4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jna3hso4) = log_gamZ_;

    jA = mosaic.jnh4no3;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4no3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4no3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4cl, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4cl) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jna2so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jna2so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnano3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnano3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnacl, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnacl) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jcano3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jcano3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jcacl2, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jcacl2) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jhno3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jhno3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jhcl, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jhcl) = log_gamZ_;

    jA = mosaic.jnh4cl;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4no3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4no3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4cl, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4cl) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jna2so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jna2so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnano3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnano3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnacl, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnacl) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jcano3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jcano3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jcacl2, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jcacl2) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jhno3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jhno3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jhcl, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jhcl) = log_gamZ_;

    jA = mosaic.jna2so4;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4no3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4no3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4cl, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4cl) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jna2so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jna2so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnano3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnano3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnacl, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnacl) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jcano3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jcano3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jcacl2, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jcacl2) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jhno3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jhno3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jhcl, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jhcl) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jh2so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jh2so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4hso4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4hso4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jlvcite, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jlvcite) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnahso4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnahso4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jna3hso4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jna3hso4) = log_gamZ_;

    jA = mosaic.jnano3;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4no3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4no3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4cl, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4cl) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jna2so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jna2so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnano3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnano3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnacl, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnacl) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jcano3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jcano3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jcacl2, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jcacl2) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jhno3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jhno3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jhcl, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jhcl) = log_gamZ_;

    jA = mosaic.jnacl;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4no3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4no3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4cl, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4cl) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jna2so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jna2so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnano3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnano3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnacl, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnacl) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jcano3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jcano3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jcacl2, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jcacl2) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jhno3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jhno3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jhcl, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jhcl) = log_gamZ_;

    jA = mosaic.jcano3;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4no3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4no3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4cl, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4cl) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jna2so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jna2so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnano3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnano3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnacl, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnacl) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jcano3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jcano3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jcacl2, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jcacl2) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jhno3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jhno3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jhcl, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jhcl) = log_gamZ_;

    jA = mosaic.jcacl2;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4no3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4no3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4cl, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4cl) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jna2so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jna2so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnano3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnano3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnacl, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnacl) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jcano3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jcano3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jcacl2, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jcacl2) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jhno3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jhno3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jhcl, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jhcl) = log_gamZ_;

  // sulfate-rich species
    jA = mosaic.jh2so4;
    fnlog_gamZ(mosaic, jA, mosaic.jh2so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jh2so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4hso4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4hso4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jlvcite, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jlvcite) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnahso4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnahso4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jna3hso4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jna3hso4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jna2so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jna2so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jhno3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jhno3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jhcl, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jhcl) = log_gamZ_;

    jA = mosaic.jhhso4;
    fnlog_gamZ(mosaic, jA, mosaic.jh2so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jh2so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4hso4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4hso4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jlvcite, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jlvcite) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnahso4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnahso4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jna3hso4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jna3hso4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jna2so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jna2so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jhno3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jhno3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jhcl, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jhcl) = log_gamZ_;

    jA = mosaic.jnh4hso4;
    fnlog_gamZ(mosaic, jA, mosaic.jh2so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jh2so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4hso4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4hso4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jlvcite, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jlvcite) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnahso4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnahso4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jna3hso4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jna3hso4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jna2so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jna2so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jhno3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jhno3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jhcl, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jhcl) = log_gamZ_;

    jA = mosaic.jlvcite;
    fnlog_gamZ(mosaic, jA, mosaic.jh2so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jh2so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4hso4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4hso4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jlvcite, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jlvcite) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnahso4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnahso4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jna3hso4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jna3hso4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jna2so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jna2so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jhno3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jhno3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jhcl, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jhcl) = log_gamZ_;

    jA = mosaic.jnahso4;
    fnlog_gamZ(mosaic, jA, mosaic.jh2so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jh2so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4hso4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4hso4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jlvcite, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jlvcite) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnahso4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnahso4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jna3hso4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jna3hso4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jna2so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jna2so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jhno3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jhno3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jhcl, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jhcl) = log_gamZ_;

    jA = mosaic.jna3hso4;
    fnlog_gamZ(mosaic, jA, mosaic.jh2so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jh2so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4hso4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4hso4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jlvcite, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jlvcite) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnh4so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnh4so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jnahso4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jnahso4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jna3hso4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jna3hso4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jna2so4, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jna2so4) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jhno3, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jhno3) = log_gamZ_;
    fnlog_gamZ(mosaic, jA, mosaic.jhcl, aH2O, log_gamZ_);
    log_gamZ(jA,mosaic.jhcl) = log_gamZ_;
  } // MTEM_compute_log_gamZ

  KOKKOS_INLINE_FUNCTION static
  void form_nahso4(const MosaicModelData<DeviceType>& mosaic,
                   const real_type_1d_view_type& electrolyte,
                   const real_type_1d_view_type& store) {

  electrolyte(mosaic.jnahso4) = min(store(mosaic.ina_a), store(mosaic.iso4_a));

  store(mosaic.ina_a)  = store(mosaic.ina_a)  - electrolyte(mosaic.jnahso4);
  store(mosaic.iso4_a) = store(mosaic.iso4_a) - electrolyte(mosaic.jnahso4);

  store(mosaic.ina_a)  = max(0.0, store(mosaic.ina_a));
  store(mosaic.iso4_a) = max(0.0, store(mosaic.iso4_a));
  } // form_nahso4
  
  KOKKOS_INLINE_FUNCTION static
  void form_na2so4(const MosaicModelData<DeviceType>& mosaic,
                   const real_type_1d_view_type& electrolyte,
                   const real_type_1d_view_type& store) {

    electrolyte(mosaic.jna2so4) = min(0.5*store(mosaic.ina_a), store(mosaic.iso4_a));

    store(mosaic.ina_a)  = store(mosaic.ina_a)  - 2.*electrolyte(mosaic.jna2so4);
    store(mosaic.iso4_a) = store(mosaic.iso4_a) - electrolyte(mosaic.jna2so4);

    store(mosaic.ina_a)  = max(0.0, store(mosaic.ina_a));
    store(mosaic.iso4_a) = max(0.0, store(mosaic.iso4_a));
  } // form_na2so4

  KOKKOS_INLINE_FUNCTION static
  void form_caco3(const MosaicModelData<DeviceType>& mosaic,
                  const real_type& jp,
                  const real_type_1d_view_type& electrolyte,
                  const real_type_1d_view_type& aer,
                  const real_type_1d_view_type& store) {
  
    if (jp == mosaic.jtotal || jp == mosaic.jsolid) {
      electrolyte(mosaic.jcaco3) = store(mosaic.ica_a);

      aer(mosaic.ico3_a) = electrolyte(mosaic.jcaco3); // force co3 = caco3

      store(mosaic.ica_a) = 0.0;
      store(mosaic.ico3_a) = 0.0;
    }
  } // form_caco3

  KOKKOS_INLINE_FUNCTION static
  void form_cacl2(const MosaicModelData<DeviceType>& mosaic,
                  const real_type_1d_view_type& electrolyte,
                  const real_type_1d_view_type& store) {

    electrolyte(mosaic.jcacl2) = min(store(mosaic.ica_a), 0.5*store(mosaic.icl_a));

    store(mosaic.ica_a)  = store(mosaic.ica_a)  - electrolyte(mosaic.jcacl2);
    store(mosaic.icl_a)  = store(mosaic.icl_a)  - 2.*electrolyte(mosaic.jcacl2);

    store(mosaic.ica_a)  = max(0.0, store(mosaic.ica_a));
    store(mosaic.icl_a)  = max(0.0, store(mosaic.icl_a));
  } // form_cacl2

  KOKKOS_INLINE_FUNCTION static
  void form_cano3(const MosaicModelData<DeviceType>& mosaic,
                  const real_type_1d_view_type& electrolyte,
                  const real_type_1d_view_type& store) {

    electrolyte(mosaic.jcano3) = min(store(mosaic.ica_a), 0.5*store(mosaic.ino3_a));

    store(mosaic.ica_a)  = store(mosaic.ica_a)  - electrolyte(mosaic.jcano3);
    store(mosaic.ino3_a) = store(mosaic.ino3_a) - 2.*electrolyte(mosaic.jcano3);

    store(mosaic.ica_a)  = max(0.0, store(mosaic.ica_a));
    store(mosaic.ino3_a) = max(0.0, store(mosaic.ino3_a));
  } // form_cano3

  KOKKOS_INLINE_FUNCTION static
  void form_camsa2(const MosaicModelData<DeviceType>& mosaic,
                   const real_type_1d_view_type& electrolyte,
                   const real_type_1d_view_type& store) {

    electrolyte(mosaic.jcaso4) = min(store(mosaic.ica_a), 0.5*store(mosaic.imsa_a));

    store(mosaic.ica_a)  = store(mosaic.ica_a)  - electrolyte(mosaic.jcamsa2);
    store(mosaic.imsa_a) = store(mosaic.imsa_a) - 2.*electrolyte(mosaic.jcamsa2);

    store(mosaic.ica_a)  = max(0.0, store(mosaic.ica_a));
    store(mosaic.imsa_a) = max(0.0, store(mosaic.imsa_a));
  } // form_camsa2

  KOKKOS_INLINE_FUNCTION static
  void form_caso4(const MosaicModelData<DeviceType>& mosaic,
                  const real_type_1d_view_type& electrolyte,
                  const real_type_1d_view_type& store) {

    electrolyte(mosaic.jcaso4) = min(store(mosaic.ica_a), store(mosaic.iso4_a));

    store(mosaic.ica_a)  = store(mosaic.ica_a)  - electrolyte(mosaic.jcaso4);
    store(mosaic.iso4_a) = store(mosaic.iso4_a) - electrolyte(mosaic.jcaso4);

    store(mosaic.ica_a)  = max(0.0, store(mosaic.ica_a));
    store(mosaic.iso4_a) = max(0.0, store(mosaic.iso4_a));
  } // form_caso4

  KOKKOS_INLINE_FUNCTION static
  void form_nh4msa(const MosaicModelData<DeviceType>& mosaic,
                   const real_type_1d_view_type& electrolyte,
                   const real_type_1d_view_type& store) {

    electrolyte(mosaic.jnh4msa) = min(store(mosaic.inh4_a), store(mosaic.imsa_a));

    store(mosaic.inh4_a) = store(mosaic.inh4_a) - electrolyte(mosaic.jnh4msa);
    store(mosaic.imsa_a) = store(mosaic.imsa_a) - electrolyte(mosaic.jnh4msa);

    store(mosaic.inh4_a) = max(0.0, store(mosaic.inh4_a));
    store(mosaic.imsa_a) = max(0.0, store(mosaic.imsa_a));
  } // form_nh4msa

  KOKKOS_INLINE_FUNCTION static
  void form_nh4so4(const MosaicModelData<DeviceType>& mosaic,
                   const real_type_1d_view_type& electrolyte,
                   const real_type_1d_view_type& store) {

    electrolyte(mosaic.jnh4so4) = min(0.5*store(mosaic.inh4_a), store(mosaic.iso4_a));

    store(mosaic.inh4_a)  = store(mosaic.inh4_a)  - 2.*electrolyte(mosaic.jnh4so4);
    store(mosaic.iso4_a) = store(mosaic.iso4_a) - electrolyte(mosaic.jnh4so4);

    store(mosaic.inh4_a)  = max(0.0, store(mosaic.inh4_a));
    store(mosaic.iso4_a) = max(0.0, store(mosaic.iso4_a));
  } // form_nh4so4

  KOKKOS_INLINE_FUNCTION static
  void form_nh4hso4(const MosaicModelData<DeviceType>& mosaic,
                   const real_type_1d_view_type& electrolyte,
                   const real_type_1d_view_type& store) {

    electrolyte(mosaic.jnh4hso4) = min(store(mosaic.inh4_a), store(mosaic.iso4_a));

    store(mosaic.inh4_a)  = store(mosaic.inh4_a) - electrolyte(mosaic.jnh4hso4);
    store(mosaic.iso4_a) = store(mosaic.iso4_a) - electrolyte(mosaic.jnh4hso4);

    store(mosaic.inh4_a)  = max(0.0, store(mosaic.inh4_a));
    store(mosaic.iso4_a) = max(0.0, store(mosaic.iso4_a));
  } // form_nh4so4

  KOKKOS_INLINE_FUNCTION static void
  form_nacl(const MosaicModelData<DeviceType>& mosaic,
            const real_type& jp,
            const real_type_1d_view_type& electrolyte,
            const real_type_1d_view_type& store,
            const real_type_1d_view_type& aer_curr,
            const real_type_1d_view_type& aer_solid,
            const real_type_1d_view_type& aer_liquid,
            const real_type_1d_view_type& aer_total,
            const real_type_1d_view_type& gas,
            const real_type_1d_view_type& total_species,
            real_type& tot_cl_in) {
    
    electrolyte(mosaic.jnacl) = store(mosaic.ina_a);
    
    store(mosaic.ina_a) = 0.0;
    store(mosaic.icl_a) = store(mosaic.icl_a) - electrolyte(mosaic.jnacl);
    
    if (store(mosaic.icl_a) < 0.0) { // cl deficit in aerosol. take some from gas
      aer_curr(mosaic.icl_a) = aer_curr(mosaic.icl_a) - store(mosaic.icl_a);

      if (jp != mosaic.jtotal) {
        aer_total(mosaic.icl_a) = aer_solid(mosaic.icl_a) + aer_liquid(mosaic.icl_a);
      }
      
      gas(mosaic.ihcl_g) = gas(mosaic.ihcl_g) + store(mosaic.icl_a);
      
      if (gas(mosaic.ihcl_g) < 0.0) {
        total_species(mosaic.ihcl_g) = total_species(mosaic.ihcl_g) - gas(mosaic.ihcl_g);
        tot_cl_in = tot_cl_in - gas(mosaic.ihcl_g);
      }
      
      gas(mosaic.ihcl_g) = max(0.0, gas(mosaic.ihcl_g));
      
      store(mosaic.icl_a) = 0.0;
    }
    
    store(mosaic.icl_a) = max(0.0, store(mosaic.icl_a));
  } //form_nacl
  
  KOKKOS_INLINE_FUNCTION static
  void form_namsa(const MosaicModelData<DeviceType>& mosaic,
                  const real_type_1d_view_type& electrolyte,
                  const real_type_1d_view_type& store) {

    electrolyte(mosaic.jnamsa) = min(store(mosaic.ina_a), store(mosaic.imsa_a));

    store(mosaic.ina_a)  = store(mosaic.ina_a) - electrolyte(mosaic.jnamsa);
    store(mosaic.imsa_a) = store(mosaic.imsa_a) - electrolyte(mosaic.jnamsa);

    store(mosaic.ina_a)  = max(0.0, store(mosaic.ina_a));
    store(mosaic.imsa_a) = max(0.0, store(mosaic.imsa_a));
  } // form_namsa
  
  KOKKOS_INLINE_FUNCTION static
  void form_nano3(const MosaicModelData<DeviceType>& mosaic,
                  const real_type_1d_view_type& electrolyte,
                  const real_type_1d_view_type& store) {

    electrolyte(mosaic.jnano3) = min(store(mosaic.ina_a), store(mosaic.ino3_a));

    store(mosaic.ina_a)  = store(mosaic.ina_a)  - electrolyte(mosaic.jnano3);
    store(mosaic.ino3_a) = store(mosaic.ino3_a) - electrolyte(mosaic.jnano3);

    store(mosaic.ina_a)  = max(0.0, store(mosaic.ina_a));
    store(mosaic.ino3_a) = max(0.0, store(mosaic.ino3_a));
  } // form_nano3

  KOKKOS_INLINE_FUNCTION static
  void form_hcl(const MosaicModelData<DeviceType>& mosaic,
                const real_type_1d_view_type& electrolyte,
                const real_type_1d_view_type& store) {

    electrolyte(mosaic.jhcl) = max(0.0, store(mosaic.icl_a));

    store(mosaic.icl_a) = 0.0;
  } // form_hcl

  KOKKOS_INLINE_FUNCTION static
  void form_hno3(const MosaicModelData<DeviceType>& mosaic,
                 const real_type_1d_view_type& electrolyte,
                 const real_type_1d_view_type& store) {

    electrolyte(mosaic.jhno3) = max(0.0, store(mosaic.ino3_a));

    store(mosaic.ino3_a) = 0.0;
  } // form_hno3

  KOKKOS_INLINE_FUNCTION static
  void form_msa(const MosaicModelData<DeviceType>& mosaic,
                const real_type_1d_view_type& electrolyte,
                const real_type_1d_view_type& store) {

    electrolyte(mosaic.jmsa) = max(0.0, store(mosaic.imsa_a));

    store(mosaic.imsa_a) = 0.0;
  } // form_msa

  KOKKOS_INLINE_FUNCTION static
  void form_h2so4(const MosaicModelData<DeviceType>& mosaic,
                  const real_type_1d_view_type& electrolyte,
                  const real_type_1d_view_type& store) {

    electrolyte(mosaic.jh2so4) = max(0.0, store(mosaic.iso4_a));

    store(mosaic.iso4_a) = 0.0;
  } // form_h2so4
  
  KOKKOS_INLINE_FUNCTION static
  void form_na2so4_nahso4(const MosaicModelData<DeviceType>& mosaic,
                          const real_type_1d_view_type& electrolyte,
                          const real_type_1d_view_type& store) {
      
      electrolyte(mosaic.jna2so4) = store(mosaic.ina_a) - store(mosaic.iso4_a);
      electrolyte(mosaic.jnahso4) = 2.*store(mosaic.iso4_a) - store(mosaic.ina_a);
      electrolyte(mosaic.jna2so4) = max(0.0, electrolyte(mosaic.jna2so4));
      electrolyte(mosaic.jnahso4) = max(0.0, electrolyte(mosaic.jnahso4));

      store(mosaic.ina_a)  = 0.0;
      store(mosaic.iso4_a) = 0.0;
  } // form_na2so4_nahso4
  
  KOKKOS_INLINE_FUNCTION static
  void form_lvcite_nh4hso4(const MosaicModelData<DeviceType>& mosaic,
                           const real_type_1d_view_type& electrolyte,
                           const real_type_1d_view_type& store) {

    electrolyte(mosaic.jlvcite) = store(mosaic.inh4_a) - store(mosaic.iso4_a);
    electrolyte(mosaic.jnh4hso4) = 3.*store(mosaic.iso4_a) - 2.*store(mosaic.inh4_a);
    electrolyte(mosaic.jlvcite) = max(0.0, electrolyte(mosaic.jlvcite));
    electrolyte(mosaic.jnh4hso4) = max(0.0, electrolyte(mosaic.jnh4hso4));

    store(mosaic.inh4_a) = 0.0;
    store(mosaic.iso4_a) = 0.0;
  } // form_lvcite_nh4hso4
  
  KOKKOS_INLINE_FUNCTION static
  void form_nh4so4_lvcite(const MosaicModelData<DeviceType>& mosaic,
                          const real_type_1d_view_type& electrolyte,
                          const real_type_1d_view_type& store) {

    electrolyte(mosaic.jnh4so4) = 2.*store(mosaic.inh4_a) - 3.*store(mosaic.iso4_a);
    electrolyte(mosaic.jlvcite) = 2.*store(mosaic.iso4_a) - store(mosaic.inh4_a);
    electrolyte(mosaic.jnh4so4) = max(0.0, electrolyte(mosaic.jnh4so4));
    electrolyte(mosaic.jlvcite) = max(0.0, electrolyte(mosaic.jlvcite));

    store(mosaic.inh4_a) = 0.0;
    store(mosaic.iso4_a) = 0.0;
  } // form_nh4so4_lvcite
  
  KOKKOS_INLINE_FUNCTION static
  void form_nh4no3(const MosaicModelData<DeviceType>& mosaic,
                  const real_type_1d_view_type& electrolyte,
                  const real_type_1d_view_type& store) {

    electrolyte(mosaic.jnh4no3) = min(store(mosaic.inh4_a), store(mosaic.ino3_a));

    store(mosaic.inh4_a)  = store(mosaic.inh4_a)  - electrolyte(mosaic.jnh4no3);
    store(mosaic.ino3_a) = store(mosaic.ino3_a) - electrolyte(mosaic.jnh4no3);

    store(mosaic.inh4_a)  = max(0.0, store(mosaic.inh4_a));
    store(mosaic.ino3_a) = max(0.0, store(mosaic.ino3_a));
  } // form_nh4no3
  
  KOKKOS_INLINE_FUNCTION static
  void form_nh4cl(const MosaicModelData<DeviceType>& mosaic,
             const real_type_1d_view_type& electrolyte,
             const real_type_1d_view_type& store) {

    electrolyte(mosaic.jnh4cl) = min(store(mosaic.inh4_a), store(mosaic.icl_a));

    store(mosaic.inh4_a)  = store(mosaic.inh4_a)  - electrolyte(mosaic.jnh4cl);
    store(mosaic.icl_a) = store(mosaic.icl_a) - electrolyte(mosaic.jnh4cl);

    store(mosaic.inh4_a)  = max(0.0, store(mosaic.inh4_a));
    store(mosaic.icl_a) = max(0.0, store(mosaic.icl_a));
  } // form_nh4cl

#endif
