#ifndef __TCHEM_IMPL_MOSAIC_ASTEM_FLUX_WET_HPP__
#define __TCHEM_IMPL_MOSAIC_ASTEM_FLUX_WET_HPP__

  //--------------------------------------------------------------------
  // WET PARTICLES
  //--------------------------------------------------------------------
  // part of ASTEM: computes gas-aerosol fluxes over wet aerosols
  KOKKOS_INLINE_FUNCTION static
  void ASTEM_flux_wet(const MosaicModelData<DeviceType>& mosaic,
                      const real_type_1d_view_type& gas,                 // (ngas_volatile)
                      const real_type_1d_view_type& kg,                  // (ngas_volatile)
                      const real_type_1d_view_type& kel,                 // (ngas_volatile)
                      const real_type_1d_view_type& Keq_ll,              // (nrxn_aer_ll)
                      const real_type_1d_view_type& Keq_gl,              // (nrxn_aer_gl)
                      const real_type_2d_view_type& log_gamZ,            // (nelectrolyte, nelectrolyte)
                      const real_type& Kp_nh4no3,
                      const real_type& Kp_nh4cl,
                      const real_type_1d_view_type& total_species,       // (ngas_volatile)
                      real_type& tot_cl_in,
                      const real_type& delta_nh3_max,
                      const real_type& delta_hno3_max,
                      const real_type& delta_hcl_max,
                      const real_type_1d_view_type& aer_solid,           // (naer)
                      const real_type_1d_view_type& aer_liquid,          // (naer)
                      const real_type_1d_view_type& aer_total,           // (naer)
                      const real_type_1d_view_type& electrolyte_solid,   // (nelectrolyte)
                      const real_type_1d_view_type& electrolyte_liquid,  // (nelectrolyte)
                      const real_type_1d_view_type& electrolyte_total,   // (nelectrolyte)
                      real_type& electrolyte_sum_liquid,
                      const real_type& electrolyte_sum_total,
                      const real_type_1d_view_type& epercent_liquid,     // (nelectrolyte)
                      real_type& aer_sum_liquid,
                      const real_type_1d_view_type& aer_percent_liquid,  // (naer)
                      const real_type_1d_view_type& ma,                  // (nanion)
                      const real_type_1d_view_type& mc,                  // (ncation)
                      const real_type_1d_view_type& gam,                 // (nelectrolyte)
                      const real_type_1d_view_type& activity,            // (nelectrolyte)
                      real_type& gam_ratio,
                      real_type& jaerosolstate,
                      real_type& jphase,
                      real_type& jhyst_leg,
                      real_type& aH2O_a,
                      real_type& water_a,
                      real_type& Keq_nh4no3,
                      real_type& Keq_nh4cl,
                      real_type& kelvin_nh4no3,
                      real_type& kelvin_nh4cl,
                      const real_type_1d_view_type& sfc_a,               // (ngas_volatile)
                      const real_type_1d_view_type& df_gas_s,            // (ngas_volatile)
                      const real_type_1d_view_type& df_gas_l,            // (ngas_volatile)
                      const real_type_1d_view_type& phi_volatile_s,      // (ngas_volatile)
                      const real_type_1d_view_type& phi_volatile_l,      // (ngas_volatile)
                      const real_type_1d_view_type& flux_s,              // (ngas_volatile)
                      const real_type_1d_view_type& Heff,                // (ngas_volatile)
                      const real_type_1d_view_type& integrate_solid,     // (ngas_volatile)
                      const real_type_1d_view_type& integrate_liquid,    // (ngas_volatile)
                      real_type& ieqblm_ASTEM,
                      const real_type_1d_view_type& store,               // (naer)
                      const real_type_1d_view_type& na,                  // (nanion)
                      const real_type_1d_view_type& nc,                  // (ncation)
                      const real_type_1d_view_type& xeq_a,               // (nanion)
                      const real_type_1d_view_type& xeq_c,               // (ncation)
                      const real_type_1d_view_type& na_Ma,               // (nanion)
                      const real_type_1d_view_type& nc_Mc,               // (ncation)
                      const real_type_1d_view_type& molalities,          // (nelectrolyte)
                      const real_type_1d_view_type& xmol,                // (nelectrolyte)
                      const real_type_1d_view_type& log_gam) {           // (nelectrolyte)

    ordinal_type iadjust, iadjust_intermed;
    real_type XT, g_nh3_hno3, g_nh3_hcl, a_nh4_no3, a_nh4_cl;

    XT = 0.0;
    ions_to_electrolytes(mosaic, mosaic.jliquid, aer_liquid, aer_solid, aer_liquid, aer_total,
                         store, electrolyte_liquid, electrolyte_solid, electrolyte_liquid,
                         na, nc, xeq_a, xeq_c, na_Ma, nc_Mc,
                         electrolyte_sum_liquid, epercent_liquid, aer_sum_liquid, aer_percent_liquid,
                         XT); // for water content calculation
    compute_activities(mosaic, molalities, xmol, aer_liquid, ma, mc, Keq_ll,
                       electrolyte_solid, electrolyte_liquid, electrolyte_total,
                       log_gam, log_gamZ, gam, activity,
                       jaerosolstate, jphase, jhyst_leg, aH2O_a, water_a, gam_ratio);

    if (water_a == 0.0) {
      Kokkos::abort("Water is zero in liquid phase. Stopping in ASTEM_flux_wet\n");
    }

    //-------------------------------------------------------------------
    // CASE 1: caco3 > 0 absorb acids (and indirectly degas co2)

    if (electrolyte_solid(mosaic.jcaco3) > 0.0) {
      ASTEM_flux_wet_case1(mosaic, gas, kg, Keq_ll, mc,
                           sfc_a, df_gas_s, phi_volatile_s, flux_s, integrate_solid,
                           jphase, ieqblm_ASTEM);
      return;
    }

    //-------------------------------------------------------------------
    // CASE 2: Sulfate-Rich Domain

    if (XT < 1.9999 && XT >= 0.0) {
      ASTEM_flux_wet_case2(mosaic, gas, kel, Keq_ll, Keq_gl, ma, mc, gam, gam_ratio, water_a,
                           sfc_a, df_gas_l, phi_volatile_l, Heff, integrate_liquid,
                           ieqblm_ASTEM);
      return;
    }

    //-------------------------------------------------------------------
    // no ammonia in the system: go to CASE 4 (MOSAIC: goto 10)
    const bool no_ammonia = (gas(mosaic.inh3_g) + aer_liquid(mosaic.inh4_a)) < 1.e-25;

    if (!no_ammonia) {
      //-------------------------------------------------------------------
      // CASE 3: nh4no3 and/or nh4cl maybe active
      // do some small adjustments (if needed) before deciding case 3

      iadjust = 0;          // mNO, default
      iadjust_intermed = 0; // mNO, default

      // nh4no3
      g_nh3_hno3 = gas(mosaic.inh3_g)*gas(mosaic.ihno3_g);
      a_nh4_no3  = aer_liquid(mosaic.inh4_a)*aer_liquid(mosaic.ino3_a);

      if (g_nh3_hno3 > 0.0 && a_nh4_no3 == 0.0) {
        absorb_tiny_nh4no3(mosaic, aer_solid, aer_liquid, aer_total, gas,
                           delta_nh3_max, delta_hno3_max, electrolyte_sum_total);
        iadjust = 1;          // mYES
        iadjust_intermed = 1; // mYES
      }

      if (iadjust_intermed == 1) {
        ions_to_electrolytes(mosaic, mosaic.jliquid, aer_liquid, aer_solid, aer_liquid, aer_total,
                             store, electrolyte_liquid, electrolyte_solid, electrolyte_liquid,
                             na, nc, xeq_a, xeq_c, na_Ma, nc_Mc,
                             electrolyte_sum_liquid, epercent_liquid, aer_sum_liquid, aer_percent_liquid,
                             XT); // update after adjustments
        iadjust_intermed = 0; // mNO, reset
      }

      // nh4cl
      g_nh3_hcl = gas(mosaic.inh3_g)*gas(mosaic.ihcl_g);
      a_nh4_cl  = aer_liquid(mosaic.inh4_a)*aer_liquid(mosaic.icl_a);

      if (g_nh3_hcl > 0.0 && a_nh4_cl == 0.0) {
        absorb_tiny_nh4cl(mosaic, aer_solid, aer_liquid, aer_total, gas,
                          delta_nh3_max, delta_hcl_max, electrolyte_sum_total);
        iadjust = 1;          // mYES
        iadjust_intermed = 1; // mYES
      }

      if (iadjust_intermed == 1) {
        ions_to_electrolytes(mosaic, mosaic.jliquid, aer_liquid, aer_solid, aer_liquid, aer_total,
                             store, electrolyte_liquid, electrolyte_solid, electrolyte_liquid,
                             na, nc, xeq_a, xeq_c, na_Ma, nc_Mc,
                             electrolyte_sum_liquid, epercent_liquid, aer_sum_liquid, aer_percent_liquid,
                             XT); // update after adjustments
      }

      if (iadjust == 1) {
        compute_activities(mosaic, molalities, xmol, aer_liquid, ma, mc, Keq_ll,
                           electrolyte_solid, electrolyte_liquid, electrolyte_total,
                           log_gam, log_gamZ, gam, activity,
                           jaerosolstate, jphase, jhyst_leg, aH2O_a, water_a, gam_ratio); // update after adjustments
      }

      kelvin_nh4no3 = kel(mosaic.inh3_g)*kel(mosaic.ihno3_g);
      Keq_nh4no3 = kelvin_nh4no3*activity(mosaic.jnh4no3)*Kp_nh4no3; // = [NH3]s * [HNO3]s

      kelvin_nh4cl = kel(mosaic.inh3_g)*kel(mosaic.ihcl_g);
      Keq_nh4cl = kelvin_nh4cl*activity(mosaic.jnh4cl)*Kp_nh4cl; // = [NH3]s * [HCl]s

      ASTEM_flux_wet_case3(mosaic, gas, kg, kel, Keq_ll, Keq_gl, Keq_nh4no3, Keq_nh4cl,
                           gam, gam_ratio, water_a, ma, mc, activity,
                           aer_solid, aer_liquid, aer_total, total_species, tot_cl_in,
                           electrolyte_liquid, electrolyte_total, electrolyte_sum_liquid,
                           epercent_liquid, aer_sum_liquid, aer_percent_liquid,
                           sfc_a, df_gas_l, phi_volatile_l, Heff, integrate_liquid,
                           ieqblm_ASTEM, store);
      return;
    }

    //-------------------------------------------------------------------
    // CASE 4: ammonia = 0. hno3 and hcl exchange may happen here
    // do small adjustments (if needed) before deciding case 4

    iadjust = 0;          // mNO, default
    iadjust_intermed = 0; // mNO, default

    // hno3
    if (gas(mosaic.ihno3_g) > 0.0 && aer_liquid(mosaic.ino3_a) == 0.0 &&
        aer_liquid(mosaic.icl_a) > 0.0) {
      absorb_tiny_hno3(mosaic, aer_solid, aer_liquid, aer_total, gas,
                       delta_hno3_max); // and degas tiny hcl
      iadjust = 1;          // mYES
      iadjust_intermed = 1; // mYES
    }

    if (iadjust_intermed == 1) {
      ions_to_electrolytes(mosaic, mosaic.jliquid, aer_liquid, aer_solid, aer_liquid, aer_total,
                           store, electrolyte_liquid, electrolyte_solid, electrolyte_liquid,
                           na, nc, xeq_a, xeq_c, na_Ma, nc_Mc,
                           electrolyte_sum_liquid, epercent_liquid, aer_sum_liquid, aer_percent_liquid,
                           XT); // update after adjustments
      iadjust_intermed = 0; // mNO, reset
    }

    // hcl
    if (gas(mosaic.ihcl_g) > 0.0 && aer_liquid(mosaic.icl_a) == 0.0 &&
        aer_liquid(mosaic.ino3_a) > 0.0) {
      absorb_tiny_hcl(mosaic, aer_solid, aer_liquid, aer_total, gas,
                      delta_hcl_max); // and degas tiny hno3
      iadjust = 1;          // mYES
      iadjust_intermed = 1; // mYES
    }

    if (iadjust_intermed == 1) {
      ions_to_electrolytes(mosaic, mosaic.jliquid, aer_liquid, aer_solid, aer_liquid, aer_total,
                           store, electrolyte_liquid, electrolyte_solid, electrolyte_liquid,
                           na, nc, xeq_a, xeq_c, na_Ma, nc_Mc,
                           electrolyte_sum_liquid, epercent_liquid, aer_sum_liquid, aer_percent_liquid,
                           XT); // update after adjustments
    }

    if (iadjust == 1) {
      compute_activities(mosaic, molalities, xmol, aer_liquid, ma, mc, Keq_ll,
                         electrolyte_solid, electrolyte_liquid, electrolyte_total,
                         log_gam, log_gamZ, gam, activity,
                         jaerosolstate, jphase, jhyst_leg, aH2O_a, water_a, gam_ratio); // update after adjustments
    }

    ASTEM_flux_wet_case4(mosaic, gas, kg, kel, Keq_ll, Keq_gl, gam, water_a, ma, mc,
                         sfc_a, df_gas_l, phi_volatile_l, Heff, integrate_liquid,
                         ieqblm_ASTEM);
    return;
  } // ASTEM_flux_wet

  KOKKOS_INLINE_FUNCTION static
  void ASTEM_flux_wet_case1(const MosaicModelData<DeviceType>& mosaic,
                            const real_type_1d_view_type& gas,
                            const real_type_1d_view_type& kg,
                            const real_type_1d_view_type& Keq_ll,
                            const real_type_1d_view_type& mc,
                            const real_type_1d_view_type& sfc_a,
                            const real_type_1d_view_type& df_gas_s,
                            const real_type_1d_view_type& phi_volatile_s,
                            const real_type_1d_view_type& flux_s,
                            const real_type_1d_view_type& integrate_solid,
                            real_type& jphase,
                            real_type& ieqblm_ASTEM) {

    // CASE 1: caco3 > 0 absorb acids (and indirectly degas co2)
    mc(mosaic.jc_h) = ats<real_type>::sqrt(Keq_ll(2));

    // same as dry case1
    if (gas(mosaic.ihno3_g) > 1.e-6) {
      sfc_a(mosaic.ihno3_g) = 0.0;
      df_gas_s(mosaic.ihno3_g) = gas(mosaic.ihno3_g);
      phi_volatile_s(mosaic.ihno3_g) = 1.0;
      flux_s(mosaic.ihno3_g) = kg(mosaic.ihno3_g)*df_gas_s(mosaic.ihno3_g);
      integrate_solid(mosaic.ihno3_g) = 1.0; // mYES
      jphase = mosaic.jsolid;
      ieqblm_ASTEM = 0.0; // mNO
    }

    if (gas(mosaic.ihcl_g) > 1.e-6) {
      sfc_a(mosaic.ihcl_g) = 0.0;
      df_gas_s(mosaic.ihcl_g) = gas(mosaic.ihcl_g);
      phi_volatile_s(mosaic.ihcl_g) = 1.0;
      flux_s(mosaic.ihcl_g) = kg(mosaic.ihcl_g)*df_gas_s(mosaic.ihcl_g);
      integrate_solid(mosaic.ihcl_g) = 1.0; // mYES
      jphase = mosaic.jsolid;
      ieqblm_ASTEM = 0.0; // mNO
    }

    return;
  } // ASTEM_flux_wet_case1

  KOKKOS_INLINE_FUNCTION static
  void ASTEM_flux_wet_case2(const MosaicModelData<DeviceType>& mosaic,
                            const real_type_1d_view_type& gas,
                            const real_type_1d_view_type& kel,
                            const real_type_1d_view_type& Keq_ll,
                            const real_type_1d_view_type& Keq_gl,
                            const real_type_1d_view_type& ma,
                            const real_type_1d_view_type& mc,
                            const real_type_1d_view_type& gam,
                            const real_type& gam_ratio,
                            const real_type& water_a,
                            const real_type_1d_view_type& sfc_a,
                            const real_type_1d_view_type& df_gas_l,
                            const real_type_1d_view_type& phi_volatile_l,
                            const real_type_1d_view_type& Heff,
                            const real_type_1d_view_type& integrate_liquid,
                            real_type& ieqblm_ASTEM) {

    // CASE 2: Sulfate-Rich Domain
    real_type dum_hno3, dum_hcl, dum_nh3;

    sfc_a(mosaic.inh3_g)  = kel(mosaic.inh3_g)*
                            gam_ratio*mc(mosaic.jc_nh4)*Keq_ll(2)/
                            (mc(mosaic.jc_h)*Keq_ll(1)*Keq_gl(1));

    sfc_a(mosaic.ihno3_g) = kel(mosaic.ihno3_g)*
                            mc(mosaic.jc_h)*ma(mosaic.ja_no3)*ats<real_type>::pow(gam(mosaic.jhno3), 2.)/
                            Keq_gl(2);

    sfc_a(mosaic.ihcl_g)  = kel(mosaic.ihcl_g)*
                            mc(mosaic.jc_h)*ma(mosaic.ja_cl)*ats<real_type>::pow(gam(mosaic.jhcl), 2.)/
                            Keq_gl(3);

    dum_hno3 = max(sfc_a(mosaic.ihno3_g), gas(mosaic.ihno3_g));
    dum_hcl  = max(sfc_a(mosaic.ihcl_g),  gas(mosaic.ihcl_g));
    dum_nh3  = max(sfc_a(mosaic.inh3_g),  gas(mosaic.inh3_g));

    // compute relative driving forces
    if (dum_hno3 > 0.0) {
      df_gas_l(mosaic.ihno3_g) = gas(mosaic.ihno3_g) - sfc_a(mosaic.ihno3_g);
      phi_volatile_l(mosaic.ihno3_g) = df_gas_l(mosaic.ihno3_g)/dum_hno3;
    } else {
      phi_volatile_l(mosaic.ihno3_g) = 0.0;
    }

    if (dum_hcl > 0.0) {
      df_gas_l(mosaic.ihcl_g) = gas(mosaic.ihcl_g) - sfc_a(mosaic.ihcl_g);
      phi_volatile_l(mosaic.ihcl_g) = df_gas_l(mosaic.ihcl_g)/dum_hcl;
    } else {
      phi_volatile_l(mosaic.ihcl_g) = 0.0;
    }

    if (dum_nh3 > 0.0) {
      df_gas_l(mosaic.inh3_g) = gas(mosaic.inh3_g) - sfc_a(mosaic.inh3_g);
      phi_volatile_l(mosaic.inh3_g) = df_gas_l(mosaic.inh3_g)/dum_nh3;
    } else {
      phi_volatile_l(mosaic.inh3_g) = 0.0;
    }

    // compute Heff
    if (dum_hno3 > 0.0) {
      Heff(mosaic.ihno3_g) =
        kel(mosaic.ihno3_g)*ats<real_type>::pow(gam(mosaic.jhno3), 2.)*mc(mosaic.jc_h)*1.e-9/
        (water_a*Keq_gl(2));
      integrate_liquid(mosaic.ihno3_g) = 1.0; // mYES
      ieqblm_ASTEM = 0.0; // mNO
    }

    if (dum_hcl > 0.0) {
      Heff(mosaic.ihcl_g) =
        kel(mosaic.ihcl_g)*ats<real_type>::pow(gam(mosaic.jhcl), 2.)*mc(mosaic.jc_h)*1.e-9/
        (water_a*Keq_gl(3));
      integrate_liquid(mosaic.ihcl_g) = 1.0; // mYES
      ieqblm_ASTEM = 0.0; // mNO
    }

    if (dum_nh3 > 0.0) {
      Heff(mosaic.inh3_g) =
        kel(mosaic.inh3_g)*gam_ratio*1.e-9*Keq_ll(2)/
        (water_a*mc(mosaic.jc_h)*Keq_ll(1)*Keq_gl(1));
      integrate_liquid(mosaic.inh3_g) = 1.0; // mYES
      ieqblm_ASTEM = 0.0; // mNO
    }

    return;
  } // ASTEM_flux_wet_case2

  KOKKOS_INLINE_FUNCTION static
  void ASTEM_flux_wet_case3(const MosaicModelData<DeviceType>& mosaic,
                            const real_type_1d_view_type& gas,
                            const real_type_1d_view_type& kg,
                            const real_type_1d_view_type& kel,
                            const real_type_1d_view_type& Keq_ll,
                            const real_type_1d_view_type& Keq_gl,
                            const real_type& Keq_nh4no3,
                            const real_type& Keq_nh4cl,
                            const real_type_1d_view_type& gam,
                            const real_type& gam_ratio,
                            const real_type& water_a,
                            const real_type_1d_view_type& ma,
                            const real_type_1d_view_type& mc,
                            const real_type_1d_view_type& activity,
                            const real_type_1d_view_type& aer_solid,
                            const real_type_1d_view_type& aer_liquid,
                            const real_type_1d_view_type& aer_total,
                            const real_type_1d_view_type& total_species,
                            real_type& tot_cl_in,
                            const real_type_1d_view_type& electrolyte_liquid,
                            const real_type_1d_view_type& electrolyte_total,
                            real_type& electrolyte_sum_liquid,
                            const real_type_1d_view_type& epercent_liquid,
                            real_type& aer_sum_liquid,
                            const real_type_1d_view_type& aer_percent_liquid,
                            const real_type_1d_view_type& sfc_a,
                            const real_type_1d_view_type& df_gas_l,
                            const real_type_1d_view_type& phi_volatile_l,
                            const real_type_1d_view_type& Heff,
                            const real_type_1d_view_type& integrate_liquid,
                            real_type& ieqblm_ASTEM,
                            const real_type_1d_view_type& store) {

    // CASE 3: nh4no3 and/or nh4cl maybe active
    real_type a, b, c, dum_hno3, dum_hcl, dum_nh3;

    a =   kg(mosaic.inh3_g);
    b = - kg(mosaic.inh3_g)*gas(mosaic.inh3_g)
        + kg(mosaic.ihno3_g)*gas(mosaic.ihno3_g)
        + kg(mosaic.ihcl_g)*gas(mosaic.ihcl_g);
    c = -(kg(mosaic.ihno3_g)*Keq_nh4no3 + kg(mosaic.ihcl_g)*Keq_nh4cl);

    quadratic(a, b, c, sfc_a(mosaic.inh3_g));
    sfc_a(mosaic.ihno3_g) = Keq_nh4no3/max(sfc_a(mosaic.inh3_g), 1.e-20);
    sfc_a(mosaic.ihcl_g)  = Keq_nh4cl /max(sfc_a(mosaic.inh3_g), 1.e-20);

    // diagnose mH+
    if (sfc_a(mosaic.ihno3_g) > 0.0 && ma(mosaic.ja_no3) > 0.0) {
      mc(mosaic.jc_h) = Keq_gl(2)*sfc_a(mosaic.ihno3_g)/
        (kel(mosaic.ihno3_g)*ats<real_type>::pow(gam(mosaic.jhno3), 2.)*ma(mosaic.ja_no3));
    } else if (sfc_a(mosaic.ihcl_g) > 0.0 && ma(mosaic.ja_cl) > 0.0) {
      mc(mosaic.jc_h) = Keq_gl(3)*sfc_a(mosaic.ihcl_g)/
        (kel(mosaic.ihcl_g)*ats<real_type>::pow(gam(mosaic.jhcl), 2.)*ma(mosaic.ja_cl));
    } else {
      equilibrate_acids(mosaic, Keq_ll, Keq_gl, water_a, gam,
                        aer_solid, aer_liquid, aer_total, gas, total_species, tot_cl_in,
                        electrolyte_liquid, electrolyte_total, electrolyte_sum_liquid,
                        epercent_liquid, aer_sum_liquid, aer_percent_liquid,
                        ma, mc, activity, store); // hno3 and/or hcl may be > 0 in the gas phase
      mc(mosaic.jc_h) = max(mc(mosaic.jc_h), ats<real_type>::sqrt(Keq_ll(2)));

      sfc_a(mosaic.inh3_g)  = kel(mosaic.inh3_g)*
                              gam_ratio*mc(mosaic.jc_nh4)*Keq_ll(2)/
                              (mc(mosaic.jc_h)*Keq_ll(1)*Keq_gl(1));

      sfc_a(mosaic.ihno3_g) = kel(mosaic.ihno3_g)*
                              mc(mosaic.jc_h)*ma(mosaic.ja_no3)*ats<real_type>::pow(gam(mosaic.jhno3), 2.)/
                              Keq_gl(2);

      sfc_a(mosaic.ihcl_g)  = kel(mosaic.ihcl_g)*
                              mc(mosaic.jc_h)*ma(mosaic.ja_cl)*ats<real_type>::pow(gam(mosaic.jhcl), 2.)/
                              Keq_gl(3);
    }

    // compute Heff
    dum_hno3 = max(sfc_a(mosaic.ihno3_g), gas(mosaic.ihno3_g));
    dum_hcl  = max(sfc_a(mosaic.ihcl_g),  gas(mosaic.ihcl_g));
    dum_nh3  = max(sfc_a(mosaic.inh3_g),  gas(mosaic.inh3_g));

    // compute relative driving forces
    if (dum_hno3 > 0.0) {
      df_gas_l(mosaic.ihno3_g) = gas(mosaic.ihno3_g) - sfc_a(mosaic.ihno3_g);
      phi_volatile_l(mosaic.ihno3_g) = df_gas_l(mosaic.ihno3_g)/dum_hno3;
    } else {
      phi_volatile_l(mosaic.ihno3_g) = 0.0;
    }

    if (dum_hcl > 0.0) {
      df_gas_l(mosaic.ihcl_g) = gas(mosaic.ihcl_g) - sfc_a(mosaic.ihcl_g);
      phi_volatile_l(mosaic.ihcl_g) = df_gas_l(mosaic.ihcl_g)/dum_hcl;
    } else {
      phi_volatile_l(mosaic.ihcl_g) = 0.0;
    }

    if (dum_nh3 > 0.0) {
      df_gas_l(mosaic.inh3_g) = gas(mosaic.inh3_g) - sfc_a(mosaic.inh3_g);
      phi_volatile_l(mosaic.inh3_g) = df_gas_l(mosaic.inh3_g)/dum_nh3;
    } else {
      phi_volatile_l(mosaic.inh3_g) = 0.0;
    }

    if (dum_hno3 > 0.0) {
      Heff(mosaic.ihno3_g) =
        kel(mosaic.ihno3_g)*ats<real_type>::pow(gam(mosaic.jhno3), 2.)*mc(mosaic.jc_h)*1.e-9/
        (water_a*Keq_gl(2));
      integrate_liquid(mosaic.ihno3_g) = 1.0; // mYES
      ieqblm_ASTEM = 0.0; // mNO
    }

    if (dum_hcl > 0.0) {
      Heff(mosaic.ihcl_g) =
        kel(mosaic.ihcl_g)*ats<real_type>::pow(gam(mosaic.jhcl), 2.)*mc(mosaic.jc_h)*1.e-9/
        (water_a*Keq_gl(3));
      integrate_liquid(mosaic.ihcl_g) = 1.0; // mYES
      ieqblm_ASTEM = 0.0; // mNO
    }

    if (dum_nh3 > 0.0) {
      Heff(mosaic.inh3_g) =
        kel(mosaic.inh3_g)*gam_ratio*1.e-9*Keq_ll(2)/
        (water_a*mc(mosaic.jc_h)*Keq_ll(1)*Keq_gl(1));
      integrate_liquid(mosaic.inh3_g) = 1.0; // mYES
      ieqblm_ASTEM = 0.0; // mNO
    }

    return;
  } // ASTEM_flux_wet_case3

  KOKKOS_INLINE_FUNCTION static
  void ASTEM_flux_wet_case3a(const MosaicModelData<DeviceType>& mosaic,
                             const real_type_1d_view_type& gas,
                             const real_type_1d_view_type& kg,
                             const real_type_1d_view_type& kel,
                             const real_type_1d_view_type& Keq_ll,
                             const real_type_1d_view_type& Keq_gl,
                             const real_type& Keq_nh4no3,
                             const real_type_1d_view_type& gam,
                             const real_type& gam_ratio,
                             const real_type& water_a,
                             const real_type_1d_view_type& ma,
                             const real_type_1d_view_type& mc,
                             const real_type_1d_view_type& sfc_a,
                             const real_type_1d_view_type& df_gas_l,
                             const real_type_1d_view_type& phi_volatile_l,
                             const real_type_1d_view_type& Heff,
                             const real_type_1d_view_type& integrate_liquid,
                             real_type& ieqblm_ASTEM) {

    // CASE 3a: nh4no3 (aq) is active
    real_type a, b, c, dum_hno3, dum_nh3;

    a =   kg(mosaic.inh3_g);
    b = - kg(mosaic.inh3_g)*gas(mosaic.inh3_g)
        + kg(mosaic.ihno3_g)*gas(mosaic.ihno3_g);
    c = -(kg(mosaic.ihno3_g)*Keq_nh4no3);

    quadratic(a, b, c, sfc_a(mosaic.inh3_g));
    sfc_a(mosaic.ihno3_g) = Keq_nh4no3/sfc_a(mosaic.inh3_g);

    // diagnose mH+
    if (sfc_a(mosaic.ihno3_g) > 0.0 && ma(mosaic.ja_no3) > 0.0) {
      mc(mosaic.jc_h) = Keq_gl(2)*sfc_a(mosaic.ihno3_g)/
        (kel(mosaic.ihno3_g)*ats<real_type>::pow(gam(mosaic.jhno3), 2.)*ma(mosaic.ja_no3));
    } else {
      mc(mosaic.jc_h) = ats<real_type>::sqrt(Keq_ll(2));
    }

    // compute Heff
    dum_hno3 = max(sfc_a(mosaic.ihno3_g), gas(mosaic.ihno3_g));
    dum_nh3  = max(sfc_a(mosaic.inh3_g),  gas(mosaic.inh3_g));

    // compute relative driving forces
    if (dum_hno3 > 0.0) {
      df_gas_l(mosaic.ihno3_g) = gas(mosaic.ihno3_g) - sfc_a(mosaic.ihno3_g);
      phi_volatile_l(mosaic.ihno3_g) = df_gas_l(mosaic.ihno3_g)/dum_hno3;
    } else {
      phi_volatile_l(mosaic.ihno3_g) = 0.0;
    }

    if (dum_nh3 > 0.0) {
      df_gas_l(mosaic.inh3_g) = gas(mosaic.inh3_g) - sfc_a(mosaic.inh3_g);
      phi_volatile_l(mosaic.inh3_g) = df_gas_l(mosaic.inh3_g)/dum_nh3;
    } else {
      phi_volatile_l(mosaic.inh3_g) = 0.0;
    }

    // compute Heff
    Heff(mosaic.ihno3_g) =
      kel(mosaic.ihno3_g)*ats<real_type>::pow(gam(mosaic.jhno3), 2.)*mc(mosaic.jc_h)*1.e-9/
      (water_a*Keq_gl(2));
    integrate_liquid(mosaic.ihno3_g) = 1.0; // mYES

    Heff(mosaic.inh3_g) =
      kel(mosaic.inh3_g)*gam_ratio*1.e-9*Keq_ll(2)/
      (water_a*mc(mosaic.jc_h)*Keq_ll(1)*Keq_gl(1));
    integrate_liquid(mosaic.inh3_g) = 1.0; // mYES

    ieqblm_ASTEM = 0.0; // mNO

    return;
  } // ASTEM_flux_wet_case3a

  KOKKOS_INLINE_FUNCTION static
  void ASTEM_flux_wet_case3b(const MosaicModelData<DeviceType>& mosaic,
                             const real_type_1d_view_type& gas,
                             const real_type_1d_view_type& kg,
                             const real_type_1d_view_type& kel,
                             const real_type_1d_view_type& Keq_ll,
                             const real_type_1d_view_type& Keq_gl,
                             const real_type& Keq_nh4cl,
                             const real_type_1d_view_type& gam,
                             const real_type& gam_ratio,
                             const real_type& water_a,
                             const real_type_1d_view_type& ma,
                             const real_type_1d_view_type& mc,
                             const real_type_1d_view_type& sfc_a,
                             const real_type_1d_view_type& df_gas_l,
                             const real_type_1d_view_type& phi_volatile_l,
                             const real_type_1d_view_type& Heff,
                             const real_type_1d_view_type& integrate_liquid,
                             real_type& ieqblm_ASTEM) {

    // CASE 3b: nh4cl (aq) is active
    real_type a, b, c, dum_hcl, dum_nh3;

    a =   kg(mosaic.inh3_g);
    b = - kg(mosaic.inh3_g)*gas(mosaic.inh3_g)
        + kg(mosaic.ihcl_g)*gas(mosaic.ihcl_g);
    c = -(kg(mosaic.ihcl_g)*Keq_nh4cl);

    quadratic(a, b, c, sfc_a(mosaic.inh3_g));
    sfc_a(mosaic.ihcl_g) = Keq_nh4cl/sfc_a(mosaic.inh3_g);

    // diagnose mH+
    if (sfc_a(mosaic.ihcl_g) > 0.0 && ma(mosaic.ja_cl) > 0.0) {
      mc(mosaic.jc_h) = Keq_gl(3)*sfc_a(mosaic.ihcl_g)/
        (kel(mosaic.ihcl_g)*ats<real_type>::pow(gam(mosaic.jhcl), 2.)*ma(mosaic.ja_cl));
    } else {
      mc(mosaic.jc_h) = ats<real_type>::sqrt(Keq_ll(2));
    }

    // compute Heff
    dum_hcl = max(sfc_a(mosaic.ihcl_g), gas(mosaic.ihcl_g));
    dum_nh3 = max(sfc_a(mosaic.inh3_g), gas(mosaic.inh3_g));

    // compute relative driving forces
    if (dum_hcl > 0.0) {
      df_gas_l(mosaic.ihcl_g) = gas(mosaic.ihcl_g) - sfc_a(mosaic.ihcl_g);
      phi_volatile_l(mosaic.ihcl_g) = df_gas_l(mosaic.ihcl_g)/dum_hcl;
    } else {
      phi_volatile_l(mosaic.ihcl_g) = 0.0;
    }

    if (dum_nh3 > 0.0) {
      df_gas_l(mosaic.inh3_g) = gas(mosaic.inh3_g) - sfc_a(mosaic.inh3_g);
      phi_volatile_l(mosaic.inh3_g) = df_gas_l(mosaic.inh3_g)/dum_nh3;
    } else {
      phi_volatile_l(mosaic.inh3_g) = 0.0;
    }

    // compute Heff
    Heff(mosaic.ihcl_g) =
      kel(mosaic.ihcl_g)*ats<real_type>::pow(gam(mosaic.jhcl), 2.)*mc(mosaic.jc_h)*1.e-9/
      (water_a*Keq_gl(3));
    integrate_liquid(mosaic.ihcl_g) = 1.0; // mYES

    Heff(mosaic.inh3_g) =
      kel(mosaic.inh3_g)*gam_ratio*1.e-9*Keq_ll(2)/
      (water_a*mc(mosaic.jc_h)*Keq_ll(1)*Keq_gl(1));
    integrate_liquid(mosaic.inh3_g) = 1.0; // mYES

    ieqblm_ASTEM = 0.0; // mNO

    return;
  } // ASTEM_flux_wet_case3b

  KOKKOS_INLINE_FUNCTION static
  void ASTEM_flux_wet_case4(const MosaicModelData<DeviceType>& mosaic,
                            const real_type_1d_view_type& gas,
                            const real_type_1d_view_type& kg,
                            const real_type_1d_view_type& kel,
                            const real_type_1d_view_type& Keq_ll,
                            const real_type_1d_view_type& Keq_gl,
                            const real_type_1d_view_type& gam,
                            const real_type& water_a,
                            const real_type_1d_view_type& ma,
                            const real_type_1d_view_type& mc,
                            const real_type_1d_view_type& sfc_a,
                            const real_type_1d_view_type& df_gas_l,
                            const real_type_1d_view_type& phi_volatile_l,
                            const real_type_1d_view_type& Heff,
                            const real_type_1d_view_type& integrate_liquid,
                            real_type& ieqblm_ASTEM) {

    // CASE 4: ammonia = 0. hno3 and hcl exchange may happen here
    real_type dum_numer, dum_denom, gas_eqb_ratio, dum_hno3, dum_hcl;

    dum_numer = kel(mosaic.ihno3_g)*Keq_gl(3)*ma(mosaic.ja_no3)*
                ats<real_type>::pow(gam(mosaic.jhno3), 2.);
    dum_denom = kel(mosaic.ihcl_g)*Keq_gl(2)*ma(mosaic.ja_cl)*
                ats<real_type>::pow(gam(mosaic.jhcl), 2.);

    if (dum_denom == 0.0 || dum_numer == 0.0) {
      mc(mosaic.jc_h) = ats<real_type>::sqrt(Keq_ll(2));
      return;
    }

    gas_eqb_ratio = dum_numer/dum_denom; // Ce,hno3/Ce,hcl

    // compute equilibrium surface concentrations
    sfc_a(mosaic.ihcl_g) =
      (kg(mosaic.ihno3_g)*gas(mosaic.ihno3_g) + kg(mosaic.ihcl_g)*gas(mosaic.ihcl_g))/
      (kg(mosaic.ihcl_g) + gas_eqb_ratio*kg(mosaic.ihno3_g));
    sfc_a(mosaic.ihno3_g) = gas_eqb_ratio*sfc_a(mosaic.ihcl_g);

    // diagnose mH+
    if (sfc_a(mosaic.ihno3_g) > 0.0 && ma(mosaic.ja_no3) > 0.0) {
      mc(mosaic.jc_h) = Keq_gl(2)*sfc_a(mosaic.ihno3_g)/
        (kel(mosaic.ihno3_g)*ats<real_type>::pow(gam(mosaic.jhno3), 2.)*ma(mosaic.ja_no3));
    } else if (sfc_a(mosaic.ihcl_g) > 0.0 && ma(mosaic.ja_cl) > 0.0) {
      mc(mosaic.jc_h) = Keq_gl(3)*sfc_a(mosaic.ihcl_g)/
        (kel(mosaic.ihcl_g)*ats<real_type>::pow(gam(mosaic.jhcl), 2.)*ma(mosaic.ja_cl));
    } else {
      mc(mosaic.jc_h) = ats<real_type>::sqrt(Keq_ll(2));
    }

    // compute Heff
    dum_hno3 = max(sfc_a(mosaic.ihno3_g), gas(mosaic.ihno3_g));
    dum_hcl  = max(sfc_a(mosaic.ihcl_g),  gas(mosaic.ihcl_g));

    // compute relative driving forces
    if (dum_hno3 > 0.0) {
      df_gas_l(mosaic.ihno3_g) = gas(mosaic.ihno3_g) - sfc_a(mosaic.ihno3_g);
      phi_volatile_l(mosaic.ihno3_g) = df_gas_l(mosaic.ihno3_g)/dum_hno3;
    } else {
      phi_volatile_l(mosaic.ihno3_g) = 0.0;
    }

    if (dum_hcl > 0.0) {
      df_gas_l(mosaic.ihcl_g) = gas(mosaic.ihcl_g) - sfc_a(mosaic.ihcl_g);
      phi_volatile_l(mosaic.ihcl_g) = df_gas_l(mosaic.ihcl_g)/dum_hcl;
    } else {
      phi_volatile_l(mosaic.ihcl_g) = 0.0;
    }

    // compute Heff
    Heff(mosaic.ihno3_g) =
      kel(mosaic.ihno3_g)*ats<real_type>::pow(gam(mosaic.jhno3), 2.)*mc(mosaic.jc_h)*1.e-9/
      (water_a*Keq_gl(2));
    integrate_liquid(mosaic.ihno3_g) = 1.0; // mYES

    Heff(mosaic.ihcl_g) =
      kel(mosaic.ihcl_g)*ats<real_type>::pow(gam(mosaic.jhcl), 2.)*mc(mosaic.jc_h)*1.e-9/
      (water_a*Keq_gl(3));
    integrate_liquid(mosaic.ihcl_g) = 1.0; // mYES

    ieqblm_ASTEM = 0.0; // mNO

    return;
  } // ASTEM_flux_wet_case4

#endif
