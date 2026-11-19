#ifndef __TCHEM_IMPL_MOSAIC_ASTEM_FLUX_MIX_HPP__
#define __TCHEM_IMPL_MOSAIC_ASTEM_FLUX_MIX_HPP__

  //--------------------------------------------------------------------
  // MIXED-PHASE PARTICLES
  //--------------------------------------------------------------------
  // part of ASTEM: computes gas-aerosol fluxes over mixed-phase aerosols
  KOKKOS_INLINE_FUNCTION static
  void ASTEM_flux_mix(const MosaicModelData<DeviceType>& mosaic,
                      const real_type_1d_view_type& gas,                 // (ngas_volatile)
                      const real_type_1d_view_type& kg,                  // (ngas_volatile)
                      const real_type_1d_view_type& kel,                 // (ngas_volatile)
                      const real_type_1d_view_type& Keq_ll,              // (nrxn_aer_ll)
                      const real_type_1d_view_type& Keq_gl,              // (nrxn_aer_gl)
                      const real_type_1d_view_type& Keq_sg,              // (nrxn_aer_sg)
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
                      const real_type_1d_view_type& epercent_solid,      // (nelectrolyte)
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
                      real_type& Keq_nh4no3_0,
                      real_type& Keq_nh4cl_0,
                      real_type& phi_nh4no3_s,
                      real_type& phi_nh4cl_s,
                      const real_type_1d_view_type& sfc_a,               // (ngas_volatile)
                      const real_type_1d_view_type& df_gas_s,            // (ngas_volatile)
                      const real_type_1d_view_type& df_gas_l,            // (ngas_volatile)
                      const real_type_1d_view_type& phi_volatile_s,      // (ngas_volatile)
                      const real_type_1d_view_type& phi_volatile_l,      // (ngas_volatile)
                      const real_type_1d_view_type& flux_s,              // (ngas_volatile)
                      const real_type_1d_view_type& Heff,                // (ngas_volatile)
                      const real_type_1d_view_type& integrate_solid,     // (ngas_volatile)
                      const real_type_1d_view_type& integrate_liquid,    // (ngas_volatile)
                      real_type& idry_case3a,
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
    real_type XT, g_nh3_hno3, g_nh3_hcl, a_nh4_no3, a_nh4_cl, prod_nh4cl;
    real_type volatile_cl;

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
      // MOSAIC prints "ASTEM_flux_wet" here as well
      Kokkos::abort("Water is zero in liquid phase. Stopping in ASTEM_flux_wet\n");
    }

    // reset
    Keq_nh4no3_0 = Keq_sg(0); // 6/25/2008
    Keq_nh4cl_0  = Keq_sg(1); // 6/25/2008

    //-----------------------------------------------------------------
    // MIXED CASE 1:  caco3 > 0 absorb all acids (and indirectly degas co2)

    if (epercent_solid(mosaic.jcaco3) > 0.0) {
      jphase = mosaic.jliquid;
      ASTEM_flux_wet_case1(mosaic, gas, kg, Keq_ll, mc,
                           sfc_a, df_gas_s, phi_volatile_s, flux_s, integrate_solid,
                           jphase, ieqblm_ASTEM);
      return;
    }

    //-----------------------------------------------------------------
    // MIXED CASE 2: Sulfate-Rich Domain

    if (XT < 1.9999 && XT >= 0.0) { // excess sulfate (acidic)
      jphase = mosaic.jliquid;
      ASTEM_flux_wet_case2(mosaic, gas, kel, Keq_ll, Keq_gl, ma, mc, gam, gam_ratio, water_a,
                           sfc_a, df_gas_l, phi_volatile_l, Heff, integrate_liquid,
                           ieqblm_ASTEM);
      return;
    }

    //-------------------------------------------------------------------
    // MIXED CASE 3: hno3 and hcl exchange may happen here and nh4cl may form/evaporate

    volatile_cl = electrolyte_solid(mosaic.jnacl) +
                  electrolyte_solid(mosaic.jcacl2);

    if (volatile_cl > 0.0 && gas(mosaic.ihno3_g) > 0.0) {

      ASTEM_flux_dry_case3a(mosaic, gas, kg, aer_solid,
                            sfc_a, df_gas_s, phi_volatile_s, flux_s, integrate_solid,
                            idry_case3a, ieqblm_ASTEM);

      prod_nh4cl = max((gas(mosaic.inh3_g)*gas(mosaic.ihcl_g) - Keq_sg(1)), 0.0) +
                   electrolyte_solid(mosaic.jnh4cl);

      if (prod_nh4cl > 0.0) {
        ASTEM_flux_dry_case3b(mosaic, gas, kg, Keq_sg,
                              electrolyte_solid, epercent_solid,
                              aer_solid, aer_liquid, aer_total,
                              Keq_nh4cl_0, phi_nh4cl_s,
                              sfc_a, df_gas_s, phi_volatile_s, flux_s, integrate_solid,
                              ieqblm_ASTEM);
      }

      jphase = mosaic.jsolid;

      return;
    }

    //-------------------------------------------------------------------
    // MIXED CASE 4: nh4no3 or nh4cl or both may be active

    if (electrolyte_solid(mosaic.jnh4no3) > 0.0 &&
        electrolyte_solid(mosaic.jnh4cl) > 0.0) {
      jphase = mosaic.jsolid;
      ASTEM_flux_dry_case4(mosaic, gas, kg, Keq_sg,
                           electrolyte_solid, epercent_solid,
                           aer_solid, aer_liquid, aer_total,
                           Keq_nh4no3_0, Keq_nh4cl_0, phi_nh4no3_s, phi_nh4cl_s,
                           sfc_a, df_gas_s, phi_volatile_s, flux_s, integrate_solid,
                           ieqblm_ASTEM);

      if (sfc_a(mosaic.ihno3_g) > 0.0 && ma(mosaic.ja_no3) > 0.0) {
        mc(mosaic.jc_h) = Keq_gl(2)*sfc_a(mosaic.ihno3_g)/
          (kel(mosaic.ihno3_g)*ats<real_type>::pow(gam(mosaic.jhno3), 2.)*ma(mosaic.ja_no3));
      } else if (sfc_a(mosaic.ihcl_g) > 0.0 && ma(mosaic.ja_cl) > 0.0) {
        mc(mosaic.jc_h) = Keq_gl(3)*sfc_a(mosaic.ihcl_g)/
          (kel(mosaic.ihcl_g)*ats<real_type>::pow(gam(mosaic.jhcl), 2.)*ma(mosaic.ja_cl));
      } else {
        mc(mosaic.jc_h) = ats<real_type>::sqrt(Keq_ll(2));
      }

      return;

    } else if (electrolyte_solid(mosaic.jnh4no3) > 0.0) {
      // do small adjustments for nh4cl aq
      g_nh3_hcl = gas(mosaic.inh3_g)*gas(mosaic.ihcl_g);
      a_nh4_cl  = aer_liquid(mosaic.inh4_a)*aer_liquid(mosaic.icl_a);

      iadjust = 0; // mNO, initialize
      if (g_nh3_hcl > 0.0 && a_nh4_cl == 0.0) {
        absorb_tiny_nh4cl(mosaic, aer_solid, aer_liquid, aer_total, gas,
                          delta_nh3_max, delta_hcl_max, electrolyte_sum_total);
        iadjust = 1; // mYES
      } else if (g_nh3_hcl == 0.0 && a_nh4_cl > 0.0) {
        degas_tiny_nh4cl(mosaic, aer_solid, aer_liquid, aer_total, gas,
                         electrolyte_liquid);
        iadjust = 1; // mYES
      }

      if (iadjust == 1) {
        ions_to_electrolytes(mosaic, mosaic.jliquid, aer_liquid, aer_solid, aer_liquid, aer_total,
                             store, electrolyte_liquid, electrolyte_solid, electrolyte_liquid,
                             na, nc, xeq_a, xeq_c, na_Ma, nc_Mc,
                             electrolyte_sum_liquid, epercent_liquid, aer_sum_liquid, aer_percent_liquid,
                             XT); // update after adjustments
        compute_activities(mosaic, molalities, xmol, aer_liquid, ma, mc, Keq_ll,
                           electrolyte_solid, electrolyte_liquid, electrolyte_total,
                           log_gam, log_gamZ, gam, activity,
                           jaerosolstate, jphase, jhyst_leg, aH2O_a, water_a, gam_ratio); // update after adjustments
      }

      ASTEM_flux_mix_case4a(mosaic, gas, kg, kel, Keq_ll, Keq_gl, Keq_sg, Kp_nh4cl, Keq_nh4no3_0,
                            electrolyte_solid, epercent_solid, aer_solid, aer_liquid, aer_total,
                            gam, gam_ratio, water_a, ma, mc, activity,
                            Keq_nh4cl, kelvin_nh4cl, phi_nh4no3_s,
                            sfc_a, df_gas_s, df_gas_l, phi_volatile_s, phi_volatile_l, flux_s, Heff,
                            integrate_solid, integrate_liquid,
                            jphase, ieqblm_ASTEM); // nh4no3 solid + nh4cl aq
      jphase = mosaic.jtotal;
      return;

    } else if (electrolyte_solid(mosaic.jnh4cl) > 0.0) {
      // do small adjustments for nh4no3 aq
      g_nh3_hno3 = gas(mosaic.inh3_g)*gas(mosaic.ihno3_g);
      a_nh4_no3  = aer_liquid(mosaic.inh4_a)*aer_liquid(mosaic.ino3_a);

      iadjust = 0; // mNO, initialize
      if (g_nh3_hno3 > 0.0 && a_nh4_no3 == 0.0) {
        absorb_tiny_nh4no3(mosaic, aer_solid, aer_liquid, aer_total, gas,
                           delta_nh3_max, delta_hno3_max, electrolyte_sum_total);
        iadjust = 1; // mYES
      } else if (g_nh3_hno3 == 0.0 && a_nh4_no3 > 0.0) {
        degas_tiny_nh4no3(mosaic, aer_solid, aer_liquid, aer_total, gas,
                          electrolyte_liquid);
        iadjust = 1; // mYES
      }

      if (iadjust == 1) {
        ions_to_electrolytes(mosaic, mosaic.jliquid, aer_liquid, aer_solid, aer_liquid, aer_total,
                             store, electrolyte_liquid, electrolyte_solid, electrolyte_liquid,
                             na, nc, xeq_a, xeq_c, na_Ma, nc_Mc,
                             electrolyte_sum_liquid, epercent_liquid, aer_sum_liquid, aer_percent_liquid,
                             XT); // update after adjustments
        compute_activities(mosaic, molalities, xmol, aer_liquid, ma, mc, Keq_ll,
                           electrolyte_solid, electrolyte_liquid, electrolyte_total,
                           log_gam, log_gamZ, gam, activity,
                           jaerosolstate, jphase, jhyst_leg, aH2O_a, water_a, gam_ratio); // update after adjustments
      }

      kelvin_nh4no3 = kel(mosaic.inh3_g)*kel(mosaic.ihno3_g);
      Keq_nh4no3 = kelvin_nh4no3*activity(mosaic.jnh4no3)*Kp_nh4no3; // = [NH3]s * [HNO3]s

      ASTEM_flux_mix_case4b(mosaic, gas, kg, kel, Keq_ll, Keq_gl, Keq_sg, Kp_nh4no3, Keq_nh4cl_0,
                            electrolyte_solid, epercent_solid, aer_solid, aer_liquid, aer_total,
                            gam, gam_ratio, water_a, ma, mc, activity,
                            Keq_nh4no3, kelvin_nh4no3, phi_nh4cl_s,
                            sfc_a, df_gas_s, df_gas_l, phi_volatile_s, phi_volatile_l, flux_s, Heff,
                            integrate_solid, integrate_liquid,
                            jphase, ieqblm_ASTEM); // nh4cl solid + nh4no3 aq
      jphase = mosaic.jtotal;
      return;
    }

    //-------------------------------------------------------------------
    // no ammonia in the system: go to MIXED CASE 6 (MOSAIC: goto 10)
    const bool no_ammonia = (gas(mosaic.inh3_g) + aer_liquid(mosaic.inh4_a)) < 1.e-25;

    if (!no_ammonia) {
      //-------------------------------------------------------------------
      // MIXED CASE 5: liquid nh4no3 and/or nh4cl maybe active
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

      // all adjustments done...

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
      jphase = mosaic.jliquid;

      return;
    }

    //-------------------------------------------------------------------
    // MIXED CASE 6: ammonia = 0. liquid hno3 and hcl exchange may happen here
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

    // all adjustments done...

    ASTEM_flux_wet_case4(mosaic, gas, kg, kel, Keq_ll, Keq_gl, gam, water_a, ma, mc,
                         sfc_a, df_gas_l, phi_volatile_l, Heff, integrate_liquid,
                         ieqblm_ASTEM);
    jphase = mosaic.jliquid;

    return;
  } // ASTEM_flux_mix

  //------------------------------------------------------------------
  // Mix Case 4a: NH4NO3 solid maybe active. NH4Cl aq maybe active
  KOKKOS_INLINE_FUNCTION static
  void ASTEM_flux_mix_case4a(const MosaicModelData<DeviceType>& mosaic,
                             const real_type_1d_view_type& gas,
                             const real_type_1d_view_type& kg,
                             const real_type_1d_view_type& kel,
                             const real_type_1d_view_type& Keq_ll,
                             const real_type_1d_view_type& Keq_gl,
                             const real_type_1d_view_type& Keq_sg,
                             const real_type& Kp_nh4cl,
                             const real_type& Keq_nh4no3_0,
                             const real_type_1d_view_type& electrolyte_solid,
                             const real_type_1d_view_type& epercent_solid,
                             const real_type_1d_view_type& aer_solid,
                             const real_type_1d_view_type& aer_liquid,
                             const real_type_1d_view_type& aer_total,
                             const real_type_1d_view_type& gam,
                             const real_type& gam_ratio,
                             const real_type& water_a,
                             const real_type_1d_view_type& ma,
                             const real_type_1d_view_type& mc,
                             const real_type_1d_view_type& activity,
                             real_type& Keq_nh4cl,
                             real_type& kelvin_nh4cl,
                             real_type& phi_nh4no3_s,
                             const real_type_1d_view_type& sfc_a,
                             const real_type_1d_view_type& df_gas_s,
                             const real_type_1d_view_type& df_gas_l,
                             const real_type_1d_view_type& phi_volatile_s,
                             const real_type_1d_view_type& phi_volatile_l,
                             const real_type_1d_view_type& flux_s,
                             const real_type_1d_view_type& Heff,
                             const real_type_1d_view_type& integrate_solid,
                             const real_type_1d_view_type& integrate_liquid,
                             real_type& jphase,
                             real_type& ieqblm_ASTEM) {

    ordinal_type iactive_nh4no3, iactive_nh4cl;

    // set default values for flags
    iactive_nh4no3 = 1; // mYES
    iactive_nh4cl  = 1; // mYES

    // nh4no3 (solid)
    phi_nh4no3_s = (gas(mosaic.inh3_g)*gas(mosaic.ihno3_g) - Keq_sg(0))/
                   max(gas(mosaic.inh3_g)*gas(mosaic.ihno3_g), Keq_sg(0));

    // nh4cl (liquid)
    kelvin_nh4cl = kel(mosaic.inh3_g)*kel(mosaic.ihcl_g);
    Keq_nh4cl = kelvin_nh4cl*activity(mosaic.jnh4cl)*Kp_nh4cl; // = [NH3]s * [HCl]s

    //-------------------
    // now determine if nh4no3 and/or nh4cl are active or significant
    // nh4no3 solid
    if (Kokkos::abs(phi_nh4no3_s) <= mosaic.rtol_eqb_astem) {
      iactive_nh4no3 = 0; // mNO
    } else if (gas(mosaic.inh3_g)*gas(mosaic.ihno3_g) < Keq_sg(0) &&
               epercent_solid(mosaic.jnh4no3) <= mosaic.ptol_mol_astem) {
      iactive_nh4no3 = 0; // mNO
      if (epercent_solid(mosaic.jnh4no3) > 0.0) {
        degas_solid_nh4no3(mosaic, gas, Keq_sg, electrolyte_solid,
                           aer_solid, aer_liquid, aer_total);
      }
    }

    // nh4cl aq
    if (gas(mosaic.inh3_g)*gas(mosaic.ihcl_g) == 0.0 || Keq_nh4cl == 0.0) {
      iactive_nh4cl = 0; // mNO
    }

    //---------------------------------
    if (iactive_nh4no3 == 1) {

      jphase = mosaic.jsolid;
      ASTEM_flux_dry_case4a(mosaic, gas, kg, Keq_nh4no3_0, phi_nh4no3_s,
                            sfc_a, df_gas_s, phi_volatile_s, flux_s, integrate_solid,
                            ieqblm_ASTEM); // NH4NO3 (solid)

      if (sfc_a(mosaic.ihno3_g) > 0.0 && ma(mosaic.ja_no3) > 0.0) {
        mc(mosaic.jc_h) = Keq_gl(2)*sfc_a(mosaic.ihno3_g)/
          (kel(mosaic.ihno3_g)*ats<real_type>::pow(gam(mosaic.jhno3), 2.)*ma(mosaic.ja_no3));
      } else if (sfc_a(mosaic.ihcl_g) > 0.0 && ma(mosaic.ja_cl) > 0.0) {
        mc(mosaic.jc_h) = Keq_gl(3)*sfc_a(mosaic.ihcl_g)/
          (kel(mosaic.ihcl_g)*ats<real_type>::pow(gam(mosaic.jhcl), 2.)*ma(mosaic.ja_cl));
      } else {
        mc(mosaic.jc_h) = ats<real_type>::sqrt(Keq_ll(2));
      }

    }

    if (iactive_nh4cl == 1) {

      jphase = mosaic.jliquid;
      ASTEM_flux_wet_case3b(mosaic, gas, kg, kel, Keq_ll, Keq_gl, Keq_nh4cl,
                            gam, gam_ratio, water_a, ma, mc,
                            sfc_a, df_gas_l, phi_volatile_l, Heff, integrate_liquid,
                            ieqblm_ASTEM); // NH4Cl (liquid)

      if (sfc_a(mosaic.ihcl_g) > 0.0 && ma(mosaic.ja_cl) > 0.0) {
        mc(mosaic.jc_h) = Keq_gl(3)*sfc_a(mosaic.ihcl_g)/
          (kel(mosaic.ihcl_g)*ats<real_type>::pow(gam(mosaic.jhcl), 2.)*ma(mosaic.ja_cl));
      } else {
        mc(mosaic.jc_h) = ats<real_type>::sqrt(Keq_ll(2));
      }

    }

    if (iactive_nh4cl == 1 && iactive_nh4no3 == 1) {
      jphase = mosaic.jtotal;
    }

    return;
  } // ASTEM_flux_mix_case4a

  //------------------------------------------------------------------
  // Mix Case 4b: NH4Cl solid maybe active. NH4NO3 aq may or maybe active
  KOKKOS_INLINE_FUNCTION static
  void ASTEM_flux_mix_case4b(const MosaicModelData<DeviceType>& mosaic,
                             const real_type_1d_view_type& gas,
                             const real_type_1d_view_type& kg,
                             const real_type_1d_view_type& kel,
                             const real_type_1d_view_type& Keq_ll,
                             const real_type_1d_view_type& Keq_gl,
                             const real_type_1d_view_type& Keq_sg,
                             const real_type& Kp_nh4no3,
                             const real_type& Keq_nh4cl_0,
                             const real_type_1d_view_type& electrolyte_solid,
                             const real_type_1d_view_type& epercent_solid,
                             const real_type_1d_view_type& aer_solid,
                             const real_type_1d_view_type& aer_liquid,
                             const real_type_1d_view_type& aer_total,
                             const real_type_1d_view_type& gam,
                             const real_type& gam_ratio,
                             const real_type& water_a,
                             const real_type_1d_view_type& ma,
                             const real_type_1d_view_type& mc,
                             const real_type_1d_view_type& activity,
                             real_type& Keq_nh4no3,
                             real_type& kelvin_nh4no3,
                             real_type& phi_nh4cl_s,
                             const real_type_1d_view_type& sfc_a,
                             const real_type_1d_view_type& df_gas_s,
                             const real_type_1d_view_type& df_gas_l,
                             const real_type_1d_view_type& phi_volatile_s,
                             const real_type_1d_view_type& phi_volatile_l,
                             const real_type_1d_view_type& flux_s,
                             const real_type_1d_view_type& Heff,
                             const real_type_1d_view_type& integrate_solid,
                             const real_type_1d_view_type& integrate_liquid,
                             real_type& jphase,
                             real_type& ieqblm_ASTEM) {

    ordinal_type iactive_nh4no3, iactive_nh4cl;

    // set default values for flags
    iactive_nh4cl  = 1; // mYES
    iactive_nh4no3 = 1; // mYES

    // nh4cl (solid)
    phi_nh4cl_s = (gas(mosaic.inh3_g)*gas(mosaic.ihcl_g) - Keq_sg(1))/
                  max(gas(mosaic.inh3_g)*gas(mosaic.ihcl_g), Keq_sg(1));

    // nh4no3 (liquid)
    kelvin_nh4no3 = kel(mosaic.inh3_g)*kel(mosaic.ihno3_g);
    Keq_nh4no3 = kelvin_nh4no3*activity(mosaic.jnh4no3)*Kp_nh4no3; // = [NH3]s * [HNO3]s

    //-------------------
    // now determine if nh4no3 and/or nh4cl are active or significant
    // nh4cl (solid)
    if (Kokkos::abs(phi_nh4cl_s) <= mosaic.rtol_eqb_astem) {
      iactive_nh4cl = 0; // mNO
    } else if (gas(mosaic.inh3_g)*gas(mosaic.ihcl_g) < Keq_sg(1) &&
               epercent_solid(mosaic.jnh4cl) <= mosaic.ptol_mol_astem) {
      iactive_nh4cl = 0; // mNO
      if (epercent_solid(mosaic.jnh4cl) > 0.0) {
        degas_solid_nh4cl(mosaic, gas, Keq_sg, electrolyte_solid,
                          aer_solid, aer_liquid, aer_total);
      }
    }

    // nh4no3 (liquid)
    if (gas(mosaic.inh3_g)*gas(mosaic.ihno3_g) == 0.0 || Keq_nh4no3 == 0.0) {
      iactive_nh4no3 = 0; // mNO
    }

    //---------------------------------
    if (iactive_nh4cl == 1) {

      jphase = mosaic.jsolid;
      ASTEM_flux_dry_case4b(mosaic, gas, kg, Keq_nh4cl_0, phi_nh4cl_s,
                            sfc_a, df_gas_s, phi_volatile_s, flux_s, integrate_solid,
                            ieqblm_ASTEM); // NH4Cl (solid)

      if (sfc_a(mosaic.ihcl_g) > 0.0 && ma(mosaic.ja_cl) > 0.0) {
        mc(mosaic.jc_h) = Keq_gl(3)*sfc_a(mosaic.ihcl_g)/
          (kel(mosaic.ihcl_g)*ats<real_type>::pow(gam(mosaic.jhcl), 2.)*ma(mosaic.ja_cl));
      } else if (sfc_a(mosaic.ihno3_g) > 0.0 && ma(mosaic.ja_no3) > 0.0) {
        mc(mosaic.jc_h) = Keq_gl(2)*sfc_a(mosaic.ihno3_g)/
          (kel(mosaic.ihno3_g)*ats<real_type>::pow(gam(mosaic.jhno3), 2.)*ma(mosaic.ja_no3));
      } else {
        mc(mosaic.jc_h) = ats<real_type>::sqrt(Keq_ll(2));
      }

    }

    if (iactive_nh4no3 == 1) {

      jphase = mosaic.jliquid;
      ASTEM_flux_wet_case3a(mosaic, gas, kg, kel, Keq_ll, Keq_gl, Keq_nh4no3,
                            gam, gam_ratio, water_a, ma, mc,
                            sfc_a, df_gas_l, phi_volatile_l, Heff, integrate_liquid,
                            ieqblm_ASTEM); // NH4NO3 (liquid)

      if (sfc_a(mosaic.ihno3_g) > 0.0 && ma(mosaic.ja_no3) > 0.0) {
        mc(mosaic.jc_h) = Keq_gl(2)*sfc_a(mosaic.ihno3_g)/
          (kel(mosaic.ihno3_g)*ats<real_type>::pow(gam(mosaic.jhno3), 2.)*ma(mosaic.ja_no3));
      } else {
        mc(mosaic.jc_h) = ats<real_type>::sqrt(Keq_ll(2));
      }

    }

    if (iactive_nh4cl == 1 && iactive_nh4no3 == 1) {
      jphase = mosaic.jtotal;
    }

    return;
  } // ASTEM_flux_mix_case4b

#endif
