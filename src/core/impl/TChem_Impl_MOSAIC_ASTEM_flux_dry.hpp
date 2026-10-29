#ifndef __TCHEM_IMPL_MOSAIC_ASTEM_FLUX_DRY_HPP__
#define __TCHEM_IMPL_MOSAIC_ASTEM_FLUX_DRY_HPP__

  //--------------------------------------------------------------------
  // DRY PARTICLES
  //--------------------------------------------------------------------
  // part of ASTEM: computes gas-aerosol fluxes over dry aerosols
  KOKKOS_INLINE_FUNCTION static
  void ASTEM_flux_dry(const MosaicModelData<DeviceType>& mosaic,
                      const real_type_1d_view_type& gas,                // (ngas_volatile)
                      const real_type_1d_view_type& kg,                 // (ngas_volatile)
                      const real_type_1d_view_type& Keq_sg,             // (nrxn_aer_sg)
                      const real_type& Kp_nh4no3_0,
                      const real_type& Kp_nh4cl_0,
                      const real_type_1d_view_type& electrolyte_solid,  // (nelectrolyte)
                      const real_type_1d_view_type& epercent_solid,     // (nelectrolyte)
                      const real_type_1d_view_type& aer_solid,          // (naer)
                      const real_type_1d_view_type& aer_liquid,         // (naer)
                      const real_type_1d_view_type& aer_total,          // (naer)
                      real_type& Keq_nh4no3_0,
                      real_type& Keq_nh4cl_0,
                      real_type& phi_nh4no3_s,
                      real_type& phi_nh4cl_s,
                      const real_type_1d_view_type& sfc_a,              // (ngas_volatile)
                      const real_type_1d_view_type& df_gas_s,           // (ngas_volatile)
                      const real_type_1d_view_type& phi_volatile_s,     // (ngas_volatile)
                      const real_type_1d_view_type& flux_s,             // (ngas_volatile)
                      const real_type_1d_view_type& integrate_solid,    // (ngas_volatile)
                      real_type& idry_case3a,
                      real_type& ieqblm_ASTEM) {

    // local variables
    real_type XT, prod_nh4no3, prod_nh4cl, volatile_cl;

    XT = 0.0;
    calculate_XT(mosaic, aer_solid, XT);

    //-----------------------------------------------------------------
    // CASE 1:  caco3 > 0 absorb all acids (and indirectly degas co2)

    if (electrolyte_solid(mosaic.jcaco3) > 0.0) {
      ASTEM_flux_dry_case1(mosaic, gas, kg,
                           sfc_a, df_gas_s, phi_volatile_s, flux_s, integrate_solid,
                           ieqblm_ASTEM);
      return;
    }

    //-----------------------------------------------------------------
    // CASE 2: Sulfate-Rich Domain

    if (XT < 1.9999 && XT >= 0.0) { // excess sulfate (acidic)
      ASTEM_flux_dry_case2(mosaic, gas, kg,
                           sfc_a, df_gas_s, phi_volatile_s, flux_s, integrate_solid,
                           ieqblm_ASTEM);
      return;
    }

    //-------------------------------------------------------------------
    // CASE 3: hno3 and hcl exchange may happen here and nh4cl may form/evaporate

    volatile_cl = electrolyte_solid(mosaic.jnacl) +
                  electrolyte_solid(mosaic.jcacl2);

    if (volatile_cl > 0.0 && gas(mosaic.ihno3_g) > 0.0) {

      ASTEM_flux_dry_case3a(mosaic, gas, kg, aer_solid,
                            sfc_a, df_gas_s, phi_volatile_s, flux_s, integrate_solid,
                            idry_case3a, ieqblm_ASTEM);

      Keq_nh4cl_0 = min(Kp_nh4cl_0, Keq_sg(1));
      prod_nh4cl = max((gas(mosaic.inh3_g)*gas(mosaic.ihcl_g) - Keq_nh4cl_0), 0.0) +
                   electrolyte_solid(mosaic.jnh4cl);

      if (prod_nh4cl > 0.0) {
        ASTEM_flux_dry_case3b(mosaic, gas, kg, Keq_sg,
                              electrolyte_solid, epercent_solid,
                              aer_solid, aer_liquid, aer_total,
                              Keq_nh4cl_0, phi_nh4cl_s,
                              sfc_a, df_gas_s, phi_volatile_s, flux_s, integrate_solid,
                              ieqblm_ASTEM);
      }

      return;
    }

    //-----------------------------------------------------------------
    // CASE 4: nh4no3 or nh4cl or both may be active

    Keq_nh4no3_0 = min(Kp_nh4no3_0, Keq_sg(0));
    Keq_nh4cl_0  = min(Kp_nh4cl_0,  Keq_sg(1));

    prod_nh4no3 = max((gas(mosaic.inh3_g)*gas(mosaic.ihno3_g) - Keq_nh4no3_0), 0.0) +
                  electrolyte_solid(mosaic.jnh4no3);
    prod_nh4cl  = max((gas(mosaic.inh3_g)*gas(mosaic.ihcl_g) - Keq_nh4cl_0), 0.0) +
                  electrolyte_solid(mosaic.jnh4cl);

    if (prod_nh4no3 > 0.0 || prod_nh4cl > 0.0) {
      ASTEM_flux_dry_case4(mosaic, gas, kg, Keq_sg,
                           electrolyte_solid, epercent_solid,
                           aer_solid, aer_liquid, aer_total,
                           Keq_nh4no3_0, Keq_nh4cl_0, phi_nh4no3_s, phi_nh4cl_s,
                           sfc_a, df_gas_s, phi_volatile_s, flux_s, integrate_solid,
                           ieqblm_ASTEM);
      return;
    }

    //-----------------------------------------------------------------
    return;
  } // ASTEM_flux_dry

  KOKKOS_INLINE_FUNCTION static
  void ASTEM_flux_dry_case1(const MosaicModelData<DeviceType>& mosaic,
                            const real_type_1d_view_type& gas,
                            const real_type_1d_view_type& kg,
                            const real_type_1d_view_type& sfc_a,
                            const real_type_1d_view_type& df_gas_s,
                            const real_type_1d_view_type& phi_volatile_s,
                            const real_type_1d_view_type& flux_s,
                            const real_type_1d_view_type& integrate_solid,
                            real_type& ieqblm_ASTEM) {

    // CASE 1:  caco3 > 0 absorb all acids (and indirectly degas co2)
    if (gas(mosaic.ihno3_g) > 1.e-6) {
      sfc_a(mosaic.ihno3_g) = 0.0;
      df_gas_s(mosaic.ihno3_g) = gas(mosaic.ihno3_g);
      phi_volatile_s(mosaic.ihno3_g) = 1.0;
      flux_s(mosaic.ihno3_g) = kg(mosaic.ihno3_g)*df_gas_s(mosaic.ihno3_g);
      integrate_solid(mosaic.ihno3_g) = 1.0; // mYES
      ieqblm_ASTEM = 0.0; // mNO
    }

    if (gas(mosaic.ihcl_g) > 1.e-6) {
      sfc_a(mosaic.ihcl_g) = 0.0;
      df_gas_s(mosaic.ihcl_g) = gas(mosaic.ihcl_g);
      phi_volatile_s(mosaic.ihcl_g) = 1.0;
      flux_s(mosaic.ihcl_g) = kg(mosaic.ihcl_g)*df_gas_s(mosaic.ihcl_g);
      integrate_solid(mosaic.ihcl_g) = 1.0; // mYES
      ieqblm_ASTEM = 0.0; // mNO
    }

    return;
  } // ASTEM_flux_dry_case1

  KOKKOS_INLINE_FUNCTION static
  void ASTEM_flux_dry_case2(const MosaicModelData<DeviceType>& mosaic,
                            const real_type_1d_view_type& gas,
                            const real_type_1d_view_type& kg,
                            const real_type_1d_view_type& sfc_a,
                            const real_type_1d_view_type& df_gas_s,
                            const real_type_1d_view_type& phi_volatile_s,
                            const real_type_1d_view_type& flux_s,
                            const real_type_1d_view_type& integrate_solid,
                            real_type& ieqblm_ASTEM) {

    // CASE 2: Sulfate-Rich Domain
    if (gas(mosaic.inh3_g) > 1.e-6) {
      sfc_a(mosaic.inh3_g) = 0.0;
      df_gas_s(mosaic.inh3_g) = gas(mosaic.inh3_g);
      phi_volatile_s(mosaic.inh3_g) = 1.0;
      flux_s(mosaic.inh3_g) = kg(mosaic.inh3_g)*gas(mosaic.inh3_g);
      integrate_solid(mosaic.inh3_g) = 1.0; // mYES
      ieqblm_ASTEM = 0.0; // mNO
    }

    return;
  } // ASTEM_flux_dry_case2

  //--------------------------------------------------------------------
  // Fortran: ASTEM.f90, subroutine ASTEM_flux_dry_case3a(ibin)
  // CASE 3a: hno3 and hcl exchange
  //--------------------------------------------------------------------
  KOKKOS_INLINE_FUNCTION static
  void ASTEM_flux_dry_case3a(const MosaicModelData<DeviceType>& mosaic,
                             const real_type_1d_view_type& gas,
                             const real_type_1d_view_type& kg,
                             const real_type_1d_view_type& aer_solid,
                             const real_type_1d_view_type& sfc_a,
                             const real_type_1d_view_type& df_gas_s,
                             const real_type_1d_view_type& phi_volatile_s,
                             const real_type_1d_view_type& flux_s,
                             const real_type_1d_view_type& integrate_solid,
                             real_type& idry_case3a,
                             real_type& ieqblm_ASTEM) {

    // CASE 3a: degas hcl from nacl or cacl2 by flux_s balance with hno3
    if (gas(mosaic.ihno3_g) > 1.e-6) {
      sfc_a(mosaic.ihno3_g) = 0.0;
      sfc_a(mosaic.ihcl_g)  = gas(mosaic.ihcl_g) + aer_solid(mosaic.icl_a);

      df_gas_s(mosaic.ihno3_g) = gas(mosaic.ihno3_g);
      df_gas_s(mosaic.ihcl_g)  = -aer_solid(mosaic.icl_a);

      flux_s(mosaic.ihno3_g) = kg(mosaic.ihno3_g)*gas(mosaic.ihno3_g);
      flux_s(mosaic.ihcl_g)  = -flux_s(mosaic.ihno3_g);

      phi_volatile_s(mosaic.ihno3_g) = 1.0;
      phi_volatile_s(mosaic.ihcl_g)  = df_gas_s(mosaic.ihcl_g)/sfc_a(mosaic.ihcl_g);

      integrate_solid(mosaic.ihno3_g) = 1.0; // mYES
      integrate_solid(mosaic.ihcl_g)  = 1.0; // mYES

      idry_case3a = 1.0; // mYES
      ieqblm_ASTEM = 0.0; // mNO
    }

    return;
  } // ASTEM_flux_dry_case3a

  KOKKOS_INLINE_FUNCTION static
  void ASTEM_flux_dry_case3b(const MosaicModelData<DeviceType>& mosaic,
                             const real_type_1d_view_type& gas,
                             const real_type_1d_view_type& kg,
                             const real_type_1d_view_type& Keq_sg,
                             const real_type_1d_view_type& electrolyte_solid,
                             const real_type_1d_view_type& epercent_solid,
                             const real_type_1d_view_type& aer_solid,
                             const real_type_1d_view_type& aer_liquid,
                             const real_type_1d_view_type& aer_total,
                             const real_type& Keq_nh4cl_0,
                             real_type& phi_nh4cl_s,
                             const real_type_1d_view_type& sfc_a,
                             const real_type_1d_view_type& df_gas_s,
                             const real_type_1d_view_type& phi_volatile_s,
                             const real_type_1d_view_type& flux_s,
                             const real_type_1d_view_type& integrate_solid,
                             real_type& ieqblm_ASTEM) {

    // CASE 3b: nh4cl may form/evaporate here
    ordinal_type iactive_nh4cl;
    real_type a, b, c;

    // set default values for flags
    iactive_nh4cl = 1;

    // compute relative driving force
    phi_nh4cl_s = (gas(mosaic.inh3_g)*gas(mosaic.ihcl_g) - Keq_nh4cl_0)/
                  max(gas(mosaic.inh3_g)*gas(mosaic.ihcl_g), Keq_nh4cl_0);

    // now determine if nh4cl is active or significant
    // nh4cl
    if (Kokkos::abs(phi_nh4cl_s) < mosaic.rtol_eqb_astem) {
      iactive_nh4cl = 0;
    } else if (gas(mosaic.inh3_g)*gas(mosaic.ihcl_g) < Keq_nh4cl_0 &&
               epercent_solid(mosaic.jnh4cl) <= mosaic.ptol_mol_astem) {
      iactive_nh4cl = 0;
      if (epercent_solid(mosaic.jnh4cl) > 0.0) {
        degas_solid_nh4cl(mosaic, gas, Keq_sg, electrolyte_solid,
                          aer_solid, aer_liquid, aer_total);
      }
    }

    // check the outcome
    if (iactive_nh4cl == 0) return;

    // nh4cl is active

    a =   kg(mosaic.inh3_g);
    b = - kg(mosaic.inh3_g)*gas(mosaic.inh3_g)
        + kg(mosaic.ihcl_g)*gas(mosaic.ihcl_g);
    c = -(kg(mosaic.ihcl_g)*Keq_nh4cl_0);

    quadratic(a, b, c, sfc_a(mosaic.inh3_g));
    sfc_a(mosaic.ihcl_g) = Keq_nh4cl_0/sfc_a(mosaic.inh3_g);

    df_gas_s(mosaic.ihcl_g) = gas(mosaic.ihcl_g) - sfc_a(mosaic.ihcl_g);
    df_gas_s(mosaic.inh3_g) = gas(mosaic.inh3_g) - sfc_a(mosaic.inh3_g);

    flux_s(mosaic.inh3_g) = kg(mosaic.inh3_g)*df_gas_s(mosaic.inh3_g);
    flux_s(mosaic.ihcl_g) = flux_s(mosaic.ihcl_g) + flux_s(mosaic.inh3_g);

    phi_volatile_s(mosaic.inh3_g) = phi_nh4cl_s;

    if (flux_s(mosaic.ihcl_g) > 0.0) {
      df_gas_s(mosaic.ihcl_g) = flux_s(mosaic.ihcl_g)/kg(mosaic.ihcl_g); // recompute df_gas
      phi_volatile_s(mosaic.ihcl_g) = phi_nh4cl_s;
    } else {
      sfc_a(mosaic.ihcl_g) = gas(mosaic.ihcl_g) + aer_solid(mosaic.icl_a);
      df_gas_s(mosaic.ihcl_g) = -aer_solid(mosaic.icl_a);
      phi_volatile_s(mosaic.ihcl_g) = df_gas_s(mosaic.ihcl_g)/sfc_a(mosaic.ihcl_g); // not to be used
    }

    integrate_solid(mosaic.inh3_g) = 1.0; // mYES
    integrate_solid(mosaic.ihcl_g) = 1.0; // mYES: integrate HCl with explicit euler

    ieqblm_ASTEM = 0.0; // mNO

    return;
  } // ASTEM_flux_dry_case3b

  KOKKOS_INLINE_FUNCTION static
  void ASTEM_flux_dry_case4(const MosaicModelData<DeviceType>& mosaic,
                            const real_type_1d_view_type& gas,
                            const real_type_1d_view_type& kg,
                            const real_type_1d_view_type& Keq_sg,
                            const real_type_1d_view_type& electrolyte_solid,
                            const real_type_1d_view_type& epercent_solid,
                            const real_type_1d_view_type& aer_solid,
                            const real_type_1d_view_type& aer_liquid,
                            const real_type_1d_view_type& aer_total,
                            const real_type& Keq_nh4no3_0,
                            const real_type& Keq_nh4cl_0,
                            real_type& phi_nh4no3_s,
                            real_type& phi_nh4cl_s,
                            const real_type_1d_view_type& sfc_a,
                            const real_type_1d_view_type& df_gas_s,
                            const real_type_1d_view_type& phi_volatile_s,
                            const real_type_1d_view_type& flux_s,
                            const real_type_1d_view_type& integrate_solid,
                            real_type& ieqblm_ASTEM) {

    // CASE 4: nh4no3 and/or nh4cl may be active
    // local variables
    ordinal_type iactive_nh4no3, iactive_nh4cl, iactive;

    // set default values for flags
    iactive_nh4no3 = 1;
    iactive_nh4cl  = 2;

    // compute diagnostic products and ratios
    phi_nh4no3_s = (gas(mosaic.inh3_g)*gas(mosaic.ihno3_g) - Keq_nh4no3_0)/
                   max(gas(mosaic.inh3_g)*gas(mosaic.ihno3_g), Keq_nh4no3_0);
    phi_nh4cl_s  = (gas(mosaic.inh3_g)*gas(mosaic.ihcl_g) - Keq_nh4cl_0)/
                   max(gas(mosaic.inh3_g)*gas(mosaic.ihcl_g), Keq_nh4cl_0);

    // now determine if nh4no3 and/or nh4cl are active or significant

    // nh4no3
    if (Kokkos::abs(phi_nh4no3_s) < mosaic.rtol_eqb_astem) {
      iactive_nh4no3 = 0;
    } else if (gas(mosaic.inh3_g)*gas(mosaic.ihno3_g) < Keq_nh4no3_0 &&
               epercent_solid(mosaic.jnh4no3) <= mosaic.ptol_mol_astem) {
      iactive_nh4no3 = 0;
      if (epercent_solid(mosaic.jnh4no3) > 0.0) {
        degas_solid_nh4no3(mosaic, gas, Keq_sg, electrolyte_solid,
                           aer_solid, aer_liquid, aer_total);
      }
    }

    // nh4cl
    if (Kokkos::abs(phi_nh4cl_s) < mosaic.rtol_eqb_astem) {
      iactive_nh4cl = 0;
    } else if (gas(mosaic.inh3_g)*gas(mosaic.ihcl_g) < Keq_nh4cl_0 &&
               epercent_solid(mosaic.jnh4cl) <= mosaic.ptol_mol_astem) {
      iactive_nh4cl = 0;
      if (epercent_solid(mosaic.jnh4cl) > 0.0) {
        degas_solid_nh4cl(mosaic, gas, Keq_sg, electrolyte_solid,
                          aer_solid, aer_liquid, aer_total);
      }
    }

    iactive = iactive_nh4no3 + iactive_nh4cl;

    // check the outcome
    if (iactive == 0) return;

    // goto (1,2,3), iactive
    if (iactive == 1) {
      // only nh4no3 solid is active
      ASTEM_flux_dry_case4a(mosaic, gas, kg, Keq_nh4no3_0, phi_nh4no3_s,
                            sfc_a, df_gas_s, phi_volatile_s, flux_s, integrate_solid,
                            ieqblm_ASTEM);
      return;
    } else if (iactive == 2) {
      // only nh4cl solid is active
      ASTEM_flux_dry_case4b(mosaic, gas, kg, Keq_nh4cl_0, phi_nh4cl_s,
                            sfc_a, df_gas_s, phi_volatile_s, flux_s, integrate_solid,
                            ieqblm_ASTEM);
      return;
    } else if (iactive == 3) {
      // both nh4no3 and nh4cl are active
      ASTEM_flux_dry_case4ab(mosaic, gas, kg, Keq_nh4no3_0, Keq_nh4cl_0,
                             phi_nh4no3_s, phi_nh4cl_s,
                             sfc_a, df_gas_s, phi_volatile_s, flux_s, integrate_solid,
                             ieqblm_ASTEM);
      return;
    }
  } // ASTEM_flux_dry_case4

  KOKKOS_INLINE_FUNCTION static
  void ASTEM_flux_dry_case4a(const MosaicModelData<DeviceType>& mosaic,
                             const real_type_1d_view_type& gas,
                             const real_type_1d_view_type& kg,
                             const real_type& Keq_nh4no3_0,
                             const real_type& phi_nh4no3_s,
                             const real_type_1d_view_type& sfc_a,
                             const real_type_1d_view_type& df_gas_s,
                             const real_type_1d_view_type& phi_volatile_s,
                             const real_type_1d_view_type& flux_s,
                             const real_type_1d_view_type& integrate_solid,
                             real_type& ieqblm_ASTEM) {

    // CASE 4a: nh4no3 solid is active
    // local variables
    real_type a, b, c;

    a =   kg(mosaic.inh3_g);
    b = - kg(mosaic.inh3_g)*gas(mosaic.inh3_g)
        + kg(mosaic.ihno3_g)*gas(mosaic.ihno3_g);
    c = -(kg(mosaic.ihno3_g)*Keq_nh4no3_0);

    quadratic(a, b, c, sfc_a(mosaic.inh3_g));
    sfc_a(mosaic.ihno3_g) = Keq_nh4no3_0/sfc_a(mosaic.inh3_g);

    integrate_solid(mosaic.ihno3_g) = 1.0; // mYES
    integrate_solid(mosaic.inh3_g)  = 1.0; // mYES

    df_gas_s(mosaic.ihno3_g) = gas(mosaic.ihno3_g) - sfc_a(mosaic.ihno3_g);
    df_gas_s(mosaic.inh3_g)  = gas(mosaic.inh3_g)  - sfc_a(mosaic.inh3_g);

    phi_volatile_s(mosaic.ihno3_g) = phi_nh4no3_s;
    phi_volatile_s(mosaic.inh3_g)  = phi_nh4no3_s;

    flux_s(mosaic.ihno3_g) = kg(mosaic.ihno3_g)*df_gas_s(mosaic.ihno3_g);
    flux_s(mosaic.inh3_g)  = flux_s(mosaic.ihno3_g);

    ieqblm_ASTEM = 0.0; // mNO

    return;
  } // ASTEM_flux_dry_case4a

  KOKKOS_INLINE_FUNCTION static
  void ASTEM_flux_dry_case4b(const MosaicModelData<DeviceType>& mosaic,
                             const real_type_1d_view_type& gas,
                             const real_type_1d_view_type& kg,
                             const real_type& Keq_nh4cl_0,
                             const real_type& phi_nh4cl_s,
                             const real_type_1d_view_type& sfc_a,
                             const real_type_1d_view_type& df_gas_s,
                             const real_type_1d_view_type& phi_volatile_s,
                             const real_type_1d_view_type& flux_s,
                             const real_type_1d_view_type& integrate_solid,
                             real_type& ieqblm_ASTEM) {

    // CASE 4b: nh4cl solid is active
    // local variables
    real_type a, b, c;

    a =   kg(mosaic.inh3_g);
    b = - kg(mosaic.inh3_g)*gas(mosaic.inh3_g)
        + kg(mosaic.ihcl_g)*gas(mosaic.ihcl_g);
    c = -(kg(mosaic.ihcl_g)*Keq_nh4cl_0);

    quadratic(a, b, c, sfc_a(mosaic.inh3_g));
    sfc_a(mosaic.ihcl_g) = Keq_nh4cl_0/sfc_a(mosaic.inh3_g);

    integrate_solid(mosaic.ihcl_g) = 1.0; // mYES
    integrate_solid(mosaic.inh3_g) = 1.0; // mYES

    df_gas_s(mosaic.ihcl_g) = gas(mosaic.ihcl_g) - sfc_a(mosaic.ihcl_g);
    df_gas_s(mosaic.inh3_g) = gas(mosaic.inh3_g) - sfc_a(mosaic.inh3_g);

    phi_volatile_s(mosaic.ihcl_g) = phi_nh4cl_s;
    phi_volatile_s(mosaic.inh3_g) = phi_nh4cl_s;

    flux_s(mosaic.ihcl_g) = kg(mosaic.ihcl_g)*df_gas_s(mosaic.ihcl_g);
    flux_s(mosaic.inh3_g) = flux_s(mosaic.ihcl_g);

    ieqblm_ASTEM = 0.0; // mNO

    return;
  } // ASTEM_flux_dry_case4b

  KOKKOS_INLINE_FUNCTION static
  void ASTEM_flux_dry_case4ab(const MosaicModelData<DeviceType>& mosaic,
                              const real_type_1d_view_type& gas,
                              const real_type_1d_view_type& kg,
                              const real_type& Keq_nh4no3_0,
                              const real_type& Keq_nh4cl_0,
                              const real_type& phi_nh4no3_s,
                              const real_type& phi_nh4cl_s,
                              const real_type_1d_view_type& sfc_a,
                              const real_type_1d_view_type& df_gas_s,
                              const real_type_1d_view_type& phi_volatile_s,
                              const real_type_1d_view_type& flux_s,
                              const real_type_1d_view_type& integrate_solid,
                              real_type& ieqblm_ASTEM) {

    // CASE 4ab: both nh4no3 and nh4cl solids are active
    // local variables
    real_type flux_nh3_est, flux_nh3_max, ratio_flux;

    ASTEM_flux_dry_case4a(mosaic, gas, kg, Keq_nh4no3_0, phi_nh4no3_s,
                          sfc_a, df_gas_s, phi_volatile_s, flux_s, integrate_solid,
                          ieqblm_ASTEM);
    ASTEM_flux_dry_case4b(mosaic, gas, kg, Keq_nh4cl_0, phi_nh4cl_s,
                          sfc_a, df_gas_s, phi_volatile_s, flux_s, integrate_solid,
                          ieqblm_ASTEM);

    // estimate nh3 flux and adjust hno3 and/or hcl if necessary

    flux_nh3_est = flux_s(mosaic.ihno3_g) + flux_s(mosaic.ihcl_g);
    flux_nh3_max = kg(mosaic.inh3_g)*gas(mosaic.inh3_g);

    if (flux_nh3_est <= flux_nh3_max) {

      flux_s(mosaic.inh3_g) = flux_nh3_est; // all ok - no adjustments needed
      sfc_a(mosaic.inh3_g)  = gas(mosaic.inh3_g) - // recompute sfc_a(ihno3_g)
                              flux_s(mosaic.inh3_g)/kg(mosaic.inh3_g);
      phi_volatile_s(mosaic.inh3_g) = max(Kokkos::abs(phi_nh4no3_s),
                                          Kokkos::abs(phi_nh4cl_s));

    } else { // reduce hno3 and hcl flux_ses as necessary so that nh3 flux_s = flux_s_nh3_max

      ratio_flux             = flux_nh3_max/flux_nh3_est;
      flux_s(mosaic.inh3_g)  = flux_nh3_max;
      flux_s(mosaic.ihno3_g) = flux_s(mosaic.ihno3_g)*ratio_flux;
      flux_s(mosaic.ihcl_g)  = flux_s(mosaic.ihcl_g) *ratio_flux;

      sfc_a(mosaic.inh3_g)  = 0.0;
      sfc_a(mosaic.ihno3_g) = gas(mosaic.ihno3_g) - // recompute sfc_a(ihno3_g)
                              flux_s(mosaic.ihno3_g)/kg(mosaic.ihno3_g);
      sfc_a(mosaic.ihcl_g)  = gas(mosaic.ihcl_g) - // recompute sfc_a(ihcl_g)
                              flux_s(mosaic.ihcl_g)/kg(mosaic.ihcl_g);

      df_gas_s(mosaic.inh3_g)  = gas(mosaic.inh3_g)  - sfc_a(mosaic.inh3_g);
      df_gas_s(mosaic.ihno3_g) = gas(mosaic.ihno3_g) - sfc_a(mosaic.ihno3_g);
      df_gas_s(mosaic.ihcl_g)  = gas(mosaic.ihcl_g)  - sfc_a(mosaic.ihcl_g);

      phi_volatile_s(mosaic.inh3_g) = max(Kokkos::abs(phi_nh4no3_s),
                                          Kokkos::abs(phi_nh4cl_s));
    }

    ieqblm_ASTEM = 0.0; // mNO

    return;
  } // ASTEM_flux_dry_case4ab

#endif
