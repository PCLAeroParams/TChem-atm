#ifndef __TCHEM_IMPL_MOSAIC_MESA_HPP__
#define __TCHEM_IMPL_MOSAIC_MESA_HPP__

  KOKKOS_INLINE_FUNCTION static
  void adjust_liquid_aerosol(const MosaicModelData<DeviceType>& mosaic,
                             const real_type_1d_view_type& aer_solid,
                             const real_type_1d_view_type& aer_liquid,
                             const real_type_1d_view_type& aer_total,
                             const real_type_1d_view_type& electrolyte_solid,
                             const real_type_1d_view_type& electrolyte_liquid,
                             const real_type_1d_view_type& electrolyte_total,
                             const real_type_1d_view_type& epercent_solid,
                             const real_type_1d_view_type& epercent_liquid,
                             const real_type_1d_view_type& epercent_total,
                             real_type& jphase, real_type& jhyst_leg) {

    jphase    = mosaic.jliquid;
    jhyst_leg = mosaic.jhyst_up; // upper curve

    // Partition all electrolytes into liquid phase
    for (ordinal_type je = 0; je < mosaic.nelectrolyte; ++je) {
      electrolyte_solid(je) = 0.0;
      epercent_solid(je) = 0.0;
      electrolyte_liquid(je) = electrolyte_total(je);
      epercent_liquid(je) = epercent_total(je);
    }

    // Except these electrolytes, which always remain in the solid phase
    electrolyte_solid(mosaic.jcaco3) = electrolyte_total(mosaic.jcaco3);
    electrolyte_solid(mosaic.jcaso4) = electrolyte_total(mosaic.jcaso4);
    epercent_solid(mosaic.jcaco3) = epercent_total(mosaic.jcaco3);
    epercent_solid(mosaic.jcaso4) = epercent_total(mosaic.jcaso4);
    electrolyte_liquid(mosaic.jcaco3) = 0.0;
    electrolyte_liquid(mosaic.jcaso4) = 0.0;
    epercent_liquid(mosaic.jcaco3) = 0.0;
    epercent_liquid(mosaic.jcaso4) = 0.0;

    // Partition all the aer species into solid and liquid phases
    // Solid phase
    aer_solid(mosaic.iso4_a) = electrolyte_solid(mosaic.jcaso4);
    aer_solid(mosaic.ino3_a) = 0.0;
    aer_solid(mosaic.icl_a) = 0.0;
    aer_solid(mosaic.inh4_a) = 0.0;
    aer_solid(mosaic.ioc_a) = aer_total(mosaic.ioc_a);
    aer_solid(mosaic.imsa_a) = 0.0;
    aer_solid(mosaic.ico3_a) = aer_total(mosaic.ico3_a);
    aer_solid(mosaic.ina_a) = 0.0;
    aer_solid(mosaic.ica_a) = electrolyte_solid(mosaic.jcaco3) + electrolyte_solid(mosaic.jcaso4);
    aer_solid(mosaic.ibc_a) = aer_total(mosaic.ibc_a);
    aer_solid(mosaic.ioin_a) = aer_total(mosaic.ioin_a);
    aer_solid(mosaic.iaro1_a) = aer_total(mosaic.iaro1_a);
    aer_solid(mosaic.iaro2_a) = aer_total(mosaic.iaro2_a);
    aer_solid(mosaic.ialk1_a) = aer_total(mosaic.ialk1_a);
    aer_solid(mosaic.iole1_a) = aer_total(mosaic.iole1_a);
    aer_solid(mosaic.iapi1_a) = aer_total(mosaic.iapi1_a);
    aer_solid(mosaic.iapi2_a) = aer_total(mosaic.iapi2_a);
    aer_solid(mosaic.ilim1_a) = aer_total(mosaic.ilim1_a);
    aer_solid(mosaic.ilim2_a) = aer_total(mosaic.ilim2_a);

    // Liquid phase
    aer_liquid(mosaic.iso4_a) = max(0.0, aer_total(mosaic.iso4_a) - aer_solid(mosaic.iso4_a));
    aer_liquid(mosaic.ino3_a) = aer_total(mosaic.ino3_a);
    aer_liquid(mosaic.icl_a) = aer_total(mosaic.icl_a);
    aer_liquid(mosaic.inh4_a) = aer_total(mosaic.inh4_a);
    aer_liquid(mosaic.ioc_a) = 0.0;
    aer_liquid(mosaic.imsa_a) = aer_total(mosaic.imsa_a);
    aer_liquid(mosaic.ico3_a) = 0.0;
    aer_liquid(mosaic.ina_a) = aer_total(mosaic.ina_a);
    aer_liquid(mosaic.ica_a) = max(0.0, aer_total(mosaic.ica_a) - aer_solid(mosaic.ica_a));
    aer_liquid(mosaic.ibc_a) = 0.0;
    aer_liquid(mosaic.ioin_a) = 0.0;
    aer_liquid(mosaic.iaro1_a) = 0.0;
    aer_liquid(mosaic.iaro2_a) = 0.0;
    aer_liquid(mosaic.ialk1_a) = 0.0;
    aer_liquid(mosaic.iole1_a) = 0.0;
    aer_liquid(mosaic.iapi1_a) = 0.0;
    aer_liquid(mosaic.iapi2_a) = 0.0;
    aer_liquid(mosaic.ilim1_a) = 0.0;
    aer_liquid(mosaic.ilim2_a) = 0.0;
  } // adjust_liquid_aerosol

  KOKKOS_INLINE_FUNCTION static
  void adjust_solid_aerosol(const MosaicModelData<DeviceType>& mosaic,
                            const real_type_1d_view_type& aer_solid,
                            const real_type_1d_view_type& aer_liquid,
                            const real_type_1d_view_type& aer_total,
                            const real_type_1d_view_type& electrolyte_solid,
                            const real_type_1d_view_type& electrolyte_liquid,
                            const real_type_1d_view_type& electrolyte_total,
                            const real_type_1d_view_type& epercent_solid,
                            const real_type_1d_view_type& epercent_liquid,
                            const real_type_1d_view_type& epercent_total,
                            real_type& water_a, real_type& jphase, real_type& jhyst_leg) {

    jphase = mosaic.jsolid;
    jhyst_leg = mosaic.jhyst_lo;
    water_a = 0.0;

    // Transfer aer(total) to aer(jsolid) and set aer(liquid) to 0
    for (ordinal_type iaer = 0; iaer < mosaic.naer; iaer++) {
        aer_solid(iaer) = aer_total(iaer);
        aer_liquid(iaer) = 0.0;
    }

    // Transfer electrolyte(total) to electrolyte(solid) and set electrolyte(liquid) to 0
    for (ordinal_type je = 0; je < mosaic.nelectrolyte; je++) {
      electrolyte_liquid(je) = 0.0;
      epercent_liquid(je) = 0.0;
      electrolyte_solid(je) = electrolyte_total(je);
      epercent_solid(je) = epercent_total(je);
    }

    // Update aer(total) that may have been affected above
    aer_total(mosaic.inh4_a) = aer_solid(mosaic.inh4_a);
    aer_total(mosaic.ino3_a) = aer_solid(mosaic.ino3_a);
    aer_total(mosaic.icl_a) = aer_solid(mosaic.icl_a);

    // Update electrolyte(total)
    for (ordinal_type je = 0; je < mosaic.nelectrolyte; je++) {
      electrolyte_total(je) = electrolyte_solid(je);
      epercent_total(je) = epercent_solid(je);
    }
  } // adjust_solid_aerosol

  KOKKOS_INLINE_FUNCTION static
  void do_full_deliquescence(const MosaicModelData<DeviceType>& mosaic,
                             const real_type_1d_view_type& electrolyte_solid,
                             const real_type_1d_view_type& electrolyte_liquid,
                             const real_type_1d_view_type& electrolyte_total,
                             const real_type_1d_view_type& aer_solid,
                             const real_type_1d_view_type& aer_liquid,
                             const real_type_1d_view_type& aer_total) {

    // Partition all electrolytes to the liquid phase
    for (ordinal_type js = 0; js < mosaic.nelectrolyte; js++) {
        electrolyte_solid(js) = 0.0;
        electrolyte_liquid(js) = electrolyte_total(js);
    }

    // Except these electrolytes, which always remain in the solid phase
    electrolyte_solid(mosaic.jcaco3) = electrolyte_total(mosaic.jcaco3);
    electrolyte_solid(mosaic.jcaso4) = electrolyte_total(mosaic.jcaso4);
    electrolyte_liquid(mosaic.jcaco3) = 0.0;
    electrolyte_liquid(mosaic.jcaso4) = 0.0;

    // Partition all the generic aer species into solid and liquid phases
    // Solid phase
    aer_solid(mosaic.iso4_a) = electrolyte_solid(mosaic.jcaso4);
    aer_solid(mosaic.ino3_a) = 0.0;
    aer_solid(mosaic.icl_a) = 0.0;
    aer_solid(mosaic.inh4_a) = 0.0;
    aer_solid(mosaic.ioc_a) = aer_total(mosaic.ioc_a);
    aer_solid(mosaic.imsa_a) = 0.0;
    aer_solid(mosaic.ico3_a) = aer_total(mosaic.ico3_a);
    aer_solid(mosaic.ina_a) = 0.0;
    aer_solid(mosaic.ica_a) = electrolyte_solid(mosaic.jcaco3) + electrolyte_solid(mosaic.jcaso4);
    aer_solid(mosaic.ibc_a) = aer_total(mosaic.ibc_a);
    aer_solid(mosaic.ioin_a) = aer_total(mosaic.ioin_a);
    aer_solid(mosaic.iaro1_a) = aer_total(mosaic.iaro1_a);
    aer_solid(mosaic.iaro2_a) = aer_total(mosaic.iaro2_a);
    aer_solid(mosaic.ialk1_a) = aer_total(mosaic.ialk1_a);
    aer_solid(mosaic.iole1_a) = aer_total(mosaic.iole1_a);
    aer_solid(mosaic.iapi1_a) = aer_total(mosaic.iapi1_a);
    aer_solid(mosaic.iapi2_a) = aer_total(mosaic.iapi2_a);
    aer_solid(mosaic.ilim1_a) = aer_total(mosaic.ilim1_a);
    aer_solid(mosaic.ilim2_a) = aer_total(mosaic.ilim2_a);

    // Liquid phase
    aer_liquid(mosaic.iso4_a) = aer_total(mosaic.iso4_a) - electrolyte_solid(mosaic.jcaso4);
    aer_liquid(mosaic.ino3_a) = aer_total(mosaic.ino3_a);
    aer_liquid(mosaic.icl_a) = aer_total(mosaic.icl_a);
    aer_liquid(mosaic.inh4_a) = aer_total(mosaic.inh4_a);
    aer_liquid(mosaic.ioc_a) = 0.0;
    aer_liquid(mosaic.imsa_a) = aer_total(mosaic.imsa_a);
    aer_liquid(mosaic.ico3_a) = 0.0;
    aer_liquid(mosaic.ina_a) = aer_total(mosaic.ina_a);
    aer_liquid(mosaic.ica_a) = electrolyte_total(mosaic.jcano3) + electrolyte_total(mosaic.jcacl2);
    aer_liquid(mosaic.ibc_a) = 0.0;
    aer_liquid(mosaic.ioin_a) = 0.0;
    aer_liquid(mosaic.iaro1_a) = 0.0;
    aer_liquid(mosaic.iaro2_a) = 0.0;
    aer_liquid(mosaic.ialk1_a) = 0.0;
    aer_liquid(mosaic.iole1_a) = 0.0;
    aer_liquid(mosaic.iapi1_a) = 0.0;
    aer_liquid(mosaic.iapi2_a) = 0.0;
    aer_liquid(mosaic.ilim1_a) = 0.0;
    aer_liquid(mosaic.ilim2_a) = 0.0;
  } // do_full_deliquescence

  KOKKOS_INLINE_FUNCTION static
  void electrolytes_to_ions(const MosaicModelData<DeviceType>& mosaic,
                            const real_type_1d_view_type& aer,
                            const real_type_1d_view_type& electrolyte,
                            real_type& aer_sum,
                            const real_type_1d_view_type& aer_percent) {
  
    aer(mosaic.iso4_a) = electrolyte(mosaic.jcaso4)  +
                         electrolyte(mosaic.jna2so4) +
                      2.*electrolyte(mosaic.jna3hso4)+
                         electrolyte(mosaic.jnahso4) +
                         electrolyte(mosaic.jnh4so4) +
                      2.*electrolyte(mosaic.jlvcite) +
                         electrolyte(mosaic.jnh4hso4)+
                         electrolyte(mosaic.jh2so4);

    aer(mosaic.ino3_a) = electrolyte(mosaic.jnano3) +
                      2.*electrolyte(mosaic.jcano3) +
                         electrolyte(mosaic.jnh4no3)+
                         electrolyte(mosaic.jhno3);

    aer(mosaic.icl_a)  = electrolyte(mosaic.jnacl) +
                      2.*electrolyte(mosaic.jcacl2)+
                         electrolyte(mosaic.jnh4cl)+
                         electrolyte(mosaic.jhcl);

    aer(mosaic.imsa_a) = electrolyte(mosaic.jnh4msa)+
                         electrolyte(mosaic.jnamsa) +
                      2.*electrolyte(mosaic.jcamsa2)+
                         electrolyte(mosaic.jmsa);

    aer(mosaic.ico3_a) = electrolyte(mosaic.jcaco3);

    aer(mosaic.ica_a)  = electrolyte(mosaic.jcaso4) +
                         electrolyte(mosaic.jcano3) +
                         electrolyte(mosaic.jcacl2) +
                         electrolyte(mosaic.jcaco3) +
                         electrolyte(mosaic.jcamsa2);
    
    aer(mosaic.ina_a)  = electrolyte(mosaic.jnano3)  +
                         electrolyte(mosaic.jnacl)   +
                      2.*electrolyte(mosaic.jna2so4) +
                      3.*electrolyte(mosaic.jna3hso4)+
                         electrolyte(mosaic.jnahso4) +
                         electrolyte(mosaic.jnamsa);

    aer(mosaic.inh4_a) = electrolyte(mosaic.jnh4no3) +
                         electrolyte(mosaic.jnh4cl)  +
                      2.*electrolyte(mosaic.jnh4so4) + 
                      3.*electrolyte(mosaic.jlvcite) +
                         electrolyte(mosaic.jnh4hso4)+
                         electrolyte(mosaic.jnh4msa);

    real_type sum_dum = aer(mosaic.ica_a) +
                        aer(mosaic.ina_a) +
                        aer(mosaic.inh4_a)+
                        aer(mosaic.iso4_a)+
                        aer(mosaic.ino3_a)+
                        aer(mosaic.icl_a) +
                        aer(mosaic.imsa_a)+
                        aer(mosaic.ico3_a);

    if(sum_dum == 0.0) {
      sum_dum = 1.0;
    }

    aer_sum = sum_dum;

    aer_percent(mosaic.ica_a) = 100.*aer(mosaic.ica_a)/sum_dum;
    aer_percent(mosaic.ina_a) = 100.*aer(mosaic.ina_a)/sum_dum;
    aer_percent(mosaic.inh4_a)= 100.*aer(mosaic.inh4_a)/sum_dum;
    aer_percent(mosaic.iso4_a)= 100.*aer(mosaic.iso4_a)/sum_dum;
    aer_percent(mosaic.ino3_a)= 100.*aer(mosaic.ino3_a)/sum_dum;
    aer_percent(mosaic.icl_a) = 100.*aer(mosaic.icl_a)/sum_dum;
    aer_percent(mosaic.imsa_a)= 100.*aer(mosaic.imsa_a)/sum_dum;
    aer_percent(mosaic.ico3_a)= 100.*aer(mosaic.ico3_a)/sum_dum;
  } // electrolytes_to_ions

  KOKKOS_INLINE_FUNCTION static
  void ions_to_electrolytes(const MosaicModelData<DeviceType>& mosaic,
                            const real_type& jp,
                            const real_type_1d_view_type& aer_curr,
                            const real_type_1d_view_type& aer_solid,
                            const real_type_1d_view_type& aer_liquid,
                            const real_type_1d_view_type& aer_total,
                            const real_type_1d_view_type& store,
                            const real_type_1d_view_type& electrolyte_curr,
                            const real_type_1d_view_type& electrolyte_solid,
                            const real_type_1d_view_type& electrolyte_liquid,
                            const real_type_1d_view_type& na,
                            const real_type_1d_view_type& nc,
                            const real_type_1d_view_type& xeq_a,
                            const real_type_1d_view_type& xeq_c,
                            const real_type_1d_view_type& na_Ma,
                            const real_type_1d_view_type& nc_Mc,
                            real_type& electrolyte_sum,
                            const real_type_1d_view_type& epercent,
                            real_type& aer_sum,
                            const real_type_1d_view_type& aer_percent) {

    // Note: this function should only be called on liquid-phase or total-phase
    // transfers caso4 and caco3 from liquid to solid phase

    auto za = mosaic.za.template view<DeviceType>();
    auto zc = mosaic.zc.template view<DeviceType>();
    auto mw_a = mosaic.mw_a.template view<DeviceType>();
    auto mw_c = mosaic.mw_c.template view<DeviceType>();
    auto mw_electrolyte = mosaic.mw_electrolyte.template view<DeviceType>();

    // remove negative concentrations if any 
    for (int i = 0; i < mosaic.naer; i++) {
      aer_curr(i) = max(0.0, aer_curr(i));
    }

    // transfer caso4 from liquid to solid phase (caco3 should not be present here)
    store(mosaic.ica_a)  = aer_curr(mosaic.ica_a);
    store(mosaic.iso4_a) = aer_curr(mosaic.iso4_a);

    form_caso4(mosaic, electrolyte_curr, store);

    if (jp == mosaic.jliquid) { // transfer caso4 from liquid to solid phase
      aer_liquid(mosaic.ica_a) = aer_liquid(mosaic.ica_a) - electrolyte_liquid(mosaic.jcaso4);
      aer_liquid(mosaic.iso4_a) = aer_liquid(mosaic.iso4_a) - electrolyte_liquid(mosaic.jcaso4);
      aer_solid(mosaic.ica_a)  = aer_solid(mosaic.ica_a)  + electrolyte_liquid(mosaic.jcaso4);
      aer_solid(mosaic.iso4_a) = aer_solid(mosaic.iso4_a) + electrolyte_liquid(mosaic.jcaso4);
      electrolyte_solid(mosaic.jcaso4) = electrolyte_solid(mosaic.jcaso4) + electrolyte_liquid(mosaic.jcaso4);
      electrolyte_liquid(mosaic.jcaso4) = 0.0;
    }

    real_type XT = 0.0;
    calculate_XT(mosaic, aer_liquid, XT);

    ordinal_type icase = 0;
    if (XT >= 1.9999 || XT < 0.0) {
      icase = 1; // near neutral (acidity is caused by HCl and/or HNO3)
    } else {
      icase = 2; // acidic (acidity is caused by excess H2SO4)
    }

    // initialize to zero
    for (int je = 0; je < mosaic.nelectrolyte; je++) {
      electrolyte_curr(je) = 0.0;
    }

    // initialize moles of ions depending on the sulfate domain
    if (icase == 1) {
      na(mosaic.ja_hso4) = 0.0;
      na(mosaic.ja_so4) = aer_curr(mosaic.iso4_a);
      na(mosaic.ja_no3) = aer_curr(mosaic.ino3_a);
      na(mosaic.ja_cl)  = aer_curr(mosaic.icl_a);
      na(mosaic.ja_msa) = aer_curr(mosaic.imsa_a);

      nc(mosaic.jc_ca)  = aer_curr(mosaic.ica_a);
      nc(mosaic.jc_na)  = aer_curr(mosaic.ina_a);
      nc(mosaic.jc_nh4) = aer_curr(mosaic.inh4_a);

      real_type cat_net = (
          (2.0*na(mosaic.ja_so4)+na(mosaic.ja_no3)+na(mosaic.ja_cl)+na(mosaic.ja_msa)) -
          (2.0*nc(mosaic.jc_ca) +nc(mosaic.jc_nh4)+nc(mosaic.jc_na)) );

      if (cat_net < 0.0) {
        nc(mosaic.jc_h) = 0.0;
      } else { // cat_net must be 0.0 or positive
        nc(mosaic.jc_h) = cat_net;
      }

      // now compute equivalent fractions
      real_type sum_naza = 0.0;
      for (int ja = 0; ja < mosaic.nanion; ja++) {
        sum_naza += na(ja)*za(ja);
      }

      real_type sum_nczc = 0.0;
      for (int jc = 0; jc < mosaic.ncation; jc++) {
        sum_nczc += nc(jc)*zc(jc);
      }

      if (sum_naza == 0.0 || sum_nczc == 0.0) {
        return;
      }

      for (int ja = 0; ja < mosaic.nanion; ja++) {
        xeq_a(ja) = na(ja)*za(ja)/sum_naza;
      }

      for (int jc = 0; jc < mosaic.ncation; jc++) {
        xeq_c(jc) = nc(jc)*zc(jc)/sum_nczc;
      }

      na_Ma(mosaic.ja_so4) = na(mosaic.ja_so4) * mw_a(mosaic.ja_so4);
      na_Ma(mosaic.ja_no3) = na(mosaic.ja_no3) * mw_a(mosaic.ja_no3);
      na_Ma(mosaic.ja_cl)  = na(mosaic.ja_cl)  * mw_a(mosaic.ja_cl);
      na_Ma(mosaic.ja_msa) = na(mosaic.ja_msa) * mw_a(mosaic.ja_msa);
      na_Ma(mosaic.ja_hso4)= na(mosaic.ja_hso4)* mw_a(mosaic.ja_hso4);

      nc_Mc(mosaic.jc_ca)  = nc(mosaic.jc_ca) * mw_c(mosaic.jc_ca);
      nc_Mc(mosaic.jc_na)  = nc(mosaic.jc_na) * mw_c(mosaic.jc_na);
      nc_Mc(mosaic.jc_nh4) = nc(mosaic.jc_nh4)* mw_c(mosaic.jc_nh4);
      nc_Mc(mosaic.jc_h)   = nc(mosaic.jc_h)  * mw_c(mosaic.jc_h);

      // now compute electrolyte moles
      if(xeq_c(mosaic.jc_na) > 0.0 && xeq_a(mosaic.ja_so4) > 0.0) {
        electrolyte_curr(mosaic.jna2so4) = (xeq_c(mosaic.jc_na) *na_Ma(mosaic.ja_so4) +
                                            xeq_a(mosaic.ja_so4)*nc_Mc(mosaic.jc_na))/
                                            mw_electrolyte(mosaic.jna2so4);
      }

      electrolyte_curr(mosaic.jnahso4) = 0.0;

      if(xeq_c(mosaic.jc_na) > 0.0 && xeq_a(mosaic.ja_msa) > 0.0) {
        electrolyte_curr(mosaic.jnamsa)  = (xeq_c(mosaic.jc_na) *na_Ma(mosaic.ja_msa) +
                                            xeq_a(mosaic.ja_msa)*nc_Mc(mosaic.jc_na))/
                                            mw_electrolyte(mosaic.jnamsa);
      }

      if(xeq_c(mosaic.jc_na) > 0.0 && xeq_a(mosaic.ja_no3) > 0.0) {
        electrolyte_curr(mosaic.jnano3)  = (xeq_c(mosaic.jc_na) *na_Ma(mosaic.ja_no3) +
                                            xeq_a(mosaic.ja_no3)*nc_Mc(mosaic.jc_na))/
                                            mw_electrolyte(mosaic.jnano3);
      }

      if(xeq_c(mosaic.jc_na) > 0.0 && xeq_a(mosaic.ja_cl) > 0.0) {
        electrolyte_curr(mosaic.jnacl)   = (xeq_c(mosaic.jc_na) *na_Ma(mosaic.ja_cl) +
                                            xeq_a(mosaic.ja_cl) *nc_Mc(mosaic.jc_na))/
                                            mw_electrolyte(mosaic.jnacl);
      }

      if(xeq_c(mosaic.jc_nh4) > 0.0 && xeq_a(mosaic.ja_so4) > 0.0) {
        electrolyte_curr(mosaic.jnh4so4) = (xeq_c(mosaic.jc_nh4)*na_Ma(mosaic.ja_so4) +
                                            xeq_a(mosaic.ja_so4)*nc_Mc(mosaic.jc_nh4))/
                                            mw_electrolyte(mosaic.jnh4so4);
      }

      electrolyte_curr(mosaic.jnh4hso4)= 0.0;

      if(xeq_c(mosaic.jc_nh4) > 0.0 && xeq_a(mosaic.ja_msa) > 0.0) {
        electrolyte_curr(mosaic.jnh4msa) = (xeq_c(mosaic.jc_nh4)*na_Ma(mosaic.ja_msa) +
                                            xeq_a(mosaic.ja_msa)*nc_Mc(mosaic.jc_nh4))/
                                            mw_electrolyte(mosaic.jnh4msa);
      }

      if(xeq_c(mosaic.jc_nh4) > 0.0 && xeq_a(mosaic.ja_no3) > 0.0) {
        electrolyte_curr(mosaic.jnh4no3) = (xeq_c(mosaic.jc_nh4)*na_Ma(mosaic.ja_no3) +
                                            xeq_a(mosaic.ja_no3)*nc_Mc(mosaic.jc_nh4))/
                                            mw_electrolyte(mosaic.jnh4no3);
      }

      if(xeq_c(mosaic.jc_nh4) > 0.0 && xeq_a(mosaic.ja_cl) > 0.0) {
        electrolyte_curr(mosaic.jnh4cl)  = (xeq_c(mosaic.jc_nh4)*na_Ma(mosaic.ja_cl) +
                                            xeq_a(mosaic.ja_cl) *nc_Mc(mosaic.jc_nh4))/
                                            mw_electrolyte(mosaic.jnh4cl);
      }

      if(xeq_c(mosaic.jc_ca) > 0.0 && xeq_a(mosaic.ja_no3) > 0.0) {
        electrolyte_curr(mosaic.jcano3) = (xeq_c(mosaic.jc_ca) *na_Ma(mosaic.ja_no3) +
                                           xeq_a(mosaic.ja_no3)*nc_Mc(mosaic.jc_ca))/
                                           mw_electrolyte(mosaic.jcano3);
      }

      if(xeq_c(mosaic.jc_ca) > 0.0 && xeq_a(mosaic.ja_cl) > 0.0) {
        electrolyte_curr(mosaic.jcacl2)  = (xeq_c(mosaic.jc_ca) *na_Ma(mosaic.ja_cl) +
                                            xeq_a(mosaic.ja_cl) *nc_Mc(mosaic.jc_ca))/
                                            mw_electrolyte(mosaic.jcacl2);
      }

      if(xeq_c(mosaic.jc_ca) > 0.0 && xeq_a(mosaic.ja_msa) > 0.0) {
        electrolyte_curr(mosaic.jcamsa2) = (xeq_c(mosaic.jc_ca) *na_Ma(mosaic.ja_msa) +
                                            xeq_a(mosaic.ja_msa)*nc_Mc(mosaic.jc_ca))/
                                            mw_electrolyte(mosaic.jcamsa2);
      }

      electrolyte_curr(mosaic.jh2so4) = 0.0;

      if(xeq_c(mosaic.jc_h) > 0.0 && xeq_a(mosaic.ja_no3) > 0.0) {
        electrolyte_curr(mosaic.jhno3) = (xeq_c(mosaic.jc_h)  *na_Ma(mosaic.ja_no3) +
                                          xeq_a(mosaic.ja_no3)*nc_Mc(mosaic.jc_h))/
                                          mw_electrolyte(mosaic.jhno3);
      }

      if(xeq_c(mosaic.jc_h) > 0.0 && xeq_a(mosaic.ja_cl) > 0.0) {
        electrolyte_curr(mosaic.jhcl) = (xeq_c(mosaic.jc_h) *na_Ma(mosaic.ja_cl) +
                                         xeq_a(mosaic.ja_cl)*nc_Mc(mosaic.jc_h))/
                                         mw_electrolyte(mosaic.jhcl);
      }

      if(xeq_c(mosaic.jc_h) > 0.0 && xeq_a(mosaic.ja_msa) > 0.0) {
        electrolyte_curr(mosaic.jmsa) = (xeq_c(mosaic.jc_h)  *na_Ma(mosaic.ja_msa) +
                                         xeq_a(mosaic.ja_msa)*nc_Mc(mosaic.jc_h))/
                                         mw_electrolyte(mosaic.jmsa);
      }
    } else if (icase == 2) { // XT < 2: SULFATE RICH DOMAIN
      store(mosaic.imsa_a) = aer_curr(mosaic.imsa_a);
      store(mosaic.ica_a)  = aer_curr(mosaic.ica_a);

      form_camsa2(mosaic, electrolyte_curr, store);

      real_type sum_na_nh4 = aer_curr(mosaic.ina_a) + aer_curr(mosaic.inh4_a);

      real_type f_na, f_nh4 = 0.0;
      if(sum_na_nh4 > 0.0) {
        f_na  = aer_curr(mosaic.ina_a)/sum_na_nh4;
        f_nh4 = aer_curr(mosaic.inh4_a)/sum_na_nh4;
      }

      // first form msa electrolytes
      real_type rem_na, rem_nh4 = 0.0;
      if(sum_na_nh4 > store(mosaic.imsa_a)) {
        electrolyte_curr(mosaic.jnamsa)  = f_na *store(mosaic.imsa_a);
        electrolyte_curr(mosaic.jnh4msa) = f_nh4*store(mosaic.imsa_a);
        rem_na = aer_curr(mosaic.ina_a) - electrolyte_curr(mosaic.jnamsa);  // remaining na
        rem_nh4= aer_curr(mosaic.inh4_a)- electrolyte_curr(mosaic.jnh4msa); // remaining nh4
      } else {
        electrolyte_curr(mosaic.jnamsa)  = aer_curr(mosaic.ina_a);
        electrolyte_curr(mosaic.jnh4msa) = aer_curr(mosaic.inh4_a);
        electrolyte_curr(mosaic.jmsa)    = store(mosaic.imsa_a) - sum_na_nh4;
      }


      // recompute XT
      if(aer_curr(mosaic.iso4_a) > 0.0) {
        XT = (rem_nh4 + rem_na)/aer_curr(mosaic.iso4_a);
      } else {
        real_type sum_dum = 0.0;
        for (int je = 0; je < mosaic.nelectrolyte; je++) {
          sum_dum += electrolyte_curr(je);
        }

        if (sum_dum == 0.0) {
          sum_dum = 1.0;
        }
        electrolyte_sum = sum_dum;

        for (int je = 0; je < mosaic.nelectrolyte; je++) {
          epercent(je) = 100.0*electrolyte_curr(je)/sum_dum;
        }

        sum_dum = aer_curr(mosaic.ica_a) +
                  aer_curr(mosaic.ina_a) +
                  aer_curr(mosaic.inh4_a)+
                  aer_curr(mosaic.iso4_a)+
                  aer_curr(mosaic.ino3_a)+
                  aer_curr(mosaic.icl_a) +
                  aer_curr(mosaic.imsa_a)+
                  aer_curr(mosaic.ico3_a);

        if(sum_dum == 0.0) {
          sum_dum = 1.0;
        }
        aer_sum = sum_dum;

        aer_percent(mosaic.ica_a) = 100.0*aer_curr(mosaic.ica_a)/sum_dum;
        aer_percent(mosaic.ina_a) = 100.0*aer_curr(mosaic.ina_a)/sum_dum;
        aer_percent(mosaic.inh4_a)= 100.0*aer_curr(mosaic.inh4_a)/sum_dum;
        aer_percent(mosaic.iso4_a)= 100.0*aer_curr(mosaic.iso4_a)/sum_dum;
        aer_percent(mosaic.ino3_a)= 100.0*aer_curr(mosaic.ino3_a)/sum_dum;
        aer_percent(mosaic.icl_a) = 100.0*aer_curr(mosaic.icl_a)/sum_dum;
        aer_percent(mosaic.imsa_a)= 100.0*aer_curr(mosaic.imsa_a)/sum_dum;
        aer_percent(mosaic.ico3_a)= 100.0*aer_curr(mosaic.ico3_a)/sum_dum;
        return;
      }

      real_type xh, xb, xl, xs = 0.0;
      if (XT <= 1.0) {	// h2so4 + bisulfate
        xh = 1.0 - XT;
        xb = XT;
        electrolyte_curr(mosaic.jh2so4)   = xh*aer_curr(mosaic.iso4_a);
        electrolyte_curr(mosaic.jnh4hso4) = xb*f_nh4*aer_curr(mosaic.iso4_a);
        electrolyte_curr(mosaic.jnahso4)  = xb*f_na *aer_curr(mosaic.iso4_a);
      } else if (XT <= 1.5) {	// bisulfate + letovicite
        xb = 3.0 - 2.0*XT;
        xl = XT - 1.0;
        electrolyte_curr(mosaic.jnh4hso4) = xb*f_nh4*aer_curr(mosaic.iso4_a);
        electrolyte_curr(mosaic.jnahso4)  = xb*f_na *aer_curr(mosaic.iso4_a);
        electrolyte_curr(mosaic.jlvcite)  = xl*f_nh4*aer_curr(mosaic.iso4_a);
        electrolyte_curr(mosaic.jna3hso4) = xl*f_na *aer_curr(mosaic.iso4_a);
      } else { // letovicite + sulfate
        xl = 2.0 - XT;
        xs = 2.0*XT - 3.0;
        electrolyte_curr(mosaic.jlvcite)  = xl*f_nh4*aer_curr(mosaic.iso4_a);
        electrolyte_curr(mosaic.jna3hso4) = xl*f_na *aer_curr(mosaic.iso4_a);
        electrolyte_curr(mosaic.jnh4so4)  = xs*f_nh4*aer_curr(mosaic.iso4_a);
        electrolyte_curr(mosaic.jna2so4)  = xs*f_na *aer_curr(mosaic.iso4_a);
      }

      electrolyte_curr(mosaic.jhno3) = aer_curr(mosaic.ino3_a);
      electrolyte_curr(mosaic.jhcl)  = aer_curr(mosaic.icl_a);
    }

    real_type sum_dum = 0.0;
    for (int je = 0; je < mosaic.nelectrolyte; je++) {
      sum_dum += electrolyte_curr(je);
    }

    if (sum_dum == 0.0) {
      sum_dum = 1.0;
    }
    electrolyte_sum = sum_dum;

    for (int je = 0; je < mosaic.nelectrolyte; je++) {
      epercent(je) = 100.0*electrolyte_curr(je)/sum_dum;
    }

    sum_dum = aer_curr(mosaic.ica_a) +
              aer_curr(mosaic.ina_a) +
              aer_curr(mosaic.inh4_a)+
              aer_curr(mosaic.iso4_a)+
              aer_curr(mosaic.ino3_a)+
              aer_curr(mosaic.icl_a) +
              aer_curr(mosaic.imsa_a)+
              aer_curr(mosaic.ico3_a);

    if(sum_dum == 0.0) {
      sum_dum = 1.0;
    }
    aer_sum = sum_dum;

    aer_percent(mosaic.ica_a) = 100.0*aer_curr(mosaic.ica_a)/sum_dum;
    aer_percent(mosaic.ina_a) = 100.0*aer_curr(mosaic.ina_a)/sum_dum;
    aer_percent(mosaic.inh4_a)= 100.0*aer_curr(mosaic.inh4_a)/sum_dum;
    aer_percent(mosaic.iso4_a)= 100.0*aer_curr(mosaic.iso4_a)/sum_dum;
    aer_percent(mosaic.ino3_a)= 100.0*aer_curr(mosaic.ino3_a)/sum_dum;
    aer_percent(mosaic.icl_a) = 100.0*aer_curr(mosaic.icl_a)/sum_dum;
    aer_percent(mosaic.imsa_a)= 100.0*aer_curr(mosaic.imsa_a)/sum_dum;
    aer_percent(mosaic.ico3_a)= 100.0*aer_curr(mosaic.ico3_a)/sum_dum;

    // Final consistency enforcement
    if (jp == mosaic.jtotal) {
      for (ordinal_type i = 0; i < mosaic.naer; ++i) {
        aer_total(i) = aer_curr(i);
      }
    }
  } // ions_to_electrolytes

  KOKKOS_INLINE_FUNCTION static
  void MESA_dissolve_small_salt(const MosaicModelData<DeviceType>& mosaic,
                                const ordinal_type& js,
                                const real_type_1d_view_type& aer_liquid,
                                const real_type_1d_view_type& aer_solid,
                                const real_type_1d_view_type& electrolyte_solid) {

    auto sum_nh4_solid = [&]() {
      return electrolyte_solid(mosaic.jnh4no3)  +
             electrolyte_solid(mosaic.jnh4cl)   +
         2.0*electrolyte_solid(mosaic.jnh4so4)  +
         3.0*electrolyte_solid(mosaic.jlvcite)  +
             electrolyte_solid(mosaic.jnh4hso4) +
             electrolyte_solid(mosaic.jnh4msa);
    };

    auto sum_so4_solid = [&]() {
      return electrolyte_solid(mosaic.jcaso4)   +
             electrolyte_solid(mosaic.jna2so4)  +
         2.0*electrolyte_solid(mosaic.jna3hso4) +
             electrolyte_solid(mosaic.jnahso4)  +
             electrolyte_solid(mosaic.jnh4so4)  +
         2.0*electrolyte_solid(mosaic.jlvcite)  +
             electrolyte_solid(mosaic.jnh4hso4) +
             electrolyte_solid(mosaic.jh2so4);
    };

    auto sum_na_solid = [&]() {
      return electrolyte_solid(mosaic.jnano3)   +
             electrolyte_solid(mosaic.jnacl)    +
         2.0*electrolyte_solid(mosaic.jna2so4)  +
         3.0*electrolyte_solid(mosaic.jna3hso4) +
             electrolyte_solid(mosaic.jnahso4)  +
             electrolyte_solid(mosaic.jnamsa);
    };

    auto sum_no3_solid = [&]() {
      return electrolyte_solid(mosaic.jnano3)  +
         2.0*electrolyte_solid(mosaic.jcano3)  +
             electrolyte_solid(mosaic.jnh4no3) +
             electrolyte_solid(mosaic.jhno3);
    };

    auto sum_cl_solid = [&]() {
      return electrolyte_solid(mosaic.jnacl)   +
         2.0*electrolyte_solid(mosaic.jcacl2)  +
             electrolyte_solid(mosaic.jnh4cl)  +
             electrolyte_solid(mosaic.jhcl);
    };

    auto sum_ca_solid = [&]() {
      return electrolyte_solid(mosaic.jcaso4)  +
             electrolyte_solid(mosaic.jcano3)  +
             electrolyte_solid(mosaic.jcacl2)  +
             electrolyte_solid(mosaic.jcaco3)  +
             electrolyte_solid(mosaic.jcamsa2);
    };

    if (js == mosaic.jnh4so4) {
      aer_liquid(mosaic.inh4_a) += 2.0*electrolyte_solid(js);
      aer_liquid(mosaic.iso4_a) +=     electrolyte_solid(js);
      electrolyte_solid(js) = 0.0;
      aer_solid(mosaic.inh4_a) = sum_nh4_solid();
      aer_solid(mosaic.iso4_a) = sum_so4_solid();

    } else if (js == mosaic.jlvcite) {
      aer_liquid(mosaic.inh4_a) += 3.0*electrolyte_solid(js);
      aer_liquid(mosaic.iso4_a) += 2.0*electrolyte_solid(js);
      electrolyte_solid(js) = 0.0;
      aer_solid(mosaic.inh4_a) = sum_nh4_solid();
      aer_solid(mosaic.iso4_a) = sum_so4_solid();

    } else if (js == mosaic.jnh4hso4) {
      aer_liquid(mosaic.inh4_a) += electrolyte_solid(js);
      aer_liquid(mosaic.iso4_a) += electrolyte_solid(js);
      electrolyte_solid(js) = 0.0;
      aer_solid(mosaic.inh4_a) = sum_nh4_solid();
      aer_solid(mosaic.iso4_a) = sum_so4_solid();

    } else if (js == mosaic.jna2so4) {
      aer_liquid(mosaic.ina_a)  += 2.0*electrolyte_solid(js);
      aer_liquid(mosaic.iso4_a) +=     electrolyte_solid(js);
      electrolyte_solid(js) = 0.0;
      aer_solid(mosaic.ina_a)  = sum_na_solid();
      aer_solid(mosaic.iso4_a) = sum_so4_solid();

    } else if (js == mosaic.jna3hso4) {
      aer_liquid(mosaic.ina_a)  += 3.0*electrolyte_solid(js);
      aer_liquid(mosaic.iso4_a) += 2.0*electrolyte_solid(js);
      electrolyte_solid(js) = 0.0;
      aer_solid(mosaic.ina_a)  = sum_na_solid();
      aer_solid(mosaic.iso4_a) = sum_so4_solid();

    } else if (js == mosaic.jnahso4) {
      aer_liquid(mosaic.ina_a)  += electrolyte_solid(js);
      aer_liquid(mosaic.iso4_a) += electrolyte_solid(js);
      electrolyte_solid(js) = 0.0;
      aer_solid(mosaic.ina_a)  = sum_na_solid();
      aer_solid(mosaic.iso4_a) = sum_so4_solid();

    } else if (js == mosaic.jnh4no3) {
      aer_liquid(mosaic.inh4_a) += electrolyte_solid(js);
      aer_liquid(mosaic.ino3_a) += electrolyte_solid(js);
      electrolyte_solid(js) = 0.0;
      aer_solid(mosaic.inh4_a) = sum_nh4_solid();
      aer_solid(mosaic.ino3_a) = sum_no3_solid();

    } else if (js == mosaic.jnh4cl) {
      aer_liquid(mosaic.inh4_a) += electrolyte_solid(js);
      aer_liquid(mosaic.icl_a)  += electrolyte_solid(js);
      electrolyte_solid(js) = 0.0;
      aer_solid(mosaic.inh4_a) = sum_nh4_solid();
      aer_solid(mosaic.icl_a)  = sum_cl_solid();

    } else if (js == mosaic.jnano3) {
      aer_liquid(mosaic.ina_a)  += electrolyte_solid(js);
      aer_liquid(mosaic.ino3_a) += electrolyte_solid(js);
      electrolyte_solid(js) = 0.0;
      aer_solid(mosaic.ina_a)  = sum_na_solid();
      aer_solid(mosaic.ino3_a) = sum_no3_solid();

    } else if (js == mosaic.jnacl) {
      aer_liquid(mosaic.ina_a) += electrolyte_solid(js);
      aer_liquid(mosaic.icl_a) += electrolyte_solid(js);
      electrolyte_solid(js) = 0.0;
      aer_solid(mosaic.ina_a) = sum_na_solid();
      aer_solid(mosaic.icl_a) = sum_cl_solid();

    } else if (js == mosaic.jcano3) {
      aer_liquid(mosaic.ica_a)  +=     electrolyte_solid(js);
      aer_liquid(mosaic.ino3_a) += 2.0*electrolyte_solid(js);
      electrolyte_solid(js) = 0.0;
      aer_solid(mosaic.ica_a)  = sum_ca_solid();
      aer_solid(mosaic.ino3_a) = sum_no3_solid();

    } else if (js == mosaic.jcacl2) {
      aer_liquid(mosaic.ica_a) +=     electrolyte_solid(js);
      aer_liquid(mosaic.icl_a) += 2.0*electrolyte_solid(js);
      electrolyte_solid(js) = 0.0;
      aer_solid(mosaic.ica_a) = sum_ca_solid();
      aer_solid(mosaic.icl_a) = sum_cl_solid();
    }

  } // MESA_dissolve_small_salt

  KOKKOS_INLINE_FUNCTION static
  void MESA_estimate_eleliquid(
      const MosaicModelData<DeviceType>& mosaic,
      const real_type_1d_view_type& aer_liquid,   // aer(:,jliquid,ibin)
      const real_type_1d_view_type& eleliquid,    // nelectrolyte
      real_type& electrolyte_sum,                 // electrolyte_sum(jliquid,ibin)
      const real_type_1d_view_type& epercent,     // epercent(:,jliquid,ibin)
      const real_type_1d_view_type& na,           // nanion
      const real_type_1d_view_type& nc,           // ncation
      const real_type_1d_view_type& xeq_a,        // nanion
      const real_type_1d_view_type& xeq_c,        // ncation
      const real_type_1d_view_type& na_Ma,        // nanion
      const real_type_1d_view_type& nc_Mc,        // ncation
      const real_type_1d_view_type& store) {      // naer

    auto za             = mosaic.za.template view<DeviceType>();
    auto zc             = mosaic.zc.template view<DeviceType>();
    auto mw_a           = mosaic.mw_a.template view<DeviceType>();
    auto mw_c           = mosaic.mw_c.template view<DeviceType>();
    auto mw_electrolyte = mosaic.mw_electrolyte.template view<DeviceType>();

    // remove negative concentrations, if any
    for (ordinal_type iaer = 0; iaer < mosaic.naer; ++iaer) {
      aer_liquid(iaer) = max(0.0, aer_liquid(iaer));
    }

    // calculate sulfate ratio
    real_type XT = 0;
    calculate_XT(mosaic, aer_liquid, XT);

    ordinal_type icase;
    if (XT >= 2.0 || XT < 0.0) {
      icase = 1; // near neutral: acidity from HCl and/or HNO3
    } else {
      icase = 2; // acidic: acidity from excess SO4
    }

    for (ordinal_type je = 0; je < mosaic.nelectrolyte; ++je) {
      eleliquid(je) = 0.0;
    }

    if (icase == 1) {
      // SULFATE-POOR DOMAIN

      na(mosaic.ja_hso4) = 0.0;
      na(mosaic.ja_so4)  = aer_liquid(mosaic.iso4_a);
      na(mosaic.ja_no3)  = aer_liquid(mosaic.ino3_a);
      na(mosaic.ja_cl)   = aer_liquid(mosaic.icl_a);
      na(mosaic.ja_msa)  = aer_liquid(mosaic.imsa_a);

      nc(mosaic.jc_h)   = 0.0;
      nc(mosaic.jc_ca)  = aer_liquid(mosaic.ica_a);
      nc(mosaic.jc_na)  = aer_liquid(mosaic.ina_a);
      nc(mosaic.jc_nh4) = aer_liquid(mosaic.inh4_a);

      real_type cat_net =
          (2.0*na(mosaic.ja_so4) + na(mosaic.ja_no3) +
           na(mosaic.ja_cl)      + na(mosaic.ja_msa)) -
          (nc(mosaic.jc_h) + 2.0*nc(mosaic.jc_ca) +
           nc(mosaic.jc_nh4)    + nc(mosaic.jc_na));

      if (cat_net < 0.0) {
        nc(mosaic.jc_h) = 0.0;
      } else {
        nc(mosaic.jc_h) = cat_net;
      }

      real_type sum_naza = 0.0;
      for (ordinal_type ja = 0; ja < mosaic.nanion; ++ja) {
        sum_naza += na(ja) * za(ja);
      }

      real_type sum_nczc = 0.0;
      for (ordinal_type jc = 0; jc < mosaic.ncation; ++jc) {
        sum_nczc += nc(jc) * zc(jc);
      }

      if (sum_naza == 0.0 || sum_nczc == 0.0) {
        return;
      }

      for (ordinal_type ja = 0; ja < mosaic.nanion; ++ja) {
        xeq_a(ja) = na(ja) * za(ja) / sum_naza;
      }
      for (ordinal_type jc = 0; jc < mosaic.ncation; ++jc) {
        xeq_c(jc) = nc(jc) * zc(jc) / sum_nczc;
      }

      na_Ma(mosaic.ja_so4)  = na(mosaic.ja_so4)  * mw_a(mosaic.ja_so4);
      na_Ma(mosaic.ja_no3)  = na(mosaic.ja_no3)  * mw_a(mosaic.ja_no3);
      na_Ma(mosaic.ja_cl)   = na(mosaic.ja_cl)   * mw_a(mosaic.ja_cl);
      na_Ma(mosaic.ja_hso4) = na(mosaic.ja_hso4) * mw_a(mosaic.ja_hso4);
      na_Ma(mosaic.ja_msa)  = na(mosaic.ja_msa)  * mw_a(mosaic.ja_msa);

      nc_Mc(mosaic.jc_ca)  = nc(mosaic.jc_ca)  * mw_c(mosaic.jc_ca);
      nc_Mc(mosaic.jc_na)  = nc(mosaic.jc_na)  * mw_c(mosaic.jc_na);
      nc_Mc(mosaic.jc_nh4) = nc(mosaic.jc_nh4) * mw_c(mosaic.jc_nh4);
      nc_Mc(mosaic.jc_h)   = nc(mosaic.jc_h)   * mw_c(mosaic.jc_h);

      eleliquid(mosaic.jna2so4)  = (xeq_c(mosaic.jc_na)  * na_Ma(mosaic.ja_so4)  +
                                    xeq_a(mosaic.ja_so4)  * nc_Mc(mosaic.jc_na))  /
                                    mw_electrolyte(mosaic.jna2so4);

      eleliquid(mosaic.jnahso4)  = (xeq_c(mosaic.jc_na)  * na_Ma(mosaic.ja_hso4) +
                                    xeq_a(mosaic.ja_hso4) * nc_Mc(mosaic.jc_na))  /
                                    mw_electrolyte(mosaic.jnahso4);

      eleliquid(mosaic.jnamsa)   = (xeq_c(mosaic.jc_na)  * na_Ma(mosaic.ja_msa)  +
                                    xeq_a(mosaic.ja_msa)  * nc_Mc(mosaic.jc_na))  /
                                    mw_electrolyte(mosaic.jnamsa);

      eleliquid(mosaic.jnano3)   = (xeq_c(mosaic.jc_na)  * na_Ma(mosaic.ja_no3)  +
                                    xeq_a(mosaic.ja_no3)  * nc_Mc(mosaic.jc_na))  /
                                    mw_electrolyte(mosaic.jnano3);

      eleliquid(mosaic.jnacl)    = (xeq_c(mosaic.jc_na)  * na_Ma(mosaic.ja_cl)   +
                                    xeq_a(mosaic.ja_cl)   * nc_Mc(mosaic.jc_na))  /
                                    mw_electrolyte(mosaic.jnacl);

      eleliquid(mosaic.jnh4so4)  = (xeq_c(mosaic.jc_nh4) * na_Ma(mosaic.ja_so4)  +
                                    xeq_a(mosaic.ja_so4)  * nc_Mc(mosaic.jc_nh4)) /
                                    mw_electrolyte(mosaic.jnh4so4);

      eleliquid(mosaic.jnh4hso4) = (xeq_c(mosaic.jc_nh4) * na_Ma(mosaic.ja_hso4) +
                                    xeq_a(mosaic.ja_hso4) * nc_Mc(mosaic.jc_nh4)) /
                                    mw_electrolyte(mosaic.jnh4hso4);

      eleliquid(mosaic.jnh4msa)  = (xeq_c(mosaic.jc_nh4) * na_Ma(mosaic.ja_msa)  +
                                    xeq_a(mosaic.ja_msa)  * nc_Mc(mosaic.jc_nh4)) /
                                    mw_electrolyte(mosaic.jnh4msa);

      eleliquid(mosaic.jnh4no3)  = (xeq_c(mosaic.jc_nh4) * na_Ma(mosaic.ja_no3)  +
                                    xeq_a(mosaic.ja_no3)  * nc_Mc(mosaic.jc_nh4)) /
                                    mw_electrolyte(mosaic.jnh4no3);

      eleliquid(mosaic.jnh4cl)   = (xeq_c(mosaic.jc_nh4) * na_Ma(mosaic.ja_cl)   +
                                    xeq_a(mosaic.ja_cl)   * nc_Mc(mosaic.jc_nh4)) /
                                    mw_electrolyte(mosaic.jnh4cl);

      eleliquid(mosaic.jcamsa2)  = (xeq_c(mosaic.jc_ca)  * na_Ma(mosaic.ja_msa)  +
                                    xeq_a(mosaic.ja_msa)  * nc_Mc(mosaic.jc_ca))  /
                                    mw_electrolyte(mosaic.jcamsa2);

      eleliquid(mosaic.jcano3)   = (xeq_c(mosaic.jc_ca)  * na_Ma(mosaic.ja_no3)  +
                                    xeq_a(mosaic.ja_no3)  * nc_Mc(mosaic.jc_ca))  /
                                    mw_electrolyte(mosaic.jcano3);

      eleliquid(mosaic.jcacl2)   = (xeq_c(mosaic.jc_ca)  * na_Ma(mosaic.ja_cl)   +
                                    xeq_a(mosaic.ja_cl)   * nc_Mc(mosaic.jc_ca))  /
                                    mw_electrolyte(mosaic.jcacl2);

      eleliquid(mosaic.jh2so4)   = (xeq_c(mosaic.jc_h)   * na_Ma(mosaic.ja_hso4) +
                                    xeq_a(mosaic.ja_hso4) * nc_Mc(mosaic.jc_h))   /
                                    mw_electrolyte(mosaic.jh2so4);

      eleliquid(mosaic.jhno3)    = (xeq_c(mosaic.jc_h)   * na_Ma(mosaic.ja_no3)  +
                                    xeq_a(mosaic.ja_no3)  * nc_Mc(mosaic.jc_h))   /
                                    mw_electrolyte(mosaic.jhno3);

      eleliquid(mosaic.jhcl)     = (xeq_c(mosaic.jc_h)   * na_Ma(mosaic.ja_cl)   +
                                    xeq_a(mosaic.ja_cl)   * nc_Mc(mosaic.jc_h))   /
                                    mw_electrolyte(mosaic.jhcl);

      eleliquid(mosaic.jmsa)     = (xeq_c(mosaic.jc_h)   * na_Ma(mosaic.ja_msa)  +
                                    xeq_a(mosaic.ja_msa)  * nc_Mc(mosaic.jc_h))   /
                                    mw_electrolyte(mosaic.jmsa);

    } else { // SULFATE-RICH DOMAIN

      for (ordinal_type i = 0; i < mosaic.naer; ++i) store(i) = 0.0;

      store(mosaic.iso4_a) = aer_liquid(mosaic.iso4_a);
      store(mosaic.imsa_a) = aer_liquid(mosaic.imsa_a);
      store(mosaic.inh4_a) = aer_liquid(mosaic.inh4_a);
      store(mosaic.ina_a)  = aer_liquid(mosaic.ina_a);
      store(mosaic.ica_a)  = aer_liquid(mosaic.ica_a);

      eleliquid(mosaic.jcamsa2) = min(store(mosaic.ica_a), 0.5*store(mosaic.imsa_a));
      store(mosaic.ica_a)  = max(0.0, store(mosaic.ica_a)  - eleliquid(mosaic.jcamsa2));
      store(mosaic.imsa_a) = max(0.0, store(mosaic.imsa_a) - 2.0*eleliquid(mosaic.jcamsa2));

      real_type sum_na_nh4 = store(mosaic.ina_a) + store(mosaic.inh4_a);
      real_type f_nh4, f_na;
      if (sum_na_nh4 > 0.0) {
        f_nh4 = store(mosaic.inh4_a) / sum_na_nh4;
        f_na  = store(mosaic.ina_a)  / sum_na_nh4;
      } else {
        f_nh4 = 0.0;
        f_na  = 0.0;
      }

      // form MSA electrolytes first
      if (sum_na_nh4 > store(mosaic.imsa_a)) {
        eleliquid(mosaic.jnh4msa) = f_nh4 * store(mosaic.imsa_a);
        eleliquid(mosaic.jnamsa)  = f_na  * store(mosaic.imsa_a);
        store(mosaic.inh4_a) -= eleliquid(mosaic.jnh4msa);
        store(mosaic.ina_a)  -= eleliquid(mosaic.jnamsa);
      } else {
        eleliquid(mosaic.jnh4msa) = store(mosaic.inh4_a);
        eleliquid(mosaic.jnamsa)  = store(mosaic.ina_a);
        eleliquid(mosaic.jmsa)    = store(mosaic.imsa_a) - sum_na_nh4;
        store(mosaic.inh4_a) = 0.0;
        store(mosaic.ina_a)  = 0.0;
      }

      if (store(mosaic.iso4_a) != 0.0) {

        real_type XT_d  = XT;
        real_type XNa_d = 1.0 + 0.5*store(mosaic.ina_a) / store(mosaic.iso4_a);
        real_type xdum  = store(mosaic.iso4_a) - store(mosaic.inh4_a);

        real_type dum    = 2.0*store(mosaic.iso4_a) - store(mosaic.ina_a);
        real_type XNH4_d;
        if (store(mosaic.inh4_a) > 0.0 && dum > 0.0) {
          XNH4_d = 2.0*store(mosaic.inh4_a) /
                   (2.0*store(mosaic.iso4_a) - store(mosaic.ina_a));
        } else {
          XNH4_d = 0.0;
        }

        if (store(mosaic.inh4_a) > 0.0) {

          if (XT_d >= XNa_d) {
            eleliquid(mosaic.jna2so4) = 0.5*store(mosaic.ina_a);

            if (XNH4_d >= 5.0/3.0) {
              eleliquid(mosaic.jnh4so4) = 1.5*store(mosaic.ina_a)
                                         - 3.0*xdum - store(mosaic.inh4_a);
              eleliquid(mosaic.jlvcite) = 2.0*xdum + store(mosaic.inh4_a)
                                         - store(mosaic.ina_a);
            } else if (XNH4_d >= 1.5) {
              eleliquid(mosaic.jnh4so4) = store(mosaic.inh4_a) / 5.0;
              eleliquid(mosaic.jlvcite) = store(mosaic.inh4_a) / 5.0;
            } else if (XNH4_d >= 1.0) {
              eleliquid(mosaic.jnh4so4)  = store(mosaic.inh4_a) / 6.0;
              eleliquid(mosaic.jlvcite)  = store(mosaic.inh4_a) / 6.0;
              eleliquid(mosaic.jnh4hso4) = store(mosaic.inh4_a) / 6.0;
            }

          } else if (XT_d > 1.0) {
            eleliquid(mosaic.jnh4so4)  = store(mosaic.inh4_a) / 6.0;
            eleliquid(mosaic.jlvcite)  = store(mosaic.inh4_a) / 6.0;
            eleliquid(mosaic.jnh4hso4) = store(mosaic.inh4_a) / 6.0;
            eleliquid(mosaic.jna2so4)  = store(mosaic.ina_a)  / 3.0;
            eleliquid(mosaic.jnahso4)  = store(mosaic.ina_a)  / 3.0;

          } else { // XT_d <= 1.0
            eleliquid(mosaic.jna2so4)  = store(mosaic.ina_a)  / 4.0;
            eleliquid(mosaic.jnahso4)  = store(mosaic.ina_a)  / 2.0;
            eleliquid(mosaic.jlvcite)  = store(mosaic.inh4_a) / 6.0;
            eleliquid(mosaic.jnh4hso4) = store(mosaic.inh4_a) / 2.0;
          }

        } else { // inh4_a == 0

          if (XT_d > 1.0) {
            eleliquid(mosaic.jna2so4) = store(mosaic.ina_a)  - store(mosaic.iso4_a);
            eleliquid(mosaic.jnahso4) = 2.0*store(mosaic.iso4_a) - store(mosaic.ina_a);
          } else {
            eleliquid(mosaic.jna2so4) = store(mosaic.ina_a) / 4.0;
            eleliquid(mosaic.jnahso4) = store(mosaic.ina_a) / 2.0;
          }

        }

      }
    }

    real_type sum_dum = 0.0;
    for (ordinal_type je = 0; je < mosaic.nelectrolyte; ++je) {
      sum_dum += eleliquid(je);
    }

    electrolyte_sum = sum_dum;

    if (sum_dum == 0.0) sum_dum = 1.0;
    for (ordinal_type je = 0; je < mosaic.nelectrolyte; ++je) {
      epercent(je) = 100.0 * eleliquid(je) / sum_dum;
    }

  } // MESA_estimate_eleliquid

  KOKKOS_INLINE_FUNCTION static
  void MESA_flux_salt(const MosaicModelData<DeviceType>& mosaic,
                      const real_type_1d_view_type& aer_liquid,
                      const real_type_1d_view_type& aer_solid,
                      const real_type_1d_view_type& aer_total,
                      const real_type_1d_view_type& electrolyte_solid,
                      const real_type_1d_view_type& electrolyte_liquid,
                      const real_type_1d_view_type& electrolyte_total,
                      const real_type_1d_view_type& Keq_sl,
                      const real_type_1d_view_type& Keq_ll,
                      const real_type_1d_view_type& store,
                      const real_type_1d_view_type& na,
                      const real_type_1d_view_type& nc,
                      const real_type_1d_view_type& xeq_a,
                      const real_type_1d_view_type& xeq_c,
                      const real_type_1d_view_type& na_Ma,
                      const real_type_1d_view_type& nc_Mc,
                      real_type& electrolyte_sum_liq,
                      const real_type_1d_view_type& epercent_liquid,
                      real_type& aer_sum,
                      const real_type_1d_view_type& aer_percent,
                      const real_type_1d_view_type& molalities,
                      const real_type_1d_view_type& xmol,
                      const real_type_1d_view_type& ma,
                      const real_type_1d_view_type& mc,
                      const real_type_1d_view_type& log_gam,
                      const real_type_2d_view_type& log_gamZ,
                      const real_type_1d_view_type& gam,
                      const real_type_1d_view_type& activity,
                      const real_type_1d_view_type& eleliquid,
                      real_type& jaerosolstate,
                      real_type& jphase,
                      real_type& jhyst_leg,
                      real_type& aH2O_a,
                      const real_type_1d_view_type& jsalt_present,
                      const real_type_1d_view_type& flux_sl,
                      const real_type_1d_view_type& phi_salt,
                      const real_type_1d_view_type& sat_ratio,
                      const real_type_1d_view_type& frac_salt_solid,
                      const real_type_1d_view_type& frac_salt_liq,
                      real_type& electrolyte_sum_solid) {

    // compute activitie and water content
    ions_to_electrolytes(mosaic, mosaic.jliquid,
                         aer_liquid, aer_solid, aer_liquid, aer_total,
                         store, electrolyte_liquid, electrolyte_solid, electrolyte_liquid,
                         na, nc, xeq_a, xeq_c, na_Ma, nc_Mc,
                         electrolyte_sum_liq, epercent_liquid, aer_sum, aer_percent);

    compute_activities(mosaic, molalities, xmol, aer_liquid, ma, mc,
                       Keq_ll, electrolyte_solid, electrolyte_liquid, electrolyte_total,
                       log_gam, log_gamZ, gam, activity,
                       jaerosolstate, jphase, jhyst_leg, aH2O_a);
    real_type water_a = 0.0;
    aerosol_water(mosaic,electrolyte_liquid,aH2O_a,molalities,jaerosolstate,jphase,jhyst_leg,water_a);

    activity(mosaic.jna3hso4) = 0.0;

    if (water_a <= 0.0) {
      for (ordinal_type js = 0; js < mosaic.nsalt; ++js) {
        flux_sl(js) = 0.0;
      }
      return;
    }

    MESA_estimate_eleliquid(mosaic, aer_liquid, eleliquid, electrolyte_sum_liq,
                            epercent_liquid, na, nc, xeq_a, xeq_c, na_Ma, nc_Mc, store);

    const real_type calcium = aer_liquid(mosaic.ica_a);

    // calculate % electrolyte composition in the solid and liquid phases
    real_type sum_salt = 0.0;
    for (ordinal_type js = 0; js < mosaic.nsalt; ++js) {
      sum_salt += electrolyte_solid(js);
    }
    electrolyte_sum_solid = sum_salt;
    if (sum_salt == 0.0) sum_salt = 1.0;
    for (ordinal_type js = 0; js < mosaic.nsalt; ++js) {
      frac_salt_solid(js) = electrolyte_solid(js) / sum_salt;
      frac_salt_liq(js)   = epercent_liquid(js) / 100.0;
    }

    // compute salt fluxes
    for (ordinal_type js = 0; js < mosaic.nsalt; ++js) {
      // compute new saturation ratio
      sat_ratio(js) = activity(js) / Keq_sl(js);
      // compute relative driving force
      phi_salt(js)  = (sat_ratio(js) - 1.0) / max(sat_ratio(js), 1.0);

      // check if too little solid-phase salt is trying to dissolve
      if (sat_ratio(js) < 1.0 &&
          frac_salt_solid(js) < 0.01 &&
          frac_salt_solid(js) > 0.0) {
        MESA_dissolve_small_salt(mosaic, js, aer_liquid, aer_solid, electrolyte_solid);
        MESA_estimate_eleliquid(mosaic, aer_liquid, eleliquid, electrolyte_sum_liq,
                                epercent_liquid, na, nc, xeq_a, xeq_c, na_Ma, nc_Mc, store);
        sat_ratio(js) = activity(js) / Keq_sl(js);
      }

      // compute flux
      flux_sl(js) = sat_ratio(js) - 1.0;

      // apply Heaviside function
      if ((sat_ratio(js) < 1.0 && electrolyte_solid(js) == 0.0) ||
          (calcium > 0.0 && frac_salt_liq(js) < 0.01) ||
          (calcium > 0.0 && jsalt_present(js) == 0)) {
        flux_sl(js)  = 0.0;
        phi_salt(js) = 0.0;
      }
    }

    // force cacl2 and cano3 fluxes to zero
    sat_ratio(mosaic.jcano3) = 1.0;
    phi_salt(mosaic.jcano3)  = 0.0;
    flux_sl(mosaic.jcano3)   = 0.0;

    sat_ratio(mosaic.jcacl2) = 1.0;
    phi_salt(mosaic.jcacl2)  = 0.0;
    flux_sl(mosaic.jcacl2)   = 0.0;

  } // MESA_flux_salt

#endif
