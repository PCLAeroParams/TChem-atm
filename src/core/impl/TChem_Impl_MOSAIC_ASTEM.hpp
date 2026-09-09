#ifndef __TCHEM_IMPL_MOSAIC_ASTEM_HPP__
#define __TCHEM_IMPL_MOSAIC_ASTEM_HPP__

  KOKKOS_INLINE_FUNCTION static
  void aerosol_water_up(const MosaicModelData<DeviceType>& mosaic,
                        const real_type_1d_view_type& electrolyte_total,
                        real_type& aerosol_water) {

    real_type dum = 0.0;

    // TODO: change this for loop to a parallel_reduce
    for (ordinal_type je = 0; je < mosaic.nsalt + 4; je++) {
      real_type molality = 0.0;
      bin_molality_60(mosaic, je, molality);
      dum += 1.e-9*electrolyte_total(je)/molality;
    }

    aerosol_water = dum;
  } // aerosol_water_up

  KOKKOS_INLINE_FUNCTION static
  void degas_nh3(const MosaicModelData<DeviceType>& mosaic,
                 const real_type& jp,
                 const real_type_1d_view_type& electrolyte,
                 const real_type_1d_view_type& store,
                 const real_type_1d_view_type& aer_curr,
                 const real_type_1d_view_type& aer_solid,
                 const real_type_1d_view_type& aer_liquid,
                 const real_type_1d_view_type& aer_total,
                 const real_type_1d_view_type& gas) {

    store(mosaic.inh4_a) = max(0.0, store(mosaic.inh4_a));

    gas(mosaic.inh3_g) = gas(mosaic.inh3_g) + store(mosaic.inh4_a);
  
    aer_curr(mosaic.inh4_a) = (aer_curr(mosaic.inh4_a)) - (store(mosaic.inh4_a));
    aer_curr(mosaic.inh4_a) = max(0.0, aer_curr(mosaic.inh4_a));

  // also do it for jtotal
    if(jp != mosaic.jtotal) {
      aer_total(mosaic.inh4_a) = aer_solid(mosaic.inh4_a) + aer_liquid(mosaic.inh4_a);
    }

    store(mosaic.inh4_a) = 0.0;
  } // degas_nh3

  KOKKOS_INLINE_FUNCTION static
  void degas_hcl(const MosaicModelData<DeviceType>& mosaic,
                 const real_type& jp,
                 const real_type_1d_view_type& electrolyte,
                 const real_type_1d_view_type& store,
                 const real_type_1d_view_type& aer_curr,
                 const real_type_1d_view_type& aer_solid,
                 const real_type_1d_view_type& aer_liquid,
                 const real_type_1d_view_type& aer_total,
                 const real_type_1d_view_type& gas) {

    store(mosaic.icl_a) = max(0.0, store(mosaic.icl_a));

    gas(mosaic.ihcl_g) = gas(mosaic.ihcl_g) + store(mosaic.icl_a);

    aer_curr(mosaic.icl_a) = (aer_curr(mosaic.icl_a)) - (store(mosaic.icl_a));
    aer_curr(mosaic.icl_a) = max(0.0, aer_curr(mosaic.icl_a));

    // also do it for jtotal
    if(jp != mosaic.jtotal) {
      aer_curr(mosaic.icl_a) = aer_solid(mosaic.icl_a) + aer_liquid(mosaic.icl_a);
    }

    electrolyte(mosaic.jhcl) = 0.0;
    store(mosaic.icl_a) = 0.0;
  } // degas_hcl

  KOKKOS_INLINE_FUNCTION static
  void degas_hno3(const MosaicModelData<DeviceType>& mosaic,
                  const real_type& jp,
                  const real_type_1d_view_type& electrolyte,
                  const real_type_1d_view_type& store,
                  const real_type_1d_view_type& aer_curr,
                  const real_type_1d_view_type& aer_solid,
                  const real_type_1d_view_type& aer_liquid,
                  const real_type_1d_view_type& aer_total,
                  const real_type_1d_view_type& gas) {
    
    store(mosaic.ino3_a) = max(0.0, store(mosaic.ino3_a));
    
    gas(mosaic.ihno3_g) = gas(mosaic.ihno3_g) + store(mosaic.ino3_a);

    aer_curr(mosaic.ino3_a) = (aer_curr(mosaic.ino3_a)) - (store(mosaic.ino3_a));
    aer_curr(mosaic.ino3_a) = max(0.0, aer_curr(mosaic.ino3_a));

    // also do it for jtotal
    if(jp != mosaic.jtotal) {
      aer_total(mosaic.ino3_a) = aer_solid(mosaic.ino3_a) + aer_liquid(mosaic.ino3_a);
    }

    electrolyte(mosaic.jhno3) = 0.0;
    store(mosaic.ino3_a) = 0.0;
  } // degas_hno3

  KOKKOS_INLINE_FUNCTION static
  void conform_electrolytes(const MosaicModelData<DeviceType>& mosaic,
                            const real_type& jp,
                            const real_type_1d_view_type& aer_curr,
                            const real_type_1d_view_type& aer_solid,
                            const real_type_1d_view_type& aer_liquid,
                            const real_type_1d_view_type& aer_total,
                            const real_type_1d_view_type& gas,
                            const real_type_1d_view_type& store,
                            const real_type_1d_view_type& electrolyte,
                            const real_type_1d_view_type& total_species,
                            real_type& tot_cl_in,
                            real_type& electrolyte_sum,
                            const real_type_1d_view_type& epercent,
                            real_type& aer_sum,
                            const real_type_1d_view_type& aer_percent) {

    // remove negative concentrations if any 
    for (int i = 0; i < mosaic.naer; i++) {
      aer_curr(i) = max(0.0, aer_curr(i));
    }

    real_type XT = 0.0;
    calculate_XT(mosaic, aer_liquid, XT);

    ordinal_type iXT_case = 0;
    if (XT >= 1.9999 || XT < 0.0) {
      iXT_case = 1; // near neutral (acidity is caused by HCl and/or HNO3)
    } else {
      iXT_case = 2; // acidic (acidity is caused by excess H2SO4)
    }

    // initialize 
    // put aer(*) into store(*)
    store(mosaic.iso4_a) = aer_curr(mosaic.iso4_a);
    store(mosaic.ino3_a) = aer_curr(mosaic.ino3_a);
    store(mosaic.icl_a)  = aer_curr(mosaic.icl_a);
    store(mosaic.imsa_a) = aer_curr(mosaic.imsa_a);
    store(mosaic.ico3_a) = aer_curr(mosaic.ico3_a);
    store(mosaic.inh4_a) = aer_curr(mosaic.inh4_a);
    store(mosaic.ina_a)  = aer_curr(mosaic.ina_a);
    store(mosaic.ica_a)  = aer_curr(mosaic.ica_a);

    for (int j = 0; j < mosaic.nelectrolyte; j++) {
      electrolyte(j) = 0.0;
    }

    if (iXT_case == 1) {
      form_caso4(mosaic, electrolyte, store);
      form_camsa2(mosaic, electrolyte, store);
      form_na2so4(mosaic, electrolyte, store);
      form_namsa(mosaic, electrolyte, store);
      form_cano3(mosaic, electrolyte, store);
      form_nano3(mosaic, electrolyte, store);
      form_nacl(mosaic, jp, electrolyte, store, aer_curr, aer_solid, aer_liquid, aer_total, gas, total_species, tot_cl_in);
      form_cacl2(mosaic, electrolyte, store);
      form_caco3(mosaic, jp, electrolyte, aer_curr, store);
      form_nh4so4(mosaic, electrolyte, store);
      form_nh4msa(mosaic, electrolyte, store);
      form_nh4no3(mosaic, electrolyte, store);
      form_nh4cl(mosaic, electrolyte, store);
      form_msa(mosaic, electrolyte, store);
      degas_hno3(mosaic, jp, electrolyte, store, aer_curr, aer_solid, aer_liquid, aer_total, gas);
      degas_hcl(mosaic, jp, electrolyte, store, aer_curr, aer_solid, aer_liquid, aer_total, gas);
      degas_nh3(mosaic, jp, electrolyte, store, aer_curr, aer_solid, aer_liquid, aer_total, gas);
    } else if (iXT_case == 2) {
      form_caso4(mosaic, electrolyte, store);
      form_camsa2(mosaic, electrolyte, store);
      form_namsa(mosaic, electrolyte, store);
      form_nh4msa(mosaic, electrolyte, store);
      form_msa(mosaic, electrolyte, store);

      if (store(mosaic.iso4_a) == 0.0) {
        degas_hno3(mosaic, jp, electrolyte, store, aer_curr, aer_solid, aer_liquid, aer_total, gas);
        degas_hcl(mosaic, jp, electrolyte, store, aer_curr, aer_solid, aer_liquid, aer_total, gas);
        degas_nh3(mosaic, jp, electrolyte, store, aer_curr, aer_solid, aer_liquid, aer_total, gas);
      } else {
        real_type XT_prime = (store(mosaic.ina_a) + store(mosaic.inh4_a)) / store(mosaic.iso4_a);
        real_type XNa_prime = 0.5*store(mosaic.ina_a) / store(mosaic.iso4_a) + 1.0;
        if (XT_prime >= XNa_prime) {
          form_na2so4(mosaic, electrolyte, store);
          real_type XNH4_prime = 0.0;
          if (store(mosaic.iso4_a) > 1e-15) {
            XNH4_prime = store(mosaic.inh4_a) / store(mosaic.iso4_a);
          }

          if (XNH4_prime > 1.5) {
            form_nh4so4_lvcite(mosaic, electrolyte, store);
          } else {
            form_lvcite_nh4hso4(mosaic, electrolyte, store);
          }

        } else if (XT_prime > 1) {
          form_nh4hso4(mosaic, electrolyte, store);
          form_na2so4_nahso4(mosaic, electrolyte, store);
        } else if (XT_prime < 1) {
          form_nahso4(mosaic, electrolyte, store);
          form_nh4hso4(mosaic, electrolyte, store);
          form_h2so4(mosaic, electrolyte, store);
        }
      }

      degas_hno3(mosaic, jp, electrolyte, store, aer_curr, aer_solid, aer_liquid, aer_total, gas);
      degas_hcl(mosaic, jp, electrolyte, store, aer_curr, aer_solid, aer_liquid, aer_total, gas);
      degas_nh3(mosaic, jp, electrolyte, store, aer_curr, aer_solid, aer_liquid, aer_total, gas);
    }

    electrolytes_to_ions(mosaic, aer_curr, electrolyte, aer_sum, aer_percent);

    real_type sum_dum = 0.0;
    for (int je = 0; je < mosaic.nelectrolyte; je++) {
      electrolyte(je) = max(0.0, electrolyte(je));
      sum_dum += electrolyte(je);
    }

    if (sum_dum == 0.0) {
      sum_dum = 1.0;
    }
    electrolyte_sum = sum_dum;

    for (int je = 0; je < mosaic.nelectrolyte; je++) {
      epercent(je) = 100.0*electrolyte(je)/sum_dum;
    }

    sum_dum = aer_curr(mosaic.ica_a)  +
              aer_curr(mosaic.ina_a)  +
              aer_curr(mosaic.inh4_a) +
              aer_curr(mosaic.iso4_a) +
              aer_curr(mosaic.ino3_a) +
              aer_curr(mosaic.icl_a)  +
              aer_curr(mosaic.imsa_a) +
              aer_curr(mosaic.ico3_a);
    
    if (sum_dum == 0.0) {
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
  } // conform_electrolytes

  KOKKOS_INLINE_FUNCTION static
  void form_electrolytes(const MosaicModelData<DeviceType>& mosaic,
                         const real_type& jp,
                         const real_type_1d_view_type& aer_curr,
                         const real_type_1d_view_type& aer_solid,
                         const real_type_1d_view_type& aer_liquid,
                         const real_type_1d_view_type& aer_total,
                         const real_type_1d_view_type& gas,
                         const real_type_1d_view_type& store,
                         const real_type_1d_view_type& electrolyte,
                         const real_type_1d_view_type& total_species,
                         real_type& tot_cl_in,
                         real_type& electrolyte_sum,
                         const real_type_1d_view_type& epercent,
                         real_type& aer_sum,
                         const real_type_1d_view_type& aer_percent) {
             
    // remove negative concentrations if any 
    for (int i = 0; i < mosaic.naer; i++) {
      aer_curr(i) = max(0.0, aer_curr(i));
    }

    real_type XT = 0.0;
    calculate_XT(mosaic, aer_liquid, XT);

    ordinal_type iXT_case = 0;
    if (XT >= 1.9999 || XT < 0.0) {
      iXT_case = 1; // near neutral (acidity is caused by HCl and/or HNO3)
    } else {
      iXT_case = 2; // acidic (acidity is caused by excess H2SO4)
    }

    // initialize 
    // put aer(*) into store(*)
    store(mosaic.iso4_a) = aer_curr(mosaic.iso4_a);
    store(mosaic.ino3_a) = aer_curr(mosaic.ino3_a);
    store(mosaic.icl_a)  = aer_curr(mosaic.icl_a);
    store(mosaic.imsa_a) = aer_curr(mosaic.imsa_a);
    store(mosaic.ico3_a) = aer_curr(mosaic.ico3_a);
    store(mosaic.inh4_a) = aer_curr(mosaic.inh4_a);
    store(mosaic.ina_a)  = aer_curr(mosaic.ina_a);
    store(mosaic.ica_a)  = aer_curr(mosaic.ica_a);

    for (int j = 0; j < mosaic.nelectrolyte; j++) {
      electrolyte(j) = 0.0;
    }

    // XT >= 2: sulfate deficient
    if (iXT_case == 1) {
      form_caso4(mosaic, electrolyte, store);
      form_camsa2(mosaic, electrolyte, store);
      form_na2so4(mosaic, electrolyte, store);
      form_namsa(mosaic, electrolyte, store);
      form_cano3(mosaic, electrolyte, store);
      form_nano3(mosaic, electrolyte, store);
      form_nacl(mosaic, jp, electrolyte, store, aer_curr, aer_solid, aer_liquid, aer_total, gas, total_species, tot_cl_in);
      form_cacl2(mosaic, electrolyte, store);
      form_caco3(mosaic, jp, electrolyte, aer_curr, store);
      form_nh4so4(mosaic, electrolyte, store);
      form_nh4msa(mosaic, electrolyte, store);
      form_nh4no3(mosaic, electrolyte, store);
      form_nh4cl(mosaic, electrolyte, store);
      form_msa(mosaic, electrolyte, store);

      if (jp == mosaic.jsolid) {
        degas_hno3(mosaic, jp, electrolyte, store, aer_curr, aer_solid, aer_liquid, aer_total, gas);
        degas_hcl(mosaic, jp, electrolyte, store, aer_curr, aer_solid, aer_liquid, aer_total, gas);
        degas_nh3(mosaic, jp, electrolyte, store, aer_curr, aer_solid, aer_liquid, aer_total, gas);
      } else {
        form_hno3(mosaic, electrolyte, store);
        form_hcl(mosaic, electrolyte, store);
        degas_nh3(mosaic, jp, electrolyte, store, aer_curr, aer_solid, aer_liquid, aer_total, gas);
      }

    } else if (iXT_case == 2) {
      // XT < 2: sulfate enough or sulfate excess
      form_caso4(mosaic, electrolyte, store);
      form_camsa2(mosaic, electrolyte, store);
      form_namsa(mosaic, electrolyte, store);
      form_nh4msa(mosaic, electrolyte, store);
      form_msa(mosaic, electrolyte, store);

      if (store(mosaic.iso4_a) == 0.0) {
        if (jp == mosaic.jsolid) {
          degas_hno3(mosaic, jp, electrolyte, store, aer_curr, aer_solid, aer_liquid, aer_total, gas);
          degas_hcl(mosaic, jp, electrolyte, store, aer_curr, aer_solid, aer_liquid, aer_total, gas);
          degas_nh3(mosaic, jp, electrolyte, store, aer_curr, aer_solid, aer_liquid, aer_total, gas);
        } else {
          form_hno3(mosaic, electrolyte, store);
          form_hcl(mosaic, electrolyte, store);
          degas_nh3(mosaic, jp, electrolyte, store, aer_curr, aer_solid, aer_liquid, aer_total, gas);
        }
      } else {
        real_type XT_prime = (store(mosaic.ina_a) + store(mosaic.inh4_a)) / store(mosaic.iso4_a);
        real_type XNa_prime = 0.5*store(mosaic.ina_a) / store(mosaic.iso4_a) + 1.0;
        if (XT_prime >= XNa_prime) {
          form_na2so4(mosaic, electrolyte, store);
          real_type XNH4_prime = 0.0;
          if (store(mosaic.iso4_a) > 1e-15) {
            XNH4_prime = store(mosaic.inh4_a) / store(mosaic.iso4_a);
          }

          if (XNH4_prime > 1.5) {
            form_nh4so4_lvcite(mosaic, electrolyte, store);
          } else {
            form_lvcite_nh4hso4(mosaic, electrolyte, store);
          }

        } else if (XT_prime > 1) {
          form_nh4hso4(mosaic, electrolyte, store);
          form_na2so4_nahso4(mosaic, electrolyte, store);
        } else if (XT_prime < 1) {
          form_nahso4(mosaic, electrolyte, store);
          form_nh4hso4(mosaic, electrolyte, store);
          form_h2so4(mosaic, electrolyte, store);
        }
      }

      if (jp == mosaic.jsolid) {
        degas_hno3(mosaic, jp, electrolyte, store, aer_curr, aer_solid, aer_liquid, aer_total, gas);
        degas_hcl(mosaic, jp, electrolyte, store, aer_curr, aer_solid, aer_liquid, aer_total, gas);
        degas_nh3(mosaic, jp, electrolyte, store, aer_curr, aer_solid, aer_liquid, aer_total, gas);
      } else {
        form_hno3(mosaic, electrolyte, store);
        form_hcl(mosaic, electrolyte, store);
        degas_nh3(mosaic, jp, electrolyte, store, aer_curr, aer_solid, aer_liquid, aer_total, gas);
      }
    }

    electrolytes_to_ions(mosaic, aer_curr, electrolyte, aer_sum, aer_percent);

    real_type sum_dum = 0.0;
    for (int je = 0; je < mosaic.nelectrolyte; je++) {
      electrolyte(je) = max(0.0, electrolyte(je));
      sum_dum += electrolyte(je);
    }

    if (sum_dum == 0.0) {
      sum_dum = 1.0;
    }
    electrolyte_sum = sum_dum;

    for (int je = 0; je < mosaic.nelectrolyte; je++) {
      epercent(je) = 100.0*electrolyte(je)/sum_dum;
    }

    sum_dum = aer_curr(mosaic.ica_a)  +
              aer_curr(mosaic.ina_a)  +
              aer_curr(mosaic.inh4_a) +
              aer_curr(mosaic.iso4_a) +
              aer_curr(mosaic.ino3_a) +
              aer_curr(mosaic.icl_a)  +
              aer_curr(mosaic.imsa_a) +
              aer_curr(mosaic.ico3_a);
    
    if (sum_dum == 0.0) {
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
  } // form_electrolytes

  KOKKOS_INLINE_FUNCTION static
  void absorb_tiny_nh4cl(const MosaicModelData<DeviceType>& mosaic,
                         const real_type_1d_view_type& aer_solid,
                         const real_type_1d_view_type& aer_liquid,
                         const real_type_1d_view_type& aer_total,
                         const real_type_1d_view_type& gas,
                         const real_type& delta_nh3_max,
                         const real_type& delta_hcl_max,
                         const real_type& electrolyte_sum) {

    real_type small_gas = 0.01 * min(delta_nh3_max, delta_hcl_max);
    real_type small_aer = 0.01 * electrolyte_sum;
    if (small_aer == 0.0) {
      small_aer = small_gas;
    }

    real_type small_amt = min(small_gas, small_aer);

    aer_liquid(mosaic.inh4_a) = aer_liquid(mosaic.inh4_a) + small_amt;
    aer_liquid(mosaic.icl_a)  = aer_liquid(mosaic.icl_a)  + small_amt;

    aer_total(mosaic.inh4_a) = aer_solid(mosaic.inh4_a) + aer_liquid(mosaic.inh4_a);
    aer_total(mosaic.icl_a)  = aer_solid(mosaic.icl_a)  + aer_liquid(mosaic.icl_a);

    gas(mosaic.inh3_g) = gas(mosaic.inh3_g) - small_amt;
    gas(mosaic.ihcl_g) = gas(mosaic.ihcl_g) - small_amt;
  } // absorb_tiny_nh4cl
  
  KOKKOS_INLINE_FUNCTION static
  void absorb_tiny_nh4no3(const MosaicModelData<DeviceType>& mosaic,
                          const real_type_1d_view_type& aer_solid,
                          const real_type_1d_view_type& aer_liquid,
                          const real_type_1d_view_type& aer_total,
                          const real_type_1d_view_type& gas,
                          const real_type& delta_nh3_max,
                          const real_type& delta_hno3_max,
                          const real_type& electrolyte_sum_total) {

    real_type small_gas = 0.01 * min(delta_nh3_max, delta_hno3_max);
    real_type small_aer = 0.01 * electrolyte_sum_total;
    if (small_aer == 0.0) {
      small_aer = small_gas;
    }

    real_type small_amt = min(small_gas, small_aer);

    aer_liquid(mosaic.inh4_a) = aer_liquid(mosaic.inh4_a) + small_amt;
    aer_liquid(mosaic.ino3_a) = aer_liquid(mosaic.ino3_a) + small_amt;

    aer_total(mosaic.inh4_a) = aer_solid(mosaic.inh4_a) + aer_liquid(mosaic.inh4_a);
    aer_total(mosaic.ino3_a) = aer_solid(mosaic.ino3_a) + aer_liquid(mosaic.ino3_a);

    gas(mosaic.inh3_g)  = gas(mosaic.inh3_g)  - small_amt;
    gas(mosaic.ihno3_g) = gas(mosaic.ihno3_g) - small_amt;
  } // absorb_tiny_nh4no3

  KOKKOS_INLINE_FUNCTION static
  void degas_solid_nh4no3(const MosaicModelData<DeviceType>& mosaic,
                          const real_type_1d_view_type& gas,
                          const real_type_1d_view_type& Keq_sg,
                          const real_type_1d_view_type& electrolyte_solid,
                          const real_type_1d_view_type& aer_solid,
                          const real_type_1d_view_type& aer_liquid,
                          const real_type_1d_view_type& aer_total) {

    real_type a = 1.0;
    real_type b = gas(mosaic.inh3_g) + gas(mosaic.ihno3_g);
    real_type c = gas(mosaic.inh3_g)*gas(mosaic.ihno3_g) - Keq_sg(0);
    real_type xgas = 0.0;
    quadratic(a, b, c, xgas);

    if (xgas >= electrolyte_solid(mosaic.jnh4no3)) {
      gas(mosaic.inh3_g) = gas(mosaic.inh3_g) + electrolyte_solid(mosaic.jnh4no3);
      gas(mosaic.ihno3_g)= gas(mosaic.ihno3_g) + electrolyte_solid(mosaic.jnh4no3);
      aer_solid(mosaic.inh4_a) = aer_solid(mosaic.inh4_a) - electrolyte_solid(mosaic.jnh4no3);
      aer_solid(mosaic.ino3_a) = aer_solid(mosaic.ino3_a) - electrolyte_solid(mosaic.jnh4no3);
    } else { // degas only xgas amount of nh4no3
      gas(mosaic.inh3_g) = gas(mosaic.inh3_g)  + xgas;
      gas(mosaic.ihno3_g)= gas(mosaic.ihno3_g) + xgas;
      aer_solid(mosaic.inh4_a) = aer_solid(mosaic.inh4_a) - xgas;
      aer_solid(mosaic.ino3_a) = aer_solid(mosaic.ino3_a) - xgas;
    }

    // update jtotal
    aer_total(mosaic.inh4_a) = aer_solid(mosaic.inh4_a) + aer_liquid(mosaic.inh4_a);
    aer_total(mosaic.ino3_a) = aer_solid(mosaic.ino3_a) + aer_liquid(mosaic.ino3_a);
  } // degas_solid_nh4no3

  KOKKOS_INLINE_FUNCTION static
  void degas_solid_nh4cl(const MosaicModelData<DeviceType>& mosaic,
                         const real_type_1d_view_type& gas,
                         const real_type_1d_view_type& Keq_sg,
                         const real_type_1d_view_type& electrolyte_solid,
                         const real_type_1d_view_type& aer_solid,
                         const real_type_1d_view_type& aer_liquid,
                         const real_type_1d_view_type& aer_total) {
  
    real_type a = 1.0;
    real_type b = gas(mosaic.inh3_g) + gas(mosaic.ihcl_g);
    real_type c = gas(mosaic.inh3_g)*gas(mosaic.ihcl_g) - Keq_sg(1);
    real_type xgas = 0.0;
    quadratic(a, b, c, xgas);

    if (xgas >= electrolyte_solid(mosaic.jnh4cl)) {
      gas(mosaic.inh3_g) = gas(mosaic.inh3_g) + electrolyte_solid(mosaic.jnh4cl);
      gas(mosaic.ihcl_g)= gas(mosaic.ihcl_g) + electrolyte_solid(mosaic.jnh4cl);
      aer_solid(mosaic.inh4_a) = aer_solid(mosaic.inh4_a) - electrolyte_solid(mosaic.jnh4cl);
      aer_solid(mosaic.icl_a) = aer_solid(mosaic.icl_a) - electrolyte_solid(mosaic.jnh4cl);
    } else { // degas only xgas amount of nh4cl
      gas(mosaic.inh3_g) = gas(mosaic.inh3_g)  + xgas;
      gas(mosaic.ihcl_g)= gas(mosaic.ihcl_g) + xgas;
      aer_solid(mosaic.inh4_a) = aer_solid(mosaic.inh4_a) - xgas;
      aer_solid(mosaic.icl_a) = aer_solid(mosaic.icl_a) - xgas;
    }

    // update jtotal
    aer_total(mosaic.inh4_a) = aer_solid(mosaic.inh4_a) + aer_liquid(mosaic.inh4_a);
    aer_total(mosaic.icl_a) = aer_solid(mosaic.icl_a) + aer_liquid(mosaic.icl_a);
  } // degas_solid_nh4cl

  KOKKOS_INLINE_FUNCTION static
  void aerosol_phase_state(const MosaicModelData<DeviceType>& mosaic,
                           const real_type_1d_view_type& aer_liquid,
                           const real_type_1d_view_type& aer_solid,
                           const real_type_1d_view_type& aer_total,
                           const real_type_1d_view_type& electrolyte_solid,
                           const real_type_1d_view_type& electrolyte_liquid,
                           const real_type_1d_view_type& electrolyte_total,
                           const real_type_1d_view_type& epercent_liquid,
                           const real_type_1d_view_type& epercent_solid,
                           const real_type_1d_view_type& epercent_total,
                           const real_type_1d_view_type& Keq_sl,
                           const real_type_1d_view_type& Keq_ll,
                           const real_type_2d_view_type& log_gamZ,
                           const real_type_1d_view_type& MDRH_T,
                           real_type& jaerosolstate,
                           real_type& jphase,
                           real_type& jhyst_leg,
                           real_type& aH2O_a,
                           real_type& water_a,
                           real_type& mass_dry_a,
                           real_type& vol_dry_a,
                           real_type& mass_wet_a,
                           real_type& vol_wet_a,
                           real_type& growth_factor,
                           const real_type_1d_view_type& jsalt_present,
                           const real_type_1d_view_type& hsalt,
                           const real_type_1d_view_type& phi_salt,
                           const real_type_1d_view_type& phi_salt_old,
                           const real_type_1d_view_type& phi_bar,
                           const real_type_1d_view_type& alpha_salt,
                           const real_type_1d_view_type& sat_ratio,
                           const real_type_1d_view_type& flux_sl,
                           const real_type_1d_view_type& frac_salt_solid,
                           const real_type_1d_view_type& frac_salt_liq,
                           const real_type_1d_view_type& eleliquid,
                           const real_type_1d_view_type& store,
                           const real_type_1d_view_type& na,
                           const real_type_1d_view_type& nc,
                           const real_type_1d_view_type& xeq_a,
                           const real_type_1d_view_type& xeq_c,
                           const real_type_1d_view_type& na_Ma,
                           const real_type_1d_view_type& nc_Mc,
                           const real_type_1d_view_type& aer_percent,
                           const real_type_1d_view_type& molalities,
                           const real_type_1d_view_type& xmol,
                           const real_type_1d_view_type& ma,
                           const real_type_1d_view_type& mc_view,
                           const real_type_1d_view_type& log_gam,
                           const real_type_1d_view_type& gam,
                           const real_type_1d_view_type& activity,
                           const real_type_1d_view_type& tau_p,
                           const real_type_1d_view_type& tau_d,
                           real_type& electrolyte_sum_solid,
                           real_type& electrolyte_sum_liq,
                           real_type& aer_sum,
                           real_type& kelvin,
                           const real_type_1d_view_type& kel,
                           const real_type& num_a,
                           real_type& sigma_soln,
                           real_type& DpmV,
                           real_type& volume_a,
                           real_type& water_a_up,
                           const real_type& sigma_water,
                           const real_type& T_K,
                           const real_type& RH_pc) {

    auto mw_aer_mac   = mosaic.mw_aer_mac.template view<DeviceType>();
    auto dens_aer_mac = mosaic.dens_aer_mac.template view<DeviceType>();
    auto partial_molar_vol = mosaic.partial_molar_vol.template view<DeviceType>();

    aH2O_a = RH_pc * 0.01;
    kelvin  = 1.0;
    for (ordinal_type iv = 0; iv < mosaic.ngas_volatile; ++iv) {
      kel(iv) = 1.0;
    }

    const real_type kelvin_toler = (RH_pc <= 99.0) ? 1.0e-4 : 1.0e-8;

    real_type aer_H =   2.0*aer_total(mosaic.iso4_a)
                      +     aer_total(mosaic.ino3_a)
                      +     aer_total(mosaic.icl_a)
                      +     aer_total(mosaic.imsa_a)
                      + 2.0*aer_total(mosaic.ico3_a)
                      - 2.0*aer_total(mosaic.ica_a)
                      -     aer_total(mosaic.ina_a)
                      -     aer_total(mosaic.inh4_a);
    aer_H = max(aer_H, 0.0);

    mass_dry_a = 0.0;
    vol_dry_a  = 0.0;
    for (ordinal_type iaer = 0; iaer < mosaic.naer; ++iaer) {
      mass_dry_a += aer_total(iaer) * mw_aer_mac(iaer); // ng/m^3(air)
      vol_dry_a  += aer_total(iaer) * mw_aer_mac(iaer) / dens_aer_mac(iaer); // ncc/m^3(air)
    }
    mass_dry_a = (mass_dry_a + aer_H) * 1.0e-15; // g/cc(air)
    vol_dry_a  = (vol_dry_a  + aer_H) * 1.0e-15; // cc(aer)/cc(air) or m^3/m^3(air)

    mass_wet_a = mass_dry_a + water_a * 1.0e-3; // g/cc(air)
    vol_wet_a  = vol_dry_a  + water_a * 1.0e-3; // cc(aer)/cc(air) or m^3/m^3(air)

    // water uptake at reference RH (60%)
    aerosol_water_up(mosaic, electrolyte_total, water_a_up); // for hysteresis curve determination

    for (ordinal_type iter_kelvin = 1; iter_kelvin < 100; ++iter_kelvin) {

      // binary molalities at current aH2O_a
      for (ordinal_type je = 0; je < mosaic.nelectrolyte; ++je) {
        real_type mol = 0.0;
        bin_molality(mosaic, je, aH2O_a, mol);
        molalities(je) = mol;
      }

      MESA(mosaic,
           aer_liquid, aer_solid, aer_total,
           electrolyte_solid, electrolyte_liquid, electrolyte_total,
           epercent_liquid, epercent_solid, epercent_total,
           Keq_sl, Keq_ll, log_gamZ, MDRH_T,
           jaerosolstate, jphase, jhyst_leg,
           aH2O_a, water_a,
           mass_dry_a, vol_dry_a, mass_wet_a, vol_wet_a,
           growth_factor,
           jsalt_present, hsalt, phi_salt, phi_salt_old, phi_bar,
           alpha_salt, sat_ratio, flux_sl,
           frac_salt_solid, frac_salt_liq,
           eleliquid, store,
           na, nc, xeq_a, xeq_c, na_Ma, nc_Mc,
           aer_percent, molalities, xmol, ma, mc_view,
           log_gam, gam, activity,
           tau_p, tau_d,
           electrolyte_sum_solid, electrolyte_sum_liq, aer_sum);

      if (static_cast<ordinal_type>(jaerosolstate) == mosaic.all_solid) return;

      calculate_kelvin(vol_wet_a, num_a, aH2O_a, sigma_water, T_K,
                       volume_a, DpmV, sigma_soln, kelvin);

      real_type aH2O_a_new = RH_pc*0.01 / kelvin;
      real_type rel_err = ats<real_type>::abs((aH2O_a_new - aH2O_a) / aH2O_a);
      if (rel_err <= kelvin_toler) break;

      aH2O_a = aH2O_a_new;
    }

    if (static_cast<ordinal_type>(jaerosolstate) == mosaic.all_liquid) {
      jhyst_leg = static_cast<real_type>(mosaic.jhyst_up);
    }

    // now compute kelvin effect terms for condensing species (nh3, hno3, and hcl)
    for (ordinal_type iv = 0; iv < mosaic.ngas_volatile; ++iv) {
      real_type term = 4.0*sigma_soln*partial_molar_vol(iv) / (8.3144e7*T_K*DpmV);
      kel(iv) = 1.0 + term*(1.0 + 0.5*term*(1.0 + term/3.0));
    }
  } // aerosol_phase_state

#endif
