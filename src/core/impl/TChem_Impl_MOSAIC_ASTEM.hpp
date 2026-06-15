#ifndef __TCHEM_IMPL_MOSAIC_ASTEM_HPP__
#define __TCHEM_IMPL_MOSAIC_ASTEM_HPP__

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

#endif
