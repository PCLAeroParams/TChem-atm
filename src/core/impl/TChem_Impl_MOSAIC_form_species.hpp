#ifndef __TCHEM_IMPL_MOSAIC_FORM_SPECIES_HPP__
#define __TCHEM_IMPL_MOSAIC_FORM_SPECIES_HPP__

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