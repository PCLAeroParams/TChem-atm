#ifndef __TCHEM_IMPL_MOSAIC_ASTEM_FLUX_SOA_HPP__
#define __TCHEM_IMPL_MOSAIC_ASTEM_FLUX_SOA_HPP__

  //--------------------------------------------------------------------
  // SOA module
  //--------------------------------------------------------------------
  // part of ASTEM: computes fluxes of soa species
  KOKKOS_INLINE_FUNCTION static
  void ASTEM_flux_soa(const MosaicModelData<DeviceType>& mosaic,
                      const real_type_1d_view_type& gas,              // (ngas_volatile)
                      const real_type_1d_view_type& aer_total,        // (naer)
                      const real_type_1d_view_type& sat_soa,          // (ngas_volatile)
                      const real_type_1d_view_type& sfc_a,            // (ngas_volatile)
                      const real_type_1d_view_type& df_gas_o,         // (ngas_volatile)
                      const real_type_1d_view_type& flux_o,           // (ngas_volatile)
                      const real_type_1d_view_type& phi_volatile_o,   // (ngas_volatile)
                      const real_type_1d_view_type& Heff,             // (ngas_volatile)
                      const real_type_1d_view_type& integrate_total,  // (ngas_volatile)
                      real_type& ieqblm_soa) {

    real_type dum, sum_dum, sum_soa;

    // MOSAIC also sets small_oc = 1.e-15 (ng/m^3) here, but never uses it

    // default fluxes and other stuff
    for (ordinal_type iv = mosaic.iaro1_g; iv < mosaic.ngas_volatile; ++iv) {
      sfc_a(iv)          = gas(iv);
      df_gas_o(iv)       = 0.0;
      flux_o(iv)         = 0.0;
      phi_volatile_o(iv) = 0.0;
    }

    // jp = jtotal: aer(:,jp,ibin) is aer_total, integrate(:,jp,ibin) is integrate_total

    // compute mole fractions of soa species
    sum_soa = 0.0;
    for (ordinal_type iv = mosaic.iaro1_g; iv < mosaic.ngas_volatile; ++iv) {
      sum_soa = sum_soa + aer_total(iv);
    }
    sum_soa = sum_soa + aer_total(mosaic.ioc_a)/200.0; // 200 is assumed MW of primary OC

    // check threshold concentration for SOA formation in the absence of primary OC
    if (aer_total(mosaic.ioc_a) == 0.0) {
      sum_dum = 0.0;
      for (ordinal_type iv = mosaic.iaro1_g; iv < mosaic.ngas_volatile; ++iv) {
        sum_dum = sum_dum + (gas(iv) + aer_total(iv))/sat_soa(iv);
      }

      if (sum_dum <= 1.0) { // transfer all aer to gas and quit
        for (ordinal_type iv = mosaic.iaro1_g; iv < mosaic.ngas_volatile; ++iv) {
          gas(iv)       = gas(iv) + aer_total(iv);
          aer_total(iv) = 0.0;
        }
        return;
      }

      sum_soa = max(sum_soa, 1.0e-10);

    }

    // compute Heff
    for (ordinal_type iv = mosaic.iaro1_g; iv < mosaic.ngas_volatile; ++iv) {

      Heff(iv) = sat_soa(iv)/sum_soa;
      sfc_a(iv) = aer_total(iv)*Heff(iv); // nmol/m^3
      df_gas_o(iv) = gas(iv) - sfc_a(iv);

      dum = max(sfc_a(iv), gas(iv));
      if (dum > 0.0) {
        phi_volatile_o(iv) = df_gas_o(iv)/dum;
      } else {
        phi_volatile_o(iv) = 0.0;
      }

      // check equilibrium
      if (Kokkos::abs(phi_volatile_o(iv)) <= mosaic.rtol_eqb_astem) {
        integrate_total(iv) = 0.0; // mNO
      } else {
        integrate_total(iv) = 1.0; // mYES
        ieqblm_soa = 0.0; // mNO
      }

    }

    return;
  } // ASTEM_flux_soa

#endif
