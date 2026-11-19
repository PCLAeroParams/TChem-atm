#ifndef __TCHEM_IMPL_MOSAIC_ASTEM_FLUX_SOA_HPP__
#define __TCHEM_IMPL_MOSAIC_ASTEM_FLUX_SOA_HPP__

  //--------------------------------------------------------------------
  // SOA module
  //--------------------------------------------------------------------
  // part of ASTEM: condenses secondary organic species over TSI time interval
  // mechanism adapted from SORGAM
  //
  // Every thread of the team calls this routine. Per-particle arrays are (nbin_a, n).
  // The thread-private scalars at the end of the argument list get the same value on every thread.
  template<typename MemberType>
  KOKKOS_INLINE_FUNCTION static
  void ASTEM_secondary_organics(const MemberType& member,
                                const MosaicModelData<DeviceType>& mosaic,
                                const ordinal_type& nbin_a,
                                const real_type& dtchem,
                                const real_type_1d_view_type& gas,              // (ngas_volatile)
                                const real_type_1d_view_type& total_species,    // (ngas_volatile)
                                const real_type_1d_view_type& sat_soa,          // (ngas_volatile)
                                const real_type_2d_view_type& aer_total,        // (nbin_a, naer)
                                const real_type_2d_view_type& kg,               // (nbin_a, ngas_volatile)
                                const real_type_1d_view_type& jaerosolstate,    // (nbin_a)
                                const real_type_2d_view_type& df_gas_o,         // (nbin_a, ngas_volatile)
                                const real_type_2d_view_type& flux_o,           // (nbin_a, ngas_volatile)
                                const real_type_2d_view_type& phi_volatile_o,   // (nbin_a, ngas_volatile)
                                const real_type_2d_view_type& phi_volatile_l,   // (nbin_a, ngas_volatile)
                                const real_type_2d_view_type& Heff,             // (nbin_a, ngas_volatile)
                                const real_type_2d_view_type& integrate_total,  // (nbin_a, ngas_volatile)
                                // work array, one row per particle
                                const real_type_2d_view_type& sfc_a,            // (nbin_a, ngas_volatile)
                                // thread-private
                                real_type& ieqblm_soa,
                                real_type& isteps_SOA) {

    // local variables
    // (MOSAIC also declares nsteps_max = 400 here, but never uses it)
    real_type dtmax = 0.0, t_new = 0.0, t_old = 0.0, t_out = 0.0;

    const auto row = [](const real_type_2d_view_type& v, const ordinal_type ibin) {
      return Kokkos::subview(v, ibin, Kokkos::ALL());
    };

    // initialize time
    t_old = 0.0;
    t_out = dtchem;
    isteps_SOA = 0.0;

    for (ordinal_type iv = mosaic.iaro1_g; iv < mosaic.ngas_volatile; ++iv) {
      // total_species(iv) = gas(iv) plus the sum of aer(iv,jtotal,ibin) over particles
      real_type sum_aer = 0.0;
      Kokkos::parallel_reduce(Kokkos::TeamVectorRange(member, nbin_a),
        [&](const ordinal_type ibin, real_type& sum) {
          sum = sum + aer_total(ibin, iv);
        }, sum_aer);
      Kokkos::single(Kokkos::PerTeam(member), [&]() {
        total_species(iv) = gas(iv) + sum_aer;
      });
    }
    member.team_barrier();

    // overall integration loop begins over dtchem seconds
    while (true) { // label 10
      isteps_SOA = isteps_SOA + 1.0;

      // compute new fluxes
      ieqblm_soa = 1.0; // mYES, reset to default

      // ASTEM_flux_soa can change gas, and MOSAIC's later particles see those
      // changes, so this loop visits the particles in order on one thread.
      Kokkos::single(Kokkos::PerTeam(member), [&](real_type& ieqblm) {
        for (ordinal_type ibin = 0; ibin < nbin_a; ++ibin) { // do 501
          ASTEM_flux_soa(mosaic, gas, row(aer_total, ibin), sat_soa, row(sfc_a, ibin),
                         row(df_gas_o, ibin), row(flux_o, ibin), row(phi_volatile_o, ibin),
                         row(Heff, ibin), row(integrate_total, ibin), ieqblm);
        } // 501
      }, ieqblm_soa);
      member.team_barrier();

      if (ieqblm_soa == 1.0) break; // all bins have reached equilibrium (goto 30)

      // calculate maximum possible internal time-step (label 11)
      ASTEM_dtmax_soa(member, mosaic, nbin_a, dtchem, kg, phi_volatile_l, integrate_total, dtmax);
      t_new = t_old + dtmax; // update time
      if (t_new > t_out) {   // check if the new time step is too large
        dtmax = t_out - t_old;
        t_new = t_out*1.01;
      }

      // do internal time-step (dtmax) integration
      // jp = jtotal: aer(:,jp,:) is aer_total and integrate(:,jp,:) is integrate_total

      for (ordinal_type iv = mosaic.iaro1_g; iv < mosaic.ngas_volatile; ++iv) { // do 20
        real_type sum1 = 0.0, sum2 = 0.0;

        Kokkos::parallel_reduce(Kokkos::TeamVectorRange(member, nbin_a),
          [&](const ordinal_type ibin, real_type& sum) { // do 21
            if (static_cast<ordinal_type>(jaerosolstate(ibin)) == mosaic.no_aerosol) return; // goto 21
            sum = sum + aer_total(ibin, iv)/
                  (1.0 + dtmax*kg(ibin, iv)*Heff(ibin, iv)*integrate_total(ibin, iv));
          }, sum1);
        Kokkos::parallel_reduce(Kokkos::TeamVectorRange(member, nbin_a),
          [&](const ordinal_type ibin, real_type& sum) { // do 21
            if (static_cast<ordinal_type>(jaerosolstate(ibin)) == mosaic.no_aerosol) return; // goto 21
            sum = sum + kg(ibin, iv)*integrate_total(ibin, iv)/
                  (1.0 + dtmax*kg(ibin, iv)*Heff(ibin, iv)*integrate_total(ibin, iv));
          }, sum2);

        // first update gas concentration
        Kokkos::single(Kokkos::PerTeam(member), [&]() {
          gas(iv) = (total_species(iv) - sum1)/(1.0 + dtmax*sum2);
        });
        member.team_barrier();

        // now update aer concentration in the jp phase
        Kokkos::parallel_for(Kokkos::TeamVectorRange(member, nbin_a), [&](const ordinal_type ibin) { // do 22
          if (integrate_total(ibin, iv) == 1.0) { // mYES
            aer_total(ibin, iv) = (aer_total(ibin, iv) + dtmax*kg(ibin, iv)*gas(iv))/
                                  (1.0 + dtmax*kg(ibin, iv)*Heff(ibin, iv));
          }
        });
        member.team_barrier();
      } // 20
      // sub-step integration done

      // (MOSAIC has a commented-out update of jtotal here)

      // update time
      t_old = t_new;

      if (t_new < 0.9999*t_out) continue; // goto 10

      break;
    } // end of integration (label 30)
  } // ASTEM_secondary_organics

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

  //--------------------------------------------------------------------
  // part of ASTEM: computes max time step for the soa integration
  //
  // Every thread of the team calls this routine. The loop over particles
  // (ibin) is a team reduction, whose result every thread receives.
  template<typename MemberType>
  KOKKOS_INLINE_FUNCTION static
  void ASTEM_dtmax_soa(const MemberType& member,
                       const MosaicModelData<DeviceType>& mosaic,
                       const ordinal_type& nbin_a,
                       const real_type& dtchem,
                       const real_type_2d_view_type& kg,               // (nbin_a, ngas_volatile)
                       const real_type_2d_view_type& phi_volatile_l,   // (nbin_a, ngas_volatile)
                       const real_type_2d_view_type& integrate_total,  // (nbin_a, ngas_volatile)
                       real_type& dtmax) {

    // local variables
    real_type h_gas, h_gas_i, h_sub_max, sum_kg_phi;

    h_sub_max = dtchem/6.0; // sec

    // jp = jtotal: integrate(:,jtotal,:) is integrate_total

    // GAS-SIDE
    // calculate h_gas_i and h_gas

    h_gas = 2.e16;

    for (ordinal_type iv = mosaic.iaro1_g; iv < mosaic.ngas_volatile; ++iv) { // do 6
      h_gas_i = 1.e16;
      sum_kg_phi = 0.0;
      Kokkos::parallel_reduce(Kokkos::TeamVectorRange(member, nbin_a),
        [&](const ordinal_type ibin, real_type& sum) {
          if (integrate_total(ibin, iv) == 1.0) { // mYES
            sum = sum + ats<real_type>::abs(phi_volatile_l(ibin, iv))*kg(ibin, iv);
          }
        }, sum_kg_phi);

      if (sum_kg_phi > 0.0) {
        h_gas_i = mosaic.alpha_astem/sum_kg_phi;
        h_gas = min(h_gas, h_gas_i);
      }
    } // 6

    dtmax = min(h_gas, h_sub_max);
    // (MOSAIC prints dtmax here when it is 1.0e-10 or less)
  } // ASTEM_dtmax_soa

#endif
