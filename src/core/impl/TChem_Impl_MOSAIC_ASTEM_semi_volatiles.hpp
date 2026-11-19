#ifndef __TCHEM_IMPL_MOSAIC_ASTEM_SEMI_VOLATILES_HPP__
#define __TCHEM_IMPL_MOSAIC_ASTEM_SEMI_VOLATILES_HPP__

  //--------------------------------------------------------------------
  // part of ASTEM: computes max time step for gas-aerosol integration
  //
  // Every thread of the team calls this routine. The loops over particles
  // (ibin) are spread over the team, and the sums over particles are team
  // reductions, whose results every thread receives.
  template<typename MemberType>
  KOKKOS_INLINE_FUNCTION static
  void ASTEM_calculate_dtmax(const MemberType& member,
                             const MosaicModelData<DeviceType>& mosaic,
                             const ordinal_type& nbin_a,
                             const real_type& dtchem,
                             const real_type_1d_view_type& gas,                // (ngas_volatile)
                             const real_type_2d_view_type& kg,                 // (nbin_a, ngas_volatile)
                             const real_type_2d_view_type& flux_s,             // (nbin_a, ngas_volatile)
                             const real_type_2d_view_type& df_gas_s,           // (nbin_a, ngas_volatile)
                             const real_type_2d_view_type& phi_volatile_l,     // (nbin_a, ngas_volatile)
                             const real_type_2d_view_type& integrate_liquid,   // (nbin_a, ngas_volatile)
                             const real_type_2d_view_type& electrolyte_solid,  // (nbin_a, nelectrolyte)
                             const real_type_2d_view_type& aer_solid,          // (nbin_a, naer)
                             const real_type_1d_view_type& idry_case3a,        // (nbin_a)
                             const real_type_1d_view_type& jaerosolstate,      // (nbin_a)
                             const real_type_2d_view_type& volatile_s,         // (nbin_a, ngas_volatile)
                             const real_type_1d_view_type& sum_bin_s,          // (ngas_volatile)
                             const real_type_1d_view_type& sum_vdf_s,          // (ngas_volatile)
                             const real_type_1d_view_type& sum_vol_s,          // (ngas_volatile)
                             const real_type_1d_view_type& avg_df_gas_s,       // (ngas_volatile)
                             const real_type_2d_view_type& h_s_i_m,            // (nbin_a, ngas_volatile)
                             real_type& dtmax) {

    // local variables
    real_type h_gas, h_sub_max, h_gas_i, h_gas_l, h_gas_s, sum_kg_phi, sumflux_s;

    h_sub_max = 60.0; // sec

    // GAS-SIDE

    // solid-phase
    // calculate h_gas_i and h_gas_l

    h_gas_s = 2.e16;

    for (ordinal_type iv = 1; iv < mosaic.ngas_ioa; ++iv) {
      h_gas_i = 1.e16;
      sumflux_s = 0.0;
      Kokkos::parallel_reduce(Kokkos::TeamVectorRange(member, nbin_a),
        [&](const ordinal_type ibin, real_type& sum) {
          if (flux_s(ibin, iv) > 0.0) {
            sum = sum + flux_s(ibin, iv);
          }
        }, sumflux_s);

      if (sumflux_s > 0.0) {
        h_gas_i = 0.1*gas(iv)/sumflux_s;
        h_gas_s = min(h_gas_s, h_gas_i);
      }
    }

    // liquid-phase
    // calculate h_gas_s and h_gas_l

    h_gas_l = 2.e16;

    for (ordinal_type iv = 1; iv < mosaic.ngas_ioa; ++iv) {
      h_gas_i = 1.e16;
      sum_kg_phi = 0.0;
      Kokkos::parallel_reduce(Kokkos::TeamVectorRange(member, nbin_a),
        [&](const ordinal_type ibin, real_type& sum) {
          if (integrate_liquid(ibin, iv) == 1.0) { // mYES
            sum = sum + ats<real_type>::abs(phi_volatile_l(ibin, iv))*kg(ibin, iv);
          }
        }, sum_kg_phi);

      if (sum_kg_phi > 0.0) {
        h_gas_i = mosaic.alpha_astem/sum_kg_phi;
        h_gas_l = min(h_gas_l, h_gas_i);
      }
    }

    h_gas = min(h_gas_s, h_gas_l);
    h_gas = min(h_gas, h_sub_max);

    // AEROSOL-SIDE: solid-phase

    // first load volatile_solid array
    Kokkos::parallel_for(Kokkos::TeamVectorRange(member, nbin_a), [&](const ordinal_type ibin) {
      volatile_s(ibin, mosaic.ino3_a) = electrolyte_solid(ibin, mosaic.jnh4no3);
      volatile_s(ibin, mosaic.inh4_a) = electrolyte_solid(ibin, mosaic.jnh4cl) +
                                        electrolyte_solid(ibin, mosaic.jnh4no3);

      if (idry_case3a(ibin) == 1.0) { // mYES
        volatile_s(ibin, mosaic.icl_a) = aer_solid(ibin, mosaic.icl_a);
      } else {
        volatile_s(ibin, mosaic.icl_a) = electrolyte_solid(ibin, mosaic.jnh4cl);
      }
    });
    member.team_barrier();

    // next calculate weighted avg_df_gas_s
    for (ordinal_type iv = 1; iv < mosaic.ngas_ioa; ++iv) {
      // one team reduction per sum
      real_type sum_bin = 0.0, sum_vdf = 0.0, sum_vol = 0.0;
      Kokkos::parallel_reduce(Kokkos::TeamVectorRange(member, nbin_a),
        [&](const ordinal_type ibin, real_type& sum) {
          if (flux_s(ibin, iv) < 0.0) { // aer -> gas
            sum = sum + 1.0;
          }
        }, sum_bin);
      Kokkos::parallel_reduce(Kokkos::TeamVectorRange(member, nbin_a),
        [&](const ordinal_type ibin, real_type& sum) {
          if (flux_s(ibin, iv) < 0.0) { // aer -> gas
            sum = sum + volatile_s(ibin, iv)*df_gas_s(ibin, iv);
          }
        }, sum_vdf);
      Kokkos::parallel_reduce(Kokkos::TeamVectorRange(member, nbin_a),
        [&](const ordinal_type ibin, real_type& sum) {
          if (flux_s(ibin, iv) < 0.0) { // aer -> gas
            sum = sum + volatile_s(ibin, iv);
          }
        }, sum_vol);

      Kokkos::single(Kokkos::PerTeam(member), [&]() {
        sum_bin_s(iv) = sum_bin;
        sum_vdf_s(iv) = sum_vdf;
        sum_vol_s(iv) = sum_vol;

        if (sum_vol_s(iv) > 0.0) {
          avg_df_gas_s(iv) = sum_vdf_s(iv)/sum_vol_s(iv);
        } else {
          avg_df_gas_s(iv) = 1.0; // never used, but set to 1.0 just to be safe
        }
      });
    }
    member.team_barrier();

    // calculate h_s_i_m
    Kokkos::parallel_for(Kokkos::TeamVectorRange(member, nbin_a), [&](const ordinal_type ibin) {
      if (static_cast<ordinal_type>(jaerosolstate(ibin)) == mosaic.no_aerosol) return; // goto 20

      for (ordinal_type iv = 1; iv < mosaic.ngas_ioa; ++iv) {
        if (flux_s(ibin, iv) < 0.0) { // aer -> gas
          real_type alpha = ats<real_type>::abs(avg_df_gas_s(iv))/
                            (volatile_s(ibin, iv)*sum_bin_s(iv));
          alpha = min(alpha, 1.0);

          if (idry_case3a(ibin) == 1.0) alpha = 1.0; // mYES

          h_s_i_m(ibin, iv) = -alpha*volatile_s(ibin, iv)/flux_s(ibin, iv);
        }
      }
    });
    member.team_barrier();

    dtmax = min(dtchem, h_gas);
    // (MOSAIC prints dtmax here when it is zero)
  } // ASTEM_calculate_dtmax

  //--------------------------------------------------------------------
  // part of ASTEM: updates solid-liquid partitioning after each gas-aerosol
  // mass transfer step
  KOKKOS_INLINE_FUNCTION static
  void ASTEM_update_phase_eqblm(const MosaicModelData<DeviceType>& mosaic,
                                const real_type_1d_view_type& aer_liquid,          // (naer)
                                const real_type_1d_view_type& aer_solid,           // (naer)
                                const real_type_1d_view_type& aer_total,           // (naer)
                                const real_type_1d_view_type& electrolyte_solid,   // (nelectrolyte)
                                const real_type_1d_view_type& electrolyte_liquid,  // (nelectrolyte)
                                const real_type_1d_view_type& electrolyte_total,   // (nelectrolyte)
                                const real_type_1d_view_type& epercent_liquid,     // (nelectrolyte)
                                const real_type_1d_view_type& epercent_solid,      // (nelectrolyte)
                                const real_type_1d_view_type& epercent_total,      // (nelectrolyte)
                                const real_type_1d_view_type& Keq_sl,              // (nrxn_aer_sl)
                                const real_type_1d_view_type& Keq_ll,              // (nrxn_aer_ll)
                                const real_type_2d_view_type& log_gamZ,            // (nelectrolyte, nelectrolyte)
                                const real_type_1d_view_type& MDRH_T,              // (63)
                                const real_type& aH2O,
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
                                real_type& electrolyte_sum_solid,
                                real_type& MDRH,
                                const real_type_1d_view_type& jsalt_present,       // (nsalt)
                                const real_type_1d_view_type& hsalt,               // (nsalt)
                                const real_type_1d_view_type& phi_salt,            // (nsalt)
                                const real_type_1d_view_type& phi_salt_old,        // (nsalt)
                                const real_type_1d_view_type& phi_bar,             // (nsalt)
                                const real_type_1d_view_type& alpha_salt,          // (nsalt)
                                const real_type_1d_view_type& sat_ratio,           // (nsalt)
                                const real_type_1d_view_type& flux_sl,             // (nsalt)
                                const real_type_1d_view_type& frac_salt_solid,     // (nsalt)
                                const real_type_1d_view_type& frac_salt_liq,       // (nsalt)
                                const real_type_1d_view_type& eleliquid,           // (nelectrolyte)
                                const real_type_1d_view_type& store,               // (naer)
                                const real_type_1d_view_type& na,                  // (nanion)
                                const real_type_1d_view_type& nc,                  // (ncation)
                                const real_type_1d_view_type& xeq_a,               // (nanion)
                                const real_type_1d_view_type& xeq_c,               // (ncation)
                                const real_type_1d_view_type& na_Ma,               // (nanion)
                                const real_type_1d_view_type& nc_Mc,               // (ncation)
                                const real_type_1d_view_type& aer_percent,         // (naer)
                                const real_type_1d_view_type& molalities,          // (nelectrolyte)
                                const real_type_1d_view_type& xmol,                // (nelectrolyte)
                                const real_type_1d_view_type& ma,                  // (nanion)
                                const real_type_1d_view_type& mc_view,             // (ncation)
                                const real_type_1d_view_type& log_gam,             // (nelectrolyte)
                                const real_type_1d_view_type& gam,                 // (nelectrolyte)
                                const real_type_1d_view_type& activity,            // (nelectrolyte)
                                const real_type_1d_view_type& tau_p,               // (nsalt)
                                const real_type_1d_view_type& tau_d) {             // (nsalt)

    auto jsalt_index_view = mosaic.jsalt_index.template view<DeviceType>();
    auto jsulf_poor_view  = mosaic.jsulf_poor.template view<DeviceType>();
    auto jsulf_rich_view  = mosaic.jsulf_rich.template view<DeviceType>();

    ordinal_type jdum, j_index;
    real_type XT = 0.0;

    // calculate overall sulfate ratio
    calculate_XT(mosaic, aer_total, XT); // calc updated XT

    // now diagnose MDRH
    if (!(XT < 1.0 && XT >= 0.0)) { // excess sulfate domain - no MDRH exists (goto 10)

      jdum = 0;
      for (ordinal_type js = 0; js < mosaic.nsalt; ++js) {
        jsalt_present(js) = 0.0; // default value - salt absent
        if (epercent_total(js) > mosaic.ptol_mol_astem) {
          jsalt_present(js) = 1.0; // salt present
          jdum = jdum + jsalt_index_view(js);
        }
      }

      if (jdum == 0) {
        jaerosolstate = static_cast<real_type>(mosaic.all_solid); // no significant soluble material present
        jphase        = static_cast<real_type>(mosaic.jsolid);
        adjust_solid_aerosol(mosaic, aer_solid, aer_liquid, aer_total,
                             electrolyte_solid, electrolyte_liquid, electrolyte_total,
                             epercent_solid, epercent_liquid, epercent_total,
                             water_a, jphase, jhyst_leg);
        return;
      }

      if (XT >= 2.0 || XT < 0.0) {
        j_index = jsulf_poor_view(jdum);
      } else {
        j_index = jsulf_rich_view(jdum);
      }

      MDRH = MDRH_T(j_index - 1);

      if (aH2O*100.0 < MDRH) {
        jaerosolstate = static_cast<real_type>(mosaic.all_solid);
        jphase        = static_cast<real_type>(mosaic.jsolid);
        adjust_solid_aerosol(mosaic, aer_solid, aer_liquid, aer_total,
                             electrolyte_solid, electrolyte_liquid, electrolyte_total,
                             epercent_solid, epercent_liquid, epercent_total,
                             water_a, jphase, jhyst_leg);
        return;
      }

    }

    // none of the above means it must be sub-saturated or mixed-phase (label 10)
    if (static_cast<ordinal_type>(jphase) == mosaic.jsolid) {
      do_full_deliquescence(mosaic,
                            electrolyte_solid, electrolyte_liquid, electrolyte_total,
                            aer_solid, aer_liquid, aer_total);
      MESA_PTC(mosaic,
               aer_liquid, aer_solid, aer_total,
               electrolyte_solid, electrolyte_liquid, electrolyte_total,
               epercent_liquid, epercent_solid, epercent_total,
               Keq_sl, Keq_ll, log_gamZ,
               jaerosolstate, jphase, jhyst_leg, aH2O_a, water_a,
               mass_dry_a, vol_dry_a, mass_wet_a, vol_wet_a, growth_factor,
               electrolyte_sum_solid,
               jsalt_present, hsalt, phi_salt, phi_salt_old, phi_bar,
               alpha_salt, sat_ratio, flux_sl, frac_salt_solid, frac_salt_liq,
               eleliquid,
               store, na, nc, xeq_a, xeq_c, na_Ma, nc_Mc, aer_percent,
               molalities, xmol, ma, mc_view, log_gam, gam, activity,
               tau_p, tau_d);
    } else {
      MESA_PTC(mosaic,
               aer_liquid, aer_solid, aer_total,
               electrolyte_solid, electrolyte_liquid, electrolyte_total,
               epercent_liquid, epercent_solid, epercent_total,
               Keq_sl, Keq_ll, log_gamZ,
               jaerosolstate, jphase, jhyst_leg, aH2O_a, water_a,
               mass_dry_a, vol_dry_a, mass_wet_a, vol_wet_a, growth_factor,
               electrolyte_sum_solid,
               jsalt_present, hsalt, phi_salt, phi_salt_old, phi_bar,
               alpha_salt, sat_ratio, flux_sl, frac_salt_solid, frac_salt_liq,
               eleliquid,
               store, na, nc, xeq_a, xeq_c, na_Ma, nc_Mc, aer_percent,
               molalities, xmol, ma, mc_view, log_gam, gam, activity,
               tau_p, tau_d);
    }

    return;
  } // ASTEM_update_phase_eqblm

#endif
