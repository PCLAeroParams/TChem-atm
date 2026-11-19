#ifndef __TCHEM_IMPL_MOSAIC_ASTEM_SEMI_VOLATILES_HPP__
#define __TCHEM_IMPL_MOSAIC_ASTEM_SEMI_VOLATILES_HPP__

  //--------------------------------------------------------------------
  // part of ASTEM: integrates the gas-aerosol mass transfer of the semi-volatile
  // species (HNO3, HCl and NH3) over dtchem
  //
  // Every thread of the team calls this routine. Per-particle arrays are (nbin_a, n).
  // The thread-private scalars at the end of the argument list get the same value on every thread.
  template<typename MemberType>
  KOKKOS_INLINE_FUNCTION static
  void ASTEM_semi_volatiles(const MemberType& member,
                            const MosaicModelData<DeviceType>& mosaic,
                            const ordinal_type& nbin_a,
                            const real_type& dtchem,
                            const real_type& aH2O,
                            const real_type& Kp_nh4no3,
                            const real_type& Kp_nh4cl,
                            const real_type& Kp_nh4no3_0,
                            const real_type& Kp_nh4cl_0,
                            const real_type_1d_view_type& Keq_ll,                 // (nrxn_aer_ll)
                            const real_type_1d_view_type& Keq_gl,                 // (nrxn_aer_gl)
                            const real_type_1d_view_type& Keq_sg,                 // (nrxn_aer_sg)
                            const real_type_1d_view_type& Keq_sl,                 // (nrxn_aer_sl)
                            const real_type_2d_view_type& log_gamZ,               // (nelectrolyte, nelectrolyte)
                            const real_type_1d_view_type& MDRH_T,                 // (63)
                            const real_type_1d_view_type& gas,                    // (ngas_volatile)
                            const real_type_1d_view_type& total_species,          // (ngas_volatile)
                            real_type& tot_cl_in,
                            const real_type_1d_view_type& sum_bin_s,              // (ngas_volatile)
                            const real_type_1d_view_type& sum_vdf_s,              // (ngas_volatile)
                            const real_type_1d_view_type& sum_vol_s,              // (ngas_volatile)
                            const real_type_1d_view_type& avg_df_gas_s,           // (ngas_volatile)
                            const real_type_2d_view_type& kg,                     // (nbin_a, ngas_volatile)
                            const real_type_2d_view_type& kel,                    // (nbin_a, ngas_volatile)
                            const real_type_1d_view_type& delta_nh3_max,          // (nbin_a)
                            const real_type_1d_view_type& delta_hno3_max,         // (nbin_a)
                            const real_type_1d_view_type& delta_hcl_max,          // (nbin_a)
                            const real_type_2d_view_type& aer_solid,              // (nbin_a, naer)
                            const real_type_2d_view_type& aer_liquid,             // (nbin_a, naer)
                            const real_type_2d_view_type& aer_total,              // (nbin_a, naer)
                            const real_type_2d_view_type& electrolyte_solid,      // (nbin_a, nelectrolyte)
                            const real_type_2d_view_type& electrolyte_liquid,     // (nbin_a, nelectrolyte)
                            const real_type_2d_view_type& electrolyte_total,      // (nbin_a, nelectrolyte)
                            const real_type_1d_view_type& electrolyte_sum_solid,  // (nbin_a)
                            const real_type_1d_view_type& electrolyte_sum_liquid, // (nbin_a)
                            const real_type_1d_view_type& electrolyte_sum_total,  // (nbin_a)
                            const real_type_2d_view_type& epercent_solid,         // (nbin_a, nelectrolyte)
                            const real_type_2d_view_type& epercent_liquid,        // (nbin_a, nelectrolyte)
                            const real_type_2d_view_type& epercent_total,         // (nbin_a, nelectrolyte)
                            const real_type_2d_view_type& ma,                     // (nbin_a, nanion)
                            const real_type_2d_view_type& mc,                     // (nbin_a, ncation)
                            const real_type_2d_view_type& gam,                    // (nbin_a, nelectrolyte)
                            const real_type_2d_view_type& activity,               // (nbin_a, nelectrolyte)
                            const real_type_1d_view_type& gam_ratio,              // (nbin_a)
                            const real_type_1d_view_type& jaerosolstate,          // (nbin_a)
                            const real_type_1d_view_type& jphase,                 // (nbin_a)
                            const real_type_1d_view_type& jhyst_leg,              // (nbin_a)
                            const real_type_1d_view_type& aH2O_a,                 // (nbin_a)
                            const real_type_1d_view_type& water_a,                // (nbin_a)
                            const real_type_1d_view_type& mass_dry_a,             // (nbin_a)
                            const real_type_1d_view_type& vol_dry_a,              // (nbin_a)
                            const real_type_1d_view_type& mass_wet_a,             // (nbin_a)
                            const real_type_1d_view_type& vol_wet_a,              // (nbin_a)
                            const real_type_1d_view_type& growth_factor,          // (nbin_a)
                            const real_type_1d_view_type& MDRH,                   // (nbin_a)
                            const real_type_2d_view_type& df_gas_s,               // (nbin_a, ngas_volatile)
                            const real_type_2d_view_type& df_gas_l,               // (nbin_a, ngas_volatile)
                            const real_type_2d_view_type& phi_volatile_s,         // (nbin_a, ngas_volatile)
                            const real_type_2d_view_type& phi_volatile_l,         // (nbin_a, ngas_volatile)
                            const real_type_2d_view_type& flux_s,                 // (nbin_a, ngas_volatile)
                            const real_type_2d_view_type& flux_l,                 // (nbin_a, ngas_volatile)
                            const real_type_2d_view_type& Heff,                   // (nbin_a, ngas_volatile)
                            const real_type_2d_view_type& integrate_solid,        // (nbin_a, ngas_volatile)
                            const real_type_2d_view_type& integrate_liquid,       // (nbin_a, ngas_volatile)
                            const real_type_2d_view_type& volatile_s,             // (nbin_a, ngas_volatile)
                            const real_type_2d_view_type& h_s_i_m,                // (nbin_a, ngas_volatile)
                            const real_type_1d_view_type& idry_case3a,            // (nbin_a)
                            const real_type_2d_view_type& sfc_a,                  // work (nbin_a, ngas_volatile)
                            const real_type_2d_view_type& store,                  // work (nbin_a, naer)
                            const real_type_2d_view_type& aer_percent,            // work (nbin_a, naer)
                            const real_type_2d_view_type& na,                     // work (nbin_a, nanion)
                            const real_type_2d_view_type& nc,                     // work (nbin_a, ncation)
                            const real_type_2d_view_type& xeq_a,                  // work (nbin_a, nanion)
                            const real_type_2d_view_type& xeq_c,                  // work (nbin_a, ncation)
                            const real_type_2d_view_type& na_Ma,                  // work (nbin_a, nanion)
                            const real_type_2d_view_type& nc_Mc,                  // work (nbin_a, ncation)
                            const real_type_2d_view_type& molalities,             // work (nbin_a, nelectrolyte)
                            const real_type_2d_view_type& xmol,                   // work (nbin_a, nelectrolyte)
                            const real_type_2d_view_type& log_gam,                // work (nbin_a, nelectrolyte)
                            const real_type_2d_view_type& jsalt_present,          // work (nbin_a, nsalt)
                            const real_type_2d_view_type& hsalt,                  // work (nbin_a, nsalt)
                            const real_type_2d_view_type& phi_salt,               // work (nbin_a, nsalt)
                            const real_type_2d_view_type& phi_salt_old,           // work (nbin_a, nsalt)
                            const real_type_2d_view_type& phi_bar,                // work (nbin_a, nsalt)
                            const real_type_2d_view_type& alpha_salt,             // work (nbin_a, nsalt)
                            const real_type_2d_view_type& sat_ratio,              // work (nbin_a, nsalt)
                            const real_type_2d_view_type& flux_sl,                // work (nbin_a, nsalt)
                            const real_type_2d_view_type& frac_salt_solid,        // work (nbin_a, nsalt)
                            const real_type_2d_view_type& frac_salt_liq,          // work (nbin_a, nsalt)
                            const real_type_2d_view_type& eleliquid,              // work (nbin_a, nelectrolyte)
                            const real_type_2d_view_type& tau_p,                  // work (nbin_a, nsalt)
                            const real_type_2d_view_type& tau_d,                  // work (nbin_a, nsalt)
                            const real_type_2d_view_type& gas_bin,                // work (nbin_a, ngas_volatile)
                            const real_type_2d_view_type& total_species_bin,      // work (nbin_a, ngas_volatile)
                            const real_type_2d_view_type& hcl_net,                // work (nbin_a, 3)
                            const real_type_1d_view_type& aer_sum,                // work (nbin_a)
                            real_type& cumul_steps_ASTEM,                         // thread-private
                            real_type& isteps_ASTEM_max,                          // thread-private
                            real_type& jASTEM_fail,                               // thread-private
                            real_type& ieqblm_ASTEM,                              // thread-private
                            real_type& isteps_ASTEM) {                            // thread-private

    // local variables
    // (MOSAIC module scalars that the flux routines set before they use them)
    real_type Keq_nh4no3 = 0.0, Keq_nh4cl = 0.0, kelvin_nh4no3 = 0.0, kelvin_nh4cl = 0.0;
    real_type Keq_nh4no3_0 = 0.0, Keq_nh4cl_0 = 0.0, phi_nh4no3_s = 0.0, phi_nh4cl_s = 0.0;
    real_type dtmax = 0.0, t_new = 0.0, t_old = 0.0, t_out = 0.0;

    const auto row = [](const real_type_2d_view_type& v, const ordinal_type ibin) {
      return Kokkos::subview(v, ibin, Kokkos::ALL());
    };

    const auto flux_wet_ibin = [&](const ordinal_type ibin, real_type& ieqblm) {
      ASTEM_flux_wet(mosaic, gas, row(kg, ibin), row(kel, ibin), Keq_ll, Keq_gl, log_gamZ,
                     Kp_nh4no3, Kp_nh4cl, total_species, tot_cl_in,
                     delta_nh3_max(ibin), delta_hno3_max(ibin), delta_hcl_max(ibin),
                     row(aer_solid, ibin), row(aer_liquid, ibin), row(aer_total, ibin),
                     row(electrolyte_solid, ibin), row(electrolyte_liquid, ibin), row(electrolyte_total, ibin),
                     electrolyte_sum_liquid(ibin), electrolyte_sum_total(ibin),
                     row(epercent_liquid, ibin), aer_sum(ibin), row(aer_percent, ibin),
                     row(ma, ibin), row(mc, ibin), row(gam, ibin), row(activity, ibin), gam_ratio(ibin),
                     jaerosolstate(ibin), jphase(ibin), jhyst_leg(ibin), aH2O_a(ibin), water_a(ibin),
                     Keq_nh4no3, Keq_nh4cl, kelvin_nh4no3, kelvin_nh4cl,
                     row(sfc_a, ibin), row(df_gas_s, ibin), row(df_gas_l, ibin),
                     row(phi_volatile_s, ibin), row(phi_volatile_l, ibin), row(flux_s, ibin), row(Heff, ibin),
                     row(integrate_solid, ibin), row(integrate_liquid, ibin), ieqblm,
                     row(store, ibin), row(na, ibin), row(nc, ibin), row(xeq_a, ibin), row(xeq_c, ibin),
                     row(na_Ma, ibin), row(nc_Mc, ibin), row(molalities, ibin), row(xmol, ibin), row(log_gam, ibin));
    };

    const auto flux_mix_ibin = [&](const ordinal_type ibin, real_type& ieqblm) {
      ASTEM_flux_mix(mosaic, gas, row(kg, ibin), row(kel, ibin), Keq_ll, Keq_gl, Keq_sg, log_gamZ,
                     Kp_nh4no3, Kp_nh4cl, total_species, tot_cl_in,
                     delta_nh3_max(ibin), delta_hno3_max(ibin), delta_hcl_max(ibin),
                     row(aer_solid, ibin), row(aer_liquid, ibin), row(aer_total, ibin),
                     row(electrolyte_solid, ibin), row(electrolyte_liquid, ibin), row(electrolyte_total, ibin),
                     electrolyte_sum_liquid(ibin), electrolyte_sum_total(ibin),
                     row(epercent_solid, ibin), row(epercent_liquid, ibin), aer_sum(ibin), row(aer_percent, ibin),
                     row(ma, ibin), row(mc, ibin), row(gam, ibin), row(activity, ibin), gam_ratio(ibin),
                     jaerosolstate(ibin), jphase(ibin), jhyst_leg(ibin), aH2O_a(ibin), water_a(ibin),
                     Keq_nh4no3, Keq_nh4cl, kelvin_nh4no3, kelvin_nh4cl,
                     Keq_nh4no3_0, Keq_nh4cl_0, phi_nh4no3_s, phi_nh4cl_s,
                     row(sfc_a, ibin), row(df_gas_s, ibin), row(df_gas_l, ibin),
                     row(phi_volatile_s, ibin), row(phi_volatile_l, ibin), row(flux_s, ibin), row(Heff, ibin),
                     row(integrate_solid, ibin), row(integrate_liquid, ibin), idry_case3a(ibin), ieqblm,
                     row(store, ibin), row(na, ibin), row(nc, ibin), row(xeq_a, ibin), row(xeq_c, ibin),
                     row(na_Ma, ibin), row(nc_Mc, ibin), row(molalities, ibin), row(xmol, ibin), row(log_gam, ibin));
    };

    t_old = 0.0;
    t_out = dtchem;

    // reset counters
    isteps_ASTEM = 0.0;

    // begin integration loop over dtchem seconds
    while (true) { // label 10
      isteps_ASTEM = isteps_ASTEM + 1.0;

      ieqblm_ASTEM = 1.0; // mYES, reset to default

      // The flux routines can change gas, and MOSAIC's later particles see those
      // changes, so this loop visits the particles in order on one thread.
      Kokkos::single(Kokkos::PerTeam(member), [&](real_type& ieqblm) {
        phi_nh4no3_s = 0.0;
        phi_nh4cl_s  = 0.0;

        for (ordinal_type ibin = 0; ibin < nbin_a; ++ibin) { // do 501
          idry_case3a(ibin) = 0.0; // mNO, reset to default

          for (ordinal_type iv = 0; iv < mosaic.ngas_ioa; ++iv) {
            sfc_a(ibin, iv)            = gas(iv);
            df_gas_s(ibin, iv)         = 0.0;
            df_gas_l(ibin, iv)         = 0.0;
            flux_s(ibin, iv)           = 0.0;
            flux_l(ibin, iv)           = 0.0;
            Heff(ibin, iv)             = 0.0;
            volatile_s(ibin, iv)       = 0.0;
            phi_volatile_s(ibin, iv)   = 0.0;
            phi_volatile_l(ibin, iv)   = 0.0;
            integrate_solid(ibin, iv)  = 0.0; // mNO, reset to default
            integrate_liquid(ibin, iv) = 0.0; // mNO, reset to default
          }

          const ordinal_type state = static_cast<ordinal_type>(jaerosolstate(ibin));
          if (state == mosaic.all_solid) {
            jphase(ibin) = static_cast<real_type>(mosaic.jsolid);
            ASTEM_flux_dry(mosaic, gas, row(kg, ibin), Keq_sg, Kp_nh4no3_0, Kp_nh4cl_0,
                           row(electrolyte_solid, ibin), row(epercent_solid, ibin),
                           row(aer_solid, ibin), row(aer_liquid, ibin), row(aer_total, ibin),
                           Keq_nh4no3_0, Keq_nh4cl_0, phi_nh4no3_s, phi_nh4cl_s,
                           row(sfc_a, ibin), row(df_gas_s, ibin), row(phi_volatile_s, ibin), row(flux_s, ibin),
                           row(integrate_solid, ibin), idry_case3a(ibin), ieqblm);
          } else if (state == mosaic.all_liquid) {
            jphase(ibin) = static_cast<real_type>(mosaic.jliquid);
            flux_wet_ibin(ibin, ieqblm);
          } else if (state == mosaic.mixed) {
            flux_mix_ibin(ibin, ieqblm); // jphase(ibin) will be determined in this subr.
          }
        } // 501
      }, ieqblm_ASTEM);
      member.team_barrier();

      if (ieqblm_ASTEM == 1.0) break; // all bins have reached eqblm, so quit (goto 30)

      // calculate maximum possible internal time-step (label 11)
      ASTEM_calculate_dtmax(member, mosaic, nbin_a, dtchem, gas, kg, flux_s, df_gas_s, phi_volatile_l,
                            integrate_liquid, electrolyte_solid, aer_solid, idry_case3a, jaerosolstate,
                            volatile_s, sum_bin_s, sum_vdf_s, sum_vol_s, avg_df_gas_s, h_s_i_m, dtmax);

      t_new = t_old + dtmax; // update time
      if (t_new > t_out) {   // check if the new time step is too large
        dtmax = t_out - t_old;
        t_new = t_out*1.01;
      }

      // Integrate the gas-aerosol mass transfer of the semi-volatile species
      for (ordinal_type iv = 1; iv < 4; ++iv) { // do 20 iv = 2, 4
        real_type sum1 = 0.0, sum2 = 0.0, sum3 = 0.0, sum4a = 0.0, sum4b = 0.0;

        Kokkos::parallel_reduce(Kokkos::TeamVectorRange(member, nbin_a),
          [&](const ordinal_type ibin, real_type& sum) {
            if (static_cast<ordinal_type>(jaerosolstate(ibin)) == mosaic.no_aerosol) return; // goto 21
            sum = sum + aer_liquid(ibin, iv)/
                  (1.0 + dtmax*kg(ibin, iv)*Heff(ibin, iv)*integrate_liquid(ibin, iv));
          }, sum1);
        Kokkos::parallel_reduce(Kokkos::TeamVectorRange(member, nbin_a),
          [&](const ordinal_type ibin, real_type& sum) {
            if (static_cast<ordinal_type>(jaerosolstate(ibin)) == mosaic.no_aerosol) return; // goto 21
            sum = sum + kg(ibin, iv)*integrate_liquid(ibin, iv)/
                  (1.0 + dtmax*kg(ibin, iv)*Heff(ibin, iv)*integrate_liquid(ibin, iv));
          }, sum2);
        Kokkos::parallel_reduce(Kokkos::TeamVectorRange(member, nbin_a),
          [&](const ordinal_type ibin, real_type& sum) {
            if (static_cast<ordinal_type>(jaerosolstate(ibin)) == mosaic.no_aerosol) return; // goto 21
            sum = sum + aer_solid(ibin, iv);
          }, sum3);
        Kokkos::parallel_reduce(Kokkos::TeamVectorRange(member, nbin_a),
          [&](const ordinal_type ibin, real_type& sum) {
            if (static_cast<ordinal_type>(jaerosolstate(ibin)) == mosaic.no_aerosol) return; // goto 21
            if (flux_s(ibin, iv) > 0.0) sum = sum + dtmax*flux_s(ibin, iv);
          }, sum4a);
        Kokkos::parallel_reduce(Kokkos::TeamVectorRange(member, nbin_a),
          [&](const ordinal_type ibin, real_type& sum) {
            if (static_cast<ordinal_type>(jaerosolstate(ibin)) == mosaic.no_aerosol) return; // goto 21
            if (flux_s(ibin, iv) < 0.0) sum = sum + min(h_s_i_m(ibin, iv), dtmax)*flux_s(ibin, iv);
          }, sum4b);

        // the solid-phase update of do 21, after the sums have read the old aer_solid
        Kokkos::parallel_for(Kokkos::TeamVectorRange(member, nbin_a), [&](const ordinal_type ibin) {
          if (static_cast<ordinal_type>(jaerosolstate(ibin)) == mosaic.no_aerosol) return; // goto 21
          if (flux_s(ibin, iv) > 0.0) {
            aer_solid(ibin, iv) = aer_solid(ibin, iv) + dtmax*flux_s(ibin, iv);
          } else if (flux_s(ibin, iv) < 0.0) {
            aer_solid(ibin, iv) = aer_solid(ibin, iv) + min(h_s_i_m(ibin, iv), dtmax)*flux_s(ibin, iv);
            aer_solid(ibin, iv) = max(aer_solid(ibin, iv), 0.0);
          }
        });
        member.team_barrier();

        // first update gas concentration
        Kokkos::single(Kokkos::PerTeam(member), [&]() {
          const real_type sum4 = sum4a + sum4b;
          gas(iv) = (total_species(iv) - (sum1 + sum3 + sum4))/(1.0 + dtmax*sum2);
          gas(iv) = max(gas(iv), 0.0);
        });
        member.team_barrier();

        // now update aer concentration in the liquid phase
        Kokkos::parallel_for(Kokkos::TeamVectorRange(member, nbin_a), [&](const ordinal_type ibin) { // do 22
          if (integrate_liquid(ibin, iv) == 1.0) { // mYES
            aer_liquid(ibin, iv) = (aer_liquid(ibin, iv) + dtmax*kg(ibin, iv)*gas(iv))/
                                   (1.0 + dtmax*kg(ibin, iv)*Heff(ibin, iv));
          }
        });
        member.team_barrier();
      } // 20

      // Update the particles. form_electrolytes only adds to gas, so each particle
      // works on its own copy of gas and the changes are applied in order below.
      Kokkos::parallel_for(Kokkos::TeamVectorRange(member, nbin_a), [&](const ordinal_type ibin) { // do 40
        for (ordinal_type k = 0; k < 3; ++k) hcl_net(ibin, k) = 0.0;
        if (static_cast<ordinal_type>(jaerosolstate(ibin)) == mosaic.no_aerosol) return; // goto 40

        const auto gas_p = row(gas_bin, ibin);
        const auto total_species_p = row(total_species_bin, ibin);
        for (ordinal_type iv = 0; iv < mosaic.ngas_volatile; ++iv) {
          gas_p(iv) = gas(iv);
          total_species_p(iv) = total_species(iv);
        }
        real_type tot_cl_p = tot_cl_in;

        // k-th form_electrolytes call; hcl_net holds what it added to gas(ihcl_g)
        // before form_nacl's limit of gas(ihcl_g) at zero
        const auto form_electrolytes_k = [&](const ordinal_type k, const ordinal_type jp,
                                             const real_type_2d_view_type& aer_curr,
                                             const real_type_2d_view_type& electrolyte,
                                             real_type& electrolyte_sum,
                                             const real_type_2d_view_type& epercent) {
          const real_type hcl_old = gas_p(mosaic.ihcl_g);
          const real_type tot_cl_old = tot_cl_p;
          form_electrolytes(mosaic, static_cast<real_type>(jp), row(aer_curr, ibin),
                            row(aer_solid, ibin), row(aer_liquid, ibin), row(aer_total, ibin),
                            gas_p, row(store, ibin), row(electrolyte, ibin), total_species_p, tot_cl_p,
                            electrolyte_sum, row(epercent, ibin), aer_sum(ibin), row(aer_percent, ibin));
          hcl_net(ibin, k) = gas_p(mosaic.ihcl_g) - hcl_old - (tot_cl_p - tot_cl_old);
        };

        const ordinal_type phase = static_cast<ordinal_type>(jphase(ibin));
        if (phase == mosaic.jsolid) {
          form_electrolytes_k(0, mosaic.jsolid, aer_solid, electrolyte_solid, electrolyte_sum_solid(ibin), epercent_solid); // degas excess nh3 (if present)
        } else if (phase == mosaic.jliquid) {
          form_electrolytes_k(0, mosaic.jliquid, aer_liquid, electrolyte_liquid, electrolyte_sum_liquid(ibin), epercent_liquid); // degas excess nh3 (if present)
        } else if (phase == mosaic.jtotal) {
          form_electrolytes_k(0, mosaic.jsolid, aer_solid, electrolyte_solid, electrolyte_sum_solid(ibin), epercent_solid); // degas excess nh3 (if present)
          form_electrolytes_k(1, mosaic.jliquid, aer_liquid, electrolyte_liquid, electrolyte_sum_liquid(ibin), epercent_liquid); // degas excess nh3 (if present)
        }

        // now update jtotal
        for (ordinal_type iv = 1; iv < mosaic.ngas_ioa; ++iv) {
          aer_total(ibin, iv) = aer_solid(ibin, iv) + aer_liquid(ibin, iv);
        }

        form_electrolytes_k(2, mosaic.jtotal, aer_total, electrolyte_total, electrolyte_sum_total(ibin), epercent_total); // for MDRH diagnosis

        // update internal phase equilibrium
        if (static_cast<ordinal_type>(jhyst_leg(ibin)) == mosaic.jhyst_lo) {
          ASTEM_update_phase_eqblm(mosaic, row(aer_liquid, ibin), row(aer_solid, ibin), row(aer_total, ibin),
                                   row(electrolyte_solid, ibin), row(electrolyte_liquid, ibin), row(electrolyte_total, ibin),
                                   row(epercent_liquid, ibin), row(epercent_solid, ibin), row(epercent_total, ibin),
                                   Keq_sl, Keq_ll, log_gamZ, MDRH_T, aH2O,
                                   jaerosolstate(ibin), jphase(ibin), jhyst_leg(ibin), aH2O_a(ibin), water_a(ibin),
                                   mass_dry_a(ibin), vol_dry_a(ibin), mass_wet_a(ibin), vol_wet_a(ibin), growth_factor(ibin),
                                   electrolyte_sum_solid(ibin), MDRH(ibin),
                                   row(jsalt_present, ibin), row(hsalt, ibin), row(phi_salt, ibin), row(phi_salt_old, ibin),
                                   row(phi_bar, ibin), row(alpha_salt, ibin), row(sat_ratio, ibin), row(flux_sl, ibin),
                                   row(frac_salt_solid, ibin), row(frac_salt_liq, ibin), row(eleliquid, ibin),
                                   row(store, ibin), row(na, ibin), row(nc, ibin), row(xeq_a, ibin), row(xeq_c, ibin),
                                   row(na_Ma, ibin), row(nc_Mc, ibin), row(aer_percent, ibin), row(molalities, ibin), row(xmol, ibin),
                                   row(ma, ibin), row(mc, ibin), row(log_gam, ibin), row(gam, ibin), row(activity, ibin),
                                   row(tau_p, ibin), row(tau_d, ibin));
        } else {
          do_full_deliquescence(mosaic, row(electrolyte_solid, ibin), row(electrolyte_liquid, ibin), row(electrolyte_total, ibin),
                                row(aer_solid, ibin), row(aer_liquid, ibin), row(aer_total, ibin)); // simply do liquid <-- total
        }
      }); // 40
      member.team_barrier();

      // apply the particles' changes to gas in MOSAIC's particle order
      Kokkos::single(Kokkos::PerTeam(member), [&]() {
        const real_type hno3_old = gas(mosaic.ihno3_g);
        const real_type nh3_old  = gas(mosaic.inh3_g);
        for (ordinal_type ibin = 0; ibin < nbin_a; ++ibin) {
          if (static_cast<ordinal_type>(jaerosolstate(ibin)) == mosaic.no_aerosol) continue;
          gas(mosaic.ihno3_g) = gas(mosaic.ihno3_g) + gas_bin(ibin, mosaic.ihno3_g) - hno3_old;
          gas(mosaic.inh3_g)  = gas(mosaic.inh3_g)  + gas_bin(ibin, mosaic.inh3_g)  - nh3_old;
          for (ordinal_type k = 0; k < 3; ++k) {
            gas(mosaic.ihcl_g) = gas(mosaic.ihcl_g) + hcl_net(ibin, k);
            if (gas(mosaic.ihcl_g) < 0.0) { // as in form_nacl
              total_species(mosaic.ihcl_g) = total_species(mosaic.ihcl_g) - gas(mosaic.ihcl_g);
              tot_cl_in = tot_cl_in - gas(mosaic.ihcl_g);
              gas(mosaic.ihcl_g) = 0.0;
            }
          }
        }
      });
      member.team_barrier();

      t_old = t_new;

      if (isteps_ASTEM >= static_cast<real_type>(mosaic.nmax_astem)) {
        jASTEM_fail = jASTEM_fail + 1.0;
        break; // goto 30
      } else if (t_new < t_out) {
        continue; // goto 10
      }

      // check if end of dtchem reached
      if (t_new < 0.9999*t_out) continue; // goto 10

      break;
    } // end of overall integration loop over dtchem seconds (label 30)

    cumul_steps_ASTEM = cumul_steps_ASTEM + isteps_ASTEM;
    isteps_ASTEM_max  = max(isteps_ASTEM_max, isteps_ASTEM);

    // call subs to calculate fluxes over mixed-phase particles to update H+ ions,
    // which were wiped off during update_phase_eqblm
    Kokkos::single(Kokkos::PerTeam(member), [&](real_type& ieqblm) {
      for (ordinal_type ibin = 0; ibin < nbin_a; ++ibin) {
        if (static_cast<ordinal_type>(jaerosolstate(ibin)) == mosaic.mixed) {
          if (electrolyte_solid(ibin, mosaic.jnh4no3) > 0.0 ||
              electrolyte_solid(ibin, mosaic.jnh4cl)  > 0.0) {
            flux_mix_ibin(ibin, ieqblm); // jphase(ibin) will be determined in this subr.
          } else {
            jphase(ibin) = static_cast<real_type>(mosaic.jliquid);
            flux_wet_ibin(ibin, ieqblm);
          }
        }
      }
    }, ieqblm_ASTEM);
    member.team_barrier();
  } // ASTEM_semi_volatiles

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
