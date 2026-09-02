#ifndef __TCHEM_IMPL_MOSAIC_ACTIVITIES_HPP__
#define __TCHEM_IMPL_MOSAIC_ACTIVITIES_HPP__

KOKKOS_INLINE_FUNCTION static
  void compute_activities(const MosaicModelData<DeviceType>& mosaic,
                          const real_type_1d_view_type& molalities,
                          const real_type_1d_view_type& xmol,
                          const real_type_1d_view_type& aer_liquid,
                          const real_type_1d_view_type& ma,
                          const real_type_1d_view_type& mc,
                          const real_type_1d_view_type& Keq_ll,
                          const real_type_1d_view_type& electrolyte_solid,
                          const real_type_1d_view_type& electrolyte_liquid,
                          const real_type_1d_view_type& electrolyte_total,
                          const real_type_1d_view_type& log_gam,
                          const real_type_2d_view_type& log_gamZ,
                          const real_type_1d_view_type& gam,
                          const real_type_1d_view_type& activity,
                          real_type& jaerosolstate,
                          real_type& jphase,
                          real_type& jhyst_leg,
                          real_type& aH2O_a,
                          real_type& water_a) {

    // local variables
    real_type a_c, sum_elec, gam_ratio, mSULF = 0.0;
    ordinal_type jA = 0;

    // get aerosol water activity
    aerosol_water(mosaic, electrolyte_liquid, aH2O_a, molalities, jaerosolstate, jphase, jhyst_leg, water_a);

    if (water_a == 0.0) {
      return;
    }

    // get sulfate ratio to determine regime
    real_type XT = 0.0;
    calculate_XT(mosaic, aer_liquid, XT);

    if (XT > 2.0 || XT < 0.0) {
      // SULFATE POOR: fully dissociated electrolytes

      // anion molalities (mol / kg water)
      ma(mosaic.ja_so4)  = 1.e-9 * aer_liquid(mosaic.iso4_a) / water_a;
      ma(mosaic.ja_hso4) = 0.0;
      ma(mosaic.ja_no3)  = 1.e-9 * aer_liquid(mosaic.ino3_a) / water_a;
      ma(mosaic.ja_cl)   = 1.e-9 * aer_liquid(mosaic.icl_a)  / water_a;
      ma(mosaic.ja_msa)  = 1.e-9 * aer_liquid(mosaic.imsa_a) / water_a;

      // cation molalities (mol / kg water)
      mc(mosaic.jc_ca)  = 1.e-9 * aer_liquid(mosaic.ica_a)  / water_a;
      mc(mosaic.jc_nh4) = 1.e-9 * aer_liquid(mosaic.inh4_a) / water_a;
      mc(mosaic.jc_na)  = 1.e-9 * aer_liquid(mosaic.ina_a)  / water_a;
      a_c               = (
                          (2. * ma(mosaic.ja_so4)  +
                                ma(mosaic.ja_no3)  +
                                ma(mosaic.ja_cl)   +
                                ma(mosaic.ja_msa)) -
                          (2. * mc(mosaic.jc_ca)   +
                                mc(mosaic.jc_nh4)  +
                                mc(mosaic.jc_na))  );
      mc(mosaic.jc_h)   = 0.5 * ( (a_c) +
                                  (ats<real_type>::sqrt(a_c*a_c + 4. * Keq_ll(2))) );

      if (mc(mosaic.jc_h) == 0.0) {
        mc(mosaic.jc_h) = 1.e-10;
      }

      sum_elec = 2. * electrolyte_liquid(mosaic.jnh4no3) +
                 2. * electrolyte_liquid(mosaic.jnh4cl)  +
                 3. * electrolyte_liquid(mosaic.jnh4so4) +
                 3. * electrolyte_liquid(mosaic.jna2so4) +
                 2. * electrolyte_liquid(mosaic.jnano3)  +
                 2. * electrolyte_liquid(mosaic.jnacl)   +
                 3. * electrolyte_liquid(mosaic.jcano3)  +
                 3. * electrolyte_liquid(mosaic.jcacl2)  +
                 2. * electrolyte_liquid(mosaic.jhno3)   +
                 2. * electrolyte_liquid(mosaic.jhcl);

      if (sum_elec == 0.0) {
        for (ordinal_type jA = 0; jA < mosaic.nelectrolyte; jA++) {
          gam(jA) = 1.0;
        }
        // FIXME: duplicated code to avoid goto statement
        gam(mosaic.jlvcite)  = 1.0;
        gam(mosaic.jnh4hso4) = 1.0;
        gam(mosaic.jnh4msa)  = 1.0;
        gam(mosaic.jna3hso4) = 1.0;
        gam(mosaic.jnamsa)   = 1.0;
        gam(mosaic.jcamsa2)  = 1.0;

        activity(mosaic.jlvcite)  = 0.0;
        activity(mosaic.jnh4hso4) = 0.0;

        activity(mosaic.jnh4msa) = mc(mosaic.jc_nh4) * ma(mosaic.ja_msa) *
                                   gam(mosaic.jnh4msa) * gam(mosaic.jnh4msa);

        activity(mosaic.jna3hso4) = 0.0;
        activity(mosaic.jnahso4)  = 0.0;

        activity(mosaic.jnamsa)  = mc(mosaic.jc_na) * ma(mosaic.ja_msa) *
                                   gam(mosaic.jnamsa) * gam(mosaic.jnamsa);
        activity(mosaic.jcamsa2) = mc(mosaic.jc_ca) * ma(mosaic.ja_msa) *
                                   gam(mosaic.jcamsa2) * gam(mosaic.jcamsa2) * gam(mosaic.jcamsa2);

        gam_ratio = gam(mosaic.jnh4no3) * gam(mosaic.jnh4no3) /
                    gam(mosaic.jhno3)   * gam(mosaic.jhno3);

        return;
      }

      // ionic mole fractions
      xmol(mosaic.jnh4no3) = 2. * electrolyte_liquid(mosaic.jnh4no3) / sum_elec;
      xmol(mosaic.jnh4cl)  = 2. * electrolyte_liquid(mosaic.jnh4cl)  / sum_elec;
      xmol(mosaic.jnh4so4) = 3. * electrolyte_liquid(mosaic.jnh4so4) / sum_elec;
      xmol(mosaic.jna2so4) = 3. * electrolyte_liquid(mosaic.jna2so4) / sum_elec;
      xmol(mosaic.jnano3)  = 2. * electrolyte_liquid(mosaic.jnano3)  / sum_elec;
      xmol(mosaic.jnacl)   = 2. * electrolyte_liquid(mosaic.jnacl)   / sum_elec;
      xmol(mosaic.jcano3)  = 3. * electrolyte_liquid(mosaic.jcano3)  / sum_elec;
      xmol(mosaic.jcacl2)  = 3. * electrolyte_liquid(mosaic.jcacl2)  / sum_elec;
      xmol(mosaic.jhno3)   = 2. * electrolyte_liquid(mosaic.jhno3)   / sum_elec;
      xmol(mosaic.jhcl)    = 2. * electrolyte_liquid(mosaic.jhcl)    / sum_elec;

      jA = mosaic.jnh4so4;
      if (xmol(jA) > 0.0) {
        log_gam(jA) = xmol(mosaic.jnh4no3) * log_gamZ(jA,mosaic.jnh4no3) +
                      xmol(mosaic.jnh4cl)  * log_gamZ(jA,mosaic.jnh4cl)  +
                      xmol(mosaic.jnh4so4) * log_gamZ(jA,mosaic.jnh4so4) +
                      xmol(mosaic.jna2so4) * log_gamZ(jA,mosaic.jna2so4) +
                      xmol(mosaic.jnano3)  * log_gamZ(jA,mosaic.jnano3)  +
                      xmol(mosaic.jnacl)   * log_gamZ(jA,mosaic.jnacl)   +
                      xmol(mosaic.jcano3)  * log_gamZ(jA,mosaic.jcano3)  +
                      xmol(mosaic.jcacl2)  * log_gamZ(jA,mosaic.jcacl2)  +
                      xmol(mosaic.jhno3)   * log_gamZ(jA,mosaic.jhno3)   +
                      xmol(mosaic.jhcl)    * log_gamZ(jA,mosaic.jhcl);
        gam(jA) = ats<real_type>::pow(10., log_gam(jA));
        activity(jA) = ats<real_type>::pow(mc(mosaic.jc_nh4),2.) *
                      ma(mosaic.ja_so4) *
                      ats<real_type>::pow(gam(jA),3.);
      }

      jA = mosaic.jnh4no3;
      if (xmol(jA) > 0.0) {
        log_gam(jA) = xmol(mosaic.jnh4no3) * log_gamZ(jA,mosaic.jnh4no3) +
                      xmol(mosaic.jnh4cl)  * log_gamZ(jA,mosaic.jnh4cl)  +
                      xmol(mosaic.jnh4so4) * log_gamZ(jA,mosaic.jnh4so4) +
                      xmol(mosaic.jna2so4) * log_gamZ(jA,mosaic.jna2so4) +
                      xmol(mosaic.jnano3)  * log_gamZ(jA,mosaic.jnano3)  +
                      xmol(mosaic.jnacl)   * log_gamZ(jA,mosaic.jnacl)   +
                      xmol(mosaic.jcano3)  * log_gamZ(jA,mosaic.jcano3)  +
                      xmol(mosaic.jcacl2)  * log_gamZ(jA,mosaic.jcacl2)  +
                      xmol(mosaic.jhno3)   * log_gamZ(jA,mosaic.jhno3)   +
                      xmol(mosaic.jhcl)    * log_gamZ(jA,mosaic.jhcl);
        gam(jA) = ats<real_type>::pow(10., log_gam(jA));
        activity(jA) = mc(mosaic.jc_nh4) *
                      ma(mosaic.ja_no3) *
                      ats<real_type>::pow(gam(jA),2.);
      }

      jA = mosaic.jnh4cl;
      if (xmol(jA) > 0.0) {
        log_gam(jA) = xmol(mosaic.jnh4no3) * log_gamZ(jA,mosaic.jnh4no3) +
                      xmol(mosaic.jnh4cl)  * log_gamZ(jA,mosaic.jnh4cl)  +
                      xmol(mosaic.jnh4so4) * log_gamZ(jA,mosaic.jnh4so4) +
                      xmol(mosaic.jna2so4) * log_gamZ(jA,mosaic.jna2so4) +
                      xmol(mosaic.jnano3)  * log_gamZ(jA,mosaic.jnano3)  +
                      xmol(mosaic.jnacl)   * log_gamZ(jA,mosaic.jnacl)   +
                      xmol(mosaic.jcano3)  * log_gamZ(jA,mosaic.jcano3)  +
                      xmol(mosaic.jcacl2)  * log_gamZ(jA,mosaic.jcacl2)  +
                      xmol(mosaic.jhno3)   * log_gamZ(jA,mosaic.jhno3)   +
                      xmol(mosaic.jhcl)    * log_gamZ(jA,mosaic.jhcl);
        gam(jA) = ats<real_type>::pow(10., log_gam(jA));
        activity(jA) = mc(mosaic.jc_nh4) *
                      ma(mosaic.ja_cl) *
                      ats<real_type>::pow(gam(jA),2.);
      }

      jA = mosaic.jna2so4;
      if (xmol(jA) > 0.0) {
        log_gam(jA) = xmol(mosaic.jnh4no3) * log_gamZ(jA,mosaic.jnh4no3) +
                      xmol(mosaic.jnh4cl)  * log_gamZ(jA,mosaic.jnh4cl)  +
                      xmol(mosaic.jnh4so4) * log_gamZ(jA,mosaic.jnh4so4) +
                      xmol(mosaic.jna2so4) * log_gamZ(jA,mosaic.jna2so4) +
                      xmol(mosaic.jnano3)  * log_gamZ(jA,mosaic.jnano3)  +
                      xmol(mosaic.jnacl)   * log_gamZ(jA,mosaic.jnacl)   +
                      xmol(mosaic.jcano3)  * log_gamZ(jA,mosaic.jcano3)  +
                      xmol(mosaic.jcacl2)  * log_gamZ(jA,mosaic.jcacl2)  +
                      xmol(mosaic.jhno3)   * log_gamZ(jA,mosaic.jhno3)   +
                      xmol(mosaic.jhcl)    * log_gamZ(jA,mosaic.jhcl);
        gam(jA) = ats<real_type>::pow(10., log_gam(jA));
        activity(jA) = ats<real_type>::pow(mc(mosaic.jc_na),2.) *
                      ma(mosaic.ja_so4) *
                      ats<real_type>::pow(gam(jA),3.);
      }

      jA = mosaic.jnano3;
      if (xmol(jA) > 0.0) {
        log_gam(jA) = xmol(mosaic.jnh4no3) * log_gamZ(jA,mosaic.jnh4no3) +
                      xmol(mosaic.jnh4cl)  * log_gamZ(jA,mosaic.jnh4cl)  +
                      xmol(mosaic.jnh4so4) * log_gamZ(jA,mosaic.jnh4so4) +
                      xmol(mosaic.jna2so4) * log_gamZ(jA,mosaic.jna2so4) +
                      xmol(mosaic.jnano3)  * log_gamZ(jA,mosaic.jnano3)  +
                      xmol(mosaic.jnacl)   * log_gamZ(jA,mosaic.jnacl)   +
                      xmol(mosaic.jcano3)  * log_gamZ(jA,mosaic.jcano3)  +
                      xmol(mosaic.jcacl2)  * log_gamZ(jA,mosaic.jcacl2)  +
                      xmol(mosaic.jhno3)   * log_gamZ(jA,mosaic.jhno3)   +
                      xmol(mosaic.jhcl)    * log_gamZ(jA,mosaic.jhcl);
        gam(jA) = ats<real_type>::pow(10., log_gam(jA));
        activity(jA) = mc(mosaic.jc_na) *
                      ma(mosaic.ja_no3) *
                      ats<real_type>::pow(gam(jA),2.);
      }

      jA = mosaic.jnacl;
      if (xmol(jA) > 0.0) {
        log_gam(jA) = xmol(mosaic.jnh4no3) * log_gamZ(jA,mosaic.jnh4no3) +
                      xmol(mosaic.jnh4cl)  * log_gamZ(jA,mosaic.jnh4cl)  +
                      xmol(mosaic.jnh4so4) * log_gamZ(jA,mosaic.jnh4so4) +
                      xmol(mosaic.jna2so4) * log_gamZ(jA,mosaic.jna2so4) +
                      xmol(mosaic.jnano3)  * log_gamZ(jA,mosaic.jnano3)  +
                      xmol(mosaic.jnacl)   * log_gamZ(jA,mosaic.jnacl)   +
                      xmol(mosaic.jcano3)  * log_gamZ(jA,mosaic.jcano3)  +
                      xmol(mosaic.jcacl2)  * log_gamZ(jA,mosaic.jcacl2)  +
                      xmol(mosaic.jhno3)   * log_gamZ(jA,mosaic.jhno3)   +
                      xmol(mosaic.jhcl)    * log_gamZ(jA,mosaic.jhcl);
        gam(jA) = ats<real_type>::pow(10., log_gam(jA));
        activity(jA) = mc(mosaic.jc_na) *
                      ma(mosaic.ja_cl) *
                      ats<real_type>::pow(gam(jA),2.);
      }

  // Note: these are commented out in MOSAIC also.
  //  jA = mosaic.jcano3;
  //  if (xmol(jA) > 0.) {
  //    gam(jA) = 1.;
  //    activity(jA) = 1.;
  //  }

  //  jA = mosaic.jacl2;
  //  if (xmol(jA) > 0.) {
  //    gam(jA) = 1.;
  //    activity(jA) = 1.;
  //  }

      jA = mosaic.jcano3;
      if (xmol(jA) > 0.0) {
        log_gam(jA) = xmol(mosaic.jnh4no3) * log_gamZ(jA,mosaic.jnh4no3) +
                      xmol(mosaic.jnh4cl)  * log_gamZ(jA,mosaic.jnh4cl)  +
                      xmol(mosaic.jnh4so4) * log_gamZ(jA,mosaic.jnh4so4) +
                      xmol(mosaic.jna2so4) * log_gamZ(jA,mosaic.jna2so4) +
                      xmol(mosaic.jnano3)  * log_gamZ(jA,mosaic.jnano3)  +
                      xmol(mosaic.jnacl)   * log_gamZ(jA,mosaic.jnacl)   +
                      xmol(mosaic.jcano3)  * log_gamZ(jA,mosaic.jcano3)  +
                      xmol(mosaic.jcacl2)  * log_gamZ(jA,mosaic.jcacl2)  +
                      xmol(mosaic.jhno3)   * log_gamZ(jA,mosaic.jhno3)   +
                      xmol(mosaic.jhcl)    * log_gamZ(jA,mosaic.jhcl);
        gam(jA) = ats<real_type>::pow(10., log_gam(jA));
        activity(jA) = mc(mosaic.jc_ca) *
                      ats<real_type>::pow(ma(mosaic.ja_no3),2.) *
                      ats<real_type>::pow(gam(jA),3.);
      }

      jA = mosaic.jcacl2;
      if (xmol(jA) > 0.0) {
        log_gam(jA) = xmol(mosaic.jnh4no3) * log_gamZ(jA,mosaic.jnh4no3) +
                      xmol(mosaic.jnh4cl)  * log_gamZ(jA,mosaic.jnh4cl)  +
                      xmol(mosaic.jnh4so4) * log_gamZ(jA,mosaic.jnh4so4) +
                      xmol(mosaic.jna2so4) * log_gamZ(jA,mosaic.jna2so4) +
                      xmol(mosaic.jnano3)  * log_gamZ(jA,mosaic.jnano3)  +
                      xmol(mosaic.jnacl)   * log_gamZ(jA,mosaic.jnacl)   +
                      xmol(mosaic.jcano3)  * log_gamZ(jA,mosaic.jcano3)  +
                      xmol(mosaic.jcacl2)  * log_gamZ(jA,mosaic.jcacl2)  +
                      xmol(mosaic.jhno3)   * log_gamZ(jA,mosaic.jhno3)   +
                      xmol(mosaic.jhcl)    * log_gamZ(jA,mosaic.jhcl);
        gam(jA) = ats<real_type>::pow(10., log_gam(jA));
        activity(jA) = mc(mosaic.jc_ca) *
                      ats<real_type>::pow(ma(mosaic.ja_cl),2.) *
                      ats<real_type>::pow(gam(jA),3.);
      }

      jA = mosaic.jhno3;
      if (xmol(jA) > 0.0) {
        log_gam(jA) = xmol(mosaic.jnh4no3) * log_gamZ(jA,mosaic.jnh4no3) +
                      xmol(mosaic.jnh4cl)  * log_gamZ(jA,mosaic.jnh4cl)  +
                      xmol(mosaic.jnh4so4) * log_gamZ(jA,mosaic.jnh4so4) +
                      xmol(mosaic.jna2so4) * log_gamZ(jA,mosaic.jna2so4) +
                      xmol(mosaic.jnano3)  * log_gamZ(jA,mosaic.jnano3)  +
                      xmol(mosaic.jnacl)   * log_gamZ(jA,mosaic.jnacl)   +
                      xmol(mosaic.jcano3)  * log_gamZ(jA,mosaic.jcano3)  +
                      xmol(mosaic.jcacl2)  * log_gamZ(jA,mosaic.jcacl2)  +
                      xmol(mosaic.jhno3)   * log_gamZ(jA,mosaic.jhno3)   +
                      xmol(mosaic.jhcl)    * log_gamZ(jA,mosaic.jhcl);
        gam(jA) = ats<real_type>::pow(10., log_gam(jA));
        activity(jA) = mc(mosaic.jc_h) *
                      ma(mosaic.ja_no3) *
                      ats<real_type>::pow(gam(jA),2.);
      }

      jA = mosaic.jhcl;
      if (xmol(jA) > 0.0) {
        log_gam(jA) = xmol(mosaic.jnh4no3) * log_gamZ(jA,mosaic.jnh4no3) +
                      xmol(mosaic.jnh4cl)  * log_gamZ(jA,mosaic.jnh4cl)  +
                      xmol(mosaic.jnh4so4) * log_gamZ(jA,mosaic.jnh4so4) +
                      xmol(mosaic.jna2so4) * log_gamZ(jA,mosaic.jna2so4) +
                      xmol(mosaic.jnano3)  * log_gamZ(jA,mosaic.jnano3)  +
                      xmol(mosaic.jnacl)   * log_gamZ(jA,mosaic.jnacl)   +
                      xmol(mosaic.jcano3)  * log_gamZ(jA,mosaic.jcano3)  +
                      xmol(mosaic.jcacl2)  * log_gamZ(jA,mosaic.jcacl2)  +
                      xmol(mosaic.jhno3)   * log_gamZ(jA,mosaic.jhno3)   +
                      xmol(mosaic.jhcl)    * log_gamZ(jA,mosaic.jhcl);
        gam(jA) = ats<real_type>::pow(10., log_gam(jA));
        activity(jA) = mc(mosaic.jc_h) *
                      ma(mosaic.ja_cl) *
                      ats<real_type>::pow(gam(jA),2.);
      }

      // FIXME: duplicated code to avoid goto statement
      gam(mosaic.jlvcite)  = 1.0;
      gam(mosaic.jnh4hso4) = 1.0;
      gam(mosaic.jnh4msa)  = 1.0;
      gam(mosaic.jna3hso4) = 1.0;
      gam(mosaic.jnamsa)   = 1.0;
      gam(mosaic.jcamsa2)  = 1.0;

      activity(mosaic.jlvcite)  = 0.0;
      activity(mosaic.jnh4hso4) = 0.0;

      activity(mosaic.jnh4msa) = mc(mosaic.jc_nh4) * ma(mosaic.ja_msa) *
                                gam(mosaic.jnh4msa) * gam(mosaic.jnh4msa);

      activity(mosaic.jna3hso4) = 0.0;
      activity(mosaic.jnahso4)  = 0.0;

      activity(mosaic.jnamsa)  = mc(mosaic.jc_na) * ma(mosaic.ja_msa) *
                                gam(mosaic.jnamsa) * gam(mosaic.jnamsa);
      activity(mosaic.jcamsa2) = mc(mosaic.jc_ca) * ma(mosaic.ja_msa) *
                                gam(mosaic.jcamsa2) * gam(mosaic.jcamsa2) * gam(mosaic.jcamsa2);

      gam_ratio = gam(mosaic.jnh4no3) * gam(mosaic.jnh4no3) /
                  gam(mosaic.jhno3)   * gam(mosaic.jhno3);

    // end SULFATE POOR regime
    } else {
      // SULFATE RICH: solve for SO4= and HSO4- ions

      sum_elec = 3. * electrolyte_liquid(mosaic.jh2so4)   +
                 2. * electrolyte_liquid(mosaic.jnh4hso4) +
                 5. * electrolyte_liquid(mosaic.jlvcite)  +
                 3. * electrolyte_liquid(mosaic.jnh4so4)  +
                 2. * electrolyte_liquid(mosaic.jnahso4)  +
                 5. * electrolyte_liquid(mosaic.jna3hso4) +
                 3. * electrolyte_liquid(mosaic.jna2so4)  +
                 2. * electrolyte_liquid(mosaic.jhno3)    +
                 2. * electrolyte_liquid(mosaic.jhcl);

      if (sum_elec == 0.0) {
        for (ordinal_type jA = 0; jA < mosaic.nelectrolyte; jA++) {
          gam(jA) = 1.0;
        }
        // FIXME: duplicated code to avoid goto statement
        gam(mosaic.jnh4no3) = 1.0;
        gam(mosaic.jnh4cl)  = 1.0;
        gam(mosaic.jnano3)  = 1.0;
        gam(mosaic.jnacl)   = 1.0;
        gam(mosaic.jcano3)  = 1.0;
        gam(mosaic.jcacl2)  = 1.0;
        gam(mosaic.jnh4msa) = 1.0;
        gam(mosaic.jnamsa)  = 1.0;
        gam(mosaic.jcamsa2) = 1.0;

      // compute equilibrium pH
      // cation molalities (mol / kg water)
        mc(mosaic.jc_ca)  = 1.e-9 * aer_liquid(mosaic.ica_a)  / water_a;
        mc(mosaic.jc_nh4) = 1.e-9 * aer_liquid(mosaic.inh4_a) / water_a;
        mc(mosaic.jc_na)  = 1.e-9 * aer_liquid(mosaic.ina_a)  / water_a;

      // anion molalities (mol / kg water)
        mSULF              = 1.e-9 * aer_liquid(mosaic.iso4_a) / water_a;
        ma(mosaic.ja_hso4) = 0.0;
        ma(mosaic.ja_so4)  = 0.0;
        ma(mosaic.ja_no3)  = 1.e-9 * aer_liquid(mosaic.ino3_a) / water_a;
        ma(mosaic.ja_cl)   = 1.e-9 * aer_liquid(mosaic.icl_a)  / water_a;
        ma(mosaic.ja_msa)  = 1.e-9 * aer_liquid(mosaic.imsa_a) / water_a;

        real_type dumK, c_bal, aq, bq, cq;

        gam_ratio = ats<real_type>::pow(gam(mosaic.jnh4hso4),2.) /
                    ats<real_type>::pow(gam(mosaic.jhhso4),2.);
        dumK = Keq_ll(0) * ats<real_type>::pow(gam(mosaic.jhhso4),2.) /
                           ats<real_type>::pow(gam(mosaic.jh2so4),3.);

        c_bal = mc(mosaic.jc_nh4) + mc(mosaic.jc_na) + 2. * mc(mosaic.jc_ca) -
                ma(mosaic.ja_no3) - ma(mosaic.ja_cl) - mSULF - ma(mosaic.ja_msa);

        aq = 1.0;
        bq = dumK + c_bal;
        cq = dumK * (c_bal -mSULF);

        real_type xq = 0.0;
      //--quadratic solution
        if (bq != 0.0) {
          xq = 4. * (1./bq) * (cq/bq);

        } else {
          xq = 1.e6;
        }

        real_type quad = 0.0;
        if(ats<real_type>::abs(xq) < 1.e-6) {
          real_type dum = xq * (0.5 + xq*(0.125 + xq*0.0625));
          quad  = (-0.5*bq/aq)*dum;
          if (quad < 0.0) {
            quad = -bq/aq - quad;
          }
        } else {
          quad = 0.5 * (-bq + ats<real_type>::sqrt(bq*bq - 4.*cq));
        }

      //--end of quadratic solution
        mc(mosaic.jc_h) = quad > 1.e-7 ? quad : 1.e-7;
        ma(mosaic.ja_so4)  = mSULF * dumK / (mc(mosaic.jc_h) + dumK);
        ma(mosaic.ja_hso4) = mSULF - ma(mosaic.ja_so4);

        activity(mosaic.jcamsa2) = mc(mosaic.jc_ca) *
                                  ats<real_type>::pow(ma(mosaic.ja_msa),2.) *
                                  ats<real_type>::pow(gam(mosaic.jcamsa2),3.);

        activity(mosaic.jnh4so4) = ats<real_type>::pow(mc(mosaic.jc_nh4),2.) *
                                  ma(mosaic.ja_so4) *
                                  ats<real_type>::pow(gam(mosaic.jnh4so4),3.);

        activity(mosaic.jlvcite) = ats<real_type>::pow(mc(mosaic.jc_nh4),3.) *
                                  ma(mosaic.ja_hso4) *
                                  ma(mosaic.ja_so4) *
                                  ats<real_type>::pow(gam(mosaic.jlvcite),5.);

        activity(mosaic.jnh4hso4) = mc(mosaic.jc_nh4) *
                                    ma(mosaic.ja_hso4) *
                                    ats<real_type>::pow(gam(mosaic.jnh4hso4),2.);

        activity(mosaic.jnh4msa) = mc(mosaic.jc_nh4) *
                                  ma(mosaic.ja_msa) *
                                  ats<real_type>::pow(gam(mosaic.jnh4msa),2.);

        activity(mosaic.jna2so4) = ats<real_type>::pow(mc(mosaic.jc_na),2.) *
                                  ma(mosaic.ja_so4) *
                                  ats<real_type>::pow(gam(mosaic.jna2so4),3.);

        activity(mosaic.jnahso4) = mc(mosaic.jc_na) *
                                  ma(mosaic.ja_hso4) *
                                  ats<real_type>::pow(gam(mosaic.jnahso4),2.);

        activity(mosaic.jnamsa)  = mc(mosaic.jc_na) *
                                  ma(mosaic.ja_msa) *
                                  ats<real_type>::pow(gam(mosaic.jnamsa),2.);

  // Note: these lines are also commented out in MOSAIC
  //      activity(jna3hso4,ibin)= mc(jc_na,ibin)**3 * ma(ja_hso4,ibin) *
  //     &                         ma(ja_so4,ibin) * gam(jna3hso4,ibin)**5

        activity(mosaic.jna3hso4) = 0.0;

        activity(mosaic.jhno3) = mc(mosaic.jc_h) *
                                ma(mosaic.ja_no3) *
                                ats<real_type>::pow(gam(mosaic.jhno3),2.);

        activity(mosaic.jhcl)  = mc(mosaic.jc_h) *
                                ma(mosaic.ja_cl) *
                                ats<real_type>::pow(gam(mosaic.jhcl),2.);

        activity(mosaic.jmsa)  = mc(mosaic.jc_h) *
                                ma(mosaic.ja_msa) *
                                ats<real_type>::pow(gam(mosaic.jmsa),2.);

  // sulfate-poor species
        activity(mosaic.jnh4no3) = 0.0;
        activity(mosaic.jnh4cl)  = 0.0;
        activity(mosaic.jnano3)  = 0.0;
        activity(mosaic.jnacl)   = 0.0;
        activity(mosaic.jcano3)  = 0.0;
        activity(mosaic.jcacl2)  = 0.0;
      }

      xmol(mosaic.jh2so4)   = 3. * electrolyte_liquid(mosaic.jh2so4)   / sum_elec;
      xmol(mosaic.jnh4hso4) = 2. * electrolyte_liquid(mosaic.jnh4hso4) / sum_elec;
      xmol(mosaic.jlvcite)  = 5. * electrolyte_liquid(mosaic.jlvcite)  / sum_elec;
      xmol(mosaic.jnh4so4)  = 3. * electrolyte_liquid(mosaic.jnh4so4)  / sum_elec;
      xmol(mosaic.jnahso4)  = 2. * electrolyte_liquid(mosaic.jnahso4)  / sum_elec;
      xmol(mosaic.jna3hso4) = 5. * electrolyte_liquid(mosaic.jna3hso4) / sum_elec;
      xmol(mosaic.jna2so4)  = 3. * electrolyte_liquid(mosaic.jna2so4)  / sum_elec;
      xmol(mosaic.jhno3)    = 2. * electrolyte_liquid(mosaic.jhno3)    / sum_elec;
      xmol(mosaic.jhcl)     = 2. * electrolyte_liquid(mosaic.jhcl)     / sum_elec;

    // 2H.SO4
      jA = mosaic.jh2so4;
      log_gam(jA) = xmol(mosaic.jh2so4)   * log_gamZ(jA,mosaic.jh2so4)   +
                    xmol(mosaic.jnh4hso4) * log_gamZ(jA,mosaic.jnh4hso4) +
                    xmol(mosaic.jlvcite)  * log_gamZ(jA,mosaic.jlvcite)  +
                    xmol(mosaic.jnh4so4)  * log_gamZ(jA,mosaic.jnh4so4)  +
                    xmol(mosaic.jnahso4)  * log_gamZ(jA,mosaic.jnahso4)  +
                    xmol(mosaic.jna3hso4) * log_gamZ(jA,mosaic.jna3hso4) +
                    xmol(mosaic.jna2so4)  * log_gamZ(jA,mosaic.jna2so4)  +
                    xmol(mosaic.jhno3)    * log_gamZ(jA,mosaic.jhno3)    +
                    xmol(mosaic.jhcl)     * log_gamZ(jA,mosaic.jhcl);
      gam(jA) = ats<real_type>::pow(10.,log_gam(jA));

    // H.HSO4
      jA = mosaic.jhhso4;
      log_gam(jA) = xmol(mosaic.jh2so4)   * log_gamZ(jA,mosaic.jh2so4)   +
                    xmol(mosaic.jnh4hso4) * log_gamZ(jA,mosaic.jnh4hso4) +
                    xmol(mosaic.jlvcite)  * log_gamZ(jA,mosaic.jlvcite)  +
                    xmol(mosaic.jnh4so4)  * log_gamZ(jA,mosaic.jnh4so4)  +
                    xmol(mosaic.jnahso4)  * log_gamZ(jA,mosaic.jnahso4)  +
                    xmol(mosaic.jna3hso4) * log_gamZ(jA,mosaic.jna3hso4) +
                    xmol(mosaic.jna2so4)  * log_gamZ(jA,mosaic.jna2so4)  +
                    xmol(mosaic.jhno3)    * log_gamZ(jA,mosaic.jhno3)    +
                    xmol(mosaic.jhcl)     * log_gamZ(jA,mosaic.jhcl);
      gam(jA) = ats<real_type>::pow(10.,log_gam(jA));


    // NH4HSO4
      jA = mosaic.jnh4hso4;
      log_gam(jA) = xmol(mosaic.jh2so4)   * log_gamZ(jA,mosaic.jh2so4)   +
                    xmol(mosaic.jnh4hso4) * log_gamZ(jA,mosaic.jnh4hso4) +
                    xmol(mosaic.jlvcite)  * log_gamZ(jA,mosaic.jlvcite)  +
                    xmol(mosaic.jnh4so4)  * log_gamZ(jA,mosaic.jnh4so4)  +
                    xmol(mosaic.jnahso4)  * log_gamZ(jA,mosaic.jnahso4)  +
                    xmol(mosaic.jna3hso4) * log_gamZ(jA,mosaic.jna3hso4) +
                    xmol(mosaic.jna2so4)  * log_gamZ(jA,mosaic.jna2so4)  +
                    xmol(mosaic.jhno3)    * log_gamZ(jA,mosaic.jhno3)    +
                    xmol(mosaic.jhcl)     * log_gamZ(jA,mosaic.jhcl);
      gam(jA) = ats<real_type>::pow(10.,log_gam(jA));


    // LETOVICITE
      jA = mosaic.jlvcite;
      log_gam(jA) = xmol(mosaic.jh2so4)   * log_gamZ(jA,mosaic.jh2so4)   +
                    xmol(mosaic.jnh4hso4) * log_gamZ(jA,mosaic.jnh4hso4) +
                    xmol(mosaic.jlvcite)  * log_gamZ(jA,mosaic.jlvcite)  +
                    xmol(mosaic.jnh4so4)  * log_gamZ(jA,mosaic.jnh4so4)  +
                    xmol(mosaic.jnahso4)  * log_gamZ(jA,mosaic.jnahso4)  +
                    xmol(mosaic.jna3hso4) * log_gamZ(jA,mosaic.jna3hso4) +
                    xmol(mosaic.jna2so4)  * log_gamZ(jA,mosaic.jna2so4)  +
                    xmol(mosaic.jhno3)    * log_gamZ(jA,mosaic.jhno3)    +
                    xmol(mosaic.jhcl)     * log_gamZ(jA,mosaic.jhcl);
      gam(jA) = ats<real_type>::pow(10.,log_gam(jA));


    // (NH4)2SO4
      jA = mosaic.jnh4so4;
      log_gam(jA) = xmol(mosaic.jh2so4)   * log_gamZ(jA,mosaic.jh2so4)   +
                    xmol(mosaic.jnh4hso4) * log_gamZ(jA,mosaic.jnh4hso4) +
                    xmol(mosaic.jlvcite)  * log_gamZ(jA,mosaic.jlvcite)  +
                    xmol(mosaic.jnh4so4)  * log_gamZ(jA,mosaic.jnh4so4)  +
                    xmol(mosaic.jnahso4)  * log_gamZ(jA,mosaic.jnahso4)  +
                    xmol(mosaic.jna3hso4) * log_gamZ(jA,mosaic.jna3hso4) +
                    xmol(mosaic.jna2so4)  * log_gamZ(jA,mosaic.jna2so4)  +
                    xmol(mosaic.jhno3)    * log_gamZ(jA,mosaic.jhno3)    +
                    xmol(mosaic.jhcl)     * log_gamZ(jA,mosaic.jhcl);
      gam(jA) = ats<real_type>::pow(10.,log_gam(jA));


    // NaHSO4
      jA = mosaic.jnahso4;
      log_gam(jA) = xmol(mosaic.jh2so4)   * log_gamZ(jA,mosaic.jh2so4)   +
                    xmol(mosaic.jnh4hso4) * log_gamZ(jA,mosaic.jnh4hso4) +
                    xmol(mosaic.jlvcite)  * log_gamZ(jA,mosaic.jlvcite)  +
                    xmol(mosaic.jnh4so4)  * log_gamZ(jA,mosaic.jnh4so4)  +
                    xmol(mosaic.jnahso4)  * log_gamZ(jA,mosaic.jnahso4)  +
                    xmol(mosaic.jna3hso4) * log_gamZ(jA,mosaic.jna3hso4) +
                    xmol(mosaic.jna2so4)  * log_gamZ(jA,mosaic.jna2so4)  +
                    xmol(mosaic.jhno3)    * log_gamZ(jA,mosaic.jhno3)    +
                    xmol(mosaic.jhcl)     * log_gamZ(jA,mosaic.jhcl);
      gam(jA) = ats<real_type>::pow(10.,log_gam(jA));


 // Na3H(SO4)2
      jA = mosaic.jna3hso4;
 // Note: commented out in MOSAIC also.
 //      log_gam(jA) = xmol(jh2so4)  *log_gamZ(jA,jh2so4)  +
 //     &              xmol(jnh4hso4)*log_gamZ(jA,jnh4hso4)+
 //     &              xmol(jlvcite) *log_gamZ(jA,jlvcite) +
 //     &              xmol(jnh4so4) *log_gamZ(jA,jnh4so4) +
 //     &              xmol(jnahso4) *log_gamZ(jA,jnahso4) +
 //     &              xmol(jna3hso4)*log_gamZ(jA,jna3hso4)+
 //     &              xmol(jna2so4) *log_gamZ(jA,jna2so4) +
 //     &              xmol(jhno3)   *log_gamZ(jA,jhno3)   +
 //     &              xmol(jhcl)    *log_gamZ(jA,jhcl)
 //      gam(jA,ibin) = 10.**log_gam(jA)
      gam(jA) = 1.0;


      // Na2SO4
      jA = mosaic.jna2so4;
      log_gam(jA) = xmol(mosaic.jh2so4)   * log_gamZ(jA,mosaic.jh2so4)   +
                    xmol(mosaic.jnh4hso4) * log_gamZ(jA,mosaic.jnh4hso4) +
                    xmol(mosaic.jlvcite)  * log_gamZ(jA,mosaic.jlvcite)  +
                    xmol(mosaic.jnh4so4)  * log_gamZ(jA,mosaic.jnh4so4)  +
                    xmol(mosaic.jnahso4)  * log_gamZ(jA,mosaic.jnahso4)  +
                    xmol(mosaic.jna3hso4) * log_gamZ(jA,mosaic.jna3hso4) +
                    xmol(mosaic.jna2so4)  * log_gamZ(jA,mosaic.jna2so4)  +
                    xmol(mosaic.jhno3)    * log_gamZ(jA,mosaic.jhno3)    +
                    xmol(mosaic.jhcl)     * log_gamZ(jA,mosaic.jhcl);
      gam(jA) = ats<real_type>::pow(10.,log_gam(jA));


      // HNO3
      jA = mosaic.jhno3;
      log_gam(jA) = xmol(mosaic.jh2so4)   * log_gamZ(jA,mosaic.jh2so4)   +
                    xmol(mosaic.jnh4hso4) * log_gamZ(jA,mosaic.jnh4hso4) +
                    xmol(mosaic.jlvcite)  * log_gamZ(jA,mosaic.jlvcite)  +
                    xmol(mosaic.jnh4so4)  * log_gamZ(jA,mosaic.jnh4so4)  +
                    xmol(mosaic.jnahso4)  * log_gamZ(jA,mosaic.jnahso4)  +
                    xmol(mosaic.jna3hso4) * log_gamZ(jA,mosaic.jna3hso4) +
                    xmol(mosaic.jna2so4)  * log_gamZ(jA,mosaic.jna2so4)  +
                    xmol(mosaic.jhno3)    * log_gamZ(jA,mosaic.jhno3)    +
                    xmol(mosaic.jhcl)     * log_gamZ(jA,mosaic.jhcl);
      gam(jA) = ats<real_type>::pow(10.,log_gam(jA));


      // HCl
      jA = mosaic.jhcl;
      log_gam(jA) = xmol(mosaic.jh2so4)   * log_gamZ(jA,mosaic.jh2so4)   +
                    xmol(mosaic.jnh4hso4) * log_gamZ(jA,mosaic.jnh4hso4) +
                    xmol(mosaic.jlvcite)  * log_gamZ(jA,mosaic.jlvcite)  +
                    xmol(mosaic.jnh4so4)  * log_gamZ(jA,mosaic.jnh4so4)  +
                    xmol(mosaic.jnahso4)  * log_gamZ(jA,mosaic.jnahso4)  +
                    xmol(mosaic.jna3hso4) * log_gamZ(jA,mosaic.jna3hso4) +
                    xmol(mosaic.jna2so4)  * log_gamZ(jA,mosaic.jna2so4)  +
                    xmol(mosaic.jhno3)    * log_gamZ(jA,mosaic.jhno3)    +
                    xmol(mosaic.jhcl)     * log_gamZ(jA,mosaic.jhcl);
      gam(jA) = ats<real_type>::pow(10.,log_gam(jA));

      // FIXME: duplicated code to avoid goto statement
      gam(mosaic.jnh4no3) = 1.0;
      gam(mosaic.jnh4cl)  = 1.0;
      gam(mosaic.jnano3)  = 1.0;
      gam(mosaic.jnacl)   = 1.0;
      gam(mosaic.jcano3)  = 1.0;
      gam(mosaic.jcacl2)  = 1.0;
      gam(mosaic.jnh4msa) = 1.0;
      gam(mosaic.jnamsa)  = 1.0;
      gam(mosaic.jcamsa2) = 1.0;

      // compute equilibrium pH
      // cation molalities (mol / kg water)
      mc(mosaic.jc_ca)  = 1.e-9 * aer_liquid(mosaic.ica_a)  / water_a;
      mc(mosaic.jc_nh4) = 1.e-9 * aer_liquid(mosaic.inh4_a) / water_a;
      mc(mosaic.jc_na)  = 1.e-9 * aer_liquid(mosaic.ina_a)  / water_a;

      // anion molalities (mol / kg water)
      mSULF              = 1.e-9 * aer_liquid(mosaic.iso4_a) / water_a;
      ma(mosaic.ja_hso4) = 0.0;
      ma(mosaic.ja_so4)  = 0.0;
      ma(mosaic.ja_no3)  = 1.e-9 * aer_liquid(mosaic.ino3_a) / water_a;
      ma(mosaic.ja_cl)   = 1.e-9 * aer_liquid(mosaic.icl_a)  / water_a;
      ma(mosaic.ja_msa)  = 1.e-9 * aer_liquid(mosaic.imsa_a) / water_a;

      real_type dumK, c_bal, aq, bq, cq;

      gam_ratio = ats<real_type>::pow(gam(mosaic.jnh4hso4),2.) /
                  ats<real_type>::pow(gam(mosaic.jhhso4),2.);
      dumK = Keq_ll(0) * ats<real_type>::pow(gam(mosaic.jhhso4),2.) /
                         ats<real_type>::pow(gam(mosaic.jh2so4),3.);


      c_bal = mc(mosaic.jc_nh4) + mc(mosaic.jc_na) + 2. * mc(mosaic.jc_ca) -
              ma(mosaic.ja_no3) - ma(mosaic.ja_cl) - mSULF - ma(mosaic.ja_msa);

      aq = 1.0;
      bq = dumK + c_bal;
      cq = dumK * (c_bal -mSULF);

      //--quadratic solution
      real_type xq = 0.0;
      if (bq != 0.0) {
        xq = 4. * (1./bq) * (cq/bq);
      } else {
        xq = 1.e6;
      }

      real_type quad = 0.0;
      if(ats<real_type>::abs(xq) < 1.e-6) {
        real_type dum = xq * (0.5 + xq*(0.125 + xq*0.0625));
        quad  = (-0.5*bq/aq)*dum;
        if (quad < 0.0) {
          quad = -bq/aq - quad;
        }
      } else {
        quad = 0.5 * (-bq + ats<real_type>::sqrt(bq*bq - 4.*cq));
      }

      //--end of quadratic solution
      mc(mosaic.jc_h) = quad > 1.e-7 ? quad : 1.e-7;
      ma(mosaic.ja_so4)  = mSULF * dumK / (mc(mosaic.jc_h) + dumK);
      ma(mosaic.ja_hso4) = mSULF - ma(mosaic.ja_so4);

      activity(mosaic.jcamsa2) = mc(mosaic.jc_ca) *
                                ats<real_type>::pow(ma(mosaic.ja_msa),2.) *
                                ats<real_type>::pow(gam(mosaic.jcamsa2),3.);

      activity(mosaic.jnh4so4) = ats<real_type>::pow(mc(mosaic.jc_nh4),2.) *
                                ma(mosaic.ja_so4) *
                                ats<real_type>::pow(gam(mosaic.jnh4so4),3.);

      activity(mosaic.jlvcite) = ats<real_type>::pow(mc(mosaic.jc_nh4),3.) *
                                ma(mosaic.ja_hso4) *
                                ma(mosaic.ja_so4) *
                                ats<real_type>::pow(gam(mosaic.jlvcite),5.);

      activity(mosaic.jnh4hso4) = mc(mosaic.jc_nh4) *
                                  ma(mosaic.ja_hso4) *
                                  ats<real_type>::pow(gam(mosaic.jnh4hso4),2.);

      activity(mosaic.jnh4msa) = mc(mosaic.jc_nh4) *
                                ma(mosaic.ja_msa) *
                                ats<real_type>::pow(gam(mosaic.jnh4msa),2.);

      activity(mosaic.jna2so4) = ats<real_type>::pow(mc(mosaic.jc_na),2.) *
                                ma(mosaic.ja_so4) *
                                ats<real_type>::pow(gam(mosaic.jna2so4),3.);

      activity(mosaic.jnahso4) = mc(mosaic.jc_na) *
                                ma(mosaic.ja_hso4) *
                                ats<real_type>::pow(gam(mosaic.jnahso4),2.);

      activity(mosaic.jnamsa)  = mc(mosaic.jc_na) *
                                ma(mosaic.ja_msa) *
                                ats<real_type>::pow(gam(mosaic.jnamsa),2.);

  // Note: these lines are also commented out in MOSAIC.
  //      activity(jna3hso4,ibin)= mc(jc_na,ibin)**3 * ma(ja_hso4,ibin) *
  //     &                         ma(ja_so4,ibin) * gam(jna3hso4,ibin)**5

      activity(mosaic.jna3hso4) = 0.0;

      activity(mosaic.jhno3) = mc(mosaic.jc_h) *
                              ma(mosaic.ja_no3) *
                              ats<real_type>::pow(gam(mosaic.jhno3),2.);

      activity(mosaic.jhcl)  = mc(mosaic.jc_h) *
                              ma(mosaic.ja_cl) *
                              ats<real_type>::pow(gam(mosaic.jhcl),2.);

      activity(mosaic.jmsa)  = mc(mosaic.jc_h) *
                              ma(mosaic.ja_msa) *
                              ats<real_type>::pow(gam(mosaic.jmsa),2.);

  // sulfate-poor species
      activity(mosaic.jnh4no3) = 0.0;
      activity(mosaic.jnh4cl)  = 0.0;
      activity(mosaic.jnano3)  = 0.0;
      activity(mosaic.jnacl)   = 0.0;
      activity(mosaic.jcano3)  = 0.0;
      activity(mosaic.jcacl2)  = 0.0;
    }
  } // compute_activities

# endif