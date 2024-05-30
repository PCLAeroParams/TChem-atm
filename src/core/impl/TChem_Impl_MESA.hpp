#ifndef __TCHEM_IMPL_MESA_HPP__
#define __TCHEM_IMPL_MESA_HPP__

/* MESA: Multicomponent Equilibrium Solver for Aerosols.
 Computes equilibrum solid and liquid phases by integrating
 pseudo-transient dissolution and precipitation reactions

 author: Rahul A. Zaveri */

#include "TChem_KineticModelData.hpp"
#include "TChem_Util.hpp"
// #define TCHEM_ENABLE_SERIAL_TEST_OUTPUT
namespace TChem {
namespace Impl {

struct MosaicModelData {
  public:
    const ordinal_type iso4_a = 0;   // <-> ih2so4_g
    const ordinal_type ino3_a = 1;   // <-> ihno3_g
    const ordinal_type icl_a = 2;    // <-> ihcl_g
    const ordinal_type inh4_a = 3;   // <-> inh3_g
    const ordinal_type imsa_a = 4;   // <-> imsa_g
    const ordinal_type iaro1_a = 5;  // <-> iaro1_g
    const ordinal_type iaro2_a = 6;  // <-> iaro2_g
    const ordinal_type ialk1_a = 7;  // <-> ialk1_g
    const ordinal_type iole1_a = 8;  // <-> iole1_g
    const ordinal_type iapi1_a = 9;  // <-> iapi1_g
    const ordinal_type iapi2_a = 10; // <-> iapi2_g
    const ordinal_type ilim1_a = 11; // <-> ilim1_g
    const ordinal_type ilim2_a = 12; // <-> ilim2_g
    const ordinal_type ico3_a = 13;  // <-> ico2_g
    const ordinal_type ina_a = 14;
    const ordinal_type ica_a = 15;
    const ordinal_type ioin_a = 16;
    const ordinal_type ioc_a = 17;
    const ordinal_type ibc_a = 18;

    //electrolyte indices (used for water content calculations)
    // these indices are order sensitive


    const ordinal_type jnh4so4 = 0; // soluble
    const ordinal_type jlvcite = 1; // soluble
    const ordinal_type jnh4hso4 = 2; // soluble
    const ordinal_type jnh4msa = 3; // soluble: new
    const ordinal_type jnh4no3 = 4; // soluble
    const ordinal_type jnh4cl = 5; // soluble
    const ordinal_type jna2so4 = 6; // soluble
    const ordinal_type jna3hso4 = 7; // soluble
    const ordinal_type jnahso4 = 8; // soluble
    const ordinal_type jnamsa = 9; // soluble: new
    const ordinal_type jnano3 = 10; // soluble
    const ordinal_type jnacl = 11; // soluble
    const ordinal_type jcano3 = 12; // soluble
    const ordinal_type jcacl2 = 13; // soluble
    const ordinal_type jcamsa2 = 14; // soluble     nsalt
    const ordinal_type jh2so4 = 15; // soluble
    const ordinal_type jmsa = 16; // soluble
    const ordinal_type jhno3 = 17; // soluble
    const ordinal_type jhcl = 18; // soluble
    const ordinal_type jhhso4 = 19; // soluble
    const ordinal_type jcaso4 = 20; // insoluble
    const ordinal_type jcaco3 = 21; // insoluble
    const ordinal_type joc = 22; // insoluble - part of naercomp
    const ordinal_type jbc = 23; // insoluble - part of naercomp
    const ordinal_type join = 24; // insoluble - part of naercomp
    const ordinal_type jaro1 = 25; // insoluble - part of naercomp
    const ordinal_type jaro2 = 26; // insoluble - part of naercomp
    const ordinal_type jalk1 = 27; // insoluble - part of naercomp
    const ordinal_type jole1 = 28; // insoluble - part of naercomp
    const ordinal_type japi1 = 29; // insoluble - part of naercomp
    const ordinal_type japi2 = 30; // insoluble - part of naercomp
    const ordinal_type jlim1 = 31; // insoluble - part of naercomp
    const ordinal_type jlim2 = 32; // insoluble - part of naercomp
    const ordinal_type jh2o = 33; // water - part of naercomp

    const ordinal_type naer=1;
    const ordinal_type nelectrolyte=1;

    // FIXME: update values and make enum
    const ordinal_type jliquid=1;
    const ordinal_type jsolid=2;

    const ordinal_type all_solid=3;
    const ordinal_type all_liquid=3;

    const ordinal_type jhyst_up=0;
    const ordinal_type jhyst_lo=0;

    const ordinal_type nsoluble=1;
 };


template<typename ValueType, typename DeviceType>
struct MESA
{
  using value_type = ValueType;
  using device_type = DeviceType;
  using scalar_type = typename ats<value_type>::scalar_type;

  using real_type = scalar_type;
  /// sacado is value type
  using value_type_1d_view_type = Tines::value_type_1d_view<value_type,device_type>;
  using real_type_1d_view_type = Tines::value_type_1d_view<real_type,device_type>;

  using aerosol_model_data_type= AerosolModel_ConstData<device_type>;

  KOKKOS_INLINE_FUNCTION static ordinal_type getWorkSpaceSize()
  {
    ordinal_type workspace_size=0;
    return workspace_size;
  }

  KOKKOS_INLINE_FUNCTION static
  void calculate_XT(const real_type_1d_view_type& aer,
                    const MosaicModelData& mosaic,
                    real_type &XT
                  ) {
  //aer nmol/m^3
  // I remove ibin
    if ((aer(mosaic.iso4_a) + aer(mosaic.imsa_a)) > 0.0) {
        XT = (aer(mosaic.inh4_a) + aer(mosaic.ina_a) +
             2.0 * aer(mosaic.ica_a)) /
            (aer(mosaic.iso4_a) + 0.5 * aer(mosaic.imsa_a));
    } else {
        XT = -1.0;
    }
}// calculate_XT


/*
 this subroutine completely deliquesces an aerosol and partitions
 all the soluble electrolytes into the liquid phase and insoluble
 ones into the solid phase. It also calculates the corresponding
 aer(js,liquid,ibin) and aer(js,jsolid,ibin) generic species
 concentrations

 author: Rahul A. Zaveri
 update: jan 2005
-----------------------------------------------------------------------
*/

KOKKOS_INLINE_FUNCTION static
void adjust_liquid_aerosol(const MosaicModelData& mosaic,
                           const real_type_1d_view_type& electrolyte_solid,
                           const real_type_1d_view_type& electrolyte_liquid,
                           const real_type_1d_view_type& electrolyte_total,
                           const real_type_1d_view_type& epercent_solid,
                           const real_type_1d_view_type& epercent_liquid,
                           const real_type_1d_view_type& epercent_total,
                           const real_type_1d_view_type& aer_solid,
                           const real_type_1d_view_type& aer_liquid,
                           const real_type_1d_view_type& aer_total,
                           real_type& jphase,
                           real_type& jhyst_leg) {

    jphase    = mosaic.jliquid;
    jhyst_leg = mosaic.jhyst_up;// upper curve

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
}

/*
! called when aerosol bin is completely solid.
!
! author: Rahul A. Zaveri
! update: jan 2005
*/

KOKKOS_INLINE_FUNCTION static
void adjust_solid_aerosol(const MosaicModelData& mosaic,
                          const real_type_1d_view_type& aer_solid,
                          const real_type_1d_view_type& aer_liquid,
                          const real_type_1d_view_type& aer_total,
                          const real_type_1d_view_type& electrolyte_solid,
                          const real_type_1d_view_type& electrolyte_liquid,
                          const real_type_1d_view_type& electrolyte_total,
                          const real_type_1d_view_type& epercent_solid,
                          const real_type_1d_view_type& epercent_liquid,
                          const real_type_1d_view_type& epercent_total,
                          real_type& water_a, real_type&jphase, real_type&jhyst_leg ) {
    // Assuming naer, nelectrolyte are accessible or part of the Mosaic object

    // Since water_a is now an input parameter, we can directly modify it.
    water_a = 0.0;
    jphase = mosaic.jsolid;
    jhyst_leg = mosaic.jhyst_lo;// lower curve

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
}// adjust_solid_aerosol


  /*
! called when aerosol bin is completely liquid.
!
! author: Rahul A. Zaveri
! update: jan 2005
  */
  KOKKOS_INLINE_FUNCTION static void
  do_full_deliquescence(const MosaicModelData& mosaic,
                           const real_type_1d_view_type& electrolyte_solid,
                           const real_type_1d_view_type& electrolyte_liquid,
                           const real_type_1d_view_type& electrolyte_total,
                           const real_type_1d_view_type& aer_solid,
                           const real_type_1d_view_type& aer_liquid,
                           const real_type_1d_view_type& aer_total) {
    // Partition all electrolytes ordinal_typeo liquid phase
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
}

  KOKKOS_INLINE_FUNCTION static void
  compute_activities()
  {

  }

KOKKOS_INLINE_FUNCTION static void
  aerosol_water_total()
    {

  }

  template<typename MemberType>
  KOKKOS_INLINE_FUNCTION static void team_invoke(
    const MemberType& member,
    const real_type& t,
    const real_type& p,
    const real_type_1d_view_type& aer_solid,
    const real_type_1d_view_type& aer_liquid,
    const real_type_1d_view_type& aer_total,
    const real_type_1d_view_type& electrolyte_solid,
    const real_type_1d_view_type& electrolyte_liquid,
    const real_type_1d_view_type& electrolyte_total,
    const real_type_1d_view_type& epercent_solid,
    const real_type_1d_view_type& epercent_liquid,
    const real_type_1d_view_type& epercent_total,
    const MosaicModelData& mosaic)
    {

    // total is 3
    // FIXME: aero uses total: sum of solid and liquid.
    // we may sum solid and liquid here.
    // aero can be a view of 1D liquid first and then solid.
    real_type XT=0.0;
    calculate_XT(aer_total,
                 mosaic,
                 XT);

    const real_type CRH = 0.35;

    // FIXME
    real_type ptol_mol_astem =1.0;
    // FIXME are thses inputs?
    real_type jaerosolstate, jphase, jhyst_leg, water_a, aH2O_a=0;
    real_type mass_dry_a, vol_dry_a, growth_factor, vol_wet_a, mass_wet_a;

    if (aH2O_a < CRH &&
        (XT > 1.0 || XT < 0.0) &&
        epercent_total(mosaic.jcano3) <= ptol_mol_astem &&
        epercent_total(mosaic.jcacl2) <= ptol_mol_astem) {
        jaerosolstate = mosaic.all_solid;
        jphase = mosaic.jsolid;
        jhyst_leg = mosaic.jhyst_lo;
        adjust_solid_aerosol(mosaic,
                     aer_solid, aer_liquid, aer_total,
                     electrolyte_solid, electrolyte_liquid, electrolyte_total,
                     epercent_solid, epercent_liquid, epercent_total,
                     water_a, jphase, jhyst_leg);

    return;
    }
    // step 2: Check for supersaturation/metastable state
    // if(water_a_hyst(ibin)*0.0 .gt. 0.5*water_a_up(ibin))then ! 3-D
    if (jhyst_leg == mosaic.jhyst_up) { // BOX
      do_full_deliquescence(mosaic,
                      electrolyte_solid, electrolyte_liquid, electrolyte_total,
                      aer_solid, aer_liquid, aer_total);

      // Assuming electrolyte is a 3D vector or an array
      real_type sum_soluble = 0.0;
      // FIXME nsoluble
      for (ordinal_type js = 0; js < mosaic.nsoluble; ++js) {
        sum_soluble += electrolyte_total(js);
      }

      const real_type solids = electrolyte_total(mosaic.jcaso4) +
                       electrolyte_total(mosaic.jcaco3) +
                       aer_total(mosaic.ioin_a);
#if 1
      if (sum_soluble < 1.e-15 && solids > 0.0) {
        jaerosolstate = mosaic.all_solid; // no soluble material present
        jphase = mosaic.jsolid;
        adjust_solid_aerosol(mosaic,
                     aer_solid, aer_liquid, aer_total,
                     electrolyte_solid, electrolyte_liquid, electrolyte_total,
                     epercent_solid, epercent_liquid, epercent_total,
                     water_a, jphase, jhyst_leg);

        // New wet mass and wet volume
        mass_wet_a = mass_dry_a + water_a * 1.e-3; // g/cc(air)
        vol_wet_a = vol_dry_a + water_a * 1.e-3; // cc(aer)/cc(air) or m^3/m^3(air)
        growth_factor = mass_wet_a / mass_dry_a; // mass growth factor
        return;
    } else if (sum_soluble > 0.0 && solids == 0.0) {
#if 1
        jaerosolstate = mosaic.all_liquid;
        jhyst_leg = mosaic.jhyst_up;
        jphase = mosaic.jliquid;
        // FIXME: aerosol_water_total is an input
        water_a = aerosol_water_total();

        if (water_a < 0.0) {
            jaerosolstate = mosaic.all_solid; // no soluble material present
            jphase = mosaic.jsolid;
            jhyst_leg = mosaic.jhyst_lo;
            adjust_solid_aerosol(mosaic,
                     aer_solid, aer_liquid, aer_total,
                     electrolyte_solid, electrolyte_liquid, electrolyte_total,
                     epercent_solid, epercent_liquid, epercent_total,
                     water_a, jphase, jhyst_leg);
        } else {
            adjust_liquid_aerosol(mosaic,
                      electrolyte_solid, electrolyte_liquid, electrolyte_total,
                      epercent_solid, epercent_liquid, epercent_total,
                      aer_solid, aer_liquid, aer_total,
                      jphase, jhyst_leg);
            compute_activities();
        }

        // New wet mass and wet volume
        mass_wet_a = mass_dry_a + water_a * 1.e-3; // g/cc(air)
        vol_wet_a = vol_dry_a + water_a * 1.e-3; // cc(aer)/cc(air) or m^3/m^3(air)
        growth_factor = mass_wet_a / mass_dry_a; // mass growth factor
#endif
        return;
    }
#endif
    }


    }

};

} // namespace Impl
} // namespace TChem

#endif
