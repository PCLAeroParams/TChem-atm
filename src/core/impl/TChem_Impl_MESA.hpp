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

  struct Mosaic {
    const int iso4_a = 0;   // <-> ih2so4_g
    const int ino3_a = 1;   // <-> ihno3_g
    const int icl_a = 2;    // <-> ihcl_g
    const int inh4_a = 3;   // <-> inh3_g
    const int imsa_a = 4;   // <-> imsa_g
    const int iaro1_a = 5;  // <-> iaro1_g
    const int iaro2_a = 6;  // <-> iaro2_g
    const int ialk1_a = 7;  // <-> ialk1_g
    const int iole1_a = 8;  // <-> iole1_g
    const int iapi1_a = 9;  // <-> iapi1_g
    const int iapi2_a = 10; // <-> iapi2_g
    const int ilim1_a = 11; // <-> ilim1_g
    const int ilim2_a = 12; // <-> ilim2_g
    const int ico3_a = 13;  // <-> ico2_g
    const int ina_a = 14;
    const int ica_a = 15;
    const int ioin_a = 16;
    const int ioc_a = 17;
    const int ibc_a = 18;

    // FIXME: need to update indices. 
    static const int jcaco3;
    static const int jcaso4;
    static const int iso4_a;
    static const int ino3_a;
    static const int icl_a;


    int naer;
    int nelectrolyte;
    int inh4_a; 
    int ino3_a;
    int icl_a;

  };


  KOKKOS_INLINE_FUNCTION static 
  void calculate_XT(const real_type_1d_view_type& aer,
                    const Mosaic& mosaic, 
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
 aer(js,jliquid,ibin) and aer(js,jsolid,ibin) generic species
 concentrations

 author: Rahul A. Zaveri
 update: jan 2005
-----------------------------------------------------------------------
*/
KOKKOS_INLINE_FUNCTION static 
void adjust_liquid_aerosol(const MosaicIndices& mosaic,
                           const real_type_1d_view_type& electrolyte_jsolid, 
                           const real_type_1d_view_type& electrolyte_jliquid, 
                           const real_type_1d_view_type& electrolyte_jtotal, 
                           const real_type_1d_view_type& epercent_jsolid, 
                           const real_type_1d_view_type& epercent_jliquid, 
                           const real_type_1d_view_type& epercent_jtotal, 
                           const real_type_1d_view_type& aer_jsolid, 
                           const real_type_1d_view_type& aer_jliquid, 
                           const real_type_1d_view_type& aer_jtotal
                           real_type& jphase, 
                           real_type& jhyst_leg) {

    jphase    = jliquid
    jhyst_leg = jhyst_up;// upper curve

    // Partition all electrolytes into liquid phase
    for (int je = 0; je < mosaic.nelectrolyte; ++je) {
        electrolyte_jsolid(je) = 0.0;
        epercent_jsolid(je) = 0.0;
        electrolyte_jliquid(je) = electrolyte_jtotal(je);
        epercent_jliquid(je) = epercent_jtotal(je);
    }

    // Except these electrolytes, which always remain in the solid phase
    electrolyte_jsolid(mosaic.jcaco3) = electrolyte_jtotal(mosaic.jcaco3);
    electrolyte_jsolid(mosaic.jcaso4) = electrolyte_jtotal(mosaic.jcaso4);
    epercent_jsolid(mosaic.jcaco3) = epercent_jtotal(mosaic.jcaco3);
    epercent_jsolid(mosaic.jcaso4) = epercent_jtotal(mosaic.jcaso4);
    electrolyte_jliquid(mosaic.jcaco3) = 0.0;
    electrolyte_jliquid(mosaic.jcaso4) = 0.0;
    epercent_jliquid(mosaic.jcaco3) = 0.0;
    epercent_jliquid(mosaic.jcaso4) = 0.0;

    // Partition all the aer species into solid and liquid phases
    // Solid phase
    aer_jsolid(mosaic.iso4_a) = electrolyte_jsolid(mosaic.jcaso4);
    aer_jsolid(mosaic.ino3_a) = 0.0;
    aer_jsolid(mosaic.icl_a) = 0.0;
    aer_jsolid(mosaic.inh4_a) = 0.0;
    aer_jsolid(mosaic.ioc_a) = aer_jtotal(mosaic.ioc_a);
    aer_jsolid(mosaic.imsa_a) = 0.0;
    aer_jsolid(mosaic.ico3_a) = aer_jtotal(mosaic.ico3_a);
    aer_jsolid(mosaic.ina_a) = 0.0;
    aer_jsolid(mosaic.ica_a) = electrolyte_jsolid(mosaic.jcaco3) + electrolyte_jsolid(mosaic.jcaso4);
    aer_jsolid(mosaic.ibc_a) = aer_jtotal(mosaic.ibc_a);
    aer_jsolid(mosaic.ioin_a) = aer_jtotal(mosaic.ioin_a);
    aer_jsolid(mosaic.iaro1_a) = aer_jtotal(mosaic.iaro1_a);
    aer_jsolid(mosaic.iaro2_a) = aer_jtotal(mosaic.iaro2_a);
    aer_jsolid(mosaic.ialk1_a) = aer_jtotal(mosaic.ialk1_a);
    aer_jsolid(mosaic.iole1_a) = aer_jtotal(mosaic.iole1_a);
    aer_jsolid(mosaic.iapi1_a) = aer_jtotal(mosaic.iapi1_a);
    aer_jsolid(mosaic.iapi2_a) = aer_jtotal(mosaic.iapi2_a);
    aer_jsolid(mosaic.ilim1_a) = aer_jtotal(mosaic.ilim1_a);
    aer_jsolid(mosaic.ilim2_a) = aer_jtotal(mosaic.ilim2_a);

    // Liquid phase
    aer_jliquid(mosaic.iso4_a) = max(0.0, aer_jtotal(mosaic.iso4_a) - aer_jsolid(mosaic.iso4_a));
    aer_jliquid(mosaic.ino3_a) = aer_jtotal(mosaic.ino3_a);
    aer_jliquid(mosaic.icl_a) = aer_jtotal(mosaic.icl_a);
    aer_jliquid(mosaic.inh4_a) = aer_jtotal(mosaic.inh4_a);
    aer_jliquid(mosaic.ioc_a) = 0.0;
    aer_jliquid(mosaic.imsa_a) = aer_jtotal(mosaic.imsa_a);
    aer_jliquid(mosaic.ico3_a) = 0.0;
    aer_jliquid(mosaic.ina_a) = aer_jtotal(mosaic.ina_a);
    aer_jliquid(mosaic.ica_a) = max(0.0, aer_jtotal(mosaic.ica_a) - aer_jsolid(mosaic.ica_a));
    aer_jliquid(mosaic.ibc_a) = 0.0;
    aer_jliquid(mosaic.ioin_a) = 0.0;
    aer_jliquid(mosaic.iaro1_a) = 0.0;
    aer_jliquid(mosaic.iaro2_a) = 0.0;
    aer_jliquid(mosaic.ialk1_a) = 0.0;
    aer_jliquid(mosaic.iole1_a) = 0.0;
    aer_jliquid(mosaic.iapi1_a) = 0.0;
    aer_jliquid(mosaic.iapi2_a) = 0.0;
    aer_jliquid(mosaic.ilim1_a) = 0.0;
    aer_jliquid(mosaic.ilim2_a) = 0.0;
}

/*
! called when aerosol bin is completely solid.
!
! author: Rahul A. Zaveri
! update: jan 2005
*/

KOKKOS_INLINE_FUNCTION static 
void adjust_solid_aerosol(const Mosaic& mosaic, 
                          const real_type_1d_view_type& aer_jsolid,
                          const real_type_1d_view_type& aer_jliquid,
                          const real_type_1d_view_type& aer_jtotal, 
                          const real_type_1d_view_type& electrolyte_jsolid,
                          const real_type_1d_view_type& electrolyte_jliquid,
                          const real_type_1d_view_type& electrolyte_jtotal,
                          const real_type_1d_view_type& epercent_jsolid,
                          const real_type_1d_view_type& epercent_jliquid,
                          const real_type_1d_view_type& epercent_jtotal,
                          real_type& water_a,real_type&jphase, real_type&jhyst_leg ) {
    // Assuming naer, nelectrolyte are accessible or part of the Mosaic object

    // Since water_a is now an input parameter, we can directly modify it.
    water_a = 0.0;
    jphase = jsolid
    jhyst_leg = jhyst_lo;// lower curve

    // Transfer aer(jtotal) to aer(jsolid) and set aer(jliquid) to 0
    for (int iaer = 0; iaer < mosaic.naer; iaer++) {
        aer_jsolid(iaer) = aer_jtotal(iaer);
        aer_jliquid(iaer) = 0.0;
    }

    // Transfer electrolyte(jtotal) to electrolyte(jsolid) and set electrolyte(jliquid) to 0
    for (int je = 0; je < mosaic.nelectrolyte; je++) {
        electrolyte_jliquid(je) = 0.0;
        epercent_jliquid(je) = 0.0;
        electrolyte_jsolid(je) = electrolyte_jtotal(je);
        epercent_jsolid(je) = epercent_jtotal(je);
    }

    // Update aer(jtotal) that may have been affected above
    aer_jtotal(mosaic.inh4_a) = aer_jsolid(mosaic.inh4_a);
    aer_jtotal(mosaic.ino3_a) = aer_jsolid(mosaic.ino3_a);
    aer_jtotal(mosaic.icl_a) = aer_jsolid(mosaic.icl_a);

    // Update electrolyte(jtotal)
    for (int je = 0; je < mosaic.nelectrolyte; je++) {
        electrolyte_jtotal(je) = electrolyte_jsolid(je);
        epercent_jtotal(je) = epercent_jsolid(je);
    }
}// adjust_solid_aerosol


  /* 
! called when aerosol bin is completely liquid.
!
! author: Rahul A. Zaveri
! update: jan 2005
  */
  KOKKOS_INLINE_FUNCTION static void
void do_full_deliquescence(const Mosaic& mosaic, 
                           const real_type_1d_view_type& electrolyte_jsolid,
                           const real_type_1d_view_type& electrolyte_jliquid,
                           const real_type_1d_view_type& electrolyte_jtotal,
                           const real_type_1d_view_type& aer_jsolid,
                           const real_type_1d_view_type& aer_jliquid, 
                           const real_type_1d_view_type& aer_jtotal) {
    // Partition all electrolytes into liquid phase
    for (int js = 0; js < mosaic.nelectrolyte; js++) {
        electrolyte_jsolid(js) = 0.0;
        electrolyte_jliquid(js) = electrolyte_jtotal(js);
    }

    // Except these electrolytes, which always remain in the solid phase
    electrolyte_jsolid(mosaic.jcaco3) = electrolyte_jtotal(mosaic.jcaco3);
    electrolyte_jsolid(mosaic.jcaso4) = electrolyte_jtotal(mosaic.jcaso4);
    electrolyte_jliquid(mosaic.jcaco3) = 0.0;
    electrolyte_jliquid(mosaic.jcaso4) = 0.0;

    // Partition all the generic aer species into solid and liquid phases
    // Solid phase
    aer_jsolid(mosaic.iso4_a) = electrolyte_jsolid(mosaic.jcaso4);
    aer_jsolid(mosaic.ino3_a) = 0.0;
    aer_jsolid(mosaic.icl_a) = 0.0;
    aer_jsolid(mosaic.inh4_a) = 0.0;
    aer_jsolid(mosaic.ioc_a) = aer_jtotal(mosaic.ioc_a);
    aer_jsolid(mosaic.imsa_a) = 0.0;
    aer_jsolid(mosaic.ico3_a) = aer_jtotal(mosaic.ico3_a);
    aer_jsolid(mosaic.ina_a) = 0.0;
    aer_jsolid(mosaic.ica_a) = electrolyte_jsolid(mosaic.jcaco3) + electrolyte_jsolid(mosaic.jcaso4);
    aer_jsolid(mosaic.ibc_a) = aer_jtotal(mosaic.ibc_a);
    aer_jsolid(mosaic.ioin_a) = aer_jtotal(mosaic.ioin_a);
    aer_jsolid(mosaic.iaro1_a) = aer_jtotal(mosaic.iaro1_a);
    aer_jsolid(mosaic.iaro2_a) = aer_jtotal(mosaic.iaro2_a);
    aer_jsolid(mosaic.ialk1_a) = aer_jtotal(mosaic.ialk1_a);
    aer_jsolid(mosaic.iole1_a) = aer_jtotal(mosaic.iole1_a);
    aer_jsolid(mosaic.iapi1_a) = aer_jtotal(mosaic.iapi1_a);
    aer_jsolid(mosaic.iapi2_a) = aer_jtotal(mosaic.iapi2_a);
    aer_jsolid(mosaic.ilim1_a) = aer_jtotal(mosaic.ilim1_a);
    aer_jsolid(mosaic.ilim2_a) = aer_jtotal(mosaic.ilim2_a);

    // Liquid phase
    aer_jliquid(mosaic.iso4_a) = aer_jtotal(mosaic.iso4_a) - electrolyte_jsolid(mosaic.jcaso4);
    aer_jliquid(mosaic.ino3_a) = aer_jtotal(mosaic.ino3_a);
    aer_jliquid(mosaic.icl_a) = aer_jtotal(mosaic.icl_a);
    aer_jliquid(mosaic.inh4_a) = aer_jtotal(mosaic.inh4_a);
    aer_jliquid(mosaic.ioc_a) = 0.0;
    aer_jliquid(mosaic.imsa_a) = aer_jtotal(mosaic.imsa_a);
    aer_jliquid(mosaic.ico3_a) = 0.0;
    aer_jliquid(mosaic.ina_a) = aer_jtotal(mosaic.ina_a);
    aer_jliquid(mosaic.ica_a) = electrolyte_jtotal(mosaic.jcano3) + electrolyte_jtotal(mosaic.jcacl2);
    aer_jliquid(mosaic.ibc_a) = 0.0;
    aer_jliquid(mosaic.ioin_a) = 0.0;
    aer_jliquid(mosaic.iaro1_a) = 0.0;
    aer_jliquid(mosaic.iaro2_a) = 0.0;
    aer_jliquid(mosaic.ialk1_a) = 0.0;
    aer_jliquid(mosaic.iole1_a) = 0.0;
    aer_jliquid(mosaic.iapi1_a) = 0.0;
    aer_jliquid(mosaic.iapi2_a) = 0.0;
    aer_jliquid(mosaic.ilim1_a) = 0.0;
    aer_jliquid(mosaic.ilim2_a) = 0.0;
}

  KOKKOS_INLINE_FUNCTION static void
  compute_activities()
  {
    
  }


  template<typename MemberType>
  KOKKOS_INLINE_FUNCTION static void team_invoke(
    const MemberType& member,
    const real_type& t,
    const real_type& p,
    const real_type_1d_view_type& aer_total,
    )
    {

    // FIXME: 
    const Mosaic& mosaic{}; 
    // jtotal is 3
    // FIXME: aero uses jtotal: sum of solid and liquid. 
    // we may sum solid and liquid here. 
    // aero can be a view of 1D liquid first and then solid. 
    real_type XT=0.0;
    calculate_XT(aer_total,
                 mosaic, 
                 XT); 

    const real_type CRH = 0.35;

    if (aH2O_a < CRH &&
        (XT > 1.0 || XT < 0.0) &&
        epercent(jcano3) <= ptol_mol_astem &&
        epercent(jcacl2) <= ptol_mol_astem) {
        jaerosolstate = all_solid;
        jphase = jsolid;
        jhyst_leg = jhyst_lo;
        adjust_solid_aerosol(mosaic, 
                     aer_jsolid, aer_jliquid, aer_jtotal, 
                     electrolyte_jsolid, electrolyte_jliquid, electrolyte_jtotal,
                     epercent_jsolid, epercent_jliquid, epercent_jtotal,
                     water_a, jphase, jhyst_leg);

    // step 2: Check for supersaturation/metastable state
    // if(water_a_hyst(ibin)*0.0 .gt. 0.5*water_a_up(ibin))then ! 3-D
    if (jhyst_leg == jhyst_up) { // BOX
      do_full_deliquescence(mosaic, 
                      electrolyte_jsolid, electrolyte_jliquid, electrolyte_jtotal,
                      aer_jsolid, aer_jliquid, aer_jtotal);

      // Assuming electrolyte is a 3D vector or an array
      real_type sum_soluble = 0.0;
      // FIXME nsoluble
      for (int js = 0; js < nsoluble; ++js) {
        sum_soluble += electrolyte_total(js);
      }

      const real_type solids = electrolyte_total(jcaso4) +
                       electrolyte_total(jcaco3) +
                       aer_total(ioin_a);

      if (sum_soluble < 1.e-15 && solids > 0.0) {
        jaerosolstate = all_solid; // no soluble material present
        jphase = jsolid;
        adjust_solid_aerosol(mosaic, 
                     aer_jsolid, aer_jliquid, aer_jtotal, 
                     electrolyte_jsolid, electrolyte_jliquid, electrolyte_jtotal,
                     epercent_jsolid, epercent_jliquid, epercent_jtotal,
                     water_a, jphase, jhyst_leg);

        // New wet mass and wet volume
        mass_wet_a = mass_dry_a + water_a * 1.e-3; // g/cc(air)
        vol_wet_a = vol_dry_a + water_a * 1.e-3; // cc(aer)/cc(air) or m^3/m^3(air)
        growth_factor = mass_wet_a / mass_dry_a; // mass growth factor

        return;
    } else if (sum_soluble > 0.0 && solids == 0.0) {
        jaerosolstate = all_liquid;
        jhyst_leg = jhyst_up;
        jphase = jliquid;
        // FIXME: aerosol_water_total is an input
        water_a = aerosol_water_total;

        if (water_a < 0.0) {
            jaerosolstate = all_solid; // no soluble material present
            jphase = jsolid;
            jhyst_leg = jhyst_lo;
            adjust_solid_aerosol(mosaic, 
                     aer_jsolid, aer_jliquid, aer_jtotal, 
                     electrolyte_jsolid, electrolyte_jliquid, electrolyte_jtotal,
                     epercent_jsolid, epercent_jliquid, epercent_jtotal,
                     water_a, jphase, jhyst_leg);
        } else {
            adjust_liquid_aerosol(mosaic,
                      electrolyte_jsolid, electrolyte_jliquid, electrolyte_jtotal,
                      epercent_jsolid, epercent_jliquid, epercent_jtotal,
                      aer_jsolid, aer_jliquid, aer_jtotal,
                      jphase, jhyst_leg);
            compute_activities();
        }

        // New wet mass and wet volume
        mass_wet_a = mass_dry_a + water_a * 1.e-3; // g/cc(air)
        vol_wet_a = vol_dry_a + water_a * 1.e-3; // cc(aer)/cc(air) or m^3/m^3(air)
        growth_factor = mass_wet_a / mass_dry_a; // mass growth factor

        return;
    }
    }// end step 2


    }


    }
};    	
} // namespace Impl
} // namespace TChem

#endif
