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
! called when aerosol bin is completely solid.
!
! author: Rahul A. Zaveri
! update: jan 2005
*/
KOKKOS_INLINE_FUNCTION static 
void adjust_solid_aerosol(
                          // int jsolid,
                          // int jhyst_lo,
                          const real_type_1d_view_type& aer_total,
                          const real_type_1d_view_type& aer_solid,
                          real_type_1d_view_type& aer_liquid,
                          
                          const real_type_1d_view_type& electrolyte_total,
                          const real_type_1d_view_type& electrolyte_solid,
                          real_type_1d_view_type& electrolyte_liquid,
                          const real_type_1d_view_type& epercent_total,
                          const real_type_1d_view_type& epercent_solid,
                          const Mosaic& mosaic
                          )
{
    // Set phase and hysteresis leg to solid
    jphase = jsolid;
    jhyst_leg = jhyst_lo; // lower curve
    water_a = 0.0;

    int naer = mosaic.naer;
    int nelectrolyte = mosaic.nelectrolyte;
    int inh4_a = mosaic.inh4_a; 
    int ino3_a = mosaic.ino3_a;
    int icl_a = mosaic.icl_a;

    // Transfer aer(jtotal) to aer(jsolid) and set aer(jliquid) to 0.0
    for (int iaer = 0; iaer < naer; ++iaer) {
        aer_solid(iaer) = aer_total(iaer);
        aer_liquid(iaer) = 0.0;
    }

    // Transfer electrolyte(jtotal) to electrolyte(jsolid) and set electrolyte(jliquid) and epercent(jliquid) to 0.0
    for (int je = 0; je < nelectrolyte; ++je) {
        electrolyte_liquid(je) = 0.0;
        epercent_liquid(je) = 0.0;
        electrolyte_solid(je) = electrolyte_total(je);
        epercent_solid(je) = epercent_total(je);
    }

    // Update aer(jtotal) that may have been affected above
    aer_total(inh4_a) = aer_solid(inh4_a);
    aer_total(ino3_a) = aer_solid(ino3_a);
    aer_total(icl_a) = aer_solid(icl_a);

    // Update electrolyte(jtotal)
    for (int je = 0; je < nelectrolyte; ++je) {
        electrolyte_total(je) = electrolyte_solid(je);
        epercent_total(je) = epercent_solid(je);
    }
}

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
  KOKKOS_INLINE_FUNCTION static void
  do_full_deliquescence()
  {

  }
  /* 
! called when aerosol bin is completely liquid.
!
! author: Rahul A. Zaveri
! update: jan 2005
  */
  KOKKOS_INLINE_FUNCTION static void
  adjust_liquid_aerosol()
  {

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
        adjust_solid_aerosol(
    // jsolid,
    // jhyst_lo,
    aer_total,
    aer_solid,
    aer_liquid,
    electrolyte_total,
    electrolyte_solid,
    electrolyte_liquid,
    epercent_total,
    epercent_solid,
    mosaic);


    // step 2: Check for supersaturation/metastable state
    // if(water_a_hyst(ibin)*0.0 .gt. 0.5*water_a_up(ibin))then ! 3-D
    if (jhyst_leg == jhyst_up) { // BOX
      do_full_deliquescence();

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
        adjust_solid_aerosol(
    // jsolid,
    // jhyst_lo,
    aer_total,
    aer_solid,
    aer_liquid,
    electrolyte_total,
    electrolyte_solid,
    electrolyte_liquid,
    epercent_total,
    epercent_solid,
    mosaic);

        // New wet mass and wet volume
        mass_wet_a = mass_dry_a + water_a * 1.e-3; // g/cc(air)
        vol_wet_a = vol_dry_a + water_a * 1.e-3; // cc(aer)/cc(air) or m^3/m^3(air)
        growth_factor = mass_wet_a / mass_dry_a; // mass growth factor

        return;
    } else if (sum_soluble > 0.0 && solids == 0.0) {
        jaerosolstate = all_liquid;
        jhyst_leg = jhyst_up;
        jphase = jliquid;
        water_a = aerosol_water_total();

        if (water_a < 0.0) {
            jaerosolstate = all_solid; // no soluble material present
            jphase = jsolid;
            jhyst_leg = jhyst_lo;
            adjust_solid_aerosol(
    // jsolid,
    // jhyst_lo,
    aer_total,
    aer_solid,
    aer_liquid,
    electrolyte_total,
    electrolyte_solid,
    electrolyte_liquid,
    epercent_total,
    epercent_solid,
    mosaic);
        } else {
            adjust_liquid_aerosol();
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
