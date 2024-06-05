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

    const ordinal_type naer=19;
    const ordinal_type nelectrolyte=22;
    const ordinal_type nsalt=15;
    const ordinal_type nsoluble=20;
    const ordinal_type ncation=4;
    const ordinal_type nanion=5;

    // FIXME: update values and make enum
    const ordinal_type jsolid=1;
    const ordinal_type jliquid=2;
    const ordinal_type jtotal=3;

    const ordinal_type no_aerosol=0;
    const ordinal_type all_solid=1;
    const ordinal_type all_liquid=2;
    const ordinal_type mixed=3;

    const ordinal_type soluble=1;
    const ordinal_type insoluble=2;

    const ordinal_type jhyst_up=1;
    const ordinal_type jhyst_lo=0;

    const ordinal_type nsoluble=20;

    // polynomial coefficients for binary molality (used in ZSR)
       // for aw < 0.97
    const real_type_2d_view_type a_zsr("a_zsr", 6, nelectrolyte);
    const real_type_1d_view_type aw_min("aw_min", nelectrolyte);

    a_zsr(0, jnh4hso4) =   1.30894;
    a_zsr(1, jnh4hso4) =  -7.09922;
    a_zsr(2, jnh4hso4) =  20.62831;
    a_zsr(3, jnh4hso4) = -32.19965;
    a_zsr(4, jnh4hso4) =  25.17026;
    a_zsr(5, jnh4hso4) =  -7.81632;
    aw_min(jnh4hso4)   =       0.1;

    a_zsr(0, jlvcite)  =   1.10725;
    a_zsr(1, jlvcite)  =  -5.17978;
    a_zsr(2, jlvcite)  =  12.29534;
    a_zsr(3, jlvcite)  = -16.32545;
    a_zsr(4, jlvcite)  =  11.29274;
    a_zsr(5, jlvcite)  =  -3.19164;
    aw_min(jlvcite)    =       0.1;

    a_zsr(0, jnh4hso4) =  1.15510;
    a_zsr(1, jnh4hso4) = -3.20815;
    a_zsr(2, jnh4hso4) =  2.71141;
    a_zsr(3, jnh4hso4) =  2.01155;
    a_zsr(4, jnh4hso4) = -4.71014;
    a_zsr(5, jnh4hso4) =  2.04616;
    aw_min(jnh4hso4)   =      0.1;

    a_zsr(0, jnh4msa) =  1.15510;
    a_zsr(1, jnh4msa) = -3.20815;
    a_zsr(2, jnh4msa) =  2.71141;
    a_zsr(3, jnh4msa) =  2.01155;
    a_zsr(4, jnh4msa) = -4.71014;
    a_zsr(5, jnh4msa) =  2.04616;
    aw_min(jnh4msa)   =      0.1;

    a_zsr(0, jnh4no3) =   0.43507;
    a_zsr(1, jnh4no3) =   6.38220;
    a_zsr(2, jnh4no3) = -30.19797;
    a_zsr(3, jnh4no3) =  53.36470;
    a_zsr(4, jnh4no3) = -43.44203;
    a_zsr(5, jnh4no3) =  13.46158;
    aw_min(jnh4no3)   =       0.1;

    a_zsr(0, jnh4cl) =  0.45309;
    a_zsr(1, jnh4cl) =  2.65606;
    a_zsr(2, jnh4cl) = -14.7730;
    a_zsr(3, jnh4cl) =  26.2936;
    a_zsr(4, jnh4cl) = -20.5735;
    a_zsr(5, jnh4cl) =  5.94255;
    aw_min(jnh4cl)   =      0.1;

    a_zsr(0, jnacl) =  0.42922;
    a_zsr(1, jnacl) = -1.17718;
    a_zsr(2, jnacl) =  2.80208;
    a_zsr(3, jnacl) = -4.51097;
    a_zsr(4, jnacl) =  3.76963;
    a_zsr(5, jnacl) = -1.31359;
    aw_min(jnacl)   =      0.1;

    a_zsr(0, jnano3) =   1.34966;
    a_zsr(1, jnano3) =  -5.20116;
    a_zsr(2, jnano3) =  11.49011;
    a_zsr(3, jnano3) = -14.41380;
    a_zsr(4, jnano3) =   9.07037;
    a_zsr(5, jnano3) =  -2.29769;
    aw_min(jnano3)   =       0.1;

    a_zsr(0, jna2so4) =  0.39888;
    a_zsr(1, jna2so4) = -1.27150;
    a_zsr(2, jna2so4) =  3.42792;
    a_zsr(3, jna2so4) = -5.92632;
    a_zsr(4, jna2so4) =  5.33351;
    a_zsr(5, jna2so4) = -1.96541;
    aw_min(jna2so4)   =      0.1;

    a_zsr(0, jna3hso4) =  0.31480;
    a_zsr(1, jna3hso4) = -1.01087;
    a_zsr(2, jna3hso4) =  2.44029;
    a_zsr(3, jna3hso4) = -3.66095;
    a_zsr(4, jna3hso4) =  2.77632;
    a_zsr(5, jna3hso4) = -0.86058;
    aw_min(jna3hso4)   =      0.1;

    a_zsr(0, jnahso4)  =   0.62764;
    a_zsr(1, jnahso4)  =  -1.63520;
    a_zsr(2, jnahso4)  =   4.62531;
    a_zsr(3, jnahso4)  = -10.06925;
    a_zsr(4, jnahso4)  =  10.33547;
    a_zsr(5, jnahso4)  =  -3.88729;
    aw_min(jnahso4)    =       0.1;

    a_zsr(0, jnamsa)  =   0.62764;
    a_zsr(1, jnamsa)  =  -1.63520;
    a_zsr(2, jnamsa)  =   4.62531;
    a_zsr(3, jnamsa)  = -10.06925;
    a_zsr(4, jnamsa)  =  10.33547;
    a_zsr(5, jnamsa)  =  -3.88729;
    aw_min(jnamsa)    =       0.1;

    a_zsr(0, jcano3) =  0.38895;
    a_zsr(1, jcano3) = -1.16013;
    a_zsr(2, jcano3) =  2.16819;
    a_zsr(3, jcano3) = -2.23079;
    a_zsr(4, jcano3) =  1.00268;
    a_zsr(5, jcano3) = -0.16923;
    aw_min(jcano3)   =      0.1;

    a_zsr(0, jcacl2) =  0.29891;
    a_zsr(1, jcacl2) = -1.31104;
    a_zsr(2, jcacl2) =  3.68759;
    a_zsr(3, jcacl2) = -5.81708;
    a_zsr(4, jcacl2) =  4.67520;
    a_zsr(5, jcacl2) = -1.53223;
    aw_min(jcacl2)   =      0.1;

    a_zsr(0, jh2so4) =  0.32751;
    a_zsr(1, jh2so4) = -1.00692;
    a_zsr(2, jh2so4) =  2.59750;
    a_zsr(3, jh2so4) = -4.40014;
    a_zsr(4, jh2so4) =  3.88212;
    a_zsr(5, jh2so4) = -1.39916;
    aw_min(jh2so4)   =      0.1;

    a_zsr(0, jmsa) =  0.32751;
    a_zsr(1, jmsa) = -1.00692;
    a_zsr(2, jmsa) =  2.59750;
    a_zsr(3, jmsa) = -4.40014;
    a_zsr(4, jmsa) =  3.88212;
    a_zsr(5, jmsa) = -1.39916;
    aw_min(jmsa)   =      0.1;

    a_zsr(0, jhhso4) =  0.32751;
    a_zsr(1, jhhso4) = -1.00692;
    a_zsr(2, jhhso4) =  2.59750;
    a_zsr(3, jhhso4) = -4.40014;
    a_zsr(4, jhhso4) =  3.88212;
    a_zsr(5, jhhso4) = -1.39916;
    aw_min(jhhso4)   =     1.0;

    a_zsr(0, jhno3) =   0.75876;
    a_zsr(1, jhno3) =  -3.31529;
    a_zsr(2, jhno3) =   9.26392;
    a_zsr(3, jhno3) = -14.89799;
    a_zsr(4, jhno3) =  12.08781;
    a_zsr(5, jhno3) =  -3.89958;
    aw_min(jhno3)   =       0.1;

    a_zsr(0, jhcl) =  0.31133;
    a_zsr(1, jhcl) = -0.79688;
    a_zsr(2, jhcl) =  1.93995;
    a_zsr(3, jhcl) = -3.31582;
    a_zsr(4, jhcl) =  2.93513;
    a_zsr(5, jhcl) = -1.07268;
    aw_min(jhcl)   =      0.1;

    a_zsr(0, jcaso4)  =  0.0;
    a_zsr(1, jcaso4)  =  0.0;
    a_zsr(2, jcaso4)  =  0.0;
    a_zsr(3, jcaso4)  =  0.0;
    a_zsr(4, jcaso4)  =  0.0;
    a_zsr(5, jcaso4)  =  0.0;
    aw_min(jcaso4)    =  1.0;

    a_zsr(0, jcamsa2) =  0.38895;
    a_zsr(1, jcamsa2) = -1.16013;
    a_zsr(2, jcamsa2) =  2.16819;
    a_zsr(3, jcamsa2) = -2.23079;
    a_zsr(4, jcamsa2) =  1.00268;
    a_zsr(5, jcamsa2) = -0.16923;
    aw_min(jcamsa2)   =      0.1;

    a_zsr(0, jcaco3)  =  0.0;
    a_zsr(1, jcaco3)  =  0.0;
    a_zsr(2, jcaco3)  =  0.0;
    a_zsr(3, jcaco3)  =  0.0;
    a_zsr(4, jcaco3)  =  0.0;
    a_zsr(5, jcaco3)  =  0.0;
    aw_min(jcaco3)    =  1.0;

       // for aw => 0.97 to 0.999999
    const real_type_1d_view_type b_zsr("b_zsr", nelectrolyte);

    b_zsr(jnh4so4)  = 28.0811;
    b_zsr(jlvcite)  = 14.7178;
    b_zsr(jnh4hso4) = 29.4779;
    b_zsr(jnh4msa)  = 29.4779;
    b_zsr(jnh4no3)  = 33.4049;
    b_zsr(jnh4cl)   = 30.8888;
    b_zsr(jnacl)    = 29.8375;
    b_zsr(jnano3)   = 32.2756;
    b_zsr(jna2so4)  = 27.6889;
    b_zsr(jna3hso4) = 14.2184;
    b_zsr(jnahso4)  = 28.3367;
    b_zsr(jnamsa)   = 28.3367;
    b_zsr(jcano3)   = 18.3661;
    b_zsr(jcacl2)   = 20.8792;
    b_zsr(jh2so4)   = 26.7347;
    b_zsr(jhhso4)   = 26.7347;
    b_zsr(jhno3)    = 28.8257;
    b_zsr(jhcl)     = 27.7108;
    b_zsr(jmsa)     = 26.7347;
    b_zsr(jcaso4)   = 0.0;
    b_zsr(jcamsa2)  = 18.3661;
    b_zsr(jcaco3)   = 0.0;

    // local aerosol ions
    // cations
    const ordinal_type jc_h   = 1;
    const ordinal_type jc_nh4 = 2;
    const ordinal_type jc_na  = 3;
    const ordinal_type jc_ca  = 4;

    // anions
    const ordinal_type ja_hso4 = 1;
    const ordinal_type ja_so4  = 2;
    const ordinal_type ja_no3  = 3;
    const ordinal_type ja_cl   = 4;
    const ordinal_type ja_msa  = 5;
    // const ordinal_type ja_co3  = 6; (note: this is commented out in MOSAIC as well)

    // Temperature-dependent thermodynamic parameters
    // liquid-liquid
    const ordinal_type nrxn_aer_ll = 3;
    const real_type_1d_view_type Keq_298_ll("Keq_298_ll", nrxn_aer_ll);
    Keq_298_ll(0) = 1.0502e-2; // HSO4- <=> SO4= + H+
    Keq_298_ll(1) = 1.805e-5;  // NH3(l) + H2O = NH4+ + OH-
    Keq_298_ll(2) = 1.01e-14;  // H2O(l) <=> H+ + OH-
    const real_type_1d_view_type Keq_a_ll("Keq_a_ll", nrxn_aer_ll);
    Keq_a_ll(0) =   8.85;
    Keq_a_ll(1) =  -1.50;
    Keq_a_ll(2) = -22.52;
    const real_type_1d_view_type Keq_b_ll("Keq_b_ll", nrxn_aer_ll);
    Keq_b_ll(0) = 25.14;
    Keq_b_ll(1) = 26.92;
    Keq_b_ll(2) = 26.92;


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

  KOKKOS_INLINE_FUNCTION static
  void compute_activities(const MosaicModelData& mosaic,
                          const real_type_1d_view_type& aer_liquid,
                          const real_type_1d_view_type& ma,
                          const real_type_1d_view_type& mc,
                          const real_type_1d_view_type& electrolyte_solid,
                          const real_type_1d_view_type& electrolyte_liquid,
                          const real_type_1d_view_type& electrolyte_total,
                          const real_type_1d_view_type& xmol,
                          const real_type_1d_view_type& log_gam,
                          const real_type_2d_view_type& log_gamZ,
                          const real_type_1d_view_type& gam,
                          const real_type_1d_view_type& activity,
                          real_type T_K) {

    // get aerosol water activity
    real_type water_a;
    // FIXME: these need to be set beforehand
    ordinal_type jaerosolstate, jphase, jhyst_leg;
    aerosol_water(mosaic, electrolyte_liquid, water_a, jaerosolstate, jphase, jhyst_leg);

    if (water_a == 0.0) {
      return;
    }

    // get sulfate ratio to determine regime
    real_type XT = 0.0;
    calcuate_XT(aer_liquid, mosaic, XT);

    if (XT > 2.0 || XT < 0.0) {
      // SULFATE POOR: fully dissociated electrolytes
      real_type a_c, Keq_ll;

      // anion molalities (mol / kg water)
      ma(mosaic.ja_so4)  = 1.e-9 * aer_liquid(mosaic.iso4_a) / water_a;
      ma(mosaic.ja_hso4) = 0.0;
      ma(mosaic.ja_no3)  = 1.e-9 * aer_liquid(mosaic.ino3_a) / water_a;
      ma(mosaic.ja_cl)   = 1.e-9 * aer_liquid(mosaic.icl_a)  / water_a;
      ma(mosaic.ja_msa)  = 1.e-9 * aer_liquid(mosaic.imsa_a) / water_a;

      // cation molalities (mol / kg water)
      mc(mosaic.jc_ca)  = 1.e-9 * aer_liquid(mosaic.ica_a)  / water_a;
      mc(mosaic.jc_nh4) = 1.e-9 * aer_liquid(mosaic.inh4_a) / water_a;
      mc(mosaic.jc_na)  = 1.e-9 * aer_liquid(mosaic.na_a)   / water_a;
      a_c               = (
                          (2. * ma(mosaic.ja_so4)  +
                                ma(mosaic.ja_no3)  +
                                ma(mosaic.ja_cl)   +
                                ma(mosaic.ja_msa)) -
                          (2. * mc(mosaic.jc_ca)   +
                                mc(mosaic.jc_nh4)  +
                                mc(mosaic.jc_na))  );
      // FIXME: consider adding update_thermodynamic_constants instead 
      // of computing equil. constants directly 
      fn_Keq(mosaic.Keq_ll_298(2), mosaic.Keq_a_ll(2), mosiac.Keq_b_ll(2), T_K, Keq_ll);
      mc(mosaic.jc_h)   = 0.5 * ( (a_c) +
                                  (sqrt(a_c*a_c + 4. * Keq_ll)) ); 
    }
  }

  KOKKOS_INLINE_FUNCTION static
  void bin_molality(const MosaicModelData& mosaic,
                    const ordinal_type& je,
                    real_type& bin_molality) {

    real_type aw, xm;
    // FIXME: aH20_a should be set to the relative humidity
    real_type aH2O_a = 0.0;

    aw = max(aH20_a, mosaic.aw_min(je));
    aw = min(aw, 0.999999);

    if (aw < 0.97) {
      xm = mosaic.a_zsr(0, je) +
           aw * (mosaic.a_zsr(1, je) +
           aw * (mosaic.a_zsr(2, je) +
           aw * (mosaic.a_zsr(3, je) +
           aw * (mosaic.a_zsr(4, je) +
           aw *  mosaic.a_zsr(5, je) ))));
      bin_molality = 55.509 * xm / (1.0 - xm);
    } else {
      bin_molality = -1.0 * mosaic.b_zsr(je) * log(aw);
    }
  }

  // ZSR method (water uptake)
  KOKKOS_INLINE_FUNCTION static
  void aerosol_water(const MosaicModelData& mosaic,
                           const real_type_1d_view_type& electrolyte,
                           real_type& aerosol_water,
                           ordinal_type& jaerosolstate,
                           ordinal_type& jphase,
                           ordinal_type& jhyst_leg) {

    real_type_1d_view_type molality0("molality_0", mosaic.nelectrolyte);
    real_type bin_molality_0;
    for (ordinal_type je = 0; je < mosaic.nelectrolyte; je++) {
      bin_molality(mosaic, je, bin_molality_0);
      molality0(je) = bin_molality_0;
    }

    real_type dum = 0.0;
    for (ordinal_type je = 0; je < (mosaic.nsalt + 4), je++) {
      dum = dum + electrolyte_total(je) / molality0(je);
    }

    aerosol_water = dum * 1.e-9;
    if (aerosol_water <= 0.0) {
      jaerosolstate = mosaic.all_solid;
      jphase = mosaic.jsolid;
      jhyst_leg = mosaic.jhyst_lo;
    }
  }

  KOKKOS_INLINE_FUNCTION static
  void fn_Keq(const real_type Keq_298,
              const real_type a,
              const real_type b,
              const real_type& T,
              real_type& Keq) {

    real_type tt;
    tt = 298.15 / T;

    Keq = Keq_298 * exp(a * (tt - 1.0) + b * (1.0 + log(tt) - tt));
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
