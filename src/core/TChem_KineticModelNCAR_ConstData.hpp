/* =====================================================================================
TChem-atm version 1.0
Copyright (2024) NTESS
https://github.com/sandialabs/TChem-atm

Copyright 2024 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains
certain rights in this software.

This file is part of TChem-atm. TChem-atm is open source software: you can redistribute it
and/or modify it under the terms of BSD 2-Clause License
(https://opensource.org/licenses/BSD-2-Clause). A copy of the licese is also
provided under the main directory

Questions? Contact Oscar Diaz-Ibarra at <odiazib@sandia.gov>, or
           Mike Schmidt at <mjschm@sandia.gov>, or
           Cosmin Safta at <csafta@sandia.gov>

Sandia National Laboratories, New Mexico/Livermore, NM/CA, USA
===================================================================================== */
#ifndef __TCHEM_KINETIC_MODEL_NCAR_CONST_MDATA_HPP__
#define __TCHEM_KINETIC_MODEL_NCAR_CONST_MDATA_HPP__

#include "TChem_Util.hpp"
#include "TChem_KineticModelData.hpp"

namespace TChem {

  template<typename DeviceType>
  struct KineticModelNCAR_ConstData {
  public:
    using device_type = DeviceType;

    /// non const
    using real_type_0d_view_type = Tines::value_type_0d_view<real_type,device_type>;
    using real_type_1d_view_type = Tines::value_type_1d_view<real_type,device_type>;
    using real_type_2d_view_type = Tines::value_type_2d_view<real_type,device_type>;
    using real_type_3d_view_type = Tines::value_type_3d_view<real_type,device_type>;

    using ordinal_type_1d_dual_view = Tines::value_type_1d_dual_view<ordinal_type,device_type>;
    using ordinal_type_2d_dual_view = Tines::value_type_2d_dual_view<ordinal_type,device_type>;

    using ordinal_type_1d_view = Tines::value_type_1d_view<ordinal_type,device_type>;
    using ordinal_type_2d_view = Tines::value_type_2d_view<ordinal_type,device_type>;

    using string_type_1d_dual_view_type = Tines::value_type_1d_dual_view<char [LENGTHOFSPECNAME + 1],device_type>;
    using string_type_1d_view_type = Tines::value_type_1d_view<char [LENGTHOFSPECNAME + 1],device_type>;

    using arrhenius_reaction_type_1d_view_type =  Tines::value_type_1d_view<ArrheniusReactionType,device_type>;

    using cmaq_h2o2_reaction_type_1d_view_type =  Tines::value_type_1d_view<CMAQ_H2O2ReactionType,device_type>;

    using jpl_reaction_type_1d_view_type =  Tines::value_type_1d_view<JPL_ReactionType,device_type>;

    using r_jpl_arrhenius_reaction_type_1d_view_type =  Tines::value_type_1d_view<R_JPL_ArrheniusReactionType,device_type>;

    using emission_source_type_1d_view_type =  Tines::value_type_1d_view<EMISSION_SourceType,device_type>;

    using adjust_reaction_type_1d_view_type =  Tines::value_type_1d_view<adjust_reactionType,device_type>;

    using prod_O1D_type_1d_view_type =  Tines::value_type_1d_view<prod_O1DType,device_type>;

    //// const views
    using kmcd_ordinal_type_1d_view = ConstUnmanaged<ordinal_type_1d_view>;
    using kmcd_ordinal_type_2d_view = ConstUnmanaged<ordinal_type_2d_view>;

    using kmcd_real_type_1d_view = ConstUnmanaged<real_type_1d_view_type>;
    using kmcd_real_type_2d_view = ConstUnmanaged<real_type_2d_view_type>;
    using kmcd_real_type_3d_view = ConstUnmanaged<real_type_3d_view_type>;

    using kmcd_string_type_1d_view =ConstUnmanaged<string_type_1d_view_type>;

    using kmcd_arrhenius_reaction_type_1d_view = ConstUnmanaged<arrhenius_reaction_type_1d_view_type>;
    using kmcd_cmaq_h2o2_reaction_type_1d_view= ConstUnmanaged<cmaq_h2o2_reaction_type_1d_view_type>;

    // e3sm implementation of troe reaction
    using kmcd_jpl_reaction_type_1d_view_1d_view = ConstUnmanaged<jpl_reaction_type_1d_view_type>;

    using kmcd_r_jpl_arrhenius_reaction_type_1d_view = ConstUnmanaged<r_jpl_arrhenius_reaction_type_1d_view_type>;

    using kmcd_adjust_reaction_type_1d_view = ConstUnmanaged<adjust_reaction_type_1d_view_type>;
    using kmcd_prod_O1D_type_1d_view = ConstUnmanaged<prod_O1D_type_1d_view_type>;


    using kmcd_emission_source_type_1d_view_type  = ConstUnmanaged<emission_source_type_1d_view_type>;

    kmcd_string_type_1d_view speciesNames;
    kmcd_ordinal_type_1d_view reacNreac;
    kmcd_ordinal_type_1d_view reacNprod;
    kmcd_real_type_2d_view reacArhenFor;
    kmcd_real_type_2d_view reacNuki;
    kmcd_ordinal_type_2d_view reacSidx;
    kmcd_arrhenius_reaction_type_1d_view ArrheniusCoef;
    kmcd_cmaq_h2o2_reaction_type_1d_view CMAQ_H2O2Coef;
    kmcd_jpl_reaction_type_1d_view_1d_view JPL_Coef;
    kmcd_r_jpl_arrhenius_reaction_type_1d_view R_JPL_ArrheniusCoef;
    kmcd_adjust_reaction_type_1d_view adjust_reaction;
    kmcd_prod_O1D_type_1d_view prod_O1D;
    ordinal_type number_of_prod_O1D;

    kmcd_emission_source_type_1d_view_type EmissionCoef;
    ordinal_type M_index;


    ordinal_type nSpec;
    ordinal_type nReac;
    kmcd_ordinal_type_1d_view reacPfal;
    kmcd_real_type_2d_view reacPpar;
    ordinal_type nConstSpec;
    real_type CONV_PPM;
    kmcd_ordinal_type_1d_view convExponent;
  };

  template<typename SpT>
  KineticModelNCAR_ConstData<SpT> createNCAR_KineticModelConstData(const KineticModelData & kmd) {
    KineticModelNCAR_ConstData<SpT> data;
    // forward arrhenius coefficients
    data.speciesNames = kmd.sNames_.template view<SpT>();
    data.reacNuki = kmd.reacNuki_.template view<SpT>();
    data.reacSidx = kmd.reacSidx_.template view<SpT>();
    data.reacArhenFor = kmd.reacArhenFor_.template view<SpT>();
    data.reacNreac = kmd.reacNreac_.template view<SpT>();
    data.reacNprod = kmd.reacNprod_.template view<SpT>();
    data.nSpec = kmd.nSpec_;
    data.nReac = kmd.nReac_;
    // troe type
    data.reacPpar = kmd.reacPpar_.template view<SpT>();
    data.reacPfal = kmd.reacPfal_.template view<SpT>();
    // arrhenius pressure parameters
    data.ArrheniusCoef = kmd.ArrheniusCoef_.template view<SpT>();
    // cmaq h2o2 special type
    data.CMAQ_H2O2Coef = kmd.CMAQ_H2O2Coef_.template view<SpT>();
    data.EmissionCoef = kmd.EmissionCoef_.template view<SpT>();
    // species that are assumed constant like tracers
    data.nConstSpec = kmd.nConstSpec_;
    data.CONV_PPM = kmd.CONV_PPM_;
    data.convExponent = kmd.convExponent_.template view<SpT>();
    //
    data.JPL_Coef = kmd.JPL_Coef_.template view<SpT>();
    data.R_JPL_ArrheniusCoef = kmd.R_JPL_ArrheniusCoef_.template view<SpT>();
    data.adjust_reaction =kmd.adjust_reaction_.template view<SpT>();
    data.prod_O1D = kmd.prod_O1D_.template view<SpT>();
    data.number_of_prod_O1D = kmd.number_of_prod_O1D_;
    data.M_index =kmd.M_index_;

    return data;
  }

  template<typename DT>
  static inline
  Kokkos::View<KineticModelNCAR_ConstData<DT>*,DT>
  createNCAR_KineticModelConstData(const kmd_type_1d_view_host kmds) {
    Kokkos::View<KineticModelNCAR_ConstData<DT>*,DT>
      r_val(do_not_init_tag("KMCD::NCAR const data objects"),
	    kmds.extent(0));
    auto r_val_host = Kokkos::create_mirror_view(r_val);
    Kokkos::parallel_for
      (Kokkos::RangePolicy<host_exec_space>(0, kmds.extent(0)),
       [=](const int i) {
	r_val_host(i) = createNCAR_KineticModelConstData<DT>(kmds(i));
      });
    Kokkos::deep_copy(r_val, r_val_host);
    return r_val;
  }

} // namespace TChem
#endif
