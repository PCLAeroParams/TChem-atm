#ifndef __TCHEM_KINETIC_AEROSOL_MODEL_DATA_HPP__
#define __TCHEM_KINETIC_AEROSOL_MODEL_DATA__HPP__

#include "TChem_Util.hpp"

#include "yaml-cpp/yaml.h"
#include <iostream>

namespace TChem {
  struct AerosolModelData {
  public:

  // number of species per aerosol
  ordinal_type nSpec_;
  ordinal_type nSpec_gas_;
  // number of particles
  ordinal_type nParticles_;
  // aerosol molecular weigths and density
  real_type_1d_dual_view molecular_weigths_, aerosol_density_;

  public:

    AerosolModelData(const std::string& mechfile, std::ostream& echofile);
    AerosolModelData(const std::string& mechfile);

      /// constructor and destructor
  AerosolModelData() = default;
  AerosolModelData(const AerosolModelData& b) = default;
  ~AerosolModelData() = default;

  void initFile(const std::string &mechfile, std::ostream& echofile);
  ordinal_type initChem(YAML::Node& doc, std::ostream& echofile);
  };

  template<typename DeviceType>
  struct AerosolModel_ConstData {
  public:
    using device_type = DeviceType;

    using real_type_1d_view_type = Tines::value_type_1d_view<real_type,device_type>;

    using amcd_real_type_1d_view = ConstUnmanaged<real_type_1d_view_type>;

    ordinal_type GAS_SPEC_;
    ordinal_type AERO_SPEC_i_phase;
    real_type DIFF_COEFF_;
    real_type MW_;
    //
    real_type N_star;
    bool compute_alpha;
    real_type B1_;
    real_type B2_;
    real_type B3_;
    real_type B4_;

    ordinal_type nSpec;
    ordinal_type nSpec_gas;
    ordinal_type nParticles;

    amcd_real_type_1d_view molecular_weigths;
    amcd_real_type_1d_view aerosol_density;

   };

   template<typename SpT>
  AerosolModel_ConstData<SpT> create_AerosolModelConstData(const AerosolModelData & amd) {
    AerosolModel_ConstData<SpT> data;
        // FIXME: input parameters
    data.GAS_SPEC_=0;
    data.AERO_SPEC_i_phase=1;
    // data.num_of_phases =1;
        // from gas species
    data.DIFF_COEFF_=0.95E-05; // m2 s-1
    data.MW_=0.04607; //
    // FIXME: inputs
    data.N_star=2.55; //only gas species
    // reaction info
    data.compute_alpha=true;
    data.B1_= -1.97E+03;
    data.B2_=2.91E+00;
    data.B3_=1.96E-03;
    data.B4_=-4.96E-01;
    //
    data.nSpec_gas=amd.nSpec_gas_;
    data.nSpec=amd.nSpec_;
    data.nParticles=amd.nParticles_;
    data.molecular_weigths = amd.molecular_weigths_.template view<SpT>();
    data.aerosol_density = amd.aerosol_density_.template view<SpT>();

  return data;
  }

} // namespace TChem
#endif
