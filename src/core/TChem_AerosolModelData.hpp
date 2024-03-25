#ifndef __TCHEM_AEROSOL_MODEL_DATA_HPP__
#define __TCHEM_AEROSOL_MODEL_DATA_HPP__

#include "TChem_Util.hpp"
#include "TChem_ReactionTypes.hpp"
#include "yaml-cpp/yaml.h"
#include <iostream>

namespace TChem {
  struct AerosolModelData {
  public:

  // number of species per aerosol
  ordinal_type nSpec_;
  ordinal_type nSpec_gas_;
  ordinal_type nConstSpec_gas_;
  // number of particles
  ordinal_type nParticles_;
  // aerosol molecular weigths and density
  real_type_1d_dual_view molecular_weigths_, aerosol_density_;
  simplo_phase_transfer_type_1d_dual_view simpol_params_;

  // only use this map in host
  std::map<std::string, int> aerosol_sp_name_idx_;

  public:

    AerosolModelData(const std::string& mechfile, std::ostream& echofile);
    AerosolModelData(const std::string& mechfile);

      /// constructor and destructor
  AerosolModelData() = default;
  AerosolModelData(const AerosolModelData& b) = default;
  ~AerosolModelData() = default;

  void initFile(const std::string &mechfile, std::ostream& echofile);
  ordinal_type initChem(YAML::Node& doc, std::ostream& echofile);
  void scenarioConditionParticles(const std::string &mechfile,
                                  const ordinal_type nBatch,
                                  real_type_2d_view_host& num_concentration,
                                  real_type_2d_view_host& state);
  };

  template<typename DeviceType>
  struct AerosolModel_ConstData {
  public:
    using device_type = DeviceType;

    using real_type_1d_view_type = Tines::value_type_1d_view<real_type,device_type>;

    using amcd_real_type_1d_view = ConstUnmanaged<real_type_1d_view_type>;

    using simplo_phase_transfer_type_1d_view_type =  Tines::value_type_1d_view<SIMPOL_PhaseTransferType,device_type>;
    using amcd_simplo_phase_transfer_type_1d_view = ConstUnmanaged<simplo_phase_transfer_type_1d_view_type>;
    amcd_simplo_phase_transfer_type_1d_view simpol_params;

    ordinal_type nSpec;
    ordinal_type nSpec_gas;
    ordinal_type nParticles;

    amcd_real_type_1d_view molecular_weigths;
    amcd_real_type_1d_view aerosol_density;

   };

   template<typename SpT>
  AerosolModel_ConstData<SpT> create_AerosolModelConstData(const AerosolModelData & amd) {
    AerosolModel_ConstData<SpT> data;

    data.nSpec_gas=amd.nSpec_gas_;
    data.nSpec=amd.nSpec_;
    data.nParticles=amd.nParticles_;
    data.molecular_weigths = amd.molecular_weigths_.template view<SpT>();
    data.aerosol_density = amd.aerosol_density_.template view<SpT>();
    data.simpol_params = amd.simpol_params_.template view<SpT>();

  return data;
  }

} // namespace TChem
#endif
