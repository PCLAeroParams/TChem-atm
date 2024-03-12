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
  // number of particles
  ordinal_type nParticles_;


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


   };

   template<typename SpT>
  AerosolModel_ConstData<SpT> create_AerosolModelConstData(const AerosolModelData & kmd) {
    AerosolModel_ConstData<SpT> data;

  return data;
  }

} // namespace TChem
#endif
