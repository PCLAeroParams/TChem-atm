#include "TChem_AerosolModelData.hpp"

#include <cstdarg>
#include <algorithm>

namespace TChem {


AerosolModelData::AerosolModelData(const std::string &mechfile,
 std::ostream& echofile) {
    initFile(mechfile, echofile);
}

AerosolModelData::AerosolModelData(const std::string &mechfile) {
    std::ofstream echofile;
    echofile.open("kmod.echo");
    initFile(mechfile,echofile);
    echofile.close();
}

void AerosolModelData::initFile(const std::string &mechfile,
                                    std::ostream& echofile){

  YAML::Node doc = YAML::LoadFile(mechfile);
  // FIXME: add error checking in yaml parser
  if (doc["NCAR-version"]) {
    if (verboseEnabled) {
      printf("Using TChem parser for atmospheric chemistry\n");
    }
    initChem(doc, echofile);
  } else {
    TCHEM_CHECK_ERROR(true,"we do not have a parser for this file." );
    // exit
  }
}

int AerosolModelData::initChem(YAML::Node &root, std::ostream& echofile) {
    // implimentation of parser
    nSpec_=2;
    nSpec_gas_=1;
    nParticles_=2;
    molecular_weigths_ = real_type_1d_dual_view(do_not_init_tag("AMD::molecular_weigths"), nSpec_);
    auto molecular_weigths_host = molecular_weigths_.view_host();
    molecular_weigths_host(0)=0.04607;
    molecular_weigths_host(1)=0.01801;

    aerosol_density_= real_type_1d_dual_view(do_not_init_tag("AMD::aerosol_density_"), nSpec_);
    auto aerosol_density_host = aerosol_density_.view_host();
    aerosol_density_host(0)=1e3;
    aerosol_density_host(1)=1e3;

    molecular_weigths_.modify_host();
    aerosol_density_.modify_host();

    molecular_weigths_.sync_device();
    molecular_weigths_.sync_device();

    return 0;
}





} // namespace TChem
