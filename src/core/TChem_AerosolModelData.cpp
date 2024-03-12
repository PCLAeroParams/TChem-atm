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
    return 0;
}





} // namespace TChem
