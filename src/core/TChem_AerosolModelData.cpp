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
    nParticles_=5;
    molecular_weigths_ = real_type_1d_dual_view(do_not_init_tag("AMD::molecular_weigths"), nSpec_);
    auto molecular_weigths_host = molecular_weigths_.view_host();
    molecular_weigths_host(0)=0.04607;
    molecular_weigths_host(1)=0.01801;

    aerosol_density_= real_type_1d_dual_view(do_not_init_tag("AMD::aerosol_density_"), nSpec_);
    auto aerosol_density_host = aerosol_density_.view_host();
    aerosol_density_host(0)=1e3;
    aerosol_density_host(1)=1e3;

    //
    ordinal_type nSimpol_tran=1;

    simpol_params_ = simplo_phase_transfer_type_1d_dual_view(do_not_init_tag("AMD::simpol_params_"), nSimpol_tran);
    const auto simpol_params_host = simpol_params_.view_host();

    SIMPOL_PhaseTransferType simpol_params_host_at_i;
    // auto simpol_params_host_at_i = simpol_params_host(0);
    simpol_params_host_at_i.aerosol_sp_index=1;
    simpol_params_host_at_i.B1= -1.97E+03;
    simpol_params_host_at_i.B2=2.91E+00;
    simpol_params_host_at_i.B3=1.96E-03;
    simpol_params_host_at_i.B4=-4.96E-01;

    simpol_params_host_at_i.gas_sp_index=0;
    simpol_params_host_at_i.diffusion_coeff=0.95E-05; // m2 s-;
    simpol_params_host_at_i.N_star=2.55;
    simpol_params_host_at_i.compute_alpha=true;
    simpol_params_host_at_i.molecular_weight=0.04607;
    simpol_params_host(0)=simpol_params_host_at_i;

    molecular_weigths_.modify_host();
    aerosol_density_.modify_host();
    simpol_params_.modify_host();

    molecular_weigths_.sync_device();
    molecular_weigths_.sync_device();
    simpol_params_.sync_device();

    return 0;
}





} // namespace TChem
