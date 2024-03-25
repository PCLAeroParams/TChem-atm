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
    // implementation of parser
    // FIXME: I need to get the number of gas species from gas mechanism
    nSpec_gas_=1;
    nConstSpec_gas_=1;
    // FIXME: I maybe need to add a map of gas species as an input in initChem
    // 1. let's make a map of aerosol species and gas species.
    // std::map<std::string, int> aerosol_sp_name_idx_;
    int i_aero_sp=0;
    // 2. FIXME: we maybe need to created a list with gas species from gas mechanism
    std::map<std::string, int> gas_sp_name_idx;
    std::vector<YAML::Node> gas_sp_info;
    int i_gas_sp=0;
    // 3. get molecular weitghs and density of aerosol_species
    std::vector<real_type> mw_aerosol_sp, density_aero_sp;
    for (auto const &item : root["NCAR-version"] ) {
      auto type =item["type"].as<std::string>();
      if (type=="CHEM_SPEC"){
        const auto sp_name = item["name"].as<std::string>();
        if (item["phase"]){
           auto phase_name =item["phase"].as<std::string>();
           if (phase_name == "AEROSOL" ){
            // std::cout <<"sp_name " << sp_name<< item["phase"]<<"\n";
             aerosol_sp_name_idx_.insert(std::pair<std::string, int>(sp_name, i_aero_sp));
             i_aero_sp++;
            //  std::cout <<"molecular weight" << item["molecular weight [kg mol-1]"]<<"\n";
             mw_aerosol_sp.push_back(item["molecular weight [kg mol-1]"].as<real_type>());
             density_aero_sp.push_back(item["density [kg m-3]"].as<real_type>());
           }
        } else
        {
          gas_sp_name_idx.insert(std::pair<std::string, int>(sp_name, i_gas_sp));
          gas_sp_info.push_back(item);
          i_gas_sp++;
        }// phase
      } // if CHEM_SPEC

      if (type=="AERO_REP_SINGLE_PARTICLE"){
        nParticles_= item["maximum computational particles"].as<int>();
      }// AERO_REP_SINGLE_PARTICLE

    } // item loop

    std::vector<SIMPOL_PhaseTransferType> simpol_info;
    for (auto const &item : root["NCAR-version"]) {
      auto type =item["type"].as<std::string>();
      if (type=="MECHANISM"){
       auto reactions = item["reactions"];
       auto reaction_type = reactions["type"].as<std::string>();
       if (reaction_type=="SIMPOL_PHASE_TRANSFER")
       {
        SIMPOL_PhaseTransferType simpol_info_at;
        simpol_info_at.B1 = reactions["B"][0].as<real_type>();
        simpol_info_at.B2 = reactions["B"][1].as<real_type>();
        simpol_info_at.B3 = reactions["B"][2].as<real_type>();
        simpol_info_at.B4 = reactions["B"][3].as<real_type>();
        // getting aerosol and gas species index.
        const auto aero_sp = reactions["aerosol-phase species"].as<std::string>();
         auto it = aerosol_sp_name_idx_.find(aero_sp);
         if (it != aerosol_sp_name_idx_.end()) {
          // aerosol species are place after gas species.
          //Note: that we do not use number of particles here.
          simpol_info_at.aerosol_sp_index = it->second + nSpec_gas_;
          } else {
          printf("species does not exit %s in reactants of SIMPOL_PHASE_TRANSFER reaction No \n", aero_sp);
          TCHEM_CHECK_ERROR(true,"Yaml : Error when interpreting kinetic model" );
         }

         const auto gas_sp = reactions["gas-phase species"].as<std::string>();
         if (it != gas_sp_name_idx.end()) {
          //
          simpol_info_at.gas_sp_index = it->second;
          } else {
          printf("species does not exit %s in reactants of SIMPOL_PHASE_TRANSFER reaction No \n", gas_sp);
          TCHEM_CHECK_ERROR(true,"Yaml : Error when interpreting kinetic model" );
         }
         simpol_info.push_back(simpol_info_at);
       }
      }
    }// item loop


    // number of aerosol species
    nSpec_= density_aero_sp.size();

    ordinal_type nSimpol_tran=simpol_info.size();
    // update simpol
    simpol_params_ = simplo_phase_transfer_type_1d_dual_view(do_not_init_tag("AMD::simpol_params_"), nSimpol_tran);
    const auto simpol_params_host = simpol_params_.view_host();

    for (ordinal_type isimpol = 0; isimpol < nSimpol_tran; isimpol++)
    {
      // SIMPOL_PhaseTransferType simpol_params_host_at_i;
      auto simpol_params_host_at_i = simpol_info[isimpol];
      // FIXME: Assuming that gas info is provided in the conf yaml
      // If gas species exit in gas mechanism, we maybe have an issue.
      // It is possible that is information is given in gas mechanism conf file.
      // get gas properties for simpol transfer.
      auto gas_sp_at_i = gas_sp_info[isimpol];

      simpol_params_host_at_i.diffusion_coeff=gas_sp_at_i["diffusion coeff [m2 s-1]"].as<real_type>(); // m2 s-;
      if (gas_sp_at_i["N star"]) {
        simpol_params_host_at_i.N_star=gas_sp_at_i["N star"].as<real_type>();
        simpol_params_host_at_i.compute_alpha=true;
      } else {
        simpol_params_host_at_i.compute_alpha=false;
        simpol_params_host_at_i.N_star=999;
      }
      simpol_params_host_at_i.molecular_weight=gas_sp_at_i["molecular weight [kg mol-1]"].as<real_type>();
      // save simpol_params_host_at_i in kokkos array
      simpol_params_host(isimpol)=simpol_params_host_at_i;
    } // nsimpol


    molecular_weigths_ = real_type_1d_dual_view(do_not_init_tag("AMD::molecular_weigths"), nSpec_);
    auto molecular_weigths_host = molecular_weigths_.view_host();
    aerosol_density_= real_type_1d_dual_view(do_not_init_tag("AMD::aerosol_density_"), nSpec_);
    auto aerosol_density_host = aerosol_density_.view_host();
    for (int i = 0; i < nSpec_; i++)
    {
      molecular_weigths_host(i) = mw_aerosol_sp[i];
      aerosol_density_host(i) = density_aero_sp[i];
    }

    molecular_weigths_.modify_host();
    aerosol_density_.modify_host();
    simpol_params_.modify_host();

    molecular_weigths_.sync_device();
    molecular_weigths_.sync_device();
    simpol_params_.sync_device();

    return 0;
}

void AerosolModelData::scenarioConditionParticles(const std::string &mechfile,
                                                  const ordinal_type nBatch,
                                                  real_type_2d_view_host& num_concentration,
                                                  real_type_2d_view_host& state)
{
  YAML::Node root = YAML::LoadFile(mechfile);

  if (root["particles"]){

    auto num_conc_info = root["particles"]["num_concentration"]["initial_value"];
    num_concentration = real_type_2d_view_host("num_concentration",nBatch, nParticles_ );
    // loop over batch
    for (ordinal_type ibatch = 0; ibatch < nBatch; ibatch++)
    {
      // Note: we assume all particles have equal num_concentration
      const auto value = num_conc_info[ibatch].as<real_type>();
      auto num_con_ibatch = Kokkos::subview(num_concentration, ibatch, Kokkos::ALL());
      Kokkos::deep_copy(num_con_ibatch,value);
    }

    // read initial conditions for part concentraion
    auto initial_state = root["particles"]["initial_state"];
    // density, pressure, temperature, active gas species, and invariant gas species
    const ordinal_type initial_index = 3 +nSpec_gas_+nConstSpec_gas_;
    for (auto const& sp_cond : initial_state) {
      auto species_name = sp_cond.first.as<std::string>();
      auto it = aerosol_sp_name_idx_.find(species_name);
      if (it == aerosol_sp_name_idx_.end()) {
          printf("species does not exit %s in reactants of SIMPOL_PHASE_TRANSFER reaction No \n", species_name);
          TCHEM_CHECK_ERROR(true,"Yaml : Error when interpreting kinetic model" );
      }
      const ordinal_type species_idx = it->second;
      std::cout << species_idx << "\n";
      auto species_info = initial_state[species_name]["initial_value"];

      for (ordinal_type ibatch = 0; ibatch < nBatch; ibatch++)
      {
       const auto value = species_info[ibatch].as<real_type>();
       auto state_ibatch = Kokkos::subview(state, ibatch, Kokkos::ALL());
       auto num_con_ibatch = Kokkos::subview(num_concentration, ibatch, Kokkos::ALL());
       for (int ipart = 0; ipart < nParticles_; ipart++)
      {
      state_ibatch(initial_index+nSpec_*ipart+species_idx) = value/num_con_ibatch(ipart);
    }
    } // batches

  }





  }

}




} // namespace TChem
