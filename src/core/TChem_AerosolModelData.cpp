#include "TChem_AerosolModelData.hpp"

#include <cstdarg>
#include <algorithm>

namespace TChem {


AerosolModelData::AerosolModelData(const std::string &mechfile,
                                   const KineticModelData& kmd,
 std::ostream& echofile) {
    setGasParameters(kmd);
    initFile(mechfile, echofile);
}

AerosolModelData::AerosolModelData(const std::string &mechfile,
                                   const KineticModelData& kmd) {
    std::ofstream echofile;
    echofile.open("kmod.echo");
    setGasParameters(kmd);
    initFile(mechfile ,echofile);
    echofile.close();
}

void AerosolModelData::setGasParameters(const KineticModelData& kmd){
  // FIXME: I maybe need to add a map of gas species as an input in initChem
  const ordinal_type n_active_species = kmd.nSpec_- kmd.nConstSpec_;
  nSpec_gas_=n_active_species;
  nConstSpec_gas_=kmd.nConstSpec_;//
  // maps from kmd has order of gas species.
  gas_sp_name_idx_ = kmd.species_indx_;
  is_gas_parameters_set_=true;
  printf("[AerosolModelData::setGasParameters] Setting gas parameters... nSpec_gas_ %d  nConstSpec_gas_%d \n", nSpec_gas_, nConstSpec_gas_);
 }
void AerosolModelData::initFile(const std::string &mechfile,
                                std::ostream& echofile){

  YAML::Node doc = YAML::LoadFile(mechfile);
  // FIXME: add error checking in yaml parser
  if (doc["NCAR-version"]) {
    if (verboseEnabled) {
      printf("[AerosolModelData::initFile] Using TChem parser for atmospheric chemistry\n");
    }
    initChem(doc, echofile);
  } else {
    TCHEM_CHECK_ERROR(true,"we do not have a parser for this file." );
    // exit
  }
}

int AerosolModelData::initChem(YAML::Node &root,
                               std::ostream& echofile) {

    TCHEM_CHECK_ERROR(!is_gas_parameters_set_,"[AerosolModelData::initChem] Error: gas parameters are not set use: setGasParameters(n_active_gas) ." );

    //gas_idx_sp_name : given index return species name
    std::map<int, std::string> gas_idx_sp_name;
    for (std::map<std::string, int>::iterator
         i = gas_sp_name_idx_.begin();
         i != gas_sp_name_idx_.end(); ++i)
    gas_idx_sp_name[i->second] = i->first;

    // 1. let's make a map of aerosol species and gas species.
    // given the name of species return the YAML::NODE
    std::map<std::string, YAML::Node> gas_sp_info;
    // 2. get molecular weights and density of aerosol_species
    std::vector<real_type> mw_aerosol_sp, density_aero_sp;
    int i_aero_sp=0;
    // loops over species, only make map from aerosol species.
    // we assume that map of gas species was previously created in kmd.
    for (auto const &item : root["NCAR-version"] ) {
      auto type =item["type"].as<std::string>();
      if (type=="CHEM_SPEC"){
        // std::cout <<"sp_name " << item["name"]<<"\n";
        const auto sp_name = item["name"].as<std::string>();
        if (item["phase"]){
           auto phase_name =item["phase"].as<std::string>();
          //  std::cout <<"sp_name " << sp_name<< item["phase"]<<"\n";
           if (phase_name == "AEROSOL" ){
            // std::cout <<"sp_name " << sp_name<< item["phase"]<<"\n";
             aerosol_sp_name_idx_[sp_name] = i_aero_sp;
             i_aero_sp++;
            //  std::cout <<"molecular weight" << item["molecular weight [kg mol-1]"]<<"\n";
             mw_aerosol_sp.push_back(item["molecular weight [kg mol-1]"].as<real_type>());
             density_aero_sp.push_back(item["density [kg m-3]"].as<real_type>());
           }
        } else
        {
          // save YAML::NODE, we need it in simpol.
          // order of species can be different in gas conf file and aero mech file.
          gas_sp_info.insert(std::pair<std::string, YAML::Node>(sp_name,item));
        }// phase
      } // if CHEM_SPEC

    // get number of particles
      if (type=="AERO_REP_SINGLE_PARTICLE"){
        nParticles_= item["maximum computational particles"].as<int>();
      }// AERO_REP_SINGLE_PARTICLE

    } // item loop
    printf("[AerosolModelData::initChem] Done with species...\n");
    std::vector<SIMPOL_PhaseTransferType> simpol_info;
    std::vector<aerosol_ion_pair_type> ion_pair_vec;
    aerosol_water_type aerowater_model;
    for (auto const &item : root["NCAR-version"]) {
      auto type =item["type"].as<std::string>();
      if (type=="MECHANISM"){
        // loop over reactions
       auto reactions = item["reactions"];
       for (auto const &ireac : reactions) {
        // std::cout <<"reactions " << ireac<<"\n";
         auto reaction_type = ireac["type"].as<std::string>();
        //  std::cout <<"reactions " << reaction_type<<"\n";
         if (reaction_type=="SIMPOL_PHASE_TRANSFER")
         {
          SIMPOL_PhaseTransferType simpol_info_at;

          simpol_info_at.B1 = ireac["B"][0].as<real_type>();
          simpol_info_at.B2 = ireac["B"][1].as<real_type>();
          simpol_info_at.B3 = ireac["B"][2].as<real_type>();
          simpol_info_at.B4 = ireac["B"][3].as<real_type>();
          // getting aerosol and gas species index.
          const auto aero_sp = ireac["aerosol-phase species"].as<std::string>();
          // std::cout <<"aero_sp " << aero_sp<<"\n";
          auto it_aero = aerosol_sp_name_idx_.find(aero_sp);
          if (it_aero != aerosol_sp_name_idx_.end()) {
           // aerosol species are place after gas species.
           //Note: that we do not use number of particles here.
           simpol_info_at.aerosol_sp_index = it_aero->second + nSpec_gas_;
           } else {
           printf("[AerosolModelData::initChem] species does not exit %s in reactants of SIMPOL_PHASE_TRANSFER reaction No \n", aero_sp);
           TCHEM_CHECK_ERROR(true,"Yaml : Error when interpreting kinetic model" );
          }

           const auto gas_sp = ireac["gas-phase species"].as<std::string>();
           auto it_gas = gas_sp_name_idx_.find(gas_sp);
           if (it_gas != gas_sp_name_idx_.end()) {
            //
            simpol_info_at.gas_sp_index = it_gas->second;
           } else {
            printf("[AerosolModelData::initChem] Error : species does not exit %s in reactants of SIMPOL_PHASE_TRANSFER\n", gas_sp);
            TCHEM_CHECK_ERROR(true,"Yaml : Error when interpreting aerosol model" );
          }
          simpol_info.push_back(simpol_info_at);
         } // if SIMPOL_PHASE_TRANSFER 
         else {
          printf("[AerosolModelData::initChem] Warning : TChem (amd) did not parse this reaction type %s \n", reaction_type.c_str());
       }
      } // ireac loop
    } // if MECHANISM

    // Extract attributes for aerosol water sub model
    if (type=="SUB_MODEL_ZSR_AEROSOL_WATER"){
       printf("[AerosolModelData::initChem] Parsing aerosol water sub model\n"); 
       
       aerowater_model.name = item["name"].as<std::string>();
       aerowater_model.aero_phase = item["aerosol phase"].as<std::string>();
       aerowater_model.gas_phase_water = item["gas-phase water"].as<std::string>();
       aerowater_model.aerosol_phase_water = item["aerosol-phase water"].as<std::string>();

       auto ion_pairs = item["ion pairs"];
       // loop over ion pair and extract aerosol water calc. attribs
       for (auto it = ion_pairs.begin(); it != ion_pairs.end(); ++it){
        // each entry in ion pairs is a key value pair
        auto i_ionpair_name = it->first.as<std::string>();
        auto i_ionpair_data = it->second;
        printf("[AerosolModelData::initChem] Ion pair: %s\n", i_ionpair_name.c_str());
        auto calc_type = i_ionpair_data["type"].as<std::string>();
        printf("[AerosolModelData::initChem] ..Calc type: %s\n", calc_type.c_str());

        if (calc_type=="JACOBSON"){
              IonPair jacobson_ionpair;
              jacobson_ionpair.calc_type = calc_type;
              
              // parse ions, assign ion names and ion quantities
              auto ions = i_ionpair_data["ions"];
              for (auto i_ion = ions.begin(); i_ion != ions.end(); ++i_ion){
                auto i_ion_name = i_ion->first.as<std::string>();
                auto i_ion_data = i_ion->second;
                jacobson_ionpair.ions.push_back(i_ion_name);
                printf("[AerosolModelData::initChem] ..Ion: %s\n", i_ion_name.c_str());

                // if qty in i_ion_data then use value, otherwise set qty to 1
                int i_ion_qty = 1;
                if (i_ion_data["qty"]){
                  i_ion_qty = i_ion_data["qty"].as<int>(); // type alias for int?
                }
                printf("[AerosolModelData::initChem] ....Ion qty: %d\n", i_ion_qty);
                //jacobson_ionpair.ion_quantities.push_back(i_ion_qty);
                // NOTE: in Python prototyping using dictionary (key=ion name, value=qty integer)
                jacobson_ionpair.ion_quantities_dict.push_back(std::make_pair(i_ion_name, i_ion_qty));

                // assign ion type (anion or cation)
                char charge = i_ion_name.back();
                if (charge == 'p'){
                  jacobson_ionpair.jacobson_cation = i_ion_name;
                } 
                if (charge == 'm'){
                  jacobson_ionpair.jacobson_anion = i_ion_name;
                }
                // possibly raise an exception if an unexpected charge character is encountered?

              } // i_ion loop

              printf("[AerosolModelData::initChem] ..Jacobson anion: %s\n", jacobson_ionpair.jacobson_anion.c_str());
              printf("[AerosolModelData::initChem] ..Jacobson cation: %s\n", jacobson_ionpair.jacobson_cation.c_str());

              // assign jacobson molality polynomial coefficients (Y_j)
              auto jacobson_coeffs = i_ionpair_data["Y_j"];
              for (auto const &i_coeff : jacobson_coeffs) {
                real_type y_j = i_coeff.as<real_type>();
                jacobson_ionpair.jacobson_Y_j.push_back(y_j);
                printf("[AerosolModelData::initChem] ..Jacobson Y_j: %f\n", y_j);
              } // Y_j jacobson molality coefficient loop

              jacobson_ionpair.jacobson_low_RH = i_ionpair_data["low RH"].as<real_type>();
              printf("[AerosolModelData::initChem] ..Jacobson low RH = %f\n", jacobson_ionpair.jacobson_low_RH);

              // push back ionpair attributes to vector
              ion_pair_vec.push_back(jacobson_ionpair);

          } // calc_type JACOBSON

          if (calc_type=="EQSAM"){
              IonPair eqsam_ionpair;
              eqsam_ionpair.calc_type = calc_type;

              // parse ions, assign ion names and ion quantities
              auto ions = i_ionpair_data["ions"];
              for (auto i_ion = ions.begin(); i_ion != ions.end(); ++i_ion){
                auto i_ion_name = i_ion->first.as<std::string>();
                auto i_ion_data = i_ion->second;
                eqsam_ionpair.ions.push_back(i_ion_name);
                printf("[AerosolModelData::initChem] ..Ion: %s\n", i_ion_name.c_str());

                // if qty in i_ion_data then use value, otherwise set qty to 1
                int i_ion_qty = 1;
                if (i_ion_data["qty"]){
                  i_ion_qty = i_ion_data["qty"].as<int>(); // type alias for int?
                }
                printf("[AerosolModelData::initChem] ....Ion qty: %d\n", i_ion_qty);
                //eqsam_ionpair.ion_quantities.push_back(i_ion_qty);
                // NOTE: in Python prototyping using dictionary (key=ion name, value=qty integer)
                eqsam_ionpair.ion_quantities_dict.push_back(std::make_pair(i_ion_name, i_ion_qty));
              }// i_ion loop

              eqsam_ionpair.eqsam_NW = i_ionpair_data["NW"].as<real_type>();
              eqsam_ionpair.eqsam_MW = i_ionpair_data["MW"].as<real_type>();
              eqsam_ionpair.eqsam_ZW = i_ionpair_data["ZW"].as<real_type>();

              printf("[AerosolModelData::initChem] ..EQSAM NW = %f\n", eqsam_ionpair.eqsam_NW );
              printf("[AerosolModelData::initChem] ..EQSAM MW = %f\n", eqsam_ionpair.eqsam_MW );
              printf("[AerosolModelData::initChem] ..EQSAM ZW = %f\n", eqsam_ionpair.eqsam_ZW );

              // push back ionpair attributes to vector
              ion_pair_vec.push_back(eqsam_ionpair);

          } // calc_type EQSAM
          // TODO: Raise key error if calc_type other than jacobson or EQSAM encountered

  
       } // ion pair (it) loop

       aerowater_model.ion_pair_vec = ion_pair_vec;

    } // if SUB_MODEL_ZSR_AEROSOL_WATER
    } // item loop
    printf("[AerosolModelData::initChem] Done with reactions, mechanisms, sub-models...\n");

    // number of aerosol species
    nSpec_= density_aero_sp.size();

    nSimpol_tran_=simpol_info.size();
    auto nAerosolWater_ionpairs_ = aerowater_model.ion_pair_vec.size();
    printf("[AerosolModelData::initChem] Number of Simpol phase transfer %d \n",nSimpol_tran_ );
    printf("[AerosolModelData::initChem] Number of aerosol water ion pairs %d \n",nAerosolWater_ionpairs_ );
    // update simpol
    simpol_params_ = simpol_phase_transfer_type_1d_dual_view(do_not_init_tag("AMD::simpol_params_"), nSimpol_tran_);
    const auto simpol_params_host = simpol_params_.view_host();

    for (ordinal_type isimpol = 0; isimpol < nSimpol_tran_; isimpol++)
    {
      // SIMPOL_PhaseTransferType simpol_params_host_at_i;
      auto simpol_params_host_at_i = simpol_info[isimpol];
      // get name of gas species
      auto sp_name = gas_idx_sp_name[simpol_params_host_at_i.gas_sp_index];
      // get YAML::NODE using sp_name gas species
      auto gas_sp_at_i = gas_sp_info[sp_name];
      // std::cout <<"gas_sp_index " << simpol_params_host_at_i.gas_sp_index<<"\n";
      // std::cout <<"sp_name " << sp_name<<"\n";
      // std::cout <<"gas_sp_at_i " << gas_sp_at_i<<"\n";
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

    // NOTE-SF: Should I do something similar for the ZSR data structure that I create? I dont really understand 
    // the purpose of creating a 'host' kokkos object and when that is necessary

    molecular_weights_ = real_type_1d_dual_view(do_not_init_tag("AMD::molecular_weights"), nSpec_);
    auto molecular_weights_host = molecular_weights_.view_host();
    aerosol_density_= real_type_1d_dual_view(do_not_init_tag("AMD::aerosol_density_"), nSpec_);
    auto aerosol_density_host = aerosol_density_.view_host();
    for (int i = 0; i < nSpec_; i++)
    {
      molecular_weights_host(i) = mw_aerosol_sp[i];
      aerosol_density_host(i) = density_aero_sp[i];
    }

    molecular_weights_.modify_host();
    aerosol_density_.modify_host();
    simpol_params_.modify_host();

    molecular_weights_.sync_device();
    aerosol_density_.sync_device();
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
          printf("[AerosolModelData::scenarioConditionParticles] species does not exit %s in reactants of SIMPOL_PHASE_TRANSFER reaction No \n", species_name);
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
