#ifndef __TCHEM_IMPL_SINGLEPARTICLEAEROSOLWATER_HPP__
#define __TCHEM_IMPL_SINGLEPARTICLEAEROSOLWATER_HPP__

#include "TChem_Util.hpp"

namespace TChem {
namespace Impl {

  template<typename ValueType, typename DeviceType>
struct AerosolWater_SingleParticle
{
    using value_type = ValueType;
    using device_type = DeviceType;
    using scalar_type = typename ats<value_type>::scalar_type;

    using real_type = scalar_type;
    /// sacado is value type
    using value_type_1d_view_type = Tines::value_type_1d_view<value_type,device_type>;
    using real_type_1d_view_type = Tines::value_type_1d_view<real_type,device_type>;
    using ordinal_type_1d_view_type = Tines::value_type_1d_view<ordinal_type,device_type>;
    using aerosol_model_data_type= AerosolModel_ConstData<device_type>;

    using const_real_type_1d_view_type = Tines::value_type_1d_view<const real_type,device_type>;

    template<typename MemberType>
    KOKKOS_INLINE_FUNCTION static
    void team_invoke(const MemberType& member,
        ordinal_type i_part,
        real_type number_conc,
        real_type state[],  // TODO: figure out how the state vector actually works
        const AerosolModelData amd
        //const aerosol_ion_pair_type ionpair,
        //const std::map<std::string, int> ion_idx_map,
        //const std::map<std::string, real_type> atomic_weights  // TODO: pack into ionpair objects
        )
    {   
        // NOTE: Cannot use std::string with cuda 
        // Declare and initialize variables for later use
        
        int anion_idx;
        int cation_idx;
        int qty;
        std::string atom_name;
        std::string ion_name;
        real_type ion_molecular_weight;
        real_type cation;
        real_type anion;
        real_type e_ac;
        real_type e_aa;
        real_type ALPHA = -100.0;
        real_type ion_conc;
        real_type water_content;
        real_type NW;
        real_type MW;
        real_type ZW;
        real_type MW_H20 = 999.7351; // 55.51*18.01
        ordinal_type ion_idx;

        real_type rh = state[amd.aerosol_sp_name_idx_.size()];
        
        // note that you have to convert std::string to a C-string for printf since printf is a C function
        //printf("[AerosolWater_SingleParticle::team_invoke] calculation type %s\n", ionpair.calc_type.c_str());
        auto aerowater_model = amd.aerowater_model;

        int n_ionpairs = aerowater_model.ion_pair_vec.size();
        for (int j=0; j<n_ionpairs; j++){
            auto i_ionpair = aerowater_model.ion_pair_vec[j];
            //printf("[AerosolWater_SingleParticle::team_invoke] %s\n", i_ionpair.calc_type.c_str());
            real_type effective_aw = 0.0;
            real_type molality = 0.0;
            if (i_ionpair.calc_type == "JACOBSON") {
                //printf("[AerosolWater_SingleParticle::team_invoke] in jacobson evaluation\n");

                effective_aw = max(rh, i_ionpair.jacobson_low_RH);
                
                // calculate molality from power series 
                for (int i_order=0; i_order < i_ionpair.jacobson_Y_j.size(); i_order++) {
                    molality += i_ionpair.jacobson_Y_j[i_order]*std::pow(effective_aw, i_order);
                }
                molality = std::pow(molality, 2.0);
                //printf("[AerosolWater_SingleParticle::team_invoke] Jacobson molality %f\n", molality);
                
                anion_idx = amd.aerosol_sp_name_idx_.at(i_ionpair.jacobson_anion);
                cation_idx = amd.aerosol_sp_name_idx_.at(i_ionpair.jacobson_cation);
                //printf("[AerosolWater_SingleParticle::team_invoke] Jacobson anion index %d\n", anion_idx);
                //printf("[AerosolWater_SingleParticle::team_invoke] Jacobson cation index %d\n", cation_idx);
                
                // Get molecular weights
                qty = i_ionpair.ion_quantities_map.at(i_ionpair.jacobson_anion);
                //printf("jacobson anion quantity: %d\n", qty);
                ion_molecular_weight = qty*i_ionpair.ion_molecular_weights.at(i_ionpair.jacobson_anion);
                //printf("[AerosolWater_SingleParticle::team_invoke] jacobson anion molec. weight %f\n", ion_molecular_weight);
                anion = state[anion_idx]/(ion_molecular_weight*1000.0);

                qty = i_ionpair.ion_quantities_map.at(i_ionpair.jacobson_cation);
                //printf("jacobson cation quantity: %d\n", qty);
                ion_molecular_weight = qty*i_ionpair.ion_molecular_weights.at(i_ionpair.jacobson_cation);
                //printf("[AerosolWater_SingleParticle::team_invoke] jacobson cation molec. weight %f\n", ion_molecular_weight);
                cation = state[cation_idx]/(ion_molecular_weight*1000.0);
            
                // soft-max function to allow for smooth transition between cation/anion saturation
                e_ac = exp(ALPHA*cation);
                e_aa = exp(ALPHA*anion);
                ion_conc = (cation*e_ac + anion*e_aa)/(e_ac + e_aa);
                water_content = (ion_conc/molality)*1000.0;
                //printf("[AerosolWater_SingleParticle::team_invoke] (Jacobson) water content %f\n", water_content);
                state[amd.aerosol_sp_name_idx_.at("H2O_aq")] += water_content;
                
            } // calc_type == JACOBSON
            if (i_ionpair.calc_type == "EQSAM") {
                //printf("[AerosolWater_SingleParticle::team_invoke] in eqsam evaluation\n");

                // avoid div by zero in EQSAM calc
                if (rh > 0.99){
                    effective_aw = 0.99;
                }
                else if (rh < 0.001){
                    effective_aw = 0.001;
                }
                else {
                    effective_aw = rh;
                }

                NW = i_ionpair.eqsam_NW;
                MW = 1000.0*i_ionpair.eqsam_MW;
                ZW = i_ionpair.eqsam_ZW;
                molality = std::pow((NW*(MW_H20/MW)*(1.0/effective_aw-1.0)), ZW);
                //printf("[AerosolWater_SingleParticle::team_invoke] EQSAM molality %f\n", molality);
                
                for (int i=0; i < i_ionpair.ions.size(); i++) {
                    ion_name = i_ionpair.ions.at(i);
                    
                    qty = i_ionpair.ion_quantities_map.at(ion_name);
                    //printf("EQSAM ion quantity: %d\n", qty);
                    ion_molecular_weight = qty*i_ionpair.ion_molecular_weights.at(ion_name);
                    //printf("[AerosolWater_SingleParticle::team_invoke] EQSAM ion molec. weight %f\n", ion_molecular_weight);

                    //printf("[AerosolWater_SingleParticle::team_invoke] ion molecular weight %f\n", ion_molecular_weight);
                    ion_conc = state[amd.aerosol_sp_name_idx_.at(ion_name)]/ion_molecular_weight;
                    water_content = ion_conc/molality;
                    
                    state[amd.aerosol_sp_name_idx_.at("H2O_aq")] += water_content;

                }
                //printf("[AerosolWater_SingleParticle::team_invoke] (EQSAM) water content %f\n", water_content);
                    
        } // calc_type == EQSAM
        } // loop over ion_pair structs
    } // TeamInvoke

    }; // struct AerosolWater_SingleParticle
} // namespace Impl
} // namespace TChem
#endif