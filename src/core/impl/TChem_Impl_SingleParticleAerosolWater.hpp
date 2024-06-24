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
        real_type state[],  // TODO: figure out how the state vector actually works
        const aerosol_ion_pair_type ionpair,
        const std::map<std::string, int> ion_idx_map,
        const std::map<std::string, real_type> atomic_weights  // TODO: pack into ionpair objects
        )
    {   
        // NOTE: Cannot use std::string with cuda 
        // Declare and initialize variables for later use
        real_type effective_aw = 0.0;
        real_type molality = 0.0;
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

        real_type rh = state[2];

        // note that you have to convert std::string to a C-string for printf since printf is a C function
        //printf("[AerosolWater_SingleParticle::team_invoke] calculation type %s\n", ionpair.calc_type.c_str());

        if (ionpair.calc_type == "JACOBSON") {
            //printf("[AerosolWater_SingleParticle::team_invoke] in jacobson evaluation\n");

            effective_aw = max(rh, ionpair.jacobson_low_RH);

            // calculate molality from power series 
            for (int i_order; i_order < ionpair.jacobson_Y_j.size(); i_order++) {
                molality += ionpair.jacobson_Y_j[i_order]*std::pow(effective_aw, i_order);
            }
            molality = std::pow(molality, 2.0);
            //printf("[AerosolWater_SingleParticle::team_invoke] molality %f\n", molality);

            // calculate ion conc for anion and cation
            anion_idx = ion_idx_map.at(ionpair.jacobson_anion); 
            cation_idx = ion_idx_map.at(ionpair.jacobson_cation);
            
            // calculate the molecular weight (since some ions are diatomic like Cl2, need the number of atoms)
            qty = ionpair.ion_quantities.at(cation_idx);
            atom_name = ionpair.jacobson_cation.substr(0, ionpair.jacobson_cation.find('_'));
            ion_molecular_weight = qty*atomic_weights.at(atom_name);
            cation = state[cation_idx]/(ion_molecular_weight*1000.0);
            //printf("[AerosolWater_SingleParticle::team_invoke] cation index %d\n", cation_idx);
            //printf("[AerosolWater_SingleParticle::team_invoke] atom name %s\n", atom_name.c_str());
            //printf("[AerosolWater_SingleParticle::team_invoke] ion molecular weight %f\n", ion_molecular_weight);
            //printf("[AerosolWater_SingleParticle::team_invoke] cation conc. %e\n", cation);

            // calculate the molecular weight (since some ions are diatomic like Cl2, need the number of atoms)
            qty = ionpair.ion_quantities.at(anion_idx);
            atom_name = ionpair.jacobson_anion.substr(0, ionpair.jacobson_anion.find('_'));
            ion_molecular_weight = qty*atomic_weights.at(atom_name);
            anion = state[anion_idx]/(ion_molecular_weight*1000.0);
            //printf("[AerosolWater_SingleParticle::team_invoke] anion index %d\n", anion_idx);
            //printf("[AerosolWater_SingleParticle::team_invoke] atom name %s\n", atom_name.c_str());
            //printf("[AerosolWater_SingleParticle::team_invoke] ion molecular weight %f\n", ion_molecular_weight);
            //printf("[AerosolWater_SingleParticle::team_invoke] anion conc. %e\n", anion);

            // soft-max function to allow for smooth transition between cation/anion saturation
            e_ac = exp(ALPHA*cation);
            e_aa = exp(ALPHA*anion);
            ion_conc = (cation*e_ac + anion*e_aa)/(e_ac + e_aa);
            water_content = (ion_conc/molality)*1000.0;
            //printf("[AerosolWater_SingleParticle::team_invoke] water content %f\n", water_content);
            state[3] += water_content;

        } // calc_type == JACOBSON
        if (ionpair.calc_type == "EQSAM") {
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

                NW = ionpair.eqsam_NW;
                MW = 1000.0*ionpair.eqsam_MW;
                ZW = ionpair.eqsam_ZW;
                molality = std::pow((NW*(MW_H20/MW)*(1.0/effective_aw-1.0)), ZW);
                //printf("[AerosolWater_SingleParticle::team_invoke] molality %f\n", molality);

                for (int i; i < ionpair.ions.size(); i++) {
                    ion_name = ionpair.ions.at(i);
                    
                    atom_name = ion_name.substr(0, ion_name.find('_'));
                    //printf("atom name %s\n", atom_name.c_str());

                    ion_idx = ion_idx_map.at(ion_name);
                    //printf("ion idx %d\n", ion_idx);
                    qty = ionpair.ion_quantities.at(i);
                    ion_molecular_weight = qty*atomic_weights.at(atom_name);
                    //printf("[AerosolWater_SingleParticle::team_invoke] ion molecular weight %f\n", ion_molecular_weight);
                    ion_conc = state[ion_idx]/ion_molecular_weight;
                    water_content = ion_conc/molality;
                    //printf("[AerosolWater_SingleParticle::team_invoke] water content %f\n", water_content);
                    state[3] += water_content;

                }
                //printf("[AerosolWater_SingleParticle::team_invoke] water content %f\n", state[3]);
                    
        } // calc_type == EQSAM

    } // TeamInvoke

    }; // struct AerosolWater_SingleParticle
} // namespace Impl
} // namespace TChem
#endif