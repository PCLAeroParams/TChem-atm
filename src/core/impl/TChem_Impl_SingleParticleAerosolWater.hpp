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
        const ordinal_type i_part,
        //real_type number_conc,
        const real_type_1d_view_type& number_conc,
        const value_type_1d_view_type& state,             // NOTE: should this be const?
        const aerosol_model_data_type& amcd,
        const ordinal_type rh_idx,
        const ordinal_type aqueous_water_idx
        )
    {   
        // Declare and initialize variables for later use        
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
        // NOTE: should this be type static constexpr?
        static constexpr real_type MW_H20 = H2O_MOLALITY*H2O_MW; // constants defined in TChem_Util.hpp 
        ordinal_type ion_idx;
        ordinal_type anion_idx;
        ordinal_type cation_idx;
        ordinal_type ion_species_idx;
        ordinal_type qty;

        real_type rh = state(rh_idx);
        
        // note that you have to convert std::string to a C-string for printf since printf is a C function
        //printf("[AerosolWater_SingleParticle::team_invoke] calculation type %s\n", ionpair.calc_type.c_str());

        auto aerowater_params = amcd.aerowater_params;
        //auto aerowater_model = amd.aerowater_model;

        //int n_ionpairs = aerowater_model.ion_pair_vec.size();
        int n_ionpairs = amcd.nAerosolWater_ionpairs; // note have to use extent instead of size since this a kokkos view
        //printf("Number of elements in aerowater_params_: %d\n", n_ionpairs);
        for (int j=0; j<n_ionpairs; j++){
            //auto i_ionpair = aerowater_model.ion_pair_vec[j];
            auto i_ionpair = aerowater_params(j);
            //printf("[AerosolWater_SingleParticle::team_invoke] %s\n", i_ionpair.calc_type.c_str());
            real_type effective_aw = 0.0;
            real_type molality = 0.0;

            switch( i_ionpair.calc_type ){
                case JACOBSON:
                    //printf("[AerosolWater_SingleParticle::team_invoke] in jacobson evaluation\n");
                    //effective_aw = ats<real_type>::max(rh, i_ionpair.jacobson_low_rh);
                    
                    // if rh larger than jacobson low rh, return rh, otherwise return jacobson low rh
                    effective_aw = rh > i_ionpair.jacobson_low_rh ? rh : i_ionpair.jacobson_low_rh;

                    // calculate molality from power series 
                    for (int i_order=0; i_order < N_JACOBSON_COEFFS; i_order++) {
                        molality += i_ionpair.jacobson_y_j[i_order]*ats<real_type>::pow(effective_aw, i_order);
                    }
                    molality = ats<real_type>::pow(molality, 2.0);
                    //printf("[AerosolWater_SingleParticle::team_invoke] Jacobson molality %f\n", molality);

                    // calculate cation molecular weight
                    qty = i_ionpair.ion_quantities[i_ionpair.jacobson_cation];
                    //printf("jacobson cation quantity: %d\n", qty); 
                    ion_molecular_weight = qty*i_ionpair.ion_molec_weight[i_ionpair.jacobson_cation];
                    //printf("[AerosolWater_SingleParticle::team_invoke] jacobson cation molec. weight %f\n", ion_molecular_weight);

                    cation_idx = i_ionpair.ion_species_index[i_ionpair.jacobson_cation]; // index of the cation in the state vector
                    cation = state(amcd.nSpec_gas + amcd.nSpec*i_part + cation_idx)/(ion_molecular_weight*1000.0);

                    // calculate anion molecular weight
                    qty = i_ionpair.ion_quantities[i_ionpair.jacobson_anion];
                    //printf("jacobson anion quantity: %d\n", qty); 
                    ion_molecular_weight = qty*i_ionpair.ion_molec_weight[i_ionpair.jacobson_anion];
                    //printf("[AerosolWater_SingleParticle::team_invoke] jacobson anion molec. weight %f\n", ion_molecular_weight);

                    anion_idx = i_ionpair.ion_species_index[i_ionpair.jacobson_anion]; // index of the anion in the state vector
                    anion = state(amcd.nSpec_gas + amcd.nSpec*i_part + anion_idx)/(ion_molecular_weight*1000.0);

                    // soft-max function to allow for smooth transition between cation/anion saturation
                    e_ac = ats<real_type>::exp(ALPHA*cation);
                    e_aa = ats<real_type>::exp(ALPHA*anion);
                    ion_conc = (cation*e_ac + anion*e_aa)/(e_ac + e_aa);
                    water_content = (ion_conc/molality)*1000.0;
                    //printf("[AerosolWater_SingleParticle::team_invoke] (Jacobson) water content %f\n", water_content);
                    state(amcd.nSpec_gas + amcd.nSpec*i_part + aqueous_water_idx) += water_content;

                    break;

                case EQSAM:
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

                    NW = i_ionpair.eqsam_nw;
                    MW = 1000.0*i_ionpair.eqsam_mw;
                    ZW = i_ionpair.eqsam_zw;
                    molality = ats<real_type>::pow((NW*(MW_H20/MW)*(1.0/effective_aw-1.0)), ZW);
                    //printf("[AerosolWater_SingleParticle::team_invoke] (EQSAM) molality %f\n", molality);

                    for (ordinal_type ion_idx = 0; ion_idx < i_ionpair.num_ions; ion_idx++) {
                        
                        qty = i_ionpair.ion_quantities[ion_idx];
                        //printf("EQSAM ion quantity: %d\n", qty);
                        ion_molecular_weight = qty*i_ionpair.ion_molec_weight[ion_idx];
                        //printf("[AerosolWater_SingleParticle::team_invoke] (EQSAM) ion molec. weight %f\n", ion_molecular_weight);

                        //printf("[AerosolWater_SingleParticle::team_invoke] ion molecular weight %f\n", ion_molecular_weight);
                        ion_species_idx = i_ionpair.ion_species_index[ion_idx]; // index of the ion in the state vector

                        ion_conc = state(amcd.nSpec_gas + amcd.nSpec*i_part + ion_species_idx)/ion_molecular_weight;
                        water_content = ion_conc/molality;
                        
                        state(amcd.nSpec_gas + amcd.nSpec*i_part + aqueous_water_idx) += water_content;

                    }
                    //printf("[AerosolWater_SingleParticle::team_invoke] (EQSAM) water content %f\n", water_content);

                    break;
            } // switch calc_type

        } // loop over ion_pair structs
    } // TeamInvoke

    }; // struct AerosolWater_SingleParticle
} // namespace Impl
} // namespace TChem
#endif