inputs=${TCHEM_INSTALL_PATH}/examples/runs/atmospheric_chemistry/CB05CL_AE5_w_simpolSOA
inputfile_particles=${inputs}/scenario_conditions_particle.yaml
chemfile=${inputs}/config_full_gas.yaml
aerofile=${inputs}/mechanism_aero.yaml

export scenario_n_inputs="--chemfile=$chemfile \
                          --aerofile=$aerofile \
                          --inputfile_particles=$inputfile_particles"