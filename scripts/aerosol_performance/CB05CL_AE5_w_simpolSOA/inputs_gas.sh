inputs=${TCHEM_HOST_INSTALL_PATH}/examples/runs/atmospheric_chemistry/CB05CL_AE5_w_simpolSOA
inputfile_particles=${inputs}/scenario_conditions_particle.yaml
chemfile=${inputs}/config_full_gas.yaml
# for gas maximum computational particles: 0
aerofile=mechanism_aero_gas.yaml

export scenario_n_inputs="--chemfile=$chemfile \
                          --aerofile=$aerofile \
			              --write-time-profiles=false \
                          --max-time-iterations=10000 \
                          --tend=10 \
                          --inputfile_particles=$inputfile_particles"
