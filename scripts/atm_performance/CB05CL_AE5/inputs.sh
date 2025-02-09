inputs=${TCHEM_INSTALL_PATH}/examples/runs/atmospheric_chemistry/CB05CL_AE5
chemfile=${inputs}/config_full_gas.yaml

export scenario_n_inputs="--chemfile=$chemfile \
                          --max-time-iterations=5 \
                          --write-time-profiles=false \
                          --inputfile_particles=$inputfile_particles"
