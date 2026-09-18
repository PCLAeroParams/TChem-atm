exec=$TCHEM_INSTALL_PATH/examples/TChem_AerosolChemistry_CVODE_K.x

run_this="$exec --chemfile=config_full_gas.yaml \
        --aerofile=mechanism_aero.yaml \
        --inputfile_particles=scenario_conditions_particle.yaml \
        --outputfile=cvode_bgmr_cb05_batched.dat \
        --outputfile_times=cvode_bgmr_cb05_times_batched.json \
        --team_thread_size=8 \
        --vector_thread_size=32 \
        --solver_type=2 \
        --bgmr_max_iter=10 \
        --use_cloned_samples=true \
        --batch_size=10000 \
        --number_of_particles=1000 \
        --rtol-time=1e-8 \
        --atol-time=1e-20 \
        --tbeg=0 \
        --tend=10 \
        --dtmin=1 \
        --max-time-iterations=10000 \
        --write-time-profiles=false \
        --verbose=true"

echo $run_this
eval $run_this