exec=$TCHEM_INSTALL_PATH/examples/TChem_AerosolChemistry_KokkosKernels.x

run_this="$exec --chemfile=config_gas.yaml \
	  --aerofile=test_SIMPOL_phase_transfer.yaml \
          --inputfile_particles=scenario_conditions_particle.yaml \
	  --outputfile=full_gas_kk.dat \
	  --time-iterations-per-interval=10 \
	  --tol-time=1e-3 \
          --atol-time=1e-12 \
	  --dtmin=10 \
          --dtmax=10\
          --tend=1000 \
          --atol-newton=1e-12 \
          --rtol-newton=1e-8 \
          --max-newton-iterations=20 \
          --max-time-iterations=200000"

echo $run_this
eval $run_this
