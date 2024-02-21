#exec=$TCHEM_INSTALL_PATH/examples/TChem_AtmosphericChemistryE3SM.explicit_euler.x
#input=uci_explicit_mech.yaml
#inputFile="input_condition_explicit_part_multi_col.yaml"
#inputFile="input_conditions.yaml"
run_this="$exec --chemfile=$input \
          --inputfile=$inputFile \
          --outputfile="${output_dir}/full_gas_stratosphere.dat" \
          --outputfile_times="${output_dir}/wall_times.json" \
          --time-iterations-per-interval=100 \
          --tol-time=1e-6 \
          --dtmin=1800 \
          --dtmax=1800 \
          --tend=1800 \
          --atol-newton=1e-18 \
          --rtol-newton=1e-8 \
          --max-newton-iterations=20 \
          --max-time-iterations=20000"

echo $run_this
eval $run_this
