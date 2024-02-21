#exec=$TCHEM_INSTALL_PATH/examples/TChem_NetProductionRates.x
run_this="$exec --chemfile=$input \
          --inputfile=$inputFile \
          --outputfile=${output_dir}/net_production_rates_stratosphere.dat" 

echo $run_this
eval $run_this
