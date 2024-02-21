# TCHEM_INSTALL_PATH=/Users/odiazib/Documents/CODE/gitlab-ex/tchem_atm/code/ATM/RELEASE/install/tchem_atm
# exec=$TCHEM_INSTALL_PATH/examples/TChem_RateofProgress.x
# input=uci_v2_test3.yaml
# inputFile="input_conditions_multi_col.yaml"
run_this="$exec --chemfile=$input \
          --inputfile=$inputFile \
          --outputfile="${output_dir}/rate_of_progress.dat""

echo $run_this
eval $run_this
