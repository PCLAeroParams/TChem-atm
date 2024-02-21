# TCHEM_INSTALL_PATH=/home/odiazib/TChem_atm/tchem_atm/code/HOST/RELEASE/install/tchem_atm
# input=uci_v2_test3.yaml
# inputFile="input_conditions_multi_col.yaml"
# thread_size=160
# exec=$TCHEM_INSTALL_PATH/examples/TChem_AtmosphericChemistryE3SM_CVODE.x
run_this="OMP_NUM_THREADS=$thread_size OMP_PLACES=threads  OMP_PROC_BIND=close $exec \
  --chemfile=$input \
  --inputfile=$inputFile \
  --use-cvode=true \
  --outputfile_times="${output_dir}/wall_times_cvode.json" \
  --outputfile="${output_dir}/full_gas_cvode.dat" \
  --time-iterations-per-interval=1 \
  --jacobian-interval=4 \
  --tol-time=1e-3 \
  --atol-time=1e-3 \
  --dtmin=1800 \
  --dtmax=1800 \
  --tend=1800 \
  --atol-newton=1e-3 \
  --rtol-newton=1e-3 \
  --max-newton-iterations=20 \
  --max-time-iterations=20000"

echo $run_this
eval $run_this
#          --atol-newton=1e-18 \
#          --rtol-newton=1e-8 \
