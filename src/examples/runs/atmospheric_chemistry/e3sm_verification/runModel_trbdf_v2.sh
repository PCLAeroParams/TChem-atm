# TCHEM_INSTALL_PATH=/home/odiazib/TChem_atm/tchem_atm/code/HOST/RELEASE/install/tchem_atm
# exec=$TCHEM_INSTALL_PATH/examples/TChem_AtmosphericChemistryE3SM.x
# input=uci_v2_test3.yaml
# inputFile="input_conditions_multi_col.yaml"
# thread_size=160
run_this="OMP_NUM_THREADS=$thread_size OMP_PLACES=threads  OMP_PROC_BIND=close $exec \
    --chemfile=$input \
    --inputfile=$inputFile \
    --outputfile="${output_dir}/full_gas_trbdf.dat" \
    --outputfile_times="${output_dir}/wall_times_trbdf.json" \
    --time-iterations-per-interval=1 \
    --tol-time=1e-3 \
    --dtmin=1800 \
    --dtmax=1800 \
    --tend=1800 \
    --atol-newton=1e-13 \
    --rtol-newton=1e-8 \
    --max-newton-iterations=20 \
    --max-time-iterations=20000"

echo $run_this
eval $run_this
#          --atol-newton=1e-18 \
#          --rtol-newton=1e-8 \
