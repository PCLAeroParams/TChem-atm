exec=$TCHEM_INSTALL_PATH/examples/TChem_AtmosphericChemistry.x

run_this="$exec --chemfile=config_full_gas.yaml \
          --outputfile=full_gas.dat \
          --time-iterations-per-interval=10 \
          --tol-time=1e-3 \
          --dtmin=1e-20 \
          --dtmax=10 \
          --tend=600\
          --atol-newton=1e-12 \
          --rtol-newton=1e-8 \
          --max-newton-iterations=20 \
          --max-time-iterations=20000"
echo $run_this
eval $run_this
