#measure wall times for the following mechanisms
sacado_flags=(no_sacado)

# Aerosol timings
#dirs=(CB05CL_AE5_w_simpolSOA)
#export experiment_base=TChem_AerosolChemistry

# Gas timings
dirs=(CB05CL_AE5)
export experiment_base=TChem_GasChemistry

cpu_solver_strList=(
  # "trbdf"
  # "kokkoskernels"
  # "cvode"
  #"rhss"
)
cpu_exe_string=(
  # ""
  # "_KokkosKernels"
  # "_CVODE"
  #"_RHSs"
)
gpu_solver_strList=(
 #"trbdf"
  #"kokkoskernels"
  "rhss"
)
gpu_exe_string=(
  #""                   # note need to comment out when running Kokkos Kernels
  #"_KokkosKernels"
  "_RHSs"
)

rhss_params(){
  export numerical_params=""
}
trbdf_params(){
  t_iterPerInt=1
  min_dt='1e-3'
  max_dt='1'
  atol_t='1e-12'
  tol_time='1e-3'
  tend='1'
  export numerical_params="--tol-time=$tol_time \
                    --time-iterations-per-interval=$t_iterPerInt \
                    --dtmin=$min_dt \
                    --dtmax=$max_dt \
                    --atol-time=${atol_t} \
                    --tend=$tend \
                    --atol-newton=1e-18 \
                    --rtol-newton=1e-8 \
                    --max-newton-iterations=20 "

}

kokkoskernels_params(){
  t_iterPerInt=1
  min_dt='1e-3'
  max_dt='1'
  atol_t='1e-12'
  tol_time='1e-3'
  tend='1'
  export numerical_params="--tol-time=$tol_time \
                    --time-iterations-per-interval=$t_iterPerInt \
                    --dtmin=$min_dt \
                    --dtmax=$max_dt \
                    --atol-time=${atol_t} \
                    --tend=$tend "
}

cvode_params(){
  t_iterPerInt=10
  min_dt='1e-20'
  max_dt='1'
  atol_t='1e-12'
  tol_time='1e-3'
  tend='1'
  export numerical_params="--tol-time=$tol_time \
                    --time-iterations-per-interval=$t_iterPerInt \
                    --dtmin=$min_dt \
                    --dtmax=$max_dt \
                    --use-cvode=true \
                    --atol-time=${atol_t} \
                    --tend=$tend "
}

exe=runCUDA.sh

for i in "${!gpu_solver_strList[@]}"; do
  export experiment_name="${gpu_solver_strList[i]}"
  export experiment_suffix="${gpu_exe_string[i]}"
  params_fxn="${experiment_name}_params"
  printf "setting variables for %s case (GPU)\n" "${experiment_name}"
  ${params_fxn}
  for dir in ${dirs[@]}; do
    for sacado in ${sacado_flags[@]}; do
      if [ "${experiment_name}" == "expEuler" ] && [ "${sacado}" == "sacado" ]; then
        continue
      fi
      export sflag=$sacado
      run_this="cd $dir;./$exe;cd -"
      echo $run_this
      eval $run_this
    done
  done
  echo "end of outer loop"
done
