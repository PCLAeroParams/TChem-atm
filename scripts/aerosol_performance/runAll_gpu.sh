#measure wall times for the following mechanisms
dirs=(CB05CL_AE5_w_simpolSOA)
sacado_flags=(no_sacado)
export experiment_base=TChem_AerosolChemistry

source ./solver_conf.sh
gpu_solver_strList=(
  #"rhss"
  "rhs"
  #"jac"
  #"sundials_gmres"
)
gpu_exe_string=(
  #"_RHSs"
  "_RHSs"
  #"_RHSs"
  #"_CVODE_K"
)

exe_strList=(
"runCUDA.sh"
#"runCUDAmap.sh"
#"runCUDA_best.sh"
#"runCUDA_max_team.sh"
)

nRepeat=10

for exe in "${exe_strList[@]}"; do
for i in "${!gpu_solver_strList[@]}"; do
  # Set the base experiment name and suffix from the lists
  base_experiment_name="${gpu_solver_strList[i]}"
  export experiment_suffix="${gpu_exe_string[i]}"

  # Define the parameters function based on the base name
  params_fxn="${base_experiment_name}_params"

  # add for loop for nRepeat from 1 to nRepeat (inclusive)
  for ((j=1; j<=nRepeat; j++)); do
    # Set the initial experiment name to the base name for this iteration
    export experiment_name="${base_experiment_name}"

    # if nRepeat != 1 then modify experiment name (append "_run[j]" where j is the for loop iteration index)
    if [ "$nRepeat" -ne 1 ]; then
      export experiment_name="${experiment_name}_run${j}"
    fi

    printf "setting variables for %s case (GPU)\n" "${experiment_name}"
    echo exe ${exe}
    ${params_fxn}

    for dir in ${dirs[@]}; do
      for sacado in ${sacado_flags[@]}; do
        # This check uses the base name to ensure it is consistent for all repeats
        if [ "${base_experiment_name}" == "expEuler" ] && [ "${sacado}" == "sacado" ]; then
          continue
        fi
        export sflag=$sacado
        run_this="cd $dir;./$exe;cd -"
        echo $run_this
        eval $run_this
      done
    done
    echo "end of outer loop"

  # end scope of nRepeat for loop
  done
done
done
