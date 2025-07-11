#measure wall times for the following mechanisms
dirs=(CB05CL_AE5_w_simpolSOA)
sacado_flags=(no_sacado)
export experiment_base=TChem_AerosolChemistry

source ./solver_conf.sh
gpu_solver_strList=(
  "sundials_dense"
)
gpu_exe_string=(
  "_CVODE_K"
)

exe_strList=(
"runCUDA_dense.sh" 
#"runCUDAmap_dense.sh"
)
for exe in "${exe_strList[@]}"; do
for i in "${!gpu_solver_strList[@]}"; do
  export experiment_name="${gpu_solver_strList[i]}"
  export experiment_suffix="${gpu_exe_string[i]}"
  params_fxn="${experiment_name}_params"
  printf "setting variables for %s case (GPU)\n" "${experiment_name}"
  echo exe ${exe}
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
done
