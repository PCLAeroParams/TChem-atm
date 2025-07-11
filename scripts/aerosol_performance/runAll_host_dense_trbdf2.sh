#measure wall times for the following mechanisms
dirs=(CB05CL_AE5_w_simpolSOA)
sacado_flags=(no_sacado)
export experiment_base=TChem_AerosolChemistry
source ./solver_conf.sh

cpu_solver_strList=(
  "trbdf"
)
cpu_exe_string=(
  ""
)

exe=runHOST_dense.sh
for i in "${!cpu_solver_strList[@]}"; do
  export experiment_name="${cpu_solver_strList[i]}"
  export experiment_suffix="${cpu_exe_string[i]}"
  params_fxn="${experiment_name}_params"
  printf "setting variables for %s case (CPU)\n" "${experiment_name}"
  ${params_fxn}
  for dir in ${dirs[@]}; do
    for sacado in ${sacado_flags[@]}; do
      if [ "${experiment_name}" == "expEuler" ] && [ "${sacado}" == "sacado" ]; then
        continue
      elif [ "${experiment_name}" == "cvode" ] && [ "${sacado}" == "sacado" ]; then
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
