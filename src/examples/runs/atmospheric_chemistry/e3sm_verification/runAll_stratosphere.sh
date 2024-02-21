#! /bin/bash

#MJS_INSTALL_PATH=${HOME}/tchem_atm/code/CPU/RELEASE/install
#INSTALL_PATH=$MJS_INSTALL_PATH
# NOTE: use sacado or no_sacado?
#TCHEM_INSTALL_PATH=${INSTALL_PATH}/sacado/tchem_atm/examples
#TCHEM_INSTALL_PATH=/ascldap/users/odiazib/Documents/Oscar/CODE/gitlab-ex/TChem_QTI/ATM/RELEASE/install/tchem_atm
TCHEM_INSTALL_PATH=/Users/odiazib/Documents/CODE/gitlab-ex/tchem_atm/code/ATM/RELEASE/install/tchem_atm/examples

ROOT=/Users/odiazib/Documents/CODE/gitlab-ex/tchem_atm/verification/e3sm_verification
sacado_flags=(no_sacado sacado)
experiment_base=TChem_AtmosphericChemistryE3SM
export input=uci_explicit_mech.yaml
infile="input_condition_explicit_part_multi_col.yaml"
input_root=$ROOT"/inputs"
output_root=$ROOT"/outputs/tchem"
export thread_size=16


multicol_col=(
  "LA"
  "BRW"
)

case_list=(
  "runModel_stratosphere.sh"
  "runNetProductionRates_stratosphere.sh"
)

case_str=(
  "stratoSolver"
  "stratoNetProductionRates"
)

exe_list=(
  "TChem_AtmosphericChemistryE3SM.explicit_euler.x"
  "TChem_NetProductionRates.x"
)

for col in ${multicol_col[@]}; do
  for i_case in ${!case_list[@]}; do
    export inputFile="${input_root}/${col}/${infile}"
    export exec=${TCHEM_INSTALL_PATH}/${exe_list[i_case]}
    export output_dir="${output_root}/${col}/${case_str[i_case]}"
    mkdir -p $output_dir
    echo "running ${case_list[i_case]} for col = ${col}"
    echo "input file = ${inputFile}"
    echo "exe = ${exec}"
    ./"${case_list[i_case]}"
  done
done

# ./runModel.sh
# ./runModel_trbdf.sh
# ./runModel_sundials.sh
# ./runNetProductionRates.sh
# ./runRateofProgress.sh
# ./runModel_stratosphere.sh
# ./runNetProductionRates_stratosphere.sh
