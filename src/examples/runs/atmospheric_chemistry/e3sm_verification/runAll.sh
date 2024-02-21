#! /bin/bash

TCHEM_INSTALL_PATH=/Users/odiazib/Documents/CODE/gitlab-ex/tchem_atm/code/ATM/RELEASE/install/tchem_atm/examples

ROOT=/Users/odiazib/Documents/CODE/gitlab-ex/tchem_atm/verification/e3sm_verification
sacado_flags=(no_sacado sacado)
experiment_base=TChem_AtmosphericChemistryE3SM
export input=uciv3020624.yaml
infile="input_conditions_multi_col.yaml"
input_root=$ROOT"/inputs"
output_root=$ROOT"/outputs/tchem"
export thread_size=16


multicol_col=(
  "LA"
  "BRW"
)

case_list=(
  "runModel.sh"
  "runModel_trbdf.sh"
  "runModel_trbdf_v2.sh"
  "runModel_sundials.sh"
  "runNetProductionRates.sh"
  "runRateofProgress.sh"
)

case_str=(
  "impEuler"
  "trbdf"
  "trbdf2"
  "cvode"
  "netProductionRates"
  "rateofProgress"
)

exe_list=(
  "TChem_AtmosphericChemistryE3SM.implicit_euler.x"
  "TChem_AtmosphericChemistryE3SM.x"
  "TChem_AtmosphericChemistryE3SM.x"
  "TChem_AtmosphericChemistryE3SM_CVODE.x"
  "TChem_NetProductionRates.x"
  "TChem_RateofProgress.x"
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
