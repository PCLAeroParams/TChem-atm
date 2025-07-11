export DEVICE=CPU
export experiment_exe_root=${experiment_base}${experiment_suffix}
export sacado_flag=$sflag
echo "sacado flag = ${sacado_flag}"
source ../loadCPU.sh
source ./inputs_gas.sh

export experiment_name="${experiment_name}-${sacado_flag}"
echo "experiment name is ${experiment_name}"
echo "experiment exe is ${exec}"
use_cloned_samples=true
# change to true if you want to output the reaction rates outputs.
verbose=false
#we will save outputs in this directory
tchem_outputs=HOST_gas
mkdir -p ${tchem_outputs}/${experiment_name}
nbatch=(1)
Nthread=(128)
nparticles=(0)
source ../runThisCPU.sh
