export DEVICE=CPU
export experiment_exe_root=${experiment_base}${experiment_suffix}
export sacado_flag=$sflag
echo "sacado flag = ${sacado_flag}"
source ../loadCPU.sh
source ./inputs.sh

export experiment_name="${experiment_name}-${sacado_flag}"
echo "experiment name is ${experiment_name}"
echo "experiment exe is ${exec}"
use_cloned_samples=true
# change to true if you want to output the reaction rates outputs.
verbose=false
#we will save outputs in this directory
tchem_outputs=HOST
mkdir -p ${tchem_outputs}/${experiment_name}
nbatch=(1 1000)
Nthread=(128)
nparticles=(1 2 4 8 18 37 78 162 335)
source ../runThisCPU.sh
