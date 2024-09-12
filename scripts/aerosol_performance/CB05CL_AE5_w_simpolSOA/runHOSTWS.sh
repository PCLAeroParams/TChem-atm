export DEVICE=CPU
export experiment_exe_root=${experiment_base}${experiment_suffix}
export sacado_flag=$sflag
echo "sacado flag = ${sacado_flag}"
source ../loadCPU_WS.sh
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
#20 100 1000 10000 100000 200000 500000 1000000 2000000 3000000 4500000
nbatch=(1)
# 144 216 288)
# 360 432 504 576 648 720 792 864 936 1008)
Nthread=(52 104)
nparticles=(1 10 50 100)
# 500 1000 5000 10000
source ../runThisCPU.sh
