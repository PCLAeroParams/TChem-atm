export DEVICE=GPU
export experiment_exe_root=${experiment_base}${experiment_suffix}
export sacado_flag=$sflag
echo "sacado flag = ${sacado_flag}"
source ../loadGPU_WS.sh
source ./inputs.sh

export experiment_name="${experiment_name}-${sacado_flag}"
echo "experiment name is ${experiment_name}"
echo "experiment exe is ${exec}"
use_cloned_samples=true
# change to true if you want to output the reaction rates outputs.
verbose=false
#we will save outputs in this directory
tchem_outputs=CUDA
mkdir -p ${tchem_outputs}/${experiment_name}
nbatch=(72 144 216)
# 288 360 432 504 576 648 720 792 864 936 1008)
# let's use for now 1 and 1 for team and vector size; however, we need to tune up these parameters.
vector_size=(1)
team_size=(1)
nparticles=(1 10 100 1000)
source ../runThisGPU.sh
