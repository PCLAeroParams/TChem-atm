#!/bin/bash
#SBATCH --job-name="TChem_arrow"
#SBATCH --partition=ghx4
#SBATCH --mem=10G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --constraint="projects"
#SBATCH --gpus-per-node=1
#SBATCH --gpu-bind=closest
#SBATCH --account=bdoi-dtai-gh
#SBATCH -t 02:00:00

echo -e "job $SLURM_JOBID is starting on `hostname`\n\n"

export OMP_NUM_THREADS=1
export OMP_PLACES=threads
export OMP_PROC_BIND=close

exec=$TCHEM_DEVICE_INSTALL_PATH/examples/TChem_ArrowSolver.x

run_this="$exec"

eval $run_this &> output.txt
