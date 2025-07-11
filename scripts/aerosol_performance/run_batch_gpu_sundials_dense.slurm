#!/bin/bash
#SBATCH --job-name="gpu-sun-dense"
#SBATCH --partition=ghx4
#SBATCH --mem=80G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --constraint="projects"
#SBATCH --gpus-per-node=1
#SBATCH --gpu-bind=closest
#SBATCH --account=bdoi-dtai-gh
#SBATCH -t 24:00:00

#module load cudatoolkit/12.4
#module load nvidia/24.5
#module load craype-accel-nvidia80

echo -e "job $SLURM_JOBID is starting on `hostname`\n\n"
mkdir -p out

bash runAll_gpu_sundials_dense.sh &> "./out/GPU-timings-slurm-${SLURM_JOB_ID}.out"