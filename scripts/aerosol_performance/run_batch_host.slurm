#!/bin/bash
#SBATCH --job-name="TChem_timings"
#SBATCH --partition=cpu
#SBATCH --mem=30G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --constraint="projects"
#SBATCH --account=bdoi-dtai-gh
#SBATCH -t 24:00:00

#module load cudatoolkit/12.4
#module load nvidia/24.5
#module load craype-accel-nvidia80

echo -e "job $SLURM_JOBID is starting on `hostname`\n\n"
mkdir -p out

bash runAll_host.sh &> "./out/CPU-timings-slurm-${SLURM_JOB_ID}.out"