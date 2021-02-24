#!/bin/bash
#SBATCH --qos=debug
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --constraint=haswell
#SBATCH --time=30:00

### Call this with a command line such as 'sbatch --array=0-49 born_job.sh'

source ~/.bashrc
cd ~/ucapnje/born_raytrace/nersc
python run_pkd_grav.py ${SLURM_ARRAY_TASK_ID}