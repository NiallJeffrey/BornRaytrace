#!/bin/bash
#SBATCH --qos=regular
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --constraint=knl
#SBATCH --time=300:00

### Call this with a command line such as 'sbatch --array=2-285 born_job.sh'

source ~/.bashrc
cd ~/ucapnje/born_raytrace/nersc
python run_pkd_grav.py ${SLURM_ARRAY_TASK_ID} >> log${SLURM_ARRAY_TASK_ID}
