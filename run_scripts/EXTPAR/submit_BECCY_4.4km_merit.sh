#!/bin/bash -l
#SBATCH --job-name="run_extpar"
#SBATCH --account="pr133"
#SBATCH --time=00:58:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --partition=normal
#SBATCH --constraint=gpu
#SBATCH --hint=nomultithread

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export CRAY_CUDA_MPS=1
export HDF5_USE_FILE_LOCKING=FALSE

module load daint-gpu
module load CDO
module load netcdf4-python/1.5.8-CrayGNU-21.09

cd $SCRATCH/EXTPAR/extpar/run_scripts
srun ./extpar_BECCY_4.4km_merit.sh
