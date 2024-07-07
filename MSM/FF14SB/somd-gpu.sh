#!/bin/bash

####################################
#     ARIS slurm script template   #
#                                  #
# Submit script: sbatch filename   #
#                                  #
####################################

#SBATCH --job-name=soMD    # Job name
#SBATCH --output=log.%j.out # Stdout (%j expands to jobId)
#SBATCH --error=log.%j.err # Stderr (%j expands to jobId)
#SBATCH --ntasks=1     # Number of tasks(processes)
#SBATCH --nodes=1     # Number of nodes requested
#SBATCH --ntasks-per-node=1     # Tasks per node
#SBATCH --cpus-per-task=1     # Threads per task
#SBATCH --time=48:00:00   # walltime
#SBATCH --mem=28G   # memory per NODE
#SBATCH --partition=gpu    # Partition
#SBATCH --gres=gpu:1
#SBATCH --account=hpce3007   # Replace with your system project

export I_MPI_FABRICS=shm:dapl

if [ x$SLURM_CPUS_PER_TASK == x ]; then
  export OMP_NUM_THREADS=1
else
  export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
fi


## LOAD MODULES ##
module purge        # clean up loaded modules
module load gnu/6.4.0
module load cuda/9.2.148
module load sire/2018.2.0

export SIRE_DONT_PHONEHOME=1
export OPENMM_CPU_THREADS=$OMP_NUM_THREADS
srun somd -C ../scripts/sim_10.cfg -c sim_restart_9.s3 -t ../input/SYSTEM.top  -p CUDA
