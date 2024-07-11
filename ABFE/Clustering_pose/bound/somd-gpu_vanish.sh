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
#SBATCH --array=0-15

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

lamvals=( 0.000 0.050 0.100 0.150 0.200 0.250 0.300 0.350 0.400 0.450 0.500 0.550 0.600 0.700 0.850 1.000 )
lam=${lamvals[SLURM_ARRAY_TASK_ID]}

if [ "$SLURM_ARRAY_TASK_ID" -eq "0" ]
then
    cat ../input/distres >> ../input/sim.cfg
fi

sleep 5 

echo "lambda is: " $lam

mkdir lambda-$lam
cd lambda-$lam

srun somd-freenrg -C ../../input/sim.cfg -l $lam -p CUDA
cd ..

wait

if [ "$SLURM_ARRAY_TASK_ID" -eq "15" ]
then
    sbatch ../ljcor.sh
fi


if [ "$SLURM_ARRAY_TASK_ID" -eq "15" ]
then
    sleep 60
    sbatch ../mbar.sh
fi

