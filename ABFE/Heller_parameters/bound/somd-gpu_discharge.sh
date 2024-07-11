#!/bin/bash
#SBATCH -o somd-array-gpu-%A.%a.out
#SBATCH -p GTX
#SBATCH -n 1
#SBATCH --gres=gpu:1
#SBATCH --time 48:00:00
#SBATCH --array=0-8

module load cuda/7.5
module load sire/17.1.0_no_avx
echo "CUDA DEVICES:" $CUDA_VISIBLE_DEVICES

lamvals=( 0.000 0.125 0.250 0.375 0.500 0.625 0.750 0.875 1.000 )
lam=${lamvals[SLURM_ARRAY_TASK_ID]}

if [ "$SLURM_ARRAY_TASK_ID" -eq "0" ]
then
    cat ../input/distres >> ../input/sim.cfg
fi

sleep 5 


echo "lambda is: " $lam

mkdir lambda-$lam
cd lambda-$lam

export OPENMM_PLUGIN_DIR=/home/julien/sire.app/lib/plugins/

srun somd-freenrg -C ../../input/sim.cfg -l $lam -p CUDA
cd ..

wait


if [ "$SLURM_ARRAY_TASK_ID" -eq "8" ]
then
    sleep 60
    sbatch ../mbar.sh
fi

