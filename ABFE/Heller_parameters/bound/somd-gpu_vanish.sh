#!/bin/bash
#SBATCH -o somd-array-gpu-%A.%a.out
#SBATCH -p GTX
#SBATCH -n 1
#SBATCH --gres=gpu:1
#SBATCH --time 48:00:00
#SBATCH --array=0-15

module load cuda/7.5
module load  sire/17.1.0_no_avx
echo "CUDA DEVICES:" $CUDA_VISIBLE_DEVICES

lamvals=( 0.000 0.050 0.100 0.150 0.200 0.250 0.300 0.350 0.400 0.450 0.500 0.550 0.600 0.850 1.000 )
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

if [ "$SLURM_ARRAY_TASK_ID" -eq "11" ]
then
    sbatch ../ljcor.sh
fi

if [ "$SLURM_ARRAY_TASK_ID" -eq "11" ]
then
    sleep 60
    sbatch ../mbar.sh
fi
