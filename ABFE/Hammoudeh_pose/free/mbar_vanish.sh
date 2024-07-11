#!/bin/bash
#SBATCH -o analyse-free-nrg-%A.%a.out
#SBATCH -p serial
#SBATCH -n 1
#SBATCH --time 00:05:00
module load sire/16.1.0_no_avx
srun analyse_freenrg_mbar -i lambda-0.000/simfile.dat lambda-0.100/simfile.dat lambda-0.200/simfile.dat lambda-0.300/simfile.dat lambda-0.400/simfile.dat lambda-0.500/simfile.dat lambda-0.600/simfile.dat lambda-0.700/simfile.dat lambda-0.800/simfile.dat lambda-0.900/simfile.dat lambda-0.950/simfile.dat lambda-1.000/simfile.dat --lam 0.000 0.100 0.200 0.300 0.400 0.500 0.600 0.700 0.800 0.900 0.950 1.000 --temperature 298.0 -o mbar.pmf --subsampling percentage --percentage 82 > freenrg-MBAR.dat

