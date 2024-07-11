#!/bin/bash

cd lambda-1.000
cp ../../input/SYSTEM.top ../../input/SYSTEM.prmtop

srun standardstatecorrection -C  ../../input/sim.cfg -t ../../input/SYSTEM.prmtop -r  traj000000001.dcd -s 10  > ../freenrg-STD.dat

wait

cd ../

