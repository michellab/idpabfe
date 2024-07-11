#!/bin/bash

export OPENMM_PLUGIN_DIR=/home/michalis/sire.app/lib/plugins

cd lambda-0.000
~/sire.app/bin/lj-tailcorrection -C ../../input/sim.cfg -l 0.00 -r traj000000001.dcd -s 100 1> ../freenrg-LJCOR-lam-0.000.dat 2> /dev/null

wait

cd ..
cd lambda-1.000
~/sire.app/bin/lj-tailcorrection -C ../../input/sim.cfg -l 1.00 -r traj000000001.dcd -s 100 1> ../freenrg-LJCOR-lam-1.000.dat 2> /dev/null

wait
cd ..

wait

# utility script to get final LJ correction term
python parselj.py freenrg-LJCOR-lam-0.000.dat freenrg-LJCOR-lam-1.000.dat > freenrg-LJCOR.dat
