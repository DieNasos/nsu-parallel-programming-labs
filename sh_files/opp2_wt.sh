#!/bin/bash
#PBS -l select=1:ncpus=1:mpiprocs=1
#PBS -l walltime=00:03:00
./opp2_wt/opp2_wt_file 1024 3 2
