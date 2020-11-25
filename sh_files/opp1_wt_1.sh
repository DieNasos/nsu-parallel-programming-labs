#!/bin/bash
#PBS -l select=1:ncpus=1:mpiprocs=1
#PBS -l walltime=00:03:00
./opp1_wt/opp1_wt_file 8192 3 1
