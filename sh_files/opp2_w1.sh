#!/bin/bash
#PBS -l select=2:ncpus=8:mpiprocs=8
#PBS -l walltime=00:05:00
cd $PBS_O_WORKDIR
MPI_NP=$(wc -l $PBS_NODEFILE | awk '{ print $1 }')
./opp2_w1/opp2_w1_file 16 1024 3 2
