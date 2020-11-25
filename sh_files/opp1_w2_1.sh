#!/bin/bash
#PBS -l select=1:ncpus=2:mpiprocs=2
#PBS -l walltime=00:03:00
cd $PBS_O_WORKDIR
MPI_NP=$(wc -l $PBS_NODEFILE | awk '{ print $1 }')
mpirun -hostfile $PBS_NODEFILE -np $MPI_NP ./opp1_w1/opp1_w1_file 16384 3 1
