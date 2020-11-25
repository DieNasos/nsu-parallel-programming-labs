#!/bin/bash
#PBS -l select=2:ncpus=8:mpiprocs=8
#PBS -l walltime=00:05:00
cd $PBS_O_WORKDIR
MPI_NP=$(wc -l $PBS_NODEFILE | awk '{ print $1 }')
mpirun -hostfile $PBS_NODEFILE -np $MPI_NP ./opp1_w1/opp1_w1_file_2 160 1 2
