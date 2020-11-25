#!/bin/bash
#PBS -l select=1:ncpus=2:mpiprocs=2
#PBS -l walltime=00:10:00
cd $PBS_O_WORKDIR
MPI_NP=$(wc -l $PBS_NODEFILE | awk '{ print $1 }')
mpirun -hostfile $PBS_NODEFILE -np $MPI_NP ./opp5/opp5_file 1
