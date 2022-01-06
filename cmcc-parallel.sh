#!/bin/bash
#PBS -l select=1:ncpus=2:ompthreads=1:mem=4gb
# set max execution time
#PBS -l walltime=0:03:00
# imposta la coda di esecuzione
#PBS -q short_cpuQ
module load gcc91
module load mpich-3.2.1--gcc-9.1.0
module load netcdf-4.7.0--gcc-9.1.0
module load hdf5-1.10.5--gcc-9.1.0
#export OMP_NUM_THREADS=5
mpirun -n 2 NetCDF-Project/cmcc-parallel

