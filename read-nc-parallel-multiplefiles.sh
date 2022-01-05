#!/bin/bash
#PBS -l select=20:ncpus=10:mem=4gb
# set max execution time
#PBS -l walltime=0:02:00
# imposta la coda di esecuzione
#PBS -q short_cpuQ
module load gcc91
module load mpich-3.2.1--gcc-9.1.0
module load netcdf-4.7.0--gcc-9.1.0
module load hdf5-1.10.5--gcc-9.1.0
export OMP_NUM_THREADS=5
mpirun -n 20 NetCDF-Project/read-nc-parallel-multiplefiles