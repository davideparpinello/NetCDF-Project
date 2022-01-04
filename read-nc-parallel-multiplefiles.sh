#!/bin/bash
#PBS -l select=5:ncpus=10:mem=4gb
# set max execution time
#PBS -l walltime=0:02:00
# imposta la coda di esecuzione
#PBS -q short_cpuQ
module load mpich-3.2
module load netcdf-4.7.0--gcc-9.1.0
module load hdf5-1.10.5--gcc-9.1.0
mpirun.actual -n 1 NetCDF-Project/read-nc-parallel-multiplefiles