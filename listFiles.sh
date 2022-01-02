#!/bin/bash
#PBS -l select=1:ncpus=1:mem=2gb
# set max execution time
#PBS -l walltime=0:05:00
# imposta la coda di esecuzione
#PBS -q short_cpuQ
./NetCDF-Project/listFiles