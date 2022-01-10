# NetCDF-Project
Repository for the High Performance Computing for Data Science's Project, parallel analysis of NetCDF files.

## Set-up
### Dataset folders

The code is configured to work with the datasets' folders copied locally in the project folder (under the user home directory). 
Instead, to use the datasets under the shared folder of the cluster, please change both rows 22 inside the C source files.

### Modules

It is necessary to load the following modules to let the code work properly:

```shell
module load gcc91
module load mpich-3.2.1--gcc-9.1.0
module load netcdf-4.7.0--gcc-9.1.0
module load hdf5-1.10.5--gcc-9.1.0
```

## Compile

To compile the C source files, use the following commands:

For CMCC dataset:

```shell
mpicc -g -Wall -o cmcc-parallel cmcc-parallel.c -I/apps/netCDF4.7.0/include -L/apps/netCDF4.7.0/lib -lnetcdf -lm -ldl -lz -lcurl -std=gnu99 -fopenmp
```

For MPI-ESM dataset:

```shell
mpicc -fopenmp -g -Wall -o mpi-esm-parallel mpi-esm-parallel.c -I/apps/netCDF4.7.0/include -L/apps/netCDF4.7.0/lib -lnetcdf -lm -ldl -lz -lcurl -std=gnu99
```

## Run
Submit the jobs to the cluster queue using the two `.sh` files

```shell
qsub cmcc-parallel.sh
qsub mpi-esm-parallel.sh
```

## Output
The shell output files will contain the execution and elaboration times. Instead, the files `cmcc-med-final.nc` and `esm-med-final.nc`, should be copied locally on your PC and analyzed with Panoply or similar.