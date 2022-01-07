//to compile: mpicc -g -Wall -o cmcc-parallel cmcc-parallel.c -I/apps/netCDF4.7.0/include -L/apps/netCDF4.7.0/lib -lnetcdf -lm -ldl -lz -lcurl -std=gnu99 -fopenmp

// -Wall enables all compiler's warning messages
// -g default debug information
// -o  xxxxxx output name
// yyyyy.c input name
// -I, -L, -lnetcdf etc. , flags for NetCdf linking

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <netcdf.h>
#include <sys/time.h>
#include <math.h>
#include <dirent.h>
#include <mpi.h>
#include <omp.h>

/* This is the name of the NetCDF file we will write. */
#define WR_FILE_NAME "NetCDF-Project/cmcc-med-final.nc"
/* This is the name of the folder containing all the datasets that we will read */
#define FILES_FOLDER "NetCDF-Project/CMCC-CM2-SR5_historical/"

/* We are reading 3D data, a 1 x 192 x 288 lvl-lat-lon grid, with 9125 timesteps of data. */
#define NDIMS 3
#define NLAT 192
#define NLON 288
#define LAT_NAME "lat"
#define LON_NAME "lon"
#define REC_NAME "time"

/* THe number of dimensions that the output file will have */
#define NDIMSWR 2

/* Names of things. */
#define PREC_NAME "pr"
#define UNITS "units"
#define DEGREES_EAST "degrees_east"
#define DEGREES_NORTH "degrees_north"

/* For the units attributes. */
#define UNITS "units"
//#define PREC_UNITS "kg m-2 s-1"
#define PREC_UNITS "mm/day"
#define LAT_UNITS "degrees_north"
#define LON_UNITS "degrees_east"
#define MAX_ATT_LEN 80

/* Handle errors by printing an error message and exiting with a
  * non-zero status. */
#define ERR(e)                                 \
    {                                          \
        printf("Error: %s\n", nc_strerror(e)); \
        return 2;                              \
    }

/* function to measure time */
float time_diff(struct timeval *start, struct timeval *end)
{
    return (end->tv_sec - start->tv_sec) + 1e-6 * (end->tv_usec - start->tv_usec);
}

int main(int argc, char *argv[])
{
    /* MPI  inizialization */
    MPI_Init(&argc, &argv);

    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    // Print off a hello world message with processor name
    printf("greetings:  %s, rank %d out of %d processors\n",
           processor_name, rank, size);

    /* reading files from specific directory containing NetCDF matrices */
    struct dirent *file; // Pointer for directory entry

    // opendir() returns a pointer of DIR type.
    DIR *directory = opendir(FILES_FOLDER);

    if (directory == NULL) // opendir returns NULL if couldn't open directory
    {
        printf("Could not open current directory");
        return 0;
    }

    int file_count = 0;

    while ((file = readdir(directory)) != NULL)
    {
        if (file->d_type == DT_REG)
        {                 /* If the entry is a regular file */
            file_count++; //increase the counter containing the total number of files inside the dir
        }
    }

    closedir(directory);

    /* array containing all the names of the different files inside the directory */
    char *filesList[file_count];

    int counter = 0;

    directory = opendir(FILES_FOLDER);

    while ((file = readdir(directory)) != NULL)
        if (file->d_type == DT_REG)
        {
            char name[250]; //variable containing the file path
            strcpy(name, FILES_FOLDER);
            strcat(name, file->d_name);                            // append the name of file onto the folder path
            filesList[counter] = (char *)malloc(strlen(name) + 1); //allocating memory for filesList individual elements.
            strncpy(filesList[counter], name, strlen(name));
            counter++;
        }
    closedir(directory);

    /* variables for measuring time */
    struct timeval starttime;
    struct timeval starttime2;
    struct timeval starttime3;
    struct timeval endtime;
    struct timeval endtime2;
    double elapsed_time = 0;
    double total_time = 0;
    gettimeofday(&starttime2, NULL); //start timer of rank 0

    /*array of matrices containing all the different output average matrices calculated for each file*/
    float *final_averages = malloc(file_count * NLAT * NLON * sizeof(float));

    /* These program variables hold the latitudes and longitudes. */
    float lats[NLAT], lons[NLON];

    /* This array will hold the number of records for each iteration, it will be used to calculate the final averages */
    int nrec_array[file_count];

    double total_time_reading;
    double total_time_elaboration;

    /* loop for each file (numbered by counter); every process will open the same file and perform the reading and the calculation; 
        the next iteration the next file and so on */
    for (counter = 0; counter < file_count; counter++)
    {
        int nrec;
        int ncid, prec_varid;
        int lat_varid, lon_varid;

        /* variable containing the name of the current file to be opened and read */
        char filename[250];

        /* copy the name of the current file, contained in the array filesList, to the variable filename */
        strcpy(filename, filesList[counter]);

        /* Last file has a duration of 15 years and not 25 */
        if (strcmp(filename, "NetCDF-Project/CMCC-CM2-SR5_historical/pr_day_CMCC-CM2-SR5_historical_r1i1p1f1_gn_20000101-20141231.nc") != 0)
            nrec = 9125;
        else
            nrec = 5475;

        nrec_array[counter] = nrec;

        /* setup of NetCDF reading */

        /* The start and count arrays will tell the netCDF library where to read our data. */
        size_t start[NDIMS], count[NDIMS];

        /* Program variables to hold the data we will read. We will only need enough space to hold one timestep of data; one record. */
        float prec_in[NLAT][NLON];

        /* Loop indexes. */
        int rec = 0;

        /* Error handling. */
        int retval;

        /* Open the file. */
        if ((retval = nc_open(filename, NC_NOWRITE, &ncid)))
            ERR(retval);

        /* Get the varids of the latitude and longitude coordinate variables. */
        if ((retval = nc_inq_varid(ncid, LAT_NAME, &lat_varid)))
            ERR(retval);
        if ((retval = nc_inq_varid(ncid, LON_NAME, &lon_varid)))
            ERR(retval);

        /* Read the coordinate variable data. */
        if ((retval = nc_get_var_float(ncid, lat_varid, &lats[0])))
            ERR(retval);
        if ((retval = nc_get_var_float(ncid, lon_varid, &lons[0])))
            ERR(retval);

        /* Get the varid of the precipitation netCDF variable. */
        if ((retval = nc_inq_varid(ncid, PREC_NAME, &prec_varid)))
            ERR(retval);

        /* Read the data. Since we know the contents of the file we know that the data arrays in this program are the correct size to
            hold one timestep. */

        count[0] = 1;
        count[1] = NLAT;
        count[2] = NLON;
        start[1] = 0;
        start[2] = 0;

        /* end of setup of NetCDF reading */

        double elapsed_time_reading = 0; //variable containing the local sum of all elapsed time for reading
        double elapsed_time_matrix = 0;  //variable containing the local sum of all elapsed time for writing the sum matrix

        /* loop iteration variables */
        int i, k;

        /* sum matrix */
        float sum[NLAT][NLON] = {{0}};

        /* if 25 processes, 365 days elaborated by each proc */
        int days_per_proc = ceil((double)nrec / size);

        int limit = (rank + 1) * days_per_proc;

        if (limit > nrec)
        {
            limit = nrec;
        }

        int threadnumb;
        /* Get the total number of active threads in the program */
        #pragma omp parallel {
            threadnumb = omp_get_num_threads();
        }
        /* only master process writes general info about current iteration */
        if (rank == 0)
        {
            printf("Elaborated file: %s (counter: %d) \n", filename, counter);
            printf("Number of processes: %d (days elaborated by each process: %d)\n", size, days_per_proc);
            printf("Number of threads: %d \n", threadnumb);
            gettimeofday(&starttime3, NULL); //start timer of rank0
        }

        /* Read and check one record at a time.
            if 125 processes, rank 0 from day 0 to day 72, rank 1 from 73 to 145 ... */
        for (rec = rank * days_per_proc; rec < limit; rec++)
        {
            gettimeofday(&starttime, NULL); //start reading timer

            start[0] = rec; //starting day to elaborate

            /*read the information */
            if ((retval = nc_get_vara_float(ncid, prec_varid, start, count, &prec_in[0][0])))
                ERR(retval);

            /*end of reading timer and calculation of elapsed time */
            gettimeofday(&endtime, NULL);
            elapsed_time = time_diff(&starttime, &endtime);
            elapsed_time_reading += elapsed_time;

            /* start elaboration timer */
            gettimeofday(&starttime, NULL);

/* population of sum matrix */
#pragma omp parallel for collapse(2) private(i, k) reduction(+ \
                                                             : sum) schedule(guided)
            for (i = 0; i < NLAT; i++)
            {
                for (k = 0; k < NLON; k++)
                {
                    sum[i][k] += prec_in[i][k];
                }
            }

            /*end of elaboration timer and calculation of elapsed time */
            gettimeofday(&endtime, NULL);
            elapsed_time = time_diff(&starttime, &endtime);
            elapsed_time_matrix += elapsed_time;
        } /* next record */

        //pointer to sum matrix, to consider the matrix as an array.
        int *array_matr = (int *)sum;

        double sum_readtimes;        // variable containing the output of reduce operation, collecting all reading times of different processes
        double sum_elaborationtimes; // variable containing the output of reduce operation, collecting all elaboration times of different processes

        /* three MPI_Reduce, to collect the local contibutes of reading timer, elaboration timer, and most importantly the local sum matrix. */
        MPI_Reduce(&elapsed_time_reading, &sum_readtimes, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&elapsed_time_matrix, &sum_elaborationtimes, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

        gettimeofday(&starttime, NULL); // start communication timer
        MPI_Reduce(array_matr, &final_averages[counter * NLAT * NLON], NLAT * NLON, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
        /*end of communication timer and calculation of elapsed time */
        gettimeofday(&endtime, NULL);
        elapsed_time = time_diff(&starttime, &endtime);

        if (rank == 0)
        {
            total_time_reading += sum_readtimes / size;
            total_time_elaboration += sum_elaborationtimes / size;
            printf("Average reading time per process: %.7f\n", sum_readtimes / size);
            printf("Average elaboration of sum matrix per process: %.7f\n\n", sum_elaborationtimes / size);
            gettimeofday(&endtime, NULL);
            elapsed_time = time_diff(&starttime3, &endtime);
            printf("Elaboration time of master process: %.7f\n", elapsed_time);
            printf("*********************************************\n\n");
        }

        if ((retval = nc_close(ncid)))
            ERR(retval);
    }

    if (rank == 0)
    {
        // Master average matrix to calculate the final average from the final_averages matrix
        float *master_average = calloc(NLAT * NLON, sizeof(float));

        // From final_averages matrix, calculation of the master average
        int i, j;
        for (i = 0; i < file_count; i++)
        {
            float *average = final_averages + i * (NLAT * NLON);
            for (j = 0; j < NLAT * NLON; j++)
            {
                master_average[j] += average[j] * 86400 / nrec_array[i] / file_count;
            }
        }

        /* Writing of average matrix into a new nc file */
        int ncid_wr, prec_varid_wr;
        int lat_varid_wr, lon_varid_wr, lon_dimid, lat_dimid;

        int dimids_wr[NDIMSWR];

        /* Error handling. */
        int retval;

        /* Create the file. */
        if ((retval = nc_create(WR_FILE_NAME, NC_CLOBBER, &ncid_wr)))
            ERR(retval);

        /* Define the dimensions. */
        if ((retval = nc_def_dim(ncid_wr, LAT_NAME, NLAT, &lat_dimid)))
            ERR(retval);
        if ((retval = nc_def_dim(ncid_wr, LON_NAME, NLON, &lon_dimid)))
            ERR(retval);

        /* Define coordinate netCDF variables. They will hold the coordinate information, that is, the latitudes and longitudes. 
            A varid is returned for each.*/
        if ((retval = nc_def_var(ncid_wr, LAT_NAME, NC_FLOAT, 1, &lat_dimid, &lat_varid_wr)))
            ERR(retval);
        if ((retval = nc_def_var(ncid_wr, LON_NAME, NC_FLOAT, 1, &lon_dimid, &lon_varid_wr)))
            ERR(retval);

        /* Define units attributes for coordinate vars. This attaches a text attribute to each of the coordinate variables, containing
        the units. Note that we are not writing a trailing NULL, just "units", because the reading program may be fortran which does
        not use null-terminated strings. In general it is up to the reading C program to ensure that it puts null-terminators on strings where necessary.*/
        if ((retval = nc_put_att_text(ncid_wr, lat_varid_wr, UNITS, strlen(DEGREES_NORTH), DEGREES_NORTH)))
            ERR(retval);
        if ((retval = nc_put_att_text(ncid_wr, lon_varid_wr, UNITS, strlen(DEGREES_EAST), DEGREES_EAST)))
            ERR(retval);

        /* Define the netCDF variables. The dimids array is used to pass the dimids of the dimensions of the variables.*/
        dimids_wr[0] = lat_dimid;
        dimids_wr[1] = lon_dimid;
        if ((retval = nc_def_var(ncid_wr, PREC_NAME, NC_FLOAT, NDIMSWR, dimids_wr, &prec_varid_wr)))
            ERR(retval);

        /* Define units attributes for vars. */
        if ((retval = nc_put_att_text(ncid_wr, prec_varid_wr, UNITS, strlen(PREC_UNITS), PREC_UNITS)))
            ERR(retval);

        /* End define mode. */
        if ((retval = nc_enddef(ncid_wr)))
            ERR(retval);

        gettimeofday(&starttime, NULL);
        /* Write the coordinate variable data. This will put the latitudes and longitudes of our data grid into the netCDF file. */
        if ((retval = nc_put_var_float(ncid_wr, lat_varid_wr, &lats[0])))
            ERR(retval);
        if ((retval = nc_put_var_float(ncid_wr, lon_varid_wr, &lons[0])))
            ERR(retval);

        /* Write the pretend data. This will write our surface pressure and surface temperature data. The arrays of data are the same size
        as the netCDF variables we have defined. */
        if ((retval = nc_put_var_float(ncid_wr, prec_varid_wr, master_average)))
            ERR(retval);

        /* Close the file. */
        if ((retval = nc_close(ncid_wr)))
            ERR(retval);

        gettimeofday(&endtime, NULL);
        elapsed_time = time_diff(&starttime, &endtime);
        printf("Writing NetCDFfile misuration time : %.7f\n\n", elapsed_time);

        printf("Total average reading time: %.7f\n\n", total_time_reading / file_count);
        printf("Total average elaboration time: %.7f\n\n", total_time_elaboration / file_count);
        gettimeofday(&endtime2, NULL);
        total_time = time_diff(&starttime2, &endtime2);
        printf("Total elaboration time of master process: %.7f\n\n", total_time);
    }
    MPI_Finalize();
    return 0;
}