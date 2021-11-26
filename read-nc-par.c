/* Copyright 2006-2011 University Corporation for Atmospheric
 Research/Unidata. See COPYRIGHT file for conditions of use. */
#include <stdio.h>
#include <string.h>
#include <netcdf.h>
#include <sys/time.h>

/* This is the name of the data file we will read. */
#define FILE_NAME "/shares/HPC4DataScience/pta/CMCC-CM2-SR5_historical/pr_day_CMCC-CM2-SR5_historical_r1i1p1f1_gn_19750101-19991231.nc"
#define WR_FILE_NAME "NetCDF-Project/cmcc-med.nc"

/* We are reading 3D data, a 1 x 192 x 288 lvl-lat-lon grid, with 9125
    timesteps of data. */
#define NDIMS 3
#define NLAT 192
#define NLON 288
#define LAT_NAME "lat"
#define LON_NAME "lon"
#define NREC 5475
#define REC_NAME "time"

#define NDIMSWR 2

/* Names of things. */
#define PREC_NAME "pr"
#define UNITS "units"
#define DEGREES_EAST "degrees_east"
#define DEGREES_NORTH "degrees_north"

/* For the units attributes. */
#define UNITS "units"
#define PREC_UNITS "kg m-2 s-1"
#define LAT_UNITS "degrees_north"
#define LON_UNITS "degrees_east"
#define MAX_ATT_LEN 80

/* Handle errors by printing an error message and exiting with a
  * non-zero status. */
#define ERR(e)                               \
   {                                         \
      printf("Error: %s\n", nc_strerror(e)); \
      return 2;                              \
   }

float time_diff(struct timeval *start, struct timeval *end)
{
    return (end->tv_sec - start->tv_sec) + 1e-6*(end->tv_usec - start->tv_usec);
}

int main()
{
   int ncid, omode, ncid_wr, prec_varid, prec_varid_wr;
   int lat_varid, lon_varid, lat_varid_wr, lon_varid_wr, lon_dimid, lat_dimid;

   int dimids_wr[NDIMSWR];

   omode = NC_NOWRITE;

   /* The start and count arrays will tell the netCDF library where to
       read our data. */
   size_t start[NDIMS], count[NDIMS];

   /* Program variables to hold the data we will read. We will only
       need enough space to hold one timestep of data; one record. */
   float prec_in[NLAT][NLON];

   /* These program variables hold the latitudes and longitudes. */
   float lats[NLAT], lons[NLON];

   /* Loop indexes. */
   int rec = 0;

   /* Error handling. */
   int retval;

   /* Open the file. */
   if ((retval = nc_open_par(FILE_NAME, omode, MPI_COMM_WORLD, MPI_INFO_NULL, &ncid)))
      ERR(retval);

   /* Create the file. */
   if ((retval = nc_create(WR_FILE_NAME, NC_CLOBBER, &ncid_wr)))
      ERR(retval);

   /* Define the dimensions. */
   if ((retval = nc_def_dim(ncid_wr, LAT_NAME, NLAT, &lat_dimid)))
      ERR(retval);
   if ((retval = nc_def_dim(ncid_wr, LON_NAME, NLON, &lon_dimid)))
      ERR(retval);

   /* Define coordinate netCDF variables. They will hold the
      coordinate information, that is, the latitudes and longitudes. A
      varid is returned for each.*/
   if ((retval = nc_def_var(ncid_wr, LAT_NAME, NC_FLOAT, 1, &lat_dimid, 
			    &lat_varid_wr)))
      ERR(retval);
   if ((retval = nc_def_var(ncid_wr, LON_NAME, NC_FLOAT, 1, &lon_dimid, 
			    &lon_varid_wr)))
      ERR(retval);

   /* Define units attributes for coordinate vars. This attaches a
      text attribute to each of the coordinate variables, containing
      the units. Note that we are not writing a trailing NULL, just
      "units", because the reading program may be fortran which does
      not use null-terminated strings. In general it is up to the
      reading C program to ensure that it puts null-terminators on
      strings where necessary.*/
   if ((retval = nc_put_att_text(ncid_wr, lat_varid_wr, UNITS, 
				 strlen(DEGREES_NORTH), DEGREES_NORTH)))
      ERR(retval);
   if ((retval = nc_put_att_text(ncid_wr, lon_varid_wr, UNITS, 
				 strlen(DEGREES_EAST), DEGREES_EAST)))
      ERR(retval);

   /* Define the netCDF variables. The dimids array is used to pass
      the dimids of the dimensions of the variables.*/
   dimids_wr[0] = lat_dimid;
   dimids_wr[1] = lon_dimid;
   if ((retval = nc_def_var(ncid_wr, PREC_NAME, NC_FLOAT, NDIMSWR, 
			    dimids_wr, &prec_varid_wr)))
      ERR(retval);

      /* Define units attributes for vars. */
   if ((retval = nc_put_att_text(ncid_wr, prec_varid_wr, UNITS, 
				 strlen(PREC_UNITS), PREC_UNITS)))
      ERR(retval);

   /* End define mode. */
   if ((retval = nc_enddef(ncid_wr)))
      ERR(retval);

   /* Get the varids of the latitude and longitude coordinate
     * variables. */
   if ((retval = nc_inq_varid(ncid, LAT_NAME, &lat_varid)))
      ERR(retval);
   if ((retval = nc_inq_varid(ncid, LON_NAME, &lon_varid)))
      ERR(retval);

   /* Read the coordinate variable data. */
   if ((retval = nc_get_var_float(ncid, lat_varid, &lats[0])))
      ERR(retval);
   if ((retval = nc_get_var_float(ncid, lon_varid, &lons[0])))
      ERR(retval);

   /* Get the varids of the pressure and temperature netCDF
     * variables. */
   if ((retval = nc_inq_varid(ncid, PREC_NAME, &prec_varid)))
      ERR(retval);

   /* Read the data. Since we know the contents of the file we know
     * that the data arrays in this program are the correct size to
     * hold one timestep. */
   count[0] = 1;
   count[1] = NLAT;
   count[2] = NLON;
   start[1] = 0;
   start[2] = 0;

   struct timeval starttime;
   struct timeval endtime;
   double elapsed_time;

   int i, k;

   float somma[NLAT][NLON] = {{0}};

   gettimeofday(&starttime, NULL);
   /* Read and check one record at a time. */
   
   for (rec = 0; rec < NREC; rec++)
   {
      start[0] = rec;

      if ((retval = nc_get_vara_float(ncid, prec_varid, start,
                                      count, &prec_in[0][0])))
         ERR(retval);

      for(i = 0; i < NLAT; i++) {
         for (k = 0; k < NLON; k++) {
            somma[i][k] += prec_in[i][k];
         }
      }

   } /* next record */

   printf("Here\n\n");
   
   
   for(i = 0; i < NLAT; i++) {
      for (k = 0; k < NLON; k++) {
         somma[i][k] = (somma[i][k] / NREC) * 86400;
      }
   }

   

   gettimeofday(&endtime, NULL);
   elapsed_time = time_diff(&starttime, &endtime);
   printf("Tempo medio: %.7f\n\n", elapsed_time);

   printf("Here\n\n");
   for(i = 0; i < NLAT; i++) {
      for (k = 0; k < NLON; k++) {
         printf("%.7f \t", somma[i][k]);
      }
      fflush(stdout); 
      printf("\n");
   }

   printf("Here\n\n");
   

   /* Write the coordinate variable data. This will put the latitudes
      and longitudes of our data grid into the netCDF file. */
   if ((retval = nc_put_var_float(ncid_wr, lat_varid_wr, &lats[0])))
      ERR(retval);
   if ((retval = nc_put_var_float(ncid_wr, lon_varid_wr, &lons[0])))
      ERR(retval);

   printf("Here\n\n");

   /* Write the pretend data. This will write our surface pressure and
      surface temperature data. The arrays of data are the same size
      as the netCDF variables we have defined. */
   if ((retval = nc_put_var_float(ncid_wr, prec_varid_wr, &somma[0][0])))
      ERR(retval);

   /* Close the file. */
   if ((retval = nc_close(ncid)))
      ERR(retval);

      /* Close the file. */
   if ((retval = nc_close(ncid_wr)))
      ERR(retval);

   return 0;
}