/* Copyright 2006-2011 University Corporation for Atmospheric
 Research/Unidata. See COPYRIGHT file for conditions of use. */
#include <stdio.h>
#include <string.h>
#include <netcdf.h>
#include <time.h>

/* This is the name of the data file we will read. */
#define FILE_NAME "NetCDF-Project/data.nc"

/* We are reading 3D data, a 1 x 192 x 288 lvl-lat-lon grid, with 9125
    timesteps of data. */
#define NDIMS 3
#define NLAT 192
#define NLON 288
#define LAT_NAME "lat"
#define LON_NAME "lon"
#define NREC 5475
#define REC_NAME "time"

/* Names of things. */
#define PREC_NAME "pr"
#define UNITS "units"
#define DEGREES_EAST "degrees_east"
#define DEGREES_NORTH "degrees_north"

/* These are used to calculate the values we expect to find. */
#define SAMPLE_PRECIPITATION 0.0
#define START_LAT 25.0
#define START_LON -125.0

/* For the units attributes. */
#define UNITS "units"
#define PRES_UNITS "kg m-2 s-1"
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

int main()
{
   int ncid, prec_varid;
   int lat_varid, lon_varid;

   time_t start_t, end_t;
   double diff_t;

   time(&start_t);

   /* The start and count arrays will tell the netCDF library where to
       read our data. */
   size_t start[NDIMS], count[NDIMS];

   /* Program variables to hold the data we will read. We will only
       need enough space to hold one timestep of data; one record. */
   float prec_in[NLAT][NLON];

   /* These program variables hold the latitudes and longitudes. */
   float lats[NLAT], lons[NLON];

   /* Loop indexes. */
   int lat, lon, rec = 0;

   /* Error handling. */
   int retval;

   /* Open the file. */
   if ((retval = nc_open(FILE_NAME, NC_NOWRITE, &ncid)))
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

   /* Check the coordinate variable data. */
   /* for (lat = 0; lat < NLAT; lat++)
       if (lats[lat] != START_LAT + 5.*lat)
      return 2;
    for (lon = 0; lon < NLON; lon++)
       if (lons[lon] != START_LON + 5.*lon)
      return 2; */

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
   start[1] = 3;
   start[2] = 150;

   float data, somma = 0;

   /* Read and check one record at a time. */
   for (rec = 0; rec < NREC; rec++)
   {
      start[0] = rec;

      if ((retval = nc_get_var1_float(ncid, prec_varid, start, &data)))
         ERR(retval);

      printf("Day: %d Data %.7f\n", rec, data);

      somma += data;

   } /* next record */

   somma = somma/NREC;
   printf("\n\n--- MEDIA: %.7f\n", somma);

   time(&end_t);

   diff_t = difftime(end_t, start_t);

   printf("Execution time = %.3f\n", diff_t);

   /* Close the file. */
   if ((retval = nc_close(ncid)))
      ERR(retval);

   

   //printf("*** SUM = %.3f\n", sum);
   return 0;
}