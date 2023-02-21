/******************************************************************************
 *
 * File:        reader_h8.c        
 *
 * Created:     18/01/2017
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description: Reader for pre-preprocessed SST from Himawari-8.
 *
 * Revisions:  
 *
 *****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "ncw.h"
#include "definitions.h"
#include "utils.h"
#include "obsprm.h"
#include "model.h"
#include "grid.h"
#include "observations.h"
#include "prep_utils.h"
#include "allreaders.h"

#define ERRORSTD_DEF 0.5

/**
 */
void reader_h8(char* fname, int fid, obsmeta* meta, grid* g, observations* obs)
{
    int ksurf = grid_getsurflayerid(g);
    int ncid;
    int is1d = -1;
    char llfname[MAXSTRLEN] = "";
    int dimid_ni, dimid_nj;
    size_t ni, nj;
    int varid_lon, varid_lat, varid_sst, varid_timemin, varid_timemax;
    double* lon;
    double* lat;
    double* sst;
    double* time;
    double* time2;
    char tunits[MAXSTRLEN];
    double tunits_multiple, tunits_offset;
    int i, nobs_read;

    ncw_open(fname, NC_NOWRITE, &ncid);
    if (ncw_dim_exists(ncid, "x") && ncw_dim_exists(ncid, "y")) {
        is1d = 0;
        enkf_printf("        structure = 2D\n");
    } else if (ncw_var_exists(ncid, "index")) {
        is1d = 1;
        enkf_printf("        structure = 1D\n");
    } else
        enkf_quit("%s: structure not recognised (1d or 2d)", fname);
    ncw_close(ncid);

    for (i = 0; i < meta->npars; ++i)
        if (strcasecmp(meta->pars[i].name, "LLFNAME") == 0)
            strncpy(llfname, meta->pars[i].value, MAXSTRLEN - 1);
        else
            enkf_quit("unknown PARAMETER \"%s\"\n", meta->pars[i].name);
    if (!is1d) {
        if (strlen(llfname) == 0)
            enkf_quit("%s: reader_h8: mandatory parameter \"LLFNAME\" not specified for a 2d input file\n", meta->prmfname);
    }

    if (!is1d) {
        ncw_open(llfname, NC_NOWRITE, &ncid);
        ncw_inq_dimid(ncid, "x", &dimid_ni);
        ncw_inq_dimlen(ncid, dimid_ni, &ni);
        ncw_inq_dimid(ncid, "y", &dimid_nj);
        ncw_inq_dimlen(ncid, dimid_nj, &nj);
        enkf_printf("        (ni, nj) = (%u, %u)\n", ni, nj);
        ni *= nj;

        lon = malloc(ni * sizeof(double));
        lat = malloc(ni * sizeof(double));
        ncw_inq_varid(ncid, "lon", &varid_lon);
        ncw_inq_varid(ncid, "lat", &varid_lat);
        ncw_get_var_double(ncid, varid_lon, lon);
        ncw_get_var_double(ncid, varid_lat, lat);
        ncw_close(ncid);
    }

    ncw_open(fname, NC_NOWRITE, &ncid);

    if (is1d) {
        ncw_inq_dimid(ncid, "nobs", &dimid_ni);
        ncw_inq_dimlen(ncid, dimid_ni, &ni);
        enkf_printf("        nobs total = %u\n", ni);
        lon = malloc(ni * sizeof(double));
        lat = malloc(ni * sizeof(double));
        ncw_inq_varid(ncid, "lon", &varid_lon);
        ncw_inq_varid(ncid, "lat", &varid_lat);
        ncw_get_var_double(ncid, varid_lon, lon);
        ncw_get_var_double(ncid, varid_lat, lat);
    }

    ncw_inq_varid(ncid, "sst", &varid_sst);
    sst = malloc(ni * sizeof(double));
    ncw_get_var_double(ncid, varid_sst, sst);

    ncw_inq_varid(ncid, "time_min", &varid_timemin);
    time = malloc(ni * sizeof(double));
    ncw_get_var_double(ncid, varid_timemin, time);

    ncw_inq_varid(ncid, "time_max", &varid_timemax);
    time2 = malloc(ni * sizeof(double));
    ncw_get_var_double(ncid, varid_timemax, time2);

    ncw_get_att_text(ncid, varid_timemin, "units", tunits);
    ncw_close(ncid);
    tunits_convert(tunits, &tunits_multiple, &tunits_offset);

    nobs_read = 0;
    for (i = 0; i < (int) ni; ++i) {
        observation* o;

        if (!isfinite(sst[i]) || fabs(sst[i]) > MAXOBSVAL)
            continue;

        nobs_read++;
        obs_checkalloc(obs);
        o = &obs->data[obs->nobs];

        o->product = st_findindexbystring(obs->products, meta->product);
        assert(o->product >= 0);
        o->type = obstype_getid(obs->nobstypes, obs->obstypes, meta->type, 1);
        o->instrument = st_add_ifabsent(obs->instruments, "H8", -1);
        o->id = obs->nobs;
        o->fid = fid;
        o->batch = 0;
        /*
         * add 0.17 for for the difference between skin and foundation
         * temperature
         */
        o->value = sst[i] - 273.15 + 0.17;
        o->estd = ERRORSTD_DEF;
        o->lon = lon[i];
        o->lat = lat[i];
        o->depth = 0.0;
        o->fk = (double) ksurf;
        o->status = grid_xy2fij(g, o->lon, o->lat, o->fij);
        if (!obs->allobs && o->status == STATUS_OUTSIDEGRID)
            continue;
        o->model_depth = NAN;   /* set in obs_add() */
        o->time = (time[i] + time2[i]) / 2.0 * tunits_multiple + tunits_offset;
        o->aux = -1;

        obs->nobs++;
    }
    enkf_printf("        nobs = %d\n", nobs_read);

    free(lon);
    free(lat);
    free(sst);
    free(time);
    free(time2);
}

/**
 */
void reader_h8_describe(void)
{
    enkf_printf("\n  Reader \"h8\" reads SST data from Himawari-8 imager as pre-processed\n\
in-house by BoM.\n\
\n\
  There are a number of parameters that must (marked below with \"++\"), can\n\
  (\"+\"), or may (\"-\") be specified in the corresponding section of the\n\
  observation data parameter file. The names in brackets represent the default\n\
  names checked in the abscence of the entry for the parameter. Each parameter\n\
  needs to be entered as follows:\n\
    PARAMETER <name> = <value> ...\n\
\n\
  Parameters specific to the reader:\n\
    - LLFNAME (++)\n\
        file with geographic coordinates of the data grid\n");
    describe_commonreaderparams();
}
