/******************************************************************************
 *
 * File:        reader_amsr2.c        
 *
 * Created:     30/01/2015
 *
 * Author:      Paul Sandery
 *              Bureau of Meteorology
 *
 * Description: AMSR-2 reader, created from reader_navo.c.
 *
 * Revisions:   07/08/2018   Xinmei Huang 
 *                Added parameter ADDBIAS for reverting subtraction of sses_bias
 *                during pre-processing.
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

#define ADDBIAS_DEF 0

/**
 */
void reader_amsr2(char* fname, int fid, obsmeta* meta, grid* g, observations* obs)
{
    int ksurf = grid_getsurflayerid(g);
    int addbias = ADDBIAS_DEF;
    int ncid;
    int dimid_nobs;
    size_t nobs_local;
    int varid_lon, varid_lat, varid_sst, varid_sstb, varid_error, varid_time;
    double* lon = NULL;
    double* lat = NULL;
    double* sst = NULL;
    double* sstb = NULL;
    double* error_std = NULL;
    double* time = NULL;
    int year, month, day;
    char tunits[MAXSTRLEN];
    size_t tunits_len;
    double tunits_multiple, tunits_offset;
    char* basename = NULL;
    int i;

    for (i = 0; i < meta->npars; ++i) {
        if (strcasecmp(meta->pars[i].name, "ADDBIAS") == 0)
            addbias = (istrue(meta->pars[i].value)) ? 1 : 0;
        else
            enkf_quit("unknown PARAMETER \"%s\"\n", meta->pars[i].name);
    }
    enkf_printf("        ADDBIAS = %s\n", (addbias) ? "YES" : "NO");

    basename = strrchr(fname, '/');
    if (basename == NULL)
        basename = fname;
    else
        basename += 1;

    ncw_open(fname, NC_NOWRITE, &ncid);

    ncw_inq_dimid(ncid, "nobs", &dimid_nobs);
    ncw_inq_dimlen(ncid, dimid_nobs, &nobs_local);

    if (nobs_local == 0) {
        enkf_printf("        nobs = %u\n", (unsigned int) nobs_local);
        ncw_close(ncid);
        return;
    }

    ncw_inq_varid(ncid, "lon", &varid_lon);
    lon = malloc(nobs_local * sizeof(double));
    ncw_get_var_double(ncid, varid_lon, lon);

    ncw_inq_varid(ncid, "lat", &varid_lat);
    lat = malloc(nobs_local * sizeof(double));
    ncw_get_var_double(ncid, varid_lat, lat);

    ncw_inq_varid(ncid, "sst", &varid_sst);
    sst = malloc(nobs_local * sizeof(double));
    ncw_get_var_double(ncid, varid_sst, sst);

    if (addbias) {
        if (ncw_var_exists(ncid, "sses_bias")) {
            ncw_inq_varid(ncid, "sses_bias", &varid_sstb);
            sstb = malloc(nobs_local * sizeof(double));
            ncw_get_var_double(ncid, varid_sstb, sstb);
        } else
            addbias = 0;
    }

    ncw_inq_varid(ncid, "error", &varid_error);
    error_std = malloc(nobs_local * sizeof(double));
    ncw_get_var_double(ncid, varid_error, error_std);

    ncw_inq_varid(ncid, "age", &varid_time);
    time = malloc(nobs_local * sizeof(double));
    ncw_get_var_double(ncid, varid_time, time);
    ncw_inq_attlen(ncid, varid_time, "units", &tunits_len);
    ncw_get_att_text(ncid, varid_time, "units", tunits);
    basename[14] = 0;
    if (!str2int(&basename[12], &day))
        enkf_quit("AMSR2 reader: could not convert file name \"%s\" to date", fname);
    basename[12] = 0;
    if (!str2int(&basename[10], &month))
        enkf_quit("AMSR2 reader: could not convert file name \"%s\" to date", fname);
    basename[10] = 0;
    if (!str2int(&basename[6], &year))
        enkf_quit("AMSR2 reader: could not convert file name \"%s\" to date", fname);
    snprintf(&tunits[tunits_len], MAXSTRLEN - tunits_len, " since %4d-%02d-%02d", year, month, day);

    ncw_close(ncid);

    tunits_convert(tunits, &tunits_multiple, &tunits_offset);

    for (i = 0; i < (int) nobs_local; ++i) {
        observation* o;

        obs_checkalloc(obs);
        o = &obs->data[obs->nobs];

        o->product = st_findindexbystring(obs->products, meta->product);
        assert(o->product >= 0);
        o->type = obstype_getid(obs->nobstypes, obs->obstypes, meta->type, 1);
        o->instrument = st_add_ifabsent(obs->instruments, "AMSR-2", -1);
        o->id = obs->nobs;
        o->id_orig = i;
        o->fid = fid;
        o->batch = 0;
        o->value = (addbias) ? sst[i] + sstb[i] : sst[i];
        o->estd = error_std[i];
        o->lon = lon[i];
        o->lat = lat[i];
        o->depth = 0.0;
        o->fk = (double) ksurf;
        o->status = grid_xy2fij(g, o->lon, o->lat, o->fij);
        if (!obs->allobs && o->status == STATUS_OUTSIDEGRID)
            continue;
        o->model_depth = NAN;   /* set in obs_add() */
        o->time = time[i] * tunits_multiple + tunits_offset;
        o->aux = -1;

        obs->nobs++;
    }
    enkf_printf("        nobs = %d\n", nobs_local);

    free(lon);
    free(lat);
    free(sst);
    free(error_std);
    free(time);
    if (addbias)
        free(sstb);
}

/**
 */
void reader_amsr2_describe(void)
{
    enkf_printf("\n  Reader \"amsr2\" reads SST data from AMSR2 radiometer as pre-processed\n\
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
    - ADDBIAS (-)\n\
        reverses bias correction\n");
    describe_commonreaderparams();
}
