/******************************************************************************
 *
 * File:        reader_windsat.c        
 *
 * Created:     12/2012
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description:
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
#include "nan.h"
#include "stringtable.h"
#include "kdtree.h"
#include "definitions.h"
#include "utils.h"
#include "enkfprm.h"
#include "obsmeta.h"
#include "model.h"
#include "grid.h"
#include "observations.h"
#include "prep.h"
#include "allreaders.h"

void reader_windsat_standard(char* fname, int fid, obsmeta* meta, model* m, observations* obs)
{
    int ncid;
    int dimid_nobs;
    size_t nobs_local;
    int varid_lon, varid_lat, varid_sst, varid_error, varid_time;
    double* lon;
    double* lat;
    double* sst;
    double* error_std;
    double* time;
    int year, month, day;
    char tunits[MAXSTRLEN];
    size_t tunits_len;
    double tunits_multiple, tunits_offset;
    char* basename;
    int i;

    basename = strrchr(fname, '/');
    if (basename == NULL)
        basename = fname;
    else
        basename += 1;

    ncw_open(fname, NC_NOWRITE, &ncid);

    ncw_inq_dimid(fname, ncid, "nobs", &dimid_nobs);
    ncw_inq_dimlen(fname, ncid, dimid_nobs, &nobs_local);
    enkf_printf("        nobs = %u\n", (unsigned int) nobs_local);

    if (nobs_local == 0) {
        ncw_close(fname, ncid);
        return;
    }

    ncw_inq_varid(fname, ncid, "lon", &varid_lon);
    lon = malloc(nobs_local * sizeof(double));
    ncw_get_var_double(fname, ncid, varid_lon, lon);

    ncw_inq_varid(fname, ncid, "lat", &varid_lat);
    lat = malloc(nobs_local * sizeof(double));
    ncw_get_var_double(fname, ncid, varid_lat, lat);

    ncw_inq_varid(fname, ncid, "sst", &varid_sst);
    sst = malloc(nobs_local * sizeof(double));
    ncw_get_var_double(fname, ncid, varid_sst, sst);

    ncw_inq_varid(fname, ncid, "error", &varid_error);
    error_std = malloc(nobs_local * sizeof(double));
    ncw_get_var_double(fname, ncid, varid_error, error_std);

    ncw_inq_varid(fname, ncid, "age", &varid_time);
    time = malloc(nobs_local * sizeof(double));
    ncw_get_var_double(fname, ncid, varid_time, time);
    ncw_inq_attlen(fname, ncid, varid_time, "units", &tunits_len);
    ncw_get_att_text(fname, ncid, varid_time, "units", tunits);
    basename[16] = 0;
    if (!str2int(&basename[14], &day))
        enkf_quit("WindSat reader: could not convert file name \"%s\" to date", fname);
    basename[14] = 0;
    if (!str2int(&basename[12], &month))
        enkf_quit("WindSat reader: could not convert file name \"%s\" to date", fname);
    basename[12] = 0;
    if (!str2int(&basename[8], &year))
        enkf_quit("WindSat reader: could not convert file name \"%s\" to date", fname);
    sprintf(&tunits[tunits_len], " since %4d-%02d-%02d", year, month, day);

    ncw_close(fname, ncid);

    tunits_convert(tunits, &tunits_multiple, &tunits_offset);

    for (i = 0; i < (int) nobs_local; ++i)
        if (time[i] != 0.0)
            break;
    if (i == (int) nobs_local)
        for (i = 0; i < (int) nobs_local; ++i)
            time[i] = 0.5;

    for (i = 0; i < (int) nobs_local; ++i) {
        observation* o;
        obstype* ot;

        obs_checkalloc(obs);
        o = &obs->data[obs->nobs];

        o->product = st_findindexbystring(obs->products, meta->product);
        assert(o->product >= 0);
        o->type = obstype_getid(obs->nobstypes, obs->obstypes, meta->type);
        assert(o->type >= 0);
        ot = &obs->obstypes[o->type];
        o->instrument = st_add_ifabscent(obs->instruments, "WindSat", -1);
        o->id = obs->nobs;
        o->fid = fid;
        o->batch = 0;
        o->value = sst[i];
        o->std = error_std[i];
        o->lon = lon[i];
        o->lat = lat[i];
        o->depth = 0.0;
        o->status = model_xy2fij(m, o->lon, o->lat, &o->fi, &o->fj);
        if (!obs->allobs && o->status == STATUS_OUTSIDEGRID)
            continue;
        if ((o->status == STATUS_OK) && (o->lon <= ot->xmin || o->lon >= ot->xmax || o->lat <= ot->ymin || o->lat >= ot->ymax || o->depth <= ot->zmin || o->depth >= ot->zmax))
            o->status = STATUS_OUTSIDEOBSDOMAIN;
        o->fk = 0.0;
        o->date = time[i] * tunits_multiple + tunits_offset;
        o->aux = -1;

        obs->nobs++;
    }

    free(lon);
    free(lat);
    free(sst);
    free(error_std);
    free(time);
}
