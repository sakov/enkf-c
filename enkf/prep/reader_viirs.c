/******************************************************************************
 *
 * File:        reader_viirs.c        
 *
 * Created:     08/03/2017
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description: Reader for pre-preprocessed (L3C) SST from VIIRS. Parameter
 *              KIND allows one to assimilate NIGHTTIME, WINDY or ALL obs.
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
#include "ncutils.h"
#include "definitions.h"
#include "utils.h"
#include "obsprm.h"
#include "model.h"
#include "grid.h"
#include "observations.h"
#include "prep_utils.h"
#include "allreaders.h"

#define KIND_NIGHTTIME (1 << 0)
#define KIND_WINDY     (1 << 1)
#define KIND_ALL       (KIND_NIGHTTIME | KIND_WINDY)
#define KIND_DEF       KIND_ALL

/**
 */
void reader_viirs_standard(char* fname, int fid, obsmeta* meta, grid* g, observations* obs)
{
    int ksurf = grid_getsurflayerid(g);
    int ncid;
    int ndim;
    int dimid_ni, dimid_nj;
    size_t ni, nj, n;
    int varid_sst, varid_lon, varid_lat, varid_npoints, varid_std, varid_estd, varid_kind, varid_time;
    float* lon;
    float* lat;
    float* sst;
    float* std;
    float* estd;
    short* npoints;
    double* time;
    unsigned char* kind = NULL;
    unsigned char kind_mask = KIND_DEF;

    char tunits[MAXSTRLEN];
    double tunits_multiple, tunits_offset;
    int i, nobs;

#if defined(DEBUG)
    int varid_id;
    int* id;
#endif

    for (i = 0; i < meta->npars; ++i) {
        if (strcasecmp(meta->pars[i].name, "KIND") == 0) {
            int kind_value;

            if (strcasecmp(meta->pars[i].value, "nighttime"))
                kind_value = KIND_NIGHTTIME;
            else if (strcasecmp(meta->pars[i].value, "windy"))
                kind_value = KIND_WINDY;
            else if (!str2int(meta->pars[i].value, &kind_value))
                enkf_quit("%s: can not convert KIND = \"%s\" to int\n", meta->prmfname, meta->pars[i].value);
            if (kind_value < 0 || kind_value > KIND_ALL)
                enkf_printf("KIND: value = %d is outside allowed range [0,%d]\n", kind_value, KIND_ALL);
            kind_mask = (unsigned char) kind_value;
            enkf_printf("        KIND = %d\n", kind_value);
        } else
            enkf_quit("unknown PARAMETER \"%s\"\n", meta->pars[i].name);
    }

    ncw_open(fname, NC_NOWRITE, &ncid);
    ncw_inq_varid(ncid, "sst", &varid_sst);
    ncw_inq_varndims(ncid, varid_sst, &ndim);

    if (ndim != 1 && ndim != 2)
        enkf_quit("%s: # dimensions = %d (1 or 2 expected)", fname, ndim);
    if (ndim == 2) {
        ncw_inq_dimid(ncid, "lon", &dimid_ni);
        ncw_inq_dimlen(ncid, dimid_ni, &ni);
        ncw_inq_dimid(ncid, "lat", &dimid_nj);
        ncw_inq_dimlen(ncid, dimid_nj, &nj);
        enkf_printf("        (ni, nj) = (%u, %u)\n", ni, nj);
        lon = malloc(ni * sizeof(float));
        lat = malloc(nj * sizeof(float));

        n = ni * nj;
    } else {
        int dimid_nobs;

        ncw_inq_dimid(ncid, "nobs", &dimid_nobs);
        ncw_inq_dimlen(ncid, dimid_nobs, &n);
        lon = malloc(n * sizeof(float));
        lat = malloc(n * sizeof(float));

#if defined(DEBUG)
        id = malloc(n * sizeof(int));
        ncw_inq_varid(ncid, "id", &varid_id);
        ncw_get_var_int(ncid, varid_id, id);
#endif
    }

    /*
     * lon/lat
     */
    ncw_inq_varid(ncid, "lon", &varid_lon);
    ncw_inq_varid(ncid, "lat", &varid_lat);
    ncu_readvarfloat(ncid, varid_lon, (ndim == 1) ? n : ni, lon);
    ncu_readvarfloat(ncid, varid_lat, (ndim == 1) ? n : nj, lat);

    /*
     * sst
     */
    sst = malloc(n * sizeof(float));
    ncu_readvarfloat(ncid, varid_sst, n, sst);

    /*
     * std
     */
    ncw_inq_varid(ncid, "std", &varid_std);
    std = malloc(n * sizeof(float));
    ncu_readvarfloat(ncid, varid_std, n, std);

    /*
     * error_std
     */
    ncw_inq_varid(ncid, "error_std", &varid_estd);
    estd = malloc(n * sizeof(float));
    ncu_readvarfloat(ncid, varid_estd, n, estd);

    /*
     * npoints
     */
    ncw_inq_varid(ncid, "npoints", &varid_npoints);
    npoints = malloc(n * sizeof(short));
    ncw_get_var_short(ncid, varid_npoints, npoints);

    if (kind_mask != KIND_ALL) {
        kind = malloc(n);
        ncw_inq_varid(ncid, "kind", &varid_kind);
        ncw_get_var_uchar(ncid, varid_kind, kind);
    }

    ncw_inq_varid(ncid, "time", &varid_time);
    time = malloc(n * sizeof(double));
    ncu_readvardouble(ncid, varid_time, n, time);
    ncw_get_att_text(ncid, varid_time, "units", tunits);

    ncw_close(ncid);

    tunits_convert(tunits, &tunits_multiple, &tunits_offset);

    nobs = 0;
    for (i = 0; i < (int) n; ++i) {
        observation* o;

        if (npoints[i] == 0 || isnan(sst[i]) || isnan(std[i]) || isnan(estd[i]) || isnan(time[i]) || (kind != NULL && kind[i] != kind_mask))
            continue;

        nobs++;
        obs_checkalloc(obs);
        o = &obs->data[obs->nobs];

        o->product = st_findindexbystring(obs->products, meta->product);
        assert(o->product >= 0);
        o->type = obstype_getid(obs->nobstypes, obs->obstypes, meta->type, 1);
        o->instrument = st_add_ifabsent(obs->instruments, "VIIRS", -1);
        o->id = obs->nobs;
        o->fid = fid;
        o->batch = 0;
        o->value = sst[i] - 273.15;
        o->estd = (std[i] > estd[i]) ? std[i] : estd[i];
        if (ndim == 2) {
            o->lon = lon[i % ni];
            o->lat = lat[i / ni];
        } else {
            o->lon = lon[i];
            o->lat = lat[i];
        }
        o->depth = 0.0;
        o->fk = (double) ksurf;
        o->status = grid_xy2fij(g, o->lon, o->lat, &o->fi, &o->fj);
        if (!obs->allobs && o->status == STATUS_OUTSIDEGRID)
            continue;
        o->model_depth = NAN;   /* set in obs_add() */
        o->time = time[i] * tunits_multiple + tunits_offset;
        o->aux = -1;

        obs->nobs++;
    }
    enkf_printf("        nobs = %d\n", nobs);

#if defined(DEBUG)
    if (ndim == 1)
        free(id);
#endif
    free(lon);
    free(lat);
    if (kind != NULL)
        free(kind);
    free(sst);
    free(std);
    free(estd);
    free(npoints);
    free(time);
}
