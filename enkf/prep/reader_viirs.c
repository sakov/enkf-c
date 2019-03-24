/******************************************************************************
 *
 * File:        reader_viirs.c        
 *
 * Created:     08/03/2017
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description: Reader for pre-preprocessed (L3C) SST from VIIRS.
 *                There are a number of parameters that can be specified:
 *              - VARSHIFT (-)
 *                  data offset to be added. Note that by default the data is
 *                shifted by -273.15. Because VIIRS data by NOAA represent
 *                skin temperature, we suggest set PARAMETER VARSHIFT = 0.16
 *              - MINDEPTH (-)
 *                  minimal allowed depth
 *              
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

#define KIND_NIGHTTIME (1 << 0)
#define KIND_WINDY     (1 << 1)
#define KIND_ALL       (KIND_NIGHTTIME | KIND_WINDY)

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
    double* lon;
    double lon_add_offset, lon_scale_factor;
    double* lat;
    double lat_add_offset, lat_scale_factor;
    short* sst;
    short sst_fill_value;
    double sst_add_offset, sst_scale_factor;
    short* std;
    double std_add_offset, std_scale_factor;
    short std_fill_value;
    short* estd;
    double estd_add_offset, estd_scale_factor;
    short estd_fill_value;
    short* npoints;
    short* time;
    double time_add_offset, time_scale_factor;
    short time_fill_value;
    unsigned char* kind = NULL;
    unsigned char kind_mask = KIND_ALL;

    char tunits[MAXSTRLEN];
    double tunits_multiple, tunits_offset;
    double varshift = 0.0;
    int i, nobs;

#if defined(DEBUG)
    int varid_id;
    int* id;
#endif

    for (i = 0; i < meta->npars; ++i) {
        if (strcasecmp(meta->pars[i].name, "VARSHIFT") == 0) {
            if (!str2double(meta->pars[i].value, &varshift))
                enkf_quit("%s: can not convert VARSHIFT = \"%s\" to double\n", meta->prmfname, meta->pars[i].value);
            enkf_printf("        VARSHIFT = %s\n", meta->pars[i].value);
        }
        /*
         * (MINDEPTH is handled in obs_add() )
         */
        else if (strcasecmp(meta->pars[i].name, "MINDEPTH") == 0) {
            double mindepth;

            if (!str2double(meta->pars[i].value, &mindepth))
                enkf_quit("observation prm file: can not convert MINDEPTH = \"%s\" to double\n", meta->pars[i].value);
            enkf_printf("        MINDEPTH = %.0f\n", mindepth);
            continue;
        } else if (strcasecmp(meta->pars[i].name, "KIND") == 0) {
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
        lon = malloc(ni * sizeof(double));
        lat = malloc(nj * sizeof(double));

        n = ni * nj;
    } else {
        int dimid_nobs;

        ncw_inq_dimid(ncid, "nobs", &dimid_nobs);
        ncw_inq_dimlen(ncid, dimid_nobs, &n);
        lon = malloc(n * sizeof(double));
        lat = malloc(n * sizeof(double));

#if defined(DEBUG)
        id = malloc(n * sizeof(int));
        ncw_inq_varid(ncid, "id", &varid_id);
        ncw_get_var_int(ncid, varid_id, id);
#endif
    }

    ncw_inq_varid(ncid, "lon", &varid_lon);
    ncw_inq_varid(ncid, "lat", &varid_lat);
    ncw_get_var_double(ncid, varid_lon, lon);
    ncw_get_var_double(ncid, varid_lat, lat);
    if (ndim == 1) {
        ncw_get_att_double(ncid, varid_lon, "add_offset", &lon_add_offset);
        ncw_get_att_double(ncid, varid_lon, "scale_factor", &lon_scale_factor);
        ncw_get_att_double(ncid, varid_lat, "add_offset", &lat_add_offset);
        ncw_get_att_double(ncid, varid_lat, "scale_factor", &lat_scale_factor);
    }

    sst = malloc(n * sizeof(short));
    ncw_get_var_short(ncid, varid_sst, sst);
    ncw_get_att_double(ncid, varid_sst, "add_offset", &sst_add_offset);
    ncw_get_att_double(ncid, varid_sst, "scale_factor", &sst_scale_factor);
    ncw_get_att_short(ncid, varid_sst, "_FillValue", &sst_fill_value);

    ncw_inq_varid(ncid, "std", &varid_std);
    std = malloc(n * sizeof(short));
    ncw_get_var_short(ncid, varid_std, std);
    ncw_get_att_double(ncid, varid_std, "add_offset", &std_add_offset);
    ncw_get_att_double(ncid, varid_std, "scale_factor", &std_scale_factor);
    ncw_get_att_short(ncid, varid_std, "_FillValue", &std_fill_value);

    ncw_inq_varid(ncid, "error_std", &varid_estd);
    estd = malloc(n * sizeof(short));
    ncw_get_var_short(ncid, varid_estd, estd);
    ncw_get_att_double(ncid, varid_estd, "add_offset", &estd_add_offset);
    ncw_get_att_double(ncid, varid_estd, "scale_factor", &estd_scale_factor);
    ncw_get_att_short(ncid, varid_estd, "_FillValue", &estd_fill_value);

    ncw_inq_varid(ncid, "npoints", &varid_npoints);
    npoints = malloc(n * sizeof(short));
    ncw_get_var_short(ncid, varid_npoints, npoints);

    if (kind_mask != KIND_ALL) {
        kind = malloc(n);
        ncw_inq_varid(ncid, "kind", &varid_kind);
        ncw_get_var_uchar(ncid, varid_kind, kind);
    }

    ncw_inq_varid(ncid, "time", &varid_time);
    time = malloc(n * sizeof(short));
    ncw_get_var_short(ncid, varid_time, time);
    ncw_get_att_double(ncid, varid_time, "add_offset", &time_add_offset);
    ncw_get_att_double(ncid, varid_time, "scale_factor", &time_scale_factor);
    ncw_get_att_short(ncid, varid_time, "_FillValue", &time_fill_value);
    ncw_get_att_text(ncid, varid_time, "units", tunits);

    ncw_close(ncid);

    tunits_convert(tunits, &tunits_multiple, &tunits_offset);

    nobs = 0;
    for (i = 0; i < (int) n; ++i) {
        observation* o;

        if (npoints[i] == 0 || sst[i] == sst_fill_value || std[i] == std_fill_value || estd[i] == estd_fill_value || time[i] == time_fill_value || (kind != NULL && kind[i] != kind_mask))
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
        o->value = (double) sst[i] * sst_scale_factor + sst_add_offset + varshift - 273.15;
        {
            double std1 = (double) std[i] * std_scale_factor + std_add_offset;
            double std2 = (double) estd[i] * estd_scale_factor + estd_add_offset;

            o->estd = (std1 > std2) ? std1 : std2;
        }
        if (ndim == 2) {
            o->lon = lon[i % ni];
            o->lat = lat[i / ni];
        } else {
            o->lon = (double) lon[i] * lon_scale_factor + lon_add_offset;
            o->lat = (double) lat[i] * lat_scale_factor + lat_add_offset;
        }
        o->depth = 0.0;
        o->fk = (double) ksurf;
        o->status = grid_xy2fij(g, o->lon, o->lat, &o->fi, &o->fj);
        if (!obs->allobs && o->status == STATUS_OUTSIDEGRID)
            continue;
        o->model_depth = NAN;   /* set in obs_add() */
        o->day = ((double) time[i] * time_scale_factor + time_add_offset) * tunits_multiple + tunits_offset;
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
