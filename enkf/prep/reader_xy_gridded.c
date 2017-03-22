/******************************************************************************
 *
 * File:        reader_xy_gridded.c        
 *
 * Created:     08/03/2017
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description: Generic reader for gridded surface observations.
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
#include "definitions.h"
#include "utils.h"
#include "obsmeta.h"
#include "model.h"
#include "grid.h"
#include "observations.h"
#include "prep_utils.h"

#define TYPE_DOUBLE 0
#define TYPE_SHORT 1

void reader_xy_gridded(char* fname, int fid, obsmeta* meta, model* m, observations* obs)
{
    char* varname = NULL;
    char* lonname = NULL;
    char* latname = NULL;
    char* npointsname = NULL;
    char* stdname = NULL;
    char* estdname = NULL;
    char* timename = NULL;
    int ncid;
    int ndim;

    int dimid_i, dimid_j;
    size_t ni, nj, n;
    int varid_lon, varid_lat;
    double* lon = NULL;
    double* lat = NULL;

    int varid_var = -1, varid_npoints = -1, varid_std = -1, varid_estd = -1, varid_time = -1;
    float* var = NULL;
    float var_shift = 0.0;
    float var_fill_value = NAN;
    float var_add_offset = NAN, var_scale_factor = NAN;
    short* npoints = NULL;
    float* std = NULL;
    float std_add_offset = NAN, std_scale_factor = NAN;
    float std_fill_value = NAN;
    float* estd = NULL;
    float estd_add_offset = NAN, estd_scale_factor = NAN;
    float estd_fill_value = NAN;
    int singletime;
    float* time = NULL;
    float time_add_offset = NAN, time_scale_factor = NAN;
    float time_fill_value = NAN;
    char tunits[MAXSTRLEN];
    double tunits_multiple = NaN, tunits_offset = NaN;
    int mvid;
    float** depth;
    int ktop;
    int i, nobs_read;

    for (i = 0; i < meta->npars; ++i) {
        if (strcasecmp(meta->pars[i].name, "VARNAME") == 0)
            varname = meta->pars[i].value;
        else if (strcasecmp(meta->pars[i].name, "NPOINTSNAME") == 0)
            npointsname = meta->pars[i].value;
        else if (strcasecmp(meta->pars[i].name, "LONNAME") == 0)
            lonname = meta->pars[i].value;
        else if (strcasecmp(meta->pars[i].name, "LATNAME") == 0)
            latname = meta->pars[i].value;
        else if (strcasecmp(meta->pars[i].name, "STDNAME") == 0)
            stdname = meta->pars[i].value;
        else if (strcasecmp(meta->pars[i].name, "ESTDNAME") == 0)
            estdname = meta->pars[i].value;
        else if (strcasecmp(meta->pars[i].name, "VARSHIFT") == 0) {
            if (!str2float(meta->pars[i].value, &var_shift))
                enkf_quit("observation prm file: can not convert VARSHIFT = \"%s\" to double\n", meta->pars[i].value);
            enkf_printf("        VARSHIFT = %.0f\n", var_shift);
        } else
            enkf_quit("unknown PARAMETER \"%s\"\n", meta->pars[i].name);
    }
    if (varname == NULL)
        enkf_quit("reader_xy_gridded(): %s variable name not specified", fname);

    ncw_open(fname, NC_NOWRITE, &ncid);
    ncw_inq_varid(ncid, varname, &varid_var);
    ncw_inq_varndims(ncid, varid_var, &ndim);
    if (ndim != 2)
        enkf_quit("reader_xy_gridded(): %s: # dimensions = %d (must be 2)", fname, ndim);
    ncw_check_varndims(ncid, lonname, 1);
    ncw_check_varndims(ncid, latname, 1);

    ncw_inq_varid(ncid, lonname, &varid_lon);
    ncw_inq_vardimid(ncid, varid_lon, &dimid_i);
    ncw_inq_dimlen(ncid, dimid_i, &ni);
    ncw_inq_varid(ncid, latname, &varid_lat);
    ncw_inq_vardimid(ncid, varid_lat, &dimid_j);
    ncw_inq_dimlen(ncid, dimid_j, &nj);
    enkf_printf("        (ni, nj) = (%u, %u)\n", ni, nj);
    n = ni * nj;

    lon = malloc(ni * sizeof(double));
    lat = malloc(nj * sizeof(double));
    ncw_get_var_double(ncid, varid_lon, lon);
    ncw_get_var_double(ncid, varid_lat, lat);

    var = malloc(n * sizeof(float));
    ncw_get_var_float(ncid, varid_var, var);
    ncw_get_att_float(ncid, varid_var, "add_offset", &var_add_offset);
    ncw_get_att_float(ncid, varid_var, "scale_factor", &var_scale_factor);
    ncw_get_att_float(ncid, varid_var, "_FillValue", &var_fill_value);

    if (npointsname != NULL)
        ncw_inq_varid(ncid, npointsname, &varid_npoints);
    else if (ncw_var_exists(ncid, "npoints"))
        ncw_inq_varid(ncid, "npoints", &varid_npoints);
    if (varid_npoints >= 0) {
        npoints = malloc(n * sizeof(short));
        ncw_get_var_short(ncid, varid_npoints, npoints);
    }

    if (stdname != NULL)
        ncw_inq_varid(ncid, stdname, &varid_std);
    else if (ncw_var_exists(ncid, "std"))
        ncw_inq_varid(ncid, "std", &varid_std);
    if (varid_std >= 0) {
        std = malloc(n * sizeof(float));
        if (ncw_att_exists(ncid, varid_std, "_FillValue"))
            ncw_get_att_float(ncid, varid_std, "_FillValue", &std_fill_value);
        if (ncw_att_exists(ncid, varid_std, "add_offset")) {
            ncw_get_att_float(ncid, varid_std, "add_offset", &std_add_offset);
            ncw_get_att_float(ncid, varid_std, "scale_factor", &std_scale_factor);
        }
    }

    if (estdname != NULL)
        ncw_inq_varid(ncid, estdname, &varid_estd);
    else if (ncw_var_exists(ncid, "error_std"))
        ncw_inq_varid(ncid, "error_std", &varid_estd);
    ncw_inq_varid(ncid, "error_std", &varid_estd);
    estd = malloc(n * sizeof(float));
    if (varid_estd >= 0) {
        estd = malloc(n * sizeof(float));
        if (ncw_att_exists(ncid, varid_estd, "_FillValue"))
            ncw_get_att_float(ncid, varid_estd, "_FillValue", &estd_fill_value);
        if (ncw_att_exists(ncid, varid_estd, "add_offset")) {
            ncw_get_att_float(ncid, varid_estd, "add_offset", &estd_add_offset);
            ncw_get_att_float(ncid, varid_estd, "scale_factor", &estd_scale_factor);
        }
    }

    if (timename != NULL)
        ncw_inq_varid(ncid, timename, &varid_time);
    else if (ncw_var_exists(ncid, "time"))
        ncw_inq_varid(ncid, "time", &varid_time);
    else
        enkf_quit("reader_xy_scattered(): %s: could not find TIME variable", fname);
    {
        int timendims;
        int timedimids[NC_MAX_DIMS];
        size_t timelen = 1;

        ncw_inq_varndims(ncid, varid_time, &timendims);
        ncw_inq_vardimid(ncid, varid_time, timedimids);
        for (i = 0; i < timendims; ++i) {
            size_t dimlen;

            ncw_inq_dimlen(ncid, timedimids[i], &dimlen);
            timelen *= dimlen;
        }

        if (timelen == 1) {
            singletime = 1;
            time = malloc(sizeof(float));
        } else {
            singletime = 0;
            assert(timelen == n);
            time = malloc(n * sizeof(float));
        }
    }
    ncw_get_var_float(ncid, varid_time, time);
    if (ncw_att_exists(ncid, varid_time, "_FillValue"))
        ncw_get_att_float(ncid, varid_time, "_FillValue", &time_fill_value);
    if (ncw_att_exists(ncid, varid_time, "add_offset")) {
        ncw_get_att_float(ncid, varid_time, "add_offset", &time_add_offset);
        ncw_get_att_float(ncid, varid_time, "scale_factor", &time_scale_factor);
    }

    ncw_close(ncid);

    tunits_convert(tunits, &tunits_multiple, &tunits_offset);

    mvid = model_getvarid(m, obs->obstypes[obstype_getid(obs->nobstypes, obs->obstypes, meta->type, 1)].varnames[0], 1);
    ktop = grid_gettoplayerid(model_getvargrid(m, mvid));
    depth = model_getdepth(m, mvid, 0);

    nobs_read = 0;
    for (i = 0; i < (int) n; ++i) {
        observation* o;
        obstype* ot;

        if ((npoints != NULL && npoints[i] == 0) || var[i] == var_fill_value || (std != NULL && std[i] == std_fill_value) || estd[i] == estd_fill_value || (!singletime && time[i] == time_fill_value))
            continue;

        nobs_read++;
        obs_checkalloc(obs);
        o = &obs->data[obs->nobs];

        o->product = st_findindexbystring(obs->products, meta->product);
        assert(o->product >= 0);
        o->type = obstype_getid(obs->nobstypes, obs->obstypes, meta->type, 1);
        ot = &obs->obstypes[o->type];
        o->instrument = st_add_ifabscent(obs->instruments, "VIIRS", -1);
        o->id = obs->nobs;
        o->fid = fid;
        o->batch = 0;
        if (isnan(var_add_offset))
            o->value = (double) (var[i] + var_shift);
        else
            o->value = (double) (var[i] * var_scale_factor + var_add_offset + var_shift);
        if (std == NULL)
            o->std = 0.0;
        else {
            if (isnan(std_add_offset))
                o->std = (double) (std[i] * std_scale_factor + std_add_offset);
            else
                o->std = (double) std[i];
        }
        if (isnan(estd_add_offset)) {
            double std2 = (double) (estd[i] * estd_scale_factor + estd_add_offset);
            o->std = (o->std > std2) ? o->std : std2;
        } else
            o->std = (o->std > estd[i]) ? o->std : estd[i];
        o->lon = lon[i % ni];
        o->lat = lat[i / ni];
        o->depth = 0.0;
        o->fk = (double) ktop;
        o->status = model_xy2fij(m, mvid, o->lon, o->lat, &o->fi, &o->fj);
        if (!obs->allobs && o->status == STATUS_OUTSIDEGRID)
            continue;
        if ((o->status == STATUS_OK) && (o->lon <= ot->xmin || o->lon >= ot->xmax || o->lat <= ot->ymin || o->lat >= ot->ymax))
            o->status = STATUS_OUTSIDEOBSDOMAIN;
        o->model_depth = (depth == NULL || isnan(o->fi + o->fj)) ? NaN : depth[(int) (o->fj + 0.5)][(int) (o->fi + 0.5)];
        {
            float t = (singletime) ? time[0] : time[i];

            if (!isnan(time_add_offset))
                o->date = (double) (t * time_scale_factor + time_add_offset) * tunits_multiple + tunits_offset;
            else
                o->date = (double) t * tunits_multiple + tunits_offset;
        }
            
        o->aux = -1;

        obs->nobs++;
    }
    enkf_printf("        nobs = %d\n", nobs_read);

    free(lon);
    free(lat);
    free(var);
    if (std != NULL)
        free(std);
    free(estd);
    if (npoints != NULL)
        free(npoints);
    free(time);
}
