/******************************************************************************
 *
 * File:        reader_xyh_gridded.c        
 *
 * Created:     02/01/2018
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description: Generic reader for 3D gridded observations with hybrid vertical
 *              coordinate. Adopted from reader_xyz_gridded.c. Assumes
 *              
 *                p[k] = a[k] + b[k] * (p1 - p2),
 *
 *              where p(k) -- z coordinate (pressure) at the middle of layer k.
 *              The grid must be defined in the grid parameter file.
 *
 *                It is currently assumed that there is only one data record
 *              (3D field).
 *                Similarly to reader_xy_gridded(), there are a number of
 *              parameters that must (++) or can be specified if they differ
 *              from the default value (+). Some parameters are optional (-):
 *              - VARNAME (++)
 *              - GRIDNAME (++)
 *              - TIMENAME ("time") (+)
 *              - NPOINTSNAME ("npoints") (-)
 *                  number of collated points for each datum; used basically as
 *                  a data mask n = 0
 *              - STDNAME ("std") (-)
 *                  internal variability of the collated data
 *              - ESTDNAME ("error_std") (-)
 *                  error STD; if absent then needs to be specified externally
 *                  in the oobservation data parameter file
 *              - VARSHIFT (-)
 *                  data offset to be added
 *              - MINDEPTH (-)
 *                  minimal allowed depth
 *              - INSTRUMENT (-)
 *                  instrument string that will be used for calculating
 *                  instrument stats
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
#include "obsmeta.h"
#include "model.h"
#include "grid.h"
#include "observations.h"
#include "prep_utils.h"
#include "allreaders.h"

#define TYPE_DOUBLE 0
#define TYPE_SHORT 1

/**
 */
void reader_xyh_gridded(char* fname, int fid, obsmeta* meta, grid* gdst, observations* obs)
{
    char* varname = NULL;
    char* gridname = NULL;
    grid* gsrc = NULL;

    char* npointsname = NULL;
    char* stdname = NULL;
    char* estdname = NULL;
    char* timename = NULL;
    int ncid;
    int ndim;

    float varshift = 0.0;
    double mindepth = 0.0;
    char instrument[MAXSTRLEN];

    int ni = 0, nj = 0, nk = 0, nij = 0, nijk = 0;

    int varid_var = -1, varid_npoints = -1, varid_std = -1, varid_estd = -1, varid_time = -1;
    float* var = NULL;
    float var_fill_value = NAN;
    float var_add_offset = NAN, var_scale_factor = NAN;
    double var_estd = NAN;
    short* npoints = NULL;
    float* std = NULL;
    float std_add_offset = NAN, std_scale_factor = NAN;
    float std_fill_value = NAN;
    float* estd = NULL;
    float estd_add_offset = NAN, estd_scale_factor = NAN;
    float estd_fill_value = NAN;
    int have_time = 1;
    int singletime = -1;
    float* time = NULL;
    float time_add_offset = NAN, time_scale_factor = NAN;
    float time_fill_value = NAN;
    char tunits[MAXSTRLEN];
    double tunits_multiple = NAN, tunits_offset = NAN;
    int i, j, k, nobs_read;

    strcpy(instrument, meta->product);
    for (i = 0; i < meta->npars; ++i) {
        if (strcasecmp(meta->pars[i].name, "VARNAME") == 0)
            varname = meta->pars[i].value;
        else if (strcasecmp(meta->pars[i].name, "GRIDNAME") == 0)
            gridname = meta->pars[i].value;
        else if (strcasecmp(meta->pars[i].name, "TIMENAME") == 0)
            timename = meta->pars[i].value;
        else if (strcasecmp(meta->pars[i].name, "NPOINTSNAME") == 0)
            npointsname = meta->pars[i].value;
        else if (strcasecmp(meta->pars[i].name, "STDNAME") == 0)
            stdname = meta->pars[i].value;
        else if (strcasecmp(meta->pars[i].name, "ESTDNAME") == 0)
            estdname = meta->pars[i].value;
        else if (strcasecmp(meta->pars[i].name, "VARSHIFT") == 0) {
            if (!str2float(meta->pars[i].value, &varshift))
                enkf_quit("observation prm file: can not convert VARSHIFT = \"%s\" to float\n", meta->pars[i].value);
            enkf_printf("        VARSHIFT = %s\n", meta->pars[i].value);
        } else if (strcasecmp(meta->pars[i].name, "MINDEPTH") == 0) {
            if (!str2double(meta->pars[i].value, &mindepth))
                enkf_quit("observation prm file: can not convert MINDEPTH = \"%s\" to double\n", meta->pars[i].value);
            enkf_printf("        MINDEPTH = %f\n", mindepth);
        } else if (strcasecmp(meta->pars[i].name, "INSTRUMENT") == 0) {
            strncpy(instrument, meta->pars[i].value, MAXSTRLEN);
        } else
            enkf_quit("unknown PARAMETER \"%s\"\n", meta->pars[i].name);
    }
    if (varname == NULL)
        enkf_quit("reader_xyh_gridded(): %s: VARNAME not specified", fname);
    if (gridname == NULL)
        enkf_quit("reader_xyh_gridded(): %s: GRIDNAME not specified", fname);

    gsrc = model_getgridbyname(obs->model, gridname);

    grid_getdims(gsrc, &ni, &nj, &nk);
    nij = ni * nj;
    nijk = nij * nk;

    ncw_open(fname, NC_NOWRITE, &ncid);
    ncw_inq_varid(ncid, varname, &varid_var);
    ncw_inq_varndims(ncid, varid_var, &ndim);
    {
        int dimid[4];
        size_t dimlen[4];

        ncw_inq_vardimid(ncid, varid_var, dimid);
        for (i = 0; i < ndim; ++i)
            ncw_inq_dimlen(ncid, dimid[i], &dimlen[i]);
        if (ndim == 4) {
            if (dimlen[0] != 1)
                enkf_quit("reader_xyh_gridded(): %d records (currently only one record is allowed)", dimlen[0]);
        } else if (ndim != 3)
            enkf_quit("reader_xyh_gridded(): %s: # dimensions = %d (must be 3 or 4 with a single record)", fname, ndim);

        if (dimlen[ndim - 1] != ni || dimlen[ndim - 2] != nj || dimlen[ndim - 3] != nk)
            enkf_quit("dimension mismatch between grid \"%s\" (%d x %d x %d) and variable \"%s\" in \"%s\" (%d x %d x %d)", gridname, dimlen[ndim - 1], dimlen[ndim - 2], dimlen[ndim - 3], varname, fname, ni, nj, nk);
    }

    var = malloc(nijk * sizeof(float));
    ncw_get_var_float(ncid, varid_var, var);
    if (ncw_att_exists(ncid, varid_var, "add_offset")) {
        ncw_get_att_float(ncid, varid_var, "add_offset", &var_add_offset);
        ncw_get_att_float(ncid, varid_var, "scale_factor", &var_scale_factor);
    }
    if (ncw_att_exists(ncid, varid_var, "_FillValue"))
        ncw_get_att_float(ncid, varid_var, "_FillValue", &var_fill_value);

    if (npointsname != NULL)
        ncw_inq_varid(ncid, npointsname, &varid_npoints);
    else if (ncw_var_exists(ncid, "npoints"))
        ncw_inq_varid(ncid, "npoints", &varid_npoints);
    if (varid_npoints >= 0) {
        npoints = malloc(nijk * sizeof(short));
        ncw_get_var_short(ncid, varid_npoints, npoints);
    }

    if (stdname != NULL)
        ncw_inq_varid(ncid, stdname, &varid_std);
    else if (ncw_var_exists(ncid, "std"))
        ncw_inq_varid(ncid, "std", &varid_std);
    if (varid_std >= 0) {
        std = malloc(nijk * sizeof(float));
        ncw_get_var_float(ncid, varid_std, std);
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
    if (varid_estd >= 0) {
        estd = malloc(nijk * sizeof(float));
        ncw_get_var_float(ncid, varid_estd, estd);
        if (ncw_att_exists(ncid, varid_estd, "_FillValue"))
            ncw_get_att_float(ncid, varid_estd, "_FillValue", &estd_fill_value);
        if (ncw_att_exists(ncid, varid_estd, "add_offset")) {
            ncw_get_att_float(ncid, varid_estd, "add_offset", &estd_add_offset);
            ncw_get_att_float(ncid, varid_estd, "scale_factor", &estd_scale_factor);
        }
    }

    if (std == NULL && estd == NULL)
        if (ncw_att_exists(ncid, varid_var, "error_std")) {
            ncw_check_attlen(ncid, varid_var, "error_std", 1);
            ncw_get_att_double(ncid, varid_var, "error_std", &var_estd);
        }

    if (timename != NULL)
        ncw_inq_varid(ncid, timename, &varid_time);
    else if (ncw_var_exists(ncid, "time"))
        ncw_inq_varid(ncid, "time", &varid_time);
    else {
        enkf_printf("        reader_xyh_gridded(): %s: no TIME variable\n", fname);
        have_time = 0;
    }

    if (have_time) {
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
            assert(timelen == nijk);
            time = malloc(nijk * sizeof(float));
        }

        ncw_get_var_float(ncid, varid_time, time);
        if (ncw_att_exists(ncid, varid_time, "_FillValue"))
            ncw_get_att_float(ncid, varid_time, "_FillValue", &time_fill_value);
        if (ncw_att_exists(ncid, varid_time, "add_offset")) {
            ncw_get_att_float(ncid, varid_time, "add_offset", &time_add_offset);
            ncw_get_att_float(ncid, varid_time, "scale_factor", &time_scale_factor);
        }
        ncw_get_att_text(ncid, varid_time, "units", tunits);
        tunits_convert(tunits, &tunits_multiple, &tunits_offset);
    }

    ncw_close(ncid);

    nobs_read = 0;
    for (i = 0; i < ni; ++i) {
        for (j = 0; j < nj; ++j) {
            for (k = 0; k < nk; ++k) {
                int ii = k * nij + j * ni + i;
                observation* o;
                obstype* ot;

                if ((npoints != NULL && npoints[ii] == 0) || var[ii] == var_fill_value || (std != NULL && (std[ii] == std_fill_value || isnan(std[ii]))) || (estd != NULL && (estd[ii] == estd_fill_value || isnan(estd[ii]))) || (have_time && !singletime && (time[ii] == time_fill_value || isnan(time[ii]))))
                    continue;

                nobs_read++;
                obs_checkalloc(obs);
                o = &obs->data[obs->nobs];

                o->product = st_findindexbystring(obs->products, meta->product);
                assert(o->product >= 0);
                o->type = obstype_getid(obs->nobstypes, obs->obstypes, meta->type, 1);
                ot = &obs->obstypes[o->type];
                o->instrument = st_add_ifabsent(obs->instruments, instrument, -1);
                o->id = obs->nobs;
                o->fid = fid;
                o->batch = 0;
                if (!isnan(var_add_offset))
                    o->value = (double) (var[ii] * var_scale_factor + var_add_offset + varshift);
                else
                    o->value = (double) (var[ii] + varshift);
                if (estd == NULL)
                    o->std = var_estd;
                else {
                    if (std == NULL)
                        o->std = 0.0;
                    else {
                        if (!isnan(std_add_offset))
                            o->std = (double) (std[ii] * std_scale_factor + std_add_offset);
                        else
                            o->std = (double) std[ii];
                    }
                    if (!isnan(estd_add_offset)) {
                        double std2 = (double) (estd[ii] * estd_scale_factor + estd_add_offset);

                        o->std = (o->std > std2) ? o->std : std2;
                    } else
                        o->std = (o->std > estd[ii]) ? o->std : estd[ii];
                }
                grid_ij2xy(gsrc, i, j, &o->lon, &o->lat);
                assert(isfinite(o->lon + o->lat));
                o->status = grid_xy2fij(gdst, o->lon, o->lat, &o->fi, &o->fj);
                if (!obs->allobs && o->status == STATUS_OUTSIDEGRID)
                    continue;
                if ((o->status == STATUS_OK) && (o->lon <= ot->xmin || o->lon >= ot->xmax || o->lat <= ot->ymin || o->lat >= ot->ymax))
                    o->status = STATUS_OUTSIDEOBSDOMAIN;
                o->status = grid_fk2z(gsrc, i, j, (double) k, &o->depth);
                if (o->status == STATUS_OK)
                    o->status = grid_z2fk(gdst, o->fi, o->fj, o->depth, &o->fk);
                else
                    o->fk = NAN;
                o->model_depth = NAN;   /* set in obs_add() */
                if (have_time) {
                    float t = (singletime) ? time[0] : time[ii];

                    if (!isnan(time_add_offset))
                        o->date = (double) (t * time_scale_factor + time_add_offset) * tunits_multiple + tunits_offset;
                    else
                        o->date = (double) t* tunits_multiple + tunits_offset;
                } else
                    o->date = NAN;

                o->aux = -1;

                obs->nobs++;
            }
        }
    }
    enkf_printf("        nobs = %d\n", nobs_read);

    free(var);
    if (std != NULL)
        free(std);
    if (estd != NULL)
        free(estd);
    if (npoints != NULL)
        free(npoints);
    if (time != NULL)
        free(time);
}
