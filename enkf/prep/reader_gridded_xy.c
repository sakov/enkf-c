/******************************************************************************
 *
 * File:        reader_gridded_xy.c        
 *
 * Created:     08/03/2017
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description: Generic reader for gridded surface observations.
 *                It currently assumes the following:
 *              - there is only one data record (2D field);
 *              - longitude is the inner ("fast") coordinate of the variable.
 * Revisions:   PS 6/7/2018
 *                Added parameters QCFLAGNAME and QCFLAGVALS. The latter is
 *                supposed to contain a list of allowed flag values.
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

#define TYPE_DOUBLE 0
#define TYPE_SHORT 1

/**
 */
void reader_gridded_xy(char* fname, int fid, obsmeta* meta, grid* g, observations* obs)
{
    int ksurf = grid_getsurflayerid(g);
    char* varname = NULL;
    char* lonname = NULL;
    char* latname = NULL;
    double zvalue = NAN;
    int zvalueentered = 0;
    char* npointsname = NULL;
    char* stdname = NULL;
    char* estdname = NULL;
    char* batchname = NULL;
    char instrument[MAXSTRLEN] = "";
    int ndim_var, ndim_xy;
    size_t dimlen_var[3], dimlen_xy[2];
    int nqcflagvars = 0;
    char** qcflagvarnames = NULL;
    uint32_t* qcflagmasks = NULL;

    int instid = -1;
    int productid = -1;
    int typeid = -1;

    int ncid;
    int iscurv = -1;
    size_t ni = 0, nj = 0, nij = 0;
    double* lon = NULL;
    double* lat = NULL;
    float* var = NULL;
    double var_estd = NAN;
    short* npoints = NULL;
    float* std = NULL;
    float* estd = NULL;
    int* batch = NULL;
    uint32_t** qcflag = NULL;
    size_t ntime = 0;
    double* time = NULL;
    int varid;
    size_t i, nobs_read;

    for (i = 0; i < meta->npars; ++i) {
        if (strcasecmp(meta->pars[i].name, "VARNAME") == 0)
            varname = meta->pars[i].value;
        else if (strcasecmp(meta->pars[i].name, "NPOINTSNAME") == 0)
            npointsname = meta->pars[i].value;
        else if (strcasecmp(meta->pars[i].name, "LONNAME") == 0)
            lonname = meta->pars[i].value;
        else if (strcasecmp(meta->pars[i].name, "LATNAME") == 0)
            latname = meta->pars[i].value;
        else if (strcasecmp(meta->pars[i].name, "ZVALUE") == 0) {
            if (!str2double(meta->pars[i].value, &zvalue))
                enkf_quit("observation prm file: can not convert ZVALUE = \"%s\" to double\n", meta->pars[i].value);
            zvalueentered = 1;
        } else if (strcasecmp(meta->pars[i].name, "STDNAME") == 0)
            stdname = meta->pars[i].value;
        else if (strcasecmp(meta->pars[i].name, "ESTDNAME") == 0)
            estdname = meta->pars[i].value;
        else if (strcasecmp(meta->pars[i].name, "BATCHNAME") == 0)
            batchname = meta->pars[i].value;
        else if (strcasecmp(meta->pars[i].name, "INSTRUMENT") == 0)
            strncpy(instrument, meta->pars[i].value, MAXSTRLEN - 1);
        else if (strcasecmp(meta->pars[i].name, "TIMENAME") == 0 || strcasecmp(meta->pars[i].name, "TIMENAMES") == 0)
            /*
             * TIMENAME and TIMENAMES are dealt with separately
             */
            ;
        else if (strcasecmp(meta->pars[i].name, "QCFLAGNAME") == 0 || strcasecmp(meta->pars[i].name, "QCFLAGVARNAME") == 0 || strcasecmp(meta->pars[i].name, "QCFLAGVALS") == 0)
            /*
             * QCFLAGNAME and QCFLAGVALS are dealt with separately
             */
            ;
        else
            enkf_quit("unknown PARAMETER \"%s\"\n", meta->pars[i].name);
    }

    if (varname == NULL)
        enkf_quit("reader_xy_gridded(): %s: VARNAME not specified", fname);
    else
        enkf_printf("        VARNAME = %s\n", varname);

    ncw_open(fname, NC_NOWRITE, &ncid);

    /*
     * main variable
     */
    ncw_inq_varid(ncid, varname, &varid);
    ncw_inq_vardims(ncid, varid, 3, &ndim_var, dimlen_var);
    if (ndim_var == 3) {
        if (dimlen_var[0] != 1)
            enkf_quit("reader_xy_gridded(): %s: %s: %d records (currently only one is allowed)", fname, varname, dimlen_var[0]);
        nij = dimlen_var[1] * dimlen_var[2];
    } else if (ndim_var == 2) {
        if (ncw_var_hasunlimdim(ncid, varid))
            enkf_quit("reader_xy_gridded(): %s: %s: not enough spatial dimensions (must be 2)", fname, varname);
        nij = dimlen_var[0] * dimlen_var[1];
    } else if (ndim_var != 2)
        enkf_quit("reader_xy_gridded(): %s: %s: %d dimensions (must be either 2 or 3 with only one record)", fname, varname, ndim_var);

    var = malloc(nij * sizeof(float));
    ncu_readvarfloat(ncid, varid, nij, var);

    /*
     * longitude
     */
    lonname = get_lonname(ncid, lonname);
    if (lonname != NULL) {
        enkf_printf("        LONNAME = %s\n", lonname);
        ncw_inq_varid(ncid, lonname, &varid);
    } else
        enkf_quit("reader_xy_gridded(): %s: could not find longitude variable", fname);
    ncw_inq_vardims(ncid, varid, 2, &ndim_xy, dimlen_xy);
    if (ndim_xy == 1) {
        iscurv = 0;
        ni = dimlen_xy[0];
    } else if (ndim_xy == 2) {
        iscurv = 1;
        ni = dimlen_xy[1];
        nj = dimlen_xy[0];
    } else
        enkf_quit("reader_xy_gridded(): %s: coordinate variable \"%s\" has neither 1 or 2 dimensions", fname, lonname);

    if (iscurv == 0) {
        lon = malloc(ni * sizeof(double));
        ncu_readvardouble(ncid, varid, ni, lon);
    } else {
        lon = malloc(nij * sizeof(double));
        ncu_readvardouble(ncid, varid, nij, lon);
    }

    /*
     * latitude
     */
    latname = get_latname(ncid, latname);
    if (latname != NULL) {
        enkf_printf("        LATNAME = %s\n", latname);
        ncw_inq_varid(ncid, latname, &varid);
    } else
        enkf_quit("reader_xyz_gridded(): %s: could not find latitude variable", fname);
    if (iscurv == 0) {
        ncw_check_varndims(ncid, varid, 1);
        ncw_inq_vardims(ncid, varid, 1, NULL, &nj);
    } else
        ncw_check_vardims(ncid, varid, 2, dimlen_xy);

    enkf_printf("        (ni, nj) = (%u, %u)\n", ni, nj);
    if (ni * nj != nij)
        enkf_quit("reader_xy_gridded(): %s: dimensions of variable \"%s\" do not match coordinate dimensions", fname, varname);
    if (dimlen_var[ndim_var - 1] != ni)
        enkf_quit("reader_xy_gridded(): %s: %s: longitude must be the inner coordinate", fname, varname);

    if (iscurv == 0) {
        lat = malloc(nj * sizeof(double));
        ncu_readvardouble(ncid, varid, nj, lat);
    } else {
        lat = malloc(nij * sizeof(double));
        ncu_readvardouble(ncid, varid, nij, lat);
    }

    /*
     * npoints
     */
    varid = -1;
    if (npointsname != NULL)
        ncw_inq_varid(ncid, npointsname, &varid);
    else if (ncw_var_exists(ncid, "npoints"))
        ncw_inq_varid(ncid, "npoints", &varid);
    if (varid >= 0) {
        npoints = malloc(nij * sizeof(short));
        ncw_get_var_short(ncid, varid, npoints);
    }

    /*
     * std
     */
    varid = -1;
    if (stdname != NULL)
        ncw_inq_varid(ncid, stdname, &varid);
    else if (ncw_var_exists(ncid, "std"))
        ncw_inq_varid(ncid, "std", &varid);
    if (varid >= 0) {
        std = malloc(nij * sizeof(float));
        ncu_readvarfloat(ncid, varid, nij, std);
    }

    /*
     * estd
     */
    varid = -1;
    if (estdname != NULL)
        ncw_inq_varid(ncid, estdname, &varid);
    else if (ncw_var_exists(ncid, "error_std"))
        ncw_inq_varid(ncid, "error_std", &varid);
    if (varid >= 0) {
        estd = malloc(nij * sizeof(float));
        ncu_readvarfloat(ncid, varid, nij, estd);
    }

    if (std == NULL && estd == NULL) {
        ncw_inq_varid(ncid, varname, &varid);
        if (ncw_att_exists(ncid, varid, "error_std")) {
            ncw_check_attlen(ncid, varid, "error_std", 1);
            ncw_get_att_double(ncid, varid, "error_std", &var_estd);
        }
    }

    /*
     * batch
     */
    varid = -1;
    if (batchname != NULL)
        ncw_inq_varid(ncid, batchname, &varid);
    else if (ncw_var_exists(ncid, "batch"))
        ncw_inq_varid(ncid, "batch", &varid);
    if (varid >= 0) {
        ncw_check_varsize(ncid, varid, nij);
        batch = malloc(nij * sizeof(int));
        ncw_get_var_int(ncid, varid, batch);
    }

    /*
     * qcflags
     */
    get_qcflags(meta, &nqcflagvars, &qcflagvarnames, &qcflagmasks);
    if (nqcflagvars > 0) {
        qcflag = alloc2d(nqcflagvars, nij, sizeof(int32_t));
        for (i = 0; i < nqcflagvars; ++i) {
            ncw_inq_varid(ncid, qcflagvarnames[i], &varid);
            ncw_check_varsize(ncid, varid, nij);
            ncw_get_var_uint(ncid, varid, qcflag[i]);
        }
    }

    /*
     * time
     */
    get_time(meta, ncid, &ntime, &time);
    assert(ntime == nij || ntime <= 1);

    /*
     * instrument
     */
    if (strlen(instrument) == 0 && !get_insttag(ncid, varname, instrument))
        strncpy(instrument, meta->product, MAXSTRLEN - 1);

    ncw_close(ncid);

    instid = st_add_ifabsent(obs->instruments, instrument, -1);
    productid = st_findindexbystring(obs->products, meta->product);
    assert(productid >= 0);
    typeid = obstype_getid(obs->nobstypes, obs->obstypes, meta->type, 1);
    nobs_read = 0;
    for (i = 0; i < nij; ++i) {
        observation* o;
        int ii;

        if ((npoints != NULL && npoints[i] == 0) || isnan(var[i]) || (std != NULL && isnan(std[i])) || (estd != NULL && isnan(estd[i])) || (ntime == nij && isnan(time[i])))
            continue;
        for (ii = 0; ii < nqcflagvars; ++ii)
            if (!((1 << qcflag[ii][i]) & qcflagmasks[ii]))
                goto nextob;

        nobs_read++;
        obs_checkalloc(obs);
        o = &obs->data[obs->nobs];

        o->product = productid;
        o->type = typeid;
        o->instrument = instid;
        o->id = obs->nobs;
        o->fid = fid;
        o->batch = (batch == NULL) ? 0 : batch[i];
        o->value = (double) var[i];
        if (estd == NULL)
            o->estd = var_estd;
        else {
            if (std == NULL)
                o->estd = estd[i];
            else
                o->estd = (std[i] > estd[i]) ? std[i] : estd[i];
        }
        if (iscurv == 0) {
            o->lon = lon[i % ni];
            o->lat = lat[i / ni];
        } else {
            o->lon = lon[i];
            o->lat = lat[i];
        }
        o->status = grid_xy2fij_f(g, o->lon, o->lat, &o->fi, &o->fj);
        if (!obs->allobs && o->status == STATUS_OUTSIDEGRID)
            continue;
        if (zvalueentered) {
            o->depth = zvalue;
            if (isfinite(zvalue))
                o->status = grid_z2fk_f(g, o->fi, o->fj, o->depth, &o->fk);
            else
                o->fk = NAN;
        } else {
            o->depth = 0.0;
            o->fk = (double) ksurf;
        }
        o->model_depth = NAN;   /* set in obs_add() */
        if (ntime > 0)
            o->time = (ntime == 1) ? time[0] : time[i];
        else
            o->time = NAN;

        o->aux = -1;

        obs->nobs++;
      nextob:
        ;
    }
    enkf_printf("        nobs = %d\n", nobs_read);

    free(lon);
    free(lat);
    free(var);
    if (std != NULL)
        free(std);
    if (estd != NULL)
        free(estd);
    if (batch != NULL)
        free(batch);
    if (npoints != NULL)
        free(npoints);
    if (time != NULL)
        free(time);
    if (nqcflagvars > 0) {
        free(qcflagvarnames);
        free(qcflagmasks);
        free(qcflag);
    }
}

/**
 */
void reader_gridded_xy_describe(void)
{
   enkf_printf("\n  Generic reader \"gridded_xy\" reads 2D gridded data. It can handle either\n\
  curvilinear or geographically rectangular grids. In the latter case the\n\
  coordinates can be represented by 1D variables; longitude is assumed to be\n\
  the inner (\"fast\") coordinate. Currently there can be only 1 time record in\n\
  a data file.\n\
\n\
  There are a number of parameters that must (marked below with \"++\"), can\n\
  (\"+\"), or may (\"-\") be specified in the corresponding section of the\n\
  observation data parameter file. The names in brackets represent the default\n\
  names checked in the abscence of the entry for the parameter. Each parameter\n\
  needs to be entered as follows:\n\
    PARAMETER <name> = <value> ...\n\
\n\
  Parameters common to generic readers:\n\
    - VARNAME (++)\n\
    - TIMENAME (\"*[tT][iI][mM][eE]*\") (+)\n\
    - or TIMENAMES (when time = base_time + offset) (+)\n\
    - LONNAME (\"lon\" | \"longitude\") (+)\n\
    - LATNAME (\"lat\" | \"latitude\") (+)\n\
    - STDNAME (\"std\") (-)\n\
        dispersion of the collated data\n\
    - ESTDNAME (\"error_std\") (-)\n\
        error STD; if absent then needs to be specified in the corresponding\n\
        section of the observation data parameter file\n\
    - BATCHNAME (\"batch\") (-)\n\
        name of the variable used for batch ID (e.g. \"pass\" for SLA)\n\
    - INSTRUMENT (-)\n\
        instrument string that will be used for calculating instrument stats\n\
        (overrides the global attribute \"instrument\" in the data file)\n\
    - QCFLAGNAME (-)\n\
        name of the QC flag variable, possible values 0 <= qcflag <= 31\n\
    - QCFLAGVALS (-)\n\
        the list of allowed values of QC flag variable\n\
        Note: it is possible to have multiple entries of QCFLAGNAME and\n\
        QCFLAGVALS combination, e.g.:\n\
          PARAMETER QCFLAGNAME = TEMP_quality_control\n\
          PARAMETER QCFLAGVALS = 1\n\
          PARAMETER QCFLAGNAME = DEPTH_quality_control\n\
          PARAMETER QCFLAGVALS = 1\n\
          PARAMETER QCFLAGNAME = LONGITUDE_quality_control\n\
          PARAMETER QCFLAGVALS = 1,8\n\
          PARAMETER QCFLAGNAME = LATITUDE_quality_control\n\
          PARAMETER QCFLAGVALS = 1,8\n\
        An observation is considered valid if each of the specified flags takes\n\
        a permitted value.\n\
  Parameters specific to the reader:\n\
    - ZVALUE (0) (+)\n\
        height/depth (can be NaN)\n\
    - NPOINTSNAME (\"npoints\") (-)\n\
        number of collated points for each datum; used basically as a data mask\n\
        when n = 0\n\
  Parameters common to all readers:\n\
    - VARSHIFT (-)\n\
        data offset to be added (e.g. -273.15 to convert from K to C)\n\
    - FOOTRPINT (-)\n\
        footprint of observations in km\n\
    - MINDEPTH (-)\n\
        minimal allowed depth\n\
    - MAXDEPTH (-)\n\
        maximal allowed depth\n\
    - THIN (-)\n\
        data thinning ratio (only one out of each consequitive <THIN> values is\n\
        read\n");
}
