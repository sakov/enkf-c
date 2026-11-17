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
 *              PS 17/11/2026
 *                Factored common attribute handling into obsatts.
 *
 *****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "ncw.h"
#include "ncutils.h"
#include "global.h"
#include "utils.h"
#include "obsprm.h"
#include "model.h"
#include "grid.h"
#include "observations.h"
#include "prep_utils.h"
#include "obsatts.h"
#include "allreaders.h"

/**
 */
void reader_gridded_xy(char* fname, int fid, obssection* section, grid* g, observations* obs)
{
    int ksurf = grid_getsurflayerid(g);
    char* varname = NULL;
    char* lonname = NULL;
    char* latname = NULL;
    double zvalue = NAN;
    int zvalueentered = 0;
    char* npointsname = NULL;
    obsatts oa;
    int ndim_var, ndim_xy;
    size_t dimlen_var[3], dimlen_xy[2];

    int productid = -1;
    int typeid = -1;

    int ncid;
    int iscurv = -1;
    size_t ni = 0, nj = 0, nij = 0;
    double* lon = NULL;
    double* lat = NULL;
    float* var = NULL;
    short* npoints = NULL;
    size_t ntime = 0;
    double* time = NULL;
    int varid;
    size_t i, nobs_read;

    obsatts_init(&oa);

    for (i = 0; i < (size_t) section->npars; ++i) {
        if (strcasecmp(section->pars[i].name, "VARNAME") == 0)
            varname = section->pars[i].value;
        else if (strcasecmp(section->pars[i].name, "NPOINTSNAME") == 0)
            npointsname = section->pars[i].value;
        else if (strcasecmp(section->pars[i].name, "LONNAME") == 0)
            lonname = section->pars[i].value;
        else if (strcasecmp(section->pars[i].name, "LATNAME") == 0)
            latname = section->pars[i].value;
        else if (strcasecmp(section->pars[i].name, "ZVALUE") == 0) {
            if (!str2double(section->pars[i].value, &zvalue))
                enkf_quit("observation prm file: can not convert ZVALUE = \"%s\" to double\n", section->pars[i].value);
            zvalueentered = 1;
        } else if (!obsatts_parse_param(&oa, section->pars[i].name, section->pars[i].value))
            enkf_quit("unknown PARAMETER \"%s\"\n", section->pars[i].name);
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
            enkf_quit("reader_xy_gridded(): %s: %s: %zu records (currently only one is allowed)", fname, varname, dimlen_var[0]);
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

    enkf_printf("        (ni, nj) = (%zu, %zu)\n", ni, nj);
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
        ncw_check_varsize(ncid, varid, nij);
        npoints = malloc(nij * sizeof(short));
        ncw_get_var_short(ncid, varid, npoints);
    }

    /*
     * std, estd, batch, qcflag, instrument
     */
    obsatts_read(&oa, section, ncid, varname, nij, nij, 1, 0, 0, obs);

    /*
     * time
     */
    get_time(section, ncid, &ntime, &time);
    assert(ntime == nij || ntime == nj || ntime <= 1);
    if (ntime == nj) {
        double* tmp = calloc(nij, sizeof(double));
        size_t j, ij;

        for (j = 0, ij = 0; j < nj; ++j)
            for (i = 0; i < ni; ++i, ++ij)
                tmp[ij] = time[j];
        free(time);
        time = tmp;
        ntime = nij;
    }

    ncw_close(ncid);

    productid = st_findindexbystring(obs->products, section->product);
    assert(productid >= 0);
    typeid = obstype_getid(obs->nobstypes, obs->obstypes, section->type, 1);
    nobs_read = 0;
    for (i = 0; i < nij; ++i) {
        observation* o;

        if ((npoints != NULL && npoints[i] == 0) || isnan(var[i]) || (ntime == nij && isnan(time[i])))
            continue;
        if (!obsatts_valid(&oa, i))
            continue;

        nobs_read++;
        obs_checkalloc(obs);
        o = &obs->data[obs->nobs];

        o->product = productid;
        o->type = typeid;
        o->id = obs->nobs;
        o->id_orig = i;
        o->fid = fid;
        o->value = (double) var[i];
        obsatts_apply(&oa, i, 0, o);
        if (iscurv == 0) {
            o->lon = lon[i % ni];
            o->lat = lat[i / ni];
        } else {
            o->lon = lon[i];
            o->lat = lat[i];
        }
        o->status = grid_xy2fij(g, o->lon, o->lat, o->fij);
        if (!obs->allobs && o->status == STATUS_OUTSIDEGRID)
            continue;
        if (zvalueentered) {
            o->depth = zvalue;
            if (isfinite(zvalue))
                o->status = grid_z2fk_f(g, o->fij, o->depth, &o->fk);
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
    }
    enkf_printf("        nobs = %zu\n", nobs_read);

    free(lon);
    free(lat);
    free(var);
    if (npoints != NULL)
        free(npoints);
    if (time != NULL)
        free(time);
    obsatts_free(&oa);
}

/**
 */
void reader_gridded_xy_describe(void)
{
    enkf_printf("\n  Generic reader \"gridded_xy\" reads 2D gridded data. It can handle either\n\
  curvilinear or geographically rectangular grids. In the latter case the\n\
  coordinates can be represented by 1D variables; longitude is assumed to be\n\
  the inner (\"fast\") coordinate. Currently there can be only 1 time record in\n\
  a data file.\n");
    describe_commongenericreaderparams();
    enkf_printf("  Parameters specific to the reader:\n\
    - ZVALUE (0) (+)\n\
        height/depth (can be NaN)\n\
    - NPOINTSNAME (\"npoints\") (-)\n\
        number of collated points for each datum; used basically as a data mask\n\
        when n = 0\n");
    describe_commonreaderparams();
}
