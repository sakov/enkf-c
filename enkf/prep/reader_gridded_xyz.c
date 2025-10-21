/******************************************************************************
 *
 * File:        reader_gridded_xyz.c        
 *
 * Created:     25/08/2017
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description: Generic reader for 3D gridded observations.
 *
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
void reader_gridded_xyz(char* fname, int fid, obsmeta* meta, grid* g, observations* obs)
{
    char* varname = NULL;
    char* lonname = NULL;
    char* latname = NULL;
    char* zname = NULL;
    char* npointsname = NULL;
    char* stdname = NULL;
    char* estdname = NULL;
    char* batchname = NULL;
    char* instattname = NULL;
    char* instprefix = NULL;
    char instrument[MAXSTRLEN] = "";
    int nqcflagvars = 0;
    char** qcflagvarnames = NULL;
    uint32_t* qcflagmasks = NULL;

    int instid = -1;
    int productid = -1;
    int typeid = -1;

    int ncid;
    int iscurv = -1, zndim = -1;
    int ndim_var, ndim_xy;
    size_t dimlen_var[4], dimlen_xy[2], dimlen_z[3];
    size_t ni = 0, nj = 0, nk = 0, nij = 0, nijk = 0;
    double* lon = NULL;
    double* lat = NULL;
    float* z = NULL;

    float* var = NULL;
    double var_estd = NAN;
    short* npoints = NULL;
    float* std = NULL;
    float* estd = NULL;
    int* batch = NULL;
    int32_t** qcflag = NULL;
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
        else if (strcasecmp(meta->pars[i].name, "ZNAME") == 0)
            zname = meta->pars[i].value;
        else if (strcasecmp(meta->pars[i].name, "STDNAME") == 0)
            stdname = meta->pars[i].value;
        else if (strcasecmp(meta->pars[i].name, "ESTDNAME") == 0)
            estdname = meta->pars[i].value;
        else if (strcasecmp(meta->pars[i].name, "BATCHNAME") == 0)
            batchname = meta->pars[i].value;
        else if (strcasecmp(meta->pars[i].name, "INSTRUMENT") == 0)
            strncpy(instrument, meta->pars[i].value, MAXSTRLEN - 1);
        else if (strcasecmp(meta->pars[i].name, "TIMENAME") == 0 || strcasecmp(meta->pars[i].name, "TIMENAMES") == 0)
            /*
             * TIMENAME and TIMENAMES are dealt with later
             */
            ;
        else if (strcasecmp(meta->pars[i].name, "QCFLAGNAME") == 0 || strcasecmp(meta->pars[i].name, "QCFLAGVARNAME") == 0 || strcasecmp(meta->pars[i].name, "QCFLAGVALS") == 0)
            /*
             * QCFLAGNAME and QCFLAGVALS are dealt with later
             */
            ;
        else if (strcasecmp(meta->pars[i].name, "LOCATION_BASED_THINNING_TYPE") == 0)
            /*
             * LOCATION_BASED_THINNING_TYPE is dealt with outside
             */
            ;
        else
            enkf_quit("unknown PARAMETER \"%s\"\n", meta->pars[i].name);
    }

    if (varname == NULL)
        enkf_quit("reader_xyz_gridded(): %s: VARNAME not specified", fname);
    else
        enkf_printf("        VARNAME = %s\n", varname);

    ncw_open(fname, NC_NOWRITE, &ncid);

    /*
     * main variable
     */
    ncw_inq_varid(ncid, varname, &varid);
    ncw_inq_vardims(ncid, varid, 4, &ndim_var, dimlen_var);
    if (ndim_var == 4) {
        if (dimlen_var[0] != 1)
            enkf_quit("reader_xyz_gridded(): %s: %s: %d records (currently only one record is allowed)", fname, varname, dimlen_var[0]);
        nijk = dimlen_var[3] * dimlen_var[2] * dimlen_var[1];
        nij = dimlen_var[3] * dimlen_var[2];
    } else if (ndim_var == 3) {
        if (ncw_var_hasunlimdim(ncid, varid))
            enkf_quit("reader_xyz_gridded(): %s: %s: not enough spatial dimensions (must be 3)", fname, varname);
        nijk = dimlen_var[2] * dimlen_var[1] * dimlen_var[0];
        nij = dimlen_var[2] * dimlen_var[1];
    } else
        enkf_quit("reader_xyz_gridded(): %s: %s: %d dimensions (must be 3 or 4 with only one record)", fname, varname, ndim_var);

    var = malloc(nijk * sizeof(float));
    ncu_readvarfloat(ncid, varid, nijk, var);

    /*
     * longitude
     */
    lonname = get_lonname(ncid, lonname);
    if (lonname != NULL) {
        enkf_printf("        LONNAME = %s\n", lonname);
        ncw_inq_varid(ncid, lonname, &varid);
    } else
        enkf_quit("reader_xyz_gridded(): %s: could not find longitude variable", fname);

    ncw_inq_vardims(ncid, varid, 2, &ndim_xy, dimlen_xy);
    if (ndim_xy == 1) {
        iscurv = 0;
        ni = dimlen_xy[0];
    } else if (ndim_xy == 2) {
        iscurv = 1;
        ni = dimlen_xy[1];
        nj = dimlen_xy[0];
    } else
        enkf_quit("reader_xyz_gridded(): %s: coordinate variable \"%s\" has neither 1 or 2 dimensions", fname, lonname);

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

    if (ni * nj != nij)
        enkf_quit("reader_xyz_gridded(): %s: dimensions of variable \"%s\" do not match coordinate dimensions", fname, varname);

    if (iscurv == 0) {
        lat = malloc(nj * sizeof(double));
        ncu_readvardouble(ncid, varid, nj, lat);
    } else {
        lat = malloc(nij * sizeof(double));
        ncu_readvardouble(ncid, varid, nij, lat);
    }

    /*
     * z
     */
    zname = get_zname(ncid, zname);
    if (zname != NULL) {
        enkf_printf("        ZNAME = %s\n", zname);
        ncw_inq_varid(ncid, zname, &varid);
    } else
        enkf_quit("reader_xyz_gridded(): %s: could not find Z variable", fname);

    ncw_inq_vardims(ncid, varid, 3, &zndim, dimlen_z);
    if (zndim == 1)
        nk = dimlen_z[0];
    else if (zndim == 3)
        nk = dimlen_z[0];
    else
        enkf_quit("reader_xyz_gridded(): %s: %s (the vertical coordinate): %d-dimensional; supposed to be either 1- or 3-dimensional only", fname, zname, zndim);

    enkf_printf("        (ni, nj, nk) = (%u, %u, %u)\n", ni, nj, nk);
    if (ni * nj * nk != nijk)
        enkf_quit("reader_xyz_gridded(): %s: dimensions of variable \"%s\" do not match coordinate dimensions", fname, varname);
    if (dimlen_var[ndim_var - 1] != ni)
        enkf_quit("reader_xyz_gridded(): %s: %s: longitude must be the inner coordinate", fname, varname);

    if (zndim == 1) {
        z = malloc(nk * sizeof(float));
        ncu_readvarfloat(ncid, varid, nk, z);
    } else if (zndim == 3) {
        z = malloc(nijk * sizeof(float));
        ncu_readvarfloat(ncid, varid, nijk, z);
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
        npoints = malloc(nijk * sizeof(short));
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
        std = malloc(nijk * sizeof(float));
        ncu_readvarfloat(ncid, varid, nijk, std);
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
        estd = malloc(nijk * sizeof(float));
        ncu_readvarfloat(ncid, varid, nijk, estd);
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
        ncw_check_varsize(ncid, varid, nijk);
        batch = malloc(nijk * sizeof(int));
        ncw_get_var_int(ncid, varid, batch);
    }

    /*
     * qcflag
     */
    get_qcflags(meta, &nqcflagvars, &qcflagvarnames, &qcflagmasks);
    if (nqcflagvars > 0) {
        qcflag = alloc2d(nqcflagvars, nijk, sizeof(int32_t));
        for (i = 0; i < nqcflagvars; ++i) {
            ncw_inq_varid(ncid, qcflagvarnames[i], &varid);
            ncw_check_varsize(ncid, varid, nijk);
            ncw_get_var_int(ncid, varid, qcflag[i]);
        }
    }

    /*
     * time
     */
    get_time(meta, ncid, &ntime, &time);
    assert(ntime == nijk || ntime <= 1);

    /*
     * instrument
     */
    if (strlen(instrument) == 0) {
        if (instattname != NULL)
            ncw_get_att_text(ncid, NC_GLOBAL, instattname, instrument);
        else if (!get_insttag(ncid, varname, instrument))
            strncpy(instrument, meta->product, MAXSTRLEN - 1);
        if (strlen(instrument) > 0 && instprefix != NULL) {
            int len_p = strlen(instprefix);
            int len_i = strlen(instrument);
            int ii;

            assert(len_p + len_i < MAXSTRLEN);
            for (ii = len_i - 1; ii >= 0; --ii)
                instrument[ii + len_p] = instrument[ii];
            for (ii = 0; ii < len_p; ++ii)
                instrument[ii] = instprefix[ii];
            instrument[len_p + len_i] = 0;
        }
    }
    ncw_close(ncid);

    instid = st_add_ifabsent(obs->instruments, instrument, -1);
    productid = st_findindexbystring(obs->products, meta->product);
    assert(productid >= 0);
    typeid = obstype_getid(obs->nobstypes, obs->obstypes, meta->type, 1);
    nobs_read = 0;
    for (i = 0; i < nijk; ++i) {
        int ij = i % nij;
        observation* o;
        int ii;

        if ((npoints != NULL && npoints[i] == 0) || isnan(var[i]) || (std != NULL && isnan(std[i])) || (estd != NULL && isnan(estd[i])) || (ntime == nijk && isnan(time[i])))
            continue;
        for (ii = 0; ii < nqcflagvars; ++ii)
            if (qcflag[ii][i] < 0 || qcflag[ii][i] > 31 || !((1 << qcflag[ii][i]) & qcflagmasks[ii]))
                goto nextob;

        nobs_read++;
        obs_checkalloc(obs);
        o = &obs->data[obs->nobs];

        o->product = productid;
        assert(o->product >= 0);
        o->type = typeid;
        o->instrument = instid;
        o->id = obs->nobs;
        o->id_orig = i;
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
            o->lon = lon[ij % ni];
            o->lat = lat[ij / ni];
        } else {
            o->lon = lon[ij];
            o->lat = lat[ij];
        }
        o->status = grid_xy2fij(g, o->lon, o->lat, o->fij);
        if (!obs->allobs && o->status == STATUS_OUTSIDEGRID)
            continue;
        o->depth = (zndim == 1) ? z[i / nij] : z[i];
        if (o->status == STATUS_OK)
            o->status = grid_z2fk_f(g, o->fij, o->depth, &o->fk);
        else
            o->fk = NAN;
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
    free(z);
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
void reader_gridded_xyz_describe(void)
{
    enkf_printf("\n  Generic reader \"gridded_xyz\" reads 3D gridded data. It can handle either\n\
  curvilinear or geographically rectangular grids. In the latter case the\n\
  coordinates can be represented by 1D variables; longitude is assumed to be\n\
  the inner (\"fast\") coordinate, followed by latitude, and then depth\n\
  (height). Currently there can be only 1 time record in a data file.\n");
    describe_commongenericreaderparams();
    enkf_printf("  Parameters specific to the reader:\n\
    - NPOINTSNAME (\"npoints\") (-)\n\
        number of collated points for each datum; used basically as a data mask\n\
        when n = 0\n");
    describe_commonreaderparams();
}
