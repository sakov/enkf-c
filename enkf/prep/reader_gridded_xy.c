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

#define NADDVAR_INC 1
#define ADDVAR_ACTION_ADD 0
#define ADDVAR_ACTION_SUB 1

typedef struct {
    char* varname;
    int action;
    float* v;
} addvar;

#define TYPE_DOUBLE 0
#define TYPE_SHORT 1

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
    char* stdname = NULL;
    char* estdname = NULL;
    char* batchname = NULL;
    char* instattname = NULL;
    char* instprefix = NULL;
    char instrument[MAXSTRLEN] = "";
    int ndim_var, ndim_xy;
    size_t dimlen_var[3], dimlen_xy[2];
    int nqcflagvars = 0;
    char** qcflagvarnames = NULL;
    uint32_t* qcflagmasks = NULL;

    int naddvar = 0;
    addvar* addvars = NULL;

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
    int32_t** qcflag = NULL;
    size_t ntime = 0;
    double* time = NULL;
    int varid;
    size_t i, nobs_read;

    for (i = 0; i < section->npars; ++i) {
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
        } else if (strcasecmp(section->pars[i].name, "STDNAME") == 0)
            stdname = section->pars[i].value;
        else if (strcasecmp(section->pars[i].name, "ESTDNAME") == 0)
            estdname = section->pars[i].value;
        else if (strcasecmp(section->pars[i].name, "BATCHNAME") == 0)
            batchname = section->pars[i].value;
        else if (strcasecmp(section->pars[i].name, "INSTRUMENT") == 0)
            strncpy(instrument, section->pars[i].value, MAXSTRLEN - 1);
        else if (strcasecmp(section->pars[i].name, "TIMENAME") == 0 || strcasecmp(section->pars[i].name, "TIMENAMES") == 0)
            /*
             * TIMENAME and TIMENAMES are dealt with later
             */
            ;
        else if (strcasecmp(section->pars[i].name, "QCFLAGNAME") == 0 || strcasecmp(section->pars[i].name, "QCFLAGVARNAME") == 0 || strcasecmp(section->pars[i].name, "QCFLAGVALS") == 0)
            /*
             * QCFLAGNAME and QCFLAGVALS are dealt with later
             */
            ;
        else if (strcasecmp(section->pars[i].name, "LOCATION_BASED_THINNING_TYPE") == 0)
            /*
             * LOCATION_BASED_THINNING_TYPE is dealt with outside
             */
            ;
        else if (strcasecmp(section->pars[i].name, "ADDVAR") == 0) {
            if (naddvar % NADDVAR_INC == 0) {
                addvars = realloc(addvars, (naddvar + NADDVAR_INC) * sizeof(addvar));
                addvars[naddvar].varname = strdup(section->pars[i].value);
                addvars[naddvar].action = ADDVAR_ACTION_ADD;
                addvars[naddvar].v = NULL;
                enkf_printf("      ADDVAR = %s\n", addvars[naddvar].varname);
                naddvar++;
            }
        } else if (strcasecmp(section->pars[i].name, "SUBVAR") == 0) {
            if (naddvar % NADDVAR_INC == 0) {
                addvars = realloc(addvars, (naddvar + NADDVAR_INC) * sizeof(addvar));
                addvars[naddvar].varname = strdup(section->pars[i].value);
                addvars[naddvar].action = ADDVAR_ACTION_SUB;
                addvars[naddvar].v = NULL;
                enkf_printf("      SUBVAR = %s\n", addvars[naddvar].varname);
                naddvar++;
            }
        } else
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
     * add variables
     */
    for (i = 0; i < naddvar; ++i) {
        addvar* a = &addvars[i];

        ncw_inq_varid(ncid, a->varname, &varid);
        a->v = malloc(nij * sizeof(float));
        ncu_readvarfloat(ncid, varid, nij, a->v);
        if (a->action == ADDVAR_ACTION_SUB) {
            int ii;

            for (ii = 0; ii < nij; ++ii)
                a->v[ii] = -a->v[ii];
        }
    }

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
    get_qcflags(section, &nqcflagvars, &qcflagvarnames, &qcflagmasks);
    if (nqcflagvars > 0) {
        qcflag = alloc2d(nqcflagvars, nij, sizeof(int32_t));
        for (i = 0; i < nqcflagvars; ++i) {
            ncw_inq_varid(ncid, qcflagvarnames[i], &varid);
            ncw_check_varsize(ncid, varid, nij);
            ncw_get_var_int(ncid, varid, qcflag[i]);
        }
    }

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

    /*
     * instrument
     */
    if (strlen(instrument) == 0) {
        if (instattname != NULL)
            ncw_get_att_text(ncid, NC_GLOBAL, instattname, instrument);
        else if (!get_insttag(ncid, varname, instrument))
            strncpy(instrument, section->product, MAXSTRLEN - 1);
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
    productid = st_findindexbystring(obs->products, section->product);
    assert(productid >= 0);
    typeid = obstype_getid(obs->nobstypes, obs->obstypes, section->type, 1);
    nobs_read = 0;
    for (i = 0; i < nij; ++i) {
        observation* o;
        int ii;

        if ((npoints != NULL && npoints[i] == 0) || isnan(var[i]) || (std != NULL && isnan(std[i])) || (estd != NULL && isnan(estd[i])) || (ntime == nij && isnan(time[i])))
            continue;
        for (ii = 0; ii < nqcflagvars; ++ii)
            if (qcflag[ii][i] < 0 || qcflag[ii][i] > 31 || !((1 << qcflag[ii][i]) & qcflagmasks[ii]))
                goto nextob;

        nobs_read++;
        obs_checkalloc(obs);
        o = &obs->data[obs->nobs];

        o->product = productid;
        o->type = typeid;
        o->instrument = instid;
        o->id = obs->nobs;
        o->id_orig = i;
        o->fid = fid;
        o->batch = (batch == NULL) ? 0 : batch[i];
        o->value = (double) var[i];
        for (ii = 0; ii < naddvar; ++ii)
            o->value += addvars[ii].v[i];
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
    if (naddvar > 0) {
        for (i = 0; i < naddvar; ++i) {
            free(addvars[i].varname);
            free(addvars[i].v);
        }
        free(addvars);
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
  a data file.\n");
    describe_commongenericreaderparams();
    enkf_printf("  Parameters specific to the reader:\n\
    - ZVALUE (0) (+)\n\
        height/depth (can be NaN)\n\
    - NPOINTSNAME (\"npoints\") (-)\n\
        number of collated points for each datum; used basically as a data mask\n\
        when n = 0\n\
    - ADDVAR (-)\n\
        name of the variable to be added to the main variable (can be repeated)\n\
    - SUBVAR (-)\n\
        name of the variable to be subtracted from the main variable (can be repeated)\n");
    describe_commonreaderparams();
}
