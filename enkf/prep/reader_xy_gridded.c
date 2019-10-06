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
 *                It currently assumes the following:
 *              - there is only one data record (2D field);
 *              - longitude is the inner ("fast") coordinate of the variable.
 *                There are a number of parameters that must (++) or can be
 *              specified if they differ from the default value (+). Some
 *              parameters are optional (-):
 *              - VARNAME (++)
 *              - TIMENAME ("*[tT][iI][mM][eE]*") (+)
 *              - or TIMENAMES (when time = base_time + offset) (+)
 *              - NPOINTSNAME ("npoints") (-)
 *                  number of collated points for each datum; used basically as
 *                  a data mask n = 0
 *              - LONNAME ("lon" | "longitude") (+)
 *              - LATNAME ("lat" | "latitude") (+)
 *              - STDNAME ("std") (-)
 *                  internal variability of the collated data
 *              - ESTDNAME ("error_std") (-)
 *                  error STD; if absent then needs to be specified externally
 *                  in the observation data parameter file
 *              - VARSHIFT (-)
 *                  data offset to be added
 *              - FOOTRPINT (-)
 *                  footprint of observations in km
 *              - MINDEPTH (-)
 *                  minimal allowed depth
 *              - MAXDEPTH (-)
 *                  maximal allowed depth
 *              - INSTRUMENT (-)
 *                  instrument string that will be used for calculating
 *                  instrument stats
 *              - QCFLAGNAME (-)
 *                  name of the QC flag variable, 0 <= qcflag <= 31
 *              - QCFLAGVALS (-)
 *                  the list of allowed values of QC flag variable
 *              Note: it is possible to have multiple entries of QCFLAGNAME and
 *                QCFLAGVALS combination, e.g.:
 *                  PARAMETER QCFLAGNAME = TEMP_quality_control
 *                  PARAMETER QCFLAGVALS = 1
 *                  PARAMETER QCFLAGNAME = DEPTH_quality_control
 *                  PARAMETER QCFLAGVALS = 1
 *                  PARAMETER QCFLAGNAME = LONGITUDE_quality_control
 *                  PARAMETER QCFLAGVALS = 1,8
 *                  PARAMETER QCFLAGNAME = LATITUDE_quality_control
 *                  PARAMETER QCFLAGVALS = 1,8
 *                An observation is considered valid if each of the specified
 *                flags takes a permitted value.
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
void reader_xy_gridded(char* fname, int fid, obsmeta* meta, grid* g, observations* obs)
{
    int ksurf = grid_getsurflayerid(g);
    char* varname = NULL;
    char* lonname = NULL;
    char* latname = NULL;
    char* npointsname = NULL;
    char* stdname = NULL;
    char* estdname = NULL;
    int ndim_var, ndim_xy;
    size_t dimlen_var[3], dimlen_xy[2];
    int nqcflagvars = 0;
    char** qcflagvarnames = NULL;
    uint32_t* qcflagmasks = 0;

    int ncid;
    float varshift = 0.0;
    char instrument[MAXSTRLEN] = "";

    int iscurv = -1;
    size_t ni = 0, nj = 0, nij = 0, n_var = 0;
    int varid_lon = -1, varid_lat = -1;
    double* lon = NULL;
    double* lat = NULL;
    int varid_var = -1, varid_npoints = -1, varid_std = -1, varid_estd = -1;
    float* var = NULL;
    double var_estd = NAN;
    short* npoints = NULL;
    float* std = NULL;
    float* estd = NULL;
    uint32_t** qcflag = NULL;
    size_t ntime = 0;
    double* time = NULL;
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
        else if (strcasecmp(meta->pars[i].name, "STDNAME") == 0)
            stdname = meta->pars[i].value;
        else if (strcasecmp(meta->pars[i].name, "ESTDNAME") == 0)
            estdname = meta->pars[i].value;
        else if (strcasecmp(meta->pars[i].name, "VARSHIFT") == 0) {
            if (!str2float(meta->pars[i].value, &varshift))
                enkf_quit("%s: can not convert VARSHIFT = \"%s\" to float\n", meta->prmfname, meta->pars[i].value);
            enkf_printf("        VARSHIFT = %s\n", meta->pars[i].value);
        }
        /*
         * (FOOTPRINT, MINDEPTH and MAXDEPTH are handled in obs_add() )
         */
        else if (strcasecmp(meta->pars[i].name, "FOOTPRINT") == 0) {
            double footprint;

            if (!str2double(meta->pars[i].value, &footprint))
                enkf_quit("observation prm file: can not convert FOOTPRINT = \"%s\" to double\n", meta->pars[i].value);
            enkf_printf("        FOOTPRINT = %.0f\n", footprint);
            continue;
        } else if (strcasecmp(meta->pars[i].name, "MINDEPTH") == 0) {
            double mindepth;

            if (!str2double(meta->pars[i].value, &mindepth))
                enkf_quit("observation prm file: can not convert MINDEPTH = \"%s\" to double\n", meta->pars[i].value);
            enkf_printf("        MINDEPTH = %.0f\n", mindepth);
            continue;
        } else if (strcasecmp(meta->pars[i].name, "MAXDEPTH") == 0) {
            double maxdepth;

            if (!str2double(meta->pars[i].value, &maxdepth))
                enkf_quit("observation prm file: can not convert MAXDEPTH = \"%s\" to double\n", meta->pars[i].value);
            enkf_printf("        MAXDEPTH = %.0f\n", maxdepth);
            continue;
        } else if (strcasecmp(meta->pars[i].name, "INSTRUMENT") == 0)
            strncpy(instrument, meta->pars[i].value, MAXSTRLEN - 1);
        else if (strcasecmp(meta->pars[i].name, "TIMENAME") == 0 || strcasecmp(meta->pars[i].name, "TIMENAMES") == 0)
            /*
             * TIMENAME and TIMENAMES are dealt with separately
             */
            ;
        else if (strcasecmp(meta->pars[i].name, "QCFLAGNAME") == 0 || strcasecmp(meta->pars[i].name, "QCFLAGVALS") == 0)
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
     * main variable dims
     */
    ncw_inq_varid(ncid, varname, &varid_var);
    ncw_inq_vardims(ncid, varid_var, 3, &ndim_var, dimlen_var);
    if (ndim_var == 3) {
        if (dimlen_var[0] != 1)
            enkf_quit("reader_xy_gridded(): %s: %s: %d records (currently only one is allowed)", fname, varname, dimlen_var[0]);
        n_var = dimlen_var[1] * dimlen_var[2];
    } else if (ndim_var == 2) {
        if (ncw_var_hasunlimdim(ncid, varid_var))
            enkf_quit("reader_xy_gridded(): %s: %s: not enough spatial dimensions (must be 2)", fname, varname);
        n_var = dimlen_var[0] * dimlen_var[1];
    } else if (ndim_var != 2)
        enkf_quit("reader_xy_gridded(): %s: %s: %d dimensions (must be 2 or 3 with only one record)", fname, varname, ndim_var);

    /*
     * longitude dims
     */
    lonname = get_lonname(ncid, lonname);
    if (lonname != NULL) {
        enkf_printf("        LONNAME = %s\n", lonname);
        ncw_inq_varid(ncid, lonname, &varid_lon);
    } else
        enkf_quit("reader_xy_gridded(): %s: could not find longitude variable", fname);
    ncw_inq_vardims(ncid, varid_lon, 2, &ndim_xy, dimlen_xy);
    if (ndim_xy == 1) {
        iscurv = 0;
        ni = dimlen_xy[0];
    } else if (ndim_xy == 2) {
        iscurv = 1;
        ni = dimlen_xy[1];
        nj = dimlen_xy[0];
    } else
        enkf_quit("reader_xy_gridded(): %s: coordinate variable \"%s\" has neither 1 or 2 dimensions", fname, lonname);

    /*
     * latitude dims
     */
    latname = get_latname(ncid, latname);
    if (latname != NULL) {
        enkf_printf("        LATNAME = %s\n", latname);
        ncw_inq_varid(ncid, latname, &varid_lat);
    } else
        enkf_quit("reader_xyz_gridded(): %s: could not find latitude variable", fname);
    if (ndim_xy == 1)
        ncw_inq_vardims(ncid, varid_lat, 1, &ndim_xy, &nj);

    if (iscurv == 0) {
        ncw_check_varndims(ncid, varid_lat, 1);
        ncw_inq_vardims(ncid, varid_lat, 1, NULL, &nj);
    } else
        ncw_check_vardims(ncid, varid_lat, 2, dimlen_xy);

    enkf_printf("        (ni, nj) = (%u, %u)\n", ni, nj);
    nij = ni * nj;
    if (nij != n_var)
        enkf_quit("reader_xy_gridded(): %s: dimensions of variable \"%s\" do not match coordinate dimensions", fname, varname);
    if (dimlen_var[ndim_var - 1] != ni)
        enkf_quit("reader_xy_gridded(): %s: %s: longitude must be the inner coordinate", fname, varname);

    /*
     * lon/lat
     */
    if (iscurv == 0) {
        lon = malloc(ni * sizeof(double));
        ncu_readvardouble(ncid, varid_lon, ni, lon);

        lat = malloc(nj * sizeof(double));
        ncu_readvardouble(ncid, varid_lat, nj, lat);
    } else {
        lon = malloc(nij * sizeof(double));
        ncu_readvardouble(ncid, varid_lon, nij, lon);

        lat = malloc(nij * sizeof(double));
        ncu_readvardouble(ncid, varid_lat, nij, lat);
    }

    /*
     * main variable
     */
    var = malloc(nij * sizeof(float));
    ncu_readvarfloat(ncid, varid_var, nij, var);

    /*
     * npoints
     */
    if (npointsname != NULL)
        ncw_inq_varid(ncid, npointsname, &varid_npoints);
    else if (ncw_var_exists(ncid, "npoints"))
        ncw_inq_varid(ncid, "npoints", &varid_npoints);
    if (varid_npoints >= 0) {
        npoints = malloc(nij * sizeof(short));
        ncw_get_var_short(ncid, varid_npoints, npoints);
    }

    /*
     * std
     */
    if (stdname != NULL)
        ncw_inq_varid(ncid, stdname, &varid_std);
    else if (ncw_var_exists(ncid, "std"))
        ncw_inq_varid(ncid, "std", &varid_std);
    if (varid_std >= 0) {
        std = malloc(nij * sizeof(float));
        ncu_readvarfloat(ncid, varid_std, nij, std);
    }

    /*
     * estd
     */
    if (estdname != NULL)
        ncw_inq_varid(ncid, estdname, &varid_estd);
    else if (ncw_var_exists(ncid, "error_std"))
        ncw_inq_varid(ncid, "error_std", &varid_estd);
    if (varid_estd >= 0) {
        estd = malloc(nij * sizeof(float));
        ncu_readvarfloat(ncid, varid_estd, nij, estd);
    }

    if (std == NULL && estd == NULL) {
        if (ncw_att_exists(ncid, varid_var, "error_std")) {
            ncw_check_attlen(ncid, varid_var, "error_std", 1);
            ncw_get_att_double(ncid, varid_var, "error_std", &var_estd);
        }
    }

    /*
     * qcflags
     */
    get_qcflags(meta, &nqcflagvars, &qcflagvarnames, &qcflagmasks);
    if (nqcflagvars > 0) {
        int varid = -1;

        qcflag = alloc2d(nqcflagvars, nij, sizeof(int32_t));
        for (i = 0; i < nqcflagvars; ++i) {
            ncw_inq_varid(ncid, qcflagvarnames[i], &varid);
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

        o->product = st_findindexbystring(obs->products, meta->product);
        assert(o->product >= 0);
        o->type = obstype_getid(obs->nobstypes, obs->obstypes, meta->type, 1);
        o->instrument = st_add_ifabsent(obs->instruments, instrument, -1);
        o->id = obs->nobs;
        o->fid = fid;
        o->batch = 0;
        o->value = (double) (var[i] + varshift);
        if (estd == NULL)
            o->estd = var_estd;
        else {
            o->estd = (std == NULL) ? 0.0 : (double) std[i];
            o->estd = (o->estd > estd[i]) ? o->estd : estd[i];
        }
        if (iscurv == 0) {
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
