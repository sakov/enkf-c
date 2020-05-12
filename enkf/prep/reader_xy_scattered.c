/******************************************************************************
 *
 * File:        reader_xy_scattered.c        
 *
 * Created:     21/03/2017
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description: Generic reader for scattered surface observations.
 *                There are a number of parameters that must (++) or can be
 *              specified if they differ from the default value (+). Some
 *              parameters are optional (-):
 *              - VARNAME (++)
 *              - TIMENAME ("*[tT][iI][mM][eE]*") (+)
 *              - or TIMENAMES (when time = base_time + offset) (+)
 *              - LONNAME ("lon" | "longitude") (+)
 *              - LATNAME ("lat" | "latitude") (+)
 *              - STDNAME ("std") (-)
 *                  internal variability of the collated data
 *              - ESTDNAME ("error_std") (-)
 *                  error STD; if absent then needs to be specified externally
 *                  in the oobservation data parameter file
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
 *              - ADDVAR (-)
 *                  name of the variable to be added to the main variable
 *                  (can be repeated)
 *              - SUBVAR (-)
 *                  name of the variable to be subtracted from the main variable
 *                  (can be repeated)
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

#define NADDVAR_INC 1
#define ADDVAR_ACTION_ADD 0
#define ADDVAR_ACTION_SUB 1

typedef struct {
    char* varname;
    int action;
    double* v;
} addvar;

/**
 */
void reader_xy_scattered(char* fname, int fid, obsmeta* meta, grid* g, observations* obs)
{
    int ksurf = grid_getsurflayerid(g);
    char* varname = NULL;
    char* lonname = NULL;
    char* latname = NULL;
    char* stdname = NULL;
    char* estdname = NULL;
    char instrument[MAXSTRLEN] = "";
    int nqcflagvars = 0;
    char** qcflagvarnames = NULL;
    uint32_t* qcflagmasks = 0;
    int naddvar = 0;
    addvar* addvars = NULL;

    int ncid;
    size_t nobs;
    double* lon = NULL;
    double* lat = NULL;
    double* var = NULL;
    double var_estd = NAN;
    double* std = NULL;
    double* estd = NULL;
    uint32_t** qcflag = NULL;
    size_t ntime = 0;
    double* time = NULL;
    int varid;
    size_t i, nobs_read;

    for (i = 0; i < meta->npars; ++i) {
        if (strcasecmp(meta->pars[i].name, "VARNAME") == 0)
            varname = meta->pars[i].value;
        else if (strcasecmp(meta->pars[i].name, "LONNAME") == 0)
            lonname = meta->pars[i].value;
        else if (strcasecmp(meta->pars[i].name, "LATNAME") == 0)
            latname = meta->pars[i].value;
        else if (strcasecmp(meta->pars[i].name, "STDNAME") == 0)
            stdname = meta->pars[i].value;
        else if (strcasecmp(meta->pars[i].name, "ESTDNAME") == 0)
            estdname = meta->pars[i].value;
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
        else if (strcasecmp(meta->pars[i].name, "ADDVAR") == 0) {
            if (naddvar % NADDVAR_INC == 0) {
                addvars = realloc(addvars, (naddvar + NADDVAR_INC) * sizeof(addvar));
                addvars[naddvar].varname = strdup(meta->pars[i].value);
                addvars[naddvar].action = ADDVAR_ACTION_ADD;
                addvars[naddvar].v = NULL;
                naddvar++;
            }
        } else if (strcasecmp(meta->pars[i].name, "SUBVAR") == 0) {
            if (naddvar % NADDVAR_INC == 0) {
                addvars = realloc(addvars, (naddvar + NADDVAR_INC) * sizeof(addvar));
                addvars[naddvar].varname = strdup(meta->pars[i].value);
                addvars[naddvar].action = ADDVAR_ACTION_SUB;
                addvars[naddvar].v = NULL;
                naddvar++;
            }
        } else
            enkf_quit("unknown PARAMETER \"%s\"\n", meta->pars[i].name);
    }

    if (varname == NULL)
        enkf_quit("reader_xy_scattered(): %s: VARNAME not specified", fname);
    else
        enkf_printf("        VARNAME = %s\n", varname);

    for (i = 0; i < naddvar; ++i)
        enkf_printf("        %s = %s\n", (addvars[i].action == ADDVAR_ACTION_ADD) ? "ADDVAR" : "SUBVAR", addvars[i].varname);

    ncw_open(fname, NC_NOWRITE, &ncid);

    /*
     * main variable
     */
    ncw_inq_varid(ncid, varname, &varid);
    ncw_inq_vardims(ncid, varid, 1, NULL, &nobs);
    var = malloc(nobs * sizeof(double));
    ncu_readvardouble(ncid, varid, nobs, var);

    /*
     * add variables
     */
    for (i = 0; i < naddvar; ++i) {
        addvar* a = &addvars[i];

        ncw_inq_varid(ncid, a->varname, &varid);
        a->v = malloc(nobs * sizeof(double));
        ncu_readvardouble(ncid, varid, nobs, a->v);
        if (a->action == ADDVAR_ACTION_SUB) {
            int ii;

            for (ii = 0; ii < nobs; ++ii)
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
        ncw_check_vardims(ncid, varid, 1, &nobs);
    } else
        enkf_quit("reader_xy_scattered(): %s: could not find longitude variable", fname);
    lon = malloc(nobs * sizeof(double));
    ncu_readvardouble(ncid, varid, nobs, lon);

    /*
     * latitude
     */
    latname = get_latname(ncid, latname);
    if (latname != NULL) {
        enkf_printf("        LATNAME = %s\n", latname);
        ncw_inq_varid(ncid, latname, &varid);
        ncw_check_vardims(ncid, varid, 1, &nobs);
    } else
        enkf_quit("reader_xy_scattered(): %s: could not find latitude variable", fname);
    lat = malloc(nobs * sizeof(double));
    ncu_readvardouble(ncid, varid, nobs, lat);

    /*
     * std
     */
    varid = -1;
    if (stdname != NULL)
        ncw_inq_varid(ncid, stdname, &varid);
    else if (ncw_var_exists(ncid, "std"))
        ncw_inq_varid(ncid, "std", &varid);
    if (varid >= 0) {
        ncw_check_vardims(ncid, varid, 1, &nobs);
        std = malloc(nobs * sizeof(double));
        ncu_readvardouble(ncid, varid, nobs, std);
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
        ncw_check_vardims(ncid, varid, 1, &nobs);
        estd = malloc(nobs * sizeof(double));
        ncu_readvardouble(ncid, varid, nobs, estd);
    }

    if (std == NULL && estd == NULL) {
        ncw_inq_varid(ncid, varname, &varid);
        if (ncw_att_exists(ncid, varid, "error_std")) {
            ncw_check_attlen(ncid, varid, "error_std", 1);
            ncw_get_att_double(ncid, varid, "error_std", &var_estd);
        }
    }

    /*
     * qcflag
     */
    get_qcflags(meta, &nqcflagvars, &qcflagvarnames, &qcflagmasks);
    if (nqcflagvars > 0) {
        int varid = -1;

        qcflag = alloc2d(nqcflagvars, nobs, sizeof(int32_t));
        for (i = 0; i < nqcflagvars; ++i) {
            ncw_inq_varid(ncid, qcflagvarnames[i], &varid);
            ncw_check_vardims(ncid, varid, 1, &nobs);
            ncw_get_var_uint(ncid, varid, qcflag[i]);
        }
    }

    /*
     * time
     */
    get_time(meta, ncid, &ntime, &time);
    assert(ntime == nobs || ntime <= 1);

    /*
     * instrument
     */
    if (strlen(instrument) == 0 && !get_insttag(ncid, varname, instrument))
        strncpy(instrument, meta->product, MAXSTRLEN - 1);

    ncw_close(ncid);

    nobs_read = 0;
    for (i = 0; i < nobs; ++i) {
        observation* o;
        int ii;

        if (isnan(lon[i]) || isnan(lat[i]))
            continue;
        if (isnan(var[i]) || (std != NULL && isnan(std[i])) || (estd != NULL && isnan(estd[i])) || (ntime == nobs && isnan(time[i])))
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
        o->value = var[i];
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
        o->lon = lon[i];
        o->lat = lat[i];
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
