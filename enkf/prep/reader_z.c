/******************************************************************************
 *
 * File:        reader_z.c        
 *
 * Created:     25/06/2019
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description: Generic reader for observations in a single profile.
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

#define TYPE_DOUBLE 0
#define TYPE_SHORT 1

/**
 */
void reader_z(char* fname, int fid, obsmeta* meta, grid* g, observations* obs)
{
    char* varname = NULL;
    char* lonname = NULL;
    char* latname = NULL;
    char* zname = NULL;
    char* estdname = NULL;
    char instrument[MAXSTRLEN] = "";

    int nqcflagvars = 0;
    char** qcflagvarnames = NULL;
    uint32_t* qcflagmasks = NULL;

    int instid = -1;
    int productid = -1;
    int typeid = -1;

    int ncid;
    size_t nobs = 0;
    double lon = NAN;
    double lat = NAN;
    double* var = NULL;
    double* z = NULL;
    double var_estd = NAN;
    double* estd = NULL;
    uint32_t** qcflag = NULL;
    size_t ntime = 0;
    double* time = NULL;
    int varid;
    int i, nobs_read;

    for (i = 0; i < meta->npars; ++i) {
        if (strcasecmp(meta->pars[i].name, "VARNAME") == 0)
            varname = meta->pars[i].value;
        else if (strcasecmp(meta->pars[i].name, "LONNAME") == 0)
            lonname = meta->pars[i].value;
        else if (strcasecmp(meta->pars[i].name, "LATNAME") == 0)
            latname = meta->pars[i].value;
        else if (strcasecmp(meta->pars[i].name, "ZNAME") == 0)
            zname = meta->pars[i].value;
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
        else
            enkf_quit("unknown PARAMETER \"%s\"\n", meta->pars[i].name);
    }

    if (varname == NULL)
        enkf_quit("reader_z(): %s: VARNAME not specified", fname);
    else
        enkf_printf("        VARNAME = %s\n", varname);

    ncw_open(fname, NC_NOWRITE, &ncid);
    ncw_inq_varid(ncid, varname, &varid);
    ncw_inq_varsize(ncid, varid, &nobs);
    /*
     * main variable
     */
    var = malloc(nobs * sizeof(double));
    ncu_readvardouble(ncid, varid, nobs, var);

    /*
     * lon/lat
     */
    lonname = get_lonname(ncid, lonname);
    if (lonname != NULL) {
        enkf_printf("        LONNAME = %s\n", lonname);
        ncw_inq_varid(ncid, lonname, &varid);
    } else
        enkf_quit("reader_z(): %s: could not find longitude variable", fname);
    ncw_check_varsize(ncid, varid, 1);
    ncw_get_var_double(ncid, varid, &lon);

    latname = get_latname(ncid, latname);
    if (latname != NULL) {
        enkf_printf("        LATNAME = %s\n", latname);
        ncw_inq_varid(ncid, latname, &varid);
    } else
        enkf_quit("reader_z(): %s: could not find latitude variable", fname);
    ncw_check_varsize(ncid, varid, 1);
    ncw_get_var_double(ncid, varid, &lat);

    /*
     * z
     */
    zname = get_zname(ncid, zname);
    if (zname != NULL) {
        enkf_printf("        ZNAME = %s\n", zname);
        ncw_inq_varid(ncid, zname, &varid);
    } else
        enkf_quit("reader_z(): %s: could not find Z variable", fname);
    ncw_check_varsize(ncid, varid, nobs);
    z = malloc(nobs * sizeof(double));
    ncu_readvardouble(ncid, varid, nobs, z);

    /*
     * estd
     */
    varid = -1;
    if (estdname != NULL)
        ncw_inq_varid(ncid, estdname, &varid);
    else if (ncw_var_exists(ncid, "error_std"))
        ncw_inq_varid(ncid, "error_std", &varid);
    if (varid >= 0) {
        estd = malloc(nobs * sizeof(double));
        ncu_readvardouble(ncid, varid, nobs, estd);
    }

    if (estd == NULL) {
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

    instid = st_add_ifabsent(obs->instruments, instrument, -1);
    productid = st_findindexbystring(obs->products, meta->product);
    assert(productid >= 0);
    typeid = obstype_getid(obs->nobstypes, obs->obstypes, meta->type, 1);
    nobs_read = 0;
    for (i = 0; i < (int) nobs; ++i) {
        observation* o;
        int ii;

        if (isnan(var[i]) || (estd != NULL && isnan(estd[i])) || (ntime == nobs && isnan(time[i])))
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
        o->batch = 0;
        o->value = (double) var[i];
        if (estd == NULL)
            o->estd = var_estd;
        else
            o->estd = estd[i];
        o->lon = lon;
        o->lat = lat;
        o->depth = z[i];
        o->status = grid_xy2fij_f(g, o->lon, o->lat, &o->fi, &o->fj);
        if (!obs->allobs && o->status == STATUS_OUTSIDEGRID)
            continue;
        if (o->status == STATUS_OK)
            o->status = grid_z2fk_f(g, o->fi, o->fj, o->depth, &o->fk);
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

    free(var);
    free(z);
    if (estd != NULL)
        free(estd);
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
void reader_z_describe(void)
{
    enkf_printf("\n  Generic reader \"z\" reads data of a single profile. It currently assumes\n\
  the following:\n\
    - longitude and latitude are either 0-dimensional variables or 1-dimensional\n\
      variables of size 1;\n\
    - time is either 0-dimensional, 1-dimensional of size 1, or 1-dimensional of\n\
      full size (nz);\n\
    - profile variables are either 1-dimensional or 3-4 dimensional with dummy\n\
      dimensions of size 1 (yes, this has little sense, but some providers do\n\
      use such formats)\n\
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
    - TIMENAME (\"t\" | \"[tT]ime\" | \"TIME\") (+)\n\
    - or TIMENAMES (when time = base_time + offset) (+)\n\
    - LONNAME (\"lon\" | \"[lL]ongitude\" | \"LONGITUDE\") (+)\n\
    - LATNAME (\"lat\" | \"[lL]atitude\" | \"LATITUDE\") (+)\n\
    - ZNAME (\"z\" | \"[dD]epth\" | \"DEPTH\") (+)\n\
    - ESTDNAME (\"error_std\") (-)\n\
        error STD; if absent then needs to be specified in the corresponding\n\
        section of the observation data parameter file\n\
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
