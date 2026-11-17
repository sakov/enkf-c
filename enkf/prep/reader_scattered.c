/******************************************************************************
 *
 * File:        reader_scattered.c        
 *
 * Created:     21/03/2017
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description: Generic reader for scattered 2D or 3D observations.
 *
 * Revisions:   PS 4/7/2018
 *                Added parameters QCFLAGNAME and QCFLAGVALS. The latter is
 *                supposed to contain a list of allowed flag values.
 *              PS 16/09/2020
 *                This reader is creatred by merging reader_xy_scattered()
 *                and reader_xyz_scattered(). The functionality of 
 *                reader_xy_scattered() is achieved by specifying parameter
 *                ZVALUE (typically set to 0).
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
#include "global.h"
#include "ncutils.h"
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
void reader_scattered(char* fname, int fid, obssection* section, grid* g, observations* obs)
{
    char* varname = NULL;
    char* lonname = NULL;
    char* latname = NULL;
    char* zname = NULL;
    double zvalue = NAN;
    int zvalueentered = 0;
    obsatts oa;

    int productid = -1;
    int typeid = -1;

    int ncid;
    size_t nobs;
    double* lon = NULL;
    double* lat = NULL;
    double* z = NULL;
    double* var = NULL;
    size_t ntime = 0;
    double* time = NULL;
    int varid;
    size_t i, nobs_read;

    obsatts_init(&oa);

    for (i = 0; i < section->npars; ++i) {
        if (strcasecmp(section->pars[i].name, "VARNAME") == 0)
            varname = section->pars[i].value;
        else if (strcasecmp(section->pars[i].name, "LONNAME") == 0)
            lonname = section->pars[i].value;
        else if (strcasecmp(section->pars[i].name, "LATNAME") == 0)
            latname = section->pars[i].value;
        else if (strcasecmp(section->pars[i].name, "ZNAME") == 0) {
            if (zvalueentered)
                enkf_quit("reader_scattered(): can not simultaneously specify ZNAME and ZVALUE");
            zname = section->pars[i].value;
        } else if (strcasecmp(section->pars[i].name, "ZVALUE") == 0) {
            if (zname != NULL)
                enkf_quit("reader_scattered(): can not simultaneously specify ZNAME and ZVALUE");
            if (!str2double(section->pars[i].value, &zvalue))
                enkf_quit("observation prm file: can not convert ZVALUE = \"%s\" to double\n", section->pars[i].value);
            zvalueentered = 1;
        } else if (!obsatts_parse_param(&oa, section->pars[i].name, section->pars[i].value))
            enkf_quit("reader_scattered(): unknown PARAMETER \"%s\"\n", section->pars[i].name);
    }

    if (varname == NULL)
        enkf_quit("reader_scattered(): %s VARNAME not specified", fname);
    else
        enkf_printf("        VARNAME = %s\n", varname);

    ncw_open(fname, NC_NOWRITE, &ncid);

    /*
     * main variable
     */
    ncw_inq_varid(ncid, varname, &varid);
    ncw_inq_vardims(ncid, varid, 1, NULL, &nobs);
    var = malloc(nobs * sizeof(double));
    ncu_readvardouble(ncid, varid, nobs, var);

    /*
     * longitude
     */
    lonname = get_lonname(ncid, lonname);
    if (lonname != NULL) {
        enkf_printf("        LONNAME = %s\n", lonname);
        ncw_inq_varid(ncid, lonname, &varid);
        ncw_check_vardims(ncid, varid, 1, &nobs);
    } else
        enkf_quit("reader_scattered(): %s: could not find longitude variable", fname);
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
        enkf_quit("reader_scattered(): %s: could not find latitude variable", fname);
    lat = malloc(nobs * sizeof(double));
    ncu_readvardouble(ncid, varid, nobs, lat);

    /*
     * z
     */
    if (!zvalueentered) {
        zname = get_zname(ncid, zname);
        if (zname != NULL) {
            enkf_printf("        ZNAME = %s\n", zname);
            ncw_inq_varid(ncid, zname, &varid);
            ncw_check_vardims(ncid, varid, 1, &nobs);
        } else
            enkf_quit("reader_scattered(): %s: could not find Z variable; and Z value is not specified", fname);
        z = malloc(nobs * sizeof(double));
        ncu_readvardouble(ncid, varid, nobs, z);
    }

    /*
     * std, estd, batch, qcflag, instrument
     */
    obsatts_read(&oa, section, ncid, varname, nobs, nobs, 1, 1, 0, obs);

    /*
     * time
     */
    get_time(section, ncid, &ntime, &time);
    assert(ntime == nobs || ntime <= 1);

    ncw_close(ncid);

    productid = st_findindexbystring(obs->products, section->product);
    assert(productid >= 0);
    typeid = obstype_getid(obs->nobstypes, obs->obstypes, section->type, 1);
    nobs_read = 0;
    for (i = 0; i < nobs; ++i) {
        observation* o;

        if (isnan(lon[i]) || isnan(lat[i]) || (z != NULL && isnan(z[i])) || isnan(var[i]))
            continue;
        if ((ntime == nobs && isnan(time[i])) || !obsatts_valid(&oa, i))
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
        o->lon = lon[i];
        o->lat = lat[i];
        o->depth = (z != NULL) ? z[i] : zvalue;
        o->status = grid_xy2fij(g, o->lon, o->lat, o->fij);
        if (!obs->allobs && o->status == STATUS_OUTSIDEGRID)
            continue;
        if (o->status == STATUS_OK && isfinite(o->depth))
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
    }
    enkf_printf("        nobs = %zu\n", nobs_read);

    free(lon);
    free(lat);
    if (z != NULL)
        free(z);
    free(var);
    if (time != NULL)
        free(time);
    obsatts_free(&oa);
}

/**
 */
void reader_scattered_describe(void)
{
    enkf_printf("\n  Generic reader \"scattered\" reads 2D or 3D scattered point data.\n");
    describe_commongenericreaderparams();
    describe_commonreaderparams();
}
