/******************************************************************************
 *
 * File:        reader_gridded_xyh.c
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
 *              Unlike other readers this one requires the (hybrid) data grid
 *              to be defined in the grid parameter file.
 *
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
void reader_gridded_xyh(char* fname, int fid, obssection* section, grid* gdst, observations* obs)
{
    char* varname = NULL;
    char* gridname = NULL;
    grid* gsrc = NULL;

    char* npointsname = NULL;
    obsatts oa;

    int ni = 0, nj = 0, nk = 0;
    size_t nij = 0, nijk = 0;

    int productid = -1;
    int typeid = -1;

    int ncid;
    int ndim;
    float* var = NULL;
    short* npoints = NULL;
    size_t ntime = 0;
    double* time = NULL;
    int varid;
    size_t i, j, k, nobs_read;

    obsatts_init(&oa);

    for (i = 0; i < (size_t) section->npars; ++i) {
        if (strcasecmp(section->pars[i].name, "VARNAME") == 0)
            varname = section->pars[i].value;
        else if (strcasecmp(section->pars[i].name, "GRIDNAME") == 0)
            gridname = section->pars[i].value;
        else if (strcasecmp(section->pars[i].name, "NPOINTSNAME") == 0)
            npointsname = section->pars[i].value;
        else if (!obsatts_parse_param(&oa, section->pars[i].name, section->pars[i].value))
            enkf_quit("unknown PARAMETER \"%s\"\n", section->pars[i].name);
    }

    if (varname == NULL)
        enkf_quit("reader_xyh_gridded(): %s: VARNAME not specified", fname);
    else
        enkf_printf("        VARNAME = %s\n", varname);
    if (gridname == NULL)
        enkf_quit("reader_xyh_gridded(): %s: GRIDNAME not specified", fname);
    else
        enkf_printf("        GRIDNAME = %s\n", gridname);

    gsrc = model_getgridbyname(obs->model, gridname);

    grid_getsize(gsrc, &ni, &nj, &nk);
    nij = (size_t) ni* (size_t) nj;

    nijk = nij * (size_t) nk;

    ncw_open(fname, NC_NOWRITE, &ncid);

    ncw_inq_varid(ncid, varname, &varid);
    ncw_inq_varndims(ncid, varid, &ndim);
    if (ndim != 3 && ndim != 4)
        enkf_quit("reader_xyh_gridded(): %s: # dimensions = %d (must be 3 or 4 with a single record)", fname, ndim);
    {
        int dimid[4];
        size_t dimlen[4];

        ncw_inq_vardimid(ncid, varid, dimid);
        for (i = 0; i < (size_t) ndim; ++i)
            ncw_inq_dimlen(ncid, dimid[i], &dimlen[i]);
        if (ndim == 4 && dimlen[0] != 1)
            enkf_quit("reader_xyh_gridded(): %s: %zu records (currently only one record is allowed)", fname, dimlen[0]);

        if (dimlen[ndim - 1] != (size_t) ni || dimlen[ndim - 2] != (size_t) nj || dimlen[ndim - 3] != (size_t) nk)
            enkf_quit("dimension mismatch between grid \"%s\" (%d x %d x %d) and variable \"%s\" in \"%s\" (%zu x %zu x %zu)", gridname, ni, nj, nk, varname, fname, dimlen[ndim - 1], dimlen[ndim - 2], dimlen[ndim - 3]);
    }

    var = malloc(nijk * sizeof(float));
    ncu_readvarfloat(ncid, varid, nijk, var);

    /*
     * npoints
     */
    varid = -1;
    if (npointsname != NULL)
        ncw_inq_varid(ncid, npointsname, &varid);
    else if (ncw_var_exists(ncid, "npoints"))
        ncw_inq_varid(ncid, "npoints", &varid);
    if (varid >= 0) {
        ncw_check_varsize(ncid, varid, nijk);
        npoints = malloc(nijk * sizeof(short));
        ncw_get_var_short(ncid, varid, npoints);
    }

    /*
     * std, estd, batch, qcflag, instrument
     */
    obsatts_read(&oa, section, ncid, varname, nijk, nijk, 1, 0, 0, obs);

    /*
     * time
     */
    get_time(section, ncid, &ntime, &time);
    assert(ntime == nijk || ntime <= 1);

    ncw_close(ncid);

    productid = st_findindexbystring(obs->products, section->product);
    assert(productid >= 0);
    typeid = obstype_getid(obs->nobstypes, obs->obstypes, section->type, 1);
    nobs_read = 0;
    for (i = 0; i < (size_t) ni; ++i) {
        for (j = 0; j < (size_t) nj; ++j) {
            int ij[2] = { (int) i, (int) j };

            for (k = 0; k < (size_t) nk; ++k) {
                size_t ii = k * nij + j * (size_t) ni + i;
                observation* o;

                if ((npoints != NULL && npoints[ii] == 0) || isnan(var[ii]) || (ntime == nijk && isnan(time[ii])))
                    continue;
                if (!obsatts_valid(&oa, ii))
                    continue;

                nobs_read++;
                obs_checkalloc(obs);
                o = &obs->data[obs->nobs];

                o->product = productid;
                o->type = typeid;
                o->id = obs->nobs;
                o->id_orig = (int) ii;
                o->fid = fid;
                o->value = (double) var[ii];
                obsatts_apply(&oa, ii, 0, o);
                {
                    double lon_d, lat_d;

                    grid_ij2xy(gsrc, ij, &lon_d, &lat_d);
                    o->lon = (float) lon_d;
                    o->lat = (float) lat_d;
                }
                assert(isfinite(o->lon + o->lat));
                o->status = grid_xy2fij(gdst, o->lon, o->lat, o->fij);
                if (!obs->allobs && o->status == STATUS_OUTSIDEGRID)
                    continue;
                {
                    double depth_d;
                    int status = grid_fk2z(gsrc, ij, (double) k, &depth_d);

                    o->depth = (float) depth_d;
                    if (o->status == STATUS_OK)
                        o->status = status;
                }
                if (o->status == STATUS_OK)
                    o->status = grid_z2fk_f(gdst, o->fij, o->depth, &o->fk);
                else
                    o->fk = NAN;
                o->model_depth = NAN;   /* set in obs_add() */
                if (ntime > 0)
                    o->time = (ntime == 1) ? time[0] : time[ii];
                else
                    o->time = NAN;

                o->aux = -1;

                obs->nobs++;
            }
        }
    }
    enkf_printf("        nobs = %zu\n", nobs_read);

    free(var);
    if (npoints != NULL)
        free(npoints);
    if (time != NULL)
        free(time);
    obsatts_free(&oa);
}

/**
 */
void reader_gridded_xyh_describe(void)
{
    enkf_printf("\n  Generic reader \"gridded_xyh\" reads 3D gridded data with hybrid\n\
  vertical coordinate. It assumes\n\
\n\
    p[k] = a[k] + b[k] * (p1 - p2),\n\
\n\
  where p(k) is z coordinate (pressure) at the middle of layer k, p1 - surface\n\
  pressure, and p2 - top pressure, see\n\
  https://ametsoc.org for details. The\n\
  grid must be defined in the grid parameter file.\n");
    describe_commongenericreaderparams();
    enkf_printf("  Parameters specific to the reader:\n\
    - GRIDNAME (++)\n\
    - NPOINTSNAME (\"npoints\") (-)\n\
        number of collated points for each datum; used basically as a data mask\n\
        when n = 0\n");
    describe_commonreaderparams();
}
