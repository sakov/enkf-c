/******************************************************************************
 *
 * File:        reader_pathfinder.c        
 *
 * Created:     10/2015
 *
 * Author:      Paul Sandery
 *              Bureau of Meteorology
 *
 * Description: contains reader reader_pathfinder_standard() for files of the
 *              form y<yyyy>/m<mm>/??_d<dd>.nc with no "time" variable. Used
 *              reader_rads_standard2() as a template.
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
#include "stringtable.h"
#include "kdtree.h"
#include "definitions.h"
#include "utils.h"
#include "obsprm.h"
#include "grid.h"
#include "observations.h"
#include "prep_utils.h"
#include "allreaders.h"

/** For files of the form y<yyyy>/m<mm>/??_d<dd>.nc with no "time" variable.
 */
void reader_pathfinder_standard(char* fname, int fid, obsmeta* meta, grid* g, observations* obs)
{
    int ksurf = grid_getsurflayerid(g);
    int ncid;
    int dimid_nobs;
    size_t nobs_local;
    int varid_lon, varid_lat, varid_error, varid_sst;
    double* lon;
    double* lat;
    double* sst;
    double* error_std;
    char buf[MAXSTRLEN];
    int len;
    int year, month, day;
    double tunits_multiple, tunits_offset;
    char* basename;
    char instname[5];
    int i;

    for (i = 0; i < meta->npars; ++i)
        enkf_quit("unknown PARAMETER \"%s\"\n", meta->pars[i].name);

    ncw_open(fname, NC_NOWRITE, &ncid);
    ncw_inq_dimid(ncid, "nobs", &dimid_nobs);
    ncw_inq_dimlen(ncid, dimid_nobs, &nobs_local);
    enkf_printf("        nobs = %u\n", (unsigned int) nobs_local);

    if (nobs_local == 0) {
        ncw_close(ncid);
        return;
    }

    ncw_inq_varid(ncid, "lon", &varid_lon);
    lon = malloc(nobs_local * sizeof(double));
    ncw_get_var_double(ncid, varid_lon, lon);

    ncw_inq_varid(ncid, "lat", &varid_lat);
    lat = malloc(nobs_local * sizeof(double));
    ncw_get_var_double(ncid, varid_lat, lat);

    ncw_inq_varid(ncid, "error", &varid_error);
    error_std = malloc(nobs_local * sizeof(double));
    ncw_get_var_double(ncid, varid_error, error_std);

    ncw_inq_varid(ncid, "sst", &varid_sst);
    sst = malloc(nobs_local * sizeof(double));
    ncw_get_var_double(ncid, varid_sst, sst);

    ncw_close(ncid);

    strcpy(buf, fname);
    len = strlen(buf);
    buf[len - 3] = 0;           /* .nc */
    if (!str2int(&buf[len - 5], &day))
        enkf_quit("PATHFINDER reader: could not convert file name \"%s\" to date", fname);
    buf[len - 5] = 0;
    if (!str2int(&buf[len - 7], &month))
        enkf_quit("PATHFINDER reader: could not convert file name \"%s\" to date", fname);
    buf[len - 7] = 0;
    if (!str2int(&buf[len - 11], &year))
        enkf_quit("PATHFINDER reader: could not convert file name \"%s\" to date", fname);
    snprintf(buf, MAXSTRLEN, "days since %4d-%02d-%02d", year, month, day);

    tunits_convert(buf, &tunits_multiple, &tunits_offset);

    basename = strrchr(fname, '/');
    if (basename == NULL)
        basename = fname;
    else
        basename += 1;
    strncpy(instname, basename, 4);
    instname[4] = 0;

    for (i = 0; i < (int) nobs_local; ++i) {
        observation* o;

        obs_checkalloc(obs);
        o = &obs->data[obs->nobs];

        o->product = st_findindexbystring(obs->products, meta->product);
        assert(o->product >= 0);
        o->type = obstype_getid(obs->nobstypes, obs->obstypes, meta->type, 1);
        o->instrument = st_add_ifabsent(obs->instruments, instname, -1);
        o->id = obs->nobs;
        o->fid = fid;
        o->batch = 0;
        o->value = sst[i];
        o->estd = error_std[i];
        o->lon = lon[i];
        o->lat = lat[i];
        o->depth = 0.0;
        o->fk = (double) ksurf;
        o->status = grid_xy2fij(g, o->lon, o->lat, &o->fi, &o->fj);
        if (!obs->allobs && o->status == STATUS_OUTSIDEGRID)
            continue;
        o->model_depth = NAN;   /* set in obs_add() */
        o->day = tunits_offset + 0.5;
        o->aux = -1;

        obs->nobs++;
    }

    free(lon);
    free(lat);
    free(sst);
    free(error_std);
}
