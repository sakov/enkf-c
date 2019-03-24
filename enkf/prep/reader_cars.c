/******************************************************************************
 *
 * File:        reader_cars.c        
 *
 * Created:     12/2012
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description:
 *
 * Revisions:
 *
 *****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "ncw.h"
#include "definitions.h"
#include "utils.h"
#include "obsprm.h"
#include "grid.h"
#include "model.h"
#include "observations.h"
#include "prep_utils.h"
#include "allreaders.h"

#define EPS 1.0e-6

/**
 */
static int cmp_lonlat(const void* p1, const void* p2)
{
    double* ll1 = (double*) p1;
    double* ll2 = (double*) p2;

    if (ll1[0] > ll2[0])
        return 1;
    else if (ll1[0] < ll2[0])
        return -1;
    else if (ll1[1] > ll2[1])
        return 1;
    else if (ll1[1] < ll2[1])
        return -1;
    return 0;
}

/**
 */
void reader_cars_standard(char* fname, int fid, obsmeta* meta, grid* g, observations* obs)
{
    int ncid;
    int dimid_nprof, dimid_nz = -1;
    size_t nprof, nz;
    int varid_lon, varid_lat, varid_z, varid_type;
    int varid_v = -1;
    double* lon;
    double* lat;
    double** z;
    double** v;
    double missval, validmin, validmax;
    int* type;
    char buf[MAXSTRLEN];
    int len;
    int year, month, day;
    double tunits_multiple, tunits_offset;
    int p, i;

    for (i = 0; i < meta->npars; ++i)
        enkf_quit("unknown PARAMETER \"%s\"\n", meta->pars[i].name);

    if (meta->nestds == 0)
        enkf_quit("ERROR_STD is necessary but not specified for product \"%s\"", meta->product);

    ncw_open(fname, NC_NOWRITE, &ncid);
    ncw_inq_dimid(ncid, "nobs", &dimid_nprof);
    ncw_inq_dimlen(ncid, dimid_nprof, &nprof);
    enkf_printf("        # profiles = %u\n", (unsigned int) nprof);
    if (nprof == 0) {
        ncw_close(ncid);
        return;
    }

    if (ncw_dim_exists(ncid, "zt"))
        ncw_inq_dimid(ncid, "zt", &dimid_nz);
    else if (ncw_dim_exists(ncid, "ztd"))
        ncw_inq_dimid(ncid, "ztd", &dimid_nz);
    else
        enkf_quit("reader_cars_standard(): neither dimension \"zt\" ot \"ztd\" exist");
    ncw_inq_dimlen(ncid, dimid_nz, &nz);
    enkf_printf("        # z levels = %u\n", (unsigned int) nz);

    ncw_inq_varid(ncid, "lon", &varid_lon);
    lon = malloc(nprof * sizeof(double));
    ncw_get_var_double(ncid, varid_lon, lon);

    ncw_inq_varid(ncid, "lat", &varid_lat);
    lat = malloc(nprof * sizeof(double));
    ncw_get_var_double(ncid, varid_lat, lat);

    ncw_inq_varid(ncid, "zt", &varid_z);
    z = alloc2d(nprof, nz, sizeof(double));
    ncw_get_var_double(ncid, varid_z, z[0]);

    if (strncmp(meta->type, "TEM", 3) == 0)
        ncw_inq_varid(ncid, "temp", &varid_v);
    else if (strncmp(meta->type, "SAL", 3) == 0)
        ncw_inq_varid(ncid, "salt", &varid_v);
    else
        enkf_quit("observation type \"%s\" not handled for CARS product", meta->type);
    v = alloc2d(nprof, nz, sizeof(double));
    ncw_get_var_double(ncid, varid_v, v[0]);
    ncw_get_att_double(ncid, varid_v, "missing_value", &missval);
    ncw_get_att_double(ncid, varid_v, "valid_min", &validmin);
    ncw_get_att_double(ncid, varid_v, "valid_max", &validmax);

    ncw_inq_varid(ncid, "type", &varid_type);
    type = malloc(nprof * sizeof(int));
    ncw_get_var_int(ncid, varid_type, type);

    ncw_close(ncid);

    strcpy(buf, fname);
    len = strlen(buf);
    buf[len - 3] = 0;           /* .nc */
    if (!str2int(&buf[len - 5], &day))
        enkf_quit("CARS reader: could not convert file name \"%s\" to date", fname);
    buf[len - 17] = 0;
    if (!str2int(&buf[len - 19], &month))
        enkf_quit("CARS reader: could not convert file name \"%s\" to date", fname);
    buf[len - 21] = 0;
    if (!str2int(&buf[len - 25], &year))
        enkf_quit("CARS reader: could not convert file name \"%s\" to date", fname);
    snprintf(buf, MAXSTRLEN, "days since %4d-%02d-%02d", year, month, day);

    tunits_convert(buf, &tunits_multiple, &tunits_offset);

    for (p = 0; p < (int) nprof; ++p) {
        char inststr[MAXSTRLEN];

        if (type[p] == 11)
            strcpy(inststr, "ARGO");
        else if (type[p] == 12)
            strcpy(inststr, "TAO");
        else if (type[p] == 61)
            strcpy(inststr, "PIRATA");
        else if (type[p] == 7 || type[p] == 9 || type[p] == 13 || type[p] == 35 || type[p] == 41)
            strcpy(inststr, "CTD");
        else if (type[p] == 8 || type[p] == 17)
            strcpy(inststr, "XBT");
        else
            snprintf(inststr, MAXSTRLEN, "CARS%02u", type[p]);

        for (i = 0; i < (int) nz; ++i) {
            observation* o;

            if (fabs(v[p][i] - missval) < EPS || v[p][i] < validmin || v[p][i] > validmax)
                continue;
            if (z[p][i] < 0.0)
                continue;

            obs_checkalloc(obs);
            o = &obs->data[obs->nobs];

            o->product = st_findindexbystring(obs->products, meta->product);
            assert(o->product >= 0);
            o->type = obstype_getid(obs->nobstypes, obs->obstypes, meta->type, 1);
            o->instrument = st_add_ifabsent(obs->instruments, inststr, -1);
            o->id = obs->nobs;
            o->fid = fid;
            o->batch = p;
            o->value = v[p][i];
            o->estd = 0.0;
            o->lon = lon[p];
            o->lat = lat[p];
            o->depth = z[p][i];
            o->status = grid_xy2fij(g, o->lon, o->lat, &o->fi, &o->fj);
            if (!obs->allobs && o->status == STATUS_OUTSIDEGRID)
                break;
            if (o->status == STATUS_OK)
                o->status = grid_z2fk(g, o->fi, o->fj, o->depth, &o->fk);
            else
                o->fk = NAN;
            o->model_depth = NAN;       /* set in obs_add() */
            o->day = tunits_offset + 0.5;
            o->aux = -1;

            obs->nobs++;
        }
    }

    /*
     * get the number of unique profile locations
     */
    {
        double* lonlat = malloc(nprof * sizeof(double) * 2);
        int nunique = (nprof > 0) ? 1 : 0;
        int ii;

        for (i = 0; i < nprof; ++i) {
            lonlat[i * 2] = lon[i];
            lonlat[i * 2 + 1] = lat[i];
        }
        qsort(lonlat, nprof, sizeof(double) * 2, cmp_lonlat);
        for (i = 1, ii = 0; i < nprof; ++i) {
            if (lonlat[i * 2] == lonlat[ii * 2] && lonlat[i * 2 + 1] == lonlat[ii * 2 + 1])
                continue;
            ii = i;
            nunique++;
        }
        enkf_printf("        # unique locations = %d\n", nunique);
        free(lonlat);
    }

    free(lon);
    free(lat);
    free(v);
    free(z);
    free(type);
}
