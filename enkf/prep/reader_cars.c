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
#include "nan.h"
#include "stringtable.h"
#include "kdtree.h"
#include "definitions.h"
#include "utils.h"
#include "enkfprm.h"
#include "obsmeta.h"
#include "grid.h"
#include "model.h"
#include "observations.h"
#include "prep.h"
#include "allreaders.h"

#define EPS 1.0e-6

void reader_cars_standard(char* fname, int fid, obsmeta* meta, model* m, observations* obs)
{
    int ncid;
    int dimid_nprof, dimid_nz;
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

    if (meta->nstds == 0)
        enkf_quit("ERROR_STD is necessary but not specified for product \"%s\"", meta->product);

    ncw_open(fname, NC_NOWRITE, &ncid);

    ncw_inq_dimid(fname, ncid, "nobs", &dimid_nprof);
    ncw_inq_dimlen(fname, ncid, dimid_nprof, &nprof);
    enkf_printf("        # profiles = %u\n", (unsigned int) nprof);
    if (nprof == 0) {
        ncw_close(fname, ncid);
        return;
    }

    ncw_inq_dimid(fname, ncid, "zt", &dimid_nz);
    ncw_inq_dimlen(fname, ncid, dimid_nz, &nz);
    enkf_printf("        # z levels = %u\n", (unsigned int) nz);

    ncw_inq_varid(fname, ncid, "lon", &varid_lon);
    lon = malloc(nprof * sizeof(double));
    ncw_get_var_double(fname, ncid, varid_lon, lon);

    ncw_inq_varid(fname, ncid, "lat", &varid_lat);
    lat = malloc(nprof * sizeof(double));
    ncw_get_var_double(fname, ncid, varid_lat, lat);

    ncw_inq_varid(fname, ncid, "zt", &varid_z);
    z = alloc2d(nprof, nz, sizeof(double));
    ncw_get_var_double(fname, ncid, varid_z, z[0]);

    if (strncmp(meta->type, "TEM", 3) == 0)
        ncw_inq_varid(fname, ncid, "temp", &varid_v);
    else if (strncmp(meta->type, "SAL", 3) == 0)
        ncw_inq_varid(fname, ncid, "salt", &varid_v);
    else
        enkf_quit("observation type \"%s\" not handled for CARS product", meta->type);
    v = alloc2d(nprof, nz, sizeof(double));
    ncw_get_var_double(fname, ncid, varid_v, v[0]);
    ncw_get_att_double(fname, ncid, varid_v, "missing_value", &missval);
    ncw_get_att_double(fname, ncid, varid_v, "valid_min", &validmin);
    ncw_get_att_double(fname, ncid, varid_v, "valid_max", &validmax);

    ncw_inq_varid(fname, ncid, "type", &varid_type);
    type = malloc(nprof * sizeof(int));
    ncw_get_var_int(fname, ncid, varid_type, type);

    ncw_close(fname, ncid);

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
    sprintf(buf, "days since %4d-%02d-%02d", year, month, day);

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
            sprintf(inststr, "CARS%02u", type[p]);

        for (i = 0; i < (int) nz; ++i) {
            observation* o;
            obstype* ot;

            if (fabs(v[p][i] - missval) < EPS || v[p][i] < validmin || v[p][i] > validmax)
                continue;

            obs_checkalloc(obs);
            o = &obs->data[obs->nobs];

            o->product = st_findindexbystring(obs->products, meta->product);
            assert(o->product >= 0);
            o->type = obstype_getid(obs->nobstypes, obs->obstypes, meta->type);
            assert(o->type >= 0);
            ot = &obs->obstypes[o->type];
            o->instrument = st_add_ifabscent(obs->instruments, inststr, -1);
            o->id = obs->nobs;
            o->fid = fid;
            o->batch = p;
            o->value = v[p][i];
            o->std = 0.0;
            o->lon = lon[p];
            o->lat = lat[p];
            o->depth = z[p][i];
            o->status = model_xy2fij(m, o->lon, o->lat, &o->fi, &o->fj);
            if (!obs->allobs && o->status == STATUS_OUTSIDEGRID)
                break;
            if (o->status == STATUS_OK)
                o->status = model_z2fk(m, o->fi, o->fj, o->depth, &o->fk);
            else
                o->fk = NaN;
            if ((o->status == STATUS_OK) && (o->lon <= ot->xmin || o->lon >= ot->xmax || o->lat <= ot->ymin || o->lat >= ot->ymax || o->depth <= ot->zmin || o->depth >= ot->zmax))
                o->status = STATUS_OUTSIDEOBSDOMAIN;
            o->date = tunits_offset + 0.5;
            o->aux = -1;

            obs->nobs++;
        }
    }

    free(lon);
    free(lat);
    free2d(v);
    free2d(z);
    free(type);
}
