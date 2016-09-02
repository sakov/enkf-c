/******************************************************************************
 *
 * File:        reader_mmt.c        
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
#include <float.h>
#include "ncw.h"
#include "nan.h"
#include "definitions.h"
#include "utils.h"
#include "obsmeta.h"
#include "grid.h"
#include "model.h"
#include "observations.h"
#include "prep_utils.h"

#define EPS 1.0e-6
#define WMO_INSTSIZE 4

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

void reader_mmt_standard(char* fname, int fid, obsmeta* meta, model* m, observations* obs)
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
    double missval;
    double validmin = DBL_MAX;
    double validmax = -DBL_MAX;
    char* type;
    char buf[MAXSTRLEN];
    int len;
    int year, month, day;
    double tunits_multiple, tunits_offset;
    int mvid;
    float** depth;
    int p, i;

    for (i = 0; i < meta->npars; ++i)
        enkf_quit("unknown PARAMETER \"%s\"\n", meta->pars[i].name);

    if (meta->nstds == 0)
        enkf_quit("ERROR_STD is necessary but not specified for product \"%s\"", meta->product);

    ncw_open(fname, NC_NOWRITE, &ncid);

    ncw_inq_dimid(ncid, "N_PROF", &dimid_nprof);
    ncw_inq_dimlen(ncid, dimid_nprof, &nprof);

    ncw_inq_dimid(ncid, "N_LEVELS", &dimid_nz);
    ncw_inq_dimlen(ncid, dimid_nz, &nz);
    enkf_printf("        # profiles = %u\n", (unsigned int) nprof);
    if (nprof == 0) {
        ncw_close(ncid);
        return;
    }
    enkf_printf("        # z levels = %u\n", (unsigned int) nz);

    ncw_inq_varid(ncid, "LONGITUDE", &varid_lon);
    lon = malloc(nprof * sizeof(double));
    ncw_get_var_double(ncid, varid_lon, lon);

    ncw_inq_varid(ncid, "LATITUDE", &varid_lat);
    lat = malloc(nprof * sizeof(double));
    ncw_get_var_double(ncid, varid_lat, lat);

    ncw_inq_varid(ncid, "PRES_BLUELINK", &varid_z);
    z = alloc2d(nprof, nz, sizeof(double));
    ncw_get_var_double(ncid, varid_z, z[0]);

    if (strncmp(meta->type, "TEM", 3) == 0) {
        validmin = -2.0;
        validmax = 40.0;
        ncw_inq_varid(ncid, "TEMP_BLUELINK", &varid_v);
    } else if (strncmp(meta->type, "SAL", 3) == 0) {
        validmin = 0;
        validmax = 50.0;
        ncw_inq_varid(ncid, "PSAL_BLUELINK", &varid_v);
    } else
        enkf_quit("observation type \"%s\" not handled for MMT product", meta->type);
    v = alloc2d(nprof, nz, sizeof(double));
    ncw_get_var_double(ncid, varid_v, v[0]);
    ncw_get_att_double(ncid, varid_v, "_FillValue", &missval);

    ncw_inq_varid(ncid, "WMO_INST_TYPE", &varid_type);
    type = malloc(nprof * WMO_INSTSIZE);
    ncw_get_var_text(ncid, varid_type, type);

    ncw_close(ncid);

    strcpy(buf, fname);
    len = strlen(buf);
    buf[len - 10] = 0;          /* _mmt_qc.nc */
    if (!str2int(&buf[len - 12], &day))
        enkf_quit("MMT reader: could not convert file name \"%s\" to date", fname);
    buf[len - 12] = 0;
    if (!str2int(&buf[len - 14], &month))
        enkf_quit("MMT reader: could not convert file name \"%s\" to date", fname);
    buf[len - 14] = 0;
    if (!str2int(&buf[len - 18], &year))
        enkf_quit("MMT reader: could not convert file name \"%s\" to date", fname);
    snprintf(buf, MAXSTRLEN, "days since %4d-%02d-%02d", year, month, day);

    tunits_convert(buf, &tunits_multiple, &tunits_offset);

    mvid = model_getvarid(m, obs->obstypes[obstype_getid(obs->nobstypes, obs->obstypes, meta->type)].varnames[0], 1);
    depth = model_getdepth(m, mvid, 0);

    for (p = 0; p < (int) nprof; ++p) {
        char inststr[MAXSTRLEN];

        snprintf(inststr, MAXSTRLEN, "WMO%04u", type[p * WMO_INSTSIZE]);

        for (i = 0; i < (int) nz; ++i) {
            observation* o;
            obstype* ot;

            if (fabs(v[p][i] - missval) < EPS || v[p][i] < validmin || v[p][i] > validmax)
                continue;
            if (z[p][i] < 0.0)
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
            o->status = model_xy2fij(m, mvid, o->lon, o->lat, &o->fi, &o->fj);
            if (!obs->allobs && o->status == STATUS_OUTSIDEGRID)
                break;
            if (o->status == STATUS_OK)
                o->status = model_z2fk(m, mvid, o->fi, o->fj, o->depth, &o->fk);
            else
                o->fk = NaN;
            if ((o->status == STATUS_OK) && (o->lon <= ot->xmin || o->lon >= ot->xmax || o->lat <= ot->ymin || o->lat >= ot->ymax || o->depth <= ot->zmin || o->depth >= ot->zmax))
                o->status = STATUS_OUTSIDEOBSDOMAIN;
            o->model_depth = (depth == NULL || isnan(o->fi + o->fj)) ? NaN : depth[(int) (o->fj + 0.5)][(int) (o->fi + 0.5)];
            o->date = tunits_offset + 0.5;
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
        int ii = 0;

        for (i = 0; i < nprof; ++i) {
            lonlat[i * 2] = lon[i];
            lonlat[i + 2 + 1] = lat[i];
        }
        qsort(lonlat, nprof, sizeof(double) * 2, cmp_lonlat);
        for (i = 1; i < nprof; ++i) {
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
