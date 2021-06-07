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
#include "definitions.h"
#include "utils.h"
#include "obsprm.h"
#include "grid.h"
#include "model.h"
#include "observations.h"
#include "prep_utils.h"
#include "allreaders.h"

#define EPS 1.0e-6
#define WMO_INSTSIZE 4

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
void reader_mmt(char* fname, int fid, obsmeta* meta, grid* g, observations* obs)
{
    stringtable* st_exclude = NULL;
    int ncid;
    int dimid_nprof, dimid_nz;
    size_t nprof, nz;
    int varid;
    int* status = NULL;
    double* lon;
    double* lat;
    double** z;
    double** v;
    double* time;
    char tunits[MAXSTRLEN];
    int** qc = NULL;
    double missval;
    double validmin = DBL_MAX;
    double validmax = -DBL_MAX;
    char* type;
    double tunits_multiple, tunits_offset;
    int p, i, nobs_read;

    for (i = 0; i < meta->npars; ++i) {
        if (strcasecmp(meta->pars[i].name, "EXCLUDEINST") == 0) {
            if (st_exclude == NULL)
                st_exclude = st_create("exclude");
            st_add(st_exclude, meta->pars[i].value, -1);
        } else
            enkf_quit("unknown PARAMETER \"%s\"\n", meta->pars[i].name);
    }

    if (meta->nestds == 0)
        enkf_quit("ERROR_STD is necessary but not specified for product \"%s\"", meta->product);

    ncw_open(fname, NC_NOWRITE, &ncid);

    ncw_inq_dimid(ncid, "N_PROF", &dimid_nprof);
    ncw_inq_dimlen(ncid, dimid_nprof, &nprof);

    ncw_inq_dimid(ncid, "N_LEVELS", &dimid_nz);
    ncw_inq_dimlen(ncid, dimid_nz, &nz);
    enkf_printf("        # profiles = %u\n", (unsigned int) nprof);
    if (nprof == 0) {
        enkf_printf("        no profiles found\n");
        ncw_close(ncid);
        goto noprofiles;
    }
    status = calloc(nprof, sizeof(int));

    enkf_printf("        # z levels = %u\n", (unsigned int) nz);

    ncw_inq_varid(ncid, "LONGITUDE", &varid);
    lon = malloc(nprof * sizeof(double));
    ncw_get_var_double(ncid, varid, lon);

    ncw_inq_varid(ncid, "LATITUDE", &varid);
    lat = malloc(nprof * sizeof(double));
    ncw_get_var_double(ncid, varid, lat);

    ncw_inq_varid(ncid, "PRES_BLUELINK", &varid);
    z = alloc2d(nprof, nz, sizeof(double));
    ncw_get_var_double(ncid, varid, z[0]);

    if (strncmp(meta->type, "TEM", 3) == 0) {
        validmin = -2.0;
        validmax = 40.0;
        ncw_inq_varid(ncid, "TEMP_BLUELINK", &varid);
    } else if (strncmp(meta->type, "SAL", 3) == 0) {
        validmin = 0;
        validmax = 50.0;
        ncw_inq_varid(ncid, "PSAL_BLUELINK", &varid);
    } else
        enkf_quit("observation type \"%s\" not handled for MMT product", meta->type);
    v = alloc2d(nprof, nz, sizeof(double));
    ncw_get_var_double(ncid, varid, v[0]);
    ncw_get_att_double(ncid, varid, "_FillValue", &missval);

    ncw_inq_varid(ncid, "JULD", &varid);
    time = malloc(nprof * sizeof(double));
    ncw_get_var_double(ncid, varid, time);
    ncw_get_att_text(ncid, varid, "units", tunits);
    tunits_convert(tunits, &tunits_multiple, &tunits_offset);

    varid = -1;
    if (strncmp(meta->type, "TEM", 3) == 0) {
        if (ncw_var_exists(ncid, "TEMP_BLUELINK_QC"))
            ncw_inq_varid(ncid, "TEMP_BLUELINK_QC", &varid);
        else
            enkf_printf("        no \"TEMP_BLUELINK_QC\" flag\n");
    } else if (strncmp(meta->type, "SAL", 3) == 0) {
        if (ncw_var_exists(ncid, "PSAL_BLUELINK_QC"))
            ncw_inq_varid(ncid, "PSAL_BLUELINK_QC", &varid);
        else
            enkf_printf("        no \"PSAL_BLUELINK_QC\" flag\n");
    }
    if (varid >= 0) {
        qc = alloc2d(nprof, nz, sizeof(int));
        ncw_get_var_int(ncid, varid, qc[0]);
    }

    ncw_inq_varid(ncid, "WMO_INST_TYPE", &varid);
    type = malloc(nprof * WMO_INSTSIZE);
    ncw_get_var_text(ncid, varid, type);

    ncw_close(ncid);

    nobs_read = 0;
    for (p = 0; p < (int) nprof; ++p) {
        char inststr[MAXSTRLEN];
        int instnum;

        strncpy(inststr, &type[p * WMO_INSTSIZE], 4);
        inststr[4] = 0;
        if (!str2int(inststr, &instnum))
            instnum = 999;
        snprintf(inststr, MAXSTRLEN, "WMO%03d", instnum);

        if (st_exclude != NULL && st_findindexbystring(st_exclude, inststr) >= 0)
            continue;

        for (i = 0; i < (int) nz; ++i) {
            observation* o;

            if (fabs(v[p][i] - missval) < EPS || v[p][i] < validmin || v[p][i] > validmax)
                continue;
            if (z[p][i] < 0.0)
                continue;
            if (qc != NULL && qc[p][i] != 0)
                continue;

            nobs_read++;
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
            o->status = grid_xy2fij_f(g, o->lon, o->lat, &o->fi, &o->fj);
            if (!obs->allobs && o->status == STATUS_OUTSIDEGRID)
                break;
            if (o->status == STATUS_OK)
                o->status = grid_z2fk_f(g, o->fi, o->fj, o->depth, &o->fk);
            else
                o->fk = NAN;
            o->model_depth = NAN;       /* set in obs_add() */
            o->time = time[p] * tunits_multiple + tunits_offset;
            o->aux = -1;
            if (o->status == STATUS_OK)
                status[p] = 1;

            obs->nobs++;
        }
    }
    enkf_printf("        nobs = %d\n", nobs_read);

    /*
     * get the number of unique profile locations
     */
    {
        double* lonlat = malloc(nprof * sizeof(double) * 2);
        int ngood, ii;

        ngood = 0;
        for (i = 0; i < nprof; ++i) {
            if (status[i] == 0)
                continue;
            lonlat[ngood * 2] = lon[i];
            lonlat[ngood * 2 + 1] = lat[i];
            ngood++;
        }
        enkf_printf("        # profiles with data to process = %d\n", ngood);
        if (ngood > 1) {
            int nunique = 1;

            qsort(lonlat, ngood, sizeof(double) * 2, cmp_lonlat);
            for (i = 1, ii = 0; i < ngood; ++i) {
                if (lonlat[i * 2] == lonlat[ii * 2] && lonlat[i * 2 + 1] == lonlat[ii * 2 + 1])
                    continue;
                ii = i;
                nunique++;
            }
            enkf_printf("        # unique locations = %d\n", nunique);
        }
        free(lonlat);
    }

    free(lon);
    free(lat);
    free(status);
    free(v);
    free(z);
    if (qc != NULL)
        free(qc);
    free(type);
  noprofiles:
    if (st_exclude != NULL)
        st_destroy(st_exclude);
}

/**
 */
void reader_mmt_describe(void)
{
}
