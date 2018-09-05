/******************************************************************************
 *
 * File:        reader_cmems.c        
 *
 * Created:     5/9/2018
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description: Reader for profile data from Copernicus Marine Environment
 *              Monitoring Service. Based on the description available from
 *              http://cmems-resources.cls.fr/documents/PUM/
 *              CMEMS-INS-PUM-013-001-b.pdf.
 *              
 *              The optional parameters:
 *              - QCFLAGVALS
 *                  the list of allowed values of QC flag variables, from the
 *                  interval [0,31]; by default PARAMETER QCFLAGVALS = 1
 *              - EXCLUDEINST
 *                  WMO instrument type id (p. 27 of the above manual) to be
 *                  skipped, one entry (one line) per instrument, e.g.
 *                  PARAMETER EXCLUDEINST = 995 # marine mammals
 *                  PARAMETER EXCLUDEINST = 999 # unknown
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
#define QCFLAGVALS_DEF 2        /* 0b00000000000000000000000000000010
                                 * corresponds to QCFLAG = 1 */

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
void reader_cmems_standard(char* fname, int fid, obsmeta* meta, grid* g, observations* obs)
{
    stringtable* st_exclude = NULL;
    int ncid;
    int dimid_nprof, dimid_nz;
    size_t nprof, nz;
    int varid, varid_qc = -1;
    double* lon;
    double* lat;
    double** z;
    double** v;
    char** qcflag;
    uint32_t qcflagvals = QCFLAGVALS_DEF;
    double* time;
    double missval;
    double validmin = DBL_MAX;
    double validmax = -DBL_MAX;
    char* insttype;
    char buf[MAXSTRLEN];
    double tunits_multiple, tunits_offset;
    int p, i;

    for (i = 0; i < meta->npars; ++i) {
        if (strcasecmp(meta->pars[i].name, "EXCLUDEINST") == 0) {
            if (st_exclude == NULL)
                st_exclude = st_create("exclude");
            st_add(st_exclude, meta->pars[i].value, -1);
        } else if (strcasecmp(meta->pars[i].name, "QCFLAGVALS") == 0) {
            char seps[] = " ,";
            char* line = strdup(meta->pars[i].value);
            char* token;
            int val;

            assert(meta->pars[i].value != NULL);        /* (supposed to be
                                                         * impossible) */
            qcflagvals = 0;
            while ((token = strtok(line, seps)) != NULL) {
                if (!str2int(token, &val))
                    enkf_quit("%s: could not convert QCFLAGVALS entry \"%s\" to integer", meta->prmfname, token);
                if (val < 0 || val > 31)
                    enkf_quit("%s: QCFLAGVALS entry = %d (supposed to be in [0,31] interval", meta->prmfname, val);
                qcflagvals |= 1 << val;
                line = NULL;
            }
            free(line);
            if (qcflagvals == 0)
                enkf_quit("%s: no valid flag entries found after QCFLAGVALS\n", meta->prmfname);
        } else
            enkf_quit("unknown PARAMETER \"%s\"\n", meta->pars[i].name);
    }

    if (meta->nstds == 0)
        enkf_quit("ERROR_STD is necessary but not specified for product \"%s\"", meta->product);

    ncw_open(fname, NC_NOWRITE, &ncid);

    ncw_inq_dimid(ncid, "N_PROF", &dimid_nprof);
    ncw_inq_dimlen(ncid, dimid_nprof, &nprof);
    enkf_printf("        # profiles = %u\n", (unsigned int) nprof);

    ncw_inq_dimid(ncid, "N_LEVELS", &dimid_nz);
    ncw_inq_dimlen(ncid, dimid_nz, &nz);
    if (nprof == 0) {
        ncw_close(ncid);
        return;
    }
    enkf_printf("        # levels = %u\n", (unsigned int) nz);

    ncw_inq_varid(ncid, "LONGITUDE", &varid);
    lon = malloc(nprof * sizeof(double));
    ncw_get_var_double(ncid, varid, lon);

    ncw_inq_varid(ncid, "LATITUDE", &varid);
    lat = malloc(nprof * sizeof(double));
    ncw_get_var_double(ncid, varid, lat);

    if (ncw_var_exists(ncid, "DEPH_ADJUSTED"))
        ncw_inq_varid(ncid, "DEPH_ADJUSTED", &varid);
    else
        ncw_inq_varid(ncid, "PRES_ADJUSTED", &varid);
    z = alloc2d(nprof, nz, sizeof(double));
    ncw_get_var_double(ncid, varid, z[0]);

    if (strncmp(meta->type, "TEM", 3) == 0) {
        validmin = -2.0;
        validmax = 40.0;
        if (ncw_var_exists(ncid, "TEMP_ADJUSTED")) {
            enkf_printf("        reading TEMP_ADJUSTED\n");
            ncw_inq_varid(ncid, "TEMP_ADJUSTED", &varid);
            ncw_inq_varid(ncid, "TEMP_ADJUSTED_QC", &varid_qc);
        } else {
            ncw_close(ncid);
            goto nodata;
        }
    } else if (strncmp(meta->type, "SAL", 3) == 0) {
        validmin = 0;
        validmax = 50.0;
        if (ncw_var_exists(ncid, "PSAL_ADJUSTED")) {
            enkf_printf("        reading PSAL_ADJUSTED\n");
            ncw_inq_varid(ncid, "PSAL_ADJUSTED", &varid);
            ncw_inq_varid(ncid, "PSAL_ADJUSTED_QC", &varid_qc);
        } else {
            ncw_close(ncid);
            goto nodata;
        }
    } else
        enkf_quit("observation type \"%s\" not handled for CMEMS product", meta->type);
    v = alloc2d(nprof, nz, sizeof(double));
    ncw_get_var_double(ncid, varid, v[0]);
    ncw_get_att_double(ncid, varid, "_FillValue", &missval);
    qcflag = alloc2d(nprof, nz, sizeof(char));
    ncw_get_var_text(ncid, varid_qc, qcflag[0]);

    ncw_inq_varid(ncid, "JULD_LOCATION", &varid);
    ncw_get_att_text(ncid, varid, "units", buf);
    tunits_convert(buf, &tunits_multiple, &tunits_offset);
    time = malloc(nprof * sizeof(double));
    ncw_get_var_double(ncid, varid, time);

    ncw_inq_varid(ncid, "WMO_INST_TYPE", &varid);
    insttype = malloc(nprof * WMO_INSTSIZE);
    ncw_get_var_text(ncid, varid, insttype);

    ncw_close(ncid);

    for (p = 0; p < (int) nprof; ++p) {
        char inststr[MAXSTRLEN];
        int instnum;

        strncpy(inststr, &insttype[p * WMO_INSTSIZE], 4);
        inststr[4] = 0;
        if (!str2int(inststr, &instnum))
            instnum = 999;
        snprintf(inststr, MAXSTRLEN, "WMO%03d", instnum);

        if (st_exclude != NULL && st_findindexbystring(st_exclude, inststr) >= 0)
            continue;

        for (i = 0; i < (int) nz; ++i) {
            observation* o;
            obstype* ot;

            if (fabs(v[p][i] - missval) < EPS || v[p][i] < validmin || v[p][i] > validmax)
                continue;
            if (z[p][i] < 0.0)
                continue;
            {
                int qcflagint = (int) qcflag[p][i] - (int) '0';

                assert(qcflagint >= 0 && qcflagint <= 9);
                if (!((1 << qcflagint) & qcflagvals))
                    continue;
            }

            obs_checkalloc(obs);
            o = &obs->data[obs->nobs];

            o->product = st_findindexbystring(obs->products, meta->product);
            assert(o->product >= 0);
            o->type = obstype_getid(obs->nobstypes, obs->obstypes, meta->type, 1);
            ot = &obs->obstypes[o->type];
            o->instrument = st_add_ifabsent(obs->instruments, inststr, -1);
            o->id = obs->nobs;
            o->fid = fid;
            o->batch = p;
            o->value = v[p][i];
            o->std = 0.0;
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
            if ((o->status == STATUS_OK) && (o->lon <= ot->xmin || o->lon >= ot->xmax || o->lat <= ot->ymin || o->lat >= ot->ymax || o->depth <= ot->zmin || o->depth >= ot->zmax))
                o->status = STATUS_OUTSIDEOBSDOMAIN;
            o->model_depth = NAN;       /* set in obs_add() */
            o->date = time[p] * tunits_multiple + tunits_offset;
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

        for (i = 0, ii = 0; i < nprof; ++i) {
            lonlat[ii * 2] = lon[ii];
            lonlat[ii * 2 + 1] = lat[ii];
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

    free(v);
    free(z);
    free(qcflag);
    free(insttype);
  nodata:
    free(lon);
    free(lat);
    if (st_exclude != NULL)
        st_destroy(st_exclude);
}
