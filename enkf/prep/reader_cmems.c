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
 *              Optional parameters:
 *              - QCFLAGVALS
 *                  the list of allowed values of QC flag variables, from the
 *                  interval [0,31]; by default PARAMETER QCFLAGVALS = 1
 *              - EXCLUDEINST
 *                  WMO instrument type id (p. 27 of the above manual) to be
 *                  skipped, one entry (one line) per instrument, e.g.
 *                  PARAMETER EXCLUDEINST = WMO995 # marine mammals
 *                  PARAMETER EXCLUDEINST = WMO999 # unknown
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
#include "ncutils.h"
#include "definitions.h"
#include "utils.h"
#include "obsprm.h"
#include "grid.h"
#include "model.h"
#include "observations.h"
#include "prep_utils.h"
#include "allreaders.h"

#define WMO_INSTSIZE 4
#define QCFLAGVALS_DEF 2        /* 0b00000000000000000000000000000010
                                 * corresponds to QCFLAG = 1 */
#define QCFLAGVALMAX 9

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
static int allmissing(int ncid, char varname[])
{
    int varid;
    int ndims;
    size_t dimlen[2];
    nc_type type;
    int status = 1;
    int n, i;

    if (!ncw_var_exists(ncid, varname))
        return 1;
    ncw_inq_varid(ncid, varname, &varid);
    ncw_inq_vardims(ncid, varid, 2, &ndims, dimlen);
    n = dimlen[0] * dimlen[1];
    assert(ndims == 2);
    ncw_inq_vartype(ncid, varid, &type);
    if (type != NC_CHAR) {
        double* v = malloc(n * sizeof(double));

        ncu_readvardouble(ncid, varid, n, v);
        for (i = 0; i < n; ++i)
            if (isfinite(v[i])) {
                status = 0;
                break;
            }
        free(v);
    } else {
        char* v = malloc(n);
        int nofill;
        char fillval;

        ncw_get_var(ncid, varid, v);
        ncw_inq_var_fill(ncid, varid, &nofill, &fillval);
        if (nofill == 0)
            for (i = 0; i < n; ++i)
                if (v[i] != fillval) {
                    status = 0;
                    break;
                }
        free(v);
    }

    return status;
}

/**
 */
static int datamissing(double* v, int n)
{
    int i;

    for (i = 0; i < n; ++i)
        if (isfinite(v[i]))
            return 0;
    return 1;
}

/**
 */
void reader_cmems(char* fname, int fid, obsmeta* meta, grid* g, observations* obs)
{
    stringtable* st_exclude = NULL;
    int ncid, varid;
    int dimid_nprof, dimid_nz;
    size_t nprof, nz;
    int* status_prof;
    double* lon;
    double* lat;
    double** z;
    double** v_adj = NULL;
    char** qcflag_adj = NULL;
    double** v_nonadj = NULL;
    char** qcflag_nonadj = NULL;
    uint32_t qcflagvals = QCFLAGVALS_DEF;
    int qcflagcounts[QCFLAGVALMAX + 1];
    double* time;
    double validmin = DBL_MAX;
    double validmax = -DBL_MAX;
    char* insttype;
    char buf[MAXSTRLEN];
    double tunits_multiple, tunits_offset;
    int npexcluded;
    int p, i, nobs_read, id, nprof_adj, nprof_nonadj;

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
                if (val < 0 || val > QCFLAGVALMAX)
                    enkf_quit("%s: QCFLAGVALS entry = %d (supposed to be in [0,%d] interval", meta->prmfname, val, QCFLAGVALMAX);
                qcflagvals |= 1 << val;
                line = NULL;
            }
            free(line);
            if (qcflagvals == 0)
                enkf_quit("%s: no valid flag entries found after QCFLAGVALS\n", meta->prmfname);
            enkf_printf("        QCFLAGS used =");
            for (i = 0; i <= QCFLAGVALMAX; ++i)
                if (qcflagvals & (1 << i))
                    enkf_printf(" %d", i);
            enkf_printf("\n");
        } else
            enkf_quit("unknown PARAMETER \"%s\"\n", meta->pars[i].name);
    }
    if (st_exclude != NULL) {
        enkf_printf("        excluding instruments:");
        st_printentries(st_exclude, " ");
        enkf_printf("\n");
    }

    if (meta->nestds == 0)
        enkf_quit("ERROR_STD is necessary but not specified for product \"%s\"", meta->product);

    ncw_open(fname, NC_NOWRITE, &ncid);

    ncw_inq_dimid(ncid, "N_PROF", &dimid_nprof);
    ncw_inq_dimlen(ncid, dimid_nprof, &nprof);
    enkf_printf("        # profiles = %u\n", (unsigned int) nprof);
    if (nprof == 0) {
        enkf_printf("        no profiles found\n");
        ncw_close(ncid);
        goto noprofiles;
    }
    status_prof = calloc(nprof, sizeof(int));

    ncw_inq_varid(ncid, "LONGITUDE", &varid);
    lon = malloc(nprof * sizeof(double));
    ncw_get_var_double(ncid, varid, lon);

    ncw_inq_varid(ncid, "LATITUDE", &varid);
    lat = malloc(nprof * sizeof(double));
    ncw_get_var_double(ncid, varid, lat);

    ncw_inq_dimid(ncid, "N_LEVELS", &dimid_nz);
    ncw_inq_dimlen(ncid, dimid_nz, &nz);
    enkf_printf("        # levels = %u\n", (unsigned int) nz);
    z = alloc2d(nprof, nz, sizeof(double));
    {
        double* zall[4] = { NULL, NULL, NULL, NULL };
        char* zall_qc[4] = { NULL, NULL, NULL, NULL };
        char znames[4][20] = { "DEPH_ADJUSTED",
            "PRES_ADJUSTED",
            "DEPH",
            "PRES"
        };
        char znames_qc[4][20] = { "DEPH_ADJUSTED_QC",
            "PRES_ADJUSTED_QC",
            "DEPH_QC",
            "PRES_QC"
        };
        int j;

        for (i = 0; i < 4; ++i) {
            int ii;

            if (!ncw_var_exists(ncid, znames[i]) || !ncw_var_exists(ncid, znames_qc[i]))
                continue;
            ncw_inq_varid(ncid, znames[i], &varid);
            zall[i] = malloc(nprof * nz * sizeof(double));
            ncu_readvardouble(ncid, varid, nprof * nz, zall[i]);
            ncw_inq_varid(ncid, znames_qc[i], &varid);
            zall_qc[i] = malloc(nprof * nz * sizeof(char));
            ncw_get_var_text(ncid, varid, zall_qc[i]);

            for (ii = 0; ii < nz * nprof; ++ii)
                if ((int) zall_qc[i][ii] - (int) '0' != 1)
                    zall[i][ii] = NAN;
        }

        for (i = 0; i < nprof * nz; ++i) {
            for (j = 0; j < 4; ++j) {
                if (zall[j] == NULL)
                    continue;
                if (isfinite(zall[j][i]))
                    break;
            }
            z[0][i] = (j < 4) ? zall[j][i] : NAN;
        }
        for (i = 0; i < 4; ++i)
            if (zall[i] != NULL) {
                free(zall[i]);
                free(zall_qc[i]);
            }
    }

    if (strncmp(meta->type, "TEM", 3) == 0) {
        validmin = -2.0;
        validmax = 40.0;
        if (!allmissing(ncid, "TEMP_ADJUSTED")) {
            ncw_inq_varid(ncid, "TEMP_ADJUSTED", &varid);
            v_adj = alloc2d(nprof, nz, sizeof(double));
            ncu_readvardouble(ncid, varid, nz * nprof, v_adj[0]);
            ncw_inq_varid(ncid, "TEMP_ADJUSTED_QC", &varid);
            qcflag_adj = alloc2d(nprof, nz, sizeof(char));
            ncw_get_var_text(ncid, varid, qcflag_adj[0]);
        }
        if (!allmissing(ncid, "TEMP")) {
            ncw_inq_varid(ncid, "TEMP", &varid);
            v_nonadj = alloc2d(nprof, nz, sizeof(double));
            ncu_readvardouble(ncid, varid, nz * nprof, v_nonadj[0]);
            ncw_inq_varid(ncid, "TEMP_QC", &varid);
            qcflag_nonadj = alloc2d(nprof, nz, sizeof(char));
            ncw_get_var_text(ncid, varid, qcflag_nonadj[0]);
        }
        if (v_adj == NULL && v_nonadj == NULL) {
            enkf_printf("        no data of specified type\n");
            ncw_close(ncid);
            goto nodata;
        }
    } else if (strncmp(meta->type, "SAL", 3) == 0) {
        validmin = 0;
        validmax = 50.0;
        if (!allmissing(ncid, "PSAL_ADJUSTED")) {
            ncw_inq_varid(ncid, "PSAL_ADJUSTED", &varid);
            v_adj = alloc2d(nprof, nz, sizeof(double));
            ncu_readvardouble(ncid, varid, nz * nprof, v_adj[0]);
            ncw_inq_varid(ncid, "PSAL_ADJUSTED_QC", &varid);
            qcflag_adj = alloc2d(nprof, nz, sizeof(char));
            ncw_get_var_text(ncid, varid, qcflag_adj[0]);
        }
        if (!allmissing(ncid, "PSAL")) {
            ncw_inq_varid(ncid, "PSAL", &varid);
            v_nonadj = alloc2d(nprof, nz, sizeof(double));
            ncu_readvardouble(ncid, varid, nz * nprof, v_nonadj[0]);
            ncw_inq_varid(ncid, "PSAL_QC", &varid);
            qcflag_nonadj = alloc2d(nprof, nz, sizeof(char));
            ncw_get_var_text(ncid, varid, qcflag_nonadj[0]);
        }
        if (v_adj == NULL && v_nonadj == NULL) {
            enkf_printf("        no data of specified type\n");
            ncw_close(ncid);
            goto nodata;
        }
    } else
        enkf_quit("observation type \"%s\" not handled for CMEMS product", meta->type);

    ncw_inq_varid(ncid, "JULD", &varid);
    ncw_get_att_text(ncid, varid, "units", buf);
    tunits_convert(buf, &tunits_multiple, &tunits_offset);
    time = malloc(nprof * sizeof(double));
    ncw_get_var_double(ncid, varid, time);

    ncw_inq_varid(ncid, "WMO_INST_TYPE", &varid);
    insttype = malloc(nprof * WMO_INSTSIZE);
    ncw_get_var_text(ncid, varid, insttype);

    ncw_close(ncid);

    npexcluded = 0;
    memset(qcflagcounts, 0, (QCFLAGVALMAX + 1) * sizeof(int));
    nobs_read = 0;
    nprof_adj = 0;
    nprof_nonadj = 0;
    for (p = 0, id = 0; p < (int) nprof; ++p) {
        double* v;
        char* qcflag;
        char inststr[MAXSTRLEN];

        {
            int instnum;

            strncpy(inststr, &insttype[p * WMO_INSTSIZE], 4);
            inststr[4] = 0;
            if (!str2int(inststr, &instnum))
                instnum = 999;
            snprintf(inststr, MAXSTRLEN, "WMO%03d", instnum);

            if (st_exclude != NULL && st_findindexbystring(st_exclude, inststr) >= 0) {
                npexcluded++;
                id += nz;
                continue;
            }
        }

        {
            if (v_adj != NULL && !datamissing(v_adj[p], nz)) {
                v = v_adj[p];
                qcflag = qcflag_adj[p];
                nprof_adj++;
            } else if (v_nonadj != NULL && !datamissing(v_nonadj[p], nz)) {
                v = v_nonadj[p];
                qcflag = qcflag_nonadj[p];
                nprof_nonadj++;
            } else {
                id += nz;
                continue;
            }
        }

        for (i = 0; i < (int) nz; ++i, ++id) {
            observation* o;
            int qcflagint;

            if (isnan(v[i]) || v[i] < validmin || v[i] > validmax)
                continue;
            if (isnan(z[p][i]) || z[p][i] < 0.0)
                continue;
            {
                qcflagint = (int) qcflag[i] - (int) '0';
                if (qcflagint < 0 || qcflagint > QCFLAGVALMAX)  /* impossible 
                                                                 */
                    continue;
                if (!((1 << qcflagint) & qcflagvals))
                    continue;
            }
            qcflagcounts[qcflagint]++;
            nobs_read++;

            obs_checkalloc(obs);
            o = &obs->data[obs->nobs];

            o->product = st_findindexbystring(obs->products, meta->product);
            assert(o->product >= 0);
            o->type = obstype_getid(obs->nobstypes, obs->obstypes, meta->type, 1);
            o->instrument = st_add_ifabsent(obs->instruments, inststr, -1);
            o->id = obs->nobs;
            o->id_orig = id;
            o->fid = fid;
            o->batch = p;
            o->value = v[i];
            o->estd = 0.0;
            o->lon = lon[p];
            o->lat = lat[p];
            o->depth = z[p][i];
            o->status = grid_xy2fij(g, o->lon, o->lat, o->fij);
            if (!obs->allobs && o->status == STATUS_OUTSIDEGRID)
                break;
            if (o->status == STATUS_OK)
                o->status = grid_z2fk_f(g, o->fij, o->depth, &o->fk);
            else
                o->fk = NAN;
            o->model_depth = NAN;       /* set in obs_add() */
            o->time = time[p] * tunits_multiple + tunits_offset;
            o->aux = -1;
            if (o->status == STATUS_OK)
                status_prof[p] = 1;

            obs->nobs++;
        }
    }
    if (npexcluded > 0)
        enkf_printf("        # profiles excluded by inst. type = %d\n", npexcluded);
    if (nprof_adj > 0 && nprof_adj != nprof)
        enkf_printf("        processed %d adjusted profiles\n", nprof_adj);
    if (nprof_nonadj > 0 && nprof_nonadj != nprof)
        enkf_printf("        processed %d non-adjusted profiles\n", nprof_nonadj);

    /*
     * get the number of unique profile locations
     */
    if (nobs_read > 0) {
        double* lonlat = malloc(nprof * sizeof(double) * 2);
        int ngood, ii;

        ngood = 0;
        for (i = 0; i < nprof; ++i) {
            if (status_prof[i] == 0)
                continue;
            lonlat[ngood * 2] = lon[i];
            lonlat[ngood * 2 + 1] = lat[i];
            ngood++;
        }
        enkf_printf("        # profiles with data to process = %d\n", ngood);
        if (ngood > 1) {
            int nunique = 1;

            qsort(lonlat, ngood, sizeof(double) * 2, cmp_lonlat);
            nunique = 1;
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

    enkf_printf("        nobs = %d\n", nobs_read);
    if (nobs_read > 0) {
        enkf_printf("        # of processed obs by quality flags:\n");
        for (i = 0; i <= QCFLAGVALMAX; ++i)
            if (qcflagvals & (1 << i) && qcflagcounts[i] > 0)
                enkf_printf("          %d: %d obs\n", i, qcflagcounts[i]);
    }

    free(time);
    if (v_adj != NULL) {
        free(v_adj);
        free(qcflag_adj);
    }
    if (v_nonadj != NULL) {
        free(v_nonadj);
        free(qcflag_nonadj);
    }
    free(insttype);
  nodata:
    free(lon);
    free(lat);
    free(status_prof);
    free(z);
  noprofiles:
    if (st_exclude != NULL)
        st_destroy(st_exclude);
}

/**
 */
void reader_cmems_describe(void)
{
    enkf_printf("\n  Reader \"cmems\" reads in-situ temperature and salinity data produced by\n\
Copernicus Marine Environment Monitoring Service.\n\
\n\
  There are a number of parameters that must (marked below with \"++\"), can\n\
  (\"+\"), or may (\"-\") be specified in the corresponding section of the\n\
  observation data parameter file. The names in brackets represent the default\n\
  names checked in the abscence of the entry for the parameter. Each parameter\n\
  needs to be entered as follows:\n\
    PARAMETER <name> = <value> ...\n\
\n\
  Parameters specific to the reader:\n\
    - QCFLAGVALS (1) (-)\n\
        the list of allowed values of QC flag variable\n\
        Note: the name of the QC flag variable is set by the code after\n\
        determining the name of the principal variable\n\
    - EXCLUDEINST (-)\n\
        instrument in format \"WMO*\" to be skipped\n");
    describe_commonreaderparams();
}
