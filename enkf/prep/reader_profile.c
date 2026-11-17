/******************************************************************************
 *
 * File:        reader_profile.c        
 *
 * Created:     25/06/2019
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description: Generic reader for profile observations.
 *
 * Revisions:   08/12/2023 PS: Modified to handle multiple profiles.
 *              12/11/2025 PS: Renamed this reader from "z" to "profile".
 *              17/11/2026 PS: Factored common attribute handling into obsatts.
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
    size_t size;
    nc_type type;
    int status = 1;
    int n, i;

    if (!ncw_var_exists(ncid, varname))
        return 1;
    ncw_inq_varid(ncid, varname, &varid);
    ncw_inq_varsize(ncid, varid, &size);
    n = (int) size;
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
void reader_profile(char* fname, int fid, obssection* section, grid* g, observations* obs)
{
    char* varname = NULL;
    int varname_specified = 0;
    char* lonname = NULL;
    char* latname = NULL;
    stringtable* st_znames = NULL;
    obsatts oa;
    char** instruments = NULL;
    int* status = NULL;

    stringtable* st_exclude = NULL;

    int productid = -1;
    int typeid = -1;

    int ncid;
    size_t nz = 0;
    size_t nprof = 0;
    size_t nobs = 0;
    double* lon = NULL;
    double* lat = NULL;
    double** var = NULL;
    double** z = NULL;
    size_t ntime = 0;
    double* time = NULL;
    int varid;
    int npexcluded;
    int p, k, i, nobs_read;

    obsatts_init(&oa);
    ncw_open(fname, NC_NOWRITE, &ncid);

    for (i = 0; i < section->npars; ++i) {
        if (strcasecmp(section->pars[i].name, "VARNAME") == 0) {
            varname_specified = 1;
            if (varname == NULL && ncw_var_exists(ncid, section->pars[i].value) && !allmissing(ncid, section->pars[i].value)) {
                varname = section->pars[i].value;
            }
        } else if (strcasecmp(section->pars[i].name, "LONNAME") == 0)
            lonname = section->pars[i].value;
        else if (strcasecmp(section->pars[i].name, "LATNAME") == 0)
            latname = section->pars[i].value;
        else if (strcasecmp(section->pars[i].name, "ZNAME") == 0) {
            if (st_znames == NULL)
                st_znames = st_create("znames");
            if (ncw_var_exists(ncid, section->pars[i].value) && !allmissing(ncid, section->pars[i].value))
                st_add_ifabsent(st_znames, section->pars[i].value, -1);
        } else if (strcasecmp(section->pars[i].name, "EXCLUDEINST") == 0) {
            if (st_exclude == NULL)
                st_exclude = st_create("exclude");
            st_add(st_exclude, section->pars[i].value, -1);
        } else if (!obsatts_parse_param(&oa, section->pars[i].name, section->pars[i].value))
            enkf_quit("unknown PARAMETER \"%s\"\n", section->pars[i].name);
    }

    if (varname == NULL) {
        if (varname_specified) {
            enkf_printf("        no data of specified type\n");
            ncw_close(ncid);
            goto finish;
        }
        enkf_quit("reader_z(): %s: VARNAME not specified", fname);
    } else
        enkf_printf("        VARNAME = %s\n", varname);

    /*
     * main variable
     */
    ncw_inq_varid(ncid, varname, &varid);
    ncw_inq_varsize(ncid, varid, &nobs);
    if (nobs == 0) {
        enkf_printf("        no observations found\n");
        ncw_close(ncid);
        goto finish;
    } else if (nobs == 1) {
        nprof = 1;
        nz = 1;
    } else {
        int ndims;
        size_t dimlen[4];

        ncw_inq_vardims(ncid, varid, 4, &ndims, dimlen);
        if (ndims == 1) {
            nz = dimlen[0];
            nprof = 1;
        } else if (ndims >= 2) {
            for (i = 2; i < ndims; ++i)
                assert(dimlen[i] == 1);
            nprof = dimlen[0];
            nz = dimlen[1];
        }
        assert(nobs == nprof * nz);
    }
    enkf_printf("        # profiles = %u\n", (unsigned int) nprof);
    enkf_printf("        # levels = %u\n", (unsigned int) nz);
    var = alloc2d(nprof, nz, sizeof(double));
    ncu_readvardouble(ncid, varid, nobs, var[0]);

    /*
     * lon/lat
     */
    lonname = get_lonname(ncid, lonname);
    if (lonname != NULL) {
        enkf_printf("        LONNAME = %s\n", lonname);
        ncw_inq_varid(ncid, lonname, &varid);
    } else
        enkf_quit("reader_z(): %s: could not find longitude variable", fname);
    ncw_check_varsize(ncid, varid, nprof);
    lon = malloc(nprof * sizeof(double));
    ncw_get_var_double(ncid, varid, lon);

    latname = get_latname(ncid, latname);
    if (latname != NULL) {
        enkf_printf("        LATNAME = %s\n", latname);
        ncw_inq_varid(ncid, latname, &varid);
    } else
        enkf_quit("reader_z(): %s: could not find latitude variable", fname);
    ncw_check_varsize(ncid, varid, nprof);
    lat = malloc(nprof * sizeof(double));
    ncw_get_var_double(ncid, varid, lat);

    /*
     * z
     */
    if (st_znames == NULL || st_getsize(st_znames) == 0) {
        char* zname = NULL;

        zname = get_zname(ncid, zname);
        if (zname == NULL) {
            if (st_znames == NULL)
                enkf_quit("reader_z(): ZNAME not specified and no suitable candidate found\n");
            else {
                enkf_printf("%s: no valid Z data, skipping\n", fname);
                ncw_close(ncid);
                goto finish;
            }
        }

        if (st_znames == NULL)
            st_znames = st_create("znames");
        st_add(st_znames, zname, -1);
    }
    {
        int nznames = st_getsize(st_znames);
        double*** zall = alloc3d(nznames, nprof, nz, sizeof(double));

        z = alloc2d(nprof, nz, sizeof(double));
        for (i = 0; i < nznames; ++i) {
            ncw_inq_varid(ncid, st_findstringbyindex(st_znames, i), &varid);
            ncu_readvardouble(ncid, varid, nprof * nz, zall[i][0]);
        }
        for (p = 0; p < (int) nprof; ++p) {
            for (i = 0; i < nznames; ++i) {
                for (k = 0; k < (int) nz; ++k) {
                    if (isfinite(zall[i][p][k]))
                        break;
                }
                if (k < (int) nz) {
                    for (k = 0; k < (int) nz; ++k)
                        z[p][k] = zall[i][p][k];
                    break;
                }
            }
        }
        free(zall);
        st_destroy(st_znames);
        st_znames = NULL;
    }

    /*
     * std, estd, batch, qcflag, instrument
     */
    obsatts_read(&oa, section, ncid, varname, nobs, nprof, nz, 0, 1, obs);

    /*
     * time
     */
    get_time(section, ncid, &ntime, &time);
    assert(ntime == nobs || ntime == nprof || ntime <= 1);

    /*
     * per-profile instrument array override
     */
    if (ncw_var_exists(ncid, oa.instrument)) {
        nc_type type;
        int ndim;
        size_t dimlen[2];
        char** tmp;

        ncw_inq_varid(ncid, oa.instrument, &varid);
        ncw_inq_vartype(ncid, varid, &type);
        if (strcmp(ncw_nctype2str(type), "NC_CHAR") != 0)
            enkf_quit("reader_z(): %s: instrument array \"%s\" is not of type \"char\"", fname, oa.instrument);
        ncw_inq_vardims(ncid, varid, 2, &ndim, dimlen);
        if (ndim != 2)
            enkf_quit("reader_z(): %s: instrument array \"%s\" must have 2 dimensions", fname, oa.instrument);
        assert(dimlen[0] == nprof);
        tmp = alloc2d(dimlen[0], dimlen[1], sizeof(char));
        ncw_get_var(ncid, varid, tmp[0]);
        instruments = alloc2d(dimlen[0], dimlen[1] + 2, sizeof(char));
        for (i = 0; i < (int) dimlen[0]; ++i) {
            int ii;

            for (ii = (int) dimlen[1] - 1; ii >= 0; --ii)
                if (tmp[i][ii] != ' ')
                    break;
            if (ii >= 0)
                strncpy(instruments[i], tmp[i], ii + 1);
            else
                strncpy(instruments[i], "NA", dimlen[1] + 1);
        }
        free(tmp);
    }

    ncw_close(ncid);

    productid = st_findindexbystring(obs->products, section->product);
    assert(productid >= 0);
    typeid = obstype_getid(obs->nobstypes, obs->obstypes, section->type, 1);

    nobs_read = 0;
    npexcluded = 0;
    status = calloc(nprof, sizeof(int));
    for (p = 0, i = 0; p < (int) nprof; ++p) {
        if (instruments != NULL && st_exclude != NULL && st_findindexbystring(st_exclude, instruments[p]) >= 0) {
            npexcluded++;
            i += (int) nz;
            continue;
        }

        for (k = 0; k < (int) nz; ++k, ++i) {
            observation* o;

            if (isnan(var[p][k]) || (ntime == nobs && isnan(time[i])) || (ntime == nprof && isnan(time[p])))
                continue;
            if (!obsatts_valid(&oa, (size_t) i))
                continue;

            nobs_read++;
            obs_checkalloc(obs);
            o = &obs->data[obs->nobs];

            o->product = productid;
            o->type = typeid;
            o->id = obs->nobs;
            o->id_orig = i;
            o->fid = fid;
            o->value = (double) var[p][k];
            obsatts_apply(&oa, (size_t) i, p, o);
            if (instruments != NULL)
                o->instrument = st_add_ifabsent(obs->instruments, instruments[p], -1);
            o->lon = lon[p];
            o->lat = lat[p];
            o->depth = z[p][k];
            o->status = grid_xy2fij(g, o->lon, o->lat, o->fij);
            if (!obs->allobs && o->status == STATUS_OUTSIDEGRID)
                continue;
            if (o->status == STATUS_OK)
                o->status = grid_z2fk_f(g, o->fij, o->depth, &o->fk);
            else
                o->fk = NAN;
            o->model_depth = NAN;       /* set in obs_add() */
            if (ntime == 0)
                o->time = NAN;
            else if (ntime == 1)
                o->time = time[0];
            else if (ntime == nprof)
                o->time = time[p];
            else if (ntime == nobs)
                o->time = time[i];
            else
                enkf_quit("programming error");

            o->aux = -1;
            if (o->status == STATUS_OK)
                status[p] = 1;

            obs->nobs++;
        }
    }
    if (npexcluded > 0)
        enkf_printf("        # profiles excluded by inst. type = %d\n", npexcluded);
    /*
     * get the number of unique profile locations
     */
    if (nobs_read > 0) {
        double* lonlat = malloc(nprof * sizeof(double) * 2);
        int ngood, ii;

        ngood = 0;
        for (i = 0; i < (int) nprof; ++i) {
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

  finish:

    if (st_znames != NULL)
        st_destroy(st_znames);
    if (var != NULL)
        free(var);
    if (lon != NULL)
        free(lon);
    if (lat != NULL)
        free(lat);
    if (z != NULL)
        free(z);
    if (time != NULL)
        free(time);
    obsatts_free(&oa);
    if (instruments != NULL)
        free(instruments);
    if (st_exclude != NULL)
        st_destroy(st_exclude);
    if (status != NULL)
        free(status);
}

/**
 */
void reader_profile_describe(void)
{
    enkf_printf("\n  Generic reader \"profile\" reads profile data. It currently assumes the following:\n\
    - longitude and latitude are either 0-dimensional or 1-dimensional of size\n\
      [nprof];\n\
    - time is either 0- or 1-dimensional of size 1, or 1-dimensional of size\n\
      [nprof], or 1-dimensional of size [nprof * nz], or 2-dimensional of size\n\
      [nprof][nz];\n\
    - profile variables are either 1-dimensional of size [nz] or 2-dimensional\n\
      of size [nprof][nz]\n");
    describe_commongenericreaderparams();
    enkf_printf("  Parameters specific to the reader:\n\
    - EXCLUDEINST (-)\n\
        tag of the instrument to be skipped (can be repeated)\n");
    describe_commonreaderparams();
}
