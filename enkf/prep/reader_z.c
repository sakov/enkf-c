/******************************************************************************
 *
 * File:        reader_z.c        
 *
 * Created:     25/06/2019
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description: Generic reader for observations in a single profile.
 *
 * Revisions:   08/12/2023 PS: Modified to handle multiple profiles.
 *
 *****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "ncw.h"
#include "ncutils.h"
#include "definitions.h"
#include "utils.h"
#include "obsprm.h"
#include "model.h"
#include "grid.h"
#include "observations.h"
#include "prep_utils.h"
#include "allreaders.h"

#define TYPE_DOUBLE 0
#define TYPE_SHORT 1

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
    n = (ndims == 1) ? dimlen[0] : dimlen[0] * dimlen[1];
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
void reader_z(char* fname, int fid, obsmeta* meta, grid* g, observations* obs)
{
    char* varname = NULL;
    int varname_specified = 0;
    char* lonname = NULL;
    char* latname = NULL;
    stringtable* st_znames = NULL;
    char* estdname = NULL;
    char instrument[MAXSTRLEN] = "";
    char** instruments = NULL;
    int* status = NULL;

    stringtable* st_exclude = NULL;

    int nqcflagvars = 0;
    char** qcflagvarnames = NULL;
    uint32_t* qcflagmasks = NULL;

    int instid = -1;
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
    double var_estd = NAN;
    double** estd = NULL;
    uint32_t*** qcflag = NULL;
    size_t ntime = 0;
    double* time = NULL;
    int varid;
    int npexcluded;
    int p, k, i, nobs_read;

    ncw_open(fname, NC_NOWRITE, &ncid);

    for (i = 0; i < meta->npars; ++i) {
        if (strcasecmp(meta->pars[i].name, "VARNAME") == 0) {
            varname_specified = 1;
            if (varname == NULL && ncw_var_exists(ncid, meta->pars[i].value) && !allmissing(ncid, meta->pars[i].value)) {
                varname = meta->pars[i].value;
            }
        } else if (strcasecmp(meta->pars[i].name, "LONNAME") == 0)
            lonname = meta->pars[i].value;
        else if (strcasecmp(meta->pars[i].name, "LATNAME") == 0)
            latname = meta->pars[i].value;
        else if (strcasecmp(meta->pars[i].name, "ZNAME") == 0) {
            if (st_znames == NULL)
                st_znames = st_create("znames");
            if (ncw_var_exists(ncid, meta->pars[i].value) && !allmissing(ncid, meta->pars[i].value))
                st_add_ifabsent(st_znames, meta->pars[i].value, -1);
        } else if (strcasecmp(meta->pars[i].name, "ESTDNAME") == 0)
            estdname = meta->pars[i].value;
        else if (strcasecmp(meta->pars[i].name, "INSTRUMENT") == 0)
            strncpy(instrument, meta->pars[i].value, MAXSTRLEN - 1);
        else if (strcasecmp(meta->pars[i].name, "TIMENAME") == 0 || strcasecmp(meta->pars[i].name, "TIMENAMES") == 0)
            /*
             * TIMENAME and TIMENAMES are dealt with separately
             */
            ;
        else if (strcasecmp(meta->pars[i].name, "QCFLAGNAME") == 0 || strcasecmp(meta->pars[i].name, "QCFLAGVARNAME") == 0 || strcasecmp(meta->pars[i].name, "QCFLAGVALS") == 0)
            /*
             * QCFLAGNAME and QCFLAGVALS are dealt with separately
             */
            ;
        else if (strcasecmp(meta->pars[i].name, "EXCLUDEINST") == 0) {
            if (st_exclude == NULL)
                st_exclude = st_create("exclude");
            st_add(st_exclude, meta->pars[i].value, -1);
        } else
            enkf_quit("unknown PARAMETER \"%s\"\n", meta->pars[i].name);
    }

    if (varname == NULL) {
        if (varname_specified) {
            enkf_printf("        no data of specified type\n");
            ncw_close(ncid);
            return;
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
        for (p = 0; p < nprof; ++p) {
            for (i = 0; i < nznames; ++i) {
                for (k = 0; k < nz; ++k) {
                    if (isfinite(zall[i][p][k]))
                        break;
                }
                if (k < nz) {
                    for (k = 0; k < nz; ++k)
                        z[p][k] = zall[i][p][k];
                    break;
                }
            }
        }
        free(zall);
        st_destroy(st_znames);
    }

    /*
     * estd
     */
    varid = -1;
    if (estdname != NULL)
        ncw_inq_varid(ncid, estdname, &varid);
    else if (ncw_var_exists(ncid, "error_std"))
        ncw_inq_varid(ncid, "error_std", &varid);
    if (varid >= 0) {
        int ndims;
        size_t dimlen[2];

        estd = alloc2d(nprof, nz, sizeof(double));

        ncw_inq_vardims(ncid, varid, 2, &ndims, dimlen);
        if (ndims == 1) {
            /*
             * estd for prifile
             */
            assert(dimlen[0] == nprof);
            ncu_readvardouble(ncid, varid, nprof, estd[0]);
            for (p = 0; p < nprof; ++p)
                estd[p][0] = estd[0][p];
            for (p = 0; p < nprof; ++p)
                for (k = 1; k < nz; ++k)
                    estd[p][k] = estd[p][0];
        } else if (ndims == 2) {
            /*
             * estd for each ob.
             */
            assert(dimlen[0] == nprof && dimlen[1] == nz);
            ncu_readvardouble(ncid, varid, nobs, estd[0]);
        }
    }

    if (estd == NULL) {
        ncw_inq_varid(ncid, varname, &varid);
        if (ncw_att_exists(ncid, varid, "error_std")) {
            ncw_check_attlen(ncid, varid, "error_std", 1);
            ncw_get_att_double(ncid, varid, "error_std", &var_estd);
        }
    }

    /*
     * qcflag
     */
    get_qcflags(meta, &nqcflagvars, &qcflagvarnames, &qcflagmasks);
    if (nqcflagvars > 0) {
        qcflag = alloc3d(nqcflagvars, nprof, nz, sizeof(int32_t));
        for (i = 0; i < nqcflagvars; ++i) {
            uint32_t** qcflagi = qcflag[i];
            int ndims;
            size_t dimlen[2];
            nc_type type;

            ncw_inq_varid(ncid, qcflagvarnames[i], &varid);
            ncw_inq_vardims(ncid, varid, 2, &ndims, dimlen);
            ncw_inq_vartype(ncid, varid, &type);
            if (type != NC_CHAR) {
                if (ndims == 1) {
                    if (nprof > 1) {
                        /*
                         * flag for profile
                         */
                        assert(dimlen[0] == nprof);
                        ncw_get_var_uint(ncid, varid, qcflagi[0]);
                        for (p = 0; p < nprof; ++p)
                            qcflagi[p][0] = qcflagi[0][p];
                        for (p = 0; p < nprof; ++p)
                            for (k = 1; k < nz; ++k)
                                qcflagi[p][k] = qcflagi[p][0];
                    }
                } else if (ndims == 2) {
                    /*
                     * flag for each ob.
                     */
                    assert(dimlen[0] == nprof && dimlen[1] == nz);
                    ncw_check_varsize(ncid, varid, nobs);
                    ncw_get_var_uint(ncid, varid, qcflagi[0]);
                }
                ncw_inq_vartype(ncid, varid, &type);
            } else {
                char* buf = NULL;
                int o;

                ncw_check_varsize(ncid, varid, nobs);
                buf = calloc(nobs + 1, 1);
                ncw_get_var_text(ncid, varid, buf);
                for (o = 0; o < nobs; ++o)
                    qcflagi[0][o] = (uint32_t) buf[o] - (uint32_t) '0';
            }
        }
    }

    /*
     * time
     */
    get_time(meta, ncid, &ntime, &time);
    assert(ntime == nobs || ntime == nprof || ntime <= 1);

    /*
     * instrument
     */
    if (!ncw_var_exists(ncid, instrument)) {
        /*
         * use parameter INSTRUMENT as the instrument tag
         */
        if (strlen(instrument) == 0 && !get_insttag(ncid, varname, instrument))
            strncpy(instrument, meta->product, MAXSTRLEN - 1);
    } else {
        int varid;
        nc_type type;
        int ndim;
        size_t dimlen[2];
        char** tmp;

        ncw_inq_varid(ncid, instrument, &varid);
        ncw_inq_vartype(ncid, varid, &type);
        if (strcmp(ncw_nctype2str(type), "NC_CHAR") != 0)
            enkf_quit("reader_z(): %s: instrument array \"%s\" is not of type \"char\"", fname, instrument);
        ncw_inq_vardims(ncid, varid, 2, &ndim, dimlen);
        if (ndim != 2)
            enkf_quit("reader_z(): %s: instrument array \"%s\" must have 2 dimensions", fname, instrument);
        assert(dimlen[0] == nprof);
        tmp = alloc2d(dimlen[0], dimlen[1], sizeof(char));
        ncw_get_var(ncid, varid, tmp[0]);
        instruments = alloc2d(dimlen[0], dimlen[1] + 2, sizeof(char));
        for (i = 0; i < dimlen[0]; ++i) {
            int ii;

            for (ii = dimlen[1] - 1; ii >= 0; --ii)
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

    if (instruments == NULL)
        instid = st_add_ifabsent(obs->instruments, instrument, -1);
    productid = st_findindexbystring(obs->products, meta->product);
    assert(productid >= 0);
    typeid = obstype_getid(obs->nobstypes, obs->obstypes, meta->type, 1);

    nobs_read = 0;
    npexcluded = 0;
    status = calloc(nprof, sizeof(int));
    for (p = 0, i = 0; p < (int) nprof; ++p) {
        if (instruments != NULL && st_exclude != NULL && st_findindexbystring(st_exclude, instruments[p]) >= 0) {
            npexcluded++;
            continue;
        }

        for (k = 0; k < (int) nz; ++k, ++i) {
            observation* o;
            int ii;

            if (isnan(var[p][k]) || (estd != NULL && isnan(estd[p][k])) || (ntime == nobs && isnan(time[i])) || (ntime == nprof && isnan(time[p])))
                continue;
            for (ii = 0; ii < nqcflagvars; ++ii)
                if (!((1 << qcflag[ii][p][k]) & qcflagmasks[ii]))
                    goto nextob;

            nobs_read++;
            obs_checkalloc(obs);
            o = &obs->data[obs->nobs];

            o->product = productid;
            o->type = typeid;
            if (instruments == NULL)
                o->instrument = instid;
            else
                o->instrument = st_add_ifabsent(obs->instruments, instruments[p], -1);
            o->id = obs->nobs;
            o->fid = fid;
            o->batch = p;
            o->value = (double) var[p][k];
            if (estd == NULL)
                o->estd = var_estd;
            else
                o->estd = estd[p][k];
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
          nextob:
            ;
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

    free(var);
    free(lon);
    free(lat);
    if (z != NULL)
        free(z);
    if (estd != NULL)
        free(estd);
    if (time != NULL)
        free(time);
    if (nqcflagvars > 0) {
        free(qcflagvarnames);
        free(qcflagmasks);
        free(qcflag);
    }
    if (instruments != NULL)
        free(instruments);
    if (st_exclude != NULL)
        st_destroy(st_exclude);
    if (status != NULL)
        free(status);
}

/**
 */
void reader_z_describe(void)
{
    enkf_printf("\n  Generic reader \"z\" reads profile data. It currently assumes the following:\n\
    - longitude and latitude are either 0-dimensional or 1-dimensional of size\n\
      [nprof];\n\
    - time is either 0- or 1-dimensional of size 1, or 1-dimensional of size\n\
      [nprof], or 1-dimensional of size [nprof * nz], or 2-dimensional of size\n\
      [nprof][nz];\n\
    - profile variables are either 1-dimensional of size [nz] or 2-dimensional\n\
      of size [nprof][nz]\n\
\n\
  There are a number of parameters that must (marked below with \"++\"), can\n\
  (\"+\"), or may (\"-\") be specified in the corresponding section of the\n\
  observation data parameter file. The names in brackets represent the default\n\
  names checked in the abscence of the entry for the parameter. Each parameter\n\
  needs to be entered as follows:\n\
    PARAMETER <name> = <value> ...\n\
\n\
  Parameters common to generic readers:\n\
    - VARNAME (++)\n\
        multiple entries are possible; in the first valid entry is used\n\
    - TIMENAME (\"t\" | \"[tT]ime\" | \"TIME\") (+)\n\
    - or TIMENAMES (when time = base_time + offset) (+)\n\
    - LONNAME (\"lon\" | \"[lL]ongitude\" | \"LONGITUDE\") (+)\n\
    - LATNAME (\"lat\" | \"[lL]atitude\" | \"LATITUDE\") (+)\n\
    - ZNAME (\"z\" | \"[dD]epth\" | \"DEPTH\") (+)\n\
        multiple entries are possible; for each profile the first valid entry is\n\
        used\n\
    - ESTDNAME (\"error_std\") (-)\n\
        error STD; if absent then needs to be specified in the corresponding\n\
        section of the observation data parameter file\n\
    - INSTRUMENT (-)\n\
        either text variable [nprof][len] with instrument tag for each profile\n\
        or (if absent) the instrument tag that will be used for calculating\n\
        instrument stats (overrides the global attribute \"instrument\" in the\n\
        data file)\n\
    - QCFLAGNAME (-)\n\
        name of the QC flag variable, possible values 0 <= qcflag <= 31\n\
    - QCFLAGVALS (-)\n\
        the list of allowed values of QC flag variable\n\
        Note: it is possible to have multiple entries of QCFLAGNAME and\n\
        QCFLAGVALS combination, e.g.:\n\
          PARAMETER QCFLAGNAME = TEMP_quality_control\n\
          PARAMETER QCFLAGVALS = 1\n\
          PARAMETER QCFLAGNAME = DEPTH_quality_control\n\
          PARAMETER QCFLAGVALS = 1\n\
          PARAMETER QCFLAGNAME = LONGITUDE_quality_control\n\
          PARAMETER QCFLAGVALS = 1,8\n\
          PARAMETER QCFLAGNAME = LATITUDE_quality_control\n\
          PARAMETER QCFLAGVALS = 1,8\n\
        An observation is considered valid if each of the specified flags takes\n\
        a permitted value.\n");
    describe_commonreaderparams();
}
