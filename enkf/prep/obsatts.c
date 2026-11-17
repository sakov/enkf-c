/******************************************************************************
 *
 * File:        obsatts.c
 *
 * Created:     17/11/2026
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description:
 *
 * Revisions:   Observation reader helper
 *
 *****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include "ncw.h"
#include "ncutils.h"
#include "global.h"
#include "utils.h"
#include "observations.h"
#include "obsprm.h"
#include "prep_utils.h"
#include "allreaders.h"
#include "obsatts.h"

/**
 */
static int oa_allmissing(size_t n, const double v[])
{
    size_t i;

    for (i = 0; i < n; ++i)
        if (isfinite(v[i]))
            return 0;
    return 1;
}

/** Reads a double variable of size nprof * nz or nprof (broadcasting per
 ** profile) into a flat array v[nprof * nz].
 */
static void oa_readvar_double(int ncid, int varid, size_t nprof, size_t nz, double* v)
{
    size_t size;
    size_t nobs = nprof * nz;

    ncw_inq_varsize(ncid, varid, &size);
    if (size == nobs) {
        ncu_readvardouble(ncid, varid, size, v);
    } else if (nz > 1 && size == nprof) {
        size_t p, k;

        ncu_readvardouble(ncid, varid, size, v);
        for (p = nprof; p > 0; --p) {
            double val = v[p - 1];

            for (k = 0; k < nz; ++k)
                v[(p - 1) * nz + k] = val;
        }
    } else {
        char vname[NC_MAX_NAME + 1];

        ncw_inq_varname(ncid, varid, vname);
        if (nz > 1)
            enkf_quit("%s: %s: size = %zu; expected size = %zu (per observation) or %zu (per profile)", ncw_get_path(ncid), vname, size, nobs, nprof);
        else
            enkf_quit("%s: %s: size = %zu; expected size = %zu", ncw_get_path(ncid), vname, size, nobs);
    }
}

/** Same as oa_readvar_double(), but for int.
 */
static void oa_readvar_int(int ncid, int varid, size_t nprof, size_t nz, int* v)
{
    size_t size;
    size_t nobs = nprof * nz;

    ncw_inq_varsize(ncid, varid, &size);
    if (size == nobs) {
        ncw_get_var_int(ncid, varid, v);
    } else if (nz > 1 && size == nprof) {
        size_t p, k;

        ncw_get_var_int(ncid, varid, v);
        for (p = nprof; p > 0; --p) {
            int val = v[p - 1];

            for (k = 0; k < nz; ++k)
                v[(p - 1) * nz + k] = val;
        }
    } else {
        char vname[NC_MAX_NAME + 1];

        ncw_inq_varname(ncid, varid, vname);
        if (nz > 1)
            enkf_quit("%s: %s: size = %zu; expected size = %zu (per observation) or %zu (per profile)", ncw_get_path(ncid), vname, size, nobs, nprof);
        else
            enkf_quit("%s: %s: size = %zu; expected size = %zu", ncw_get_path(ncid), vname, size, nobs);
    }
}

/**
 */
void obsatts_init(obsatts* oa)
{
    memset(oa, 0, sizeof(obsatts));
    oa->var_estd = NAN;
    oa->instid = -1;
}

/** Parses common generic reader parameters. Returns 1 if handled, 0 otherwise.
 */
int obsatts_parse_param(obsatts* oa, const char* name, char* value)
{
    if (strcasecmp(name, "STDNAME") == 0)
        oa->stdname = value;
    else if (strcasecmp(name, "ESTDNAME") == 0)
        oa->estdname = value;
    else if (strcasecmp(name, "BATCHNAME") == 0)
        oa->batchname = value;
    else if (strcasecmp(name, "INSTATTNAME") == 0)
        oa->instattname = value;
    else if (strcasecmp(name, "INSTPREFIX") == 0)
        oa->instprefix = value;
    else if (strcasecmp(name, "INSTRUMENT") == 0) {
        strncpy(oa->instrument, value, MAXSTRLEN - 1);
        oa->instrument[MAXSTRLEN - 1] = '\0';
    } else if (strcasecmp(name, "ADDVAR") == 0 || strcasecmp(name, "SUBVAR") == 0) {
        int sign = (strcasecmp(name, "ADDVAR") == 0) ? 1 : -1;

        oa->addvars = realloc(oa->addvars, (oa->naddvar + 1) * sizeof(obsaddvar));
        oa->addvars[oa->naddvar].varname = strdup(value);
        oa->addvars[oa->naddvar].sign = sign;
        oa->addvars[oa->naddvar].v = NULL;
        enkf_printf("      %s = %s\n", (sign > 0) ? "ADDVAR" : "SUBVAR", value);
        oa->naddvar++;
    } else if (strcasecmp(name, "TIMENAME") == 0 || strcasecmp(name, "TIMENAMES") == 0 || strcasecmp(name, "QCFLAGNAME") == 0 || strcasecmp(name, "QCFLAGVARNAME") == 0 || strcasecmp(name, "QCFLAGVALS") == 0 || strcasecmp(name, "LOCATION_BASED_THINNING_TYPE") == 0);
    /*
     * handled elsewhere 
     */
    else
        return 0;

    return 1;
}

/** Reads std, estd, batch, qcflag, and instrument attributes from an open
 ** NetCDF file. For non-profile readers pass nprof = nobs, nz = 1.
 */
void obsatts_read(obsatts* oa, obssection* section, int ncid, const char* varname, size_t nobs, size_t nprof, size_t nz, int check_allmissing, int allow_char_batch_qc, observations* obs)
{
    int varid;
    int i;

    assert(nprof * nz == nobs);

    /*
     * addvars / subvars
     */
    for (i = 0; i < oa->naddvar; ++i) {
        obsaddvar* a = &oa->addvars[i];
        size_t k;

        ncw_inq_varid(ncid, a->varname, &varid);
        a->v = malloc(nobs * sizeof(double));
        oa_readvar_double(ncid, varid, nprof, nz, a->v);
        if (a->sign < 0) {
            for (k = 0; k < nobs; ++k)
                a->v[k] = -a->v[k];
        }
    }

    /*
     * std
     */
    varid = -1;
    if (oa->stdname != NULL)
        ncw_inq_varid(ncid, oa->stdname, &varid);
    else if (ncw_var_exists(ncid, "std"))
        ncw_inq_varid(ncid, "std", &varid);
    if (varid >= 0) {
        oa->std = malloc(nobs * sizeof(double));
        oa_readvar_double(ncid, varid, nprof, nz, oa->std);
        if (check_allmissing && oa_allmissing(nobs, oa->std)) {
            free(oa->std);
            oa->std = NULL;
        }
    }

    /*
     * estd
     */
    varid = -1;
    if (oa->estdname != NULL)
        ncw_inq_varid(ncid, oa->estdname, &varid);
    else if (ncw_var_exists(ncid, "error_std"))
        ncw_inq_varid(ncid, "error_std", &varid);
    if (varid >= 0) {
        oa->estd = malloc(nobs * sizeof(double));
        oa_readvar_double(ncid, varid, nprof, nz, oa->estd);
        if (check_allmissing && oa_allmissing(nobs, oa->estd)) {
            free(oa->estd);
            oa->estd = NULL;
        }
    }

    if (oa->estd == NULL) {
        ncw_inq_varid(ncid, varname, &varid);
        if (ncw_att_exists(ncid, varid, "error_std")) {
            ncw_check_attlen(ncid, varid, "error_std", 1);
            ncw_get_att_double(ncid, varid, "error_std", &oa->var_estd);
        }
    }

    /*
     * batch
     */
    varid = -1;
    if (oa->batchname != NULL)
        ncw_inq_varid(ncid, oa->batchname, &varid);
    else if (ncw_var_exists(ncid, "batch"))
        ncw_inq_varid(ncid, "batch", &varid);
    if (varid >= 0) {
        nc_type type;

        oa->batch = malloc(nobs * sizeof(int));
        ncw_inq_vartype(ncid, varid, &type);
        if (!allow_char_batch_qc || type != NC_CHAR) {
            oa_readvar_int(ncid, varid, nprof, nz, oa->batch);
        } else {
            int ndims;
            size_t dimlen[2];

            ncw_inq_vardims(ncid, varid, 2, &ndims, dimlen);
            if (ndims == 2) {
                char** data = alloc2d(dimlen[0], dimlen[1], 1);
                char* entry = malloc(dimlen[1] + 1);
                size_t p, k;

                assert(dimlen[0] == nprof);
                ncw_get_var_text(ncid, varid, data[0]);
                for (p = 0; p < nprof; ++p) {
                    int bval;

                    strncpy(entry, data[p], dimlen[1]);
                    entry[dimlen[1]] = 0;
                    if (!str2int(entry, &bval))
                        enkf_quit("%s: %s: could not convert \"%s\" to int", ncw_get_path(ncid), (oa->batchname != NULL) ? oa->batchname : "batch", entry);
                    for (k = 0; k < nz; ++k)
                        oa->batch[p * nz + k] = bval;
                }
                free(entry);
                free(data);
            } else
                enkf_quit("%s: %s: expected 2 dimensions for a batch variable of type NC_CHAR", ncw_get_path(ncid), (oa->batchname != NULL) ? oa->batchname : "batch");
        }
    }

    /*
     * qcflag
     */
    get_qcflags(section, &oa->nqcflagvars, &oa->qcflagvarnames, &oa->qcflagmasks);
    if (oa->nqcflagvars > 0) {
        oa->qcflag = alloc2d(oa->nqcflagvars, nobs, sizeof(int32_t));
        for (i = 0; i < oa->nqcflagvars; ++i) {
            nc_type type;

            ncw_inq_varid(ncid, oa->qcflagvarnames[i], &varid);
            ncw_inq_vartype(ncid, varid, &type);
            if (!allow_char_batch_qc || type != NC_CHAR) {
                oa_readvar_int(ncid, varid, nprof, nz, (int*) oa->qcflag[i]);
            } else {
                char* buf = calloc(nobs + 1, 1);
                size_t o;

                ncw_check_varsize(ncid, varid, nobs);
                ncw_get_var_text(ncid, varid, buf);
                for (o = 0; o < nobs; ++o)
                    oa->qcflag[i][o] = (int32_t) ((uint32_t) buf[o] - (uint32_t) '0');
                free(buf);
            }
        }
    }

    /*
     * instrument
     */
    if (strlen(oa->instrument) == 0) {
        if (oa->instattname != NULL)
            ncw_get_att_text_n(ncid, NC_GLOBAL, oa->instattname, oa->instrument, sizeof(oa->instrument));
        else if (!get_insttag(ncid, (char*) varname, oa->instrument))
            strncpy(oa->instrument, section->product, MAXSTRLEN - 1);
        if (strlen(oa->instrument) > 0 && oa->instprefix != NULL) {
            size_t len_p = strlen(oa->instprefix);
            size_t len_i = strlen(oa->instrument);

            assert(len_p + len_i < MAXSTRLEN);
            memmove(oa->instrument + len_p, oa->instrument, len_i + 1);
            memcpy(oa->instrument, oa->instprefix, len_p);
        }
    }
    oa->instid = st_add_ifabsent(obs->instruments, oa->instrument, -1);
}

/** Returns 1 if observation `idx` passes std/estd finiteness and all QC flag
 ** bitmasks, 0 otherwise.
 */
int obsatts_valid(const obsatts* oa, size_t idx)
{
    int i;

    if ((oa->std != NULL && isnan(oa->std[idx])) || (oa->estd != NULL && isnan(oa->estd[idx])))
        return 0;
    for (i = 0; i < oa->nqcflagvars; ++i)
        if (oa->qcflag[i][idx] < 0 || oa->qcflag[i][idx] > 31 || !((1u << oa->qcflag[i][idx]) & oa->qcflagmasks[i]))
            return 0;
    for (i = 0; i < oa->naddvar; ++i)
        if (isnan(oa->addvars[i].v[idx]))
            return 0;
    return 1;
}

/** Populates instrument, batch, and estd (applying std as lower bound) on `o`.
 */
void obsatts_apply(const obsatts* oa, size_t idx, int default_batch, observation* o)
{
    int i;

    for (i = 0; i < oa->naddvar; ++i)
        o->value += oa->addvars[i].v[idx];
    o->instrument = oa->instid;
    o->batch = (oa->batch == NULL) ? default_batch : oa->batch[idx];
    o->estd = (oa->estd == NULL) ? oa->var_estd : oa->estd[idx];
    if (oa->std != NULL && (isnan(o->estd) || oa->std[idx] > o->estd))
        o->estd = oa->std[idx];
}

/**
 */
void obsatts_free(obsatts* oa)
{
    if (oa->naddvar > 0) {
        int i;

        for (i = 0; i < oa->naddvar; ++i) {
            free(oa->addvars[i].varname);
            free(oa->addvars[i].v);
        }
        free(oa->addvars);
    }
    if (oa->std != NULL)
        free(oa->std);
    if (oa->estd != NULL)
        free(oa->estd);
    if (oa->batch != NULL)
        free(oa->batch);
    if (oa->nqcflagvars > 0) {
        free(oa->qcflagvarnames);
        free(oa->qcflagmasks);
        free(oa->qcflag);
    }
}
