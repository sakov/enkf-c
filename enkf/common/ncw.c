/******************************************************************************
 *
 * File:           ncw.c
 *  
 * Created         19/10/2000
 *  
 * Author:         Pavel Sakov
 *                 CSIRO Marine Research
 *  
 * Purpose:        Simple wrappers to netcdf library procedures for
 *                 better error messaging.
 *                 Some straightforward extensions to netcdf library 
 *                 procedures.
 *
 * Revisions:      PS 19/08/2003 Changed library name from "nc" to "ncw". 
 *                    Changed procedure prefixes from "NC_" to "ncw_".
 *
 *****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <errno.h>
#include "ncw.h"

const char ncw_version[] = "1.03";

/* This macro is substituted in error messages instead of the name of a
 * variable in cases when the name could not be found by the variable id.
 */
#define STR_UNKNOWN "???"

/* Used in ncw_find_vars()
 */
#define NALLOCATED_START 10

static void quit_def(char* format, ...);
static ncw_quit_fn quit = quit_def;

static void quit_def(char* format, ...)
{
    va_list args;

    fflush(stdout);

    fprintf(stderr, "\n  error: libncw: ");
    va_start(args, format);
    vfprintf(stderr, format, args);
    va_end(args);
    exit(1);
}

static void _ncw_inq_varname(const char fname[], int ncid, int varid, char varname[])
{
    int status;

    if (varid == -1)
        strcpy(varname, "NC_GLOBAL");
    else if ((status = nc_inq_varname(ncid, varid, varname)) != NC_NOERR)
        quit("\"%s\": nc_inq_varname(): failed for varid = %d: %s\n", fname, varid, nc_strerror(status));
}

static void _ncw_inq_varid(const char fname[], int ncid, const char varname[], int* varid)
{
    int status;

    if (strcmp(varname, "NC_GLOBAL") == 0)
        *varid = -1;
    else if ((status = nc_inq_varid(ncid, varname, varid)) != NC_NOERR)
        quit("\"%s\": nc_inq_varid(): failed for varid = %d: %s\n", fname, varid, nc_strerror(status));
}

#define STRBUFSIZE 1024

/* Prints array of integers to a string. E.g., {1,2,5} will be printed as
 * "(1,2,5)".
 */
static char* int2str(int n, const int v[])
{
    int i;
    char all[STRBUFSIZE] = "(";
    char next[STRBUFSIZE];

    for (i = 0; i < n; ++i) {
        if (i < n - 1)
            snprintf(next, STRBUFSIZE, "%d,", v[i]);
        else
            snprintf(next, STRBUFSIZE, "%d", v[i]);
        assert(strlen(all) + strlen(next) + 2 < STRBUFSIZE);
        strcat(all, next);
    }
    assert(strlen(all) + 2 < STRBUFSIZE);
    strcat(all, ")");

    return strdup(all);
}

static char* uint2str(int n, const unsigned int v[])
{
    int i;
    char all[STRBUFSIZE] = "(";
    char next[STRBUFSIZE];

    for (i = 0; i < n; ++i) {
        if (i < n - 1)
            snprintf(next, STRBUFSIZE, "%d,", v[i]);
        else
            snprintf(next, STRBUFSIZE, "%d", v[i]);
        assert(strlen(all) + strlen(next) + 2 < STRBUFSIZE);
        strcat(all, next);
    }
    assert(strlen(all) + 2 < STRBUFSIZE);
    strcat(all, ")");

    return strdup(all);
}

static char* double2str(int n, const double v[])
{
    int i;
    char all[STRBUFSIZE] = "(";
    char next[STRBUFSIZE];

    for (i = 0; i < n && i < 3; ++i) {
        if (i < n - 1)
            snprintf(next, STRBUFSIZE, "%.4g,", v[i]);
        else
            snprintf(next, STRBUFSIZE, "%.4g", v[i]);
        assert(strlen(all) + strlen(next) + 2 < STRBUFSIZE);
        strcat(all, next);
    }
    if (n > 3) {
        assert(strlen(all) + strlen(",...") + 1 < STRBUFSIZE);
        strcat(all, ",...");
    }
    assert(strlen(all) + strlen(")") + 1 < STRBUFSIZE);
    strcat(all, ")");

    return strdup(all);
}

static char* float2str(int n, const float v[])
{
    int i;
    char all[STRBUFSIZE] = "(";
    char next[STRBUFSIZE];

    for (i = 0; i < n && i < 3; ++i) {
        if (i < n - 1)
            snprintf(next, STRBUFSIZE, "%.4g,", v[i]);
        else
            snprintf(next, STRBUFSIZE, "%.4g", v[i]);
        assert(strlen(all) + strlen(next) + 2 < STRBUFSIZE);
        strcat(all, next);
    }
    if (n > 3) {
        assert(strlen(all) + strlen(",...") + 1 < STRBUFSIZE);
        strcat(all, ",...");
    }
    assert(strlen(all) + strlen(")") + 1 < STRBUFSIZE);
    strcat(all, ")");

    return strdup(all);
}

/* *INDENT-OFF* */
static struct nctype2str {
    nc_type type;
    char* str;
} nctypes2str[] = {
    {-1, "UNKNOWN"},
    {NC_BYTE, "NC_BYTE"},
    {NC_CHAR, "NC_CHAR"},
    {NC_SHORT, "NC_SHORT"},
    {NC_INT, "NC_INT"},
    {NC_FLOAT, "NC_FLOAT"},
    {NC_DOUBLE, "NC_DOUBLE"},
};
/* *INDENT-ON* */

void ncw_create(const char fname[], int mode, int* ncid)
{
    int status = nc_create(fname, mode, ncid);

    if (status != NC_NOERR)
        quit("\"%s\": nc_create(): failed: %s\n", fname, nc_strerror(status));
}

void ncw_open(const char fname[], int mode, int* ncid)
{
    int status = nc_open(fname, mode, ncid);

    if (status != NC_NOERR)
        quit("\"%s\": nc_open(): failed: %s\n", fname, nc_strerror(status));
}

void ncw_redef(const char fname[], int ncid)
{
    int status = nc_redef(ncid);

    if (status != NC_NOERR)
        quit("\"%s\": ncredef(): failed: %s\n", fname, nc_strerror(status));
}

void ncw_enddef(const char fname[], int ncid)
{
    int status = nc_enddef(ncid);

    if (status != NC_NOERR)
        quit("\"%s\": ncenddef(): failed: %s\n", fname, nc_strerror(status));
}

void ncw_sync(const char fname[], int ncid)
{
    int status;

    nc_enddef(ncid);

    status = ncsync(ncid);

    if (status != NC_NOERR)
        quit("\"%s\": ncsync(): failed: %s\n", fname, nc_strerror(status));
}

void ncw_close(const char fname[], int ncid)
{
    int status = nc_close(ncid);

    if (status != NC_NOERR)
        quit("\"%s\": nc_close(): failed: %s\n", fname, nc_strerror(status));
}

void ncw_inq(const char fname[], int ncid, int* ndims, int* nvars, int* natts, int* unlimdimid)
{
    int status = nc_inq(ncid, ndims, nvars, natts, unlimdimid);

    if (status != NC_NOERR)
        quit("\"%s\": nc_inq(): failed: %s\n", fname, nc_strerror(status));
}

void ncw_inq_ndims(const char fname[], int ncid, int* ndims)
{
    int status = nc_inq_ndims(ncid, ndims);

    if (status != NC_NOERR)
        quit("\"%s\": nc_inq_ndims(): failed: %s\n", fname, nc_strerror(status));
}

void ncw_inq_nvars(const char fname[], int ncid, int* nvars)
{
    int status = nc_inq_nvars(ncid, nvars);

    if (status != NC_NOERR)
        quit("\"%s\": nc_inq_nvars(): failed: %s\n", fname, nc_strerror(status));
}

void ncw_inq_natts(const char fname[], int ncid, int* natts)
{
    int status = nc_inq_natts(ncid, natts);

    if (status != NC_NOERR)
        quit("\"%s\": nc_inq_natts(): failed: %s\n", fname, nc_strerror(status));
}

void ncw_inq_unlimdimid(const char fname[], int ncid, int* unlimdimid)
{
    int status = nc_inq_unlimdim(ncid, unlimdimid);

    if (status != NC_NOERR)
        quit("\"%s\": nc_inq_unlimdim(): failed: %s\n", fname, nc_strerror(status));
}

void ncw_def_dim(const char fname[], int ncid, const char dimname[], size_t len, int* dimid)
{
    int status = nc_def_dim(ncid, dimname, len, dimid);

    if (status != NC_NOERR)
        quit("\"%s\": nc_def_dim(): failed for dimname = \"%s\", dimlen = %d: %s\n", fname, dimname, len, nc_strerror(status));
}

void ncw_inq_dimid(const char fname[], int ncid, const char dimname[], int* dimid)
{
    int status = nc_inq_dimid(ncid, dimname, dimid);

    if (status != NC_NOERR)
        quit("\"%s\": nc_inq_dimid(): failed for dimname = \"%s\": %s\n", fname, dimname, nc_strerror(status));
}

void ncw_inq_dim(const char fname[], int ncid, int dimid, char dimname[], size_t* len)
{
    int status = nc_inq_dim(ncid, dimid, dimname, len);

    if (status != NC_NOERR)
        quit("\"%s\": nc_inq_dim(): failed for dimid = %d: %s\n", fname, dimid, nc_strerror(status));
}

void ncw_inq_dimname(const char fname[], int ncid, int dimid, char dimname[])
{
    int status = nc_inq_dimname(ncid, dimid, dimname);

    if (status != NC_NOERR)
        quit("\"%s\": nc_inq_dimname(): failed for dimid = %d: %s\n", fname, dimid, nc_strerror(status));
}

void ncw_inq_dimlen(const char fname[], int ncid, int dimid, size_t* len)
{
    int status = nc_inq_dimlen(ncid, dimid, len);

    if (status != NC_NOERR) {
        char dimname[NC_MAX_NAME] = STR_UNKNOWN;

        ncw_inq_dimname(fname, ncid, dimid, dimname);
        quit("\"%s\": nc_inq_dimlen(): failed for dimid = %d (dimname = \"%s\"): %s\n", fname, dimid, dimname, nc_strerror(status));
    }
}

void ncw_rename_dim(const char fname[], int ncid, const char oldname[], const char newname[])
{
    int dimid;
    int status;

    ncw_inq_dimid(fname, ncid, oldname, &dimid);
    status = nc_rename_dim(ncid, dimid, newname);

    if (status != NC_NOERR)
        quit("\"%s\": nc_rename_dim(): failed for dimid = %d (oldname = \"%s\", newname = \"%s\"): %s\n", dimid, oldname, newname, nc_strerror(status));
}

void ncw_def_var(const char fname[], int ncid, const char varname[], nc_type xtype, int ndims, const int dimids[], int* varid)
{
    int status = nc_def_var(ncid, varname, xtype, ndims, dimids, varid);

    if (status != NC_NOERR)
        quit("\"%s\": nc_def_var(): failed for varname = \"%s\", vartype = %s, ndims = %d, dimids = %s: %s\n", fname, varname, ncw_nctype2str(xtype), ndims, uint2str(ndims, (const unsigned int*) dimids), nc_strerror(status));
}

void ncw_inq_varid(const char* fname, int ncid, const char varname[], int* varid)
{
    int status = nc_inq_varid(ncid, varname, varid);

    if (status != NC_NOERR)
        quit("\"%s\": nc_inq_varid(): failed for varname = \"%s\": %s\n", fname, varname, nc_strerror(status));
}

void ncw_inq_var(const char fname[], int ncid, int varid, char varname[], nc_type* xtype, int* ndims, int dimids[], int* natts)
{
    int status = nc_inq_var(ncid, varid, varname, xtype, ndims, dimids, natts);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        ncw_inq_varname(fname, ncid, varid, varname);
        quit("\"%s\": nc_inq_var(): failed for varid = %d (varname = \"%s\"): %s\n", fname, varid, varname, nc_strerror(status));
    }
}

void ncw_inq_varname(const char fname[], int ncid, int varid, char varname[])
{
    int status = nc_inq_varname(ncid, varid, varname);

    if (status != NC_NOERR)
        quit("\"%s\": nc_inq_varname(): failed for varid = %d: %s\n", fname, varid, nc_strerror(status));
}

void ncw_inq_vartype(const char fname[], int ncid, int varid, nc_type* xtype)
{
    int status = nc_inq_vartype(ncid, varid, xtype);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        ncw_inq_varname(fname, ncid, varid, varname);
        quit("\"%s\": nc_inq_vartype(): failed for varid = %d (varname = %s): %s\n", fname, varid, varname, nc_strerror(status));
    }
}

void ncw_inq_varndims(const char fname[], int ncid, int varid, int* ndims)
{
    int status = nc_inq_varndims(ncid, varid, ndims);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = "STR_UNKNOWN";

        ncw_inq_varname(fname, ncid, varid, varname);
        quit("\"%s\": nc_inq_varndims(): failed for varid = %d (varname = \"%s\"): %s\n", fname, varid, varname, nc_strerror(status));
    }
}

void ncw_inq_vardimid(const char fname[], int ncid, int varid, int dimids[])
{
    int status = nc_inq_vardimid(ncid, varid, dimids);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        ncw_inq_varname(fname, ncid, varid, varname);
        quit("\"%s\": ncw_inq_vardimid(): failed for varid = %d (varname = \"%s\"): %s\n", fname, varid, varname, nc_strerror(status));
    }
}

void ncw_inq_varnatts(const char fname[], int ncid, int varid, int* natts)
{
    int status = nc_inq_varnatts(ncid, varid, natts);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = "STR_UNKNOWN";

        _ncw_inq_varname(fname, ncid, varid, varname);
        quit("\"%s\": nc_inq_varnatts(): failed for varid = %d (varnatts = \"%s\"): %s\n", fname, varid, varname, nc_strerror(status));
    }
}

void ncw_rename_var(const char fname[], int ncid, const char oldname[], const char newname[])
{
    int varid;
    int status;

    ncw_inq_varid(fname, ncid, oldname, &varid);
    status = nc_rename_var(ncid, varid, newname);

    if (status != NC_NOERR)
        quit("\"%s\": nc_rename_var(): failed for varid = %d (oldname = \"%s\", newname = \"%s\"): %s\n", fname, varid, oldname, newname, nc_strerror(status));
}

void ncw_def_var_deflate(const char fname[], int ncid, int varid, int shuffle, int deflate, int deflate_level)
{
    int status = nc_def_var_deflate(ncid, varid, shuffle, deflate, deflate_level);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = "STR_UNKNOWN";

        _ncw_inq_varname(fname, ncid, varid, varname);
        quit("\"%s\": ncw_def_var_deflate(): failed for varid = %d (varname = \"%s\"): %s\n", fname, varid, varname, nc_strerror(status));
    }
}

void ncw_put_var_text(const char fname[], int ncid, int varid, const char v[])
{
    int status = nc_put_var_text(ncid, varid, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        ncw_inq_varname(fname, ncid, varid, varname);
        quit("\"%s\": nc_put_var_text(): failed for varid = %d (varname = \"%s\"): %s\n", fname, varid, varname, nc_strerror(status));
    }
}

void ncw_put_var_short(const char fname[], int ncid, int varid, const short v[])
{
    int status = nc_put_var_short(ncid, varid, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        ncw_inq_varname(fname, ncid, varid, varname);
        quit("\"%s\": nc_put_var_short(): failed for varid = %d (varname = \"%s\"): %s\n", fname, varid, varname, nc_strerror(status));
    }
}

void ncw_put_var_ushort(const char fname[], int ncid, int varid, const unsigned short v[])
{
    int status = nc_put_var_ushort(ncid, varid, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        ncw_inq_varname(fname, ncid, varid, varname);
        quit("\"%s\": nc_put_var_ushort(): failed for varid = %d (varname = \"%s\"): %s\n", fname, varid, varname, nc_strerror(status));
    }
}

void ncw_put_var_int(const char fname[], int ncid, int varid, const int v[])
{
    int status = nc_put_var_int(ncid, varid, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        ncw_inq_varname(fname, ncid, varid, varname);
        quit("\"%s\": nc_put_var_int(): failed for varid = %d (varname = \"%s\"): %s\n", fname, varid, varname, nc_strerror(status));
    }
}

void ncw_put_var_uint(const char fname[], int ncid, int varid, const unsigned int v[])
{
    int status = nc_put_var_uint(ncid, varid, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        ncw_inq_varname(fname, ncid, varid, varname);
        quit("\"%s\": nc_put_var_uint(): failed for varid = %d (varname = \"%s\"): %s\n", fname, varid, varname, nc_strerror(status));
    }
}

void ncw_put_var_float(const char fname[], int ncid, int varid, const float v[])
{
    int status = nc_put_var_float(ncid, varid, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        ncw_inq_varname(fname, ncid, varid, varname);
        quit("\"%s\": nc_put_var_float(): failed for varid = %d (varname = \"%s\"): %s\n", fname, varid, varname, nc_strerror(status));
    }
}

void ncw_put_var_double(const char fname[], int ncid, int varid, const double v[])
{
    int status = nc_put_var_double(ncid, varid, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        ncw_inq_varname(fname, ncid, varid, varname);
        quit("\"%s\": nc_put_var_double(): failed for varid = %d (varname = \"%s\"): %s\n", fname, varid, varname, nc_strerror(status));
    }
}

void ncw_get_var_text(const char fname[], int ncid, int varid, char v[])
{
    int status = nc_get_var_text(ncid, varid, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        ncw_inq_varname(fname, ncid, varid, varname);
        quit("\"%s\": nc_get_var_text(): failed for varid = %d (varname = \"%s\"): %s\n", fname, varid, varname, nc_strerror(status));
    }
}

void ncw_get_var_schar(const char fname[], int ncid, int varid, signed char v[])
{
    int status = nc_get_var_schar(ncid, varid, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        ncw_inq_varname(fname, ncid, varid, varname);
        quit("\"%s\": nc_get_var_schar(): failed for varid = %d (varname = \"%s\"): %s\n", fname, varid, varname, nc_strerror(status));
    }
}

void ncw_get_var_short(const char fname[], int ncid, int varid, short int v[])
{
    int status = nc_get_var_short(ncid, varid, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        ncw_inq_varname(fname, ncid, varid, varname);
        quit("\"%s\": nc_get_var_short(): failed for varid = %d (varname = \"%s\"): %s\n", fname, varid, varname, nc_strerror(status));
    }
}

void ncw_get_var_ushort(const char fname[], int ncid, int varid, unsigned short int v[])
{
    int status = nc_get_var_ushort(ncid, varid, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        ncw_inq_varname(fname, ncid, varid, varname);
        quit("\"%s\": nc_get_var_ushort(): failed for varid = %d (varname = \"%s\"): %s\n", fname, varid, varname, nc_strerror(status));
    }
}

void ncw_get_var_int(const char fname[], int ncid, int varid, int v[])
{
    int status = nc_get_var_int(ncid, varid, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        ncw_inq_varname(fname, ncid, varid, varname);
        quit("\"%s\": nc_get_var_int(): failed for varid = %d (varname = \"%s\"): %s\n", fname, varid, varname, nc_strerror(status));
    }
}

void ncw_get_var_uint(const char fname[], int ncid, int varid, unsigned int v[])
{
    int status = nc_get_var_uint(ncid, varid, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        ncw_inq_varname(fname, ncid, varid, varname);
        quit("\"%s\": nc_get_var_uint(): failed for varid = %d (varname = \"%s\"): %s\n", fname, varid, varname, nc_strerror(status));
    }
}

void ncw_get_var_float(const char fname[], int ncid, int varid, float v[])
{
    int status = nc_get_var_float(ncid, varid, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        ncw_inq_varname(fname, ncid, varid, varname);
        quit("\"%s\": nc_get_var_float(): failed for varid = %d (varname = \"%s\"): %s\n", fname, varid, varname, nc_strerror(status));
    }
}

void ncw_get_var_double(const char fname[], int ncid, int varid, double v[])
{
    int status = nc_get_var_double(ncid, varid, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        ncw_inq_varname(fname, ncid, varid, varname);
        quit("\"%s\": nc_get_var_double(): failed for varid = %d (varname = \"%s\"): %s\n", fname, varid, varname, nc_strerror(status));
    }
}

void ncw_get_var1_double(const char fname[], int ncid, int varid, const size_t len[], double* in)
{
    int status = nc_get_var1_double(ncid, varid, len, in);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        ncw_inq_varname(fname, ncid, varid, varname);
        quit("\"%s\": nc_get_var1_double(): failed for varid = %d (varname = \"%s\"): %s\n", fname, varid, varname, nc_strerror(status));
    }
}

void ncw_put_vara_text(const char fname[], int ncid, int varid, const size_t start[], const size_t count[], const char v[])
{
    int status = nc_put_vara_text(ncid, varid, start, count, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;
        int ndims = 0;

        ncw_inq_varname(fname, ncid, varid, varname);
        ncw_inq_varndims(fname, ncid, varid, &ndims);
        quit("\"%s\": ncw_put_vara_text(): failed for varid = %d (varname = \"%s\"), start = %s, count = %s: %s\n", fname, varid, varname, uint2str(ndims, (unsigned int*) start), uint2str(ndims, (unsigned int*) count), nc_strerror(status));
    }
}

void ncw_put_vara_short(const char fname[], int ncid, int varid, const size_t start[], const size_t count[], const short int v[])
{
    int status = nc_put_vara_short(ncid, varid, start, count, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;
        int ndims = 0;

        ncw_inq_varname(fname, ncid, varid, varname);
        ncw_inq_varndims(fname, ncid, varid, &ndims);
        quit("\"%s\": ncw_put_vara_short(): failed for varid = %d (varname = \"%s\"), start = %s, count = %s: %s\n", fname, varid, varname, uint2str(ndims, (unsigned int*) start), uint2str(ndims, (unsigned int*) count), nc_strerror(status));
    }
}

void ncw_put_vara_ushort(const char fname[], int ncid, int varid, const size_t start[], const size_t count[], const unsigned short int v[])
{
    int status = nc_put_vara_ushort(ncid, varid, start, count, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;
        int ndims = 0;

        ncw_inq_varname(fname, ncid, varid, varname);
        ncw_inq_varndims(fname, ncid, varid, &ndims);
        quit("\"%s\": ncw_put_vara_ushort(): failed for varid = %d (varname = \"%s\"), start = %s, count = %s: %s\n", fname, varid, varname, uint2str(ndims, (unsigned int*) start), uint2str(ndims, (unsigned int*) count), nc_strerror(status));
    }
}

void ncw_put_vara_int(const char fname[], int ncid, int varid, const size_t start[], const size_t count[], const int v[])
{
    int status = nc_put_vara_int(ncid, varid, start, count, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;
        int ndims = 0;

        ncw_inq_varname(fname, ncid, varid, varname);
        ncw_inq_varndims(fname, ncid, varid, &ndims);
        quit("\"%s\": ncw_put_vara_int(): failed for varid = %d (varname = \"%s\"), start = %s, count = %s: %s\n", fname, varid, varname, uint2str(ndims, (unsigned int*) start), uint2str(ndims, (unsigned int*) count), nc_strerror(status));
    }
}

void ncw_put_vara_uint(const char fname[], int ncid, int varid, const size_t start[], const size_t count[], const unsigned int v[])
{
    int status = nc_put_vara_uint(ncid, varid, start, count, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;
        int ndims = 0;

        ncw_inq_varname(fname, ncid, varid, varname);
        ncw_inq_varndims(fname, ncid, varid, &ndims);
        quit("\"%s\": ncw_put_vara_uint(): failed for varid = %d (varname = \"%s\"), start = %s, count = %s: %s\n", fname, varid, varname, uint2str(ndims, (unsigned int*) start), uint2str(ndims, (unsigned int*) count), nc_strerror(status));
    }
}

void ncw_put_vara_float(const char fname[], int ncid, int varid, const size_t start[], const size_t count[], const float v[])
{
    int status = nc_put_vara_float(ncid, varid, start, count, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;
        int ndims = 0;

        ncw_inq_varname(fname, ncid, varid, varname);
        ncw_inq_varndims(fname, ncid, varid, &ndims);
        quit("\"%s\": ncw_put_vara_float(): failed for varid = %d (varname = \"%s\"), start = %s, count = %s: %s\n", fname, varid, varname, uint2str(ndims, (unsigned int*) start), uint2str(ndims, (unsigned int*) count), nc_strerror(status));
    }
}

void ncw_put_vara_double(const char fname[], int ncid, int varid, const size_t start[], const size_t count[], const double v[])
{
    int status = nc_put_vara_double(ncid, varid, start, count, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;
        int ndims = 0;

        ncw_inq_varname(fname, ncid, varid, varname);
        ncw_inq_varndims(fname, ncid, varid, &ndims);
        quit("\"%s\": nc_put_vara_double(): failed for varid = %d (varname = \"%s\"), start = %s, count = %s: %s\n", fname, varid, varname, uint2str(ndims, (unsigned int*) start), uint2str(ndims, (unsigned int*) count), nc_strerror(status));
    }
}

void ncw_get_vara_text(const char fname[], int ncid, int varid, const size_t start[], const size_t count[], char v[])
{
    int status = nc_get_vara_text(ncid, varid, start, count, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        ncw_inq_varname(fname, ncid, varid, varname);
        quit("\"%s\": nc_get_vara_text(): failed for varid = %d (varname = \"%s\"): %s\n", fname, varid, varname, nc_strerror(status));
    }
}

void ncw_get_vara_short(const char fname[], int ncid, int varid, const size_t start[], const size_t count[], short int v[])
{
    int status = nc_get_vara_short(ncid, varid, start, count, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        ncw_inq_varname(fname, ncid, varid, varname);
        quit("\"%s\": nc_get_vara_short(): failed for varid = %d (varname = \"%s\"): %s\n", fname, varid, varname, nc_strerror(status));
    }
}

void ncw_get_vara_int(const char fname[], int ncid, int varid, const size_t start[], const size_t count[], int v[])
{
    int status = nc_get_vara_int(ncid, varid, start, count, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        ncw_inq_varname(fname, ncid, varid, varname);
        quit("\"%s\": nc_get_vara_int(): failed for varid = %d (varname = \"%s\"): %s\n", fname, varid, varname, nc_strerror(status));
    }
}

void ncw_get_vara_float(const char fname[], int ncid, int varid, const size_t start[], const size_t count[], float v[])
{
    int status = nc_get_vara_float(ncid, varid, start, count, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        ncw_inq_varname(fname, ncid, varid, varname);
        quit("\"%s\": nc_get_vara_float(): failed for varid = %d (varname = \"%s\"): %s\n", fname, varid, varname, nc_strerror(status));
    }
}

void ncw_get_vara_double(const char fname[], int ncid, int varid, const size_t start[], const size_t count[], double v[])
{
    int status = nc_get_vara_double(ncid, varid, start, count, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        ncw_inq_varname(fname, ncid, varid, varname);
        quit("\"%s\": nc_get_vara_double(): failed for varid = %d (varname = \"%s\"): %s\n", fname, varid, varname, nc_strerror(status));
    }
}

void ncw_put_att_text(const char fname[], int ncid, int varid, const char attname[], const char v[])
{
    int status = nc_put_att_text(ncid, varid, attname, strlen(v), v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        _ncw_inq_varname(fname, ncid, varid, varname);
        quit("\"%s\": nc_put_att_text(): failed for varid = %d (varname = \"%s\"), attname = \"%s\", attvalue = \"%s\": %s\n", fname, varid, varname, attname, v, nc_strerror(status));
    }
}

void ncw_put_att_int(const char fname[], int ncid, int varid, const char attname[], size_t len, const int v[])
{
    int status = nc_put_att_int(ncid, varid, attname, NC_INT, len, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        _ncw_inq_varname(fname, ncid, varid, varname);
        quit("\"%s\": nc_put_att_int(): failed for varid = %d (varname = \"%s\"), attname = \"%s\", attvalue(s) = \"%s\": %s\n", fname, varid, varname, attname, int2str(len, v), nc_strerror(status));
    }
}

void ncw_put_att_float(const char fname[], int ncid, int varid, const char attname[], size_t len, const float v[])
{
    int status = nc_put_att_float(ncid, varid, attname, NC_FLOAT, len, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = "STR_UNKNOWN";

        _ncw_inq_varname(fname, ncid, varid, varname);
        quit("\"%s\": nc_put_att_float(): failed for varid = %d (varname = \"%s\"), attname = \"%s\", attvalue = \"%s\": %s\n", fname, varid, varname, attname, float2str(len, v), nc_strerror(status));
    }

    if (strcmp(attname, "_FillValue") == 0 && varid != NC_GLOBAL) {
        nc_type xtype;

        ncw_inq_vartype(fname, ncid, varid, &xtype);
        if (xtype != NC_FLOAT) {
            char varname[NC_MAX_NAME] = "STR_UNKNOWN";

            ncw_inq_varname(fname, ncid, varid, varname);

            quit("\"%s\": ncw_put_att_float(): failed for varid = %d (varname = \"%s\"): attype = NC_FLOAT differs from vartype = %s for attname = \"_FillValue\"\n", fname, varid, varname, ncw_nctype2str(xtype));
        }
    }
}

void ncw_put_att_double(const char fname[], int ncid, int varid, const char attname[], size_t len, const double v[])
{
    int status = nc_put_att_double(ncid, varid, attname, NC_DOUBLE, len, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        _ncw_inq_varname(fname, ncid, varid, varname);
        quit("\"%s\": nc_put_att_double(): failed for varid = %d (varname = \"%s\"), attname = \"%s\", attvalue(s) = \"%s\": %s\n", fname, varid, varname, attname, double2str(len, v), nc_strerror(status));
    }

    if (strcmp(attname, "_FillValue") == 0 && varid != NC_GLOBAL) {
        nc_type xtype;

        ncw_inq_vartype(fname, ncid, varid, &xtype);
        if (xtype != NC_DOUBLE) {
            char varname[NC_MAX_NAME] = "STR_UNKNOWN";

            ncw_inq_varname(fname, ncid, varid, varname);

            quit("\"%s\": ncw_put_att_float(): fatal inconsistency for varid = %d (varname = \"%s\"): attype = NC_DOUBLE differs from vartype = %s for attname = \"_FillValue\"\n", fname, varid, varname, ncw_nctype2str(xtype));
        }
    }
}

void ncw_inq_attname(const char fname[], int ncid, int varid, int attrid, char attname[])
{
    int status = nc_inq_attname(ncid, varid, attrid, attname);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        _ncw_inq_varname(fname, ncid, varid, varname);
        quit("\"%s\": nc_inq_attname(): failed for varid = %d (varname = \"%s\"): %s\n", fname, varid, varname, nc_strerror(status));
    }
}

void ncw_inq_att(const char fname[], int ncid, int varid, const char attname[], nc_type* xtype, size_t* len)
{
    int status = nc_inq_att(ncid, varid, attname, xtype, len);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        _ncw_inq_varname(fname, ncid, varid, varname);
        quit("\"%s\": nc_inq_att(): failed for varid = %d (varname = \"%s\"), attname = \"%s\", atttype = %d: %s\n", fname, varid, varname, attname, ncw_nctype2str(*xtype), nc_strerror(status));
    }
}

void ncw_inq_attlen(const char fname[], int ncid, int varid, const char attname[], size_t* len)
{
    int status = nc_inq_attlen(ncid, varid, attname, len);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        _ncw_inq_varname(fname, ncid, varid, varname);
        quit("\"%s\": nc_inq_attlen(): failed for varid = %d (varname = \"%s\"), attname = \"%s\": %s\n", fname, varid, varname, attname, nc_strerror(status));
    }
}

void ncw_copy_att(const char fname_src[], int ncid_src, int varid_src, const char attname[], const char fname_dst[], int ncid_dst, int varid_dst)
{
    int status = nc_copy_att(ncid_src, varid_src, attname, ncid_dst, varid_dst);

    if (status != NC_NOERR) {
        char varname_src[NC_MAX_NAME] = STR_UNKNOWN;
        char varname_dst[NC_MAX_NAME] = STR_UNKNOWN;

        _ncw_inq_varname(fname_src, ncid_src, varid_src, varname_src);
        _ncw_inq_varname(fname_dst, ncid_dst, varid_dst, varname_dst);

        quit("\"%s\": -> %s:  nc_copy_att(): failed for varid_src = %d (varname = \"%s\"), attname = \"%s\", varid_dst = %d, (varname = \"%s\"): %s\n", fname_src, fname_dst, varid_src, varname_src, attname, varid_dst, varname_dst, nc_strerror(status));
    }
}

void ncw_rename_att(const char fname[], int ncid, const char varname[], const char oldname[], const char newname[])
{
    int varid;
    int status;

    _ncw_inq_varid(fname, ncid, varname, &varid);
    status = nc_rename_att(ncid, varid, oldname, newname);

    if (status != NC_NOERR)
        quit("\"%s\": nc_rename_att(): failed for varid = %d, oldname = \"%s\", newname = \"%s\": %s\n", varid, oldname, newname, nc_strerror(status));
}

void ncw_del_att(const char fname[], int ncid, int varid, const char attname[])
{
    int status = nc_del_att(ncid, varid, attname);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        _ncw_inq_varname(fname, ncid, varid, varname);
        quit("\"%s\": nc_del_att(): failed for varid = %d (varname = \"%s\"): %s\n", fname, varid, varname, nc_strerror(status));
    }
}

void ncw_get_att_text(const char fname[], int ncid, int varid, const char attname[], char v[])
{
    nc_type xtype;
    size_t len;
    int status;

    ncw_inq_att(fname, ncid, varid, attname, &xtype, &len);

    if (xtype != NC_CHAR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        _ncw_inq_varname(fname, ncid, varid, varname);
        quit("ncw_get_att_text(): failed for varid = %d (varname = \"%s\"), attname = \"%s\": the attribute is of %s type\n", fname, varid, varname, attname, ncw_nctype2str(xtype));
    }

    status = nc_get_att_text(ncid, varid, attname, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        _ncw_inq_varname(fname, ncid, varid, varname);
        quit("\"%s\": nc_get_att_text(): failed for varid = %d (varname = \"%s\"), attname = \"%s\": %s\n", fname, varid, varname, attname, nc_strerror(status));
    }
}

void ncw_get_att_short(const char fname[], int ncid, int varid, const char attname[], short int v[])
{
    nc_type xtype;
    size_t len;
    int status;

    ncw_inq_att(fname, ncid, varid, attname, &xtype, &len);

    if (xtype != NC_SHORT) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        _ncw_inq_varname(fname, ncid, varid, varname);
        quit("ncw_get_att_short(): failed for varid = %d (varname = \"%s\"), attname = \"%s\": the attribute is of %s type\n", fname, varid, varname, attname, ncw_nctype2str(xtype));
    }

    status = nc_get_att_short(ncid, varid, attname, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        _ncw_inq_varname(fname, ncid, varid, varname);
        quit("\"%s\": nc_get_att_short(): failed for varid = %d (varname = \"%s\"), attname = \"%s\": %s\n", fname, varid, varname, attname, nc_strerror(status));
    }
}

void ncw_get_att_int(const char fname[], int ncid, int varid, const char attname[], int v[])
{
    nc_type xtype;
    size_t len;
    int status;

    ncw_inq_att(fname, ncid, varid, attname, &xtype, &len);

    if (xtype != NC_INT && xtype != NC_BYTE && xtype != NC_SHORT) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        _ncw_inq_varname(fname, ncid, varid, varname);
        quit("ncw_get_att_int(): failed for varid = %d (varname = \"%s\"), attname = \"%s\": the attribute is of %s type\n", fname, varid, varname, attname, ncw_nctype2str(xtype));
    }

    status = nc_get_att_int(ncid, varid, attname, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        _ncw_inq_varname(fname, ncid, varid, varname);
        quit("\"%s\": nc_get_att_int(): failed for varid = %d (varname = \"%s\"), attname = \"%s\": %s\n", fname, varid, varname, attname, nc_strerror(status));
    }
}

void ncw_get_att_float(const char fname[], int ncid, int varid, const char attname[], float v[])
{
    nc_type xtype;
    size_t len;
    int status;

    ncw_inq_att(fname, ncid, varid, attname, &xtype, &len);

    if (xtype != NC_DOUBLE && xtype != NC_FLOAT) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        _ncw_inq_varname(fname, ncid, varid, varname);
        quit("ncw_get_att_double(): failed for varid = %d (varname = \"%s\"), attname = \"%s\": the attribute is of %s type\n", fname, varid, varname, attname, ncw_nctype2str(xtype));
    }

    status = nc_get_att_float(ncid, varid, attname, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        _ncw_inq_varname(fname, ncid, varid, varname);
        quit("\"%s\": ncw_get_att_float(): failed for varid = %d (varname = \"%s\"), attname = \"%s\": %s\n", fname, varid, varname, attname, nc_strerror(status));
    }
}

void ncw_get_att_double(const char fname[], int ncid, int varid, const char attname[], double v[])
{
    nc_type xtype;
    size_t len;
    int status;

    ncw_inq_att(fname, ncid, varid, attname, &xtype, &len);

    if (xtype != NC_DOUBLE && xtype != NC_FLOAT) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        _ncw_inq_varname(fname, ncid, varid, varname);
        quit("ncw_get_att_double(): failed for varid = %d (varname = \"%s\"), attname = \"%s\": the attribute is of %s type\n", fname, varid, varname, attname, ncw_nctype2str(xtype));
    }

    status = nc_get_att_double(ncid, varid, attname, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        _ncw_inq_varname(fname, ncid, varid, varname);
        quit("\"%s\": ncw_get_att_double(): failed for varid = %d (varname = \"%s\"), attname = \"%s\": %s\n", fname, varid, varname, attname, nc_strerror(status));
    }
}

/*
 * The following procedures do not have direct analogues in the netcdf library.
 */

int ncw_inq_nrecords(const char fname[], int ncid)
{
    int unlimdimid;
    size_t nrecords;

    ncw_inq_unlimdimid(fname, ncid, &unlimdimid);
    if (unlimdimid < 0)
        return 0;
    ncw_inq_dimlen(fname, ncid, unlimdimid, &nrecords);

    return nrecords;
}

const char* ncw_nctype2str(nc_type type)
{
    size_t i;

    for (i = 1; i < sizeof(nctypes2str) / sizeof(struct nctype2str); ++i) {
        if (type == nctypes2str[i].type)
            return nctypes2str[i].str;
    }

    return nctypes2str[0].str;
}

size_t ncw_sizeof(nc_type type)
{
    switch (type) {
    case NC_BYTE:
    case NC_CHAR:
        return sizeof(char);
    case NC_SHORT:
        return sizeof(short);
    case NC_INT:
        return sizeof(int);
    case NC_FLOAT:
        return sizeof(float);
    case NC_DOUBLE:
        return sizeof(double);
    default:
        quit("ncw_sizeof(): unknown type\n");
    }

    return UINT_MAX;
}

/** Copies all dimensions from one NetCDF file to another.
 *
 * @param fname_src Source file name
 * @param ncid_src Source file id
 * @param fname_dst Destination file name
 * @param ncid_dst Destination file id
 */
void ncw_copy_dims(const char* fname_src, int ncid_src, const char* fname_dst, int ncid_dst)
{
    int ndims;
    int unlimdimid = -1;
    int i;

    ncw_inq_ndims(fname_src, ncid_src, &ndims);
    ncw_inq_unlimdimid(fname_src, ncid_src, &unlimdimid);

    for (i = 0; i < ndims; ++i) {
        char dimname[NC_MAX_NAME] = STR_UNKNOWN;
        size_t size;
        int dimid;

        memset(dimname, 0, NC_MAX_NAME);

        ncw_inq_dim(fname_src, ncid_src, i, dimname, &size);
        if (i == unlimdimid)
            size = NC_UNLIMITED;
        if (!ncw_dim_exists(ncid_dst, dimname))
            ncw_def_dim(fname_dst, ncid_dst, dimname, size, &dimid);
    }
}

/** Copies dimension from one NetCDF file to another.
 *
 * @param fname_src Source file name
 * @param ncid_src Source file id
 * @param dimname Dimension name
 * @param fname_dst Destination file name
 * @param ncid_dst Destination file id
 */
void ncw_copy_dim(const char* fname_src, int ncid_src, const char dimname[], const char* fname_dst, int ncid_dst)
{
    int unlimdimid = -1;
    int dimid;
    size_t size;

    ncw_inq_unlimdimid(fname_src, ncid_src, &unlimdimid);
    ncw_inq_dimid(fname_src, ncid_src, dimname, &dimid);
    if (dimid == unlimdimid)
        size = NC_UNLIMITED;
    else
        ncw_inq_dimlen(fname_src, ncid_src, dimid, &size);
    ncw_def_dim(fname_dst, ncid_dst, dimname, size, &dimid);
}

/** Copies definition of a specified variable from one NetCDF file to another.
 *
 * @param fname_src Source file name
 * @param ncid_src Source file id
 * @param vid_src Variable id
 * @param fname_dest Destination file name
 * @param ncid_dest Destination file id
 * @return Variable id in the destination file
 */
int ncw_copy_vardef(const char* fname_src, int ncid_src, int vid_src, const char* fname_dst, int ncid_dst)
{
    char varname[NC_MAX_NAME] = STR_UNKNOWN;
    int vid_dst;
    nc_type type;
    int ndims;
    int dimids_src[NC_MAX_DIMS], dimids_dst[NC_MAX_DIMS];
    int natts;
    int status;
    int i;

    status = nc_redef(ncid_dst);

    ncw_inq_varname(fname_src, ncid_src, vid_src, varname);
    ncw_inq_var(fname_src, ncid_src, vid_src, NULL, &type, &ndims, dimids_src, &natts);

    for (i = 0; i < ndims; ++i) {
        char dimname[NC_MAX_NAME];
        size_t len;

        ncw_inq_dim(fname_src, ncid_src, dimids_src[i], dimname, &len);
        if (!ncw_dim_exists(ncid_dst, dimname)) {
            int unlimdimid = -1;

            ncw_inq_unlimdimid(fname_src, ncid_src, &unlimdimid);
            ncw_def_dim(fname_dst, ncid_dst, dimname, (dimids_src[i] == unlimdimid) ? NC_UNLIMITED : len, &dimids_dst[i]);
        } else
            ncw_inq_dimid(fname_dst, ncid_dst, dimname, &dimids_dst[i]);
    }

    ncw_def_var(fname_dst, ncid_dst, varname, type, ndims, dimids_dst, &vid_dst);
    ncw_copy_atts(fname_src, ncid_src, vid_src, fname_dst, ncid_dst, vid_dst);

    if (status == NC_NOERR)
        nc_enddef(ncid_dst);

    return vid_dst;
}

/** Copies data for a given variable from one file to another.
 *  Note that variable IDs may be different, but name and dimensions must be
 *  the same.
 *
 * @param fname_src Source file name
 * @param ncid_src Source file id
 * @param vid_src Variable id
 * @param fname_dest Destination file name
 * @param ncid_dest Destination file id
 */
void ncw_copy_vardata(const char* fname_src, int ncid_src, int vid_src, const char* fname_dst, int ncid_dst)
{
    char varname[NC_MAX_NAME] = STR_UNKNOWN;
    nc_type type;
    int ndims;
    int dimids[NC_MAX_DIMS];
    size_t dimlens[NC_MAX_DIMS];
    int natts;
    void* data = NULL;
    int vid_dst = -1;
    int size;
    int status;
    int i;

    status = nc_enddef(ncid_dst);

    ncw_inq_varname(fname_src, ncid_src, vid_src, varname);
    ncw_inq_var(fname_src, ncid_src, vid_src, NULL, &type, &ndims, dimids, &natts);
    for (i = 0; i < ndims; ++i)
        ncw_inq_dimlen(fname_src, ncid_src, dimids[i], &dimlens[i]);
    size = 1;
    for (i = 0; i < ndims; ++i)
        size *= dimlens[i];
    ncw_inq_varid(fname_dst, ncid_dst, varname, &vid_dst);

    size *= ncw_sizeof(type);
    data = malloc(size);
    switch (type) {
    case NC_BYTE:
    case NC_CHAR:
        ncw_get_var_text(fname_src, ncid_src, vid_src, data);
        ncw_put_var_text(fname_dst, ncid_dst, vid_dst, data);
        break;
    case NC_SHORT:
        ncw_get_var_short(fname_src, ncid_src, vid_src, data);
        ncw_put_var_short(fname_dst, ncid_dst, vid_dst, data);
        break;
    case NC_INT:
        ncw_get_var_int(fname_src, ncid_src, vid_src, data);
        ncw_put_var_int(fname_dst, ncid_dst, vid_dst, data);
        break;
    case NC_FLOAT:
        ncw_get_var_float(fname_src, ncid_src, vid_src, data);
        ncw_put_var_float(fname_dst, ncid_dst, vid_dst, data);
        break;
    case NC_DOUBLE:
        ncw_get_var_double(fname_src, ncid_src, vid_src, data);
        ncw_put_var_double(fname_dst, ncid_dst, vid_dst, data);
        break;
    default:
        quit("\"%s\": ncw_copy_vardata(): variable = \"%s\": unknown data type \"%s\"\n", ncw_nctype2str(type));
    }
    free(data);

    if (status == NC_NOERR)
        nc_redef(ncid_dst);
}

/** Copies definitiona and data for a given variable from one file to another.
 *
 * @param fname_src Source file name
 * @param ncid_src Source file id
 * @param varname Variable name
 * @param fname_dest Destination file name
 * @param ncid_dest Destination file id
 */
void ncw_copy_var(const char* fname_src, int ncid_src, const char varname[], const char* fname_dst, int ncid_dst)
{
    int vid_src = -1;

    ncw_inq_varid(fname_src, ncid_src, varname, &vid_src);
    (void) ncw_copy_vardef(fname_src, ncid_src, vid_src, fname_dst, ncid_dst);
    ncw_copy_vardata(fname_src, ncid_src, vid_src, fname_dst, ncid_dst);
}

/** Set deflation (compression) level for all variables in a file.
 *
 * @param fname NetCDF file name
 * @param ncid NetCDF file id
 * @param shuffle Shuffle flag, whatever it means
 * @param deflate Flag, turns deflation on/off
 * @param deflate_level Compression level
 */
void ncw_def_deflate(const char fname[], int ncid, int shuffle, int deflate, int deflate_level)
{
    int nv = -1;                /* total number of variables */
    int vid;

    ncw_inq_nvars(fname, ncid, &nv);
    for (vid = 0; vid < nv; ++vid) {
        int status = nc_def_var_deflate(ncid, vid, shuffle, deflate, deflate_level);

        if (status != NC_NOERR) {
            char varname[NC_MAX_NAME] = "STR_UNKNOWN";

            _ncw_inq_varname(fname, ncid, vid, varname);
            quit("\"%s\": ncw_def_deflate(): failed for \"%s\": %s\n", fname, varname, nc_strerror(status));
        }
    }
}

/** Gets the id for the first dimension found to be present in a NetCDF file
 * out of two dimensions specified by names. Useful for handling both new and
 * old data formats.
 *
 * @param fname NetCDF file name
 * @param ncid NetCDF file id
 * @param name1 The first dimension name to be tried
 * @param name2 The second dimension name to be tried
 * @param dimid Dimension id (output)
 */
void ncw_inq_dimid2(const char fname[], int ncid, const char dimname1[], const char dimname2[], int* dimid)
{
    int status1 = nc_inq_dimid(ncid, dimname1, dimid);

    if (status1 != NC_NOERR) {
        int status2 = nc_inq_dimid(ncid, dimname2, dimid);

        if (status2 != NC_NOERR) {
            quit("\"%s\": nc_inq_dimid(): failed for dimname = \"%s\": %s, and for dimname = \"%s\": %s\n", fname, dimname1, nc_strerror(status1), dimname2, nc_strerror(status2));
        }
    }
}

/** Gets the value(s) (converted to int type) of the first attribute found to
 * be present in a NetCDF file out of two specified by attribute names. Useful
 * for handling both new and old data formats.
 *
 * @param fname NetCDF file name
 * @param ncid NetCDF file id
 * @param varid Name of the variable the attrubute belongs to (use NC_GLOBAL
 *              for global attributes)
 * @param attname1 The first attribute name to be tried
 * @param attname2 The second attribute name to be tried
 * @param in Attribute value(s) (output)
 */
void ncw_get_att_int2(const char fname[], int ncid, int varid, const char attname1[], const char attname2[], int v[])
{
    nc_type xtype;
    size_t len;
    int status1 = nc_get_att_int(ncid, varid, attname1, v);

    if (status1 != NC_NOERR) {
        int status2 = nc_get_att_int(ncid, varid, attname2, v);

        if (status2 != NC_NOERR) {
            char varname[NC_MAX_NAME] = STR_UNKNOWN;

            _ncw_inq_varname(fname, ncid, varid, varname);
            quit("\"%s\": nc_get_att_int(): failed for varid = %d (varname = \"%s\"), attname = \"%s\": %s, and attname = \"%s\": %s\n", fname, varid, varname, attname1, nc_strerror(status1), attname2, nc_strerror(status2));
        }
    }

    ncw_inq_att(fname, ncid, varid, (status1 == NC_NOERR) ? attname1 : attname2, &xtype, &len);

    if (xtype != NC_INT && xtype != NC_BYTE && xtype != NC_SHORT) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        _ncw_inq_varname(fname, ncid, varid, varname);
        quit("ncw_get_att_int2(): failed for varid = %d (varname = \"%s\"), attname = \"%s\": the attribute is of %s type\n", fname, varid, varname, (status1 == NC_NOERR) ? attname1 : attname2, ncw_nctype2str(xtype));
    }
}

/** Finds all variables with specified dimensions and specified attribute.
 *
 * @param fname NetCDF file name
 * @param ncid NetCDF file id
 * @param ndims Number of dimensions (0 if any)
 * @param dims Array of dimension ids
 * @param attname Attribute name that is be set (NULL if any)
 * @param attval Attribute value (NULL if any)
 * @param nvars Return -- number of variables found
 * @param vids Return -- array of variable ids (needs to be deallocated when
 *             necessary)
 */
void ncw_find_vars(const char fname[], int ncid, int ndims, const int dims[], const char attname[], const void* attval, int* nvars, int** vids)
{
    int nv = -1;                /* total number of variables */
    int nallocated = 0;
    int vid;

    *nvars = 0;
    *vids = NULL;

    ncw_inq_nvars(fname, ncid, &nv);

    for (vid = 0; vid < nv; ++vid) {
        char varname[NC_MAX_NAME] = "STR_UNKNOWN";
        nc_type vtype = -1;
        int nd = -1;
        int dids[NC_MAX_DIMS];
        int natts = -1;
        int i;

        /*
         * check if the dimensions are right for the variable
         */
        if (ndims > 0) {
            ncw_inq_var(fname, ncid, vid, varname, &vtype, &nd, dids, &natts);

            if (ndims != nd)
                goto nextvar;   /* (not using "continue" to be consistent) */

            /*
             * check dimensions
             */
            for (i = 0; i < ndims; ++i) {
                if (dids[i] != dims[i])
                    goto nextvar;
            }
        }

        /*
         * (so far so good)
         */

        /*
         * check if the attribute is present for the variable
         */
        if (attname != NULL) {
            for (i = 0; i < natts; ++i) {
                char aname[NC_MAX_NAME] = STR_UNKNOWN;

                ncw_inq_attname(fname, ncid, vid, i, aname);

                if (strcmp(aname, attname) != 0)
                    continue;

                if (attval != NULL) {
                    size_t alen = UINT_MAX;
                    nc_type atype = -1;
                    char* aval = NULL;

                    ncw_inq_att(fname, ncid, vid, aname, &atype, &alen);
                    if (alen <= 0)
                        continue;
                    aval = calloc(alen, ncw_sizeof(atype));
                    if (aval == NULL)
                        quit("\"%s\": ncw_find_vars(): could not allocate memory for attribute = \"%s\", type = %s, length = %d, varid = %d (varname = \"%s\"): %s\n", fname, aname, ncw_nctype2str(atype), alen, vid, (vid == NC_GLOBAL) ? "NC_GLOBAL" : varname, strerror(errno));
                    ncw_get_att_text(fname, ncid, vid, attname, aval);
                    if (memcmp(attval, aval, alen * sizeof(atype)) == 0) {
                        free(aval);
                        break;
                    }
                    free(aval);
                } else
                    break;
            }
            if (i >= natts)
                goto nextvar;
        }

        /*
         * the variable matches the criteria; add the id to the id array
         */

        /*
         * expand alocated space for the id array if necessary
         */
        if (nallocated == 0) {
            nallocated = NALLOCATED_START;
            *vids = malloc(nallocated * sizeof(int));
        } else if (*nvars == nallocated) {
            nallocated *= 2;
            *vids = realloc(vids, nallocated * sizeof(int));
        }

        (*vids)[*nvars] = vid;
        (*nvars)++;

      nextvar:
        ;
    }
}

/** Finds the first (id-wise) variable dependent on record dimension only.
 *
 * @param fname NetCDF file name
 * @param ncid NetCDF file id
 * @param varid Return -- id of the variable
 */
void ncw_find_timevarid(const char fname[], int ncid, int* varid)
{
    int unlimdimid = -1;
    int nvars = -1;
    int i;

    ncw_inq_unlimdimid(fname, ncid, &unlimdimid);
    ncw_inq_nvars(fname, ncid, &nvars);

    for (i = 0; i < nvars; ++i) {
        int ndims = -1;
        int dimid = -1;

        ncw_inq_varndims(fname, ncid, i, &ndims);
        if (ndims != 1)
            continue;
        ncw_inq_vardimid(fname, ncid, i, &dimid);
        if (dimid == unlimdimid) {
            *varid = i;
            return;
        }
    }
    quit("\"%s\": ncw_find_timevar(): could not find any\n");
}

/** Checks if specified attribute exists in the file.
 *
 * @param ncid NetCDF file id
 * @param varid Variable id (NC_GLOBAL for a global attribute)
 * @param attname Attribute name
 * @return 1 if attribute exists, 0 if not
 */
int ncw_att_exists(int ncid, int varid, const char attname[])
{
    if (nc_inq_attid(ncid, varid, attname, NULL) == NC_NOERR)
        return 1;
    else
        return 0;
}

/** Checks if specified variable exists in the file.
 *
 * @param ncid NetCDF file id
 * @param varname Variable name
 * @return 1 if attribute exists, 0 if not
 */
int ncw_var_exists(int ncid, const char varname[])
{
    int varid;

    if (nc_inq_varid(ncid, varname, &varid) == NC_NOERR)
        return 1;
    else
        return 0;
}

/** Checks if specified dimension exists in the file.
 *
 * @param ncid NetCDF file id
 * @param dimname Dimension name
 * @return 1 if attribute exists, 0 if not
 */
int ncw_dim_exists(int ncid, const char dimname[])
{
    int dimid;

    if (nc_inq_dimid(ncid, dimname, &dimid) == NC_NOERR)
        return 1;
    else
        return 0;
}

/** Copies all attributes of a specified variable from one NetCDF file to
 * another.
 *
 * @param fname_src Source file name
 * @param ncid_src Source file id
 * @param varid Variable id
 * @param fname_dest Destination file name
 * @param ncid_dest Destination file id
 */
void ncw_copy_atts(const char fname_src[], int ncid_src, int varid_src, const char* fname_dst, int ncid_dst, int varid_dst)
{
    char varname_src[NC_MAX_NAME] = STR_UNKNOWN;
    char varname_dst[NC_MAX_NAME] = STR_UNKNOWN;
    int status;
    int natts;
    int i;

    _ncw_inq_varname(fname_src, ncid_src, varid_src, varname_src);
    _ncw_inq_varname(fname_dst, ncid_dst, varid_dst, varname_dst);
    ncw_inq_varnatts(fname_src, ncid_src, varid_src, &natts);

    for (i = 0; i < natts; ++i) {
        char attname[NC_MAX_NAME] = STR_UNKNOWN;

        ncw_inq_attname(fname_src, ncid_src, varid_src, i, attname);
        if ((status = nc_copy_att(ncid_src, varid_src, attname, ncid_dst, varid_dst)) != NC_NOERR)
            quit("\"%s\": -> %s:  nc_copy_att(): failed for varid_in = %d (varname = \"%s\"), attname = \"%s\", varid_dst = %d, (varname = \"%s\"): %s\n", fname_src, fname_dst, varid_src, varname_src, attname, varid_dst, varname_dst, nc_strerror(status));
    }
}

/** Defines a new variable using an existing variable as a template.
 *
 * @param fname NetCDF file name
 * @param ncid NetCDF file id
 * @param oldvarname Name of the existing variable
 * @param newvarname Name of the new variable
 */
void ncw_def_var_as(const char fname[], int ncid, const char oldvarname[], const char newvarname[])
{
    int oldvarid, newvarid;
    nc_type type;
    int ndims;
    int dimids[NC_MAX_DIMS];
    int natts;
    int i;

    ncw_inq_varid(fname, ncid, oldvarname, &oldvarid);
    ncw_inq_var(fname, ncid, oldvarid, NULL, &type, &ndims, dimids, &natts);

    ncw_def_var(fname, ncid, newvarname, type, ndims, dimids, &newvarid);

    for (i = 0; i < natts; ++i) {
        char attname[NC_MAX_NAME] = STR_UNKNOWN;

        ncw_inq_attname(fname, ncid, oldvarid, i, attname);
        ncw_copy_att(fname, ncid, oldvarid, attname, fname, ncid, newvarid);
    }
}

/** Reads one record of a variable.
 * @param fname NetCDF file name
 * @param ncid NetCDF file id
 * @param varid ID of the variable
 * @param r Record number
 * @param v The data
 */
void ncw_get_var_double_record(const char fname[], int ncid, int varid, int r, double v[])
{
    int ndims;
    int dimids[NC_MAX_DIMS];
    size_t dimlen[NC_MAX_VAR_DIMS];
    size_t start[NC_MAX_VAR_DIMS];
    int i;
    int status;

    ncw_inq_varndims(fname, ncid, varid, &ndims);
    ncw_inq_vardimid(fname, ncid, varid, dimids);
    start[0] = r;               /* this record */
    dimlen[0] = 1;              /* one record only */
    for (i = 1; i < ndims; ++i) {
        ncw_inq_dimlen(fname, ncid, dimids[i], &dimlen[i]);
        start[i] = 0;
    }
    status = nc_get_vara_double(ncid, varid, start, dimlen, v);

    if (status != 0) {
        char varname[NC_MAX_NAME] = "STR_UNKNOWN";

        ncw_inq_varname(fname, ncid, varid, varname);
        quit("\"%s\": ncw_get_var_double_record(): failed to read record (outer dim) %d for \"%s\": %s\n", fname, r, varname, nc_strerror(status));
    }
}

/** Reads one record of a variable.
 * @param fname NetCDF file name
 * @param ncid NetCDF file id
 * @param varid ID of the variable
 * @param r Record number
 * @param v The data
 */
void ncw_get_var_float_record(const char fname[], int ncid, int varid, int r, float v[])
{
    int ndims;
    int dimids[NC_MAX_DIMS];
    size_t dimlen[NC_MAX_VAR_DIMS];
    size_t start[NC_MAX_VAR_DIMS];
    int i;
    int status;

    ncw_inq_varndims(fname, ncid, varid, &ndims);
    ncw_inq_vardimid(fname, ncid, varid, dimids);
    start[0] = r;               /* this record */
    dimlen[0] = 1;              /* one record only */
    for (i = 1; i < ndims; ++i) {
        ncw_inq_dimlen(fname, ncid, dimids[i], &dimlen[i]);
        start[i] = 0;
    }
    status = nc_get_vara_float(ncid, varid, start, dimlen, v);

    if (status != 0) {
        char varname[NC_MAX_NAME] = "STR_UNKNOWN";

        ncw_inq_varname(fname, ncid, varid, varname);
        quit("\"%s\": ncw_get_var_float_record(): failed to read record (outer dim) %d for \"%s\": %s\n", fname, r, varname, nc_strerror(status));
    }
}

/** Writes one record of a variable.
 * @param fname NetCDF file name
 * @param ncid NetCDF file id
 * @param varid ID of the variable
 * @param r Record number
 * @param v The data
 */
void ncw_put_var_double_record(const char fname[], int ncid, int varid, int r, double v[])
{
    int ndims;
    int dimids[NC_MAX_DIMS];
    size_t dimlen[NC_MAX_VAR_DIMS];
    size_t start[NC_MAX_VAR_DIMS];
    int i;
    int status;

    ncw_inq_varndims(fname, ncid, varid, &ndims);
    ncw_inq_vardimid(fname, ncid, varid, dimids);
    start[0] = r;               /* this record */
    dimlen[0] = 1;              /* one record only */
    for (i = 1; i < ndims; ++i) {
        ncw_inq_dimlen(fname, ncid, dimids[i], &dimlen[i]);
        start[i] = 0;
    }
    status = nc_put_vara_double(ncid, varid, start, dimlen, v);

    if (status != 0) {
        char varname[NC_MAX_NAME] = "STR_UNKNOWN";

        ncw_inq_varname(fname, ncid, varid, varname);
        quit("\"%s\": ncw_put_var_double_record(): failed to read record (outer dim) %d for \"%s\": %s\n", fname, r, varname, nc_strerror(status));
    }
}

/** Writes one record of a variable.
 * @param fname NetCDF file name
 * @param ncid NetCDF file id
 * @param varid ID of the variable
 * @param r Record number
 * @param v The data
 */
void ncw_put_var_float_record(const char fname[], int ncid, int varid, int r, float v[])
{
    int ndims;
    int dimids[NC_MAX_DIMS];
    size_t dimlen[NC_MAX_VAR_DIMS];
    size_t start[NC_MAX_VAR_DIMS];
    int i;
    int status;

    ncw_inq_varndims(fname, ncid, varid, &ndims);
    ncw_inq_vardimid(fname, ncid, varid, dimids);
    start[0] = r;               /* this record */
    dimlen[0] = 1;              /* one record only */
    for (i = 1; i < ndims; ++i) {
        ncw_inq_dimlen(fname, ncid, dimids[i], &dimlen[i]);
        start[i] = 0;
    }
    status = nc_put_vara_float(ncid, varid, start, dimlen, v);

    if (status != 0) {
        char varname[NC_MAX_NAME] = "STR_UNKNOWN";

        ncw_inq_varname(fname, ncid, varid, varname);
        quit("\"%s\": ncw_put_var_float_record(): failed to read record (outer dim) %d for \"%s\": %s\n", fname, r, varname, nc_strerror(status));
    }
}

/** Checks that the attribute has certain type and length
 */
void ncw_check_att(const char fname[], int ncid, int varid, const char attname[], nc_type att_type, size_t att_len)
{
    nc_type type;
    size_t len;
    int status = nc_inq_att(ncid, varid, attname, &type, &len);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        _ncw_inq_varname(fname, ncid, varid, varname);
        quit("\"%s\": nc_inq_att(): failed for varid = %d (varname = \"%s\"), attname = \"%s\", atttype = %d: %s\n", fname, varid, varname, attname, ncw_nctype2str(type), nc_strerror(status));
    }

    if (type != att_type)
        quit("\"%s\": ncw_check_att(): attribute \"%s\" is supposed to have type \"%s\"; the actual type is \"%s\"\n", fname, attname, ncw_nctype2str(att_type), ncw_nctype2str(type));
    if (len != att_len)
        quit("\"%s\": ncw_check_att(): attribute \"%s\" is supposed to have length %z; the actual length is %z\n", fname, attname, att_len, len);
}

/** Check that the dimension has certain length
 */
void ncw_check_dimlen(const char fname[], int ncid, const char dimname[], size_t dimlen)
{
    size_t len;
    int dimid;
    int status;

    ncw_inq_dimid(fname, ncid, dimname, &dimid);
    status = nc_inq_dimlen(ncid, dimid, &len);

    if (status != NC_NOERR)
        quit("\"%s\": nc_inq_dimlen(): failed for dimid = %d (dimname = \"%s\"): %s\n", fname, dimid, dimname, nc_strerror(status));
    if (len != dimlen)
        quit("\"%s\": ncw_check_dimlen(): dimension \"%s\" is supposed to have length %zu; the actual length is %zu\n", fname, dimname, dimlen, len);
}

/** Sets the quit function.
 */
void ncw_set_quitfn(ncw_quit_fn quitfn)
{
    quit = quitfn;
}
