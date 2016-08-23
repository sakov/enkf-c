/******************************************************************************
 *
 * File:           ncw.c
 *  
 * Created         19/10/2000
 *  
 * Author:         Pavel Sakov
 *                 CSIRO Marine Research
 *                 Bureau of Meteorology
 *  
 * Purpose:        Simple wrappers to netcdf library procedures for
 *                 better error messaging.
 *                 Some straightforward extensions to netcdf library 
 *                 procedures.
 *
 * Revisions:      PS 19/08/2003 Changed library name from "nc" to "ncw". 
 *                    Changed procedure prefixes from "NC_" to "ncw_".
 *                 PS 24/08/2016 Version 2.00: eliminated file names from the
 *                    functions' arguments, using nc_inq_path() instead. Most
 *                    of the wrappers have now the same format as their siblings
 *                    in the NetCDF library.
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

const char ncw_version[] = "2.00";

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

static char* _ncw_get_path(int ncid)
{
    char* path = NULL;
    size_t len = -1;
    int status = nc_inq_path(ncid, &len, NULL);

    if (status != NC_NOERR)
        quit("nc_inq_path(): failed for ncid = %d: %s\n", ncid, nc_strerror(status));
    path = malloc(len + 1);
    status = nc_inq_path(ncid, NULL, path);

    return path;
}

static void _ncw_inq_varname(int ncid, int varid, char varname[])
{
    int status;

    if (varid == -1)
        strcpy(varname, "NC_GLOBAL");
    else if ((status = nc_inq_varname(ncid, varid, varname)) != NC_NOERR)
        quit("\"%s\": nc_inq_varname(): failed for varid = %d: %s\n", _ncw_get_path(ncid), varid, nc_strerror(status));
}

static void _ncw_inq_varid(int ncid, const char varname[], int* varid)
{
    int status;

    if (strcmp(varname, "NC_GLOBAL") == 0)
        *varid = -1;
    else if ((status = nc_inq_varid(ncid, varname, varid)) != NC_NOERR)
        quit("\"%s\": nc_inq_varid(): failed for varid = %d: %s\n", _ncw_get_path(ncid), varid, nc_strerror(status));
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

/** Sets the quit function.
 */
void ncw_set_quitfn(ncw_quit_fn quitfn)
{
    quit = quitfn;
}

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

void ncw_redef(int ncid)
{
    int status = nc_redef(ncid);

    if (status != NC_NOERR)
        quit("\"%s\": ncredef(): failed: %s\n", _ncw_get_path(ncid), nc_strerror(status));
}

void ncw_enddef(int ncid)
{
    int status = nc_enddef(ncid);

    if (status != NC_NOERR)
        quit("\"%s\": ncenddef(): failed: %s\n", _ncw_get_path(ncid), nc_strerror(status));
}

void ncw_sync(int ncid)
{
    int status;

    nc_enddef(ncid);

    status = ncsync(ncid);

    if (status != NC_NOERR)
        quit("\"%s\": ncsync(): failed: %s\n", _ncw_get_path(ncid), nc_strerror(status));
}

void ncw_close(int ncid)
{
    int status = nc_close(ncid);

    if (status != NC_NOERR)
        quit("\"%s\": nc_close(): failed: %s\n", _ncw_get_path(ncid), nc_strerror(status));
}

void ncw_inq(int ncid, int* ndims, int* nvars, int* natts, int* unlimdimid)
{
    int status = nc_inq(ncid, ndims, nvars, natts, unlimdimid);

    if (status != NC_NOERR)
        quit("\"%s\": nc_inq(): failed: %s\n", _ncw_get_path(ncid), nc_strerror(status));
}

void ncw_inq_ndims(int ncid, int* ndims)
{
    int status = nc_inq_ndims(ncid, ndims);

    if (status != NC_NOERR)
        quit("\"%s\": nc_inq_ndims(): failed: %s\n", _ncw_get_path(ncid), nc_strerror(status));
}

void ncw_inq_nvars(int ncid, int* nvars)
{
    int status = nc_inq_nvars(ncid, nvars);

    if (status != NC_NOERR)
        quit("\"%s\": nc_inq_nvars(): failed: %s\n", _ncw_get_path(ncid), nc_strerror(status));
}

void ncw_inq_natts(int ncid, int* natts)
{
    int status = nc_inq_natts(ncid, natts);

    if (status != NC_NOERR)
        quit("\"%s\": nc_inq_natts(): failed: %s\n", _ncw_get_path(ncid), nc_strerror(status));
}

void ncw_inq_unlimdimid(int ncid, int* unlimdimid)
{
    int status = nc_inq_unlimdim(ncid, unlimdimid);

    if (status != NC_NOERR)
        quit("\"%s\": nc_inq_unlimdim(): failed: %s\n", _ncw_get_path(ncid), nc_strerror(status));
}

void ncw_def_dim(int ncid, const char dimname[], size_t len, int* dimid)
{
    int status = nc_def_dim(ncid, dimname, len, dimid);

    if (status != NC_NOERR)
        quit("\"%s\": nc_def_dim(): failed for dimname = \"%s\", dimlen = %d: %s\n", _ncw_get_path(ncid), dimname, len, nc_strerror(status));
}

void ncw_inq_dimid(int ncid, const char dimname[], int* dimid)
{
    int status = nc_inq_dimid(ncid, dimname, dimid);

    if (status != NC_NOERR)
        quit("\"%s\": nc_inq_dimid(): failed for dimname = \"%s\": %s\n", _ncw_get_path(ncid), dimname, nc_strerror(status));
}

void ncw_inq_dim(int ncid, int dimid, char dimname[], size_t* len)
{
    int status = nc_inq_dim(ncid, dimid, dimname, len);

    if (status != NC_NOERR)
        quit("\"%s\": nc_inq_dim(): failed for dimid = %d: %s\n", _ncw_get_path(ncid), dimid, nc_strerror(status));
}

void ncw_inq_dimname(int ncid, int dimid, char dimname[])
{
    int status = nc_inq_dimname(ncid, dimid, dimname);

    if (status != NC_NOERR)
        quit("\"%s\": nc_inq_dimname(): failed for dimid = %d: %s\n", _ncw_get_path(ncid), dimid, nc_strerror(status));
}

void ncw_inq_dimlen(int ncid, int dimid, size_t* len)
{
    int status = nc_inq_dimlen(ncid, dimid, len);

    if (status != NC_NOERR) {
        char dimname[NC_MAX_NAME] = STR_UNKNOWN;

        ncw_inq_dimname(ncid, dimid, dimname);
        quit("\"%s\": nc_inq_dimlen(): failed for dimid = %d (dimname = \"%s\"): %s\n", _ncw_get_path(ncid), dimid, dimname, nc_strerror(status));
    }
}

void ncw_rename_dim(int ncid, const char oldname[], const char newname[])
{
    int dimid;
    int status;

    ncw_inq_dimid(ncid, oldname, &dimid);
    status = nc_rename_dim(ncid, dimid, newname);

    if (status != NC_NOERR)
        quit("\"%s\": nc_rename_dim(): failed for dimid = %d (oldname = \"%s\", newname = \"%s\"): %s\n", dimid, oldname, newname, nc_strerror(status));
}

void ncw_def_var(int ncid, const char varname[], nc_type xtype, int ndims, const int dimids[], int* varid)
{
    int status = nc_def_var(ncid, varname, xtype, ndims, dimids, varid);

    if (status != NC_NOERR)
        quit("\"%s\": nc_def_var(): failed for varname = \"%s\", vartype = %s, ndims = %d, dimids = %s: %s\n", _ncw_get_path(ncid), varname, ncw_nctype2str(xtype), ndims, uint2str(ndims, (const unsigned int*) dimids), nc_strerror(status));
}

void ncw_inq_varid(int ncid, const char varname[], int* varid)
{
    int status = nc_inq_varid(ncid, varname, varid);

    if (status != NC_NOERR)
        quit("\"%s\": nc_inq_varid(): failed for varname = \"%s\": %s\n", _ncw_get_path(ncid), varname, nc_strerror(status));
}

void ncw_inq_var(int ncid, int varid, char varname[], nc_type* xtype, int* ndims, int dimids[], int* natts)
{
    int status = nc_inq_var(ncid, varid, varname, xtype, ndims, dimids, natts);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        ncw_inq_varname(ncid, varid, varname);
        quit("\"%s\": nc_inq_var(): failed for varid = %d (varname = \"%s\"): %s\n", _ncw_get_path(ncid), varid, varname, nc_strerror(status));
    }
}

void ncw_inq_varname(int ncid, int varid, char varname[])
{
    int status = nc_inq_varname(ncid, varid, varname);

    if (status != NC_NOERR)
        quit("\"%s\": nc_inq_varname(): failed for varid = %d: %s\n", _ncw_get_path(ncid), varid, nc_strerror(status));
}

void ncw_inq_vartype(int ncid, int varid, nc_type* xtype)
{
    int status = nc_inq_vartype(ncid, varid, xtype);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        ncw_inq_varname(ncid, varid, varname);
        quit("\"%s\": nc_inq_vartype(): failed for varid = %d (varname = %s): %s\n", _ncw_get_path(ncid), varid, varname, nc_strerror(status));
    }
}

void ncw_inq_varndims(int ncid, int varid, int* ndims)
{
    int status = nc_inq_varndims(ncid, varid, ndims);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = "STR_UNKNOWN";

        ncw_inq_varname(ncid, varid, varname);
        quit("\"%s\": nc_inq_varndims(): failed for varid = %d (varname = \"%s\"): %s\n", _ncw_get_path(ncid), varid, varname, nc_strerror(status));
    }
}

void ncw_inq_vardimid(int ncid, int varid, int dimids[])
{
    int status = nc_inq_vardimid(ncid, varid, dimids);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        ncw_inq_varname(ncid, varid, varname);
        quit("\"%s\": ncw_inq_vardimid(): failed for varid = %d (varname = \"%s\"): %s\n", _ncw_get_path(ncid), varid, varname, nc_strerror(status));
    }
}

void ncw_inq_varnatts(int ncid, int varid, int* natts)
{
    int status = nc_inq_varnatts(ncid, varid, natts);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = "STR_UNKNOWN";

        _ncw_inq_varname(ncid, varid, varname);
        quit("\"%s\": nc_inq_varnatts(): failed for varid = %d (varnatts = \"%s\"): %s\n", _ncw_get_path(ncid), varid, varname, nc_strerror(status));
    }
}

void ncw_rename_var(int ncid, const char oldname[], const char newname[])
{
    int varid;
    int status;

    ncw_inq_varid(ncid, oldname, &varid);
    status = nc_rename_var(ncid, varid, newname);

    if (status != NC_NOERR)
        quit("\"%s\": nc_rename_var(): failed for varid = %d (oldname = \"%s\", newname = \"%s\"): %s\n", _ncw_get_path(ncid), varid, oldname, newname, nc_strerror(status));
}

void ncw_def_var_deflate(int ncid, int varid, int shuffle, int deflate, int deflate_level)
{
    int status = nc_def_var_deflate(ncid, varid, shuffle, deflate, deflate_level);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = "STR_UNKNOWN";

        _ncw_inq_varname(ncid, varid, varname);
        quit("\"%s\": ncw_def_var_deflate(): failed for varid = %d (varname = \"%s\"): %s\n", _ncw_get_path(ncid), varid, varname, nc_strerror(status));
    }
}

void ncw_put_var_text(int ncid, int varid, const char v[])
{
    int status = nc_put_var_text(ncid, varid, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        ncw_inq_varname(ncid, varid, varname);
        quit("\"%s\": nc_put_var_text(): failed for varid = %d (varname = \"%s\"): %s\n", _ncw_get_path(ncid), varid, varname, nc_strerror(status));
    }
}

void ncw_put_var_short(int ncid, int varid, const short v[])
{
    int status = nc_put_var_short(ncid, varid, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        ncw_inq_varname(ncid, varid, varname);
        quit("\"%s\": nc_put_var_short(): failed for varid = %d (varname = \"%s\"): %s\n", _ncw_get_path(ncid), varid, varname, nc_strerror(status));
    }
}

void ncw_put_var_ushort(int ncid, int varid, const unsigned short v[])
{
    int status = nc_put_var_ushort(ncid, varid, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        ncw_inq_varname(ncid, varid, varname);
        quit("\"%s\": nc_put_var_ushort(): failed for varid = %d (varname = \"%s\"): %s\n", _ncw_get_path(ncid), varid, varname, nc_strerror(status));
    }
}

void ncw_put_var_int(int ncid, int varid, const int v[])
{
    int status = nc_put_var_int(ncid, varid, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        ncw_inq_varname(ncid, varid, varname);
        quit("\"%s\": nc_put_var_int(): failed for varid = %d (varname = \"%s\"): %s\n", _ncw_get_path(ncid), varid, varname, nc_strerror(status));
    }
}

void ncw_put_var_uint(int ncid, int varid, const unsigned int v[])
{
    int status = nc_put_var_uint(ncid, varid, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        ncw_inq_varname(ncid, varid, varname);
        quit("\"%s\": nc_put_var_uint(): failed for varid = %d (varname = \"%s\"): %s\n", _ncw_get_path(ncid), varid, varname, nc_strerror(status));
    }
}

void ncw_put_var_float(int ncid, int varid, const float v[])
{
    int status = nc_put_var_float(ncid, varid, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        ncw_inq_varname(ncid, varid, varname);
        quit("\"%s\": nc_put_var_float(): failed for varid = %d (varname = \"%s\"): %s\n", _ncw_get_path(ncid), varid, varname, nc_strerror(status));
    }
}

void ncw_put_var_double(int ncid, int varid, const double v[])
{
    int status = nc_put_var_double(ncid, varid, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        ncw_inq_varname(ncid, varid, varname);
        quit("\"%s\": nc_put_var_double(): failed for varid = %d (varname = \"%s\"): %s\n", _ncw_get_path(ncid), varid, varname, nc_strerror(status));
    }
}

void ncw_get_var_text(int ncid, int varid, char v[])
{
    int status = nc_get_var_text(ncid, varid, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        ncw_inq_varname(ncid, varid, varname);
        quit("\"%s\": nc_get_var_text(): failed for varid = %d (varname = \"%s\"): %s\n", _ncw_get_path(ncid), varid, varname, nc_strerror(status));
    }
}

void ncw_get_var_schar(int ncid, int varid, signed char v[])
{
    int status = nc_get_var_schar(ncid, varid, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        ncw_inq_varname(ncid, varid, varname);
        quit("\"%s\": nc_get_var_schar(): failed for varid = %d (varname = \"%s\"): %s\n", _ncw_get_path(ncid), varid, varname, nc_strerror(status));
    }
}

void ncw_get_var_short(int ncid, int varid, short int v[])
{
    int status = nc_get_var_short(ncid, varid, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        ncw_inq_varname(ncid, varid, varname);
        quit("\"%s\": nc_get_var_short(): failed for varid = %d (varname = \"%s\"): %s\n", _ncw_get_path(ncid), varid, varname, nc_strerror(status));
    }
}

void ncw_get_var_ushort(int ncid, int varid, unsigned short int v[])
{
    int status = nc_get_var_ushort(ncid, varid, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        ncw_inq_varname(ncid, varid, varname);
        quit("\"%s\": nc_get_var_ushort(): failed for varid = %d (varname = \"%s\"): %s\n", _ncw_get_path(ncid), varid, varname, nc_strerror(status));
    }
}

void ncw_get_var_int(int ncid, int varid, int v[])
{
    int status = nc_get_var_int(ncid, varid, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        ncw_inq_varname(ncid, varid, varname);
        quit("\"%s\": nc_get_var_int(): failed for varid = %d (varname = \"%s\"): %s\n", _ncw_get_path(ncid), varid, varname, nc_strerror(status));
    }
}

void ncw_get_var_uint(int ncid, int varid, unsigned int v[])
{
    int status = nc_get_var_uint(ncid, varid, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        ncw_inq_varname(ncid, varid, varname);
        quit("\"%s\": nc_get_var_uint(): failed for varid = %d (varname = \"%s\"): %s\n", _ncw_get_path(ncid), varid, varname, nc_strerror(status));
    }
}

void ncw_get_var_float(int ncid, int varid, float v[])
{
    int status = nc_get_var_float(ncid, varid, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        ncw_inq_varname(ncid, varid, varname);
        quit("\"%s\": nc_get_var_float(): failed for varid = %d (varname = \"%s\"): %s\n", _ncw_get_path(ncid), varid, varname, nc_strerror(status));
    }
}

void ncw_get_var_double(int ncid, int varid, double v[])
{
    int status = nc_get_var_double(ncid, varid, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        ncw_inq_varname(ncid, varid, varname);
        quit("\"%s\": nc_get_var_double(): failed for varid = %d (varname = \"%s\"): %s\n", _ncw_get_path(ncid), varid, varname, nc_strerror(status));
    }
}

void ncw_get_var1_double(int ncid, int varid, const size_t len[], double* in)
{
    int status = nc_get_var1_double(ncid, varid, len, in);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        ncw_inq_varname(ncid, varid, varname);
        quit("\"%s\": nc_get_var1_double(): failed for varid = %d (varname = \"%s\"): %s\n", _ncw_get_path(ncid), varid, varname, nc_strerror(status));
    }
}

void ncw_put_vara_text(int ncid, int varid, const size_t start[], const size_t count[], const char v[])
{
    int status = nc_put_vara_text(ncid, varid, start, count, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;
        int ndims = 0;

        ncw_inq_varname(ncid, varid, varname);
        ncw_inq_varndims(ncid, varid, &ndims);
        quit("\"%s\": ncw_put_vara_text(): failed for varid = %d (varname = \"%s\"), start = %s, count = %s: %s\n", _ncw_get_path(ncid), varid, varname, uint2str(ndims, (unsigned int*) start), uint2str(ndims, (unsigned int*) count), nc_strerror(status));
    }
}

void ncw_put_vara_short(int ncid, int varid, const size_t start[], const size_t count[], const short int v[])
{
    int status = nc_put_vara_short(ncid, varid, start, count, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;
        int ndims = 0;

        ncw_inq_varname(ncid, varid, varname);
        ncw_inq_varndims(ncid, varid, &ndims);
        quit("\"%s\": ncw_put_vara_short(): failed for varid = %d (varname = \"%s\"), start = %s, count = %s: %s\n", _ncw_get_path(ncid), varid, varname, uint2str(ndims, (unsigned int*) start), uint2str(ndims, (unsigned int*) count), nc_strerror(status));
    }
}

void ncw_put_vara_ushort(int ncid, int varid, const size_t start[], const size_t count[], const unsigned short int v[])
{
    int status = nc_put_vara_ushort(ncid, varid, start, count, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;
        int ndims = 0;

        ncw_inq_varname(ncid, varid, varname);
        ncw_inq_varndims(ncid, varid, &ndims);
        quit("\"%s\": ncw_put_vara_ushort(): failed for varid = %d (varname = \"%s\"), start = %s, count = %s: %s\n", _ncw_get_path(ncid), varid, varname, uint2str(ndims, (unsigned int*) start), uint2str(ndims, (unsigned int*) count), nc_strerror(status));
    }
}

void ncw_put_vara_int(int ncid, int varid, const size_t start[], const size_t count[], const int v[])
{
    int status = nc_put_vara_int(ncid, varid, start, count, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;
        int ndims = 0;

        ncw_inq_varname(ncid, varid, varname);
        ncw_inq_varndims(ncid, varid, &ndims);
        quit("\"%s\": ncw_put_vara_int(): failed for varid = %d (varname = \"%s\"), start = %s, count = %s: %s\n", _ncw_get_path(ncid), varid, varname, uint2str(ndims, (unsigned int*) start), uint2str(ndims, (unsigned int*) count), nc_strerror(status));
    }
}

void ncw_put_vara_uint(int ncid, int varid, const size_t start[], const size_t count[], const unsigned int v[])
{
    int status = nc_put_vara_uint(ncid, varid, start, count, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;
        int ndims = 0;

        ncw_inq_varname(ncid, varid, varname);
        ncw_inq_varndims(ncid, varid, &ndims);
        quit("\"%s\": ncw_put_vara_uint(): failed for varid = %d (varname = \"%s\"), start = %s, count = %s: %s\n", _ncw_get_path(ncid), varid, varname, uint2str(ndims, (unsigned int*) start), uint2str(ndims, (unsigned int*) count), nc_strerror(status));
    }
}

void ncw_put_vara_float(int ncid, int varid, const size_t start[], const size_t count[], const float v[])
{
    int status = nc_put_vara_float(ncid, varid, start, count, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;
        int ndims = 0;

        ncw_inq_varname(ncid, varid, varname);
        ncw_inq_varndims(ncid, varid, &ndims);
        quit("\"%s\": ncw_put_vara_float(): failed for varid = %d (varname = \"%s\"), start = %s, count = %s: %s\n", _ncw_get_path(ncid), varid, varname, uint2str(ndims, (unsigned int*) start), uint2str(ndims, (unsigned int*) count), nc_strerror(status));
    }
}

void ncw_put_vara_double(int ncid, int varid, const size_t start[], const size_t count[], const double v[])
{
    int status = nc_put_vara_double(ncid, varid, start, count, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;
        int ndims = 0;

        ncw_inq_varname(ncid, varid, varname);
        ncw_inq_varndims(ncid, varid, &ndims);
        quit("\"%s\": nc_put_vara_double(): failed for varid = %d (varname = \"%s\"), start = %s, count = %s: %s\n", _ncw_get_path(ncid), varid, varname, uint2str(ndims, (unsigned int*) start), uint2str(ndims, (unsigned int*) count), nc_strerror(status));
    }
}

void ncw_get_vara_text(int ncid, int varid, const size_t start[], const size_t count[], char v[])
{
    int status = nc_get_vara_text(ncid, varid, start, count, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        ncw_inq_varname(ncid, varid, varname);
        quit("\"%s\": nc_get_vara_text(): failed for varid = %d (varname = \"%s\"): %s\n", _ncw_get_path(ncid), varid, varname, nc_strerror(status));
    }
}

void ncw_get_vara_short(int ncid, int varid, const size_t start[], const size_t count[], short int v[])
{
    int status = nc_get_vara_short(ncid, varid, start, count, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        ncw_inq_varname(ncid, varid, varname);
        quit("\"%s\": nc_get_vara_short(): failed for varid = %d (varname = \"%s\"): %s\n", _ncw_get_path(ncid), varid, varname, nc_strerror(status));
    }
}

void ncw_get_vara_int(int ncid, int varid, const size_t start[], const size_t count[], int v[])
{
    int status = nc_get_vara_int(ncid, varid, start, count, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        ncw_inq_varname(ncid, varid, varname);
        quit("\"%s\": nc_get_vara_int(): failed for varid = %d (varname = \"%s\"): %s\n", _ncw_get_path(ncid), varid, varname, nc_strerror(status));
    }
}

void ncw_get_vara_float(int ncid, int varid, const size_t start[], const size_t count[], float v[])
{
    int status = nc_get_vara_float(ncid, varid, start, count, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        ncw_inq_varname(ncid, varid, varname);
        quit("\"%s\": nc_get_vara_float(): failed for varid = %d (varname = \"%s\"): %s\n", _ncw_get_path(ncid), varid, varname, nc_strerror(status));
    }
}

void ncw_get_vara_double(int ncid, int varid, const size_t start[], const size_t count[], double v[])
{
    int status = nc_get_vara_double(ncid, varid, start, count, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        ncw_inq_varname(ncid, varid, varname);
        quit("\"%s\": nc_get_vara_double(): failed for varid = %d (varname = \"%s\"): %s\n", _ncw_get_path(ncid), varid, varname, nc_strerror(status));
    }
}

void ncw_put_att_text(int ncid, int varid, const char attname[], const char v[])
{
    int status = nc_put_att_text(ncid, varid, attname, strlen(v), v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        _ncw_inq_varname(ncid, varid, varname);
        quit("\"%s\": nc_put_att_text(): failed for varid = %d (varname = \"%s\"), attname = \"%s\", attvalue = \"%s\": %s\n", _ncw_get_path(ncid), varid, varname, attname, v, nc_strerror(status));
    }
}

void ncw_put_att_int(int ncid, int varid, const char attname[], size_t len, const int v[])
{
    int status = nc_put_att_int(ncid, varid, attname, NC_INT, len, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        _ncw_inq_varname(ncid, varid, varname);
        quit("\"%s\": nc_put_att_int(): failed for varid = %d (varname = \"%s\"), attname = \"%s\", attvalue(s) = \"%s\": %s\n", _ncw_get_path(ncid), varid, varname, attname, int2str(len, v), nc_strerror(status));
    }
}

void ncw_put_att_float(int ncid, int varid, const char attname[], size_t len, const float v[])
{
    int status = nc_put_att_float(ncid, varid, attname, NC_FLOAT, len, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = "STR_UNKNOWN";

        _ncw_inq_varname(ncid, varid, varname);
        quit("\"%s\": nc_put_att_float(): failed for varid = %d (varname = \"%s\"), attname = \"%s\", attvalue = \"%s\": %s\n", _ncw_get_path(ncid), varid, varname, attname, float2str(len, v), nc_strerror(status));
    }

    if (strcmp(attname, "_FillValue") == 0 && varid != NC_GLOBAL) {
        nc_type xtype;

        ncw_inq_vartype(ncid, varid, &xtype);
        if (xtype != NC_FLOAT) {
            char varname[NC_MAX_NAME] = "STR_UNKNOWN";

            ncw_inq_varname(ncid, varid, varname);

            quit("\"%s\": ncw_put_att_float(): failed for varid = %d (varname = \"%s\"): attype = NC_FLOAT differs from vartype = %s for attname = \"_FillValue\"\n", _ncw_get_path(ncid), varid, varname, ncw_nctype2str(xtype));
        }
    }
}

void ncw_put_att_double(int ncid, int varid, const char attname[], size_t len, const double v[])
{
    int status = nc_put_att_double(ncid, varid, attname, NC_DOUBLE, len, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        _ncw_inq_varname(ncid, varid, varname);
        quit("\"%s\": nc_put_att_double(): failed for varid = %d (varname = \"%s\"), attname = \"%s\", attvalue(s) = \"%s\": %s\n", _ncw_get_path(ncid), varid, varname, attname, double2str(len, v), nc_strerror(status));
    }

    if (strcmp(attname, "_FillValue") == 0 && varid != NC_GLOBAL) {
        nc_type xtype;

        ncw_inq_vartype(ncid, varid, &xtype);
        if (xtype != NC_DOUBLE) {
            char varname[NC_MAX_NAME] = "STR_UNKNOWN";

            ncw_inq_varname(ncid, varid, varname);

            quit("\"%s\": ncw_put_att_float(): fatal inconsistency for varid = %d (varname = \"%s\"): attype = NC_DOUBLE differs from vartype = %s for attname = \"_FillValue\"\n", _ncw_get_path(ncid), varid, varname, ncw_nctype2str(xtype));
        }
    }
}

void ncw_inq_attname(int ncid, int varid, int attrid, char attname[])
{
    int status = nc_inq_attname(ncid, varid, attrid, attname);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        _ncw_inq_varname(ncid, varid, varname);
        quit("\"%s\": nc_inq_attname(): failed for varid = %d (varname = \"%s\"): %s\n", _ncw_get_path(ncid), varid, varname, nc_strerror(status));
    }
}

void ncw_inq_att(int ncid, int varid, const char attname[], nc_type* xtype, size_t* len)
{
    int status = nc_inq_att(ncid, varid, attname, xtype, len);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        _ncw_inq_varname(ncid, varid, varname);
        quit("\"%s\": nc_inq_att(): failed for varid = %d (varname = \"%s\"), attname = \"%s\", atttype = %d: %s\n", _ncw_get_path(ncid), varid, varname, attname, ncw_nctype2str(*xtype), nc_strerror(status));
    }
}

void ncw_inq_attlen(int ncid, int varid, const char attname[], size_t* len)
{
    int status = nc_inq_attlen(ncid, varid, attname, len);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        _ncw_inq_varname(ncid, varid, varname);
        quit("\"%s\": nc_inq_attlen(): failed for varid = %d (varname = \"%s\"), attname = \"%s\": %s\n", _ncw_get_path(ncid), varid, varname, attname, nc_strerror(status));
    }
}

void ncw_copy_att(int ncid_src, int varid_src, const char attname[], int ncid_dst, int varid_dst)
{
    int status = nc_copy_att(ncid_src, varid_src, attname, ncid_dst, varid_dst);

    if (status != NC_NOERR) {
        char varname_src[NC_MAX_NAME] = STR_UNKNOWN;
        char varname_dst[NC_MAX_NAME] = STR_UNKNOWN;

        _ncw_inq_varname(ncid_src, varid_src, varname_src);
        _ncw_inq_varname(ncid_dst, varid_dst, varname_dst);

        quit("\"%s\": -> %s:  nc_copy_att(): failed for varid_src = %d (varname = \"%s\"), attname = \"%s\", varid_dst = %d, (varname = \"%s\"): %s\n", _ncw_get_path(ncid_src), _ncw_get_path(ncid_dst), varid_src, varname_src, attname, varid_dst, varname_dst, nc_strerror(status));
    }
}

void ncw_rename_att(int ncid, const char varname[], const char oldname[], const char newname[])
{
    int varid;
    int status;

    _ncw_inq_varid(ncid, varname, &varid);
    status = nc_rename_att(ncid, varid, oldname, newname);

    if (status != NC_NOERR)
        quit("\"%s\": nc_rename_att(): failed for varid = %d, oldname = \"%s\", newname = \"%s\": %s\n",  _ncw_get_path(ncid), varid, oldname, newname, nc_strerror(status));
}

void ncw_del_att(int ncid, int varid, const char attname[])
{
    int status = nc_del_att(ncid, varid, attname);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        _ncw_inq_varname(ncid, varid, varname);
        quit("\"%s\": nc_del_att(): failed for varid = %d (varname = \"%s\"): %s\n", _ncw_get_path(ncid), varid, varname, nc_strerror(status));
    }
}

void ncw_get_att_text(int ncid, int varid, const char attname[], char v[])
{
    int status = nc_get_att_text(ncid, varid, attname, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        _ncw_inq_varname(ncid, varid, varname);
        quit("\"%s\": nc_get_att_text(): failed for varid = %d (varname = \"%s\"), attname = \"%s\": %s\n", _ncw_get_path(ncid), varid, varname, attname, nc_strerror(status));
    }
}

void ncw_get_att_short(int ncid, int varid, const char attname[], short int v[])
{
    int status = nc_get_att_short(ncid, varid, attname, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        _ncw_inq_varname(ncid, varid, varname);
        quit("\"%s\": nc_get_att_short(): failed for varid = %d (varname = \"%s\"), attname = \"%s\": %s\n", _ncw_get_path(ncid), varid, varname, attname, nc_strerror(status));
    }
}

void ncw_get_att_int(int ncid, int varid, const char attname[], int v[])
{
    int status = nc_get_att_int(ncid, varid, attname, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        _ncw_inq_varname(ncid, varid, varname);
        quit("\"%s\": nc_get_att_int(): failed for varid = %d (varname = \"%s\"), attname = \"%s\": %s\n", _ncw_get_path(ncid), varid, varname, attname, nc_strerror(status));
    }
}

void ncw_get_att_float(int ncid, int varid, const char attname[], float v[])
{
    int status = nc_get_att_float(ncid, varid, attname, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        _ncw_inq_varname(ncid, varid, varname);
        quit("\"%s\": ncw_get_att_float(): failed for varid = %d (varname = \"%s\"), attname = \"%s\": %s\n", _ncw_get_path(ncid), varid, varname, attname, nc_strerror(status));
    }
}

void ncw_get_att_double(int ncid, int varid, const char attname[], double v[])
{
    int status = nc_get_att_double(ncid, varid, attname, v);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        _ncw_inq_varname(ncid, varid, varname);
        quit("\"%s\": ncw_get_att_double(): failed for varid = %d (varname = \"%s\"), attname = \"%s\": %s\n", _ncw_get_path(ncid), varid, varname, attname, nc_strerror(status));
    }
}

/*
 * The following procedures do not have direct analogues in the netcdf library.
 */

int ncw_inq_nrecords(int ncid)
{
    int unlimdimid;
    size_t nrecords;

    ncw_inq_unlimdimid(ncid, &unlimdimid);
    if (unlimdimid < 0)
        return 0;
    ncw_inq_dimlen(ncid, unlimdimid, &nrecords);

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
 * @param ncid_src Source file id
 * @param ncid_dst Destination file id
 */
void ncw_copy_dims(int ncid_src, int ncid_dst)
{
    int ndims;
    int unlimdimid = -1;
    int i;

    ncw_inq_ndims(ncid_src, &ndims);
    ncw_inq_unlimdimid(ncid_src, &unlimdimid);

    for (i = 0; i < ndims; ++i) {
        char dimname[NC_MAX_NAME] = STR_UNKNOWN;
        size_t size;
        int dimid;

        memset(dimname, 0, NC_MAX_NAME);

        ncw_inq_dim(ncid_src, i, dimname, &size);
        if (i == unlimdimid)
            size = NC_UNLIMITED;
        if (!ncw_dim_exists(ncid_dst, dimname))
            ncw_def_dim(ncid_dst, dimname, size, &dimid);
    }
}

/** Copies dimension from one NetCDF file to another.
 *
 * @param ncid_src Source file id
 * @param dimname Dimension name
 * @param ncid_dst Destination file id
 */
void ncw_copy_dim(int ncid_src, const char dimname[], int ncid_dst)
{
    int unlimdimid = -1;
    int dimid;
    size_t size;

    ncw_inq_unlimdimid(ncid_src, &unlimdimid);
    ncw_inq_dimid(ncid_src, dimname, &dimid);
    if (dimid == unlimdimid)
        size = NC_UNLIMITED;
    else
        ncw_inq_dimlen(ncid_src, dimid, &size);
    ncw_def_dim(ncid_dst, dimname, size, &dimid);
}

/** Copies definition of a specified variable from one NetCDF file to another.
 *
 * @param ncid_src Source file id
 * @param vid_src Variable id
 * @param ncid_dst Destination file id
 * @return Variable id in the destination file
 */
int ncw_copy_vardef(int ncid_src, int vid_src, int ncid_dst)
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

    ncw_inq_varname(ncid_src, vid_src, varname);
    ncw_inq_var(ncid_src, vid_src, NULL, &type, &ndims, dimids_src, &natts);

    for (i = 0; i < ndims; ++i) {
        char dimname[NC_MAX_NAME];
        size_t len;

        ncw_inq_dim(ncid_src, dimids_src[i], dimname, &len);
        if (!ncw_dim_exists(ncid_dst, dimname)) {
            int unlimdimid = -1;

            ncw_inq_unlimdimid(ncid_src, &unlimdimid);
            ncw_def_dim(ncid_dst, dimname, (dimids_src[i] == unlimdimid) ? NC_UNLIMITED : len, &dimids_dst[i]);
        } else
            ncw_inq_dimid(ncid_dst, dimname, &dimids_dst[i]);
    }

    ncw_def_var(ncid_dst, varname, type, ndims, dimids_dst, &vid_dst);
    ncw_copy_atts(ncid_src, vid_src, ncid_dst, vid_dst);

    if (status == NC_NOERR)
        nc_enddef(ncid_dst);

    return vid_dst;
}

/** Copies data for a given variable from one file to another.
 *  Note that variable IDs may be different, but name and dimensions must be
 *  the same.
 *
 * @param ncid_src Source file id
 * @param vid_src Variable id
 * @param ncid_dst Destination file id
 */
void ncw_copy_vardata(int ncid_src, int vid_src, int ncid_dst)
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

    ncw_inq_varname(ncid_src, vid_src, varname);
    ncw_inq_var(ncid_src, vid_src, NULL, &type, &ndims, dimids, &natts);
    for (i = 0; i < ndims; ++i)
        ncw_inq_dimlen(ncid_src, dimids[i], &dimlens[i]);
    size = 1;
    for (i = 0; i < ndims; ++i)
        size *= dimlens[i];
    ncw_inq_varid(ncid_dst, varname, &vid_dst);

    size *= ncw_sizeof(type);
    data = malloc(size);
    switch (type) {
    case NC_BYTE:
    case NC_CHAR:
        ncw_get_var_text(ncid_src, vid_src, data);
        ncw_put_var_text(ncid_dst, vid_dst, data);
        break;
    case NC_SHORT:
        ncw_get_var_short(ncid_src, vid_src, data);
        ncw_put_var_short(ncid_dst, vid_dst, data);
        break;
    case NC_INT:
        ncw_get_var_int(ncid_src, vid_src, data);
        ncw_put_var_int(ncid_dst, vid_dst, data);
        break;
    case NC_FLOAT:
        ncw_get_var_float(ncid_src, vid_src, data);
        ncw_put_var_float(ncid_dst, vid_dst, data);
        break;
    case NC_DOUBLE:
        ncw_get_var_double(ncid_src, vid_src, data);
        ncw_put_var_double(ncid_dst, vid_dst, data);
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
 * @param ncid_src Source file id
 * @param varname Variable name
 * @param ncid_dst Destination file id
 */
void ncw_copy_var(int ncid_src, const char varname[], int ncid_dst)
{
    int vid_src = -1;

    ncw_inq_varid(ncid_src, varname, &vid_src);
    (void) ncw_copy_vardef(ncid_src, vid_src, ncid_dst);
    ncw_copy_vardata(ncid_src, vid_src, ncid_dst);
}

/** Set deflation (compression) level for all variables in a file.
 *
 * @param ncid NetCDF file id
 * @param shuffle Shuffle flag, whatever it means
 * @param deflate Flag, turns deflation on/off
 * @param deflate_level Compression level
 */
void ncw_def_deflate(int ncid, int shuffle, int deflate, int deflate_level)
{
    int nv = -1;                /* total number of variables */
    int vid;

    ncw_inq_nvars(ncid, &nv);
    for (vid = 0; vid < nv; ++vid) {
        int status = nc_def_var_deflate(ncid, vid, shuffle, deflate, deflate_level);

        if (status != NC_NOERR) {
            char varname[NC_MAX_NAME] = "STR_UNKNOWN";

            _ncw_inq_varname(ncid, vid, varname);
            quit("\"%s\": ncw_def_deflate(): failed for \"%s\": %s\n", _ncw_get_path(ncid), varname, nc_strerror(status));
        }
    }
}

/** Gets the id for the first dimension found to be present in a NetCDF file
 * out of two dimensions specified by names. Useful for handling both new and
 * old data formats.
 *
 * @param ncid NetCDF file id
 * @param name1 The first dimension name to be tried
 * @param name2 The second dimension name to be tried
 * @param dimid Dimension id (output)
 */
void ncw_inq_dimid2(int ncid, const char dimname1[], const char dimname2[], int* dimid)
{
    int status1 = nc_inq_dimid(ncid, dimname1, dimid);

    if (status1 != NC_NOERR) {
        int status2 = nc_inq_dimid(ncid, dimname2, dimid);

        if (status2 != NC_NOERR) {
            quit("\"%s\": nc_inq_dimid(): failed for dimname = \"%s\": %s, and for dimname = \"%s\": %s\n", _ncw_get_path(ncid), dimname1, nc_strerror(status1), dimname2, nc_strerror(status2));
        }
    }
}

/** Gets the value(s) (converted to int type) of the first attribute found to
 * be present in a NetCDF file out of two specified by attribute names. Useful
 * for handling both new and old data formats.
 *
 * @param ncid NetCDF file id
 * @param varid Name of the variable the attrubute belongs to (use NC_GLOBAL
 *              for global attributes)
 * @param attname1 The first attribute name to be tried
 * @param attname2 The second attribute name to be tried
 * @param in Attribute value(s) (output)
 */
void ncw_get_att_int2(int ncid, int varid, const char attname1[], const char attname2[], int v[])
{
    nc_type xtype;
    size_t len;
    int status1 = nc_get_att_int(ncid, varid, attname1, v);

    if (status1 != NC_NOERR) {
        int status2 = nc_get_att_int(ncid, varid, attname2, v);

        if (status2 != NC_NOERR) {
            char varname[NC_MAX_NAME] = STR_UNKNOWN;

            _ncw_inq_varname(ncid, varid, varname);
            quit("\"%s\": nc_get_att_int(): failed for varid = %d (varname = \"%s\"), attname = \"%s\": %s, and attname = \"%s\": %s\n", _ncw_get_path(ncid), varid, varname, attname1, nc_strerror(status1), attname2, nc_strerror(status2));
        }
    }

    ncw_inq_att(ncid, varid, (status1 == NC_NOERR) ? attname1 : attname2, &xtype, &len);

    if (xtype != NC_INT && xtype != NC_BYTE && xtype != NC_SHORT) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        _ncw_inq_varname(ncid, varid, varname);
        quit("ncw_get_att_int2(): failed for varid = %d (varname = \"%s\"), attname = \"%s\": the attribute is of %s type\n", _ncw_get_path(ncid), varid, varname, (status1 == NC_NOERR) ? attname1 : attname2, ncw_nctype2str(xtype));
    }
}

/** Finds all variables with specified dimensions and specified attribute.
 *
 * @param ncid NetCDF file id
 * @param ndims Number of dimensions (0 if any)
 * @param dims Array of dimension ids
 * @param attname Attribute name that is be set (NULL if any)
 * @param attval Attribute value (NULL if any)
 * @param nvars Return -- number of variables found
 * @param vids Return -- array of variable ids (needs to be deallocated when
 *             necessary)
 */
void ncw_find_vars(int ncid, int ndims, const int dims[], const char attname[], const void* attval, int* nvars, int** vids)
{
    int nv = -1;                /* total number of variables */
    int nallocated = 0;
    int vid;

    *nvars = 0;
    *vids = NULL;

    ncw_inq_nvars(ncid, &nv);

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
            ncw_inq_var(ncid, vid, varname, &vtype, &nd, dids, &natts);

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

                ncw_inq_attname(ncid, vid, i, aname);

                if (strcmp(aname, attname) != 0)
                    continue;

                if (attval != NULL) {
                    size_t alen = UINT_MAX;
                    nc_type atype = -1;
                    char* aval = NULL;

                    ncw_inq_att(ncid, vid, aname, &atype, &alen);
                    if (alen <= 0)
                        continue;
                    aval = calloc(alen, ncw_sizeof(atype));
                    if (aval == NULL)
                        quit("\"%s\": ncw_find_vars(): could not allocate memory for attribute = \"%s\", type = %s, length = %d, varid = %d (varname = \"%s\"): %s\n", _ncw_get_path(ncid), aname, ncw_nctype2str(atype), alen, vid, (vid == NC_GLOBAL) ? "NC_GLOBAL" : varname, strerror(errno));
                    ncw_get_att_text(ncid, vid, attname, aval);
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
 * @param ncid NetCDF file id
 * @param varid Return -- id of the variable
 */
void ncw_find_timevarid(int ncid, int* varid)
{
    int unlimdimid = -1;
    int nvars = -1;
    int i;

    ncw_inq_unlimdimid(ncid, &unlimdimid);
    ncw_inq_nvars(ncid, &nvars);

    for (i = 0; i < nvars; ++i) {
        int ndims = -1;
        int dimid = -1;

        ncw_inq_varndims(ncid, i, &ndims);
        if (ndims != 1)
            continue;
        ncw_inq_vardimid(ncid, i, &dimid);
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
 * @param ncid_src Source file id
 * @param varid Variable id
 * @param ncid_dst Destination file id
 */
void ncw_copy_atts(int ncid_src, int varid_src, int ncid_dst, int varid_dst)
{
    char varname_src[NC_MAX_NAME] = STR_UNKNOWN;
    char varname_dst[NC_MAX_NAME] = STR_UNKNOWN;
    int status;
    int natts;
    int i;

    _ncw_inq_varname(ncid_src, varid_src, varname_src);
    _ncw_inq_varname(ncid_dst, varid_dst, varname_dst);
    ncw_inq_varnatts(ncid_src, varid_src, &natts);

    for (i = 0; i < natts; ++i) {
        char attname[NC_MAX_NAME] = STR_UNKNOWN;

        ncw_inq_attname(ncid_src, varid_src, i, attname);
        if ((status = nc_copy_att(ncid_src, varid_src, attname, ncid_dst, varid_dst)) != NC_NOERR)
            quit("\"%s\": -> %s:  nc_copy_att(): failed for varid_in = %d (varname = \"%s\"), attname = \"%s\", varid_dst = %d, (varname = \"%s\"): %s\n", _ncw_get_path(ncid_src), _ncw_get_path(ncid_dst), varid_src, varname_src, attname, varid_dst, varname_dst, nc_strerror(status));
    }
}

/** Defines a new variable using an existing variable as a template.
 *
 * @param ncid NetCDF file id
 * @param oldvarname Name of the existing variable
 * @param newvarname Name of the new variable
 */
void ncw_def_var_as(int ncid, const char oldvarname[], const char newvarname[])
{
    int oldvarid, newvarid;
    nc_type type;
    int ndims;
    int dimids[NC_MAX_DIMS];
    int natts;
    int i;

    ncw_inq_varid(ncid, oldvarname, &oldvarid);
    ncw_inq_var(ncid, oldvarid, NULL, &type, &ndims, dimids, &natts);

    ncw_def_var(ncid, newvarname, type, ndims, dimids, &newvarid);

    for (i = 0; i < natts; ++i) {
        char attname[NC_MAX_NAME] = STR_UNKNOWN;

        ncw_inq_attname(ncid, oldvarid, i, attname);
        ncw_copy_att(ncid, oldvarid, attname, ncid, newvarid);
    }
}

/** Reads one record of a variable.
 * @param ncid NetCDF file id
 * @param varid ID of the variable
 * @param r Record number
 * @param v The data
 */
void ncw_get_var_double_record(int ncid, int varid, int r, double v[])
{
    int ndims;
    int dimids[NC_MAX_DIMS];
    size_t dimlen[NC_MAX_VAR_DIMS];
    size_t start[NC_MAX_VAR_DIMS];
    int i;
    int status;

    ncw_inq_varndims(ncid, varid, &ndims);
    ncw_inq_vardimid(ncid, varid, dimids);
    start[0] = r;               /* this record */
    dimlen[0] = 1;              /* one record only */
    for (i = 1; i < ndims; ++i) {
        ncw_inq_dimlen(ncid, dimids[i], &dimlen[i]);
        start[i] = 0;
    }
    status = nc_get_vara_double(ncid, varid, start, dimlen, v);

    if (status != 0) {
        char varname[NC_MAX_NAME] = "STR_UNKNOWN";

        ncw_inq_varname(ncid, varid, varname);
        quit("\"%s\": ncw_get_var_double_record(): failed to read record (outer dim) %d for \"%s\": %s\n", _ncw_get_path(ncid), r, varname, nc_strerror(status));
    }
}

/** Reads one record of a variable.
 * @param ncid NetCDF file id
 * @param varid ID of the variable
 * @param r Record number
 * @param v The data
 */
void ncw_get_var_float_record(int ncid, int varid, int r, float v[])
{
    int ndims;
    int dimids[NC_MAX_DIMS];
    size_t dimlen[NC_MAX_VAR_DIMS];
    size_t start[NC_MAX_VAR_DIMS];
    int i;
    int status;

    ncw_inq_varndims(ncid, varid, &ndims);
    ncw_inq_vardimid(ncid, varid, dimids);
    start[0] = r;               /* this record */
    dimlen[0] = 1;              /* one record only */
    for (i = 1; i < ndims; ++i) {
        ncw_inq_dimlen(ncid, dimids[i], &dimlen[i]);
        start[i] = 0;
    }
    status = nc_get_vara_float(ncid, varid, start, dimlen, v);

    if (status != 0) {
        char varname[NC_MAX_NAME] = "STR_UNKNOWN";

        ncw_inq_varname(ncid, varid, varname);
        quit("\"%s\": ncw_get_var_float_record(): failed to read record (outer dim) %d for \"%s\": %s\n", _ncw_get_path(ncid), r, varname, nc_strerror(status));
    }
}

/** Writes one record of a variable.
 * @param ncid NetCDF file id
 * @param varid ID of the variable
 * @param r Record number
 * @param v The data
 */
void ncw_put_var_double_record(int ncid, int varid, int r, double v[])
{
    int ndims;
    int dimids[NC_MAX_DIMS];
    size_t dimlen[NC_MAX_VAR_DIMS];
    size_t start[NC_MAX_VAR_DIMS];
    int i;
    int status;

    ncw_inq_varndims(ncid, varid, &ndims);
    ncw_inq_vardimid(ncid, varid, dimids);
    start[0] = r;               /* this record */
    dimlen[0] = 1;              /* one record only */
    for (i = 1; i < ndims; ++i) {
        ncw_inq_dimlen(ncid, dimids[i], &dimlen[i]);
        start[i] = 0;
    }
    status = nc_put_vara_double(ncid, varid, start, dimlen, v);

    if (status != 0) {
        char varname[NC_MAX_NAME] = "STR_UNKNOWN";

        ncw_inq_varname(ncid, varid, varname);
        quit("\"%s\": ncw_put_var_double_record(): failed to read record (outer dim) %d for \"%s\": %s\n", _ncw_get_path(ncid), r, varname, nc_strerror(status));
    }
}

/** Writes one record of a variable.
 * @param ncid NetCDF file id
 * @param varid ID of the variable
 * @param r Record number
 * @param v The data
 */
void ncw_put_var_float_record(int ncid, int varid, int r, float v[])
{
    int ndims;
    int dimids[NC_MAX_DIMS];
    size_t dimlen[NC_MAX_VAR_DIMS];
    size_t start[NC_MAX_VAR_DIMS];
    int i;
    int status;

    ncw_inq_varndims(ncid, varid, &ndims);
    ncw_inq_vardimid(ncid, varid, dimids);
    start[0] = r;               /* this record */
    dimlen[0] = 1;              /* one record only */
    for (i = 1; i < ndims; ++i) {
        ncw_inq_dimlen(ncid, dimids[i], &dimlen[i]);
        start[i] = 0;
    }
    status = nc_put_vara_float(ncid, varid, start, dimlen, v);

    if (status != 0) {
        char varname[NC_MAX_NAME] = "STR_UNKNOWN";

        ncw_inq_varname(ncid, varid, varname);
        quit("\"%s\": ncw_put_var_float_record(): failed to read record (outer dim) %d for \"%s\": %s\n", _ncw_get_path(ncid), r, varname, nc_strerror(status));
    }
}

/** Checks that the attribute has certain type and length
 */
void ncw_check_att(int ncid, int varid, const char attname[], size_t att_len)
{
    nc_type type;
    size_t len;
    int status = nc_inq_att(ncid, varid, attname, &type, &len);

    if (status != NC_NOERR) {
        char varname[NC_MAX_NAME] = STR_UNKNOWN;

        _ncw_inq_varname(ncid, varid, varname);
        quit("\"%s\": nc_inq_att(): failed for varid = %d (varname = \"%s\"), attname = \"%s\", atttype = %d: %s\n", _ncw_get_path(ncid), varid, varname, attname, ncw_nctype2str(type), nc_strerror(status));
    }

    if (len != att_len)
        quit("\"%s\": ncw_check_att(): attribute \"%s\" is supposed to have length %z; the actual length is %z\n", _ncw_get_path(ncid), attname, att_len, len);
}

/** Check that the dimension has certain length
 */
void ncw_check_dimlen(int ncid, const char dimname[], size_t dimlen)
{
    size_t len;
    int dimid;
    int status;

    ncw_inq_dimid(ncid, dimname, &dimid);
    status = nc_inq_dimlen(ncid, dimid, &len);

    if (status != NC_NOERR)
        quit("\"%s\": nc_inq_dimlen(): failed for dimid = %d (dimname = \"%s\"): %s\n", _ncw_get_path(ncid), dimid, dimname, nc_strerror(status));
    if (len != dimlen)
        quit("\"%s\": ncw_check_dimlen(): dimension \"%s\" is supposed to have length %zu; the actual length is %zu\n", _ncw_get_path(ncid), dimname, dimlen, len);
}
