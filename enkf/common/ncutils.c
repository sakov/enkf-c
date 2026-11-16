/******************************************************************************
 *
 * File:        ncutils.c        
 *
 * Created:     6/9/2019
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description: Intermediate level NetCDF r/w procedures. Among other things
 *              these procedures are supposed to handle the following
 *              attributes: _FillValue, missing_value, valid_range, valid_min,
 *              valid_max, add_offset, and scale_factor.
 *
 * Revisions:   PS 20250531 Addded argument "r" (record) to ncu_readfield().
 *
 *****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include <stdint.h>
#include <assert.h>
#include "ncw.h"
#include "ncutils.h"

static void ncu_quit_fn_def(char* format, ...);
static ncu_quit_fn quit = ncu_quit_fn_def;

/**
 */
static void ncu_quit_fn_def(char* format, ...)
{
    va_list args;

    fflush(stdout);

    fprintf(stderr, "\n\n  error: ncu: ");
    va_start(args, format);
    vfprintf(stderr, format, args);
    va_end(args);
    fprintf(stderr, "\n\n");
    exit(1);
}

/**
 */
void ncu_set_quitfn(ncu_quit_fn quitfn)
{
    quit = quitfn;
}

/** Converts a value of the specified NetCDF type to double.
 */
static double ncu_todouble(nc_type type, const void* p)
{
    switch (type) {
    case NC_BYTE:
        return (double) ((const int8_t*) p)[0];
    case NC_CHAR:
        return (double) ((const uint8_t *) p)[0];
    case NC_UBYTE:
        return (double) ((const uint8_t *) p)[0];
    case NC_SHORT:
        return (double) ((const int16_t*) p)[0];
    case NC_USHORT:
        return (double) ((const uint16_t*) p)[0];
    case NC_INT:
        return (double) ((const int32_t*) p)[0];
    case NC_UINT:
        return (double) ((const uint32_t*) p)[0];
    case NC_INT64:
        return (double) ((const int64_t*) p)[0];
    case NC_UINT64:
        return (double) ((const uint64_t *) p)[0];
    case NC_FLOAT:
        return (double) ((const float*) p)[0];
    case NC_DOUBLE:
        return (double) ((const double*) p)[0];
    default:
        quit("ncu_todouble(): unexpected NetCDF type %d", (int) type);
    }
    return NAN;
}

typedef union {
    int8_t i8;
    uint8_t u8;
    int16_t i16;
    uint16_t u16;
    int32_t i32;
    uint32_t u32;
    int64_t i64;
    uint64_t u64;
    float f32;
    double f64;
} ncu_val;

/** Applies _FillValue, missing_value, valid_min, valid_max, and valid_range
 ** masks from raw typed buffer `vv` onto either a float[] or double[] array.
 */
static void ncu_apply_masks(int ncid, int varid, nc_type vartype, size_t n, const void* vv, void* out, int is_double)
{
    int has_fill = 0, has_miss = 0, has_vmin = 0, has_vmax = 0, has_vrange = 0;
    ncu_val fill_buf = { 0 }, miss_buf = {
    0};
    ncu_val vmin_buf = { 0 }, vmax_buf = {
    0}, vrange_buf[2] = { {
    0}, {
    0}};

    if (ncw_att_exists2(ncid, varid, "_FillValue")) {
        ncw_check_attlen(ncid, varid, "_FillValue", 1);
        ncw_get_att(ncid, varid, "_FillValue", &fill_buf);
        has_fill = 1;
    } else {
        int nofill = 1;

        ncw_inq_var_fill(ncid, varid, &nofill, &fill_buf);
        has_fill = !nofill;
    }

    if (ncw_att_exists2(ncid, varid, "missing_value")) {
        ncw_check_attlen(ncid, varid, "missing_value", 1);
        ncw_get_att(ncid, varid, "missing_value", &miss_buf);
        has_miss = 1;
    }
    if (ncw_att_exists2(ncid, varid, "valid_min")) {
        ncw_check_attlen(ncid, varid, "valid_min", 1);
        ncw_get_att(ncid, varid, "valid_min", &vmin_buf);
        has_vmin = 1;
    }
    if (ncw_att_exists2(ncid, varid, "valid_max")) {
        ncw_check_attlen(ncid, varid, "valid_max", 1);
        ncw_get_att(ncid, varid, "valid_max", &vmax_buf);
        has_vmax = 1;
    }
    if (ncw_att_exists2(ncid, varid, "valid_range")) {
        ncw_check_attlen(ncid, varid, "valid_range", 2);
        ncw_get_att(ncid, varid, "valid_range", vrange_buf);
        has_vrange = 1;
    }

    if (!has_fill && !has_miss && !has_vmin && !has_vmax && !has_vrange)
        return;

#define APPLY_MASK_LOOP(CTYPE, FIELD) do {                                     \
        const CTYPE* raw = (const CTYPE*) vv;                                  \
        CTYPE fill   = fill_buf.FIELD;                                         \
        CTYPE miss   = miss_buf.FIELD;                                         \
        CTYPE vmin   = vmin_buf.FIELD;                                         \
        CTYPE vmax   = vmax_buf.FIELD;                                         \
        CTYPE vr_min = vrange_buf[0].FIELD;                                    \
        CTYPE vr_max = vrange_buf[1].FIELD;                                    \
        for (size_t i = 0; i < n; ++i) {                                       \
            CTYPE val = raw[i];                                                \
            int invalid = (has_fill   && val == fill) ||                       \
                          (has_miss   && val == miss) ||                       \
                          (has_vmin   && val < vmin)  ||                       \
                          (has_vmax   && val > vmax)  ||                       \
                          (has_vrange && (val < vr_min || val > vr_max));      \
            if (invalid) {                                                     \
                if (is_double) ((double*) out)[i] = NAN;                       \
                else           ((float*)  out)[i] = NAN;                       \
            }                                                                  \
        }                                                                      \
    } while (0)

    switch (vartype) {
    case NC_BYTE:
        APPLY_MASK_LOOP(int8_t, i8);

        break;
    case NC_CHAR:
    case NC_UBYTE:
        APPLY_MASK_LOOP(uint8_t, u8);
        break;
    case NC_SHORT:
        APPLY_MASK_LOOP(int16_t, i16);

        break;
    case NC_USHORT:
        APPLY_MASK_LOOP(uint16_t, u16);

        break;
    case NC_INT:
        APPLY_MASK_LOOP(int32_t, i32);

        break;
    case NC_UINT:
        APPLY_MASK_LOOP(uint32_t, u32);

        break;
    case NC_INT64:
        APPLY_MASK_LOOP(int64_t, i64);

        break;
    case NC_UINT64:
        APPLY_MASK_LOOP(uint64_t, u64);
        break;
    case NC_FLOAT:
        APPLY_MASK_LOOP(float, f32);

        break;
    case NC_DOUBLE:
        APPLY_MASK_LOOP(double, f64);

        break;
    default:
        quit("ncu_apply_masks(): unexpected NetCDF type %d", (int) vartype);
    }

#undef APPLY_MASK_LOOP
}

/** Allocates raw buffer if needed, reads raw values, applies masks, and unpacks
 ** scale_factor / add_offset onto an already-populated float[] or double[]
 ** array.
 */
static void ncu_postprocess_read(int ncid, int varid, const size_t start[], const size_t count[], size_t n, void* out, int is_double, int is_full_var)
{
    nc_type vartype;
    nc_type target_type = is_double ? NC_DOUBLE : NC_FLOAT;
    int needs_raw = ncw_att_exists2(ncid, varid, "_FillValue") || ncw_att_exists2(ncid, varid, "missing_value") || ncw_att_exists2(ncid, varid, "valid_range") || ncw_att_exists2(ncid, varid, "valid_min") || ncw_att_exists2(ncid, varid, "valid_max");

    ncw_inq_vartype(ncid, varid, &vartype);

    /*
     * Check implicit NetCDF fill value even if _FillValue attribute is
     * absent 
     */
    if (!needs_raw && !ncw_att_exists2(ncid, varid, "_FillValue")) {
        int nofill = 1;
        ncu_val dummy[2];

        ncw_inq_var_fill(ncid, varid, &nofill, dummy);
        if (!nofill)
            needs_raw = 1;
    }

    if (needs_raw) {
        void* vv = out;

        if (vartype != target_type) {
            int typesize = ncw_sizeof(vartype);

            vv = malloc(n * (size_t) typesize);
            if (vv == NULL)
                quit("ncutils: could not allocate memory (%zu bytes)", n * (size_t) typesize);
            if (is_full_var)
                ncw_get_var(ncid, varid, vv);
            else
                ncw_get_vara(ncid, varid, start, count, vv);
        }

        ncu_apply_masks(ncid, varid, vartype, n, vv, out, is_double);

        if (vv != out)
            free(vv);
    }

    if (ncw_att_exists(ncid, varid, "scale_factor")) {
        ncw_check_attlen(ncid, varid, "scale_factor", 1);
        if (is_double) {
            double sf;

            ncw_get_att_double(ncid, varid, "scale_factor", &sf);
            for (size_t i = 0; i < n; ++i)
                ((double*) out)[i] *= sf;
        } else {
            float sf;

            ncw_get_att_float(ncid, varid, "scale_factor", &sf);
            for (size_t i = 0; i < n; ++i)
                ((float*) out)[i] *= sf;
        }
    }

    if (ncw_att_exists(ncid, varid, "add_offset")) {
        ncw_check_attlen(ncid, varid, "add_offset", 1);
        if (is_double) {
            double ao;

            ncw_get_att_double(ncid, varid, "add_offset", &ao);
            for (size_t i = 0; i < n; ++i)
                ((double*) out)[i] += ao;
        } else {
            float ao;

            ncw_get_att_float(ncid, varid, "add_offset", &ao);
            for (size_t i = 0; i < n; ++i)
                ((float*) out)[i] += ao;
        }
    }
}

/** Rounds values to be written to a NetCDF variable of an integer type
 ** (otherwise the NetCDF library truncates them).
 */
static void ncu_roundforint(nc_type type, size_t n, float* v)
{
    size_t i;

    if (type == NC_FLOAT || type == NC_DOUBLE || type == NC_CHAR)
        return;
    for (i = 0; i < n; ++i)
        v[i] = roundf(v[i]);
}

/** Applies inverse add_offset/scale_factor packing, rounds for integer types,
 ** and clamps to valid_min/valid_max/valid_range before writing.
 ** Leaves NaNs intact so ncu_put_vara() can write exact fill values.
 */
static void ncu_preprocess_write(int ncid, int varid, size_t n, float* v)
{
    nc_type vartype;
    float attval[2];
    size_t i;

    if (ncw_att_exists(ncid, varid, "add_offset")) {
        float add_offset;

        ncw_check_attlen(ncid, varid, "add_offset", 1);
        ncw_get_att_float(ncid, varid, "add_offset", &add_offset);
        for (i = 0; i < n; ++i)
            v[i] -= add_offset;
    }

    if (ncw_att_exists(ncid, varid, "scale_factor")) {
        float scale_factor;

        ncw_check_attlen(ncid, varid, "scale_factor", 1);
        ncw_get_att_float(ncid, varid, "scale_factor", &scale_factor);
        for (i = 0; i < n; ++i)
            v[i] /= scale_factor;
    }

    ncw_inq_vartype(ncid, varid, &vartype);
    ncu_roundforint(vartype, n, v);

    if (ncw_att_exists2(ncid, varid, "valid_min")) {
        ncw_check_attlen(ncid, varid, "valid_min", 1);
        ncw_get_att_float(ncid, varid, "valid_min", attval);
        for (i = 0; i < n; ++i)
            if (v[i] < attval[0])
                v[i] = attval[0];
    }
    if (ncw_att_exists2(ncid, varid, "valid_max")) {
        ncw_check_attlen(ncid, varid, "valid_max", 1);
        ncw_get_att_float(ncid, varid, "valid_max", attval);
        for (i = 0; i < n; ++i)
            if (v[i] > attval[0])
                v[i] = attval[0];
    }
    if (ncw_att_exists2(ncid, varid, "valid_range")) {
        ncw_check_attlen(ncid, varid, "valid_range", 2);
        ncw_get_att_float(ncid, varid, "valid_range", attval);
        for (i = 0; i < n; ++i)
            if (v[i] < attval[0])
                v[i] = attval[0];
            else if (v[i] > attval[1])
                v[i] = attval[1];
    }
}

/** Resolves the fill value for a variable (checking _FillValue, implicit
 ** NetCDF fill, and missing_value in order). Returns NAN if none is set.
 */
static double ncu_get_write_fill(int ncid, int varid, nc_type vartype)
{
    if (ncw_att_exists2(ncid, varid, "_FillValue")) {
        double fill = NAN;

        ncw_check_attlen(ncid, varid, "_FillValue", 1);
        ncw_get_att_double(ncid, varid, "_FillValue", &fill);
        return fill;
    } else {
        int nofill = 1;
        ncu_val buf[2];

        ncw_inq_var_fill(ncid, varid, &nofill, buf);
        if (!nofill)
            return ncu_todouble(vartype, buf);
    }

    if (ncw_att_exists2(ncid, varid, "missing_value")) {
        double miss = NAN;

        ncw_check_attlen(ncid, varid, "missing_value", 1);
        ncw_get_att_double(ncid, varid, "missing_value", &miss);
        return miss;
    }

    return NAN;
}

/** Writes a float array to a NetCDF variable, replacing NaNs with the exact
 ** fill value. Uses double for 32/64-bit integer types so large integer fill
 ** values are never rounded by float precision.
 */
static void ncu_put_vara(int ncid, int varid, const size_t start[], const size_t count[], size_t n, float* v)
{
    nc_type vartype;
    double fill;
    size_t i;

    ncw_inq_vartype(ncid, varid, &vartype);
    fill = ncu_get_write_fill(ncid, varid, vartype);

    if (vartype == NC_INT || vartype == NC_UINT || vartype == NC_INT64 || vartype == NC_UINT64) {
        double* dv = malloc(n * sizeof(double));

        if (dv == NULL)
            quit("ncu_put_vara(): could not allocate memory (%zu bytes)", n * sizeof(double));
        for (i = 0; i < n; ++i)
            dv[i] = isnan(v[i]) ? fill : (double) v[i];
        ncw_put_vara_double(ncid, varid, start, count, dv);
        free(dv);
    } else {
        if (!isnan(fill)) {
            float ffill = (float) fill;

            for (i = 0; i < n; ++i)
                if (isnan(v[i]))
                    v[i] = ffill;
        }
        ncw_put_vara_float(ncid, varid, start, count, v);
    }
}

/**
 */
int ncu_getnlevels(char fname[], char varname[], int isstructured)
{
    int ncid;
    int varid;
    int ndims;
    size_t dimlen[4];
    int hasrecorddim;

    ncw_open(fname, NC_NOWRITE, &ncid);
    ncw_inq_varid(ncid, varname, &varid);
    ncw_inq_vardims(ncid, varid, 4, &ndims, dimlen);
    hasrecorddim = (ncw_var_hasunlimdim(ncid, varid) || dimlen[0] == 1);
    ncw_close(ncid);

    if (isstructured) {
        if (ndims > 4)
            quit("%s: %s: do not know how to handle more than 4-dimensional variables\n", fname, varname);
        if (ndims == 4) {
            if (!hasrecorddim)
                quit("%s: %s: expect an unlimited dimension to be present for a 4-dimensional variable\n", fname, varname);
            return (int) dimlen[1];
        }
        if (ndims == 3)
            return (hasrecorddim) ? 1 : (int) dimlen[0];
        if (ndims == 2)
            return (hasrecorddim) ? 0 : 1;
    } else {
        if (ndims > 3)
            quit("%s: %s: EnKF-C does not know how to handle more than 3-dimensional variables on unstructured grids\n", fname, varname);
        if (ndims == 3) {
            if (!hasrecorddim)
                quit("%s: %s: expect an unlimited dimension to be present for a 3-dimensional variable on unstructured grid\n", fname, varname);
            return (int) dimlen[1];
        }
        if (ndims == 2)
            return (hasrecorddim) ? 1 : (int) dimlen[0];
        if (ndims == 1)
            return (hasrecorddim) ? 0 : 1;
    }

    return 0;
}

/** Reads one horizontal field (layer) for a variable from a NetCDF file.
 ** Verifies that the field dimensions are ni x nj.
 */
void ncu_readfield(char fname[], char varname[], int r, int k, int ni, int nj, int nk, float* v)
{
    ncu_readfield_part(fname, varname, r, k, 0, -1, ni, nj, nk, v);
}

/** Reads part of one horizontal field (layer) from a NetCDF file.
 */
void ncu_readfield_part(char fname[], char varname[], int r, int k, int j1, int j2, int ni, int nj, int nk, float* v)
{
    int ncid;
    int varid;
    int ndims;
    size_t dimlen[4];
    size_t start[4], count[4];
    size_t i, n;
    int hasrecorddim;

    if (j1 <= 0 && j2 < 0) {
        /*
         * read full field
         */
        if (nj > 0) {
            j1 = 0;
            j2 = nj - 1;
        } else {
            j1 = 0;
            j2 = ni - 1;
        }
    }

    if (j2 < j1)
        quit("ncu_readfield_part(): %s: %s: impossible range request: j1 = %d, j2 = %d\n", fname, varname, j1, j2);
    if (j1 < 0 || (nj > 0 && j2 > nj - 1))
        quit("ncu_readfield_part(): %s: %s: impossible range request: j1 = %d, j2 = %d; nj = %d\n", fname, varname, j1, j2, nj);
    else if (nj <= 0 && j2 > ni - 1)
        quit("ncu_readfield_part(): %s: %s: impossible range request: i1 = %d, ji = %d; ni = %d\n", fname, varname, j1, j2, ni);

    ncw_open(fname, NC_NOWRITE, &ncid);
    ncw_inq_varid(ncid, varname, &varid);
    ncw_inq_vardims(ncid, varid, 4, &ndims, dimlen);
    hasrecorddim = ncw_var_hasunlimdim(ncid, varid);
    if (hasrecorddim && dimlen[0] == 0)
        quit("ncu_readfield_part(): \"%s\": %s: empty record dimension", fname, varname);
    if (hasrecorddim && r < 0) {
        if (dimlen[0] > 1)
            quit("ncu_readfield_part(): \"%s\": %s: not supposed to have multiple records", fname, varname);
        else
            r = 0;
    }

    if (ndims == 4) {
        /*
         * the first dimension is considered as record dimension, regardless of
         * whether it is unlimited or not
         */
        if (r < 0) {
            if (dimlen[0] > 1)
                quit("ncu_readfield_part(): \"%s\": %s: not supposed to have multiple records", fname, varname);
            r = 0;
        }
        if (r >= dimlen[0])
            quit("ncu_readfield_part(): \"%s\": %s: can not read record %d (nr = %zu)", fname, varname, r, dimlen[0]);
        if (nj == 0)            /* unstructured grid */
            quit("ncu_readfield_part(): \"%s\": %s: expected positive \"j\" dimension for a 4-dimensional variable", fname, varname);
        start[0] = r;
        if (dimlen[1] != nk) {
            if (dimlen[1] != 1)
                quit("ncu_readfield_part(): \"%s\": vertical dimension of variable \"%s\" (nk = %zu) does not match grid dimension (nk = %d)", fname, varname, dimlen[1], nk);
            else
                k = 0;          /* ignore k */
        }
        if (dimlen[3] != ni || dimlen[2] != nj)
            quit("ncu_readfield_part(): %s: %s: horizontal dimensions (ni = %zu, nj = %zu) do not match grid dimensions (ni = %d, nj = %d)", fname, varname, dimlen[3], dimlen[2], ni, nj);
        start[1] = k;
        start[2] = j1;
        start[3] = 0;
        count[0] = 1;
        count[1] = 1;
        count[2] = j2 - j1 + 1;
        count[3] = dimlen[3];
        assert(count[2] <= dimlen[2]);
    } else if (ndims == 3) {
        if (nj > 0) {
            if (dimlen[2] != ni || dimlen[1] != nj)
                quit("ncu_readfield_part(): %s: %s: horizontal dimensions (ni = %zu, nj = %zu) do not match grid dimensions (ni = %d, nj = %d)", fname, varname, dimlen[2], dimlen[1], ni, nj);
            if (!hasrecorddim) {
                if (dimlen[0] != nk && !(dimlen[0] == 1 && (k <= 0 || k == nk - 1)))
                    quit("ncu_readfield_part(): %s: %s: vertical dimension (nk = %zu) does not match grid dimension (nk = %d)", fname, varname, dimlen[0], nk);
                start[0] = (dimlen[0] == 1) ? 0 : k;
                start[1] = j1;
                start[2] = 0;
                count[0] = 1;
                count[1] = j2 - j1 + 1;
                count[2] = dimlen[2];
                assert(count[1] <= dimlen[1]);
            } else {
                /*
                 * ignore k in this case
                 */
                if (r < 0) {
                    if (dimlen[0] > 1)
                        quit("ncu_readfield_part(): %s: %s: not supposed to have multiple records", fname, varname);
                    r = 0;
                }
                if (r >= dimlen[0])
                    quit("ncu_readfield_part(): %s: %s: can not read record %d (nr = %zu)", fname, varname, r, dimlen[0]);
                start[0] = r;
                start[1] = j1;
                start[2] = 0;
                count[0] = 1;
                count[1] = j2 - j1 + 1;
                count[2] = dimlen[2];
                assert(count[1] <= dimlen[1]);
            }
        } else {
            if (r < 0) {
                if (dimlen[0] > 1)
                    quit("ncu_readfield_part(): %s: %s: not supposed to have multiple records", fname, varname);
                r = 0;
            }
            if (r >= dimlen[0])
                quit("ncu_readfield_part(): %s: %s: can not read record %d (nr = %zu)", fname, varname, r, dimlen[0]);
            start[0] = r;
            if (dimlen[1] != nk) {
                if (dimlen[1] != 1)
                    quit("ncu_readfield_part(): %s: %s: vertical dimension (nk = %zu) does not match grid dimension (nk = %d)", fname, varname, dimlen[1], nk);
                else
                    k = 0;      /* ignore k */
            }
            if (dimlen[2] != ni)
                quit("ncu_readfield_part(): %s: %s: horizontal dimension (ni = %zu) does not match grid dimension (ni = %d)", fname, varname, dimlen[2], ni);
            start[1] = k;
            start[2] = j1;
            count[0] = 1;
            count[1] = 1;
            count[2] = j2 - j1 + 1;
            assert(count[2] <= dimlen[2]);
        }
    } else if (ndims == 2) {
        if (nj > 0) {
            if (hasrecorddim)
                quit("ncu_readfield_part(): %s: %s: can not read layer from a 1D variable", fname, varname);
            /*
             * ignore k
             */
            start[0] = j1;
            start[1] = 0;
            count[0] = j2 - j1 + 1;
            count[1] = dimlen[1];
            assert(count[0] <= dimlen[0]);
            if (dimlen[1] != ni || dimlen[0] != nj)
                quit("ncu_readfield_part(): %s: %s: horizontal dimensions (ni = %zu, nj = %zu) do not match grid dimensions (ni = %d, nj = %d)", fname, varname, dimlen[1], dimlen[0], ni, nj);
        } else {
            if (!hasrecorddim) {
                if (dimlen[0] != nk && !(dimlen[0] == 1 && (k <= 0 || k == nk - 1)))
                    quit("ncu_readfield_part(): %s: %s: vertical dimension (nk = %zu) does not match grid dimension (nk = %d)", fname, varname, dimlen[0], nk);
                start[0] = (dimlen[0] == 1) ? 0 : k;
                start[1] = j1;
                count[0] = 1;
                count[1] = j2 - j1 + 1;
                assert(count[1] <= dimlen[1]);
            } else {
                /*
                 * ignore k in this case
                 */
                start[0] = r;
                start[1] = j1;
                count[0] = 1;
                count[1] = j2 - j1 + 1;
                assert(count[1] <= dimlen[1]);
            }
            if (dimlen[1] != ni)
                quit("ncu_readfield_part(): %s: %s: horizontal dimension (ni = %zu) does not match grid dimension (ni = %d)", fname, varname, dimlen[1], ni);
        }
    } else if (ndims == 1) {
        if (nj > 0) {
            quit("ncu_readfield_part(): %s: %s: can not read 2D field: # of dimensions = %d", fname, varname, ndims);
        } else {
            if (hasrecorddim)
                quit("ncu_readfield_part(): %s: %s: can not read layer from a 0D variable", fname, varname);
            /*
             * ignore k in this case
             */
            start[0] = j1;
            count[0] = j2 - j1 + 1;
            assert(count[0] <= dimlen[0]);
        }
    } else
        quit("ncu_readfield_part(): %s: %s: can not read 2D field: # of dimensions = %d", fname, varname, ndims);

    ncw_get_vara_float_fixerange(ncid, varid, start, count, v);

    n = 1;
    for (i = 0; i < ndims; ++i)
        n *= count[i];

    ncu_postprocess_read(ncid, varid, start, count, n, v, 0, 0);
    ncw_close(ncid);
}

/** Writes one horizontal field (layer) for a variable to a NetCDF file.
 */
void ncu_writefield(char fname[], char varname[], int k, int ni, int nj, int nk, float* v)
{
    int ncid;
    int varid;
    int ndims;
    size_t dimlen[4];
    size_t start[4], count[4];
    size_t i, n;
    int hasrecorddim;

    ncw_open(fname, NC_WRITE, &ncid);
    ncw_inq_varid(ncid, varname, &varid);
    ncw_inq_vardims(ncid, varid, 4, &ndims, dimlen);
    hasrecorddim = ncw_var_hasunlimdim(ncid, varid);

    if (ndims == 4) {
        if (nj == 0)
            quit("ncu_writefield(): %s: %s: expected positive \"j\" dimension for a 4-dimensional variable\n", fname, varname);
        if (!hasrecorddim && dimlen[0] != 1)
            quit("ncu_writefield(): %s: %s: for a 4-dimensional variable expected the first dimension to be either unlimited or of length 1\n", fname, varname);
        start[0] = (dimlen[0] == 0) ? 0 : dimlen[0] - 1;
        if (dimlen[1] != nk)
            quit("ncu_writefield(): %s: %s: vertical dimension (nk = %zu) does not match grid dimension (nk = %d)", fname, varname, dimlen[1], nk);
        start[1] = k;
        start[2] = 0;
        start[3] = 0;
        count[0] = 1;
        count[1] = 1;
        count[2] = dimlen[2];
        count[3] = dimlen[3];
        if (dimlen[3] != ni || dimlen[2] != nj)
            quit("ncu_writefield(): %s: %s: horizontal dimensions (ni = %zu, nj = %zu) do not match grid dimensions (ni = %d, nj = %d)", fname, varname, dimlen[3], dimlen[2], ni, nj);
    } else if (ndims == 3) {
        if (nj > 0) {
            if (!hasrecorddim) {
                if (dimlen[0] != nk && !(dimlen[0] == 1 && (k == 0 || k == nk - 1)))
                    quit("ncu_writefield(): %s: %s: vertical dimension (nk = %zu) does not match grid dimension (nk = %d)", fname, varname, dimlen[0], nk);
                start[0] = (dimlen[0] == 1) ? 0 : k;
                start[1] = 0;
                start[2] = 0;
                count[0] = 1;
                count[1] = dimlen[1];
                count[2] = dimlen[2];
            } else {
                /*
                 * 2D variable, ignore k
                 */
                start[0] = (dimlen[0] == 0) ? 0 : dimlen[0] - 1;
                start[1] = 0;
                start[2] = 0;
                count[0] = 1;
                count[1] = dimlen[1];
                count[2] = dimlen[2];
            }
            if (dimlen[2] != ni || dimlen[1] != nj)
                quit("ncu_writefield(): %s: %s: horizontal dimensions (ni = %zu, nj = %zu) do not match grid dimensions (ni = %d, nj = %d)", fname, varname, dimlen[2], dimlen[1], ni, nj);
        } else {
            if (!hasrecorddim && dimlen[0] != 1)
                quit("ncu_writefield(): %s: %s: for a 3-dimensional variable on unstructured grid expected the first dimension to be either unlimited or of length 1\n", fname, varname);
            start[0] = (dimlen[0] == 0) ? 0 : dimlen[0] - 1;
            if (dimlen[1] != nk) {
                if (dimlen[1] != 1)
                    quit("ncu_writefield(): %s: %s: vertical dimension (nk = %zu) does not match grid dimension (nk = %d)", fname, varname, dimlen[1], nk);
                else
                    k = 0;      /* ignore k */
            }
            if (dimlen[2] != ni)
                quit("ncu_writefield(): %s: %s: horizontal dimension (ni = %zu) does not match grid dimension (ni = %d)", fname, varname, dimlen[2], ni);
            start[1] = k;
            start[2] = 0;
            count[0] = 1;
            count[1] = 1;
            count[2] = dimlen[2];
        }
    } else if (ndims == 2) {
        if (nj > 0) {
            if (hasrecorddim)
                quit("ncu_writefield(): %s: %s: can not write a layer to a 1D variable", fname, varname);
            /*
             * ignore k
             */
            start[0] = 0;
            start[1] = 0;
            count[0] = dimlen[0];
            count[1] = dimlen[1];
            if (dimlen[1] != ni || dimlen[0] != nj)
                quit("ncu_writefield(): %s: %s: horizontal dimensions (ni = %zu, nj = %zu) do not match grid dimensions (ni = %d, nj = %d)", fname, varname, dimlen[1], dimlen[0], ni, nj);
        } else {
            if (!hasrecorddim) {
                if (dimlen[0] != nk && !(dimlen[0] == 1 && (k == 0 || k == nk - 1)))
                    quit("ncu_writefield(): %s: %s: vertical dimension (nk = %zu) does not match grid dimension (nk = %d)", fname, varname, dimlen[0], nk);
                start[0] = (dimlen[0] == 1) ? 0 : k;
                start[1] = 0;
                count[0] = 1;
                count[1] = dimlen[1];
            } else {
                /*
                 * ignore k in this case
                 */
                start[0] = (dimlen[0] == 0) ? 0 : dimlen[0] - 1;
                start[1] = 0;
                count[0] = 1;
                count[1] = dimlen[1];
            }
            if (dimlen[1] != ni)
                quit("ncu_writefield(): %s: %s: horizontal dimension (ni = %zu) does not match grid dimension (ni = %d)", fname, varname, dimlen[1], ni);
        }
    } else if (ndims == 1) {
        if (nj > 0) {
            quit("ncu_writefield(): %s: %s: can not write 2D field: # of dimensions = %d", fname, varname, ndims);
        } else {
            if (hasrecorddim)
                quit("ncu_writefield(): %s: %s: can not write layer to a 0D variable", fname, varname);
            /*
             * ignore k in this case
             */
            start[0] = 0;
            count[0] = dimlen[0];
        }
    } else
        quit("ncu_writefield(): %s: %s: can not write 2D field: # of dimensions = %d", fname, varname, ndims);

    n = 1;
    for (i = 0; i < ndims; ++i)
        n *= count[i];

    ncu_preprocess_write(ncid, varid, n, v);
    ncu_put_vara(ncid, varid, start, count, n, v);
    ncw_close(ncid);
}

/** Writes one row of a horizontal field (layer) for a variable to a NetCDF
 *  file.
 */
void ncu_writerow(char fname[], char varname[], int k, int j, int ni, int nj, int nk, float* v)
{
    int ncid;
    int varid;
    int ndims;
    size_t dimlen[4];
    size_t start[4], count[4];
    size_t n;
    int hasrecorddim;

    ncw_open(fname, NC_WRITE, &ncid);
    ncw_inq_varid(ncid, varname, &varid);
    ncw_inq_vardims(ncid, varid, 4, &ndims, dimlen);
    hasrecorddim = ncw_var_hasunlimdim(ncid, varid);

    if (ndims == 4) {
        if (nj == 0)
            quit("ncu_writerow(): \"%s\": %s: expected positive \"j\" dimension for a 4-dimensional variable\n", fname, varname);
        if (!hasrecorddim && dimlen[0] != 1)
            quit("%s: %s: for a 4-dimensional variable expected the first dimension to be either unlimited or of length 1\n", fname, varname);
        start[0] = (dimlen[0] == 0) ? 0 : dimlen[0] - 1;
        if (dimlen[1] != nk)
            quit("ncu_writerow(): \"%s\": vertical dimension of variable \"%s\" (nk = %zu) does not match grid dimension (nk = %d)", fname, varname, dimlen[1], nk);
        start[1] = k;
        start[2] = j;
        start[3] = 0;
        count[0] = 1;
        count[1] = 1;
        count[2] = 1;
        count[3] = dimlen[3];
        if (dimlen[3] != ni || dimlen[2] != nj)
            quit("ncu_writerow(): \"%s\": horizontal dimensions of variable \"%s\" (ni = %zu, nj = %zu) do not match grid dimensions (ni = %d, nj = %d)", fname, varname, dimlen[3], dimlen[2], ni, nj);
    } else if (ndims == 3) {
        if (nj > 0) {
            if (!hasrecorddim) {
                if (dimlen[0] != nk && !(dimlen[0] == 1 && (k == 0 || k == nk - 1)))
                    quit("ncu_writerow(): \"%s\": vertical dimension of variable \"%s\" (nk = %zu) does not match grid dimension (nk = %d)", fname, varname, dimlen[0], nk);
                start[0] = (dimlen[0] == 1) ? 0 : k;
                start[1] = j;
                start[2] = 0;
                count[0] = 1;
                count[1] = 1;
                count[2] = dimlen[2];
            } else {
                /*
                 * 2D variable, ignore k
                 */
                start[0] = (dimlen[0] == 0) ? 0 : dimlen[0] - 1;
                start[1] = j;
                start[2] = 0;
                count[0] = 1;
                count[1] = 1;
                count[2] = dimlen[2];
            }
            if (dimlen[2] != ni || dimlen[1] != nj)
                quit("ncu_writerow(): \"%s\": horizontal dimensions of variable \"%s\" (ni = %zu, nj = %zu) do not match grid dimensions (ni = %d, nj = %d)", fname, varname, dimlen[2], dimlen[1], ni, nj);
        } else {
            if (!hasrecorddim && dimlen[0] != 1)
                quit("ncu_writerow(): \"%s\": %s: for a 3-dimensional variable on unstructured horizontal grid expected the first dimension to be either unlimited or of length 1\n", fname, varname);
            start[0] = (dimlen[0] == 0) ? 0 : dimlen[0] - 1;
            if (dimlen[1] != nk) {
                if (dimlen[1] != 1)
                    quit("ncu_writerow(): \"%s\": vertical dimension of variable \"%s\" (nk = %zu) does not match grid dimension (nk = %d)", fname, varname, dimlen[1], nk);
                else
                    k = 0;      /* ignore k */
            }
            if (dimlen[2] != ni)
                quit("ncu_writerow(): \"%s\": horizontal dimension of variable \"%s\" (ni = %zu) does not match grid dimension (ni = %d)", fname, varname, dimlen[2], ni);
            start[1] = k;
            start[2] = j;
            count[0] = 1;
            count[1] = 1;
            count[2] = 1;
        }
    } else if (ndims == 2) {
        if (nj > 0) {
            if (hasrecorddim)
                quit("%s: can not write a row %d to a 1D variable \"%s\"", fname, j, varname);
            /*
             * ignore k
             */
            start[0] = j;
            start[1] = 0;
            count[0] = 1;
            count[1] = dimlen[1];
            if (dimlen[1] != ni || dimlen[0] != nj)
                quit("ncu_writerow(): \"%s\": horizontal dimensions of variable \"%s\" (ni = %zu, nj = %zu) do not match grid dimensions (ni = %d, nj = %d)", fname, varname, dimlen[1], dimlen[0], ni, nj);
        } else {
            if (!hasrecorddim) {
                if (dimlen[0] != nk && !(dimlen[0] == 1 && (k == 0 || k == nk - 1)))
                    quit("ncu_writerow(): \"%s\": vertical dimension of variable \"%s\" (nk = %zu) does not match grid dimension (nk = %d)", fname, varname, dimlen[0], nk);
                start[0] = (dimlen[0] == 1) ? 0 : k;
                start[1] = j;
                count[0] = 1;
                count[1] = 1;
            } else {
                /*
                 * ignore k in this case
                 */
                start[0] = (dimlen[0] == 0) ? 0 : dimlen[0] - 1;
                start[1] = j;
                count[0] = 1;
                count[1] = 1;
            }
            if (dimlen[1] != ni)
                quit("ncu_writerow(): \"%s\": horizontal dimension of variable \"%s\" (ni = %zu) does not match grid dimension (ni = %d)", fname, varname, dimlen[1], ni);
        }
    } else if (ndims == 1) {
        if (nj > 0) {
            quit("ncu_writerow(): %s: can not write 2D field for \"%s\": # of dimensions = %d", fname, varname, ndims);
        } else {
            if (hasrecorddim)
                quit("ncu_writerow(): %s: can not write layer to a 0D variable \"%s\"", fname, varname);
            /*
             * ignore k in this case
             */
            start[0] = j;
            count[0] = 1;
        }
    } else
        quit("%s: can not write a row to 2D field for \"%s\": # of dimensions = %d", fname, varname, ndims);

    n = count[ndims - 1];

    ncu_preprocess_write(ncid, varid, n, v);
    ncu_put_vara(ncid, varid, start, count, n, v);
    ncw_close(ncid);
}

/**
 */
void ncu_read3dfield(char* fname, char* varname, int r, int ni, int nj, int nk, float* v)
{
    int ncid;
    int varid;
    int ndims;
    size_t dimlen[4];
    size_t start[4], count[4];
    size_t i, n;
    int hasrecorddim;

    ncw_open(fname, NC_NOWRITE, &ncid);
    ncw_inq_varid(ncid, varname, &varid);
    ncw_inq_vardims(ncid, varid, 4, &ndims, dimlen);
    hasrecorddim = ncw_var_hasunlimdim(ncid, varid);
    if (hasrecorddim && dimlen[0] == 0)
        quit("ncu_read3dfield(): %s: %s: empty record dimension", fname, varname);

    if (ndims == 4) {
        if (r < 0) {
            if (dimlen[0] > 1)
                quit("ncu_read3dfield(): \"%s\": %s: not supposed to have multiple records", fname, varname);
            r = 0;
        }
        if (r >= dimlen[0])
            quit("ncu_read3dfield(): \"%s\": %s: can not read record %d (nr = %zu)", fname, varname, r, dimlen[0]);
        if (nj == 0)            /* unstructured grid */
            quit("ncu_read3dfield(): \"%s\": %s: too many dimensions for an unstructured grid", fname, varname);
        start[0] = r;
        start[1] = 0;
        start[2] = 0;
        start[3] = 0;
        count[0] = 1;
        count[1] = dimlen[1];
        count[2] = dimlen[2];
        count[3] = dimlen[3];
        if (dimlen[3] != ni || dimlen[2] != nj || dimlen[1] != nk)
            quit("ncu_read3dfield(): %s: %s: horizontal dimensions (ni = %zu, nj = %zu, nk = %zu) do not match grid dimensions (ni = %d, nj = %d, nk = %d)", fname, varname, dimlen[3], dimlen[2], dimlen[1], ni, nj, nk);
    } else if (ndims == 3) {
        if (nj > 0) {
            if (hasrecorddim)
                quit("ncu_read3dfield(): %s: %s: can not read 3D field because the variable is only 2D", fname, varname);
            if (r > 0)
                quit("ncu_read3dfield(): %s: %s: can not read 3D record %d because the variable has only 3 dimensions", fname, varname, r);
            if (dimlen[2] != ni || dimlen[1] != nj || dimlen[0] != nk)
                quit("ncu_read3dfield(): %s: %s: horizontal dimensions (ni = %zu, nj = %zu, nk = %zu) do not match grid dimensions (ni = %d, nj = %d, nk = %d)", fname, varname, dimlen[2], dimlen[1], dimlen[0], ni, nj, nk);
            start[0] = 0;
            start[1] = 0;
            start[2] = 0;
            count[0] = dimlen[0];
            count[1] = dimlen[1];
            count[2] = dimlen[2];
        } else {
            if (r < 0) {
                if (dimlen[0] > 1)
                    quit("ncu_read3dfield(): \"%s\": %s: (unstructured grid): not supposed to have multiple records", fname, varname);
                r = 0;
            }
            if (r >= dimlen[0])
                quit("ncu_read3dfield(): \"%s\": %s: (unstructured grid): can not read record %d (nr = %zu)", fname, varname, r, dimlen[0]);
            start[0] = r;
            start[1] = 0;
            start[2] = 0;
            count[0] = 1;
            count[1] = dimlen[1];
            count[2] = dimlen[2];
        }
    } else if (ndims == 2) {
        if (hasrecorddim)
            quit("ncu_read3dfield(): %s: %s: can not read 3D field: too few dimensions", fname, varname);
        if (r > 0)
            quit("ncu_read3dfield(): %s: %s: can not read 3D record %d: too few dimensions", fname, varname, r);
        if (nj > 0 && nk > 0)
            quit("ncu_read3dfield(): %s: %s: can not read 3D field: too few dimensions", fname, varname);
        start[0] = 0;
        start[1] = 0;
        count[0] = dimlen[0];
        count[1] = dimlen[1];
    } else if (ndims == 1) {
        if (hasrecorddim)
            quit("ncu_read3dfield(): %s: %s: can not read 3D field: too few dimensions", fname, varname);
        if (r > 0)
            quit("ncu_read3dfield(): %s: %s: can not read 3D record %d: too few dimensions", fname, varname, r);
        if (nj > 0)
            quit("ncu_read3dfield(): %s: %s: can not read 3D field: # of dimensions = %d", fname, varname, ndims);
        start[0] = 0;
        count[0] = dimlen[0];
    }

    ncw_get_vara_float(ncid, varid, start, count, v);

    n = 1;
    for (i = 0; i < ndims; ++i)
        n *= count[i];

    ncu_postprocess_read(ncid, varid, start, count, n, v, 0, 0);
    ncw_close(ncid);
}

/** Tries to determine the number of physical dimensions of a variable.
 */
int ncu_getnD(char fname[], char varname[])
{
    int ncid;
    int varid;
    int ndims;
    size_t dimlen[4];
    int i;

    ncw_open(fname, NC_NOWRITE, &ncid);
    ncw_inq_varid(ncid, varname, &varid);
    ncw_inq_vardims(ncid, varid, 4, &ndims, dimlen);
    if (ndims > 4)
        quit("%s: %s: do not know how to treat a %d-dimensional variable", fname, varname, ndims);

    for (i = ncw_var_hasunlimdim(ncid, varid) ? 1 : 0; i < ndims; ++i)
        if (dimlen[i] > 1)
            break;
    ncw_close(ncid);

    return ndims - i;
}

/**
 */
void ncu_readvarfloat(int ncid, int varid, size_t n, float v[])
{
    ncw_check_varsize(ncid, varid, n);
    ncw_get_var_float_fixerange(ncid, varid, v);
    ncu_postprocess_read(ncid, varid, NULL, NULL, n, v, 0, 1);
}

/**
 */
void ncu_readvardouble(int ncid, int varid, size_t n, double v[])
{
    ncw_check_varsize(ncid, varid, n);
    ncw_get_var_double(ncid, varid, v);
    ncu_postprocess_read(ncid, varid, NULL, NULL, n, v, 1, 1);
}
