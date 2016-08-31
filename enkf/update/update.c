/******************************************************************************
 *
 * File:        update.c        
 *
 * Created:     11/2013
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description: The code in this file structurally belongs to dasystem.c,
 *              and is put in a separate file just to break dasystem.c in
 *              smaller parts.
 *
 * Revisions:
 *
 *****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <math.h>
#include "definitions.h"
#include "nan.h"
#include "utils.h"
#include "distribute.h"
#include "grid.h"
#include "lapack.h"
#include "dasystem.h"
#include "pointlog.h"

#define EPSF 1.0e-6f

/**
 */
static void getfieldfname(char* dir, char* prefix, char* varname, int level, char* fname)
{
    snprintf(fname, MAXSTRLEN, "%s/%s_%s-%03d.nc", dir, prefix, varname, level);
}

/**
 */
static void das_writeinflation(dasystem* das, field* f, int j, float* v)
{
    assert(das->mode == MODE_ENKF);

    if (das->updatespec & UPDATE_DIRECTWRITE)
        writerow(FNAME_INFLATION, f->varname, f->level, j, v);
    else {
        char fname[MAXSTRLEN];
        int ncid, varid;

        getfieldfname(das->ensdir, "inflation", f->varname, f->level, fname);

        if (j == 0) {
            int ni, nj;
            int dimids[2];

            model_getvardims(das->m, f->varid, &ni, &nj, NULL);

            ncw_create(fname, NC_CLOBBER | NETCDF_FORMAT, &ncid);
            ncw_def_dim(ncid, "nj", nj, &dimids[0]);
            ncw_def_dim(ncid, "ni", ni, &dimids[1]);
            ncw_def_var(ncid, f->varname, NC_FLOAT, 2, dimids, NULL);
            ncw_enddef(ncid);
        } else {
            ncw_open(fname, NC_WRITE, &ncid);
            ncw_inq_varid(ncid, f->varname, &varid);
        }
        writerow(fname, f->varname, 0, j, v);
        ncw_close(ncid);
    }
}

/** Updates `nfield' fields read into `fieldbuffer' with `X5'. Applies
 * variable-dependent inflation to each field.
 * 
 */
static void das_updatefields(dasystem* das, int nfields, void** fieldbuffer, field fields[])
{
    model* m = das->m;
    int nmem = das->nmem;

    void* grid = model_getvargrid(m, fields[0].varid);
    int** nlevels = grid_getnumlevels(grid);
    int periodic_i = grid_isperiodic_x(grid);
    int writeinflation = das->updatespec & UPDATE_DOINFLATION;

    /*
     * X5
     */
    char fname_X5[MAXSTRLEN];
    int ncid;
    int varid;
    int dimids[3];
    size_t dimlens[3];
    size_t start[3], count[3];
    float** X5jj = NULL;
    float** X5jj1 = NULL;
    float** X5jj2 = NULL;
    float** X5j = NULL;

    float* v_f = NULL;          /* v_f = E_f(i, :) */
    float* v_a = NULL;          /* v_a = E_a(i, :) */
    float* infl = NULL;

    int mni, mnj;
    int ni, nj;
    int i, j;
    int jj, stepj, ii, stepi;
    int e, fid;

    assert(das->mode == MODE_ENKF);

    grid_getdims(grid, &mni, &mnj, NULL);

    das_getfname_X5(das, grid, fname_X5);

    ncw_open(fname_X5, NC_NOWRITE, &ncid);
    ncw_inq_varid(ncid, "X5", &varid);
    ncw_inq_vardimid(ncid, varid, dimids);
    for (i = 0; i < 3; ++i)
        ncw_inq_dimlen(ncid, dimids[i], &dimlens[i]);
    nj = dimlens[0];
    ni = dimlens[1];
    assert((int) dimlens[2] == nmem * nmem);

    start[0] = 0;
    start[1] = 0;
    start[2] = 0;
    count[0] = 1;
    count[1] = ni;
    count[2] = nmem * nmem;

    X5j = alloc2d(mni, nmem * nmem, sizeof(float));
    if (das->stride > 1) {
        X5jj = alloc2d(ni, nmem * nmem, sizeof(float));
        X5jj1 = alloc2d(ni, nmem * nmem, sizeof(float));
        X5jj2 = alloc2d(ni, nmem * nmem, sizeof(float));

        ncw_get_vara_float(ncid, varid, start, count, X5jj2[0]);
    }

    v_f = malloc(nmem * sizeof(float));
    v_a = malloc(nmem * sizeof(float));
    if (writeinflation)
        infl = malloc(mni * sizeof(float));

    /*
     * jj, ii are the indices of the subsampled grid; i, j are the indices
     * of the model grid 
     */
    for (jj = 0, j = 0; jj < nj; ++jj) {
        for (stepj = 0; stepj < das->stride && j < mnj; ++stepj, ++j) {
            if (das->stride == 1) {
                /*
                 * no interpolation necessary; simply read the ETMs for the
                 * j-th row from disk 
                 */
                start[0] = j;
                ncw_get_vara_float(ncid, varid, start, count, X5j[0]);
            } else {
                /*
                 * the following code interpolates the ETM back to the
                 * original grid, first by j, and then by i 
                 */
                if (stepj == 0) {
                    memcpy(X5jj1[0], X5jj2[0], ni * nmem * nmem * sizeof(float));
                    memcpy(X5jj[0], X5jj2[0], ni * nmem * nmem * sizeof(float));
                    if (jj < nj - 1) {
                        start[0] = (jj + 1) % nj;
                        ncw_get_vara_float(ncid, varid, start, count, X5jj2[0]);
                    }
                } else {
                    float weight2 = (float) stepj / (float) das->stride;
                    float weight1 = 1.0f - weight2;

                    for (ii = 0; ii < ni; ++ii) {
                        float* X5jjii = X5jj[ii];
                        float* X5jj1ii = X5jj1[ii];
                        float* X5jj2ii = X5jj2[ii];

                        for (e = 0; e < nmem * nmem; ++e)
                            X5jjii[e] = X5jj1ii[e] * weight1 + X5jj2ii[e] * weight2;
                    }
                }

                for (ii = 0, i = 0; ii < ni; ++ii) {
                    for (stepi = 0; stepi < das->stride && i < mni; ++stepi, ++i) {
                        if (stepi == 0) {
                            memcpy(X5j[i], X5jj[ii], nmem * nmem * sizeof(float));
                        } else {
                            float weight2 = (float) stepi / (float) das->stride;
                            float weight1 = 1.0f - weight2;
                            float* X5jjii1 = X5jj[ii];
                            float* X5ji = X5j[i];
                            float* X5jjii2;

                            if (ii < ni - 1)
                                X5jjii2 = X5jj[ii + 1];
                            else
                                X5jjii2 = X5jj[(periodic_i) ? (ii + 1) % ni : ii];
                            for (e = 0; e < nmem * nmem; ++e)
                                X5ji[e] = X5jjii1[e] * weight1 + X5jjii2[e] * weight2;
                        }
                    }
                }
            }                   /* stride != 1 */

            /*
             * (at this stage X5j should contain the array of X5 matrices
             * for the j-th row of the grid) 
             */

            /*
             * update the j-th row of the fields 
             */
            for (fid = 0; fid < nfields; ++fid) {
                field* f = &fields[fid];
                float*** vvv = (float***) fieldbuffer[fid];
                char do_T = 'T';
                float alpha = 1.0f;
                int inc = 1;
                float beta = 0.0f;
                float inflation0 = NaN;
                double inf_ratio = NaN;

                model_getvarinflation(m, f->varid, &inflation0, &inf_ratio);
                if (writeinflation)
                    memset(infl, 0, mni * sizeof(float));

                for (i = 0; i < mni; ++i) {
                    float inflation = inflation0;
                    double v1_a = 0.0;

                    /*
                     * For now we assume that layers are counted down from the
                     * surface. (This is not so in ROMS, but there the number
                     * of layers is always the same.) This will be easy to
                     * modify as soon as we encounter a Z model with layers
                     * counted up from the bottom.
                     */
                    if (f->level >= nlevels[j][i])
                        continue;
                    /*
                     * assume that if |value| > MAXOBSVAL, then it is filled
                     * with the missing value 
                     */
                    if (fabsf(vvv[0][j][i]) > (float) MAXOBSVAL)
                        continue;
                    /*
                     * (it would be straightforward to compare with the
                     * actual missing value, provided that it is NC_FLOAT;
                     * otherwise it may be a bit tiresome to handle all
                     * variations) 
                     */

                    /*
                     * write the actual (interpolated) transform  matrix to the
                     * pointlog for this (i,j) pair (if it exists)
                     */
                    {
                        unsigned short key[4] = { i, j, grid_getid(grid), 0 };
                        pointlog* plog = NULL;

                        if (f->id == 0 && das->nplogs > 0 && (plog = ht_find(das->ht_plogs, key)) != NULL)
                            plog_writeactualtransform(das, plog->id, X5j[i]);
                    }

                    for (e = 0; e < nmem; ++e)
                        v_f[e] = vvv[e][j][i];

                    /*
                     * E_a(i, :) = E_f(i, :) * X5 
                     */
                    sgemv_(&do_T, &nmem, &nmem, &alpha, X5j[i], &nmem, v_f, &inc, &beta, v_a, &inc);

                    for (e = 0; e < nmem; ++e)
                        v1_a += v_a[e];
                    v1_a /= (double) nmem;

                    if (!isnan(inf_ratio)) {
                        double v1_f = 0.0;
                        double v2_f = 0.0;
                        double v2_a = 0.0;
                        double var_a, var_f;

                        for (e = 0; e < nmem; ++e) {
                            double ve = (double) v_f[e];

                            v1_f += ve;
                            v2_f += ve * ve;
                        }
                        v1_f /= (double) nmem;
                        var_f = v2_f / (double) nmem - v1_f * v1_f;

                        for (e = 0; e < nmem; ++e) {
                            double ve = (double) v_a[e];

                            v2_a += ve * ve;
                        }
                        var_a = v2_a / (double) nmem - v1_a * v1_a;

                        if (var_a > 0) {
                            /*
                             * (Normal case.) Limit inflation by inf_ratio of
                             * the magnitude of spread reduction.
                             */
                            inflation = (float) (sqrt(var_f / var_a) * inf_ratio + 1.0 - inf_ratio);
                            if (inflation >= inflation0)
                                inflation = inflation0;
                        }
                    }

                    /*
                     * (Do not inflate if inflation is about 1 or less than 1.)
                     */
                    if (inflation - 1.0f > EPSF)
                        for (e = 0; e < nmem; ++e)
                            v_a[e] = (v_a[e] - (float) v1_a) * inflation + (float) v1_a;

                    if (!(das->updatespec & UPDATE_OUTPUTINC))
                        for (e = 0; e < nmem; ++e)
                            vvv[e][j][i] = v_a[e];
                    else
                        for (e = 0; e < nmem; ++e)
                            vvv[e][j][i] = v_a[e] - v_f[e];

                    if (writeinflation)
                        infl[i] = inflation;
                }               /* for i */
                if (writeinflation)
                    das_writeinflation(das, f, j, infl);
            }                   /* for fid */
        }                       /* for stepj */
    }                           /* for jj */

    ncw_close(ncid);
    free(X5j);
    if (das->stride > 1) {
        free(X5jj);
        free(X5jj1);
        free(X5jj2);
    }
    if (writeinflation)
        free(infl);
    free(v_a);
    free(v_f);

    /*
     * "randomise" ("propagate") fields if required
     */
    for (fid = 0; fid < nfields; ++fid) {
        field* f = &fields[fid];
        float*** vvv = (float***) fieldbuffer[fid];

        if (!isnan(model_getvardeflation(m, f->varid))) {
            int e;

            for (e = 0; e < das->nmem; ++e)
                model_randomisefield(m, f->varid, vvv[e]);
        }
    }
}

/** Updates `nfield' fields read into `fieldbuffer' with `w'. 
 */
static void das_updatebg(dasystem* das, int nfields, void** fieldbuffer, field fields[])
{
    model* m = das->m;
    int nmem = das->nmem;

    void* grid = model_getvargrid(m, fields[0].varid);
    int** nlevels = grid_getnumlevels(grid);
    int periodic_i = grid_isperiodic_x(grid);

    char fname_w[MAXSTRLEN];
    int ncid;
    int varid;
    int dimids[3];
    size_t dimlens[3];
    size_t start[3], count[3];
    float** wjj = NULL;
    float** wjj1 = NULL;
    float** wjj2 = NULL;
    float** wj;

    int mni, mnj;
    int i, j, ni, nj;
    int jj, stepj, ii, stepi;
    int e, fid;

    /*
     * the following code for interpolation of w essentially coincides with
     * that in das_updatefields()
     */

    assert(das->mode == MODE_ENOI);

    grid_getdims(grid, &mni, &mnj, NULL);

    das_getfname_w(das, grid, fname_w);

    ncw_open(fname_w, NC_NOWRITE, &ncid);
    ncw_inq_varid(ncid, "w", &varid);
    ncw_inq_vardimid(ncid, varid, dimids);
    for (i = 0; i < 3; ++i)
        ncw_inq_dimlen(ncid, dimids[i], &dimlens[i]);
    ni = dimlens[1];
    nj = dimlens[0];
    assert((int) dimlens[2] == nmem);

    start[0] = 0;
    start[1] = 0;
    start[2] = 0;
    count[0] = 1;
    count[1] = ni;
    count[2] = nmem;

    wj = alloc2d(mni, nmem, sizeof(float));
    if (das->stride > 1) {
        wjj = alloc2d(ni, nmem, sizeof(float));
        wjj1 = alloc2d(ni, nmem, sizeof(float));
        wjj2 = alloc2d(ni, nmem, sizeof(float));

        ncw_get_vara_float(ncid, varid, start, count, wjj2[0]);
    }
    /*
     * jj, ii are the indices of the subsampled grid; i, j are the indices
     * of the model grid 
     */
    for (jj = 0, j = 0; jj < nj; ++jj) {
        for (stepj = 0; stepj < das->stride && j < mnj; ++stepj, ++j) {
            if (das->stride == 1) {
                /*
                 * no interpolation necessary; simply read the ETMs for the
                 * j-th row from disk 
                 */
                start[0] = j;
                ncw_get_vara_float(ncid, varid, start, count, wj[0]);
            } else {
                /*
                 * the following code interpolates the ETM back to the
                 * original grid, first by j, and then by i 
                 */
                if (stepj == 0) {
                    memcpy(wjj1[0], wjj2[0], ni * nmem * sizeof(float));
                    memcpy(wjj[0], wjj2[0], ni * nmem * sizeof(float));
                    if (jj < nj - 1) {
                        start[0] = (jj + 1) % nj;
                        ncw_get_vara_float(ncid, varid, start, count, wjj2[0]);
                    }
                } else {
                    float weight2 = (float) stepj / das->stride;
                    float weight1 = (float) 1.0f - weight2;

                    for (ii = 0; ii < ni; ++ii) {
                        float* wjjii = wjj[ii];
                        float* wjj1ii = wjj1[ii];
                        float* wjj2ii = wjj2[ii];

                        for (e = 0; e < nmem; ++e)
                            wjjii[e] = wjj1ii[e] * weight1 + wjj2ii[e] * weight2;
                    }
                }

                for (ii = 0, i = 0; ii < ni; ++ii) {
                    for (stepi = 0; stepi < das->stride && i < mni; ++stepi, ++i) {
                        if (stepi == 0)
                            memcpy(wj[i], wjj[ii], nmem * sizeof(float));
                        else {
                            float weight2 = (float) stepi / das->stride;
                            float weight1 = (float) 1.0f - weight2;
                            float* wjjii1 = wjj[ii];
                            float* wji = wj[i];
                            float* wjjii2;

                            if (ii < ni - 1)
                                wjjii2 = wjj[ii + 1];
                            else
                                wjjii2 = wjj[(periodic_i) ? (ii + 1) % ni : ii];

                            for (e = 0; e < nmem; ++e)
                                wji[e] = wjjii1[e] * weight1 + wjjii2[e] * weight2;
                        }
                    }
                }
            }                   /* stride != 1 */

            /*
             * (at this stage wj should contain the array of w vectors for
             * the j-th row of the grid) 
             */

            for (fid = 0; fid < nfields; ++fid) {
                field* f = &fields[fid];
                float*** vvv = (float***) fieldbuffer[fid];

                for (i = 0; i < mni; ++i) {
                    float xmean = 0.0f;

                    /*
                     * For now we assume that layers are counted down from the
                     * surface. (This is not so in ROMS, but there the number
                     * of layers is always the same.) This will be easy to
                     * modify as soon as we encounter a Z model with layers
                     * counted up from the bottom.
                     */
                    if (f->level >= nlevels[j][i])
                        continue;
                    /*
                     * assume that if |value| > MAXOBSVAL, then it is filled
                     * with the missing value 
                     */
                    if (fabsf(vvv[nmem][j][i]) > (float) MAXOBSVAL)
                        continue;
                    if (fabsf(vvv[0][j][i]) > (float) MAXOBSVAL)
                        continue;

                    /*
                     * for fid = 0 write the actual (interpolated) weight
                     * vector to the pointlog for this (i,j) pair (if it exists)
                     */
                    if (f->id == 0) {
                        unsigned short key[4] = { i, j, grid_getid(grid), 0 };
                        int id;

                        if (f->id == 0 && das->nplogs > 0 && (id = ht_findid(das->ht_plogs, key)) >= 0)
                            plog_writeactualtransform(das, id, wj[i]);
                    }

                    for (e = 0; e < nmem; ++e)
                        xmean += vvv[e][j][i];
                    xmean /= (float) nmem;

                    /*
                     * (the case das->updatespec & UPDATE_OUTPUTINC > 0 is
                     * handled by setting vvv[nmem][][] to zero in das_update())
                     */
                    for (e = 0; e < nmem; ++e)
                        vvv[nmem][j][i] += (vvv[e][j][i] - xmean) * wj[i][e];
                }
            }
        }                       /* for stepj */
    }                           /* for jj */

    ncw_close(ncid);
    free(wj);
    if (das->stride > 1) {
        free(wjj);
        free(wjj1);
        free(wjj2);
    }
    das->s_mode = S_MODE_HE_a;
}

/**
 */
static void das_writefields_direct(dasystem* das, int nfields, void** fieldbuffer, field fields[])
{
    int i, e;

    assert(das->mode == MODE_ENKF);

    if (!(das->updatespec & UPDATE_SEPARATEOUTPUT)) {
        for (i = 0; i < nfields; ++i) {
            field* f = &fields[i];
            char varname[NC_MAX_NAME];

            strncpy(varname, f->varname, NC_MAX_NAME);
            if (!(das->updatespec & UPDATE_OUTPUTINC))
                strncat(varname, "_an", NC_MAX_NAME);
            else
                strncat(varname, "_inc", NC_MAX_NAME);

            for (e = 0; e < das->nmem; ++e) {
                char fname[MAXSTRLEN];

                model_getmemberfname(das->m, das->ensdir, f->varname, e + 1, fname);
                model_writefield(das->m, fname, INT_MAX, varname, f->level, ((float***) fieldbuffer[i])[e][0]);
            }
        }
    } else {
        for (i = 0; i < nfields; ++i) {
            field* f = &fields[i];

            for (e = 0; e < das->nmem; ++e) {
                char fname[MAXSTRLEN];

                model_getmemberfname(das->m, das->ensdir, f->varname, e + 1, fname);
                if (!(das->updatespec & UPDATE_OUTPUTINC))
                    strncat(fname, ".analysis", MAXSTRLEN);
                else
                    strncat(fname, ".increment", MAXSTRLEN);
                model_writefield(das->m, fname, INT_MAX, f->varname, f->level, ((float***) fieldbuffer[i])[e][0]);
            }
        }
    }
}

/**
 */
static void das_writefields_toassemble(dasystem* das, int nfields, void** fieldbuffer, field fields[])
{
    char fname[MAXSTRLEN];
    int ni, nj;
    int i;

    model_getvardims(das->m, fields[0].varid, &ni, &nj, NULL);

    for (i = 0; i < nfields; ++i) {
        field* f = &fields[i];

        getfieldfname(das->ensdir, "ens", f->varname, f->level, fname);

        if (!file_exists(fname)) {
            int ncid;
            int dimids[3];
            int vid;

            ncw_create(fname, NC_CLOBBER | NETCDF_FORMAT, &ncid);
            ncw_def_dim(ncid, "m", das->nmem, &dimids[0]);
            ncw_def_dim(ncid, "nj", nj, &dimids[1]);
            ncw_def_dim(ncid, "ni", ni, &dimids[2]);
            ncw_def_var(ncid, f->varname, NC_FLOAT, 3, dimids, &vid);
            ncw_enddef(ncid);
            ncw_put_var_float(ncid, vid, ((float***) fieldbuffer[i])[0][0]);
            ncw_close(ncid);
        } else {
            int ncid;
            int vid;

            ncw_open(fname, NC_WRITE, &ncid);
            ncw_inq_varid(ncid, f->varname, &vid);
            ncw_put_var_float(ncid, vid, ((float***) fieldbuffer[i])[0][0]);
            ncw_close(ncid);
        }
    }
}

/** Writes `nfields' ensemble fields from `fieldbuffer' to disk.
 */
static void das_writefields(dasystem* das, int nfields, void** fieldbuffer, field fields[])
{
    if (das->updatespec & UPDATE_DIRECTWRITE)
        das_writefields_direct(das, nfields, fieldbuffer, fields);
    else
        das_writefields_toassemble(das, nfields, fieldbuffer, fields);
}

static void das_writebg_direct(dasystem* das, int nfields, void** fieldbuffer, field fields[])
{
    model* m = das->m;
    int ni, nj;
    int i;

    assert(das->mode == MODE_ENOI);

    model_getvardims(m, fields[0].varid, &ni, &nj, NULL);

    if (!(das->updatespec & UPDATE_SEPARATEOUTPUT)) {
        for (i = 0; i < nfields; ++i) {
            field* f = &fields[i];
            char varname[NC_MAX_NAME];
            char fname[MAXSTRLEN];

            model_getbgfname(m, das->bgdir, f->varname, fname);
            strncpy(varname, f->varname, NC_MAX_NAME);
            if (!(das->updatespec & UPDATE_OUTPUTINC))
                strncat(varname, "_an", NC_MAX_NAME);
            else
                strncat(varname, "_inc", NC_MAX_NAME);
            model_writefield(m, fname, INT_MAX, varname, f->level, ((float***) fieldbuffer[i])[das->nmem][0]);
        }
    } else {
        for (i = 0; i < nfields; ++i) {
            field* f = &fields[i];
            char fname[MAXSTRLEN];

            model_getbgfname(m, das->bgdir, f->varname, fname);
            if (!(das->updatespec & UPDATE_OUTPUTINC))
                strncat(fname, ".analysis", MAXSTRLEN);
            else
                strncat(fname, ".increment", MAXSTRLEN);
            model_writefield(m, fname, INT_MAX, f->varname, f->level, ((float***) fieldbuffer[i])[das->nmem][0]);
        }
    }
}

/**
 */
static void das_writebg_toassemble(dasystem* das, int nfields, void** fieldbuffer, field fields[])
{
    char fname[MAXSTRLEN];
    int ni, nj;
    int i;

    model_getvardims(das->m, fields[0].varid, &ni, &nj, NULL);

    for (i = 0; i < nfields; ++i) {
        field* f = &fields[i];

        getfieldfname(das->bgdir, "bg", f->varname, f->level, fname);

        if (!file_exists(fname)) {
            int ncid;
            int dimids[2];
            int vid;

            ncw_create(fname, NC_CLOBBER | NETCDF_FORMAT, &ncid);
            ncw_def_dim(ncid, "nj", nj, &dimids[0]);
            ncw_def_dim(ncid, "ni", ni, &dimids[1]);
            ncw_def_var(ncid, f->varname, NC_FLOAT, 2, dimids, &vid);
            ncw_enddef(ncid);
            ncw_put_var_float(ncid, vid, ((float***) fieldbuffer[i])[das->nmem][0]);
            ncw_close(ncid);
        } else {
            writefield(fname, f->varname, 0, ((float***) fieldbuffer[i])[das->nmem][0]);
        }
    }
}

/** Writes `nfield' fields from `fieldbuffer' to disk.
 */
static void das_writebg(dasystem* das, int nfields, void** fieldbuffer, field fields[])
{
    if (das->updatespec & UPDATE_DIRECTWRITE)
        das_writebg_direct(das, nfields, fieldbuffer, fields);
    else
        das_writebg_toassemble(das, nfields, fieldbuffer, fields);
}

/** Allocates disk space for ensemble spread.
 */
static void das_allocatespread(dasystem* das, char fname[])
{
    model* m = das->m;
    int nvar = model_getnvar(m);
    char fname_src[MAXSTRLEN];
    int ncid, ncid_src;
    int vid;

    if (rank != 0)
        return;

    if (file_exists(fname))
        return;

    ncw_create(fname, NC_CLOBBER | NETCDF_FORMAT, &ncid);
    for (vid = 0; vid < nvar; ++vid) {
        char* varname_src = model_getvarname(m, vid);
        int varid_src;

        model_getmemberfname(m, das->ensdir, varname_src, 1, fname_src);
        ncw_open(fname_src, NC_NOWRITE, &ncid_src);
        ncw_copy_dims(ncid_src, ncid);
        ncw_inq_varid(ncid_src, varname_src, &varid_src);
        ncw_copy_vardef(ncid_src, varid_src, ncid);
        if (das->mode == MODE_ENKF) {
            char varname_dst[NC_MAX_NAME];

            strcpy(varname_dst, varname_src);
            strncat(varname_dst, "_an", NC_MAX_NAME);
            ncw_def_var_as(ncid, varname_src, varname_dst);
        }
        ncw_close(ncid_src);
    }
    ncw_close(ncid);
}

/** Allocates disk space for inflation magnitudes.
 */
static void das_allocateinflation(dasystem* das, char fname[])
{
    model* m = das->m;
    int nvar = model_getnvar(m);
    char fname_src[MAXSTRLEN];
    int ncid, ncid_src;
    int vid;

    if (rank != 0)
        return;

    if (file_exists(fname))
        return;

    ncw_create(fname, NC_CLOBBER | NETCDF_FORMAT, &ncid);
    for (vid = 0; vid < nvar; ++vid) {
        char* varname_src = model_getvarname(m, vid);
        int varid_src;

        model_getmemberfname(m, das->ensdir, varname_src, 1, fname_src);
        ncw_open(fname_src, NC_NOWRITE, &ncid_src);
        ncw_copy_dims(ncid_src, ncid);
        ncw_inq_varid(ncid_src, varname_src, &varid_src);
        ncw_copy_vardef(ncid_src, varid_src, ncid);
        ncw_close(ncid_src);
    }
    ncw_close(ncid);
}

/**
 */
static void das_writespread(dasystem* das, int nfields, void** fieldbuffer, field fields[], int isanalysis)
{
    char fname[MAXSTRLEN];
    model* m = das->m;
    int ni, nj;
    int fid, e, i, nv;
    double* v1 = NULL;
    double* v2 = NULL;
    float*** v_src = NULL;

    if (das->updatespec & UPDATE_DIRECTWRITE)
        strcpy(fname, FNAME_SPREAD);

    model_getvardims(m, fields[0].varid, &ni, &nj, NULL);
    nv = ni * nj;
    v1 = malloc(nv * sizeof(double));
    v2 = malloc(nv * sizeof(double));

    for (fid = 0; fid < nfields; ++fid) {
        field* f = &fields[fid];
        char varname[NC_MAX_NAME];

        v_src = (float***) fieldbuffer[fid];
        memset(v1, 0, nv * sizeof(double));
        memset(v2, 0, nv * sizeof(double));

        for (e = 0; e < das->nmem; ++e) {
            float* v = v_src[e][0];

            for (i = 0; i < nv; ++i) {
                v1[i] += v[i];
                v2[i] += v[i] * v[i];
            }
        }

        for (i = 0; i < nv; ++i) {
            v1[i] /= (double) das->nmem;
            v2[i] = v2[i] / (double) das->nmem - v1[i] * v1[i];
            v2[i] = (v2[i] < 0.0) ? 0.0 : sqrt(v2[i]);
            if (fabs(v2[i]) > (double) MAXOBSVAL)
                v2[i] = NaN;
        }

        strncpy(varname, f->varname, NC_MAX_NAME);
        if (isanalysis)
            strncat(varname, "_an", MAXSTRLEN);

        if (!(das->updatespec & UPDATE_DIRECTWRITE)) {  /* create file for *
                                                         * this field */
            int ncid, vid;
            int dimids[2];

            getfieldfname(das->mode == MODE_ENKF ? das->ensdir : das->bgdir, "spread", varname, f->level, fname);

            if (!file_exists(fname)) {
                ncw_create(fname, NC_CLOBBER | NETCDF_FORMAT, &ncid);
                ncw_def_dim(ncid, "nj", nj, &dimids[0]);
                ncw_def_dim(ncid, "ni", ni, &dimids[1]);
                ncw_def_var(ncid, varname, NC_FLOAT, 2, dimids, &vid);
                ncw_enddef(ncid);
            } else {
                ncw_open(fname, NC_WRITE, &ncid);
                if (!ncw_var_exists(ncid, varname)) {
                    ncw_redef(ncid);
                    ncw_inq_dimid(ncid, "nj", &dimids[0]);
                    ncw_inq_dimid(ncid, "ni", &dimids[1]);
                    ncw_def_var(ncid, varname, NC_FLOAT, 2, dimids, NULL);
                    ncw_enddef(ncid);
                }
                ncw_inq_varid(ncid, varname, &vid);
            }
            ncw_put_var_double(ncid, vid, v2);
            ncw_close(ncid);
        } else {
            float* v = calloc(nv, sizeof(float));

            for (i = 0; i < nv; ++i)
                v[i] = (float) v2[i];

            model_writefield(m, fname, INT_MAX, varname, f->level, v);
            free(v);
        }
    }

    free(v1);
    free(v2);
}

/**
 */
static void das_assemblemembers(dasystem* das)
{
    model* m = das->m;
    int nvar = model_getnvar(m);
    int i, e;

#if defined(MPI)
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    distribute_iterations(0, das->nmem - 1, nprocesses, rank, "    ");

    for (i = 0; i < nvar; ++i) {
        char* varname = model_getvarname(m, i);
        char varname_dst[NC_MAX_NAME];
        char fname_dst[MAXSTRLEN];
        int nlev, k;
        int ni, nj;
        float* v = NULL;

        enkf_printf("    %s:", varname);
        model_getmemberfname(m, das->ensdir, varname, 1, fname_dst);
        nlev = getnlevels(fname_dst, varname);
        strncpy(varname_dst, varname, NC_MAX_NAME);

        model_getvardims(m, i, &ni, &nj, NULL);
        v = malloc(ni * nj * sizeof(float));

        for (e = my_first_iteration; e <= my_last_iteration; ++e) {
            model_getmemberfname(m, das->ensdir, varname, e + 1, fname_dst);
            if (das->updatespec & UPDATE_SEPARATEOUTPUT) {
                if (!(das->updatespec & UPDATE_OUTPUTINC))
                    strncat(fname_dst, ".analysis", MAXSTRLEN);
                else
                    strncat(fname_dst, ".increment", MAXSTRLEN);
            } else {
                if (!(das->updatespec & UPDATE_OUTPUTINC))
                    strncat(varname_dst, "_an", NC_MAX_NAME);
                else
                    strncat(varname_dst, "_inc", NC_MAX_NAME);
            }

            for (k = 0; k < nlev; ++k) {
                char fname_src[MAXSTRLEN];
                int ncid_src, vid_src;
                size_t start[3] = { e, 0, 0 };
                size_t count[3] = { 1, nj, ni };

                getfieldfname(das->ensdir, "ens", varname, k, fname_src);
                ncw_open(fname_src, NC_NOWRITE, &ncid_src);
                ncw_inq_varid(ncid_src, varname, &vid_src);
                ncw_get_vara_float(ncid_src, vid_src, start, count, v);
                ncw_close(ncid_src);

                model_writefield(m, fname_dst, INT_MAX, varname_dst, k, v);
            }
            enkf_printf(".");
        }
        free(v);
        enkf_printf("\n");
    }

#if defined(MPI)
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    /*
     * remove tiles 
     */
    if (!(das->updatespec & UPDATE_LEAVETILES) && rank == 0) {
        for (i = 0; i < nvar; ++i) {
            char* varname = model_getvarname(m, i);
            char fname[MAXSTRLEN];
            int nlev, k;

            model_getmemberfname(m, das->ensdir, varname, 1, fname);
            nlev = getnlevels(fname, varname);
            for (k = 0; k < nlev; ++k) {
                getfieldfname(das->ensdir, "ens", varname, k, fname);
                file_delete(fname);
            }
        }
    }
}

/**
 */
static void das_assemblebg(dasystem* das)
{
    model* m = das->m;
    int nvar = model_getnvar(m);
    int i;

#if defined(MPI)
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    if (rank > 0)
        return;

    for (i = 0; i < nvar; ++i) {
        char* varname = model_getvarname(m, i);
        char varname_dst[NC_MAX_NAME];
        char fname_dst[MAXSTRLEN];
        int nlev, k;
        int ni, nj;
        float* v = NULL;

        enkf_printf("    %s:", varname);
        model_getbgfname(m, das->bgdir, varname, fname_dst);
        nlev = getnlevels(fname_dst, varname);
        strncpy(varname_dst, varname, NC_MAX_NAME);

        model_getvardims(m, i, &ni, &nj, NULL);
        v = malloc(ni * nj * sizeof(float));

        if (das->updatespec & UPDATE_SEPARATEOUTPUT) {
            if (!(das->updatespec & UPDATE_OUTPUTINC))
                strncat(fname_dst, ".analysis", MAXSTRLEN);
            else
                strncat(fname_dst, ".increment", MAXSTRLEN);
        } else {
            if (!(das->updatespec & UPDATE_OUTPUTINC))
                strncat(varname_dst, "_an", NC_MAX_NAME);
            else
                strncat(varname_dst, "_inc", NC_MAX_NAME);
        }

        for (k = 0; k < nlev; ++k) {
            char fname_src[MAXSTRLEN];
            int ncid_src, vid_src;

            getfieldfname(das->bgdir, "bg", varname, k, fname_src);
            ncw_open(fname_src, NC_NOWRITE, &ncid_src);
            ncw_inq_varid(ncid_src, varname, &vid_src);
            ncw_get_var_float(ncid_src, vid_src, v);
            ncw_close(ncid_src);

            model_writefield(m, fname_dst, INT_MAX, varname_dst, k, v);
            if (!(das->updatespec & UPDATE_LEAVETILES))
                file_delete(fname_src);

            enkf_printf(".");
        }
        free(v);
        enkf_printf("\n");
    }
}

/**
 */
static void das_assemblespread(dasystem* das)
{
    model* m = das->m;
    int nvar = model_getnvar(m);
    int i;

    for (i = 0; i < nvar; ++i) {
        char* varname = model_getvarname(m, i);
        char varname_an[NC_MAX_NAME];
        int nlev, k;
        int ni, nj;
        float* v = NULL;

        enkf_printf("    %s:", varname);
        nlev = getnlevels(FNAME_SPREAD, varname);
        if (das->mode == MODE_ENKF) {
            strncpy(varname_an, varname, NC_MAX_NAME);
            strncat(varname_an, "_an", NC_MAX_NAME);
        }

        model_getvardims(m, i, &ni, &nj, NULL);
        v = malloc(ni * nj * sizeof(float));

        for (k = 0; k < nlev; ++k) {
            char fname_src[MAXSTRLEN];
            int ncid_src, vid;

            getfieldfname(das->mode == MODE_ENKF ? das->ensdir : das->bgdir, "spread", varname, k, fname_src);
            ncw_open(fname_src, NC_NOWRITE, &ncid_src);

            ncw_inq_varid(ncid_src, varname, &vid);
            ncw_get_var_float(ncid_src, vid, v);
            ncw_close(ncid_src);
            file_delete(fname_src);

            model_writefield(m, FNAME_SPREAD, INT_MAX, varname, k, v);

            if (das->mode == MODE_ENKF) {
                getfieldfname(das->ensdir, "spread", varname_an, k, fname_src);
                ncw_open(fname_src, NC_NOWRITE, &ncid_src);
                ncw_inq_varid(ncid_src, varname_an, &vid);
                ncw_get_var_float(ncid_src, vid, v);
                ncw_close(ncid_src);
                file_delete(fname_src);

                model_writefield(m, FNAME_SPREAD, INT_MAX, varname_an, k, v);
            }

            enkf_printf(".");
        }
        free(v);
        enkf_printf("\n");
    }
}

/**
 */
static void das_assembleinflation(dasystem* das)
{
    model* m = das->m;
    int nvar = model_getnvar(m);
    int i;

    assert(das->mode == MODE_ENKF);

    for (i = 0; i < nvar; ++i) {
        char* varname = model_getvarname(m, i);
        int nlev, k;
        int ni, nj;
        float* v = NULL;

        enkf_printf("    %s:", varname);
        nlev = getnlevels(FNAME_INFLATION, varname);

        model_getvardims(m, i, &ni, &nj, NULL);
        v = malloc(ni * nj * sizeof(float));

        for (k = 0; k < nlev; ++k) {
            char fname_src[MAXSTRLEN];
            int ncid_src, vid;

            getfieldfname(das->ensdir, "inflation", varname, k, fname_src);
            ncw_open(fname_src, NC_NOWRITE, &ncid_src);

            ncw_inq_varid(ncid_src, varname, &vid);
            ncw_get_var_float(ncid_src, vid, v);
            ncw_close(ncid_src);
            model_writefield(m, FNAME_INFLATION, INT_MAX, varname, k, v);
            file_delete(fname_src);

            enkf_printf(".");
        }
        free(v);
        enkf_printf("\n");
    }
}

/** Updates ensemble/background by using calculated transform 
 * matrices/coefficients.
 */
void das_update(dasystem* das)
{
    model* m = das->m;
    int ngrid = model_getngrid(m);
    int nvar = model_getnvar(m);
    int gid;
    int i, e;

    if (das->nmem <= 0)
        das_setnmem(das);
    enkf_printf("    %d members\n", das->nmem);

    if (das->updatespec & UPDATE_DOSPREAD && rank == 0) {
        enkf_printf("    allocating disk space for spread:");
        das_allocatespread(das, FNAME_SPREAD);
        enkf_printf("\n");
    }

    if (das->updatespec & UPDATE_DOINFLATION && rank == 0) {
        enkf_printf("    allocating disk space for inflation:");
        das_allocateinflation(das, FNAME_INFLATION);
        enkf_printf("\n");
    }

    if (das->updatespec & UPDATE_DOPLOGS && rank == 0 && das->nplogs > 0) {
        enkf_printf("    defining state variables in point logs:");
        plog_definestatevars(das);
        enkf_printf("\n");
    }

    if (das->updatespec & UPDATE_DOFIELDS) {
        if (das->mode == MODE_ENKF) {
            distribute_iterations(0, das->nmem - 1, nprocesses, rank, "    ");
#if defined(MPI)
            MPI_Barrier(MPI_COMM_WORLD);
#endif
            enkf_printtime("    ");
            enkf_printf("    allocating disk space for analysis:");
            enkf_flush();
#if defined(MPI)
            MPI_Barrier(MPI_COMM_WORLD);
#endif
            if (!(das->updatespec & UPDATE_SEPARATEOUTPUT)) {
                for (i = 0; i < nvar; ++i) {
                    for (e = my_first_iteration; e <= my_last_iteration; ++e) {
                        char* varname_src = model_getvarname(m, i);
                        char fname[MAXSTRLEN];
                        int ncid;
                        char varname_dst[NC_MAX_NAME];

                        strncpy(varname_dst, varname_src, NC_MAX_NAME);
                        if (!(das->updatespec & UPDATE_OUTPUTINC))
                            strncat(varname_dst, "_an", NC_MAX_NAME);
                        else
                            strncat(varname_dst, "_inc", NC_MAX_NAME);

                        model_getmemberfname(m, das->ensdir, varname_src, e + 1, fname);
                        ncw_open(fname, NC_WRITE, &ncid);
                        if (!ncw_var_exists(ncid, varname_dst)) {
                            ncw_redef(ncid);
                            ncw_def_var_as(ncid, varname_src, varname_dst);
                        }
                        ncw_close(ncid);
                        printf(".");
                        fflush(stdout);
                    }
                }
            } else {
                for (i = 0; i < nvar; ++i) {
                    char* varname = model_getvarname(m, i);

                    for (e = my_first_iteration; e <= my_last_iteration; ++e) {
                        char fname_f[MAXSTRLEN], fname_a[MAXSTRLEN];
                        int ncid_f, ncid_a;
                        int vid_f;

                        model_getmemberfname(m, das->ensdir, varname, e + 1, fname_f);
                        ncw_open(fname_f, NC_NOWRITE, &ncid_f);

                        strncpy(fname_a, fname_f, MAXSTRLEN);
                        if (!(das->updatespec & UPDATE_OUTPUTINC))
                            strncat(fname_a, ".analysis", MAXSTRLEN);
                        else
                            strncat(fname_a, ".increment", MAXSTRLEN);
                        if (file_exists(fname_a)) {
                            ncw_open(fname_a, NC_WRITE, &ncid_a);
                            if (ncw_var_exists(ncid_a, varname)) {
                                ncw_close(ncid_a);
                                ncw_close(ncid_f);
                                continue;
                            }
                            ncw_redef(ncid_a);
                        } else {
                            ncw_create(fname_a, NC_CLOBBER | NETCDF_FORMAT, &ncid_a);
                            ncw_copy_dims(ncid_f, ncid_a);
                        }

                        ncw_inq_varid(ncid_f, varname, &vid_f);
                        ncw_copy_vardef(ncid_f, vid_f, ncid_a);
                        ncw_close(ncid_a);
                        ncw_close(ncid_f);
                        printf(".");
                        fflush(stdout);
                    }
                }
            }
#if defined(MPI)
            MPI_Barrier(MPI_COMM_WORLD);
#endif
            enkf_printf("\n");
            enkf_printtime("    ");
            enkf_flush();
        } else if (das->mode == MODE_ENOI) {
            if (rank == 0) {
                enkf_printf("    allocating disk space for the analysis:");
                enkf_printtime("    ");
                enkf_flush();

                if (!(das->updatespec & UPDATE_SEPARATEOUTPUT)) {
                    for (i = 0; i < nvar; ++i) {
                        char* varname_src = model_getvarname(m, i);
                        char fname[MAXSTRLEN];
                        int ncid;
                        char varname_dst[NC_MAX_NAME];

                        strncpy(varname_dst, varname_src, NC_MAX_NAME);
                        if (!(das->updatespec & UPDATE_OUTPUTINC))
                            strncat(varname_dst, "_an", NC_MAX_NAME);
                        else
                            strncat(varname_dst, "_inc", NC_MAX_NAME);
                        model_getbgfname(m, das->bgdir, varname_src, fname);
                        ncw_open(fname, NC_WRITE, &ncid);
                        if (!ncw_var_exists(ncid, varname_dst)) {
                            ncw_redef(ncid);
                            ncw_def_var_as(ncid, varname_src, varname_dst);
                        }
                        ncw_close(ncid);
                        printf(".");
                        fflush(stdout);
                    }
                } else {
                    for (i = 0; i < nvar; ++i) {
                        char* varname = model_getvarname(m, i);
                        char fname_f[MAXSTRLEN], fname_a[MAXSTRLEN];
                        int ncid_f, ncid_a;
                        int vid_f;

                        model_getbgfname(m, das->bgdir, varname, fname_f);
                        ncw_open(fname_f, NC_NOWRITE, &ncid_f);

                        strncpy(fname_a, fname_f, MAXSTRLEN);
                        if (!(das->updatespec & UPDATE_OUTPUTINC))
                            strncat(fname_a, ".analysis", MAXSTRLEN);
                        else
                            strncat(fname_a, ".increment", MAXSTRLEN);

                        if (file_exists(fname_a)) {
                            ncw_open(fname_a, NC_WRITE, &ncid_a);
                            if (ncw_var_exists(ncid_a, varname)) {
                                ncw_close(ncid_a);
                                ncw_close(ncid_f);
                                continue;
                            }
                            ncw_redef(ncid_a);
                        } else {
                            ncw_create(fname_a, NC_CLOBBER | NETCDF_FORMAT, &ncid_a);
                            ncw_copy_dims(ncid_f, ncid_a);
                        }

                        ncw_inq_varid(ncid_f, varname, &vid_f);
                        ncw_copy_vardef(ncid_f, vid_f, ncid_a);
                        ncw_close(ncid_a);
                        ncw_close(ncid_f);
                        printf(".");
                        fflush(stdout);
                    }
                }

                enkf_printf("\n");
                enkf_printtime("    ");
                enkf_flush();
            }
        }
    }

    for (gid = 0; gid < ngrid; ++gid) {
        void* grid = model_getgridbyid(m, gid);
        int nfields = 0;
        field* fields = NULL;
        void** fieldbuffer = NULL;
        int mni, mnj;
        int i, e;

        enkf_printf("    updating fields for %s:\n", grid_getname(grid));

        grid_getdims(grid, &mni, &mnj, NULL);

        das_getfields(das, gid, &nfields, &fields);
        enkf_printf("      %d fields\n", nfields);

        if (nfields == 0)
            continue;

        distribute_iterations(0, nfields - 1, nprocesses, rank, "      ");

        fieldbuffer = malloc(das->fieldbufsize * sizeof(void*));
        if (das->mode == MODE_ENKF) {
            for (i = 0; i < das->fieldbufsize; ++i)
                fieldbuffer[i] = alloc3d(das->nmem, mnj, mni, sizeof(float));
        } else if (das->mode == MODE_ENOI) {
            for (i = 0; i < das->fieldbufsize; ++i)
                fieldbuffer[i] = alloc3d(das->nmem + 1, mnj, mni, sizeof(float));
        }

        enkf_flush();
#if defined(MPI)
        MPI_Barrier(MPI_COMM_WORLD);
#endif
        for (i = my_first_iteration; i <= my_last_iteration; ++i) {
            int bufindex = (i - my_first_iteration) % das->fieldbufsize;
            field* f = &fields[i];
            char fname[MAXSTRLEN];

            printf("      %-8s %-3d (%d: %d: %.1f%%)\n", f->varname, f->level, rank, i, 100.0 * (double) (i - my_first_iteration + 1) / (double) (my_last_iteration - my_first_iteration + 1));
            fflush(stdout);

            for (e = 0; e < das->nmem; ++e) {
                model_getmemberfname(m, das->ensdir, f->varname, e + 1, fname);
                model_readfield(das->m, fname, INT_MAX, f->varname, f->level, ((float***) fieldbuffer[bufindex])[e][0]);
            }
            if (das->mode == MODE_ENOI) {
                if (!(das->updatespec & UPDATE_OUTPUTINC) && das->updatespec & (UPDATE_DOFIELDS | UPDATE_DOPLOGS)) {
                    model_getbgfname(m, das->bgdir, f->varname, fname);
                    model_readfield(das->m, fname, INT_MAX, f->varname, f->level, ((float***) fieldbuffer[bufindex])[das->nmem][0]);
                } else
                    memset(((float***) fieldbuffer[bufindex])[das->nmem][0], 0, mni * mnj * sizeof(float));
            }

            if (bufindex == das->fieldbufsize - 1 || i == my_last_iteration) {
                /*
                 * write forecast spread
                 */
                if (das->updatespec & UPDATE_DOSPREAD)
                    das_writespread(das, bufindex + 1, fieldbuffer, &fields[i - bufindex], 0);
                /*
                 * write forecast variables to point logs
                 */
                if (das->updatespec & UPDATE_DOPLOGS && das->nplogs > 0)
                    plog_writestatevars(das, bufindex + 1, fieldbuffer, &fields[i - bufindex], 0);

                if (das->mode == MODE_ENKF) {
                    das_updatefields(das, bufindex + 1, fieldbuffer, &fields[i - bufindex]);
                    if (das->updatespec & UPDATE_DOFIELDS)
                        das_writefields(das, bufindex + 1, fieldbuffer, &fields[i - bufindex]);
                    else if (i == my_last_iteration)
                        enkf_printf("      (skip writing the fields)\n");
                } else if (das->mode == MODE_ENOI) {
                    if (das->updatespec & (UPDATE_DOFIELDS | UPDATE_DOPLOGS))
                        das_updatebg(das, bufindex + 1, fieldbuffer, &fields[i - bufindex]);
                    if (das->updatespec & UPDATE_DOFIELDS)
                        das_writebg(das, bufindex + 1, fieldbuffer, &fields[i - bufindex]);
                }

                /*
                 * write analysis spread
                 */
                if (das->updatespec & UPDATE_DOSPREAD && das->mode == MODE_ENKF)
                    das_writespread(das, bufindex + 1, fieldbuffer, &fields[i - bufindex], 1);
                /*
                 * write analysis variables to point logs
                 */
                if (das->updatespec & UPDATE_DOPLOGS && das->nplogs > 0)
                    plog_writestatevars(das, bufindex + 1, fieldbuffer, &fields[i - bufindex], 1);
            }
        }

        free(fields);
        for (i = 0; i < das->fieldbufsize; ++i)
            free(fieldbuffer[i]);
        free(fieldbuffer);

        enkf_flush();
    }                           /* for gid */

#if defined(MPI)
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    if (!(das->updatespec & UPDATE_DIRECTWRITE)) {
        enkf_printtime("  ");
        if (das->updatespec & UPDATE_DOFIELDS) {
            enkf_printf("  assembling analysis:\n");
            if (das->mode == MODE_ENKF)
                das_assemblemembers(das);
            else if (das->mode == MODE_ENOI)
                das_assemblebg(das);
        }
        if (das->updatespec & UPDATE_DOSPREAD && rank == 0) {
            enkf_printf("  assembling spread:\n");
            das_assemblespread(das);
        }
        if (das->updatespec & UPDATE_DOINFLATION && rank == 0) {
            enkf_printf("  assembling inflation:\n");
            das_assembleinflation(das);
        }
        if (das->updatespec & UPDATE_DOPLOGS && das->nplogs > 0) {
            enkf_printf("  assembling state variables in point logs:\n");
            plog_assemblestatevars(das);
        }
    }
}

/** Calculates and writes to disk 3D field of correlation coefficients between
 * surface and other layers of a 3D variable. This function is not used at the
 * moment.
 */
void das_calccorr(dasystem* das, int mvid)
{
    model* m = das->m;
    char* varname = model_getvarname(m, mvid);
    int ncid, vid;
    int dimids[2];
    float*** v = NULL;
    float*** v0 = NULL;
    double* cor = NULL;
    double* std0 = NULL;
    double* std = NULL;
    int ni, nj, nk;
    int ktop, nv;
    int e, k, i;

    model_getvardims(m, mvid, &ni, &nj, &nk);
    nv = ni * nj;
    ktop = grid_gettoplayerid(model_getvargrid(m, mvid));

    v = alloc3d(das->nmem, nj, ni, sizeof(float));
    v0 = alloc3d(das->nmem, nj, ni, sizeof(float));
    cor = calloc(nv, sizeof(double));
    std = calloc(nv, sizeof(double));
    std0 = calloc(nv, sizeof(double));

    for (e = 0; e < das->nmem; ++e) {
        char fname[MAXSTRLEN];

        model_getmemberfname(m, das->ensdir, varname, e + 1, fname);
        model_readfield(das->m, fname, INT_MAX, varname, ktop, v0[e][0]);
    }
    for (i = 0; i < nv; ++i) {
        double vmean = 0.0;

        for (e = 0; e < das->nmem; ++e)
            vmean += (double) v0[e][0][i];
        vmean /= (double) das->nmem;
        for (e = 0; e < das->nmem; ++e)
            v0[e][0][i] -= (float) vmean;
        std0[i] = 0.0;
        for (e = 0; e < das->nmem; ++e)
            std0[i] += (double) (v0[e][0][i] * v0[e][0][i]);
        std0[i] = sqrt(std0[i] / (double) (das->nmem - 1));
    }

    distribute_iterations(0, nk - 1, nprocesses, rank, "    ");
#if defined(MPI)
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    for (k = my_first_iteration; k <= my_last_iteration; ++k) {
        char fname_src[MAXSTRLEN];
        char fname_dst[MAXSTRLEN];

        for (e = 0; e < das->nmem; ++e) {
            model_getmemberfname(m, das->ensdir, varname, e + 1, fname_src);
            model_readfield(das->m, fname_src, INT_MAX, varname, k, v[e][0]);
        }
        for (i = 0; i < nv; ++i) {
            double vmean = 0.0;

            for (e = 0; e < das->nmem; ++e)
                vmean += (double) v[e][0][i];
            vmean /= (double) das->nmem;
            for (e = 0; e < das->nmem; ++e)
                v[e][0][i] -= (float) vmean;
            std[i] = 0.0;
            for (e = 0; e < das->nmem; ++e)
                std[i] += (double) (v[e][0][i] * v[e][0][i]);
            std[i] /= sqrt(std[i] / (double) (das->nmem - 1));
            cor[i] = 0.0;
            for (e = 0; e < das->nmem; ++e)
                cor[i] += (double) (v[e][0][i] * v0[e][0][i]);
            cor[i] /= (std[i] * std0[i]);
        }

        snprintf(fname_dst, MAXSTRLEN, "corr_%s-%03d.nc", varname, k);
        if (!file_exists(fname_dst)) {
            ncw_create(fname_dst, NC_CLOBBER | NETCDF_FORMAT, &ncid);
            ncw_def_dim(ncid, "nj", nj, &dimids[0]);
            ncw_def_dim(ncid, "ni", ni, &dimids[1]);
            ncw_def_var(ncid, varname, NC_FLOAT, 2, dimids, &vid);
            ncw_enddef(ncid);
        } else {
            ncw_open(fname_dst, NC_WRITE, &ncid);
            if (!ncw_var_exists(ncid, varname)) {
                ncw_redef(ncid);
                ncw_inq_dimid(ncid, "nj", &dimids[0]);
                ncw_inq_dimid(ncid, "ni", &dimids[1]);
                ncw_def_var(ncid, varname, NC_FLOAT, 2, dimids, NULL);
                ncw_enddef(ncid);
            }
            ncw_inq_varid(ncid, varname, &vid);
        }
        ncw_put_var_double(ncid, vid, cor);
        ncw_close(ncid);
    }

    free(v);
    free(v0);
    free(std);
    free(std0);
    free(cor);

#if defined(MPI)
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    if (rank == 0) {
        float* v = malloc(nv * sizeof(float));
        char fname_dst[MAXSTRLEN];
        char fname_src[MAXSTRLEN];
        int ncid_dst, ncid_src;

        snprintf(fname_dst, MAXSTRLEN, "corr_%s.nc", varname);
        if (!file_exists(fname_dst)) {
            ncw_create(fname_dst, NC_CLOBBER | NETCDF_FORMAT, &ncid_dst);

            model_getmemberfname(m, das->ensdir, varname, 1, fname_src);
            ncw_open(fname_src, NC_NOWRITE, &ncid_src);
            ncw_copy_var(ncid_src, varname, ncid_dst);
            ncw_close(ncid_dst);
        }

        for (k = 0; k < nk; ++k) {
            snprintf(fname_src, MAXSTRLEN, "corr_%s-%03d.nc", varname, k);

            ncw_open(fname_src, NC_NOWRITE, &ncid_src);
            ncw_inq_varid(ncid_src, varname, &vid);
            ncw_get_var_float(ncid_src, vid, v);
            ncw_close(ncid_src);
            file_delete(fname_src);

            model_writefield(m, fname_dst, INT_MAX, varname, k, v);
        }

        free(v);
    }
#if defined(MPI)
    MPI_Barrier(MPI_COMM_WORLD);
#endif
}
