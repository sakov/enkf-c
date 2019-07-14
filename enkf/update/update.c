/******************************************************************************
 *
 * File:        update.c        
 *
 * Created:     11/2013
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description: The code to update model fields using ensmeble transforms
 *              calculated by CALC. It also does a number of auxiliary tasks
 *              such as calculating and writing ensemble spread, inflation,
 *              point logs etc.
 *
 * Revisions:   PS 15/07/2019: Moved procedures related to spread, inflation,
 *                and vertical correlations to diag.c.
 *
 *****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <math.h>
#include <time.h>
#include "definitions.h"
#include "utils.h"
#include "distribute.h"
#include "grid.h"
#include "lapack.h"
#include "dasystem.h"
#include "pointlog.h"
#include "diags.h"

#define EPSF 1.0e-6f

/**
 */
static void addnoise(dasystem* das, int varid, float*** v)
{
    model* m = das->m;
    double deflation = model_getvardeflation(m, varid);
    double sigma = model_getvarsigma(m, varid);
    double* random = NULL;
    double sum = 0.0;
    float mult;
    int ni, nj;
    int e, i;

    assert(isfinite(deflation));

    random = malloc((das->nmem + 1) / 2 * 2 * sizeof(double));
    for (e = 0; e < (das->nmem + 1) / 2; ++e)
        get_normalpair(&random[e * 2]);
    for (e = 0; e < das->nmem; ++e)
        sum += random[e];
    sum /= (double) das->nmem;
    for (e = 0; e < das->nmem; ++e)
        random[e] -= sum;

    mult = (float) (sqrt(1.0 - deflation * deflation) * sigma);
    model_getvargridsize(m, varid, &ni, &nj, NULL);
    for (e = 0; e < das->nmem; ++e) {
        float* vv = v[e][0];

        for (i = 0; i < ni * nj; ++i)
            vv[i] = vv[i] * (float) deflation + mult * random[e];
    }

    free(random);
}

/** Updates `nfield' fields read into `fieldbuffer' with `X5'. Applies
 * variable-dependent inflation to each field.
 */
static void das_updatefields(dasystem* das, int nfields, void** fieldbuffer, field fields[])
{
    model* m = das->m;
    int nmem = das->nmem;

    void* grid = model_getvargrid(m, fields[0].varid);
    int stride = grid_getstride(grid);
    int** nlevels = grid_getnumlevels(grid);
    int surfk = grid_getsurflayerid(grid);
    int periodic_i = grid_isperiodic_i(grid);
    int writeinflation = das->updatespec & UPDATE_DOINFLATION;

    /*
     * X5
     */
    char fname_X5[MAXSTRLEN];
    int ncid;
    int varid;
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

    grid_getsize(grid, &mni, &mnj, NULL);

    das_getfname_X5(das, grid, fname_X5);

    ncw_open(fname_X5, NC_NOWRITE, &ncid);
    ncw_inq_varid(ncid, "X5", &varid);
    ncw_inq_vardims(ncid, varid, 3, NULL, dimlens);
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
    if (stride > 1) {
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
        for (stepj = 0; stepj < stride && j < mnj; ++stepj, ++j) {
            if (stride == 1) {
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
                    float weight2 = (float) stepj / (float) stride;
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
                    for (stepi = 0; stepi < stride && i < mni; ++stepi, ++i) {
                        if (stepi == 0) {
                            memcpy(X5j[i], X5jj[ii], nmem * nmem * sizeof(float));
                        } else {
                            float weight2 = (float) stepi / (float) stride;
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
                float inflation0 = NAN;
                double inf_ratio = NAN;

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
                    if (surfk == 0) {
                        if (nlevels[j][i] <= f->level) {
                            if (das->updatespec & UPDATE_OUTPUTINC)
                                for (e = 0; e < nmem; ++e)
                                    vvv[e][j][i] = 0.0f;
                            continue;
                        }
                    } else {
                        if (nlevels[j][i] <= surfk - f->level) {
                            if (das->updatespec & UPDATE_OUTPUTINC)
                                for (e = 0; e < nmem; ++e)
                                    vvv[e][j][i] = 0.0f;
                            continue;
                        }
                    }
                    /*
                     * assume that if |value| > MAXOBSVAL, then it is filled
                     * with the missing value 
                     */
                    for (e = 0; e < nmem; ++e)
                        if (!isfinite(vvv[e][j][i]) || fabsf(vvv[e][j][i]) > (float) MAXOBSVAL)
                            break;
                    if (e < nmem) {
                        if (das->updatespec & UPDATE_OUTPUTINC)
                            for (e = 0; e < nmem; ++e)
                                vvv[e][j][i] = 0.0f;
                        continue;
                    }
                    /*
                     * (it would be straightforward to compare with the
                     * actual missing value, provided that it is NC_FLOAT;
                     * otherwise it may be a bit tiresome to handle all
                     * variations) 
                     */

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
                     * (do not inflate if inflation is about 1 or less than 1)
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
    if (stride > 1) {
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
        int varid = fields[fid].varid;

        if (!isnan(model_getvardeflation(m, varid))) {
            /*
             * (for now -- for 2D variables only) 
             */
            assert(fields[fid].level == grid_getsurflayerid(grid));

            addnoise(das, varid, fieldbuffer[fid]);
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
    int stride = grid_getstride(grid);
    int** nlevels = grid_getnumlevels(grid);
    int surfk = grid_getsurflayerid(grid);
    int periodic_i = grid_isperiodic_i(grid);

    char fname_w[MAXSTRLEN];
    int ncid;
    int varid;
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

    grid_getsize(grid, &mni, &mnj, NULL);

    das_getfname_w(das, grid, fname_w);

    ncw_open(fname_w, NC_NOWRITE, &ncid);
    ncw_inq_varid(ncid, "w", &varid);
    ncw_inq_vardims(ncid, varid, 3, NULL, dimlens);
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
    if (stride > 1) {
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
        for (stepj = 0; stepj < stride && j < mnj; ++stepj, ++j) {
            if (stride == 1) {
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
                    float weight2 = (float) stepj / stride;
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
                    for (stepi = 0; stepi < stride && i < mni; ++stepi, ++i) {
                        if (stepi == 0)
                            memcpy(wj[i], wjj[ii], nmem * sizeof(float));
                        else {
                            float weight2 = (float) stepi / stride;
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

                    if (surfk == 0) {
                        if (nlevels[j][i] <= f->level) {
                            if (das->updatespec & UPDATE_OUTPUTINC)
                                vvv[nmem][j][i] = 0.0f;
                            continue;
                        }
                    } else {
                        if (nlevels[j][i] <= surfk - f->level) {
                            if (das->updatespec & UPDATE_OUTPUTINC)
                                vvv[nmem][j][i] = 0.0f;
                            continue;
                        }
                    }
                    /*
                     * assume that if |value| > MAXOBSVAL, then it is filled
                     * with the missing value 
                     */
                    if (fabsf(vvv[nmem][j][i]) > (float) MAXOBSVAL) {
                        if (das->updatespec & UPDATE_OUTPUTINC)
                            vvv[nmem][j][i] = 0.0f;
                        continue;
                    }
                    for (e = 0; e < nmem; ++e)
                        if (fabsf(vvv[e][j][i]) > (float) MAXOBSVAL)
                            break;
                    if (e < nmem) {
                        if (das->updatespec & UPDATE_OUTPUTINC)
                            vvv[nmem][j][i] = 0.0f;
                        continue;
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
    if (stride > 1) {
        free(wjj);
        free(wjj1);
        free(wjj2);
    }
}

/** Write analyses/increments directly to the output NetCDF file.
 *
 *  The implication is that different CPUs write horizontal variables to
 *  the same file in parallel. In practice this works only with NC_CLASSIC_MODEL
 *  or NC_64BIT_OFFSET, and never 100% reliably, although the failure percent
 *  can be very very low.
 */
static void das_writefields_direct(dasystem* das, int nfields, void** fieldbuffer, field fields[])
{
    int i, e;

    assert(das->mode == MODE_ENKF);

    if (!(das->updatespec & UPDATE_SEPARATEOUTPUT)) {
        for (i = 0; i < nfields; ++i) {
            field* f = &fields[i];
            char varname[NC_MAX_NAME];

            strncpy(varname, f->varname, NC_MAX_NAME - 1);
            if (!(das->updatespec & UPDATE_OUTPUTINC))
                strncat(varname, "_an", NC_MAX_NAME);
            else
                strncat(varname, "_inc", NC_MAX_NAME);

            for (e = 0; e < das->nmem; ++e) {
                char fname[MAXSTRLEN];

                model_getmemberfname(das->m, das->ensdir, f->varname, e + 1, fname);
                model_writefieldas(das->m, fname, varname, f->varname, f->level, ((float***) fieldbuffer[i])[e][0]);
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
                model_writefield(das->m, fname, f->varname, f->level, ((float***) fieldbuffer[i])[e][0]);
            }
        }
    }
}

/** This is the main procedure for parallel settings. It writes each
 * horizontal layer (in parallel) to a separate file ("tile"). The tiles are
 * "assembled" later by a single process.
 */
static void das_writefields_toassemble(dasystem* das, int nfields, void** fieldbuffer, field fields[])
{
    char fname[MAXSTRLEN];
    int ni, nj;
    int i;

    model_getvargridsize(das->m, fields[0].varid, &ni, &nj, NULL);

    for (i = 0; i < nfields; ++i) {
        field* f = &fields[i];

        getfieldfname(das->ensdir, "ens", f->varname, f->level, fname);

        if (!file_exists(fname)) {
            int ncid;
            int dimids[3];
            int vid;

            ncw_create(fname, NC_CLOBBER | das->ncformat, &ncid);
            ncw_def_dim(ncid, "m", das->nmem, &dimids[0]);
            ncw_def_dim(ncid, "nj", nj, &dimids[1]);
            ncw_def_dim(ncid, "ni", ni, &dimids[2]);
            ncw_def_var(ncid, f->varname, NC_FLOAT, 3, dimids, &vid);
            if (das->nccompression > 0)
                ncw_def_deflate(ncid, 0, 1, das->nccompression);
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

/**
 */
static void das_writebg_direct(dasystem* das, int nfields, void** fieldbuffer, field fields[])
{
    model* m = das->m;
    int ni, nj;
    int i;

    assert(das->mode == MODE_ENOI);

    model_getvargridsize(m, fields[0].varid, &ni, &nj, NULL);

    if (!(das->updatespec & UPDATE_SEPARATEOUTPUT)) {
        for (i = 0; i < nfields; ++i) {
            field* f = &fields[i];
            char varname[NC_MAX_NAME];
            char fname[MAXSTRLEN];

            model_getbgfname(m, das->bgdir, f->varname, fname);
            strncpy(varname, f->varname, NC_MAX_NAME - 1);
            if (!(das->updatespec & UPDATE_OUTPUTINC))
                strncat(varname, "_an", NC_MAX_NAME);
            else
                strncat(varname, "_inc", NC_MAX_NAME);
            model_writefield(m, fname, varname, f->level, ((float***) fieldbuffer[i])[das->nmem][0]);
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
            model_writefield(m, fname, f->varname, f->level, ((float***) fieldbuffer[i])[das->nmem][0]);
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

    model_getvargridsize(das->m, fields[0].varid, &ni, &nj, NULL);

    for (i = 0; i < nfields; ++i) {
        field* f = &fields[i];

        getfieldfname(das->bgdir, "bg", f->varname, f->level, fname);

        if (!file_exists(fname)) {
            int ncid;
            int dimids[2];
            int vid;

            ncw_create(fname, NC_CLOBBER | das->ncformat, &ncid);
            ncw_def_dim(ncid, "nj", nj, &dimids[0]);
            ncw_def_dim(ncid, "ni", ni, &dimids[1]);
            ncw_def_var(ncid, f->varname, NC_FLOAT, 2, dimids, &vid);
            if (das->nccompression > 0)
                ncw_def_deflate(ncid, 0, 1, das->nccompression);
            ncw_enddef(ncid);
            ncw_put_var_float(ncid, vid, ((float***) fieldbuffer[i])[das->nmem][0]);
            ncw_close(ncid);
        } else {
            model_writefield(das->m, fname, f->varname, f->level, ((float***) fieldbuffer[i])[das->nmem][0]);
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
        strncpy(varname_dst, varname, NC_MAX_NAME - 1);

        model_getvargridsize(m, i, &ni, &nj, NULL);
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

                if (nlev > 1)
                    getfieldfname(das->ensdir, "ens", varname, k, fname_src);
                else
                    getfieldfname(das->ensdir, "ens", varname, grid_getsurflayerid(model_getvargrid(m, i)), fname_src);
                ncw_open(fname_src, NC_NOWRITE, &ncid_src);
                ncw_inq_varid(ncid_src, varname, &vid_src);
                ncw_get_vara_float(ncid_src, vid_src, start, count, v);
                ncw_close(ncid_src);

                model_writefield(m, fname_dst, varname_dst, k, v);
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
                if (nlev > 1)
                    getfieldfname(das->ensdir, "ens", varname, k, fname);
                else
                    getfieldfname(das->ensdir, "ens", varname, grid_getsurflayerid(model_getvargrid(m, i)), fname);
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
        strncpy(varname_dst, varname, NC_MAX_NAME - 1);

        model_getvargridsize(m, i, &ni, &nj, NULL);
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

            if (nlev > 1)
                getfieldfname(das->bgdir, "bg", varname, k, fname_src);
            else
                getfieldfname(das->bgdir, "bg", varname, grid_getsurflayerid(model_getvargrid(m, i)), fname_src);
            ncw_open(fname_src, NC_NOWRITE, &ncid_src);
            ncw_inq_varid(ncid_src, varname, &vid_src);
            ncw_get_var_float(ncid_src, vid_src, v);
            ncw_close(ncid_src);

            model_writefield(m, fname_dst, varname_dst, k, v);
            if (!(das->updatespec & UPDATE_LEAVETILES))
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

#if defined(MPI)
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    fflush(NULL);

    if (das->updatespec & UPDATE_DOSPREAD && rank == 0) {
        enkf_printf("    allocating disk space for spread:");
        das_allocatespread(das, FNAME_SPREAD);
        enkf_printf("\n");
        enkf_flush();
    }

    if (das->updatespec & UPDATE_DOINFLATION && rank == 0) {
        enkf_printf("    allocating disk space for inflation:");
        das_allocateinflation(das, FNAME_INFLATION);
        enkf_printf("\n");
        enkf_flush();
    }

    if (das->updatespec & UPDATE_DOPLOGS && rank == 0) {
        enkf_printf("    defining state variables in point logs:");
        plog_definestatevars(das);
        enkf_printf("\n");
        enkf_flush();
    }

    if (rank == 0)
        dir_createifabsent(DIRNAME_TMP);
    if (das->updatespec & UPDATE_DOFIELDS) {
        if (das->mode == MODE_ENKF) {
            distribute_iterations(0, das->nmem - 1, nprocesses, rank, "    ");

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

                        strncpy(varname_dst, varname_src, NC_MAX_NAME - 1);
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
                        if (enkf_verbose) {
                            printf(".");
                            fflush(stdout);
                        }
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
                        } else
                            ncw_create(fname_a, NC_CLOBBER | das->ncformat, &ncid_a);
                        ncw_inq_varid(ncid_f, varname, &vid_f);
                        ncw_copy_vardef(ncid_f, vid_f, ncid_a);
                        if (das->nccompression > 0)
                            ncw_def_deflate(ncid_a, 0, 1, das->nccompression);
                        ncw_close(ncid_a);
                        ncw_close(ncid_f);
                        if (enkf_verbose) {
                            printf(".");
                            fflush(stdout);
                        }
                    }
                }
            }
#if defined(MPI)
            MPI_Barrier(MPI_COMM_WORLD);
#endif
            enkf_printf("\n");
        } else if (das->mode == MODE_ENOI) {
            if (rank == 0) {
                enkf_printtime("    ");
                enkf_printf("    allocating disk space for analysis:");
                enkf_flush();

                if (!(das->updatespec & UPDATE_SEPARATEOUTPUT)) {
                    for (i = 0; i < nvar; ++i) {
                        char* varname_src = model_getvarname(m, i);
                        char fname[MAXSTRLEN];
                        int ncid;
                        char varname_dst[NC_MAX_NAME];

                        strncpy(varname_dst, varname_src, NC_MAX_NAME - 1);
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
                        if (enkf_verbose) {
                            printf(".");
                            fflush(stdout);
                        }
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
                        } else
                            ncw_create(fname_a, NC_CLOBBER | das->ncformat, &ncid_a);
                        ncw_inq_varid(ncid_f, varname, &vid_f);
                        ncw_copy_vardef(ncid_f, vid_f, ncid_a);
                        if (das->nccompression > 0)
                            ncw_def_deflate(ncid_a, 0, 1, das->nccompression);
                        ncw_close(ncid_a);
                        ncw_close(ncid_f);
                        if (enkf_verbose) {
                            printf(".");
                            fflush(stdout);
                        }
                    }
                }

                enkf_printf("\n");
            }
        }
    }

    for (gid = 0; gid < ngrid; ++gid) {
        void* grid = model_getgridbyid(m, gid);
        int nfields = 0;
        field* fields = NULL;
        void** fieldbuffer = NULL;
        int mni, mnj;
        int fid, i;

        enkf_printf("    processing fields for %s:\n", grid_getname(grid));

        enkf_printtime("      ");

        grid_getsize(grid, &mni, &mnj, NULL);

#if defined(MPI)
        MPI_Barrier(MPI_COMM_WORLD);
#endif
        das_getfields(das, gid, &nfields, &fields);
        enkf_printf("      %d fields\n", nfields);

        if (nfields == 0)
            continue;

        distribute_iterations(0, nfields - 1, nprocesses, rank, "      ");

        if (my_first_iteration > my_last_iteration)
            continue;

        fieldbuffer = malloc(das->fieldbufsize * sizeof(void*));
        if (das->mode == MODE_ENKF) {
            for (i = 0; i < das->fieldbufsize; ++i)
                fieldbuffer[i] = alloc3d(das->nmem, mnj, mni, sizeof(float));
        } else if (das->mode == MODE_ENOI) {
            for (i = 0; i < das->fieldbufsize; ++i)
                fieldbuffer[i] = alloc3d(das->nmem + 1, mnj, mni, sizeof(float));
        }

        for (fid = my_first_iteration; fid <= my_last_iteration; ++fid) {
            int bufid = (fid - my_first_iteration) % das->fieldbufsize;
            field* f = &fields[fid];
            char fname[MAXSTRLEN];
            int e;

            if (enkf_verbose) {
                printf("      %-8s %-3d (%d: %d: %.1f%%)\n", f->varname, f->level, rank, fid, 100.0 * (double) (fid - my_first_iteration + 1) / (double) (my_last_iteration - my_first_iteration + 1));
                fflush(stdout);
            }

            for (e = 0; e < das->nmem; ++e) {
                model_getmemberfname(m, das->ensdir, f->varname, e + 1, fname);
                model_readfield(das->m, fname, f->varname, f->level, ((float***) fieldbuffer[bufid])[e][0]);
            }

            /*
             * read the background to write it to pointlogs, regardless of
             * whether output is increment or analysis
             */
            if (das->mode == MODE_ENOI) {
                model_getbgfname(m, das->bgdir, f->varname, fname);
                model_readfield(das->m, fname, f->varname, f->level, ((float***) fieldbuffer[bufid])[das->nmem][0]);
            }

            if (bufid == das->fieldbufsize - 1 || fid == my_last_iteration) {
                /*
                 * write forecast spread
                 */
                if (das->updatespec & UPDATE_DOFORECASTSPREAD)
                    das_writespread(das, bufid + 1, fieldbuffer, &fields[fid - bufid], 0);

                /*
                 * write forecast variables to point logs
                 */
                if (das->updatespec & UPDATE_DOPLOGS)
                    plog_writestatevars(das, bufid + 1, fieldbuffer, &fields[fid - bufid], 0);

                /*
                 * now set the background to 0 if output is increment
                 */
                if (das->mode == MODE_ENOI && (das->updatespec & UPDATE_OUTPUTINC)) {
                    int ii;

                    for (ii = 0; ii <= bufid; ++ii)
                        memset(((float***) fieldbuffer[ii])[das->nmem][0], 0, mni * mnj * sizeof(float));
                }

                if (das->updatespec & (UPDATE_DOFIELDS | UPDATE_DOANALYSISSPREAD | UPDATE_DOPLOGS | UPDATE_DOINFLATION)) {
                    if (das->mode == MODE_ENKF) {
                        das_updatefields(das, bufid + 1, fieldbuffer, &fields[fid - bufid]);
                        if (das->updatespec & UPDATE_DOFIELDS)
                            das_writefields(das, bufid + 1, fieldbuffer, &fields[fid - bufid]);
                        else if (fid == my_last_iteration)
                            enkf_printf("      (skip writing the fields)\n");
                    } else if (das->mode == MODE_ENOI) {
                        if (das->updatespec & (UPDATE_DOFIELDS | UPDATE_DOPLOGS))
                            das_updatebg(das, bufid + 1, fieldbuffer, &fields[fid - bufid]);
                        if (das->updatespec & UPDATE_DOFIELDS)
                            das_writebg(das, bufid + 1, fieldbuffer, &fields[fid - bufid]);
                        else if (fid == my_last_iteration)
                            enkf_printf("      (skip writing the fields)\n");
                    }
                }

                /*
                 * write analysis spread
                 */
                if (das->updatespec & UPDATE_DOANALYSISSPREAD && das->mode == MODE_ENKF)
                    das_writespread(das, bufid + 1, fieldbuffer, &fields[fid - bufid], 1);
                /*
                 * write analysis variables to point logs
                 */
                if (das->updatespec & UPDATE_DOPLOGS)
                    plog_writestatevars(das, bufid + 1, fieldbuffer, &fields[fid - bufid], 1);
            }
        }                       /* for fid */

        for (i = 0; i < das->fieldbufsize; ++i)
            free(fieldbuffer[i]);
        free(fieldbuffer);
        free(fields);

        enkf_flush();
    }                           /* for gid */

#if defined(MPI)
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    if (!(das->updatespec & UPDATE_DIRECTWRITE)) {
        if (das->updatespec & UPDATE_DOFIELDS) {
            enkf_printtime("  ");
            enkf_printf("  assembling analysis:\n");
            if (das->mode == MODE_ENKF)
                das_assemblemembers(das);
            else if (das->mode == MODE_ENOI)
                das_assemblebg(das);
        }
        if (das->updatespec & UPDATE_DOSPREAD && rank == 0) {
            enkf_printtime("  ");
            enkf_printf("  assembling spread:\n");
            das_assemblespread(das);
        }
        if (das->updatespec & UPDATE_DOINFLATION && rank == 0) {
            enkf_printtime("  ");
            enkf_printf("  assembling inflation:\n");
            das_assembleinflation(das);
        }
        if (das->updatespec & UPDATE_DOPLOGS) {
            enkf_printtime("  ");
            enkf_printf("  assembling state variables in point logs:\n");
            plog_assemblestatevars(das);
        }
    }
#if defined(MPI)
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    if (rank == 0)
        dir_rmifexists(DIRNAME_TMP);
}
