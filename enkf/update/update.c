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
#include "ncutils.h"
#if defined(USE_MPIQUEUE)
#include "mpiqueue.h"
#endif
#include "distribute.h"
#include "grid.h"
#include "lapack.h"
#include "dasystem.h"
#include "pointlog.h"

#define EPSF 1.0e-6f
#define NTILE_INC 400

/**
 */
static void addnoise(dasystem* das, int varid, void* v)
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
    if (nj > 0) {
        for (e = 0; e < das->nmem; ++e) {
            float* vv = ((float***) v)[e][0];

            for (i = 0; i < ni * nj; ++i)
                vv[i] = vv[i] * (float) deflation + mult * random[e];
        }
    } else {
        for (e = 0; e < das->nmem; ++e) {
            float* vv = ((float**) v)[e];

            for (i = 0; i < ni; ++i)
                vv[i] = vv[i] * (float) deflation + mult * random[e];
        }
    }

    free(random);
}

/** Updates `nfield' fields read into `fieldbuffer' by applying transforms `X5'.
 ** Applies variable-dependent inflation to each field.
 */
static void das_updatefields(dasystem* das, int nfield, void** fieldbuffer, field fields[])
{
    model* m = das->m;
    int nmem = das->nmem;
    int nmem_dynamic = das->nmem_dynamic;

    void* g = model_getvargrid(m, fields[0].varid);
    int gridid = grid_getid(g);
    int stride = grid_getstride(g);
    void* nlevels = grid_getnumlevels(g);
    int surfk = grid_getsurflayerid(g);
    int periodic_i = grid_isperiodic_i(g);
    int writeinflation = das->updatespec & UPDATE_DOINFLATION;
    int isstructured = 1;

    /*
     * T
     */
    char fname[MAXSTRLEN];
    int ncid;
    int varid_T, varid_w;
    size_t dimlens[4];
    size_t start[4], count[4];
    float** Tjj = NULL;
    float** Tjj1 = NULL;
    float** Tjj2 = NULL;
    float** Tj = NULL;
    float** wjj = NULL;
    float** wjj1 = NULL;
    float** wjj2 = NULL;
    float** wj;

    float* v_f = NULL;          /* v_f = E_f(i, :) */
    float* v_a = NULL;          /* v_a = E_a(i, :) */
    void* infl = NULL;

    int mni, mnj;
    int ni, nj;
    int i, j;
    int jj, stepj, ii, stepi;
    int e, fid;

    assert(das->mode == MODE_ENKF || das->mode == MODE_HYBRID);

    grid_getsize(g, &mni, &mnj, NULL);
    /*
     * a treatment for unstructured grids
     */
    if (mnj <= 0) {
        mnj = mni;
        mni = 1;
        isstructured = 0;
    }

    das_getfname_transforms(das, gridid, fname);

    ncw_open(fname, NC_NOWRITE, &ncid);
    ncw_inq_varid(ncid, "T", &varid_T);
    ncw_inq_varid(ncid, "w", &varid_w);
    ncw_inq_vardims(ncid, varid_T, 4, NULL, dimlens);
    if (rank == 0) {
        int stride_transforms;

        ncw_get_att_int(ncid, NC_GLOBAL, "stride", &stride_transforms);
        if (stride_transforms != stride)
            enkf_quit("grid \"%s\": stride = %d is not equal to the stride of transforms = %d", grid_getname(g), stride, stride_transforms);
    }
    nj = dimlens[0];
    ni = dimlens[1];

    start[0] = fields[0].j1 / stride;
    start[1] = 0;
    start[2] = 0;
    start[3] = 0;
    count[0] = 1;
    count[1] = ni;
    count[2] = nmem_dynamic;
    count[3] = nmem;

    Tj = alloc2d(mni, nmem_dynamic * nmem, sizeof(float));
    wj = alloc2d(mni, nmem, sizeof(float));
    if (stride > 1) {
        Tjj = alloc2d(ni, nmem_dynamic * nmem, sizeof(float));
        Tjj1 = alloc2d(ni, nmem_dynamic * nmem, sizeof(float));
        Tjj2 = alloc2d(ni, nmem_dynamic * nmem, sizeof(float));
        ncw_get_vara_float(ncid, varid_T, start, count, Tjj2[0]);
        wjj = alloc2d(ni, nmem, sizeof(float));
        wjj1 = alloc2d(ni, nmem, sizeof(float));
        wjj2 = alloc2d(ni, nmem, sizeof(float));
        count[2] = nmem;
        ncw_get_vara_float(ncid, varid_w, start, count, wjj2[0]);
    }

    v_f = malloc(nmem * sizeof(float));
    v_a = malloc(nmem_dynamic * sizeof(float));
    if (writeinflation)
        infl = (isstructured) ? malloc(mni * sizeof(float)) : alloc2d(nfield, mnj, sizeof(float));

    /*
     * jj, ii are the indices of the subsampled grid; i, j are the indices
     * of the model grid 
     */
    assert(fields[0].j1 % stride == 0);
    for (jj = fields[0].j1 / stride, j = 0; jj <= fields[0].j2 / stride && jj < nj; ++jj) {
        for (stepj = 0; stepj < stride && j <= fields[0].j2 - fields[0].j1; ++stepj, ++j) {
            if (stride == 1) {
                /*
                 * no interpolation necessary; simply read the ETMs for the
                 * j-th row from disk 
                 */
                start[0] = jj;
                count[2] = nmem_dynamic;
                ncw_get_vara_float(ncid, varid_T, start, count, Tj[0]);
                count[2] = nmem;
                ncw_get_vara_float(ncid, varid_w, start, count, wj[0]);
            } else {
                /*
                 * the following code interpolates the ETM back to the
                 * original grid, first by j, and then by i 
                 */
                if (stepj == 0) {
                    memcpy(Tjj1[0], Tjj2[0], ni * nmem_dynamic * nmem * sizeof(float));
                    memcpy(Tjj[0], Tjj2[0], ni * nmem_dynamic * nmem * sizeof(float));
                    memcpy(wjj[0], wjj2[0], ni * nmem * sizeof(float));
                    memcpy(wjj1[0], wjj2[0], ni * nmem * sizeof(float));
                    if (jj < nj - 1) {
                        start[0] = (jj + 1) % nj;
                        count[2] = nmem_dynamic;
                        ncw_get_vara_float(ncid, varid_T, start, count, Tjj2[0]);
                        count[2] = nmem;
                        ncw_get_vara_float(ncid, varid_w, start, count, wjj2[0]);
                    }
                } else {
                    float weight2 = (float) stepj / (float) stride;
                    float weight1 = 1.0f - weight2;

                    for (ii = 0; ii < ni; ++ii) {
                        float* Tjjii = Tjj[ii];
                        float* Tjj1ii = Tjj1[ii];
                        float* Tjj2ii = Tjj2[ii];
                        float* wjjii = wjj[ii];
                        float* wjj1ii = wjj1[ii];
                        float* wjj2ii = wjj2[ii];

                        for (e = 0; e < nmem_dynamic * nmem; ++e)
                            Tjjii[e] = Tjj1ii[e] * weight1 + Tjj2ii[e] * weight2;
                        for (e = 0; e < nmem; ++e)
                            wjjii[e] = wjj1ii[e] * weight1 + wjj2ii[e] * weight2;
                    }
                }

                for (ii = 0, i = 0; ii < ni; ++ii) {
                    for (stepi = 0; stepi < stride && i < mni; ++stepi, ++i) {
                        if (stepi == 0) {
                            memcpy(Tj[i], Tjj[ii], nmem_dynamic * nmem * sizeof(float));
                            memcpy(wj[i], wjj[ii], nmem * sizeof(float));
                        } else {
                            float weight2 = (float) stepi / (float) stride;
                            float weight1 = 1.0f - weight2;
                            float* Tjjii1 = Tjj[ii];
                            float* Tji = Tj[i];
                            float* Tjjii2;
                            float* wjjii1 = wjj[ii];
                            float* wji = wj[i];
                            float* wjjii2;

                            if (ii < ni - 1) {
                                Tjjii2 = Tjj[ii + 1];
                                wjjii2 = wjj[ii + 1];
                            } else {
                                Tjjii2 = Tjj[(periodic_i) ? (ii + 1) % ni : ii];
                                wjjii2 = wjj[(periodic_i) ? (ii + 1) % ni : ii];
                            }
                            for (e = 0; e < nmem_dynamic * nmem; ++e)
                                Tji[e] = Tjjii1[e] * weight1 + Tjjii2[e] * weight2;
                            for (e = 0; e < nmem; ++e)
                                wji[e] = wjjii1[e] * weight1 + wjjii2[e] * weight2;
                        }
                    }
                }
            }                   /* stride != 1 */

            /*
             * (at this stage Tj should contain the array of T matrices
             * for the j-th row of the grid) 
             */

            /*
             * update the j-th row of the fields 
             */
            for (fid = 0; fid < nfield; ++fid) {
                field* f = &fields[fid];
                int applylog = model_getvarislog(m, f->varid);
                float inflation0 = NAN;
                double inf_ratio = NAN;
                float*** vvv = (isstructured) ? fieldbuffer[fid] : NULL;
                float** vv = (isstructured) ? NULL : fieldbuffer[fid];

                model_getvarinflation(m, f->varid, &inflation0, &inf_ratio);
                if (writeinflation && isstructured)
                    memset(infl, 0, mni * sizeof(float));

                for (i = 0; i < mni; ++i) {
                    float inflation = inflation0;
                    double v1_f, v1_a;

                    if (isstructured) {
                        if (surfk == 0) {
                            if (((int**) nlevels)[f->j1 + j][i] <= f->level) {
                                for (e = 0; e < nmem; ++e)
                                    vvv[e][j][i] = 0.0f;
                                continue;
                            }
                        } else {
                            if (((int**) nlevels)[f->j1 + j][i] <= ((f->issurfacevar) ? 0 : surfk - f->level)) {
                                for (e = 0; e < nmem; ++e)
                                    vvv[e][j][i] = 0.0f;
                                continue;
                            }
                        }
                    } else {
                        if (surfk == 0) {
                            if (((int*) nlevels)[f->j1 + j] <= f->level) {
                                for (e = 0; e < nmem; ++e)
                                    vv[e][j] = 0.0f;
                                continue;
                            }
                        } else {
                            if (((int*) nlevels)[f->j1 + j] <= ((f->issurfacevar) ? 0 : surfk - f->level)) {
                                for (e = 0; e < nmem; ++e)
                                    vv[e][j] = 0.0f;
                                continue;
                            }
                        }
                    }

                    /*
                     * Assume that if |value| > MAXOBSVAL, then it is filled
                     * with the missing value.
                     *
                     * (It would be straightforward to compare with the
                     * actual missing value, provided that it is NC_FLOAT;
                     * otherwise it may be a bit tiresome to handle all
                     * variations.) 
                     */
                    if (isstructured) {
                        for (e = 0; e < nmem; ++e)
                            if (!isfinite(vvv[e][j][i]) || fabsf(vvv[e][j][i]) > (float) MAXOBSVAL)
                                break;
                        if (e < nmem_dynamic) {
                            for (e = 0; e < nmem; ++e)
                                vvv[e][j][i] = 0.0f;
                            continue;
                        } else if (e < nmem) {
                            for (e = nmem_dynamic; e < nmem; ++e)
                                vvv[e][j][i] = 0.0f;
                            if (nmem_dynamic == 0)
                                continue;
                        }

                        for (e = 0; e < nmem; ++e)
                            v_f[e] = vvv[e][j][i];
                    } else {
                        for (e = 0; e < nmem; ++e)
                            if (!isfinite(vv[e][j]) || fabsf(vv[e][j]) > (float) MAXOBSVAL)
                                break;
                        if (e < nmem_dynamic) {
                            for (e = 0; e < nmem; ++e)
                                vv[e][j] = 0.0f;
                            continue;
                        } else if (e < nmem) {
                            for (e = nmem_dynamic; e < nmem; ++e)
                                vv[e][j] = 0.0f;
                            if (nmem_dynamic == 0)
                                continue;
                        }

                        for (e = 0; e < nmem; ++e)
                            v_f[e] = vv[e][j];
                    }

                    /*
                     * E^a(i, :) = E^f(i, :) * X5
                     *           = x^f(i) * 1^T + A^f(i, :) * (w * 1^T + T)
                     */

                    /*
                     * x^f
                     */
                    for (e = 0, v1_f = 0.0; e < nmem_dynamic; ++e)
                        v1_f += v_f[e];
                    v1_f /= (double) nmem_dynamic;

                    /*
                     * A^f = E^f - x^f * 1^T
                     */
                    for (e = 0; e < nmem; ++e)
                        v_f[e] -= v1_f;
                    /*
                     * A^a = A^f * T
                     */
                    {
                        char do_T = 'T';
                        float alpha = 1.0f;
                        float beta = 0.0f;
                        int inc = 1;

                        sgemv_(&do_T, &nmem, &nmem_dynamic, &alpha, Tj[i], &nmem, v_f, &inc, &beta, v_a, &inc);
                    }

                    /*
                     * dx = x^a - x^f
                     */
                    for (e = 0, v1_a = 0.0; e < nmem; ++e)
                        v1_a += v_f[e] * wj[i][e];
                    /*
                     * x^a = x^f + dx
                     */
                    v1_a += v1_f;

                    /*
                     * For hybrid systems -- scale back the ensemble anomalies.
                     *
                     * The anomalies scaled in das_sethybridensemble() are
                     * consistent to correctly factorise the hybrid covariance,
                     * but the analysed ensemble should be conditioned so that
                     * it does not change in a no-obs. case.
                     */
                    if (das->mode == MODE_HYBRID) {
                        double k_d = (nmem_dynamic > 1) ? sqrt((double) (nmem - 1) / (double) (nmem_dynamic - 1)) : 1.0;

                        for (e = 0; e < nmem_dynamic; ++e)
                            v_a[e] /= k_d;
                        for (e = 0; e < nmem_dynamic; ++e)
                            v_f[e] /= k_d;
                    }

                    /*
                     * calculate inflation
                     */
                    if (!isnan(inf_ratio)) {
                        double v2_f = 0.0;
                        double v2_a = 0.0;

                        for (e = 0; e < nmem_dynamic; ++e) {
                            v2_f += v_f[e] * v_f[e];
                            v2_a += v_a[e] * v_a[e];
                        }
                        if (v2_a > 0.0) {
                            /*
                             * (Normal case.) Limit inflation by inf_ratio of
                             * the magnitude of spread reduction.
                             */
                            inflation = (float) (sqrt(v2_f / v2_a) * inf_ratio + 1.0 - inf_ratio);
                            if (inflation >= inflation0)
                                inflation = inflation0;
                        }
                    }

                    /*
                     * (do not inflate if inflation is about 1 or less than 1)
                     */
                    if (inflation - 1.0f > EPSF)
                        for (e = 0; e < nmem_dynamic; ++e)
                            v_a[e] = v_a[e] * inflation;
                    else
                        inflation = 1.0f;

                    /*
                     * E^a = A^a + x^a * 1^T
                     */
                    for (e = 0; e < nmem_dynamic; ++e)
                        v_a[e] += v1_a;

                    /*
                     * E^f = A^f + x^f * 1^T
                     */
                    for (e = 0; e < nmem_dynamic; ++e)
                        v_f[e] += v1_f;

                    if (isstructured) {
                        if (!(das->updatespec & UPDATE_OUTPUTINC))
                            for (e = 0; e < nmem_dynamic; ++e)
                                vvv[e][j][i] = v_a[e];
                        else {
                            if (!applylog)
                                for (e = 0; e < nmem_dynamic; ++e)
                                    vvv[e][j][i] = v_a[e] - v_f[e];
                            else
                                for (e = 0; e < nmem_dynamic; ++e) {
                                    if (!isnormal(v_f[e]))
                                        vvv[e][j][i] = 0.0;
                                    else
                                        /*
                                         * it is necessary to take care that in
                                         * this case (--output-increment) the
                                         * variable is not transformed again
                                         * during writing
                                         */
                                        vvv[e][j][i] = exp10(v_a[e]) - exp10(v_f[e]);
                                }
                        }
                    } else {
                        if (!(das->updatespec & UPDATE_OUTPUTINC))
                            for (e = 0; e < nmem_dynamic; ++e)
                                vv[e][j] = v_a[e];
                        else {
                            if (!applylog)
                                for (e = 0; e < nmem_dynamic; ++e)
                                    vv[e][j] = v_a[e] - v_f[e];
                            else
                                for (e = 0; e < nmem_dynamic; ++e) {
                                    if (!isnormal(v_f[e]))
                                        vv[e][j] = 0.0;
                                    else
                                        /*
                                         * it is necessary to take care that in
                                         * this case (--output-increment) the
                                         * variable is not transformed again
                                         * during writing
                                         */
                                        vv[e][j] = exp10(v_a[e]) - exp10(v_f[e]);
                                }
                        }
                    }

                    if (writeinflation) {
                        if (isstructured)
                            ((float*) infl)[i] = inflation;
                        else
                            ((float**) infl)[fid][j] = inflation;
                    }
                }               /* for i */
                if (writeinflation && isstructured)
                    das_writeinflation(das, f, j, infl);
            }                   /* for fid */
        }                       /* for stepj */
    }                           /* for jj */
    if (writeinflation && !isstructured)
        for (fid = 0; fid < nfield; ++fid)
            das_writeinflation(das, &fields[fid], -1, ((float**) infl)[fid]);

    ncw_close(ncid);

    if (writeinflation)
        free(infl);
    free(v_a);
    free(v_f);
    if (stride > 1) {
        free(Tjj);
        free(Tjj1);
        free(Tjj2);
        free(wjj);
        free(wjj1);
        free(wjj2);
    }
    free(wj);
    free(Tj);

    /*
     * "randomise" ("propagate") fields if required
     */
    for (fid = 0; fid < nfield; ++fid) {
        int varid = fields[fid].varid;

        if (!isnan(model_getvardeflation(m, varid))) {
            /*
             * (for now -- for 2D variables only) 
             */
            addnoise(das, varid, fieldbuffer[fid]);
        }
    }
}

/** Updates `nfield' fields read into `fieldbuffer' with `w'. 
 */
static void das_updatebg(dasystem* das, int nfield, void** fieldbuffer, field fields[])
{
    model* m = das->m;
    int nmem = das->nmem;

    void* g = model_getvargrid(m, fields[0].varid);
    int gridid = grid_getid(g);
    int stride = grid_getstride(g);
    void* nlevels = grid_getnumlevels(g);
    int surfk = grid_getsurflayerid(g);
    int periodic_i = grid_isperiodic_i(g);
    int isstructured = 1;

    char fname[MAXSTRLEN];
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

    grid_getsize(g, &mni, &mnj, NULL);
    /*
     * a treatment for unstructured grids
     */
    if (mnj <= 0) {
        mnj = mni;
        mni = 1;
        isstructured = 0;
    }

    das_getfname_transforms(das, gridid, fname);

    ncw_open(fname, NC_NOWRITE, &ncid);
    ncw_inq_varid(ncid, "w", &varid);
    ncw_inq_vardims(ncid, varid, 3, NULL, dimlens);
    if (rank == 0) {
        int stride_transforms;

        ncw_get_att_int(ncid, NC_GLOBAL, "stride", &stride_transforms);
        if (stride_transforms != stride)
            enkf_quit("grid \"%s\": stride = %d is not equal to the stride of transforms = %d", grid_getname(g), stride, stride_transforms);
    }
    ni = dimlens[1];
    nj = dimlens[0];
    assert((int) dimlens[2] == nmem);

    start[0] = fields[0].j1 / stride;
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
    assert(fields[0].j1 % stride == 0);
    for (jj = fields[0].j1 / stride, j = 0; jj <= fields[0].j2 / stride && jj < nj; ++jj) {
        for (stepj = 0; stepj < stride && j <= fields[0].j2 - fields[0].j1; ++stepj, ++j) {
            if (stride == 1) {
                /*
                 * no interpolation necessary; simply read the ETMs for the
                 * j-th row from disk 
                 */
                start[0] = jj;
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

            /*
             * update the j-th row of the background
             */
            for (fid = 0; fid < nfield; ++fid) {
                field* f = &fields[fid];
                float*** vvv = (isstructured) ? fieldbuffer[fid] : NULL;
                float** vv = (isstructured) ? NULL : fieldbuffer[fid];

                for (i = 0; i < mni; ++i) {
                    float xmean = 0.0f;

                    if (isstructured) {
                        if (surfk == 0) {
                            if (((int**) nlevels)[f->j1 + j][i] <= f->level) {
                                vvv[nmem][j][i] = 0.0f;
                                continue;
                            }
                        } else {
                            if (((int**) nlevels)[f->j1 + j][i] <= ((f->issurfacevar) ? 1 : surfk - f->level)) {
                                vvv[nmem][j][i] = 0.0f;
                                continue;
                            }
                        }
                        /*
                         * assume that if |value| > MAXOBSVAL, then it is filled
                         * with the missing value 
                         */
                        if (fabsf(vvv[nmem][j][i]) > (float) MAXOBSVAL) {
                            vvv[nmem][j][i] = 0.0f;
                            continue;
                        }
                        for (e = 0; e < nmem; ++e)
                            if (fabsf(vvv[e][j][i]) > (float) MAXOBSVAL)
                                break;
                        if (e < nmem) {
                            vvv[nmem][j][i] = 0.0f;
                            continue;
                        }

                        for (e = 0; e < nmem; ++e)
                            xmean += vvv[e][j][i];
                        xmean /= (float) nmem;

                        /*
                         * (the case das->updatespec & UPDATE_OUTPUTINC > 0 is
                         * handled by setting vvv[nmem][][] to zero in
                         * das_update())
                         */
                        for (e = 0; e < nmem; ++e)
                            vvv[nmem][j][i] += (vvv[e][j][i] - xmean) * wj[i][e];
                    } else {
                        if (surfk == 0) {
                            if (((int*) nlevels)[f->j1 + j] <= f->level) {
                                vv[nmem][j] = 0.0f;
                                continue;
                            }
                        } else {
                            if (((int*) nlevels)[f->j1 + j] <= ((f->issurfacevar) ? 1 : surfk - f->level)) {
                                vv[nmem][j] = 0.0f;
                                continue;
                            }
                        }
                        /*
                         * assume that if |value| > MAXOBSVAL, then it is filled
                         * with the missing value 
                         */
                        if (vv[nmem][j] > (float) MAXOBSVAL) {
                            vv[nmem][j] = 0.0f;
                            continue;
                        }
                        for (e = 0; e < nmem; ++e)
                            if (vv[e][j] > (float) MAXOBSVAL)
                                break;
                        if (e < nmem) {
                            vv[nmem][j] = 0.0f;
                            continue;
                        }

                        for (e = 0; e < nmem; ++e)
                            xmean += vv[e][j];
                        xmean /= (float) nmem;

                        /*
                         * (the case das->updatespec & UPDATE_OUTPUTINC > 0 is
                         * handled by setting vvv[nmem][][] to zero in
                         * das_update())
                         */
                        for (e = 0; e < nmem; ++e)
                            vv[nmem][j] += (vv[e][j] - xmean) * wj[0][e];
                    }
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
static void das_writefields_direct(dasystem* das, int nfield, void** fieldbuffer, field fields[])
{
    int i, e;

    for (i = 0; i < nfield; ++i) {
        field* f = &fields[i];

        for (e = 0; e < das->nmem_dynamic; ++e) {
            char fname[MAXSTRLEN];

            das_getmemberfname(das, f->varname, e + 1, fname);
            if (!(das->updatespec & UPDATE_OUTPUTINC)) {
                strncat(fname, ".analysis", MAXSTRLEN - 1);
                model_writefield(das->m, fname, f->varname, f->level, (f->isstructured) ? ((float***) fieldbuffer[i])[e][0] : ((float**) fieldbuffer[i])[e], 0);
            } else {
                strncat(fname, ".increment", MAXSTRLEN - 1);
                model_writefield(das->m, fname, f->varname, f->level, (f->isstructured) ? ((float***) fieldbuffer[i])[e][0] : ((float**) fieldbuffer[i])[e], 1);
            }
        }
    }
}

/** This is the main writing procedure for parallel settings. It writes each
 * horizontal layer (in parallel) to a separate file ("tile"). The tiles are
 * "assembled" later by a single process.
 */
static void das_writefields_toassemble(dasystem* das, int nfield, void** fieldbuffer, field fields[])
{
    char fname[MAXSTRLEN];
    int ni;
    int i;

    model_getvargridsize(das->m, fields[0].varid, &ni, NULL, NULL);

    for (i = 0; i < nfield; ++i) {
        field* f = &fields[i];
        int ncid;
        int dimids[3];
        int vid;

        gettilefname(DIRNAME_TMP, "ens", f, fname);

        ncw_create(fname, NC_CLOBBER | das->ncformat, &ncid);
        ncw_def_dim(ncid, "m", das->nmem_dynamic, &dimids[0]);
        if (f->isstructured) {
            ncw_def_dim(ncid, "j", f->j2 - f->j1 + 1, &dimids[1]);
            ncw_def_dim(ncid, "i", ni, &dimids[2]);
            ncw_def_var(ncid, f->varname, NC_FLOAT, 3, dimids, &vid);
        } else {
            ncw_def_dim(ncid, "i", f->j2 - f->j1 + 1, &dimids[1]);
            ncw_def_var(ncid, f->varname, NC_FLOAT, 2, dimids, &vid);
        }
#if defined(DEFLATE_ALL)
        if (das->nccompression > 0)
            ncw_def_deflate(ncid, 0, 1, das->nccompression);
#endif
        ncw_enddef(ncid);
        ncw_put_var_float(ncid, vid, (f->isstructured) ? ((float***) fieldbuffer[i])[0][0] : ((float**) fieldbuffer[i])[0]);
        ncw_close(ncid);
    }
}

/** Writes `nfield' ensemble fields from `fieldbuffer' to disk.
 */
static void das_writefields(dasystem* das, int nfield, void** fieldbuffer, field fields[])
{
    if (das->updatespec & UPDATE_DIRECTWRITE)
        das_writefields_direct(das, nfield, fieldbuffer, fields);
    else
        das_writefields_toassemble(das, nfield, fieldbuffer, fields);
}

/**
 */
static void das_writebg_direct(dasystem* das, int nfield, void** fieldbuffer, field fields[])
{
    model* m = das->m;
    int ni, nj;
    int i;

    assert(das->mode == MODE_ENOI);

    model_getvargridsize(m, fields[0].varid, &ni, &nj, NULL);

    for (i = 0; i < nfield; ++i) {
        field* f = &fields[i];
        char fname[MAXSTRLEN];

        das_getbgfname(das, f->varname, fname);
        if (!(das->updatespec & UPDATE_OUTPUTINC))
            strncat(fname, ".analysis", MAXSTRLEN - 1);
        else
            strncat(fname, ".increment", MAXSTRLEN - 1);
        model_writefield(m, fname, f->varname, f->level, (f->isstructured) ? ((float***) fieldbuffer[i])[das->nmem][0] : ((float**) fieldbuffer[i])[das->nmem], 0);
    }
}

/**
 */
static void das_writebg_toassemble(dasystem* das, int nfield, void** fieldbuffer, field fields[])
{
    char fname[MAXSTRLEN];
    int ni;
    int i;

    model_getvargridsize(das->m, fields[0].varid, &ni, NULL, NULL);

    for (i = 0; i < nfield; ++i) {
        field* f = &fields[i];
        int ncid;
        int dimids[2];
        int vid;

        gettilefname(DIRNAME_TMP, "bg", f, fname);

        ncw_create(fname, NC_CLOBBER | das->ncformat, &ncid);
        if (f->isstructured) {
            ncw_def_dim(ncid, "j", f->j2 - f->j1 + 1, &dimids[0]);
            ncw_def_dim(ncid, "i", ni, &dimids[1]);
            ncw_def_var(ncid, f->varname, NC_FLOAT, 2, dimids, &vid);
        } else {
            ncw_def_dim(ncid, "i", f->j2 - f->j1 + 1, &dimids[0]);
            ncw_def_var(ncid, f->varname, NC_FLOAT, 1, dimids, &vid);
        }
#if defined(DEFLATE_ALL)
        if (das->nccompression > 0)
            ncw_def_deflate(ncid, 0, 1, das->nccompression);
#endif
        ncw_enddef(ncid);
        ncw_put_var_float(ncid, vid, (f->isstructured) ? ((float***) fieldbuffer[i])[das->nmem][0] : ((float**) fieldbuffer[i])[das->nmem]);
        ncw_close(ncid);
    }
}

/** Writes `nfield' fields from `fieldbuffer' to disk.
 */
static void das_writebg(dasystem* das, int nfield, void** fieldbuffer, field fields[])
{
    if (das->updatespec & UPDATE_DIRECTWRITE)
        das_writebg_direct(das, nfield, fieldbuffer, fields);
    else
        das_writebg_toassemble(das, nfield, fieldbuffer, fields);
}

/**
 */
static void das_preassemblemembers(dasystem* das)
{
    int nmem = das->nmem_dynamic;
    int ntile = 0;
    field* tiles = NULL;
    int nfield = 0;
    int* fids = NULL;
    int tid;

    das_getfields(das, -1, &ntile, &tiles);

    for (tid = 0, nfield = 0; tid < ntile; ++tid) {
        field* f = &tiles[tid];

        if (f->splitid == 0) {
            if (nfield % NTILE_INC == 0)
                fids = realloc(fids, (nfield + NTILE_INC) * sizeof(int));
            fids[nfield] = tid;
            nfield++;
        }
    }
    enkf_printf("    %d tiles, %d fields\n", ntile, nfield);
    enkf_printf("    ");
    enkf_flush();
#if defined(MPI)
    MPI_Barrier(MPI_COMM_WORLD);
#endif

#if defined(USE_MPIQUEUE)
    {
        mpiqueue* queue = NULL;

        queue = mpiqueue_create(MPI_COMM_WORLD, nfield);
        if (mpiqueue_getrank(queue) == 0)
            mpiqueue_manage(queue);
        else {
            while (1) {
                void* vout = NULL;
                char fname_dst[MAXSTRLEN];
                field* f0;
                int fid, splitid;
                int ni, nj;

                fid = mpiqueue_getjob(queue);
                if (fid < 0)
                    break;

                f0 = &tiles[fids[fid]]; /* the first tile of the field */
                model_getvargridsize(das->m, f0->varid, &ni, &nj, NULL);
                if (f0->isstructured)
                    vout = alloc3d(nmem, nj, ni, sizeof(float));
                else
                    vout = alloc2d(nmem, ni, sizeof(float));

                getfieldfname(DIRNAME_TMP, "ens", f0->varname, f0->level, fname_dst);
                if (file_exists(fname_dst))
                    enkf_quit("assembled file \"%s\" already exists", fname_dst);

                for (splitid = 0; splitid < das->nfieldsplit; ++splitid) {
                    field* f = &tiles[fids[fid] + splitid];
                    int njj = f->j2 - f->j1 + 1;
                    char fname_src[MAXSTRLEN];
                    int ncid_src, vid_src;
                    int e;

                    assert(f->splitid == splitid);

                    gettilefname(DIRNAME_TMP, "ens", f, fname_src);
                    ncw_open(fname_src, NC_NOWRITE, &ncid_src);
                    ncw_inq_varid(ncid_src, f->varname, &vid_src);

                    if (f->isstructured) {
                        float*** vin = alloc3d(nmem, njj, ni, sizeof(float));

                        ncw_get_var_float(ncid_src, vid_src, vin[0][0]);
                        for (e = 0; e < nmem; ++e)
                            memcpy(((float***) vout)[e][f->j1], vin[e][0], njj * ni * sizeof(float));
                        free(vin);
                    } else {
                        float** vin = alloc2d(nmem, njj, sizeof(float));

                        ncw_get_var_float(ncid_src, vid_src, vin[0]);
                        for (e = 0; e < nmem; ++e)
                            memcpy(&((float**) vout)[e][f->j1], vin[e], njj * sizeof(float));
                        free(vin);
                    }
                    ncw_close(ncid_src);
                    file_delete(fname_src);

                    printf(".");
                    fflush(stdout);
                }

                {
                    int ncid, vid, dimids[3];

                    ncw_create(fname_dst, NC_CLOBBER | das->ncformat, &ncid);
                    ncw_def_dim(ncid, "m", nmem, &dimids[0]);
                    if (f0->isstructured) {
                        ncw_def_dim(ncid, "j", nj, &dimids[1]);
                        ncw_def_dim(ncid, "i", ni, &dimids[2]);
                        ncw_def_var(ncid, f0->varname, NC_FLOAT, 3, dimids, &vid);
                    } else {
                        ncw_def_dim(ncid, "i", ni, &dimids[1]);
                        ncw_def_var(ncid, f0->varname, NC_FLOAT, 2, dimids, &vid);
                    }
#if defined(DEFLATE_ALL)
                    if (das->nccompression > 0)
                        ncw_def_deflate(ncid, 0, 1, das->nccompression);
#endif
                    ncw_enddef(ncid);
                    ncw_put_var_float(ncid, vid, (f0->isstructured) ? ((float***) vout)[0][0] : ((float**) vout)[0]);
                    ncw_close(ncid);
                }
                free(vout);
                mpiqueue_reportjob(queue, fid);
            }                   /* while (1) */
        }
        mpiqueue_destroy(queue);
#if defined(MPI)
        MPI_Barrier(MPI_COMM_WORLD);
#endif
        enkf_printf("\n");
        enkf_flush();
    }
#else                           /* !defined(USE_MPIQUEUE) */
    {
        int fid;

        distribute_iterations(0, nfield - 1, nprocesses, "    ");

        for (fid = my_first_iteration; fid <= my_last_iteration; ++fid) {
            field* f0 = &tiles[fids[fid]];      /* the first tile of the
                                                 * field */
            void* vout = NULL;
            char fname_dst[MAXSTRLEN];
            int ni, nj;
            int splitid;

            model_getvargridsize(das->m, f0->varid, &ni, &nj, NULL);
            if (f0->isstructured)
                vout = alloc3d(nmem, nj, ni, sizeof(float));
            else
                vout = alloc2d(nmem, ni, sizeof(float));

            getfieldfname(DIRNAME_TMP, "ens", f0->varname, f0->level, fname_dst);
            if (file_exists(fname_dst))
                enkf_quit("assembled file \"%s\" already exists", fname_dst);

            for (splitid = 0; splitid < das->nfieldsplit; ++splitid) {
                field* f = &tiles[fids[fid] + splitid];
                int njj = f->j2 - f->j1 + 1;
                char fname_src[MAXSTRLEN];
                int ncid_src, vid_src;
                int e;

                gettilefname(DIRNAME_TMP, "ens", f, fname_src);
                ncw_open(fname_src, NC_NOWRITE, &ncid_src);
                ncw_inq_varid(ncid_src, f->varname, &vid_src);

                if (f->isstructured) {
                    float*** vin = alloc3d(nmem, njj, ni, sizeof(float));

                    ncw_get_var_float(ncid_src, vid_src, vin[0][0]);
                    for (e = 0; e < nmem; ++e)
                        memcpy(((float***) vout)[e][f->j1], vin[e][0], njj * ni * sizeof(float));
                    free(vin);
                } else {
                    float** vin = alloc2d(nmem, njj, sizeof(float));

                    ncw_get_var_float(ncid_src, vid_src, vin[0]);
                    for (e = 0; e < nmem; ++e)
                        memcpy(&((float**) vout)[e][f->j1], vin[e], njj * sizeof(float));
                    free(vin);
                }
                ncw_close(ncid_src);
                file_delete(fname_src);

                printf(".");
                fflush(stdout);
            }
            {
                int ncid, vid, dimids[3];

                ncw_create(fname_dst, NC_CLOBBER | das->ncformat, &ncid);
                ncw_def_dim(ncid, "m", nmem, &dimids[0]);
                if (f0->isstructured) {
                    ncw_def_dim(ncid, "j", nj, &dimids[1]);
                    ncw_def_dim(ncid, "i", ni, &dimids[2]);
                    ncw_def_var(ncid, f0->varname, NC_FLOAT, 3, dimids, &vid);
                } else {
                    ncw_def_dim(ncid, "i", ni, &dimids[1]);
                    ncw_def_var(ncid, f0->varname, NC_FLOAT, 2, dimids, &vid);
                }
#if defined(DEFLATE_ALL)
                if (das->nccompression > 0)
                    ncw_def_deflate(ncid, 0, 1, das->nccompression);
#endif
                ncw_enddef(ncid);
                ncw_put_var_float(ncid, vid, (f0->isstructured) ? ((float***) vout)[0][0] : ((float**) vout)[0]);
                ncw_close(ncid);
            }
            free(vout);
        }
#if defined(MPI)
        MPI_Barrier(MPI_COMM_WORLD);
#endif
        enkf_printf("\n");
        enkf_flush();
    }
#endif

    if (ntile > 0) {
        free(fids);
        free(tiles);
    }
}

/**
 */
static void das_assemblemembers(dasystem* das)
{
    model* m = das->m;
    int nvar = model_getnvar(m);
    int nf, fid;

    fflush(stdout);
#if defined(MPI)
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    nf = nvar * das->nmem_dynamic;
    distribute_iterations(0, nf - 1, nprocesses, "    ");
    enkf_printf("    ");
    enkf_flush();
    MPI_Barrier(MPI_COMM_WORLD);

    for (fid = my_first_iteration; fid <= my_last_iteration; ++fid) {
        int e = fid % das->nmem_dynamic;
        int vid = fid / das->nmem_dynamic;
        char* varname = model_getvarname(m, vid);
        grid* g = model_getvargrid(m, vid);
        int isstructured = grid_isstructured(g);
        char varname_dst[NC_MAX_NAME];
        char fname_dst[MAXSTRLEN];
        int nlev, k;
        int ni, nj;
        float* v = NULL;

        das_getmemberfname(das, varname, e + 1, fname_dst);
        nlev = ncu_getnlevels(fname_dst, varname, isstructured);
        strncpy(varname_dst, varname, NC_MAX_NAME - 1);

        if (!(das->updatespec & UPDATE_OUTPUTINC))
            strncat(fname_dst, ".analysis", MAXSTRLEN - 1);
        else
            strncat(fname_dst, ".increment", MAXSTRLEN - 1);

        grid_getsize(g, &ni, &nj, NULL);
        if (isstructured)
            v = malloc(ni * nj * sizeof(float));
        else
            v = malloc(ni * sizeof(float));

        for (k = 0; k < nlev; ++k) {
            char fname_src[MAXSTRLEN];
            int ncid_src, vid_src;
            size_t start[3] = { e, 0, 0 };
            size_t count[3] = { 1, nj, ni };

            if (!isstructured)
                count[1] = ni;

            getfieldfname(DIRNAME_TMP, "ens", varname, k, fname_src);
            ncw_open(fname_src, NC_NOWRITE, &ncid_src);
            ncw_inq_varid(ncid_src, varname, &vid_src);
            ncw_get_vara_float(ncid_src, vid_src, start, count, v);
            ncw_close(ncid_src);

            if (!(das->updatespec & UPDATE_OUTPUTINC))
                model_writefield(m, fname_dst, varname_dst, k, v, 0);
            else
                model_writefield(m, fname_dst, varname_dst, k, v, 1);
        }
        free(v);
        enkf_printf(".");
        enkf_flush();
    }
#if defined(MPI)
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    enkf_printf("\n");
    enkf_flush();

    /*
     * remove tiles 
     */
    if (rank == 0) {
        int i;

        for (i = 0; i < nvar; ++i) {
            char* varname = model_getvarname(m, i);
            int isstructured = grid_isstructured(model_getvargrid(m, i));
            char fname[MAXSTRLEN];
            int nlev, k;

            das_getmemberfname(das, varname, 1, fname);
            nlev = ncu_getnlevels(fname, varname, isstructured);
            for (k = 0; k < nlev; ++k) {
                getfieldfname(DIRNAME_TMP, "ens", varname, k, fname);
                file_delete(fname);
            }
        }
    }
}

/**
 */
static void das_preassemblebg(dasystem* das)
{
    int ntile = 0;
    field* tiles = NULL;
    int nfield = 0;
    int* fids = NULL;
    int tid;

    das_getfields(das, -1, &ntile, &tiles);

    for (tid = 0, nfield = 0; tid < ntile; ++tid) {
        field* f = &tiles[tid];

        if (f->splitid == 0) {
            if (nfield % NTILE_INC == 0)
                fids = realloc(fids, (nfield + NTILE_INC) * sizeof(int));
            fids[nfield] = tid;
            nfield++;
        }
    }
    enkf_printf("    %d tiles, %d fields\n", ntile, nfield);
    enkf_printf("    ");
    enkf_flush();
    MPI_Barrier(MPI_COMM_WORLD);

#if defined(USE_MPIQUEUE)
    {
        mpiqueue* queue = NULL;

        queue = mpiqueue_create(MPI_COMM_WORLD, nfield);
        if (mpiqueue_getrank(queue) == 0)
            mpiqueue_manage(queue);
        else {
            while (1) {
                void* vout = NULL;
                char fname_dst[MAXSTRLEN];
                field* f0;
                int fid, splitid;
                int ni, nj;

                fid = mpiqueue_getjob(queue);
                if (fid < 0)
                    break;

                f0 = &tiles[fids[fid]]; /* the first tile of the field */
                model_getvargridsize(das->m, f0->varid, &ni, &nj, NULL);
                if (f0->isstructured)
                    vout = alloc2d(nj, ni, sizeof(float));
                else
                    vout = calloc(ni, sizeof(float));

                getfieldfname(DIRNAME_TMP, "bg", f0->varname, f0->level, fname_dst);
                if (file_exists(fname_dst))
                    enkf_quit("assembled file \"%s\" already exists", fname_dst);

                for (splitid = 0; splitid < das->nfieldsplit; ++splitid) {
                    field* f = &tiles[fids[fid] + splitid];
                    int njj = f->j2 - f->j1 + 1;
                    char fname_src[MAXSTRLEN];
                    int ncid_src, vid_src;

                    gettilefname(DIRNAME_TMP, "bg", f, fname_src);
                    ncw_open(fname_src, NC_NOWRITE, &ncid_src);
                    ncw_inq_varid(ncid_src, f->varname, &vid_src);

                    if (f->isstructured) {
                        float** vin = alloc2d(njj, ni, sizeof(float));

                        ncw_get_var_float(ncid_src, vid_src, vin[0]);
                        memcpy(((float**) vout)[f->j1], vin[0], njj * ni * sizeof(float));
                        free(vin);
                    } else {
                        float* vin = calloc(njj, sizeof(float));

                        ncw_get_var_float(ncid_src, vid_src, vin);
                        memcpy(&((float*) vout)[f->j1], vin, njj * sizeof(float));
                        free(vin);
                    }
                    ncw_close(ncid_src);
                    file_delete(fname_src);

                    printf(".");
                    fflush(stdout);
                }

                {
                    int ncid, vid, dimids[2];

                    ncw_create(fname_dst, NC_CLOBBER | das->ncformat, &ncid);
                    if (f0->isstructured) {
                        ncw_def_dim(ncid, "j", nj, &dimids[0]);
                        ncw_def_dim(ncid, "i", ni, &dimids[1]);
                        ncw_def_var(ncid, f0->varname, NC_FLOAT, 2, dimids, &vid);
                    } else {
                        ncw_def_dim(ncid, "i", ni, &dimids[0]);
                        ncw_def_var(ncid, f0->varname, NC_FLOAT, 1, dimids, &vid);
                    }
#if defined(DEFLATE_ALL)
                    if (das->nccompression > 0)
                        ncw_def_deflate(ncid, 0, 1, das->nccompression);
#endif
                    ncw_enddef(ncid);
                    ncw_put_var_float(ncid, vid, (f0->isstructured) ? ((float**) vout)[0] : vout);
                    ncw_close(ncid);
                }
                free(vout);
                mpiqueue_reportjob(queue, fid);
            }                   /* while (1) */
        }
        mpiqueue_destroy(queue);

#if defined(MPI)
        MPI_Barrier(MPI_COMM_WORLD);
#endif
        enkf_printf("\n");
        enkf_flush();
    }
#else                           /* !defined(USE_MPUQUEUE) */
    {
        int fid;

        distribute_iterations(0, nfield - 1, nprocesses, "    ");

        for (fid = my_first_iteration; fid <= my_last_iteration; ++fid) {
            field* f0 = &tiles[fids[fid]];      /* the first tile of the
                                                 * field */
            void* vout = NULL;
            char fname_dst[MAXSTRLEN];
            int ni, nj;
            int splitid;

            model_getvargridsize(das->m, f0->varid, &ni, &nj, NULL);
            if (f0->isstructured)
                vout = alloc2d(nj, ni, sizeof(float));
            else
                vout = calloc(ni, sizeof(float));

            getfieldfname(DIRNAME_TMP, "bg", f0->varname, f0->level, fname_dst);
            if (file_exists(fname_dst))
                enkf_quit("assembled file \"%s\" already exists", fname_dst);

            for (splitid = 0; splitid < das->nfieldsplit; ++splitid) {
                field* f = &tiles[fids[fid] + splitid];
                int njj = f->j2 - f->j1 + 1;
                char fname_src[MAXSTRLEN];
                int ncid_src, vid_src;

                gettilefname(DIRNAME_TMP, "bg", f, fname_src);
                ncw_open(fname_src, NC_NOWRITE, &ncid_src);
                ncw_inq_varid(ncid_src, f->varname, &vid_src);

                if (f->isstructured) {
                    float** vin = alloc2d(njj, ni, sizeof(float));

                    ncw_get_var_float(ncid_src, vid_src, vin[0]);
                    memcpy(((float**) vout)[f->j1], vin[0], njj * ni * sizeof(float));
                    free(vin);
                } else {
                    float* vin = calloc(njj, sizeof(float));

                    ncw_get_var_float(ncid_src, vid_src, vin);
                    memcpy(&((float*) vout)[f->j1], vin, njj * sizeof(float));
                    free(vin);
                }
                ncw_close(ncid_src);
                file_delete(fname_src);

                enkf_printf(".");
                enkf_flush();
            }

            {
                int ncid, vid, dimids[2];

                ncw_create(fname_dst, NC_CLOBBER | das->ncformat, &ncid);
                if (f0->isstructured) {
                    ncw_def_dim(ncid, "j", nj, &dimids[0]);
                    ncw_def_dim(ncid, "i", ni, &dimids[1]);
                    ncw_def_var(ncid, f0->varname, NC_FLOAT, 2, dimids, &vid);
                } else {
                    ncw_def_dim(ncid, "i", ni, &dimids[0]);
                    ncw_def_var(ncid, f0->varname, NC_FLOAT, 1, dimids, &vid);
                }
#if defined(DEFLATE_ALL)
                if (das->nccompression > 0)
                    ncw_def_deflate(ncid, 0, 1, das->nccompression);
#endif
                ncw_enddef(ncid);
                ncw_put_var_float(ncid, vid, (f0->isstructured) ? ((float**) vout)[0] : vout);
                ncw_close(ncid);
            }
            free(vout);
        }
#if defined(MPI)
        MPI_Barrier(MPI_COMM_WORLD);
#endif
        enkf_printf("\n");
        enkf_flush();
    }
#endif

    if (ntile > 0) {
        free(fids);
        free(tiles);
    }
}

/**
 */
static void das_assemblebg(dasystem* das)
{
    model* m = das->m;
    int nvar = model_getnvar(m);
    int i;

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
        strncpy(varname_dst, varname, NC_MAX_NAME - 1);

        das_getbgfname(das, varname, fname_dst);
        if (!(das->updatespec & UPDATE_OUTPUTINC))
            strncat(fname_dst, ".analysis", MAXSTRLEN - 1);
        else
            strncat(fname_dst, ".increment", MAXSTRLEN - 1);

        model_getvargridsize(m, i, &ni, &nj, NULL);
        if (nj > 0)
            v = malloc(ni * nj * sizeof(float));
        else
            v = malloc(ni * sizeof(float));

        nlev = ncu_getnlevels(fname_dst, varname, nj > 0);
        for (k = 0; k < nlev; ++k) {
            char fname_src[MAXSTRLEN];
            int ncid_src, vid_src;

            getfieldfname(DIRNAME_TMP, "bg", varname, k, fname_src);
            ncw_open(fname_src, NC_NOWRITE, &ncid_src);
            ncw_inq_varid(ncid_src, varname, &vid_src);
            ncw_get_var_float(ncid_src, vid_src, v);
            ncw_close(ncid_src);

            model_writefield(m, fname_dst, varname_dst, k, v, 0);
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
    int nvar = model_getnvar(m);

#if !defined(USE_MPIQUEUE)
    int ngrid = model_getngrid(m);
    int gid;
#endif

#if defined(MPI)
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    fflush(NULL);

    if (das->updatespec & UPDATE_DOINFLATION && rank == 0) {
        enkf_printf("    allocating disk space for inflation:\n");
        enkf_printtime("      ");
        das_allocatedst(das, ROOTNAME_INFLATION);
        enkf_flush();
        if (rank == 0 && !(das->updatespec & UPDATE_DIRECTWRITE))
            dir_createifabsent(DIRNAME_TMP);
    }

    if (das->updatespec & UPDATE_DOFIELDS) {
        int firstprint = 1;

        enkf_printf("    allocating disk space for analysis:\n");
        enkf_printtime("      ");
        if (das->mode == MODE_ENKF || das->mode == MODE_HYBRID) {
            int i, e;

            distribute_iterations(0, das->nmem_dynamic - 1, nprocesses, "      ");
#if defined(MPI)
            MPI_Barrier(MPI_COMM_WORLD);
#endif
            for (i = 0; i < nvar; ++i) {
                char* varname = model_getvarname(m, i);

                for (e = my_first_iteration; e <= my_last_iteration; ++e) {
                    char fname_f[MAXSTRLEN], fname_a[MAXSTRLEN];
                    int ncid_f, ncid_a;
                    int vid_f;

                    das_getmemberfname(das, varname, e + 1, fname_f);
                    strncpy(fname_a, fname_f, MAXSTRLEN);
                    if (!(das->updatespec & UPDATE_OUTPUTINC))
                        strncat(fname_a, ".analysis", MAXSTRLEN - 1);
                    else
                        strncat(fname_a, ".increment", MAXSTRLEN - 1);

                    if (file_exists(fname_a)) {
                        ncw_open(fname_a, NC_WRITE, &ncid_a);
                        if (ncw_var_exists(ncid_a, varname)) {
                            ncw_close(ncid_a);
                            continue;
                        }
                        ncw_redef(ncid_a);
                    } else {
                        ncw_create(fname_a, NC_CLOBBER | das->ncformat, &ncid_a);
                        das->haveanalysis = 0;
                    }

                    ncw_open(fname_f, NC_NOWRITE, &ncid_f);
                    ncw_inq_varid(ncid_f, varname, &vid_f);
                    ncw_copy_vardef(ncid_f, vid_f, ncid_a);
                    if (das->nccompression > 0)
                        ncw_def_deflate(ncid_a, 0, 1, das->nccompression);
                    ncw_close(ncid_a);
                    ncw_close(ncid_f);
                    if (firstprint) {
                        enkf_printf("      ");
                        firstprint = 0;
                    }
                    enkf_printf(".");
                    enkf_flush();
                }
            }
            if (!firstprint) {
                enkf_printf("\n");
                enkf_flush();
            }
        } else if (das->mode == MODE_ENOI) {
            int firstprint = 1;
            int i;

            if (rank == 0) {
                for (i = 0; i < nvar; ++i) {
                    char* varname = model_getvarname(m, i);
                    char fname_f[MAXSTRLEN], fname_a[MAXSTRLEN];
                    int ncid_f, ncid_a;
                    int vid_f;

                    das_getbgfname(das, varname, fname_f);
                    strncpy(fname_a, fname_f, MAXSTRLEN);
                    if ((das->updatespec & UPDATE_OUTPUTINC))
                        strncat(fname_a, ".increment", MAXSTRLEN - 1);
                    else
                        strncat(fname_a, ".analysis", MAXSTRLEN - 1);

                    if (file_exists(fname_a)) {
                        ncw_open(fname_a, NC_WRITE, &ncid_a);
                        if (ncw_var_exists(ncid_a, varname)) {
                            ncw_close(ncid_a);
                            continue;
                        }
                        ncw_redef(ncid_a);
                    } else {
                        ncw_create(fname_a, NC_CLOBBER | das->ncformat, &ncid_a);
                        das->haveanalysis = 0;
                    }

                    if (!file_exists(fname_f) && (das->updatespec & UPDATE_OUTPUTINC))
                        das_getmemberfname(das, varname, 1, fname_f);
                    ncw_open(fname_f, NC_NOWRITE, &ncid_f);
                    ncw_inq_varid(ncid_f, varname, &vid_f);
                    ncw_copy_vardef(ncid_f, vid_f, ncid_a);
                    if (das->nccompression > 0)
                        ncw_def_deflate(ncid_a, 0, 1, das->nccompression);
                    ncw_close(ncid_a);
                    ncw_close(ncid_f);
                    if (firstprint) {
                        enkf_printf("      ");
                        firstprint = 0;
                    }
                    enkf_printf(".");
                    enkf_flush();
                }
            }
#if defined(MPI)
            MPI_Barrier(MPI_COMM_WORLD);
#endif
            if (!firstprint) {
                enkf_printf("\n");
                enkf_flush();
            }
        }
    }
#if defined(MPI)
    MPI_Barrier(MPI_COMM_WORLD);
#endif
#if defined(USE_MPIQUEUE)
    /*
     * Using mpiqueue is a preferred way of updating the fields because it
     * processes all fields (on all grids) on one-by-one basis and therefore
     * reduces idle periods that occur otherwise.
     */
    if (das->updatespec & UPDATE_NEEDAN || das->updatespec & UPDATE_DOSPREAD) {
        int nfield = 0;
        field* fields = NULL;
        mpiqueue* queue = NULL;

        das_getfields(das, -1, &nfield, &fields);
        if (nfield == 0)
            enkf_quit("nothing to do: nfield = 0");
        if (nprocesses == 1)
            enkf_quit("\"mpiqueue\" can not be used on a single CPU; run on more than one CPU or recompile without -DUSE_MPIQUEUE flag");
        queue = mpiqueue_create(MPI_COMM_WORLD, nfield);
        if (das->updatespec & UPDATE_NEEDAN)
            enkf_printf("    updating %d %s using %d processes:\n", nfield, (das->nfieldsplit == 1) ? "fields" : "tiles", nprocesses);
        else
            enkf_printf("    reading %d %s using %d processes:\n", nfield, (das->nfieldsplit == 1) ? "fields" : "tiles", nprocesses);
        enkf_printtime("      ");
        enkf_flush();
        MPI_Barrier(MPI_COMM_WORLD);

        if (mpiqueue_getrank(queue) == 0)
            mpiqueue_manage(queue);
        else {
            void** fieldbuffer = calloc(das->fieldbufsize, sizeof(void*));
            field* fieldstowrite = malloc(das->fieldbufsize * sizeof(field));
            int bufid = 0;
            int gridid_storage = -1;
            int gridid_buf;
            int nij = -1;
            int i;

            while (1) {
                field* f;
                int fid, gridid_field, e;

                fid = mpiqueue_getjob(queue);
                if (fid < 0) {
                    if (bufid > 0) {
                        bufid--;
                        goto doupdate;
                    } else
                        break;
                }

                f = &fields[fid];
                gridid_field = model_getvargridid(das->m, f->varid);

                if (bufid == 0) {
                    gridid_buf = gridid_field;
                } else if (gridid_field != gridid_buf) {
                    mpiqueue_rejectjob(queue, fid);
                    bufid--;
                    goto doupdate;
                }
                /*
                 * if new grid -- allocate storage for fields
                 */
                if (gridid_buf != gridid_storage || das->nfieldsplit > 1) {
                    grid* g = model_getgridbyid(m, gridid_field);
                    int mni, mnj, isize, jsize;

                    grid_getsize(g, &mni, &mnj, NULL);
                    for (i = 0; i < das->fieldbufsize; ++i)
                        if (fieldbuffer[i] != NULL)
                            free(fieldbuffer[i]);
                    if (mnj > 0) {
                        isize = mni;
                        jsize = f->j2 - f->j1 + 1;
                        nij = isize * jsize;

                        for (i = 0; i < das->fieldbufsize; ++i) {
                            if (das->mode == MODE_ENKF)
                                fieldbuffer[i] = alloc3d(das->nmem, jsize, isize, sizeof(float));
                            else if (das->mode == MODE_ENOI)
                                fieldbuffer[i] = alloc3d(das->nmem + 1, jsize, isize, sizeof(float));
                            else if (das->mode == MODE_HYBRID)
                                /*
                                 * allocate two additional members to calculate
                                 * ensemble mean with double precision
                                 */
                                fieldbuffer[i] = alloc3d(das->nmem + 2, jsize, isize, sizeof(float));
                        }
                    } else {
                        isize = f->j2 - f->j1 + 1;
                        nij = isize;

                        for (i = 0; i < das->fieldbufsize; ++i) {
                            if (das->mode == MODE_ENKF)
                                fieldbuffer[i] = alloc2d(das->nmem, isize, sizeof(float));
                            else if (das->mode == MODE_ENOI)
                                fieldbuffer[i] = alloc2d(das->nmem + 1, isize, sizeof(float));
                            else if (das->mode == MODE_HYBRID)
                                /*
                                 * allocate two additional members to calculate
                                 * ensemble mean with double precision
                                 */
                                fieldbuffer[i] = alloc2d(das->nmem + 2, isize, sizeof(float));
                        }
                    }
                    gridid_storage = gridid_buf;
                }
                /*
                 * read ensemble of fields
                 */
                if (enkf_verbose) {
                    if (das->nfieldsplit == 1)
                        printf("      %-6s %-3d (%d: %d: %.1f%%)\n", f->varname, f->level, rank, fid, 100.0 * (double) (fid + 1) / (double) nfield);
                    else
                        printf("      %-6s %3d-%d (%d: %d: %.1f%%)\n", f->varname, f->level, f->splitid, rank, fid, 100.0 * (double) (fid + 1) / (double) nfield);
                    fflush(stdout);
                }
                for (e = 0; e < das->nmem; ++e) {
                    char fname[MAXSTRLEN];
                    int masklog = das_isstatic(das, e + 1);

                    das_getmemberfname(das, f->varname, e + 1, fname);
                    model_readfield_part(das->m, fname, f->varname, -1, f->level, f->j1, f->j2, (f->isstructured) ? ((float***) fieldbuffer[bufid])[e][0] : ((float**) fieldbuffer[bufid])[e], masklog);
                }
                if (das->mode == MODE_HYBRID) {
                    float** v = calloc(das->nmem + 1, sizeof(float*));

                    for (e = 0; e < das->nmem + 1; ++e)
                        v[e] = (f->isstructured) ? ((float***) fieldbuffer[bufid])[e][0] : ((float**) fieldbuffer[bufid])[e];

                    das_sethybridensemble(das, nij, v);
                    free(v);
                }
                /*
                 * read the background to write it to pointlogs, regardless
                 * of whether output is increment or analysis
                 */
                if (das->mode == MODE_ENOI && (((das->updatespec & UPDATE_DOFIELDS) && !(das->updatespec & UPDATE_OUTPUTINC)) || (das->updatespec & UPDATE_DOPLOGSAN))) {
                    char fname[MAXSTRLEN];

                    das_getbgfname(das, f->varname, fname);
                    model_readfield_part(das->m, fname, f->varname, -1, f->level, f->j1, f->j2, (f->isstructured) ? ((float***) fieldbuffer[bufid])[das->nmem][0] : ((float**) fieldbuffer[bufid])[das->nmem], 0);
                }

                fieldstowrite[bufid] = fields[fid];
                mpiqueue_reportjob(queue, fid);

                if (bufid < das->fieldbufsize - 1) {
                    bufid++;
                    continue;
                }

              doupdate:
                /*
                 * write forecast spread
                 */
                if (das->updatespec & UPDATE_DOSPREAD)
                    das_writespread(das, bufid + 1, fieldbuffer, fieldstowrite);
                /*
                 * set the background to 0 if output is increment
                 */
                if (das->mode == MODE_ENOI && (das->updatespec & UPDATE_OUTPUTINC)) {
                    int mni, mnj;
                    int ii;

                    model_getvargridsize(das->m, f->varid, &mni, &mnj, NULL);
                    for (ii = 0; ii <= bufid; ++ii)
                        memset((f->isstructured) ? ((float***) fieldbuffer[ii])[das->nmem][0] : ((float**) fieldbuffer[bufid])[das->nmem], 0, nij * sizeof(float));
                }
                /*
                 * update
                 */
                if (das->updatespec & UPDATE_NEEDAN) {
                    if (das->mode == MODE_ENKF || das->mode == MODE_HYBRID) {
                        das_updatefields(das, bufid + 1, fieldbuffer, fieldstowrite);
                        if (das->updatespec & UPDATE_DOFIELDS)
                            das_writefields(das, bufid + 1, fieldbuffer, fieldstowrite);
                    } else if (das->mode == MODE_ENOI) {
                        if (das->updatespec & (UPDATE_DOFIELDS | UPDATE_DOPLOGSAN))
                            das_updatebg(das, bufid + 1, fieldbuffer, fieldstowrite);
                        if (das->updatespec & UPDATE_DOFIELDS)
                            das_writebg(das, bufid + 1, fieldbuffer, fieldstowrite);
                    }
                }

                if (fid < 0)
                    break;
                bufid = 0;
            }                   /* while (1) */

            for (i = 0; i < das->fieldbufsize; ++i)
                if (fieldbuffer[i] != NULL)
                    free(fieldbuffer[i]);
            free(fieldstowrite);
            free(fieldbuffer);
        }                       /* if (worker) */
        mpiqueue_destroy(queue);
        free(fields);
    }
#else                           /* !defined(USE_MPIQUEUE) */
    if (das->updatespec & UPDATE_NEEDAN || das->updatespec & UPDATE_DOSPREAD) {
        for (gid = 0; gid < ngrid; ++gid) {
            void* g = model_getgridbyid(m, gid);
            int nfield = 0;
            field* fields = NULL;
            void** fieldbuffer = NULL;
            int nij = -1;
            int mni, mnj;
            int i, fid;

            /*
             * (just to avoid the log message below)
             */
            if (grid_getaliasid(g) >= 0)
                continue;

#if defined(MPI)
            MPI_Barrier(MPI_COMM_WORLD);
#endif
            enkf_printf("    processing fields for %s:\n", grid_getname(g));
            enkf_printtime("      ");

            das_getfields(das, gid, &nfield, &fields);
            if (nfield == 0)
                continue;

            if (das->updatespec & UPDATE_NEEDAN)
                enkf_printf("    updating %d %s using %d processes:\n", nfield, (das->nfieldsplit == 1) ? "fields" : "tiles", nprocesses);
            else
                enkf_printf("    reading %d %s using %d processes:\n", nfield, (das->nfieldsplit == 1) ? "fields" : "tiles", nprocesses);
            enkf_flush();
#if defined(MPI)
            MPI_Barrier(MPI_COMM_WORLD);
#endif

            grid_getsize(g, &mni, &mnj, NULL);

            distribute_iterations(0, nfield - 1, nprocesses, "      ");
            if (my_first_iteration > my_last_iteration)
                continue;

            fieldbuffer = calloc(das->fieldbufsize, sizeof(void*));
            for (fid = my_first_iteration; fid <= my_last_iteration; ++fid) {
                int bufid = (fid - my_first_iteration) % das->fieldbufsize;
                field* f = &fields[fid];
                char fname[MAXSTRLEN];
                int e;

                if (enkf_verbose) {
                    if (das->nfieldsplit == 1)
                        printf("      %-6s %-3d (%d: %d: %.1f%%)\n", f->varname, f->level, rank, fid, 100.0 * (double) (fid + 1) / (double) nfield);
                    else
                        printf("      %-6s %3d-%d (%d: %d: %.1f%%)\n", f->varname, f->level, f->splitid, rank, fid, 100.0 * (double) (fid + 1) / (double) nfield);
                    fflush(stdout);
                }

                if (fieldbuffer[bufid] == NULL || das->nfieldsplit > 1) {
                    if (mnj > 0) {
                        int jsize = f->j2 - f->j1 + 1;

                        nij = mni * jsize;
                        if (das->mode == MODE_ENKF) {
                            fieldbuffer[bufid] = alloc3d(das->nmem, jsize, mni, sizeof(float));
                        } else if (das->mode == MODE_ENOI) {
                            fieldbuffer[bufid] = alloc3d(das->nmem + 1, jsize, mni, sizeof(float));
                        } else if (das->mode == MODE_HYBRID) {
                            /*
                             * allocate two additional members to calculate ensemble
                             * mean with double precision
                             */
                            fieldbuffer[bufid] = alloc3d(das->nmem + 2, jsize, mni, sizeof(float));
                        }
                    } else {
                        int isize = f->j2 - f->j1 + 1;

                        nij = isize;
                        if (das->mode == MODE_ENKF) {
                            fieldbuffer[bufid] = alloc2d(das->nmem, isize, sizeof(float));
                        } else if (das->mode == MODE_ENOI) {
                            fieldbuffer[bufid] = alloc2d(das->nmem + 1, isize, sizeof(float));
                        } else if (das->mode == MODE_HYBRID) {
                            /*
                             * allocate two additional members to calculate ensemble
                             * mean with double precision
                             */
                            fieldbuffer[bufid] = alloc2d(das->nmem + 2, isize, sizeof(float));
                        }
                    }
                }
                for (e = 0; e < das->nmem; ++e) {
                    int masklog = das_isstatic(das, e + 1);

                    das_getmemberfname(das, f->varname, e + 1, fname);
                    model_readfield_part(das->m, fname, f->varname, -1, f->level, f->j1, f->j2, (f->isstructured) ? ((float***) fieldbuffer[bufid])[e][0] : ((float**) fieldbuffer[bufid])[e], masklog);
                }

                if (das->mode == MODE_HYBRID) {
                    float* v[das->nmem + 1];

                    for (e = 0; e < das->nmem + 1; ++e)
                        v[e] = (f->isstructured) ? ((float***) fieldbuffer[bufid])[e][0] : ((float**) fieldbuffer[bufid])[e];

                    das_sethybridensemble(das, nij, v);
                }

                /*
                 * read the background to write it to pointlogs, regardless of
                 * whether output is increment or analysis
                 */
                if (das->mode == MODE_ENOI && (((das->updatespec & UPDATE_DOFIELDS) && !(das->updatespec & UPDATE_OUTPUTINC)) || (das->updatespec & UPDATE_DOPLOGSAN))) {
                    das_getbgfname(das, f->varname, fname);
                    model_readfield_part(das->m, fname, f->varname, -1, f->level, f->j1, f->j2, (f->isstructured) ? ((float***) fieldbuffer[bufid])[das->nmem][0] : ((float**) fieldbuffer[bufid])[das->nmem], 0);
                }

                if (bufid == das->fieldbufsize - 1 || fid == my_last_iteration) {
                    /*
                     * write forecast spread
                     */
                    if (das->updatespec & UPDATE_DOSPREAD)
                        das_writespread(das, bufid + 1, fieldbuffer, &fields[fid - bufid]);

                    /*
                     * set the background to 0 if output is increment
                     */
                    if (das->mode == MODE_ENOI && (das->updatespec & UPDATE_OUTPUTINC)) {
                        int ii;

                        for (ii = 0; ii <= bufid; ++ii)
                            memset((f->isstructured) ? ((float***) fieldbuffer[ii])[das->nmem][0] : ((float**) fieldbuffer[bufid])[das->nmem], 0, nij * sizeof(float));
                    }

                    if (das->updatespec & UPDATE_NEEDAN) {
                        if (das->mode == MODE_ENKF || das->mode == MODE_HYBRID) {
                            das_updatefields(das, bufid + 1, fieldbuffer, &fields[fid - bufid]);
                            if (das->updatespec & UPDATE_DOFIELDS)
                                das_writefields(das, bufid + 1, fieldbuffer, &fields[fid - bufid]);
                            else if (fid == my_last_iteration)
                                enkf_printf("      (skip writing the fields)\n");
                        } else if (das->mode == MODE_ENOI) {
                            if (das->updatespec & (UPDATE_DOFIELDS | UPDATE_DOPLOGSAN))
                                das_updatebg(das, bufid + 1, fieldbuffer, &fields[fid - bufid]);
                            if (das->updatespec & UPDATE_DOFIELDS)
                                das_writebg(das, bufid + 1, fieldbuffer, &fields[fid - bufid]);
                            else if (fid == my_last_iteration)
                                enkf_printf("      (skip writing the fields)\n");
                        }
                    }
                    for (i = 0; i <= bufid; ++i) {
                        free(fieldbuffer[i]);
                        fieldbuffer[i] = NULL;
                    }
                }
            }                   /* for fid */

            free(fieldbuffer);
            free(fields);

            enkf_flush();
        }                       /* for gid */
    }
#endif                          /* USE_MPIQUEUE */

#if defined(MPI)
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    if (!(das->updatespec & UPDATE_DIRECTWRITE)) {
        /*
         * assemble restarts
         */
        if (das->updatespec & UPDATE_DOFIELDS) {
            if (das->nfieldsplit > 1) {
                if (das->mode == MODE_ENKF || das->mode == MODE_HYBRID) {
                    enkf_printf("  assembling member fields:\n");
                    enkf_printtime("    ");
                    das_preassemblemembers(das);
                } else if (das->mode == MODE_ENOI) {
                    enkf_printf("  assembling bg fields:\n");
                    enkf_printtime("    ");
                    das_preassemblebg(das);
                }
            }
            if (das->mode == MODE_ENKF || das->mode == MODE_HYBRID) {
                enkf_printf("  assembling members:\n");
                enkf_printtime("    ");
                das_assemblemembers(das);
            } else if (das->mode == MODE_ENOI) {
                enkf_printf("  assembling bg:\n");
                enkf_printtime("    ");
                das_assemblebg(das);
            }
        }
        /*
         * assemble spread
         */
        if (das->updatespec & UPDATE_DOSPREAD) {
            if (das->nfieldsplit > 1) {
                enkf_printf("  assembling spread fields:\n");
                enkf_printtime("    ");
                das_preassemble(das, ROOTNAME_SPREAD);
#if defined(MPI)
                MPI_Barrier(MPI_COMM_WORLD);
#endif
            }
            if (rank == 0) {
                enkf_printf("  assembling spread:\n");
                enkf_printtime("    ");
                das_assemble(das, ROOTNAME_SPREAD);
            }
        }
        /*
         * assemble inflation
         */
        if (das->updatespec & UPDATE_DOINFLATION) {
            if (das->nfieldsplit > 1) {
                enkf_printf("  assembling inflation fields:\n");
                enkf_printtime("    ");
                das_preassemble(das, ROOTNAME_INFLATION);
#if defined(MPI)
                MPI_Barrier(MPI_COMM_WORLD);
#endif
            }
            if (rank == 0) {
                enkf_printf("  assembling inflation:\n");
                enkf_printtime("    ");
                das_assemble(das, ROOTNAME_INFLATION);
            }
        }
    }
    if ((das->updatespec & UPDATE_DOPLOGSAN) || das->haveanalysis) {
        if (das->updatespec & UPDATE_DOPLOGSAN)
            enkf_printf("  writing analysed variables to point logs:\n");
        else
            enkf_printf("  writing analysed variables from existing analysis to point logs:\n");
        enkf_printtime("  ");
        plogs_writestatevars(das, 1);
    }
#if defined(MPI)
    MPI_Barrier(MPI_COMM_WORLD);
#endif
}
