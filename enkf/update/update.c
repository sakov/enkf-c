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

/** Updates `nfield' fields read into `fieldbuffer' with `X5'. Applies
 * variable-dependent inflation to each field.
 */
static void das_updatefields(dasystem* das, int nfields, void** fieldbuffer, field fields[])
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
    int structured = 1;

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
        structured = 0;
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

    start[0] = 0;
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
        infl = (structured) ? malloc(mni * sizeof(float)) : alloc2d(nfields, mnj, sizeof(float));

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
            for (fid = 0; fid < nfields; ++fid) {
                field* f = &fields[fid];
                int applylog = model_getvarislog(m, f->varid);
                float*** vvv = (float***) fieldbuffer[fid];
                char do_T = 'T';
                float alpha = 1.0f;
                int inc = 1;
                float beta = 0.0f;
                float inflation0 = NAN;
                double inf_ratio = NAN;

                model_getvarinflation(m, f->varid, &inflation0, &inf_ratio);
                if (writeinflation && structured)
                    memset(infl, 0, mni * sizeof(float));

                for (i = 0; i < mni; ++i) {
                    float inflation = inflation0;
                    double v1_f, v1_a;

                    /*
                     * For now we assume that layers are counted down from the
                     * surface. (This is not so in ROMS, but there the number
                     * of layers is always the same.) This will be easy to
                     * modify as soon as we encounter a Z model with layers
                     * counted up from the bottom.
                     */
                    if (structured == 1) {
                        if (surfk == 0) {
                            if (((int**) nlevels)[j][i] <= f->level) {
                                for (e = 0; e < nmem; ++e)
                                    ((float***) vvv)[e][j][i] = 0.0f;
                                continue;
                            }
                        } else {
                            if (((int**) nlevels)[j][i] <= ((f->issurfacevar) ? 0 : surfk - f->level)) {
                                for (e = 0; e < nmem; ++e)
                                    ((float***) vvv)[e][j][i] = 0.0f;
                                continue;
                            }
                        }
                    } else {
                        if (surfk == 0) {
                            if (((int*) nlevels)[j] <= f->level) {
                                for (e = 0; e < nmem; ++e)
                                    ((float**) vvv)[e][j] = 0.0f;
                                continue;
                            }
                        } else {
                            if (((int*) nlevels)[j] <= ((f->issurfacevar) ? 0 : surfk - f->level)) {
                                for (e = 0; e < nmem; ++e)
                                    ((float**) vvv)[e][j] = 0.0f;
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
                    if (structured == 1) {
                        for (e = 0; e < nmem; ++e)
                            if (!isfinite(((float***) vvv)[e][j][i]) || fabsf(((float***) vvv)[e][j][i]) > (float) MAXOBSVAL)
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
                            v_f[e] = ((float***) vvv)[e][j][i];
                    } else {
                        for (e = 0; e < nmem; ++e)
                            if (!isfinite(((float**) vvv)[e][j]) || fabsf(((float**) vvv)[e][j]) > (float) MAXOBSVAL)
                                break;
                        if (e < nmem_dynamic) {
                            for (e = 0; e < nmem; ++e)
                                ((float**) vvv)[e][j] = 0.0f;
                            continue;
                        } else if (e < nmem) {
                            for (e = nmem_dynamic; e < nmem; ++e)
                                ((float**) vvv)[e][j] = 0.0f;
                            if (nmem_dynamic == 0)
                                continue;
                        }

                        for (e = 0; e < nmem; ++e)
                            v_f[e] = ((float**) vvv)[e][j];
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
                    sgemv_(&do_T, &nmem, &nmem_dynamic, &alpha, Tj[i], &nmem, v_f, &inc, &beta, v_a, &inc);

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

                    if (structured == 1) {
                        if (!(das->updatespec & UPDATE_OUTPUTINC))
                            for (e = 0; e < nmem_dynamic; ++e)
                                ((float***) vvv)[e][j][i] = v_a[e];
                        else {
                            if (!applylog)
                                for (e = 0; e < nmem_dynamic; ++e)
                                    ((float***) vvv)[e][j][i] = v_a[e] - v_f[e];
                            else
                                for (e = 0; e < nmem_dynamic; ++e) {
                                    if (!isnormal(v_f[e]))
                                        ((float***) vvv)[e][j][i] = 0.0;
                                    else
                                        /*
                                         * it is necessary to take care that in
                                         * this case (--output-increment) the
                                         * variable is not transformed again
                                         * during writing
                                         */
                                        ((float***) vvv)[e][j][i] = exp10(v_a[e]) - exp10(v_f[e]);
                                }
                        }
                    } else {
                        if (!(das->updatespec & UPDATE_OUTPUTINC))
                            for (e = 0; e < nmem_dynamic; ++e)
                                ((float**) vvv)[e][j] = v_a[e];
                        else {
                            if (!applylog)
                                for (e = 0; e < nmem_dynamic; ++e)
                                    ((float**) vvv)[e][j] = v_a[e] - v_f[e];
                            else
                                for (e = 0; e < nmem_dynamic; ++e) {
                                    if (!isnormal(v_f[e]))
                                        ((float**) vvv)[e][j] = 0.0;
                                    else
                                        /*
                                         * it is necessary to take care that in
                                         * this case (--output-increment) the
                                         * variable is not transformed again
                                         * during writing
                                         */
                                        ((float**) vvv)[e][j] = exp10(v_a[e]) - exp10(v_f[e]);
                                }
                        }
                    }

                    if (writeinflation) {
                        if (structured)
                            ((float*) infl)[i] = inflation;
                        else
                            ((float**) infl)[fid][j] = inflation;
                    }
                }               /* for i */
                if (writeinflation && structured)
                    das_writeinflation(das, f, j, infl);
            }                   /* for fid */
        }                       /* for stepj */
    }                           /* for jj */
    if (writeinflation && !structured)
        for (fid = 0; fid < nfields; ++fid)
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
    for (fid = 0; fid < nfields; ++fid) {
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
static void das_updatebg(dasystem* das, int nfields, void** fieldbuffer, field fields[])
{
    model* m = das->m;
    int nmem = das->nmem;

    void* g = model_getvargrid(m, fields[0].varid);
    int gridid = grid_getid(g);
    int stride = grid_getstride(g);
    void* nlevels = grid_getnumlevels(g);
    int surfk = grid_getsurflayerid(g);
    int periodic_i = grid_isperiodic_i(g);
    int structured = 1;

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
        structured = 0;
    }

    das_getfname_transforms(das, gridid, fname);

    ncw_open(fname, NC_NOWRITE, &ncid);
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
                void* vvv = fieldbuffer[fid];

                for (i = 0; i < mni; ++i) {
                    float xmean = 0.0f;

                    if (structured == 1) {
                        if (surfk == 0) {
                            if (((int**) nlevels)[j][i] <= f->level) {
                                ((float***) vvv)[nmem][j][i] = 0.0f;
                                continue;
                            }
                        } else {
                            if (((int**) nlevels)[j][i] <= ((f->issurfacevar) ? 1 : surfk - f->level)) {
                                ((float***) vvv)[nmem][j][i] = 0.0f;
                                continue;
                            }
                        }
                        /*
                         * assume that if |value| > MAXOBSVAL, then it is filled
                         * with the missing value 
                         */
                        if (fabsf(((float***) vvv)[nmem][j][i]) > (float) MAXOBSVAL) {
                            ((float***) vvv)[nmem][j][i] = 0.0f;
                            continue;
                        }
                        for (e = 0; e < nmem; ++e)
                            if (fabsf(((float***) vvv)[e][j][i]) > (float) MAXOBSVAL)
                                break;
                        if (e < nmem) {
                            ((float***) vvv)[nmem][j][i] = 0.0f;
                            continue;
                        }

                        for (e = 0; e < nmem; ++e)
                            xmean += ((float***) vvv)[e][j][i];
                        xmean /= (float) nmem;

                        /*
                         * (the case das->updatespec & UPDATE_OUTPUTINC > 0 is
                         * handled by setting vvv[nmem][][] to zero in
                         * das_update())
                         */
                        for (e = 0; e < nmem; ++e)
                            ((float***) vvv)[nmem][j][i] += (((float***) vvv)[e][j][i] - xmean) * wj[i][e];
                    } else {
                        if (surfk == 0) {
                            if (((int*) nlevels)[j] <= f->level) {
                                ((float**) vvv)[nmem][j] = 0.0f;
                                continue;
                            }
                        } else {
                            if (((int*) nlevels)[j] <= ((f->issurfacevar) ? 1 : surfk - f->level)) {
                                ((float**) vvv)[nmem][j] = 0.0f;
                                continue;
                            }
                        }
                        /*
                         * assume that if |value| > MAXOBSVAL, then it is filled
                         * with the missing value 
                         */
                        if (fabsf(((float**) vvv)[nmem][j]) > (float) MAXOBSVAL) {
                            ((float**) vvv)[nmem][j] = 0.0f;
                            continue;
                        }
                        for (e = 0; e < nmem; ++e)
                            if (fabsf(((float**) vvv)[e][j]) > (float) MAXOBSVAL)
                                break;
                        if (e < nmem) {
                            ((float**) vvv)[nmem][j] = 0.0f;
                            continue;
                        }

                        for (e = 0; e < nmem; ++e)
                            xmean += ((float**) vvv)[e][j];
                        xmean /= (float) nmem;

                        /*
                         * (the case das->updatespec & UPDATE_OUTPUTINC > 0 is
                         * handled by setting vvv[nmem][][] to zero in
                         * das_update())
                         */
                        for (e = 0; e < nmem; ++e)
                            ((float**) vvv)[nmem][j] += (((float**) vvv)[e][j] - xmean) * wj[0][e];
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
static void das_writefields_direct(dasystem* das, int nfields, void** fieldbuffer, field fields[])
{
    int i, e;

    if (!(das->updatespec & UPDATE_SEPARATEOUTPUT)) {
        for (i = 0; i < nfields; ++i) {
            field* f = &fields[i];
            char varname[NC_MAX_NAME];

            strncpy(varname, f->varname, NC_MAX_NAME - 1);
            if (!(das->updatespec & UPDATE_OUTPUTINC))
                strncat(varname, "_an", NC_MAX_NAME - 1);
            else
                strncat(varname, "_inc", NC_MAX_NAME - 1);

            for (e = 0; e < das->nmem_dynamic; ++e) {
                char fname[MAXSTRLEN];

                das_getmemberfname(das, f->varname, e + 1, fname);
                model_writefieldas(das->m, fname, varname, f->varname, f->level, (f->structured) ? ((float***) fieldbuffer[i])[e][0] : ((float**) fieldbuffer[i])[e], 0);
            }
        }
    } else {
        for (i = 0; i < nfields; ++i) {
            field* f = &fields[i];

            for (e = 0; e < das->nmem_dynamic; ++e) {
                char fname[MAXSTRLEN];

                das_getmemberfname(das, f->varname, e + 1, fname);
                if (!(das->updatespec & UPDATE_OUTPUTINC)) {
                    strncat(fname, ".analysis", MAXSTRLEN - 1);
                    model_writefield(das->m, fname, f->varname, f->level, (f->structured) ? ((float***) fieldbuffer[i])[e][0] : ((float**) fieldbuffer[i])[e], 0);
                } else {
                    strncat(fname, ".increment", MAXSTRLEN - 1);
                    model_writefield(das->m, fname, f->varname, f->level, (f->structured) ? ((float***) fieldbuffer[i])[e][0] : ((float**) fieldbuffer[i])[e], 1);
                }
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
        int ncid;
        int dimids[3];
        int vid;

        getfieldfname(DIRNAME_TMP, "ens", f->varname, f->level, fname);

        ncw_create(fname, NC_CLOBBER | das->ncformat, &ncid);
        ncw_def_dim(ncid, "m", das->nmem, &dimids[0]);
        if (f->structured) {
            ncw_def_dim(ncid, "nj", nj, &dimids[1]);
            ncw_def_dim(ncid, "ni", ni, &dimids[2]);
            ncw_def_var(ncid, f->varname, NC_FLOAT, 3, dimids, &vid);
        } else {
            ncw_def_dim(ncid, "ni", ni, &dimids[1]);
            ncw_def_var(ncid, f->varname, NC_FLOAT, 2, dimids, &vid);
        }
#if defined(DEFLATE_ALL)
        if (das->nccompression > 0)
            ncw_def_deflate(ncid, 0, 1, das->nccompression);
#endif
        ncw_enddef(ncid);
        ncw_put_var_float(ncid, vid, (f->structured) ? ((float***) fieldbuffer[i])[0][0] : ((float**) fieldbuffer[i])[0]);
        ncw_close(ncid);
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

            das_getbgfname(das, f->varname, fname);
            strncpy(varname, f->varname, NC_MAX_NAME - 1);
            if (!(das->updatespec & UPDATE_OUTPUTINC))
                strncat(varname, "_an", NC_MAX_NAME - 1);
            else
                strncat(varname, "_inc", NC_MAX_NAME - 1);
            model_writefield(m, fname, varname, f->level, (f->structured) ? ((float***) fieldbuffer[i])[das->nmem][0] : ((float**) fieldbuffer[i])[das->nmem], 0);
        }
    } else {
        for (i = 0; i < nfields; ++i) {
            field* f = &fields[i];
            char fname[MAXSTRLEN];

            das_getbgfname(das, f->varname, fname);
            if (!(das->updatespec & UPDATE_OUTPUTINC))
                strncat(fname, ".analysis", MAXSTRLEN - 1);
            else
                strncat(fname, ".increment", MAXSTRLEN - 1);
            model_writefield(m, fname, f->varname, f->level, (f->structured) ? ((float***) fieldbuffer[i])[das->nmem][0] : ((float**) fieldbuffer[i])[das->nmem], 0);
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

        getfieldfname(DIRNAME_TMP, "bg", f->varname, f->level, fname);

        if (!file_exists(fname)) {
            int ncid;
            int dimids[2];
            int vid;

            ncw_create(fname, NC_CLOBBER | das->ncformat, &ncid);
            if (f->structured) {
                ncw_def_dim(ncid, "nj", nj, &dimids[0]);
                ncw_def_dim(ncid, "ni", ni, &dimids[1]);
                ncw_def_var(ncid, f->varname, NC_FLOAT, 2, dimids, &vid);
            } else {
                ncw_def_dim(ncid, "ni", ni, &dimids[0]);
                ncw_def_var(ncid, f->varname, NC_FLOAT, 1, dimids, &vid);
            }
#if defined(DEFLATE_ALL)
            if (das->nccompression > 0)
                ncw_def_deflate(ncid, 0, 1, das->nccompression);
#endif
            ncw_enddef(ncid);
            ncw_close(ncid);
        }
        /*
         * ignorelog = 1 here, as it will be handled during assembling if
         * necessary
         */
        model_writefield(das->m, fname, f->varname, f->level, (f->structured) ? ((float***) fieldbuffer[i])[das->nmem][0] : ((float**) fieldbuffer[i])[das->nmem], 1);
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
    int i;

#if defined(MPI)
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    distribute_iterations(0, das->nmem_dynamic - 1, nprocesses, rank, "    ");

    for (i = 0; i < nvar; ++i) {
        char* varname = model_getvarname(m, i);
        grid* g = model_getvargrid(m, i);
        char varname_dst[NC_MAX_NAME];
        char fname_dst[MAXSTRLEN];
        int nlev, k;
        int ni, nj;
        float* v = NULL;
        int e;

        enkf_printf("    %s:", varname);
        enkf_flush();
        das_getmemberfname(das, varname, 1, fname_dst);
        strncpy(varname_dst, varname, NC_MAX_NAME - 1);

        grid_getsize(g, &ni, &nj, NULL);
        if (nj > 0)
            v = malloc(ni * nj * sizeof(float));
        else
            v = malloc(ni * sizeof(float));

        nlev = ncu_getnlevels(fname_dst, varname, nj > 0);

        for (e = my_first_iteration; e <= my_last_iteration; ++e) {
            das_getmemberfname(das, varname, e + 1, fname_dst);
            if (das->updatespec & UPDATE_SEPARATEOUTPUT) {
                if (!(das->updatespec & UPDATE_OUTPUTINC))
                    strncat(fname_dst, ".analysis", MAXSTRLEN - 1);
                else
                    strncat(fname_dst, ".increment", MAXSTRLEN - 1);
            } else {
                if (!(das->updatespec & UPDATE_OUTPUTINC))
                    strncat(varname_dst, "_an", NC_MAX_NAME - 1);
                else
                    strncat(varname_dst, "_inc", NC_MAX_NAME - 1);
            }

            for (k = 0; k < nlev; ++k) {
                char fname_src[MAXSTRLEN];
                int ncid_src, vid_src;
                size_t start[3] = { e, 0, 0 };
                size_t count[3] = { 1, nj, ni };

                if (nj <= 0)
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
            enkf_printf(".");
            enkf_flush();
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
    if (rank == 0) {
        for (i = 0; i < nvar; ++i) {
            char* varname = model_getvarname(m, i);
            grid* g = model_getvargrid(m, i);
            char fname[MAXSTRLEN];
            int nlev, k;

            das_getmemberfname(das, varname, 1, fname);
            nlev = ncu_getnlevels(fname, varname, grid_isstructured(g));
            for (k = 0; k < nlev; ++k) {
                getfieldfname(DIRNAME_TMP, "ens", varname, k, fname);
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
        if (das->updatespec & UPDATE_SEPARATEOUTPUT) {
            if (!(das->updatespec & UPDATE_OUTPUTINC))
                strncat(fname_dst, ".analysis", MAXSTRLEN - 1);
            else
                strncat(fname_dst, ".increment", MAXSTRLEN - 1);
        } else {
            if (!(das->updatespec & UPDATE_OUTPUTINC))
                strncat(varname_dst, "_an", NC_MAX_NAME - 1);
            else
                strncat(varname_dst, "_inc", NC_MAX_NAME - 1);
        }

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

    if (das->updatespec & UPDATE_DOSPREAD && rank == 0) {
        enkf_printf("    allocating disk space for spread:");
        das_allocatespread(das, FNAME_SPREAD);
        enkf_printf("\n");
        enkf_flush();
        if (rank == 0 && !(das->updatespec & UPDATE_DIRECTWRITE))
            dir_createifabsent(DIRNAME_TMP);
    }

    if (das->updatespec & UPDATE_DOINFLATION && rank == 0) {
        enkf_printf("    allocating disk space for inflation:");
        das_allocateinflation(das, FNAME_INFLATION);
        enkf_printf("\n");
        enkf_flush();
        if (rank == 0 && !(das->updatespec & UPDATE_DIRECTWRITE))
            dir_createifabsent(DIRNAME_TMP);
    }

    if (das->updatespec & UPDATE_DOPLOGS && rank == 0) {
        enkf_printf("    defining state variables in point logs:");
        plog_definestatevars(das);
        enkf_printf("\n");
        enkf_flush();
    }

    if (rank == 0 && (das->updatespec | UPDATE_DOPLOGS))
        dir_createifabsent(DIRNAME_TMP);
    if (das->updatespec & UPDATE_DOFIELDS) {
        if (das->mode == MODE_ENKF || das->mode == MODE_HYBRID) {
            int i, e;

            distribute_iterations(0, das->nmem_dynamic - 1, nprocesses, rank, "    ");

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
                            strncat(varname_dst, "_an", NC_MAX_NAME - 1);
                        else
                            strncat(varname_dst, "_inc", NC_MAX_NAME - 1);

                        das_getmemberfname(das, varname_src, e + 1, fname);
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

                        das_getmemberfname(das, varname, e + 1, fname_f);
                        ncw_open(fname_f, NC_NOWRITE, &ncid_f);

                        strncpy(fname_a, fname_f, MAXSTRLEN);
                        if (!(das->updatespec & UPDATE_OUTPUTINC))
                            strncat(fname_a, ".analysis", MAXSTRLEN - 1);
                        else
                            strncat(fname_a, ".increment", MAXSTRLEN - 1);
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
            int i;

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
                            strncat(varname_dst, "_an", NC_MAX_NAME - 1);
                        else
                            strncat(varname_dst, "_inc", NC_MAX_NAME - 1);
                        das_getbgfname(das, varname_src, fname);
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

                        das_getbgfname(das, varname, fname_f);
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
                        } else
                            ncw_create(fname_a, NC_CLOBBER | das->ncformat, &ncid_a);

                        if (!file_exists(fname_f) && (das->updatespec & UPDATE_OUTPUTINC))
                            das_getmemberfname(das, varname, 1, fname_f);
                        ncw_open(fname_f, NC_NOWRITE, &ncid_f);
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
#if defined(USE_MPIQUEUE)
    {
        int nfields = 0;
        field* fields = NULL;
        mpiqueue* queue = NULL;

        das_getfields(das, -1, &nfields, &fields);
        if (nprocesses == 1)
            enkf_quit("\"mpiqueue\" can not be used on a single CPU; run on more than one CPU or recompile without -DUSE_MPIQUEUE flag");
        queue = mpiqueue_create(MPI_COMM_WORLD, nfields);
        enkf_printf("    updating %d fields using %d processes:\n", nfields, nprocesses);

        if (mpiqueue_getrank(queue) == 0)
            mpiqueue_manage(queue);
        else {
            void** fieldbuffer = calloc(das->fieldbufsize, sizeof(void*));
            field* fieldstowrite = malloc(das->fieldbufsize * sizeof(field));
            int bufid = 0;
            int gridid_storage = -1;
            int gridid_buf;
            int fid;
            int i, e;

            while (1) {
                field* f;
                int gridid_field;

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
                if (gridid_buf != gridid_storage) {
                    grid* g = model_getgridbyid(m, gridid_field);
                    int mni, mnj;

                    grid_getsize(g, &mni, &mnj, NULL);
                    for (i = 0; i < das->fieldbufsize; ++i)
                        if (fieldbuffer[i] != NULL)
                            free(fieldbuffer[i]);
                    if (mnj > 0) {
                        for (i = 0; i < das->fieldbufsize; ++i) {
                            if (das->mode == MODE_ENKF)
                                fieldbuffer[i] = alloc3d(das->nmem, mnj, mni, sizeof(float));
                            else if (das->mode == MODE_ENOI)
                                fieldbuffer[i] = alloc3d(das->nmem + 1, mnj, mni, sizeof(float));
                            else if (das->mode == MODE_HYBRID)
                                /*
                                 * allocate two additional members to calculate
                                 * ensemble mean with double precision
                                 */
                                fieldbuffer[i] = alloc3d(das->nmem + 2, mnj, mni, sizeof(float));
                        }
                    } else {
                        for (i = 0; i < das->fieldbufsize; ++i) {
                            if (das->mode == MODE_ENKF)
                                fieldbuffer[i] = alloc2d(das->nmem, mni, sizeof(float));
                            else if (das->mode == MODE_ENOI)
                                fieldbuffer[i] = alloc2d(das->nmem + 1, mni, sizeof(float));
                            else if (das->mode == MODE_HYBRID)
                                /*
                                 * allocate two additional members to calculate
                                 * ensemble mean with double precision
                                 */
                                fieldbuffer[i] = alloc2d(das->nmem + 2, mni, sizeof(float));
                        }
                    }
                    gridid_storage = gridid_buf;
                }
                /*
                 * read ensemble of fields
                 */
                if (enkf_verbose) {
                    printf("      %-8s %-3d (%d: %d: %.1f%%)\n", f->varname, f->level, rank, fid, 100.0 * (double) (fid + 1) / (double) nfields);
                    fflush(stdout);
                }
                for (e = 0; e < das->nmem; ++e) {
                    char fname[MAXSTRLEN];
                    int masklog = das_isstatic(das, e + 1);

                    das_getmemberfname(das, f->varname, e + 1, fname);
                    model_readfield(das->m, fname, f->varname, f->level, (f->structured) ? ((float***) fieldbuffer[bufid])[e][0] : ((float**) fieldbuffer[bufid])[e], masklog);
                }
                if (das->mode == MODE_HYBRID) {
                    float* v[das->nmem + 1];
                    int mni, mnj;

                    model_getvargridsize(das->m, f->varid, &mni, &mnj, NULL);
                    for (e = 0; e < das->nmem + 1; ++e)
                        v[e] = (f->structured) ? ((float***) fieldbuffer[bufid])[e][0] : ((float**) fieldbuffer[bufid])[e];

                    das_sethybridensemble(das, (f->structured) ? mni * mnj : mni, v);
                }
                /*
                 * read the background to write it to pointlogs, regardless
                 * of whether output is increment or analysis
                 */
                if (das->mode == MODE_ENOI && (((das->updatespec & UPDATE_DOFIELDS) && !(das->updatespec & UPDATE_OUTPUTINC)) || (das->updatespec & UPDATE_DOPLOGSAN))) {
                    char fname[MAXSTRLEN];

                    das_getbgfname(das, f->varname, fname);
                    model_readfield(das->m, fname, f->varname, f->level, (f->structured) ? ((float***) fieldbuffer[bufid])[das->nmem][0] : ((float**) fieldbuffer[bufid])[das->nmem], 0);
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
                if (das->updatespec & UPDATE_DOFORECASTSPREAD)
                    das_writespread_inupdate(das, bufid + 1, fieldbuffer, fieldstowrite, 0);
                /*
                 * write forecast variables to point logs
                 */
                if (das->updatespec & UPDATE_DOPLOGSFC)
                    plog_writestatevars(das, bufid + 1, fieldbuffer, fieldstowrite, 0);
                /*
                 * set the background to 0 if output is increment
                 */
                if (das->mode == MODE_ENOI && (das->updatespec & UPDATE_OUTPUTINC)) {
                    int mni, mnj;
                    int ii;

                    model_getvargridsize(das->m, f->varid, &mni, &mnj, NULL);
                    for (ii = 0; ii <= bufid; ++ii)
                        memset(((float***) fieldbuffer[ii])[das->nmem][0], 0, mni * mnj * sizeof(float));
                }
                if (das->updatespec & (UPDATE_DOFIELDS | UPDATE_DOANALYSISSPREAD | UPDATE_DOPLOGSAN | UPDATE_DOINFLATION)) {
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
                /*
                 * write analysis spread
                 */
                if (bufid >= 0 && das->updatespec & UPDATE_DOANALYSISSPREAD && (das->mode == MODE_ENKF || das->mode == MODE_HYBRID))
                    das_writespread_inupdate(das, bufid + 1, fieldbuffer, fieldstowrite, 1);
                /*
                 * write analysis variables to point logs
                 */
                if (bufid >= 0 && das->updatespec & UPDATE_DOPLOGSAN)
                    plog_writestatevars(das, bufid + 1, fieldbuffer, fieldstowrite, 1);

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
    }
#else                           /* !defined(USE_MPIQUEUE) */
    for (gid = 0; gid < ngrid; ++gid) {
        void* g = model_getgridbyid(m, gid);
        int nfields = 0;
        field* fields = NULL;
        void** fieldbuffer = NULL;
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

        grid_getsize(g, &mni, &mnj, NULL);

        das_getfields(das, gid, &nfields, &fields);
        enkf_printf("      %d fields\n", nfields);

        if (nfields == 0)
            continue;

        distribute_iterations(0, nfields - 1, nprocesses, rank, "      ");

        if (my_first_iteration > my_last_iteration)
            continue;

        fieldbuffer = malloc(das->fieldbufsize * sizeof(void*));
        if (mnj > 0) {
            if (das->mode == MODE_ENKF) {
                for (i = 0; i < das->fieldbufsize; ++i)
                    fieldbuffer[i] = alloc3d(das->nmem, mnj, mni, sizeof(float));
            } else if (das->mode == MODE_ENOI) {
                for (i = 0; i < das->fieldbufsize; ++i)
                    fieldbuffer[i] = alloc3d(das->nmem + 1, mnj, mni, sizeof(float));
            } else if (das->mode == MODE_HYBRID) {
                /*
                 * allocate two additional members to calculate ensemble mean
                 * with double precision
                 */
                for (i = 0; i < das->fieldbufsize; ++i)
                    fieldbuffer[i] = alloc3d(das->nmem + 2, mnj, mni, sizeof(float));
            }
        } else {
            if (das->mode == MODE_ENKF) {
                for (i = 0; i < das->fieldbufsize; ++i)
                    fieldbuffer[i] = alloc2d(das->nmem, mni, sizeof(float));
            } else if (das->mode == MODE_ENOI) {
                for (i = 0; i < das->fieldbufsize; ++i)
                    fieldbuffer[i] = alloc2d(das->nmem + 1, mni, sizeof(float));
            } else if (das->mode == MODE_HYBRID) {
                /*
                 * allocate two additional members to calculate ensemble mean
                 * with double precision
                 */
                for (i = 0; i < das->fieldbufsize; ++i)
                    fieldbuffer[i] = alloc2d(das->nmem + 2, mni, sizeof(float));
            }
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
                int masklog = das_isstatic(das, e + 1);

                das_getmemberfname(das, f->varname, e + 1, fname);
                model_readfield(das->m, fname, f->varname, f->level, (f->structured) ? ((float***) fieldbuffer[bufid])[e][0] : ((float**) fieldbuffer[bufid])[e], masklog);
            }

            if (das->mode == MODE_HYBRID) {
                float* v[das->nmem + 1];

                for (e = 0; e < das->nmem + 1; ++e)
                    v[e] = (f->structured) ? ((float***) fieldbuffer[bufid])[e][0] : ((float**) fieldbuffer[bufid])[e];

                das_sethybridensemble(das, (f->structured) ? mni * mnj : mni, v);
            }

            /*
             * read the background to write it to pointlogs, regardless of
             * whether output is increment or analysis
             */
            if (das->mode == MODE_ENOI && (((das->updatespec & UPDATE_DOFIELDS) && !(das->updatespec & UPDATE_OUTPUTINC)) || (das->updatespec & UPDATE_DOPLOGSAN))) {
                das_getbgfname(das, f->varname, fname);
                model_readfield(das->m, fname, f->varname, f->level, (f->structured) ? ((float***) fieldbuffer[bufid])[das->nmem][0] : ((float**) fieldbuffer[bufid])[das->nmem], 0);
            }

            if (bufid == das->fieldbufsize - 1 || fid == my_last_iteration) {
                /*
                 * write forecast spread
                 */
                if (das->updatespec & UPDATE_DOFORECASTSPREAD)
                    das_writespread_inupdate(das, bufid + 1, fieldbuffer, &fields[fid - bufid], 0);

                /*
                 * write forecast variables to point logs
                 */
                if (das->updatespec & UPDATE_DOPLOGSFC)
                    plog_writestatevars(das, bufid + 1, fieldbuffer, &fields[fid - bufid], 0);

                /*
                 * now set the background to 0 if output is increment
                 */
                if (das->mode == MODE_ENOI && (das->updatespec & UPDATE_OUTPUTINC)) {
                    int ii;

                    for (ii = 0; ii <= bufid; ++ii)
                        memset(((float***) fieldbuffer[ii])[das->nmem][0], 0, mni * mnj * sizeof(float));
                }

                if (das->updatespec & (UPDATE_DOFIELDS | UPDATE_DOANALYSISSPREAD | UPDATE_DOPLOGSAN | UPDATE_DOINFLATION)) {
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

                /*
                 * write analysis spread
                 */
                if (das->updatespec & UPDATE_DOANALYSISSPREAD && (das->mode == MODE_ENKF || das->mode == MODE_HYBRID))
                    das_writespread_inupdate(das, bufid + 1, fieldbuffer, &fields[fid - bufid], 1);
                /*
                 * write analysis variables to point logs
                 */
                if (das->updatespec & UPDATE_DOPLOGSAN)
                    plog_writestatevars(das, bufid + 1, fieldbuffer, &fields[fid - bufid], 1);
            }
        }                       /* for fid */

        for (i = 0; i < das->fieldbufsize; ++i)
            free(fieldbuffer[i]);
        free(fieldbuffer);
        free(fields);

        enkf_flush();
    }                           /* for gid */
#endif                          /* USE_MPIQUEUE */

#if defined(MPI)
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    if (!(das->updatespec & UPDATE_DIRECTWRITE)) {
        if (das->updatespec & UPDATE_DOFIELDS) {
            enkf_printtime("  ");
            enkf_printf("  assembling analysis:\n");
            if (das->mode == MODE_ENKF || das->mode == MODE_HYBRID)
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
}
