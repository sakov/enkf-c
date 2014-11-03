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
#include <values.h>
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

/** Updates `nfield' fields read into `fieldbuffer' with `X5'. Applies
 * variable-dependent inflation to each field.
 * 
 */
static void das_updatefields(dasystem* das, int nfields, void** fieldbuffer, field* fields)
{
    model* m = das->m;
    int** mask = model_getnumlevels(m);
    int periodic_i = grid_isperiodic_x(model_getgrid(m));
    int periodic_j = grid_isperiodic_y(model_getgrid(m));
    /*
     * X5
     */
    int ncid;
    int varid;
    int dimids[3];
    size_t dimlens[3];
    float** X5jj;
    float** X5jj1 = NULL;
    float** X5jj2 = NULL;
    float** X5j;
    size_t start[3], count[3];
    /*
     * SRF
     */
    int ncid_srf;
    int varid_srf;
    float* SRFjj;
    float* SRFjj1 = NULL;
    float* SRFjj2 = NULL;
    float* SRFj;
    size_t start_srf[2], count_srf[2];

    int ni, nj;

    int i, j, mni, mnj;
    int jj, stepj, ii, stepi;
    int e, fid;
    float* v;                   /* v = E(<some index>, :) */

    float* tmp;

    assert(das->mode == MODE_ENKF);

    model_getdims(m, &mni, &mnj, NULL);
    tmp = malloc(das->nmem * sizeof(float));
    v = malloc(das->nmem * sizeof(float));

    ncw_open(FNAME_X5, NC_NOWRITE, &ncid);
    ncw_inq_varid(FNAME_X5, ncid, "X5", &varid);
    ncw_inq_vardimid(FNAME_X5, ncid, varid, dimids);
    for (i = 0; i < 3; ++i)
        ncw_inq_dimlen(FNAME_X5, ncid, dimids[i], &dimlens[i]);
    nj = dimlens[0];
    ni = dimlens[1];
    assert((int) dimlens[2] == das->nmem * das->nmem);

    ncw_open(FNAME_ENKFSTATS, NC_NOWRITE, &ncid_srf);
    ncw_inq_varid(FNAME_ENKFSTATS, ncid_srf, "srf", &varid_srf);

    start[0] = 0;
    start[1] = 0;
    start[2] = 0;
    count[0] = 1;
    count[1] = ni;
    count[2] = das->nmem * das->nmem;
    X5j = alloc2d(mni, das->nmem * das->nmem, sizeof(float));
    if (das->stride > 1) {
        X5jj = alloc2d(ni, das->nmem * das->nmem, sizeof(float));
        X5jj1 = alloc2d(ni, das->nmem * das->nmem, sizeof(float));
        X5jj2 = alloc2d(ni, das->nmem * das->nmem, sizeof(float));
        ncw_get_vara_float(FNAME_X5, ncid, varid, start, count, X5jj2[0]);
    }

    start_srf[0] = 0;
    start_srf[1] = 0;
    count_srf[0] = 1;
    count_srf[1] = ni;
    SRFj = calloc(mni, sizeof(float));
    if (das->stride > 1) {
        SRFjj = calloc(ni, sizeof(float));
        SRFjj1 = calloc(ni, sizeof(float));
        SRFjj2 = calloc(ni, sizeof(float));
        ncw_get_vara_float(FNAME_ENKFSTATS, ncid_srf, varid_srf, start_srf, count_srf, SRFjj2);
    }

    /*
     * jj, ii are the indices of the subsampled grid; i, j are the indices of
     * the actual model grid 
     */
    for (jj = 0, j = 0; jj < nj; ++jj) {
        for (stepj = 0; stepj < das->stride && j < mnj; ++stepj, ++j) {
            if (das->stride == 1) {
                /*
                 * no interpolation necessary; simply read the ETMs for the
                 * j-th row from disk 
                 */
                start[0] = j;
                ncw_get_vara_float(FNAME_X5, ncid, varid, start, count, X5j[0]);
                start_srf[0] = j;
                ncw_get_vara_float(FNAME_ENKFSTATS, ncid_srf, varid_srf, start_srf, count_srf, SRFj);
            } else {
                /*
                 * the following code interpolates the ETM back to the
                 * original grid, first by j, and then by i 
                 */
                if (stepj == 0) {
                    memcpy(X5jj1[0], X5jj2[0], ni * das->nmem * das->nmem * sizeof(float));
                    memcpy(X5jj[0], X5jj2[0], ni * das->nmem * das->nmem * sizeof(float));
                    memcpy(SRFjj1, SRFjj2, ni * sizeof(float));
                    memcpy(SRFjj, SRFjj2, ni * sizeof(float));
                    if (jj < nj - 1 || periodic_j) {
                        start[0] = (jj + 1) % nj;
                        ncw_get_vara_float(FNAME_X5, ncid, varid, start, count, X5jj2[0]);
                        start_srf[0] = (jj + 1) % nj;
                        ncw_get_vara_float(FNAME_ENKFSTATS, ncid_srf, varid_srf, start_srf, count_srf, SRFjj2);
                    }
                } else {
                    float weight2 = (float) stepj / (float) das->stride;
                    float weight1 = 1.0f - weight2;

                    for (ii = 0; ii < ni; ++ii) {
                        float* X5jjii = X5jj[ii];
                        float* X5jj1ii = X5jj1[ii];
                        float* X5jj2ii = X5jj2[ii];

                        for (e = 0; e < das->nmem * das->nmem; ++e)
                            X5jjii[e] = X5jj1ii[e] * weight1 + X5jj2ii[e] * weight2;
                    }
                    for (ii = 0; ii < ni; ++ii)
                        SRFjj[ii] = SRFjj1[ii] * weight1 + SRFjj2[ii] * weight2;
                }

                for (ii = 0, i = 0; ii < ni; ++ii) {
                    for (stepi = 0; stepi < das->stride && i < mni; ++stepi, ++i) {
                        if (stepi == 0) {
                            memcpy(X5j[i], X5jj[ii], das->nmem * das->nmem * sizeof(float));
                            SRFj[i] = SRFjj[ii];
                        } else {
                            float weight2 = (float) stepi / (float) das->stride;
                            float weight1 = 1.0f - weight2;
                            float* X5jjii1 = X5jj[ii];
                            float* X5ji = X5j[i];
                            float* X5jjii2;
                            float SRFjjii2;

                            if (ii < ni - 1)
                                X5jjii2 = X5jj[ii + 1];
                            else
                                X5jjii2 = X5jj[(periodic_i) ? (ii + 1) % ni : ii];
                            for (e = 0; e < das->nmem * das->nmem; ++e)
                                X5ji[e] = X5jjii1[e] * weight1 + X5jjii2[e] * weight2;

                            if (ii < ni - 1)
                                SRFjjii2 = SRFjj[ii + 1];
                            else
                                SRFjjii2 = SRFjj[(periodic_i) ? (ii + 1) % ni : ii];
                            SRFj[i] = SRFjj[ii] * weight1 + SRFjjii2 * weight2;
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
                float inflation0 = model_getvarinflation(m, f->varid);

                for (i = 0; i < mni; ++i) {
                    float inflation;

                    if (inflation0 < 1.0f + SRFj[i])
                        inflation = inflation0;
                    else
                        inflation = 1.0f + SRFj[i];

                    /*
                     * assume that if |value| > MAXOBSVAL, then it is filled
                     * with the missing value 
                     */
                    if (!mask[j][i])
                        continue;
                    if (fabsf(vvv[0][j][i]) > (float) MAXOBSVAL)
                        continue;
                    /*
                     * (it would be straightforward to compare with the
                     * actual missing value, provided that it is NC_FLOAT;
                     * otherwise it may be a bit tiresome to handle all
                     * variations) 
                     */

                    for (e = 0; e < das->nmem; ++e)
                        v[e] = vvv[e][j][i];

                    /*
                     * E(i, :) = E(i, :) * X5 
                     */
                    sgemv_(&do_T, &das->nmem, &das->nmem, &alpha, X5j[i], &das->nmem, v, &inc, &beta, tmp, &inc);

                    /*
                     * applying inflation explicitely is computationally
                     * cheaper than including it into X5, I guess 
                     */
                    if (fabsf(inflation - 1.0f) > EPSF) {
                        float v_av = 0.0f;

                        for (e = 0; e < das->nmem; ++e)
                            v_av += tmp[e];
                        v_av /= (float) das->nmem;
                        for (e = 0; e < das->nmem; ++e)
                            tmp[e] = (tmp[e] - v_av) * inflation + v_av;
                    }

                    if (das->target == TARGET_ANALYSIS)
                        for (e = 0; e < das->nmem; ++e)
                            vvv[e][j][i] = tmp[e];
                    else if (das->target == TARGET_INCREMENT)
                        for (e = 0; e < das->nmem; ++e)
                            vvv[e][j][i] = tmp[e] - v[e];
                }
            }
        }                       /* for stepj */
    }                           /* for jj */

    ncw_close(FNAME_X5, ncid);
    ncw_close(FNAME_ENKFSTATS, ncid_srf);

    free(tmp);
    free(v);
    free2d(X5j);
    free(SRFj);
    if (das->stride > 1) {
        free2d(X5jj);
        free2d(X5jj1);
        free2d(X5jj2);
        free(SRFjj);
        free(SRFjj1);
        free(SRFjj2);
    }
}

/** Updates `nfield' fields read into `fieldbuffer' with `w'. 
 */
static void das_updatebg(dasystem* das, int nfields, void** fieldbuffer)
{
    model* m = das->m;
    int** mask = model_getnumlevels(m);
    int periodic_i = grid_isperiodic_x(model_getgrid(m));
    int periodic_j = grid_isperiodic_y(model_getgrid(m));
    int ncid;
    int varid;
    int dimids[3];
    size_t dimlens[3];
    float** wjj = NULL;
    float** wjj1 = NULL;
    float** wjj2 = NULL;
    float** wj;
    int i, j, ni, nj, mni, mnj;
    int jj, stepj, ii, stepi;
    size_t start[3], count[3];
    int e, f;
    float* tmp;

    assert(das->mode == MODE_ENOI);

    /*
     * the following code for interpolation of X5 essentially coincides with
     * that in das_updatefields() 
     */

    model_getdims(m, &mni, &mnj, NULL);

    ncw_open(FNAME_W, NC_NOWRITE, &ncid);
    ncw_inq_varid(FNAME_W, ncid, "w", &varid);
    ncw_inq_vardimid(FNAME_W, ncid, varid, dimids);
    for (i = 0; i < 3; ++i)
        ncw_inq_dimlen(FNAME_W, ncid, dimids[i], &dimlens[i]);
    ni = dimlens[1];
    nj = dimlens[0];

    assert((int) dimlens[2] == das->nmem);

    start[0] = 0;
    start[1] = 0;
    start[2] = 0;
    count[0] = 1;
    count[1] = ni;
    count[2] = das->nmem;
    wj = alloc2d(mni, das->nmem, sizeof(float));
    if (das->stride > 1) {
        wjj = alloc2d(ni, das->nmem, sizeof(float));
        wjj1 = alloc2d(ni, das->nmem, sizeof(float));
        wjj2 = alloc2d(ni, das->nmem, sizeof(float));
        ncw_get_vara_float(FNAME_W, ncid, varid, start, count, wjj2[0]);
    }
    tmp = malloc(das->nmem * sizeof(float));
    /*
     * jj, ii are the indices of the subsampled grid; i, j are the indices of
     * the actual model grid 
     */
    for (jj = 0, j = 0; jj < nj; ++jj) {
        for (stepj = 0; stepj < das->stride && j < mnj; ++stepj, ++j) {
            if (das->stride == 1) {
                /*
                 * no interpolation necessary; simply read the ETMs for the
                 * j-th row from disk 
                 */
                start[0] = j;
                ncw_get_vara_float(FNAME_W, ncid, varid, start, count, wj[0]);
            } else {
                /*
                 * the following code interpolates the ETM back to the
                 * original grid, first by j, and then by i 
                 */
                if (stepj == 0) {
                    memcpy(wjj1[0], wjj2[0], ni * das->nmem * sizeof(float));
                    memcpy(wjj[0], wjj2[0], ni * das->nmem * sizeof(float));
                    if (jj < nj - 1 || periodic_j) {
                        start[0] = (jj + 1) % nj;
                        ncw_get_vara_float(FNAME_W, ncid, varid, start, count, wjj2[0]);
                    }
                } else {
                    float weight2 = (float) stepj / das->stride;
                    float weight1 = (float) 1.0f - weight2;

                    for (ii = 0; ii < ni; ++ii) {
                        float* wjjii = wjj[ii];
                        float* wjj1ii = wjj1[ii];
                        float* wjj2ii = wjj2[ii];

                        for (e = 0; e < das->nmem; ++e)
                            wjjii[e] = wjj1ii[e] * weight1 + wjj2ii[e] * weight2;
                    }
                }

                for (ii = 0, i = 0; ii < ni; ++ii) {
                    for (stepi = 0; stepi < das->stride && i < mni; ++stepi, ++i) {
                        if (stepi == 0)
                            memcpy(wj[i], wjj[ii], das->nmem * sizeof(float));
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

                            for (e = 0; e < das->nmem; ++e)
                                wji[e] = wjjii1[e] * weight1 + wjjii2[e] * weight2;
                        }
                    }
                }
            }                   /* stride != 1 */

            /*
             * (at this stage wj should contain the array of b vectors for
             * the j-th row of the grid) 
             */

            for (f = 0; f < nfields; ++f) {
                float*** vvv = (float***) fieldbuffer[f];

                for (i = 0; i < mni; ++i) {
                    float xmean = 0.0f;

                    /*
                     * assume that if |value| > MAXOBSVAL, then it is filled
                     * with the missing value 
                     */
                    if (!mask[j][i])
                        continue;
                    if (fabsf(vvv[das->nmem][j][i]) > (float) MAXOBSVAL)
                        continue;
                    if (fabsf(vvv[0][j][i]) > (float) MAXOBSVAL)
                        continue;

                    for (e = 0; e < das->nmem; ++e)
                        xmean += vvv[e][j][i];
                    xmean /= (float) das->nmem;

                    /*
                     * (the case das->target = TARGET_INCREMENT is handled by
                     * setting vvv[das->nmem][][] to zero in das_update())
                     */
                    for (e = 0; e < das->nmem; ++e)
                        vvv[das->nmem][j][i] += (vvv[e][j][i] - xmean) * wj[i][e];
                }
            }
        }                       /* for stepj */
    }                           /* for jj */

    ncw_close(FNAME_W, ncid);

    free(tmp);
    free2d(wj);
    if (das->stride > 1) {
        free2d(wjj);
        free2d(wjj1);
        free2d(wjj2);
    }

    das->s_mode = S_MODE_HE_a;
}

/**
 */
static void das_writefields_direct(dasystem* das, int nfields, void** fieldbuffer, field fields[])
{
    int i, e;

    assert(das->mode == MODE_ENKF);

    if (!enkf_separateout) {
        for (i = 0; i < nfields; ++i) {
            field* f = &fields[i];
            char varname[NC_MAX_NAME];

            strncpy(varname, f->varname, NC_MAX_NAME);
            if (das->target == TARGET_ANALYSIS)
                strncat(varname, "_an", NC_MAX_NAME);
            else if (das->target == TARGET_INCREMENT)
                strncat(varname, "_inc", NC_MAX_NAME);

            for (e = 0; e < das->nmem; ++e) {
                char fname[MAXSTRLEN];

                model_getmemberfname(das->m, das->ensdir, f->varname, e + 1, fname);
                model_writefield(das->m, fname, MAXINT, varname, f->level, ((float***) fieldbuffer[i])[e][0]);
            }
        }
    } else {
        for (i = 0; i < nfields; ++i) {
            field* f = &fields[i];

            for (e = 0; e < das->nmem; ++e) {
                char fname[MAXSTRLEN];

                model_getmemberfname(das->m, das->ensdir, f->varname, e + 1, fname);
                if (das->target == TARGET_ANALYSIS)
                    strncat(fname, ".analysis", MAXSTRLEN);
                else if (das->target == TARGET_INCREMENT)
                    strncat(fname, ".increment", MAXSTRLEN);
                model_writefield(das->m, fname, MAXINT, f->varname, f->level, ((float***) fieldbuffer[i])[e][0]);
            }
        }
    }

    if (das->nplogs > 0)
        plog_writestatevars_direct(das, nfields, fieldbuffer, fields);

}

/**
 */
static void getfieldfname(char* dir, char* prefix, char* varname, int level, char* fname)
{
    snprintf(fname, MAXSTRLEN, "%s/%s_%s-%03d.nc", dir, prefix, varname, level);
}

/**
 */
static void das_writefields_toassemble(dasystem* das, int nfields, void** fieldbuffer, field fields[])
{
    char fname[MAXSTRLEN];
    int ni, nj;
    int i;

    model_getdims(das->m, &ni, &nj, NULL);

    for (i = 0; i < nfields; ++i) {
        field* f = &fields[i];

        getfieldfname(das->ensdir, "ens", f->varname, f->level, fname);

        if (!file_exists(fname)) {
            int ncid;
            int dimids[3];
            int vid;

            ncw_create(fname, NC_CLOBBER | NC_64BIT_OFFSET, &ncid);
            ncw_def_dim(fname, ncid, "m", das->nmem, &dimids[0]);
            ncw_def_dim(fname, ncid, "nj", nj, &dimids[1]);
            ncw_def_dim(fname, ncid, "ni", ni, &dimids[2]);
            ncw_def_var(fname, ncid, f->varname, NC_FLOAT, 3, dimids, &vid);
            ncw_enddef(fname, ncid);
            ncw_put_var_float(fname, ncid, vid, ((float***) fieldbuffer[i])[0][0]);
            ncw_close(fname, ncid);
        } else {
            int ncid;
            int vid;

            ncw_open(fname, NC_WRITE, &ncid);
            ncw_inq_varid(fname, ncid, f->varname, &vid);
            ncw_put_var_float(fname, ncid, vid, ((float***) fieldbuffer[i])[0][0]);
            ncw_close(fname, ncid);
        }
    }

    if (das->nplogs > 0)
        plog_writestatevars_toassemble(das, nfields, fieldbuffer, fields);
}

/** Writes `nfields' ensemble fields from `fieldbuffer' to disk.
 */
static void das_writefields(dasystem* das, int nfields, void** fieldbuffer, field fields[])
{
    if (enkf_directwrite)
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

    model_getdims(m, &ni, &nj, NULL);

    if (!enkf_separateout) {
        for (i = 0; i < nfields; ++i) {
            field* f = &fields[i];
            char varname[NC_MAX_NAME];
            char fname[MAXSTRLEN];

            model_getbgfname(m, das->bgdir, f->varname, fname);
            strncpy(varname, f->varname, NC_MAX_NAME);
            if (das->target == TARGET_ANALYSIS)
                strncat(varname, "_an", NC_MAX_NAME);
            else if (das->target == TARGET_INCREMENT)
                strncat(varname, "_inc", NC_MAX_NAME);
            model_writefield(m, fname, MAXINT, varname, f->level, ((float***) fieldbuffer[i])[das->nmem][0]);
        }
    } else {
        for (i = 0; i < nfields; ++i) {
            field* f = &fields[i];
            char fname[MAXSTRLEN];

            model_getbgfname(m, das->bgdir, f->varname, fname);
            if (das->target == TARGET_ANALYSIS)
                strncat(fname, ".analysis", MAXSTRLEN);
            else if (das->target == TARGET_INCREMENT)
                strncat(fname, ".increment", MAXSTRLEN);
            model_writefield(m, fname, MAXINT, f->varname, f->level, ((float***) fieldbuffer[i])[das->nmem][0]);
        }
    }

    if (das->nplogs > 0)
        plog_writestatevars_direct(das, nfields, fieldbuffer, fields);

}

/**
 */
static void das_writebg_toassemble(dasystem* das, int nfields, void** fieldbuffer, field fields[])
{
    char fname[MAXSTRLEN];
    int ni, nj;
    int i;

    model_getdims(das->m, &ni, &nj, NULL);

    for (i = 0; i < nfields; ++i) {
        field* f = &fields[i];

        getfieldfname(das->bgdir, "bg", f->varname, f->level, fname);

        if (!file_exists(fname)) {
            int ncid;
            int dimids[2];
            int vid;

            ncw_create(fname, NC_CLOBBER | NC_64BIT_OFFSET, &ncid);
            ncw_def_dim(fname, ncid, "nj", nj, &dimids[0]);
            ncw_def_dim(fname, ncid, "ni", ni, &dimids[1]);
            ncw_def_var(fname, ncid, f->varname, NC_FLOAT, 2, dimids, &vid);
            ncw_enddef(fname, ncid);
            ncw_put_var_float(fname, ncid, vid, ((float***) fieldbuffer[i])[das->nmem][0]);
            ncw_close(fname, ncid);
        } else {
            writefield(fname, 0, f->varname, ((float***) fieldbuffer[i])[das->nmem][0]);
        }
    }

    if (das->nplogs > 0)
        plog_writestatevars_toassemble(das, nfields, fieldbuffer, fields);
}

/** Writes `nfield' fields from `fieldbuffer' to disk.
 */
static void das_writebg(dasystem* das, int nfields, void** fieldbuffer, field fields[])
{
    if (enkf_directwrite)
        das_writebg_direct(das, nfields, fieldbuffer, fields);
    else
        das_writebg_toassemble(das, nfields, fieldbuffer, fields);
}

/**
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

    ncw_create(fname, NC_CLOBBER | NC_64BIT_OFFSET, &ncid);
    for (vid = 0; vid < nvar; ++vid) {
        char* varname_src = model_getvarname(m, vid);
        int varid_src;

        model_getmemberfname(m, das->ensdir, varname_src, 1, fname_src);
        ncw_open(fname_src, NC_NOWRITE, &ncid_src);
        ncw_copy_dims(fname_src, ncid_src, fname, ncid);
        ncw_inq_varid(fname_src, ncid_src, varname_src, &varid_src);
        ncw_copy_vardef(fname_src, ncid_src, varid_src, fname, ncid);
        if (das->mode == MODE_ENKF) {
            char varname_dst[NC_MAX_NAME];

            strcpy(varname_dst, varname_src);
            strncat(varname_dst, "_an", NC_MAX_NAME);
            ncw_def_var_as(fname, ncid, varname_src, varname_dst);
        }
        ncw_close(fname_src, ncid_src);
    }
    ncw_close(fname, ncid);
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

    if (enkf_directwrite)
        strcpy(fname, FNAME_SPREAD);

    model_getdims(m, &ni, &nj, NULL);
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
            v2[i] = (v2[i] - (v1[i] * v1[i]) / (double) das->nmem) / (double) (das->nmem - 1);
            v2[i] = (v2[i] < 0.0) ? 0.0 : sqrt(v2[i]);
            if (fabs(v2[i]) > (double) MAXOBSVAL)
                v2[i] = NaN;
        }

        strncpy(varname, f->varname, NC_MAX_NAME);
        if (isanalysis)
            strncat(varname, "_an", MAXSTRLEN);

        if (!enkf_directwrite) {        /* create file for this field */
            int ncid, vid;
            int dimids[2];

            getfieldfname(das->mode == MODE_ENKF ? das->ensdir : das->bgdir, "spread", varname, f->level, fname);

            if (!file_exists(fname)) {
                ncw_create(fname, NC_CLOBBER | NC_64BIT_OFFSET, &ncid);
                ncw_def_dim(fname, ncid, "nj", nj, &dimids[0]);
                ncw_def_dim(fname, ncid, "ni", ni, &dimids[1]);
                ncw_def_var(fname, ncid, varname, NC_FLOAT, 2, dimids, &vid);
                ncw_enddef(fname, ncid);
            } else {
                ncw_open(fname, NC_WRITE, &ncid);
                if (!ncw_var_exists(ncid, varname)) {
                    ncw_redef(fname, ncid);
                    ncw_inq_dimid(fname, ncid, "nj", &dimids[0]);
                    ncw_inq_dimid(fname, ncid, "ni", &dimids[1]);
                    ncw_def_var(fname, ncid, varname, NC_FLOAT, 2, dimids, &vid);
                    ncw_enddef(fname, ncid);
                }
                ncw_inq_varid(fname, ncid, varname, &vid);
            }
            ncw_put_var_double(fname, ncid, vid, v2);
            ncw_close(fname, ncid);
        } else {
            float* v = calloc(nv, sizeof(float));

            for (i = 0; i < nv; ++i)
                v[i] = (float) v2[i];

            model_writefield(m, fname, MAXINT, varname, f->level, v);
            free(v);
        }
    }

    free(v1);
    free(v2);
}

/**
 */
static void das_assemblemembers(dasystem* das, int leavetiles)
{
    model* m = das->m;
    int nvar = model_getnvar(m);
    float** v = NULL;
    int ni, nj;
    int i, e;

#if defined(MPI)
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    model_getdims(m, &ni, &nj, NULL);
    v = alloc2d(ni, nj, sizeof(float));

    distribute_iterations(0, das->nmem - 1, nprocesses, rank);

    for (i = 0; i < nvar; ++i) {
        char* varname = model_getvarname(m, i);
        char varname_dst[NC_MAX_NAME];
        char fname_dst[MAXSTRLEN];
        int nlev, k;

        enkf_printf("    %s:", varname);
        model_getmemberfname(m, das->ensdir, varname, 1, fname_dst);
        nlev = getnlevels(fname_dst, varname);
        strncpy(varname_dst, varname, NC_MAX_NAME);

        for (e = my_first_iteration; e <= my_last_iteration; ++e) {
            model_getmemberfname(m, das->ensdir, varname, e + 1, fname_dst);
            if (enkf_separateout) {
                if (das->target == TARGET_ANALYSIS)
                    strncat(fname_dst, ".analysis", MAXSTRLEN);
                else if (das->target == TARGET_INCREMENT)
                    strncat(fname_dst, ".increment", MAXSTRLEN);
            } else {
                if (das->target == TARGET_ANALYSIS)
                    strncat(varname_dst, "_an", NC_MAX_NAME);
                else if (das->target == TARGET_INCREMENT)
                    strncat(varname_dst, "_inc", NC_MAX_NAME);
            }

            for (k = 0; k < nlev; ++k) {
                char fname_src[MAXSTRLEN];
                int ncid_src, vid_src;
                size_t start[3] = { e, 0, 0 };
                size_t count[3] = { 1, nj, ni };

                getfieldfname(das->ensdir, "ens", varname, k, fname_src);
                ncw_open(fname_src, NC_NOWRITE, &ncid_src);
                ncw_inq_varid(fname_src, ncid_src, varname, &vid_src);
                ncw_get_vara_float(fname_src, ncid_src, vid_src, start, count, v[0]);
                ncw_close(fname_src, ncid_src);

                model_writefield(m, fname_dst, MAXINT, varname_dst, k, v[0]);
            }
            enkf_printf(".");
        }
        enkf_printf("\n");
    }

    free2d(v);

#if defined(MPI)
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    /*
     * remove tiles 
     */
    if (!leavetiles && rank == 0) {
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
static void das_assemblebg(dasystem* das, int leavetiles)
{
    model* m = das->m;
    int nvar = model_getnvar(m);
    float** v = NULL;
    int ni, nj;
    int i;

#if defined(MPI)
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    if (rank > 0)
        return;

    model_getdims(m, &ni, &nj, NULL);
    v = alloc2d(ni, nj, sizeof(float));

    for (i = 0; i < nvar; ++i) {
        char* varname = model_getvarname(m, i);
        char varname_dst[NC_MAX_NAME];
        char fname_dst[MAXSTRLEN];
        int nlev, k;

        enkf_printf("    %s:", varname);

        model_getbgfname(m, das->bgdir, varname, fname_dst);
        nlev = getnlevels(fname_dst, varname);
        strncpy(varname_dst, varname, NC_MAX_NAME);
        if (enkf_separateout) {
            if (das->target == TARGET_ANALYSIS)
                strncat(fname_dst, ".analysis", MAXSTRLEN);
            else if (das->target == TARGET_INCREMENT)
                strncat(fname_dst, ".increment", MAXSTRLEN);
        } else {
            if (das->target == TARGET_ANALYSIS)
                strncat(varname_dst, "_an", NC_MAX_NAME);
            else if (das->target == TARGET_INCREMENT)
                strncat(varname_dst, "_inc", NC_MAX_NAME);
        }

        for (k = 0; k < nlev; ++k) {
            char fname_src[MAXSTRLEN];
            int ncid_src, vid_src;

            getfieldfname(das->bgdir, "bg", varname, k, fname_src);
            ncw_open(fname_src, NC_NOWRITE, &ncid_src);
            ncw_inq_varid(fname_src, ncid_src, varname, &vid_src);
            ncw_get_var_float(fname_src, ncid_src, vid_src, v[0]);
            ncw_close(fname_src, ncid_src);

            model_writefield(m, fname_dst, MAXINT, varname_dst, k, v[0]);
            if (!leavetiles)
                file_delete(fname_src);

            enkf_printf(".");
        }
        enkf_printf("\n");
    }

    free2d(v);
}

/**
 */
static void das_assemblespread(dasystem* das)
{
    model* m = das->m;
    int nvar = model_getnvar(m);
    float** v = NULL;
    int ni, nj;
    int i;

    model_getdims(m, &ni, &nj, NULL);
    v = alloc2d(ni, nj, sizeof(float));

    for (i = 0; i < nvar; ++i) {
        char* varname = model_getvarname(m, i);
        char varname_an[NC_MAX_NAME];
        int nlev, k;

        enkf_printf("    %s:", varname);

        nlev = getnlevels(FNAME_SPREAD, varname);

        if (das->mode == MODE_ENKF) {
            strncpy(varname_an, varname, NC_MAX_NAME);
            strncat(varname_an, "_an", NC_MAX_NAME);
        }

        for (k = 0; k < nlev; ++k) {
            char fname_src[MAXSTRLEN];
            int ncid_src, vid;

            getfieldfname(das->mode == MODE_ENKF ? das->ensdir : das->bgdir, "spread", varname, k, fname_src);
            ncw_open(fname_src, NC_NOWRITE, &ncid_src);

            ncw_inq_varid(fname_src, ncid_src, varname, &vid);
            ncw_get_var_float(fname_src, ncid_src, vid, v[0]);
            ncw_close(fname_src, ncid_src);
            file_delete(fname_src);

            model_writefield(m, FNAME_SPREAD, MAXINT, varname, k, v[0]);

            if (das->mode == MODE_ENKF) {
                getfieldfname(das->ensdir, "spread", varname_an, k, fname_src);
                ncw_open(fname_src, NC_NOWRITE, &ncid_src);
                ncw_inq_varid(fname_src, ncid_src, varname_an, &vid);
                ncw_get_var_float(fname_src, ncid_src, vid, v[0]);
                ncw_close(fname_src, ncid_src);
                file_delete(fname_src);

                model_writefield(m, FNAME_SPREAD, MAXINT, varname_an, k, v[0]);
            }

            enkf_printf(".");
        }
        enkf_printf("\n");
    }

    free2d(v);
}

/** Updates ensemble/background by using calculated transform 
 * matrices/coefficients.
 */
void das_update(dasystem* das, int calcspread, int leavetiles)
{
    model* m = das->m;
    int nvar = model_getnvar(m);
    void** fieldbuffer = NULL;
    int mni, mnj;
    int i, e;

    model_getdims(m, &mni, &mnj, NULL);

    if (das->nmem <= 0)
        das_getnmem(das);
    enkf_printf("    %d members\n", das->nmem);
    das_getfields(das);
    enkf_printf("    %d fields\n", das->nfields);

    if (calcspread && rank == 0) {
        enkf_printf("    allocating disk space for spread:");
        das_allocatespread(das, FNAME_SPREAD);
        enkf_printf("\n");
    }

    if (das->nplogs > 0 && rank == 0) {
        enkf_printf("    defining state variables in point logs:");
        plog_definestatevars(das);
        enkf_printf("\n");
    }

    if (das->mode == MODE_ENKF) {
        distribute_iterations(0, das->nmem - 1, nprocesses, rank);
#if defined(MPI)
        MPI_Barrier(MPI_COMM_WORLD);
#endif
        enkf_printtime("    ");
        enkf_printf("    allocating disk space for analysis:");
        enkf_flush();
#if defined(MPI)
        MPI_Barrier(MPI_COMM_WORLD);
#endif
        if (!enkf_separateout) {
            for (i = 0; i < nvar; ++i) {
                for (e = my_first_iteration; e <= my_last_iteration; ++e) {
                    char* varname_src = model_getvarname(m, i);
                    char fname[MAXSTRLEN];
                    int ncid;
                    char varname_dst[NC_MAX_NAME];

                    strncpy(varname_dst, varname_src, NC_MAX_NAME);
                    if (das->target == TARGET_ANALYSIS)
                        strncat(varname_dst, "_an", NC_MAX_NAME);
                    else if (das->target == TARGET_INCREMENT)
                        strncat(varname_dst, "_inc", NC_MAX_NAME);

                    model_getmemberfname(m, das->ensdir, varname_src, e + 1, fname);
                    ncw_open(fname, NC_WRITE, &ncid);
                    if (!ncw_var_exists(ncid, varname_dst)) {
                        ncw_redef(fname, ncid);
                        ncw_def_var_as(fname, ncid, varname_src, varname_dst);
                    }
                    ncw_close(fname, ncid);
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
                    if (das->target == TARGET_ANALYSIS)
                        strncat(fname_a, ".analysis", MAXSTRLEN);
                    else if (das->target == TARGET_INCREMENT)
                        strncat(fname_a, ".increment", MAXSTRLEN);
                    if (file_exists(fname_a)) {
                        ncw_open(fname_a, NC_WRITE, &ncid_a);
                        if (ncw_var_exists(ncid_a, varname)) {
                            ncw_close(fname_a, ncid_a);
                            ncw_close(fname_f, ncid_f);
                            continue;
                        }
                        ncw_redef(fname_a, ncid_a);
                    } else {
                        ncw_create(fname_a, NC_CLOBBER | NC_64BIT_OFFSET, &ncid_a);
                        ncw_copy_dims(fname_f, ncid_f, fname_a, ncid_a);
                    }

                    ncw_inq_varid(fname_f, ncid_f, varname, &vid_f);
                    ncw_copy_vardef(fname_f, ncid_f, vid_f, fname_a, ncid_a);
                    ncw_close(fname_a, ncid_a);
                    ncw_close(fname_f, ncid_f);
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

            if (!enkf_separateout) {
                for (i = 0; i < nvar; ++i) {
                    char* varname_src = model_getvarname(m, i);
                    char fname[MAXSTRLEN];
                    int ncid;
                    char varname_dst[NC_MAX_NAME];

                    strncpy(varname_dst, varname_src, NC_MAX_NAME);
                    if (das->target == TARGET_ANALYSIS)
                        strncat(varname_dst, "_an", NC_MAX_NAME);
                    else if (das->target == TARGET_INCREMENT)
                        strncat(varname_dst, "_inc", NC_MAX_NAME);
                    model_getbgfname(m, das->bgdir, varname_src, fname);
                    ncw_open(fname, NC_WRITE, &ncid);
                    if (!ncw_var_exists(ncid, varname_dst)) {
                        ncw_redef(fname, ncid);
                        ncw_def_var_as(fname, ncid, varname_src, varname_dst);
                    }
                    ncw_close(fname, ncid);
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
                    if (das->target == TARGET_ANALYSIS)
                        strncat(fname_a, ".analysis", MAXSTRLEN);
                    else if (das->target == TARGET_INCREMENT)
                        strncat(fname_a, ".increment", MAXSTRLEN);

                    if (file_exists(fname_a)) {
                        ncw_open(fname_a, NC_WRITE, &ncid_a);
                        if (ncw_var_exists(ncid_a, varname)) {
                            ncw_close(fname_a, ncid_a);
                            ncw_close(fname_f, ncid_f);
                            continue;
                        }
                        ncw_redef(fname_a, ncid_a);
                    } else {
                        ncw_create(fname_a, NC_CLOBBER | NC_64BIT_OFFSET, &ncid_a);
                        ncw_copy_dims(fname_f, ncid_f, fname_a, ncid_a);
                    }

                    ncw_inq_varid(fname_f, ncid_f, varname, &vid_f);
                    ncw_copy_vardef(fname_f, ncid_f, vid_f, fname_a, ncid_a);
                    ncw_close(fname_a, ncid_a);
                    ncw_close(fname_f, ncid_f);
                    printf(".");
                    fflush(stdout);
                }
            }

            enkf_printf("\n");
            enkf_printtime("    ");
            enkf_flush();
        }
    }

    distribute_iterations(0, das->nfields - 1, nprocesses, rank);

    fieldbuffer = malloc(das->fieldbufsize * sizeof(void*));
    if (das->mode == MODE_ENKF) {
        for (i = 0; i < das->fieldbufsize; ++i)
            fieldbuffer[i] = alloc3d(das->nmem, mnj, mni, sizeof(float));
    } else if (das->mode == MODE_ENOI) {
        for (i = 0; i < das->fieldbufsize; ++i)
            fieldbuffer[i] = alloc3d(das->nmem + 1, mnj, mni, sizeof(float));
    }

    enkf_printf("    processing fields:\n");
    enkf_flush();
#if defined(MPI)
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    for (i = my_first_iteration; i <= my_last_iteration; ++i) {
        int bufindex = (i - my_first_iteration) % das->fieldbufsize;
        field* f = &das->fields[i];
        char fname[MAXSTRLEN];

        printf("      %-8s %-3d (%d: %d: %.1f%%)\n", f->varname, f->level, rank, i, 100.0 * (double) (i - my_first_iteration + 1) / (double) (my_last_iteration - my_first_iteration + 1));
        fflush(stdout);

        for (e = 0; e < das->nmem; ++e) {
            model_getmemberfname(m, das->ensdir, f->varname, e + 1, fname);
            model_readfield(das->m, fname, e + 1, MAXINT, f->varname, f->level, ((float***) fieldbuffer[bufindex])[e][0]);
        }
        if (das->mode == MODE_ENOI) {
            if (das->target == TARGET_ANALYSIS) {
                model_getbgfname(m, das->bgdir, f->varname, fname);
                model_readfield(das->m, fname, MAXINT, MAXINT, f->varname, f->level, ((float***) fieldbuffer[bufindex])[das->nmem][0]);
            } else if (das->target == TARGET_INCREMENT)
                memset(((float***) fieldbuffer[bufindex])[das->nmem][0], 0, mni * mnj * sizeof(float));
        }

        if (bufindex == das->fieldbufsize - 1 || i == my_last_iteration) {
            if (calcspread)     /* write forecast spread */
                das_writespread(das, bufindex + 1, fieldbuffer, &das->fields[i - bufindex], 0);

            if (das->mode == MODE_ENKF) {
                das_updatefields(das, bufindex + 1, fieldbuffer, &das->fields[i - bufindex]);
                das_writefields(das, bufindex + 1, fieldbuffer, &das->fields[i - bufindex]);
            } else if (das->mode == MODE_ENOI) {
                das_updatebg(das, bufindex + 1, fieldbuffer);
                das_writebg(das, bufindex + 1, fieldbuffer, &das->fields[i - bufindex]);
            }

            if (calcspread && das->mode == MODE_ENKF)   /* write analysis * * 
                                                         * spread */
                das_writespread(das, bufindex + 1, fieldbuffer, &das->fields[i - bufindex], 1);
        }
    }
    enkf_printf("\n");
    enkf_flush();

    if (!enkf_directwrite) {
        enkf_printtime("  ");
        enkf_printf("  assembling analysis:\n");
        if (das->mode == MODE_ENKF)
            das_assemblemembers(das, leavetiles);
        else if (das->mode == MODE_ENOI)
            das_assemblebg(das, leavetiles);
        if (calcspread && rank == 0) {
            enkf_printf("  assembling spread:\n");
            das_assemblespread(das);
        }
        if (das->nplogs > 0) {
            enkf_printf("  assembling state variables in point logs:\n");
            plog_assemblestatevars(das);
        }
    }

    for (i = 0; i < das->fieldbufsize; ++i)
        free3d(fieldbuffer[i]);
    free(fieldbuffer);
}
