/******************************************************************************
 *
 * File:        model2obs.c        
 *
 * Created:     12/2012
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description: H functions to calculate forecast observations from the model
 *              state, basically -- 2D and 3D linear interpolators.
 *
 * Revisions:
 *
 *****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <float.h>
#include <limits.h>
#include "ncw.h"
#include "definitions.h"
#include "utils.h"
#include "ncutils.h"
#include "observations.h"
#include "model.h"
#include "model2obs.h"

#define EPSMULT (float) 1.000001
#define EPS_IJ 0.01
#define FOOTPRINT_INC 1000

/**
 */
static void evaluate_2d_obs(model* m, observations* allobs, int nobs, int obsids[], char fname[], float** v, float out[])
{
    int otid = allobs->data[obsids[0]].type;
    int mvid = model_getvarid(m, allobs->obstypes[otid].varnames[0], 1);
    int periodic_i = grid_isperiodic_i(model_getvargrid(m, mvid));
    int** mask = model_getnumlevels(m, mvid);
    int ni, nj;
    int i;

    model_getvargridsize(m, mvid, &ni, &nj, NULL);
    for (i = 0; i < nobs; ++i) {
        int ii = obsids[i];
        observation* o = &allobs->data[ii];

        assert(o->type == otid);
        assert(isnan(out[ii]));
        if (o->footprint == 0.0)
            out[ii] = interpolate2d(o->fi, o->fj, ni, nj, v, mask, periodic_i);
        else {
            grid* g = model_getvargrid(m, mvid);
            kdtree* tree = grid_gettreeXYZ(g, 1);
            double ll[2] = { o->lon, o->lat };
            double xyz[3];
            kdset* set = NULL;
            int ncells = 0;
            size_t* ids = NULL;
            size_t id;

            ll2xyz(ll, xyz);
            set = kd_findnodeswithinrange(tree, xyz, o->footprint, 0);
            if (kdset_getsize(set) > 0) {
                while ((id = kdset_readnext(set, NULL)) != SIZE_MAX) {
                    int id_orig = kd_getnodedata(tree, id);

                    if (ncells % FOOTPRINT_INC == 0)
                        ids = realloc(ids, (ncells + FOOTPRINT_INC) * sizeof(size_t));
                    ids[ncells] = id_orig;
                    ncells++;
                }
                kdset_free(set);
                out[ii] = average2d(ncells, ids, v);
                free(ids);
            } else {
                obstype* ot = &allobs->obstypes[o->type];

                enkf_printf("\n  obs # %d: ", ii);
                obs_printob(allobs, ii);
                enkf_quit("obs # %d: it seems that the grid is too coarse to handle footprint = %f for observations of type \"%s\"; you need to set the footprint to zero (or remove the entry)", ii, o->footprint, ot->name);
            }
        }
        if (!isfinite(out[ii]) || fabs(out[ii]) > STATE_BIGNUM) {
            enkf_flush();
            enkf_printf("\n  obs # %d: ", ii);
            obs_printob(allobs, ii);
            enkf_quit("obs # %d: forecast = %.3g for \"%s\"; no point to continue", ii, out[ii], fname);
        }
    }
}

/**
 */
static void interpolate_3d_obs(model* m, observations* allobs, int nobs, int obsids[], char fname[], float*** v, float out[])
{
    int otid = allobs->data[obsids[0]].type;
    int mvid = model_getvarid(m, allobs->obstypes[otid].varnames[0], 1);
    void* grid = model_getvargrid(m, mvid);
    int ksurf = grid_getsurflayerid(grid);
    int periodic_i = grid_isperiodic_i(grid);
    int** nlevels = model_getnumlevels(m, mvid);
    int ni, nj, nk;
    int i;

    model_getvargridsize(m, mvid, &ni, &nj, &nk);
    for (i = 0; i < nobs; ++i) {
        int ii = obsids[i];
        observation* o = &allobs->data[ii];

        assert(o->type == otid);
        assert(isnan(out[ii]));
        out[ii] = interpolate3d(o->fi, o->fj, o->fk, ni, nj, nk, ksurf, v, nlevels, periodic_i);
        if (!isfinite(out[ii]) || fabs(out[ii]) > STATE_BIGNUM) {
            enkf_flush();
            enkf_printf("\n  obs # %d: ", ii);
            obs_printob(allobs, ii);
            enkf_quit("obs # %d: forecast = %.3g in \"%s\"; no point to continue", ii, out[ii], fname);
        }
    }
}

/**
 */
void H_surf_standard(dasystem* das, int nobs, int obsids[], char fname[], int mem, int t, float dst[])
{
    model* m = das->m;
    observations* allobs = das->obs;
    int otid = allobs->data[obsids[0]].type;
    obstype* ot = &allobs->obstypes[otid];
    int mvid = model_getvarid(m, ot->varnames[0], 1);
    int ksurf = grid_getsurflayerid(model_getvargrid(m, mvid));
    int ni, nj;
    float** src = NULL;
    char tag_offset[MAXSTRLEN];
    float** offset = NULL;

    assert(ot->nvar == 1);      /* should we care? */

    model_getvargridsize(m, mvid, &ni, &nj, NULL);
    src = alloc2d(nj, ni, sizeof(float));
    model_readfield(m, fname, ot->varnames[0], ksurf, src[0]);

    snprintf(tag_offset, MAXSTRLEN, "%s:OFFSET", ot->name);
    offset = model_getdata(m, tag_offset);
    if (offset != NULL) {
        int ni, nj;
        float* src0 = src[0];
        float* offset0 = offset[0];
        int i;

        assert(model_getdataalloctype(m, tag_offset) == ALLOCTYPE_2D);
        assert(mvid >= 0);
        model_getvargridsize(m, mvid, &ni, &nj, NULL);
        for (i = 0; i < ni * nj; ++i)
            src0[i] -= offset0[i];
    }

    evaluate_2d_obs(m, allobs, nobs, obsids, fname, src, dst);

    free(src);
}

/**
 */
void H_surf_biased(dasystem* das, int nobs, int obsids[], char fname[], int mem, int t, float dst[])
{
    model* m = das->m;
    observations* allobs = das->obs;
    int otid = allobs->data[obsids[0]].type;
    obstype* ot = &allobs->obstypes[otid];
    int mvid = model_getvarid(m, ot->varnames[0], 1);
    int ksurf = grid_getsurflayerid(model_getvargrid(m, mvid));
    int ni, nj, nv;
    float** src = NULL;
    float* src0 = NULL;
    char tag_offset[MAXSTRLEN];
    float** offset = NULL;
    float* bias = NULL;
    char fname2[MAXSTRLEN];
    int i;

    if (ot->nvar < 2)
        enkf_quit("%s: second variable has to be defined for the observation type when using H-function \"biased\"", ot->name);
    {
        int mvid2 = model_getvarid(m, ot->varnames[1], 1);

        if (model_getvargridid(m, mvid) != model_getvargridid(m, mvid2))
            enkf_quit("H_surf_biased(): variables \"%s\" and \"%s\" are defined on different grids", ot->varnames[0], ot->varnames[1]);
    }

    model_getvargridsize(m, mvid, &ni, &nj, NULL);
    nv = ni * nj;
    src = alloc2d(nj, ni, sizeof(float));
    src0 = src[0];

    bias = malloc(nv * sizeof(float));
    if (das->mode == MODE_ENKF)
        das_getmemberfname(das, das->ensdir, ot->varnames[1], mem, fname2);
    else if (das->mode == MODE_ENOI)
        das_getbgfname(das, das->bgdir, ot->varnames[1], fname2);
    model_readfield(m, fname2, ot->varnames[1], ksurf, bias);

    model_readfield(m, fname, ot->varnames[0], ksurf, src0);

    snprintf(tag_offset, MAXSTRLEN, "%s:OFFSET", allobs->obstypes[allobs->data[obsids[0]].type].name);
    offset = model_getdata(m, tag_offset);

    if (offset != NULL) {
        float* offset0 = offset[0];

        assert(model_getdataalloctype(m, tag_offset) == ALLOCTYPE_2D);
        for (i = 0; i < nv; ++i)
            src0[i] -= offset0[i] + bias[i];
    } else
        for (i = 0; i < nv; ++i)
            src0[i] -= bias[i];

    evaluate_2d_obs(m, allobs, nobs, obsids, fname, src, dst);

    free(bias);
    free(src);
}

/**
 */
void H_subsurf_standard(dasystem* das, int nobs, int obsids[], char fname[], int mem, int t, float dst[])
{
    model* m = das->m;
    observations* allobs = das->obs;
    int otid = allobs->data[obsids[0]].type;
    obstype* ot = &allobs->obstypes[otid];
    int mvid = model_getvarid(m, ot->varnames[0], 1);
    int ni, nj, nk;
    float*** src = NULL;
    char tag_offset[MAXSTRLEN];
    void* offset_data = NULL;

    assert(ot->nvar == 1);      /* should we care? */
    model_getvargridsize(m, mvid, &ni, &nj, &nk);
    src = alloc3d(nk, nj, ni, sizeof(float));
    if (nk > 1)
        model_read3dfield(m, fname, ot->varnames[0], src[0][0]);
    else if (nk == 1)
        model_readfield(m, fname, ot->varnames[0], 0, src[0][0]);

    snprintf(tag_offset, MAXSTRLEN, "%s:OFFSET", allobs->obstypes[allobs->data[obsids[0]].type].name);
    offset_data = model_getdata(m, tag_offset);
    if (offset_data != NULL) {
        if (model_getdataalloctype(m, tag_offset) == ALLOCTYPE_3D) {
            float* src0 = src[0][0];
            float* offset0 = ((float***) offset_data)[0][0];
            size_t ijk;

            for (ijk = 0; ijk < (size_t) ni * nj * nk; ++ijk)
                src0[ijk] -= offset0[ijk];
        } else if (model_getdataalloctype(m, tag_offset) == ALLOCTYPE_1D) {
            size_t ij, k;

            for (k = 0; k < nk; ++k) {
                float* srck = src[k][0];
                float offsetk = ((float*) offset_data)[k];

                for (ij = 0; ij < ni * nj; ++ij)
                    srck[ij] -= offsetk;
            }
        } else
            enkf_quit("obstype = %s: offset variable must be either 3D or 1D for a 3D observation type", ot->name);
    }

    interpolate_3d_obs(m, allobs, nobs, obsids, fname, src, dst);
    free(src);
}

/**
 */
static int cmp_obs_byfk(const void* p1, const void* p2, void* p)
{
    observation* obs = (observation*) p;
    observation* o1 = &obs[*(int*) p1];
    observation* o2 = &obs[*(int*) p2];

    if (o1->fk > o2->fk)
        return 1;
    if (o1->fk < o2->fk)
        return -1;
    return 0;
}

/** This is a subsurface (3D) interpolator that avoids reading the whole 3D
 ** field. Instead, it proceeds by keeping in memory two layers at a time and
 ** interpolating observations in between these layers.
 */
void H_subsurf_lowmem(dasystem* das, int nobs, int obsids[], char fname[], int mem, int t, float dst[])
{
    model* m = das->m;
    observations* allobs = das->obs;
    int otid = allobs->data[obsids[0]].type;
    obstype* ot = &allobs->obstypes[otid];
    int mvid = model_getvarid(m, ot->varnames[0], 1);
    int ni, nj, nk;
    float** src1 = NULL;
    float** src2 = NULL;
    float*** src;
    char tag_offset[MAXSTRLEN];
    void* offset_data = NULL;
    int k1, i1, i2, k1_now, k2_now;;

    assert(ot->nvar == 1);      /* should we care? */
    model_getvargridsize(m, mvid, &ni, &nj, &nk);

    if (nk == 1) {
        H_subsurf_standard(das, nobs, obsids, fname, mem, t, dst);
        return;
    }

    snprintf(tag_offset, MAXSTRLEN, "%s:OFFSET", allobs->obstypes[allobs->data[obsids[0]].type].name);
    offset_data = model_getdata(m, tag_offset);
    if (offset_data != NULL) {
        if (model_getdataalloctype(m, tag_offset) == ALLOCTYPE_3D)
            enkf_quit("obstype = %s: \"lowmem\" subsurface H-function works with 1D offsets only; use \"standard\" function for 3D offset", ot->name);
        else if (model_getdataalloctype(m, tag_offset) == ALLOCTYPE_2D)
            enkf_quit("obstype = %s:  \"lowmem\" subsurface H-function works with 1D offsets only", ot->name);
    }

    src = calloc(nk, sizeof(float**));
    src1 = alloc2d(nj, ni, sizeof(float));
    src2 = alloc2d(nj, ni, sizeof(float));

    /*
     * (we assume that obsids are not used elsewhere and thus can be reordered)
     */
    qsort_r(obsids, nobs, sizeof(int), cmp_obs_byfk, allobs->data);

    for (k1 = 0, i1 = 0, i2 = -1, k1_now = -1, k2_now = -1; k1 < nk - 1 && i2 < nobs; ++k1) {
        int k2 = k1 + 1;

        while (i2 + 1 < nobs && allobs->data[obsids[i2 + 1]].fk < (double) k2)
            i2++;

        if (i2 < i1)
            continue;

        if (k2_now == k1) {
            void* tmp = src1;

            src1 = src2;
            src2 = tmp;
            k1_now = k1;
        }
        if (k1_now != k1) {
            model_readfield(m, fname, ot->varnames[0], k1, src1[0]);
            k1_now = k1;
        }
        model_readfield(m, fname, ot->varnames[0], k2, src2[0]);
        k2_now = k2;

        /*
         * (interpolate)
         */
        src[k1] = src1;
        src[k2] = src2;
        if (offset_data != NULL) {
            int k;

            for (k = k1; k <= k2; ++k) {
                float* srck = src[k][0];
                float offsetk = ((float*) offset_data)[k];
                size_t ij;

                for (ij = 0; ij < (size_t) ni * nj; ++ij)
                    srck[ij] -= offsetk;
            }
        }
        interpolate_3d_obs(m, allobs, i2 - i1 + 1, &obsids[i1], fname, src, dst);

        i1 = i2 + 1;
    }
    free(src1);
    free(src2);
    free(src);
}

#define MLD_TRANSITION 0.1

/**
 */
static double mldtaper(double mld, double z)
{
    double v;

    v = (z / mld - 1) / MLD_TRANSITION;
    if (v < -1.0)
        return 1.0;
    else if (v > 1.0)
        return 0.0;

    return cos((v + 1) * M_PI / 2.0) * 0.5 + 0.5;
}

/** Projects surface bias into subsurface based on the mixed layer depth.
 */
void H_subsurf_wsurfbias(dasystem* das, int nobs, int obsids[], char fname[], int mem, int t, float dst[])
{
    model* m = das->m;
    observations* allobs = das->obs;
    int otid = allobs->data[obsids[0]].type;
    obstype* ot = &allobs->obstypes[otid];
    int mvid = model_getvarid(m, ot->varnames[0], 1);
    int ni, nj, nk;
    float*** src = NULL;
    int periodic_i = grid_isperiodic_i(model_getvargrid(m, mvid));
    int** mask = model_getnumlevels(m, mvid);

    float** mld = NULL;
    char tag_offset[MAXSTRLEN];
    void* offset_data = NULL;
    float** bias = NULL;
    char fname2[MAXSTRLEN];
    size_t i;

    if (ot->nvar < 2)
        enkf_quit("%s: second variable has to be defined for the observation type when using H-function \"wsurfbias\"", ot->name);
    {
        int mvid2 = model_getvarid(m, ot->varnames[1], 1);

        if (model_getvargridid(m, mvid) != model_getvargridid(m, mvid2))
            enkf_quit("H_surf_biased(): variables \"%s\" and \"%s\" are defined on different grids", ot->varnames[0], ot->varnames[1]);
    }

    model_getvargridsize(m, mvid, &ni, &nj, &nk);
    src = alloc3d(nk, nj, ni, sizeof(float));

    /*
     * this part is similar to H_subsurf_standard()
     */
    if (nk > 1)
        model_read3dfield(m, fname, ot->varnames[0], src[0][0]);
    else if (nk == 1)
        model_readfield(m, fname, ot->varnames[0], 0, src[0][0]);

    snprintf(tag_offset, MAXSTRLEN, "%s:OFFSET", ot->name);
    offset_data = model_getdata(m, tag_offset);
    if (offset_data != NULL) {
        if (model_getdataalloctype(m, tag_offset) == ALLOCTYPE_3D) {
            float* src0 = src[0][0];
            float* offset0 = ((float***) offset_data)[0][0];
            size_t ijk;

            for (ijk = 0; ijk < (size_t) ni * nj * nk; ++ijk)
                src0[ijk] -= offset0[ijk];
        } else if (model_getdataalloctype(m, tag_offset) == ALLOCTYPE_1D) {
            size_t ij, k;

            for (k = 0; k < nk; ++k) {
                float* srck = src[k][0];
                float offsetk = ((float*) offset_data)[k];

                for (ij = 0; ij < ni * nj; ++ij)
                    srck[ij] -= offsetk;
            }
        } else
            enkf_quit("obstype = %s: offset variable must be either 3D or 1D for a 3D observation type", ot->name);
    }

    interpolate_3d_obs(m, allobs, nobs, obsids, fname, src, dst);

    /*
     * now correct for the surface bias
     */
    bias = alloc2d(nj, ni, sizeof(float));
    if (das->mode == MODE_ENKF)
        das_getmemberfname(das, das->ensdir, ot->varnames[1], mem, fname2);
    else if (das->mode == MODE_ENOI)
        das_getbgfname(das, das->bgdir, ot->varnames[1], fname2);
    assert(ncu_getnD(fname2, ot->varnames[1]) == 2);
    model_readfield(m, fname2, ot->varnames[1], 0, bias[0]);

    if (isnan(ot->mld_threshold) && ot->mld_varname == NULL)
        enkf_quit("\"MLD_THRESH\" or \"MLD_VARNAME\" must be specified for observation type \"%s\" to use H function \"wsurfbias\"", ot->name);
    if (!isnan(ot->mld_threshold) && ot->mld_varname != NULL)
        enkf_quit("both \"MLD_THRESH\" and \"MLD_VARNAME\" are specified for observation type \"%s\"", ot->name);
    mld = alloc2d(nj, ni, sizeof(float));
    if (ot->mld_varname != NULL) {
        char fname_mld[MAXSTRLEN];

        if (model_getvarid(m, ot->mld_varname, 0) < 0)
            enkf_quit("\"MLD_VARNAME = %s\" for observation type \"%s\" does not exist among model variables", ot->mld_varname, ot->name);
        if (das->mode == MODE_ENKF) {
            das_getmemberfname(das, das->ensdir, ot->mld_varname, mem, fname_mld);
            model_readfield(m, fname_mld, ot->mld_varname, 0, mld[0]);
        } else if (das->mode == MODE_ENOI) {
            das_getbgfname(das, das->bgdir, ot->mld_varname, fname_mld);
            model_readfield(m, fname_mld, ot->mld_varname, 0, mld[0]);
        }
    } else {
        if (das->mode == MODE_ENKF)
            das_calcmld(das, ot, src, mld);
        else if (das->mode == MODE_ENOI) {
            char tag[MAXSTRLEN];

            snprintf(tag, MAXSTRLEN, "%s:MLD", ot->name);
            if (mem <= 0) {     /* background */
                das_calcmld(das, ot, src, mld);
                model_addorreplacedata(m, tag, mvid, ALLOCTYPE_2D, mld);
            } else              /* members */
                mld = model_getdata(m, tag);
        } else
            enkf_quit("programming error");
    }

    {
        double fi_prev = DBL_MAX;
        double fj_prev = DBL_MAX;
        double vmld = NAN, vbias = NAN;

        for (i = 0; i < nobs; ++i) {
            size_t ii = obsids[i];
            observation* o = &allobs->data[ii];

            if (fabs(fi_prev - o->fi) > EPS_IJ || fabs(fj_prev - o->fj) > EPS_IJ) {
                vmld = interpolate2d(o->fi, o->fj, ni, nj, mld, mask, periodic_i);
                vbias = interpolate2d(o->fi, o->fj, ni, nj, bias, mask, periodic_i);
                fi_prev = o->fi;
                fj_prev = o->fj;
            }

            if (!isfinite(vmld))
                continue;

            assert(isfinite(vbias));
            dst[ii] -= mldtaper(vmld, o->depth) * vbias;
        }
    }

    if (das->mode == MODE_ENKF)
        free(mld);
    free(bias);
    free(src);
}

/** Sum up the estimates over layers.
 */
void H_vertsum(dasystem* das, int nobs, int obsids[], char fname[], int mem, int t, float dst[])
{
    model* m = das->m;
    observations* allobs = das->obs;
    int otid = allobs->data[obsids[0]].type;
    obstype* ot = &allobs->obstypes[otid];
    int mvid = model_getvarid(m, ot->varnames[0], 1);
    int ni, nj, nk;
    float** src = NULL;
    float** srcsum = NULL;
    float* src0 = NULL;
    float* srcsum0 = NULL;
    char tag_offset[MAXSTRLEN];
    void* offset_data = NULL;
    size_t k, i;

    {
        int vtype = grid_getvtype(model_getvargrid(m, mvid));

        if (vtype != GRIDVTYPE_SIGMA && vtype != GRIDVTYPE_HYBRID)
            enkf_quit("obstype = %s: H-function \"vertsum\" can only be used for variables on either sigma or hybrid vertical grids", ot->name);
    }

    assert(ot->nvar == 1);      /* should we care? */
    model_getvargridsize(m, mvid, &ni, &nj, &nk);
    src = alloc2d(nj, ni, sizeof(float));
    srcsum = alloc2d(nj, ni, sizeof(float));
    src0 = src[0];
    srcsum0 = srcsum[0];

    model_readfield(m, fname, ot->varnames[0], 0, srcsum0);
    for (k = 1; k < nk; ++k) {
        model_readfield(m, fname, ot->varnames[0], k, src0);
        for (i = 0; i < (size_t) ni * nj; i++)
            srcsum0[i] += src0[i];
    }
    free(src);
    src0 = NULL;

    snprintf(tag_offset, MAXSTRLEN, "%s:OFFSET", allobs->obstypes[allobs->data[obsids[0]].type].name);
    offset_data = model_getdata(m, tag_offset);
    if (offset_data != NULL) {
        float* offset0 = ((float**) offset_data)[0];

        if (model_getdataalloctype(m, tag_offset) != ALLOCTYPE_2D)
            enkf_quit("obstype = %s: offset variable must be 2D to be used in H-function \"vertsum\"", ot->name);

        for (i = 0; i < (size_t) ni * nj; ++i)
            srcsum0[i] -= offset0[i];
    }

    evaluate_2d_obs(m, allobs, nobs, obsids, fname, srcsum, dst);
    free(srcsum);
}

/** Calculate weighted average of estimates over layers.
 */
void H_vertwavg(dasystem* das, int nobs, int obsids[], char fname[], int mem, int t, float dst[])
{
    model* m = das->m;
    observations* allobs = das->obs;
    int otid = allobs->data[obsids[0]].type;
    obstype* ot = &allobs->obstypes[otid];
    int mvid = model_getvarid(m, ot->varnames[0], 1);
    int ni, nj, nk;
    float* src = NULL;
    float** sum = NULL;
    float* w = NULL;
    float* sumw = NULL;
    char tag_offset[MAXSTRLEN];
    void* offset_data = NULL;
    size_t k, i;

    if (ot->nvar < 2)
        enkf_quit("%s: second variable has to be defined for the observation type when using H-function \"vertavg\"", ot->name);
    {
        int mvid2 = model_getvarid(m, ot->varnames[1], 1);

        if (model_getvargridid(m, mvid) != model_getvargridid(m, mvid2))
            enkf_quit("H_vertwavg(): variables \"%s\" and \"%s\" are defined on different grids", ot->varnames[0], ot->varnames[1]);
    }

    model_getvargridsize(m, mvid, &ni, &nj, &nk);
    src = calloc(nj * ni, sizeof(float));
    sum = alloc2d(nj, ni, sizeof(float));
    w = calloc(nj * ni, sizeof(float));
    sumw = calloc(nj * ni, sizeof(float));

    /*
     * get sum of weights 
     */
    model_readfield(m, fname, ot->varnames[1], 0, sumw);
    for (k = 1; k < nk; ++k) {
        model_readfield(m, fname, ot->varnames[1], k, w);
        for (i = 0; i < (size_t) ni * nj; i++)
            sumw[i] += w[i];
    }

    /*
     * calculate weighted average 
     */
    for (k = 0; k < nk; ++k) {
        float* sum0 = sum[0];

        model_readfield(m, fname, ot->varnames[0], k, src);
        model_readfield(m, fname, ot->varnames[1], k, w);
        for (i = 0; i < (size_t) ni * nj; i++) {
            double v = src[i] * w[i] / sumw[i];

            if (isfinite(v) && v > 0.0)
                sum0[i] += v;
        }
    }

    free(src);
    free(w);
    free(sumw);

    snprintf(tag_offset, MAXSTRLEN, "%s:OFFSET", allobs->obstypes[allobs->data[obsids[0]].type].name);
    offset_data = model_getdata(m, tag_offset);
    if (offset_data != NULL) {
        float* offset0 = ((float**) offset_data)[0];
        float* sum0 = sum[0];

        if (model_getdataalloctype(m, tag_offset) != ALLOCTYPE_2D)
            enkf_quit("obstype = %s: offset variable must be 2D to be used in H-function \"vertavg\"", ot->name);

        for (i = 0; i < (size_t) ni * nj; ++i)
            sum0[i] -= offset0[i];
    }

    evaluate_2d_obs(m, allobs, nobs, obsids, fname, sum, dst);
    free(sum);
}
