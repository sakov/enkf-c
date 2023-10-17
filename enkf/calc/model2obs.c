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

/**
 */
static void evaluate_2d_obs(model* m, observations* allobs, int nobs, int obsids[], char fname[], void* v, float out[])
{
    int otid = allobs->data[obsids[0]].type;
    int mvid = model_getvarid(m, allobs->obstypes[otid].varnames[0], 1);
    int i;

    for (i = 0; i < nobs; ++i) {
        int ii = obsids[i];
        observation* o = &allobs->data[ii];
        grid* g = model_getvargrid(m, mvid);

        assert(o->type == otid);
        assert(isnan(out[ii]));
        if (o->footprint == 0.0)
            out[ii] = grid_interpolate2d(g, o->fij, v);
        else {
            kdtree* tree = grid_gettreeXYZ(g, 1);
            double ll[2] = { o->lon, o->lat };
            double xyz[3];
            size_t ncells = 0;
            kdresult* results = NULL;

            ll2xyz(ll, xyz);
            kd_findnodeswithinrange(tree, xyz, o->footprint, 0, &ncells, &results);
            if (ncells > 0) {
                size_t* ids = malloc(ncells * sizeof(size_t));
                size_t iloc;

                for (iloc = 0; iloc < ncells; ++iloc)
                    ids[iloc] = kd_getnodedata(tree, results[iloc].id);
                out[ii] = average(ncells, ids, grid_isstructured(g) ? ((float**) v)[0] : v);
                free(ids);
            } else {
                obstype* ot = &allobs->obstypes[o->type];

                enkf_printf("\n  obs # %d: ", ii);
                obs_printob(allobs, ii);
                enkf_quit("obs # %d: it seems that the grid is too coarse to handle footprint = %f for observations of type \"%s\"; you need to set the footprint to zero (or remove the entry)", ii, o->footprint, ot->name);
            }
        }
        if (!isfinite(out[ii]) || fabs(out[ii]) > STATE_BIGNUM) {
            if (!skip_bad_fc_obs) {
                enkf_flush();
                enkf_verbose = -1;      /* force printing regardless of rank */
                enkf_printf("\n  obs # %d: ", ii);
                obs_printob(allobs, ii);
                enkf_quit("obs # %d: forecast = %.3g for \"%s\"; no point to continue", ii, out[ii], fname);
            } else
                o->status = STATUS_BADFC;
        }
    }
}

/**
 */
static void interpolate_3d_obs(model* m, observations* allobs, int nobs, int obsids[], char fname[], float*** v, float out[])
{
    int otid = allobs->data[obsids[0]].type;
    int mvid = model_getvarid(m, allobs->obstypes[otid].varnames[0], 1);
    void* g = model_getvargrid(m, mvid);
    int i;

    for (i = 0; i < nobs; ++i) {
        int ii = obsids[i];
        observation* o = &allobs->data[ii];

        assert(o->type == otid);
        assert(isnan(out[ii]));
        out[ii] = grid_interpolate3d(g, o->fij, o->fk, v);

        if (!isfinite(out[ii]) || fabs(out[ii]) > STATE_BIGNUM) {
            if (!skip_bad_fc_obs) {
                enkf_flush();
                enkf_verbose = -1;      /* force printing regardless of rank */
                enkf_printf("\n  obs # %d: ", ii);
                obs_printob(allobs, ii);
                enkf_quit("obs # %d: forecast = %.3g in \"%s\"; no point to continue", ii, out[ii], fname);
            } else
                o->status = STATUS_BADFC;
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
    int masklog = das_isstatic(das, mem);
    int ni, nj;
    void* src = NULL;
    char tag_offset[MAXSTRLEN];

    assert(ot->nvar == 1);      /* should we care? */

    model_getvargridsize(m, mvid, &ni, &nj, NULL);
    if (nj > 0) {
        float** offset;

        src = alloc2d(nj, ni, sizeof(float));
        model_readfield(m, fname, ot->varnames[0], ksurf, ((float**) src)[0], masklog);

        snprintf(tag_offset, MAXSTRLEN, "%s:OFFSET", ot->name);
        offset = model_getdata(m, tag_offset);
        if (offset != NULL) {
            float* src0 = ((float**) src)[0];
            float* offset0 = offset[0];
            int i;

            assert(model_getdataalloctype(m, tag_offset) == ALLOCTYPE_2D);
            for (i = 0; i < ni * nj; ++i)
                src0[i] -= offset0[i];
        }
    } else {
        float* offset;

        src = calloc(ni, sizeof(float));
        model_readfield(m, fname, ot->varnames[0], ksurf, src, masklog);

        snprintf(tag_offset, MAXSTRLEN, "%s:OFFSET", ot->name);
        offset = model_getdata(m, tag_offset);
        if (offset != NULL) {
            int i;

            assert(model_getdataalloctype(m, tag_offset) == ALLOCTYPE_1D);
            for (i = 0; i < ni; ++i)
                ((float*) src)[i] -= offset[i];
        }
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
    int masklog = das_isstatic(das, mem);
    int ni, nj, nv;
    void* src = NULL;
    float* src0 = NULL;
    char tag_offset[MAXSTRLEN];
    void* offset = NULL;
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
    if (nj > 0) {
        nv = ni * nj;
        src = alloc2d(nj, ni, sizeof(float));
        src0 = ((float**) src)[0];
    } else {
        nv = ni;
        src = calloc(ni, sizeof(float));
        src0 = src;
    }

    model_readfield(m, fname, ot->varnames[0], ksurf, src0, masklog);

    snprintf(tag_offset, MAXSTRLEN, "%s:OFFSET", allobs->obstypes[allobs->data[obsids[0]].type].name);
    offset = model_getdata(m, tag_offset);
    if (offset != NULL) {
        float* offset0 = (nj > 0) ? ((float**) offset)[0] : offset;

        assert(model_getdataalloctype(m, tag_offset) == ALLOCTYPE_2D);
        for (i = 0; i < nv; ++i)
            src0[i] -= offset0[i];
    }

    if (das->mode != MODE_HYBRID || mem <= das->nmem_dynamic) {
        bias = malloc(nv * sizeof(float));
        if (das->mode == MODE_ENKF || das->mode == MODE_HYBRID)
            das_getmemberfname(das, ot->varnames[1], mem, fname2);
        else if (das->mode == MODE_ENOI)
            das_getbgfname(das, ot->varnames[1], fname2);
        model_readfield(m, fname2, ot->varnames[1], ksurf, bias, das->mode == MODE_ENOI ? 1 : 0);

        for (i = 0; i < nv; ++i)
            src0[i] -= bias[i];
        free(bias);
    }

    evaluate_2d_obs(m, allobs, nobs, obsids, fname, src, dst);

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
    int masklog = das_isstatic(das, mem);
    int ni, nj, nk;
    void* src;
    float* src0;
    char tag_offset[MAXSTRLEN];
    void* offset_data = NULL;

    assert(ot->nvar == 1);      /* should we care? */
    model_getvargridsize(m, mvid, &ni, &nj, &nk);

    if (nj > 0) {
        src = alloc3d(nk, nj, ni, sizeof(float));
        src0 = ((float***) src)[0][0];

        if (nk > 1)
            model_read3dfield(m, fname, ot->varnames[0], src0, masklog);
        else if (nk == 1)
            model_readfield(m, fname, ot->varnames[0], 0, src0, masklog);

        snprintf(tag_offset, MAXSTRLEN, "%s:OFFSET", allobs->obstypes[allobs->data[obsids[0]].type].name);
        offset_data = model_getdata(m, tag_offset);
        if (offset_data != NULL) {
            if (model_getdataalloctype(m, tag_offset) == ALLOCTYPE_3D) {
                float* offset0 = ((float***) offset_data)[0][0];
                size_t ijk;

                for (ijk = 0; ijk < (size_t) ni * nj * nk; ++ijk)
                    src0[ijk] -= offset0[ijk];
            } else if (model_getdataalloctype(m, tag_offset) == ALLOCTYPE_1D) {
                size_t ij, k;

                for (k = 0; k < nk; ++k) {
                    float* srck = ((float***) src)[k][0];
                    float offsetk = ((float*) offset_data)[k];

                    for (ij = 0; ij < ni * nj; ++ij)
                        srck[ij] -= offsetk;
                }
            } else
                enkf_quit("obstype = %s: offset variable must be either 3D or 1D for a 3D observation type on structured grid", ot->name);
        }
    } else {
        src = alloc2d(nk, ni, sizeof(float));
        src0 = ((float**) src)[0];

        if (nk > 1)
            model_read3dfield(m, fname, ot->varnames[0], src0, masklog);
        else if (nk == 1)
            model_readfield(m, fname, ot->varnames[0], 0, src0, masklog);

        snprintf(tag_offset, MAXSTRLEN, "%s:OFFSET", allobs->obstypes[allobs->data[obsids[0]].type].name);
        offset_data = model_getdata(m, tag_offset);
        if (offset_data != NULL) {
            if (model_getdataalloctype(m, tag_offset) == ALLOCTYPE_2D) {
                float* offset0 = ((float**) offset_data)[0];
                size_t ik;

                for (ik = 0; ik < (size_t) ni * nk; ++ik)
                    src0[ik] -= offset0[ik];
            } else if (model_getdataalloctype(m, tag_offset) == ALLOCTYPE_1D) {
                size_t i, k;

                for (k = 0; k < nk; ++k) {
                    float* srck = ((float**) src)[k];
                    float offsetk = ((float*) offset_data)[k];

                    for (i = 0; i < ni; ++i)
                        srck[i] -= offsetk;
                }
            } else
                enkf_quit("obstype = %s: offset variable must be either 2D or 1D for a 3D observation type on unstructured grid", ot->name);
        }
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
    int masklog = das_isstatic(das, mem);
    void* src1 = NULL;
    void* src2 = NULL;
    void* src = NULL;
    int ni, nj, nk;
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
        if (model_getdataalloctype(m, tag_offset) == ALLOCTYPE_3D || model_getdataalloctype(m, tag_offset) == ALLOCTYPE_2D)
            enkf_quit("obstype = %s:  \"lowmem\" subsurface H-function works with 1D offsets only", ot->name);
    }

    if (nj > 0) {
        src = calloc(nk, sizeof(float**));
        src1 = alloc2d(nj, ni, sizeof(float));
        src2 = alloc2d(nj, ni, sizeof(float));
    } else {
        src = calloc(nk, sizeof(float*));
        src1 = calloc(ni, sizeof(float));
        src2 = calloc(ni, sizeof(float));
    }

    /*
     * (we assume that obsids are not used elsewhere and thus can be reordered)
     */
    qsort_r(obsids, nobs, sizeof(int), cmp_obs_byfk, allobs->data);

    for (k1 = 0, i1 = 0, i2 = -1, k1_now = -1, k2_now = -1; k1 < nk && i2 < nobs; ++k1) {
        int k2 = k1 + 1;
        int k1_isnew = 0;

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
            model_readfield(m, fname, ot->varnames[0], k1, (nj > 0) ? ((float**) src1)[0] : src1, masklog);
            k1_now = k1;
            k1_isnew = 1;
        }
        if (k2 < nk) {
            model_readfield(m, fname, ot->varnames[0], k2, (nj > 0) ? ((float**) src2)[0] : src2, masklog);
            k2_now = k2;
        }

        /*
         * (interpolate)
         */
        if (nj > 0) {
            ((float***) src)[k1] = src1;
            if (k2 < nk)
                ((float***) src)[k2] = src2;
        } else {
            ((float**) src)[k1] = src1;
            if (k2 < nk)
                ((float**) src)[k2] = src2;
        }
        if (offset_data != NULL) {
            int k;

            for (k = (k1_isnew == 1) ? k1 : k2; k <= k2 && k < nk; ++k) {
                float* srck = (nj > 0) ? ((float***) src)[k][0] : ((float**) src)[k];
                float offsetk = ((float*) offset_data)[k];
                size_t ij;

                if (nj > 0) {
                    for (ij = 0; ij < (size_t) ni * nj; ++ij)
                        srck[ij] -= offsetk;
                } else {
                    for (ij = 0; ij < (size_t) ni; ++ij)
                        srck[ij] -= offsetk;
                }
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
    grid* g = model_getvargrid(m, mvid);
    int masklog = das_isstatic(das, mem);

    float*** src = NULL;
    int ni, nj, nk;
    float** mld = NULL;
    char tag_offset[MAXSTRLEN];
    void* offset_data = NULL;
    float** bias = NULL;
    char fname2[MAXSTRLEN];
    size_t i;

    if (model_getvarislog(m, mvid))
        enkf_quit("%s: H-function \"wsurfbias\" can not be specified for model variables with APPLYLOG = 1", ot->name);

    if (ot->nvar < 2)
        enkf_quit("%s: second variable has to be defined for the observation type when using H-function \"wsurfbias\"", ot->name);
    {
        int mvid2 = model_getvarid(m, ot->varnames[1], 1);

        if (model_getvargridid(m, mvid) != model_getvargridid(m, mvid2))
            enkf_quit("H_surf_biased(): variables \"%s\" and \"%s\" are defined on different grids", ot->varnames[0], ot->varnames[1]);
    }

    model_getvargridsize(m, mvid, &ni, &nj, &nk);
    if (nj <= 0)
        enkf_quit("H-function \"wsurfbias\" is currently implemented for structured horizontal grids only");
    src = alloc3d(nk, nj, ni, sizeof(float));

    /*
     * this part is similar to H_subsurf_standard()
     */
    if (nk > 1)
        model_read3dfield(m, fname, ot->varnames[0], src[0][0], masklog);
    else if (nk == 1)
        model_readfield(m, fname, ot->varnames[0], 0, src[0][0], masklog);

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
     * now correct for surface bias
     */
    if (das->mode != MODE_HYBRID || mem <= das->nmem_dynamic) {
        bias = alloc2d(nj, ni, sizeof(float));
        if (das->mode == MODE_ENKF || das->mode == MODE_HYBRID)
            das_getmemberfname(das, ot->varnames[1], mem, fname2);
        else if (das->mode == MODE_ENOI)
            das_getbgfname(das, ot->varnames[1], fname2);
        assert(ncu_getnD(fname2, ot->varnames[1]) == 2);
        model_readfield(m, fname2, ot->varnames[1], 0, bias[0], das->mode == MODE_ENOI ? 1 : 0);

        if (isnan(ot->mld_threshold) && ot->mld_varname == NULL)
            enkf_quit("\"MLD_THRESH\" or \"MLD_VARNAME\" must be specified for observation type \"%s\" to use H function \"wsurfbias\"", ot->name);
        if (!isnan(ot->mld_threshold) && ot->mld_varname != NULL)
            enkf_quit("both \"MLD_THRESH\" and \"MLD_VARNAME\" are specified for observation type \"%s\"", ot->name);
        mld = alloc2d(nj, ni, sizeof(float));
        if (ot->mld_varname != NULL) {
            char fname_mld[MAXSTRLEN];

            if (model_getvarid(m, ot->mld_varname, 0) < 0)
                enkf_quit("\"MLD_VARNAME = %s\" for observation type \"%s\" does not exist among model variables", ot->mld_varname, ot->name);
            if (das->mode == MODE_ENKF || das->mode == MODE_HYBRID)
                das_getmemberfname(das, ot->mld_varname, mem, fname_mld);
            else if (das->mode == MODE_ENOI)
                das_getbgfname(das, ot->mld_varname, fname_mld);
            model_readfield(m, fname_mld, ot->mld_varname, 0, mld[0], das->mode == MODE_ENOI ? 1 : 0);
        } else {
            if (das->mode == MODE_ENKF || das->mode == MODE_HYBRID)
                das_calcmld(das, ot, src, mld);
            else if (das->mode == MODE_ENOI) {
                char tag[MAXSTRLEN];

                snprintf(tag, MAXSTRLEN, "%s:MLD", ot->name);
                if (mem <= 0) { /* background */
                    das_calcmld(das, ot, src, mld);
                    model_addorreplacedata(m, tag, mvid, ALLOCTYPE_2D, mld);
                } else          /* members */
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

                if (fabs(fi_prev - o->fij[0]) > EPS_IJ || fabs(fj_prev - o->fij[1]) > EPS_IJ) {
                    vmld = grid_interpolate2d(g, o->fij, mld);
                    vbias = grid_interpolate2d(g, o->fij, bias);
                    fi_prev = o->fij[0];
                    fj_prev = o->fij[1];
                }

                if (!isfinite(vmld))
                    continue;

                assert(isfinite(vbias));
                dst[ii] -= mldtaper(vmld, o->depth) * vbias;
            }
        }

        if (das->mode == MODE_ENKF || das->mode == MODE_HYBRID)
            free(mld);
        free(bias);
    }
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
    void* src = NULL;
    void* srcsum = NULL;
    float* src0 = NULL;
    float* srcsum0 = NULL;
    char tag_offset[MAXSTRLEN];
    void* offset_data = NULL;
    size_t k, i;

    if (model_getvarislog(m, mvid))
        enkf_quit("%s: H-function \"vertsum\" can not be specified for model variables with APPLYLOG = 1", ot->name);

    {
        int vtype = grid_getvtype(model_getvargrid(m, mvid));

        if (vtype != GRIDVTYPE_SIGMA && vtype != GRIDVTYPE_HYBRID)
            enkf_quit("obstype = %s: H-function \"vertsum\" can only be used for variables on either sigma or hybrid vertical grids", ot->name);
    }

    assert(ot->nvar == 1);      /* should we care? */
    model_getvargridsize(m, mvid, &ni, &nj, &nk);
    if (nj > 0) {
        src = alloc2d(nj, ni, sizeof(float));
        srcsum = alloc2d(nj, ni, sizeof(float));
        src0 = ((float**) src)[0];
        srcsum0 = ((float**) srcsum)[0];
    } else {
        src = calloc(ni, sizeof(float));
        srcsum = calloc(ni, sizeof(float));
        src0 = src;
        srcsum0 = srcsum;
    }

    model_readfield(m, fname, ot->varnames[0], 0, srcsum0, 0);
    for (k = 1; k < nk; ++k) {
        size_t nij = (nj > 0) ? ni * nj : ni;

        model_readfield(m, fname, ot->varnames[0], k, src0, 0);
        for (i = 0; i < nij; i++)
            srcsum0[i] += src0[i];
    }
    free(src);
    src0 = NULL;

    snprintf(tag_offset, MAXSTRLEN, "%s:OFFSET", allobs->obstypes[allobs->data[obsids[0]].type].name);
    offset_data = model_getdata(m, tag_offset);
    if (offset_data != NULL) {
        float* offset0 = (nj > 0) ? ((float**) offset_data)[0] : offset_data;
        size_t nij = (nj > 0) ? ni * nj : ni;

        if (nj > 0) {
            if (model_getdataalloctype(m, tag_offset) != ALLOCTYPE_2D)
                enkf_quit("obstype = %s: offset variable must be 2D to be used in H-function \"vertsum\" (for structured grids)", ot->name);
        } else {
            if (model_getdataalloctype(m, tag_offset) != ALLOCTYPE_1D)
                enkf_quit("obstype = %s: offset variable must be 1D to be used in H-function \"vertsum\" (for unstructured grids)", ot->name);
        }

        for (i = 0; i < nij; ++i)
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
    int ni, nj, nk, nij;
    float* src = NULL;
    void* sum = NULL;
    float* w = NULL;
    float* sumw = NULL;
    char tag_offset[MAXSTRLEN];
    void* offset_data = NULL;
    size_t k, i;

    if (model_getvarislog(m, mvid))
        enkf_quit("%s: H-function \"vertwavg\" can not be specified for model variables with APPLYLOG = 1", ot->name);

    if (ot->nvar < 2)
        enkf_quit("%s: second variable has to be defined for the observation type when using H-function \"vertavg\"", ot->name);
    {
        int mvid2 = model_getvarid(m, ot->varnames[1], 1);

        if (model_getvargridid(m, mvid) != model_getvargridid(m, mvid2))
            enkf_quit("H_vertwavg(): variables \"%s\" and \"%s\" are defined on different grids", ot->varnames[0], ot->varnames[1]);
    }

    model_getvargridsize(m, mvid, &ni, &nj, &nk);
    if (nj > 0) {
        nij = ni * nj;
        sum = alloc2d(nj, ni, sizeof(float));
    } else {
        nij = ni;
        sum = calloc(ni, sizeof(float));
    }
    src = calloc(nij, sizeof(float));
    w = calloc(nij, sizeof(float));
    sumw = calloc(nij, sizeof(float));

    /*
     * get sum of weights 
     */
    model_readfield(m, fname, ot->varnames[1], 0, sumw, 0);
    for (k = 1; k < nk; ++k) {
        model_readfield(m, fname, ot->varnames[1], k, w, 0);
        for (i = 0; i < (size_t) nij; i++)
            sumw[i] += w[i];
    }

    /*
     * calculate weighted average 
     */
    for (k = 0; k < nk; ++k) {
        float* sum0 = (nj > 0) ? sum : ((float**) sum)[0];

        model_readfield(m, fname, ot->varnames[0], k, src, 0);
        model_readfield(m, fname, ot->varnames[1], k, w, 0);
        for (i = 0; i < (size_t) nij; i++) {
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
        float* sum0 = (nj > 0) ? sum : ((float**) sum)[0];

        if (model_getdataalloctype(m, tag_offset) != ALLOCTYPE_2D)
            enkf_quit("obstype = %s: offset variable must be 2D to be used in H-function \"vertavg\"", ot->name);

        for (i = 0; i < (size_t) nij; ++i)
            sum0[i] -= offset0[i];
    }

    evaluate_2d_obs(m, allobs, nobs, obsids, fname, sum, dst);
    free(sum);
}
