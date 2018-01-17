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
#include "enkfprm.h"
#include "observations.h"
#include "model.h"
#include "model2obs.h"

#define EPSMULT (float) 1.000001
#define EPS_IJ 0.01

/**
 */
static void interpolate_2d_obs(model* m, observations* allobs, int nobs, int obsids[], char fname[], float** v, ENSOBSTYPE out[])
{
    int otid = allobs->data[obsids[0]].type;
    int mvid = model_getvarid(m, allobs->obstypes[otid].varnames[0], 1);
    int periodic_x = grid_isperiodic_x(model_getvargrid(m, mvid));
    int** mask = model_getnumlevels(m, mvid);
    int ni, nj;
    int i;

    model_getvardims(m, mvid, &ni, &nj, NULL);
    for (i = 0; i < nobs; ++i) {
        int ii = obsids[i];
        observation* o = &allobs->data[ii];

        assert(o->type == otid);
        assert(out[ii] == 0.0);
        out[ii] = interpolate2d(o->fi, o->fj, ni, nj, v, mask, periodic_x);
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
static void interpolate_3d_obs(model* m, observations* allobs, int nobs, int obsids[], char fname[], float*** v, ENSOBSTYPE out[])
{
    int otid = allobs->data[obsids[0]].type;
    int mvid = model_getvarid(m, allobs->obstypes[otid].varnames[0], 1);
    void* grid = model_getvargrid(m, mvid);
    int ksurf = grid_getsurflayerid(grid);
    int periodic_x = grid_isperiodic_x(grid);
    int** nlevels = model_getnumlevels(m, mvid);
    int ni, nj, nk;
    int i;

    model_getvardims(m, mvid, &ni, &nj, &nk);
    for (i = 0; i < nobs; ++i) {
        int ii = obsids[i];
        observation* o = &allobs->data[ii];

        assert(o->type == otid);
        assert(out[ii] == 0.0);
        out[ii] = interpolate3d(o->fi, o->fj, o->fk, ni, nj, nk, ksurf, v, nlevels, periodic_x);
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
void H_surf_standard(dasystem* das, int nobs, int obsids[], char fname[], int mem, int t, void* psrc, ENSOBSTYPE dst[])
{
    model* m = das->m;
    observations* allobs = das->obs;
    int otid = allobs->data[obsids[0]].type;
    obstype* ot = &allobs->obstypes[otid];
    float** src = (float**) psrc;
    char tag_offset[MAXSTRLEN];
    float** offset = NULL;
    int mvid = model_getvarid(m, ot->varnames[0], 1);
    int k = grid_getsurflayerid(model_getvargrid(m, mvid));

    assert(ot->nvar == 1);      /* should we care? */

    if (getnumlevels(fname, ot->varnames[0]) == 1)
        k = 0;
    model_readfield(m, fname, t, ot->varnames[0], k, src[0]);

    snprintf(tag_offset, MAXSTRLEN, "%s:OFFSET", ot->name);
    offset = model_getdata(m, tag_offset);
    if (offset != NULL) {
        int ni, nj;
        float* src0 = src[0];
        float* offset0 = offset[0];
        int i;

        assert(model_getdataalloctype(m, tag_offset) == ALLOCTYPE_2D);
        assert(mvid >= 0);
        model_getvardims(m, mvid, &ni, &nj, NULL);
        for (i = 0; i < ni * nj; ++i)
            src0[i] -= offset0[i];
    }

    interpolate_2d_obs(m, allobs, nobs, obsids, fname, src, dst);
}

/**
 */
void H_surf_biased(dasystem* das, int nobs, int obsids[], char fname[], int mem, int t, void* psrc, ENSOBSTYPE dst[])
{
    model* m = das->m;
    observations* allobs = das->obs;
    int otid = allobs->data[obsids[0]].type;
    obstype* ot = &allobs->obstypes[otid];
    float** src = (float**) psrc;
    float* src0 = src[0];
    char tag_offset[MAXSTRLEN];
    float** offset = NULL;
    float* bias = NULL;
    int mvid = model_getvarid(m, ot->varnames[0], 1);
    int mvid2;
    char fname2[MAXSTRLEN];
    int ni, nj, ksurf;
    int i, nv;

    if (ot->nvar < 2)
        enkf_quit("%s: second variable has to be defined for the observation type when using H-function \"biased\"", ot->name);
    mvid2 = model_getvarid(m, ot->varnames[1], 1);
    if (model_getvargridid(m, mvid) != model_getvargridid(m, mvid2))
        enkf_quit("H_surf_biased(): variables \"%s\" and \"%s\" are defined on different grids", ot->varnames[0], ot->varnames[1]);

    ksurf = grid_getsurflayerid(model_getvargrid(m, mvid));
    model_getvardims(m, mvid, &ni, &nj, NULL);
    nv = ni * nj;

    bias = malloc(nv * sizeof(float));
    if (das->mode == MODE_ENKF)
        model_getmemberfname(m, das->ensdir, ot->varnames[1], mem, fname2);
    else if (das->mode == MODE_ENOI)
        model_getbgfname(m, das->bgdir, ot->varnames[1], fname2);
    model_readfield(m, fname2, INT_MAX, ot->varnames[1], ksurf, bias);

    model_readfield(m, fname, t, ot->varnames[0], ksurf, src0);

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

    interpolate_2d_obs(m, allobs, nobs, obsids, fname, src, dst);

    free(bias);
}

/**
 */
void H_subsurf_standard(dasystem* das, int nobs, int obsids[], char fname[], int mem, int t, void* psrc, ENSOBSTYPE dst[])
{
    model* m = das->m;
    observations* allobs = das->obs;
    int otid = allobs->data[obsids[0]].type;
    obstype* ot = &allobs->obstypes[otid];
    float*** src = (float***) psrc;
    char tag_offset[MAXSTRLEN];
    float*** offset = NULL;

    assert(ot->nvar == 1);      /* should we care? */
    model_read3dfield(m, fname, t, ot->varnames[0], src[0][0]);

    snprintf(tag_offset, MAXSTRLEN, "%s:OFFSET", allobs->obstypes[allobs->data[obsids[0]].type].name);
    offset = model_getdata(m, tag_offset);
    if (offset != NULL) {
        int mvid = model_getvarid(m, ot->varnames[0], 1);
        int ni, nj, nk;
        float* src0 = src[0][0];
        float* offset0 = offset[0][0];
        int i;

        assert(model_getdataalloctype(m, tag_offset) == ALLOCTYPE_3D);
        model_getvardims(m, mvid, &ni, &nj, &nk);
        for (i = 0; i < ni * nj * nk; ++i)
            src0[i] -= offset0[i];
    }

    interpolate_3d_obs(m, allobs, nobs, obsids, fname, src, dst);
}

#define MLD_TRANSITION 0.1

static double mldtaper(double mld, double z)
{
    double v;

    assert(z >= 0.0);
    v = (z / mld - 1) / MLD_TRANSITION;
    if (v < -1.0)
        return 1.0;
    else if (v > 1.0)
        return 0.0;

    return cos((v + 1) * M_PI / 2.0) * 0.5 + 0.5;
}

/** Projects surface bias into subsurface based on the mixed layer depth.
 */
void H_subsurf_wsurfbias(dasystem* das, int nobs, int obsids[], char fname[], int mem, int t, void* psrc, ENSOBSTYPE dst[])
{
    model* m = das->m;
    observations* allobs = das->obs;
    int otid = allobs->data[obsids[0]].type;
    obstype* ot = &allobs->obstypes[otid];
    int mvid = model_getvarid(m, ot->varnames[0], 1);
    int mvid2;
    int periodic_x = grid_isperiodic_x(model_getvargrid(m, mvid));
    int** mask = model_getnumlevels(m, mvid);

    float*** src = (float***) psrc;
    float** mld = NULL;
    char tag_offset[MAXSTRLEN];
    float*** offset = NULL;
    float** bias = NULL;
    char fname2[MAXSTRLEN];
    int ni, nj, nk;
    int i;

    if (ot->nvar < 2)
        enkf_quit("%s: second variable has to be defined for the observation type when using H-function \"wsurfbias\"", ot->name);
    mvid2 = model_getvarid(m, ot->varnames[1], 1);
    if (model_getvargridid(m, mvid) != model_getvargridid(m, mvid2))
        enkf_quit("H_surf_biased(): variables \"%s\" and \"%s\" are defined on different grids", ot->varnames[0], ot->varnames[1]);
    model_getvardims(m, mvid, &ni, &nj, &nk);

    /*
     * this part is similar to H_subsurf_standard()
     */
    model_read3dfield(m, fname, t, ot->varnames[0], src[0][0]);

    snprintf(tag_offset, MAXSTRLEN, "%s:OFFSET", ot->name);
    offset = model_getdata(m, tag_offset);
    if (offset != NULL) {
        float* src0 = src[0][0];
        float* offset0 = offset[0][0];

        assert(model_getdataalloctype(m, tag_offset) == ALLOCTYPE_3D);
        for (i = 0; i < ni * nj * nk; ++i)
            src0[i] -= offset0[i];
    }

    interpolate_3d_obs(m, allobs, nobs, obsids, fname, src, dst);

    /*
     * now correct for the surface bias
     */
    bias = alloc2d(nj, ni, sizeof(float));
    if (das->mode == MODE_ENKF)
        model_getmemberfname(m, das->ensdir, ot->varnames[1], mem, fname2);
    else if (das->mode == MODE_ENOI)
        model_getbgfname(m, das->bgdir, ot->varnames[1], fname2);
    assert(!is3d(fname2, ot->varnames[1]));
    model_readfield(m, fname2, INT_MAX, ot->varnames[1], 0, bias[0]);

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
            model_getmemberfname(m, das->ensdir, ot->mld_varname, mem, fname_mld);
            model_readfield(m, fname_mld, 0, ot->mld_varname, 0, mld[0]);
        } else if (das->mode == MODE_ENOI) {
            model_getbgfname(m, das->bgdir, ot->mld_varname, fname_mld);
            model_readfield(m, fname_mld, 0, ot->mld_varname, 0, mld[0]);
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
            int ii = obsids[i];
            observation* o = &allobs->data[ii];

            if (fabs(fi_prev - o->fi) > EPS_IJ || fabs(fj_prev - o->fj) > EPS_IJ) {
                vmld = interpolate2d(o->fi, o->fj, ni, nj, mld, mask, periodic_x);
                vbias = interpolate2d(o->fi, o->fj, ni, nj, bias, mask, periodic_x);
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
}
