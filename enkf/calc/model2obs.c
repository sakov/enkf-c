/******************************************************************************
 *
 * File:        model2obs.c        
 *
 * Created:     12/2012
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description:
 *
 * Revisions:
 *
 *****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <values.h>
#include "ncw.h"
#include "nan.h"
#include "definitions.h"
#include "utils.h"
#include "enkfprm.h"
#include "observations.h"
#include "model.h"
#include "model2obs.h"

#define EPSMULT (float) 1.000001

/**
 */
static void interpolate_2d_obs(model* m, observations* allobs, int nobs, int obsids[], char fname[], float** v, ENSOBSTYPE out[])
{
    int type_prev = -1;
    int** mask = NULL;
    int ni, nj;
    int periodic_x = -1, periodic_y = -1;
    int i;

    for (i = 0; i < nobs; ++i) {
        int ii = obsids[i];
        observation* o = &allobs->data[ii];

        if (o->type != type_prev) {
            int mvid = model_getvarid(m, allobs->obstypes[o->type].varname);

            assert(mvid >= 0);
            model_getvardims(m, mvid, &ni, &nj, NULL);
            mask = model_getnumlevels(m, mvid);
            periodic_x = grid_isperiodic_x(model_getvargrid(m, mvid));
            periodic_y = grid_isperiodic_y(model_getvargrid(m, mvid));
        }

        assert(out[ii] == 0.0);
        out[ii] = interpolate2d(o->fi, o->fj, ni, nj, v, mask, periodic_x, periodic_y);
        if (!isfinite(out[ii])) {
            /*
             * the location is on land due to the round-up error after writing
             * and reading from observations.nc
             */
            if (floor(o->fi) != floor(o->fi * EPSMULT))
                o->fi *= EPSMULT;
            else if (floor(o->fi) != floor(o->fi / EPSMULT))
                o->fi /= EPSMULT;
            else if (floor(o->fj) != floor(o->fj * EPSMULT))
                o->fj *= EPSMULT;
            else if (floor(o->fj) != floor(o->fj / EPSMULT))
                o->fj /= EPSMULT;
            else {
                o->status = STATUS_ROUNDUP;
                o->value = 0.0;
                o->std = STD_BIG;
                continue;
            }
            out[ii] = interpolate2d(o->fi, o->fj, ni, nj, v, mask, periodic_x, periodic_y);
            if (!isfinite(out[ii])) {
                o->status = STATUS_ROUNDUP;
                o->value = 0.0;
                o->std = STD_BIG;
                continue;
            }
        }
        if (fabs(out[ii]) > STATE_BIGNUM)
            enkf_quit("obs # %d: forecast = %.3g for \"%s\"; no point to continue", ii, out[ii], fname);
    }
}

/**
 */
static void interpolate_3d_obs(model* m, observations* allobs, int nobs, int obsids[], char fname[], float*** v, ENSOBSTYPE out[])
{
    int** nlevels = NULL;
    int type_prev = -1;
    int ni, nj, nk;
    int periodic_x = -1, periodic_y = -1;
    int i;

    for (i = 0; i < nobs; ++i) {
        int ii = obsids[i];
        observation* o = &allobs->data[ii];

        if (o->type != type_prev) {
            int mvid = model_getvarid(m, allobs->obstypes[o->type].varname);

            assert(mvid >= 0);
            model_getvardims(m, mvid, &ni, &nj, &nk);
            nlevels = model_getnumlevels(m, mvid);
            periodic_x = grid_isperiodic_x(model_getvargrid(m, mvid));
            periodic_y = grid_isperiodic_y(model_getvargrid(m, mvid));
        }

        assert(out[ii] == 0.0);
        out[ii] = interpolate3d(o->fi, o->fj, o->fk, ni, nj, nk, v, nlevels, periodic_x, periodic_y);
        if (!isfinite(out[ii])) {
            /*
             * the location is on land due to the round-up error after writing
             * and reading from observations.nc
             */
            if (floor(o->fi) != floor(o->fi * EPSMULT))
                o->fi *= EPSMULT;
            else if (floor(o->fi) != floor(o->fi / EPSMULT))
                o->fi /= EPSMULT;
            else if (floor(o->fj) != floor(o->fj * EPSMULT))
                o->fj *= EPSMULT;
            else if (floor(o->fj) != floor(o->fj / EPSMULT))
                o->fj /= EPSMULT;
            else {
                o->status = STATUS_ROUNDUP;
                o->value = 0.0;
                o->std = STD_BIG;
                continue;
            }
            out[ii] = interpolate3d(o->fi, o->fj, o->fk, ni, nj, nk, v, nlevels, periodic_x, periodic_y);
            if (!isfinite(out[ii])) {
                o->status = STATUS_ROUNDUP;
                o->value = 0.0;
                o->std = STD_BIG;
                continue;
            }
        }
        if (fabs(out[ii]) > STATE_BIGNUM)
            enkf_quit("obs # %d: forecast = %.3g for \"%s\"; no point to continue", ii, out[ii], fname);
    }
}

/**
 */
void H_surf_standard(dasystem* das, int nobs, int obsids[], char fname[], int mem, int t, char varname[], char varname2[], void* psrc, ENSOBSTYPE dst[])
{
    model* m = das->m;
    observations* allobs = das->obs;
    float** src = (float**) psrc;
    float** offset = NULL;
    int mvid = model_getvarid(m, varname);
    int k = grid_gettoplayerid(model_getvargrid(m, mvid));

    assert(varname2 == NULL);

    model_readfield(m, fname, mem, t, varname, k, src[0]);

    offset = model_getmodeldata(m, allobs->obstypes[allobs->data[obsids[0]].type].name);
    if (offset != NULL) {
        int ni, nj;
        float* src0 = src[0];
        float* offset0 = offset[0];
        int i;

        assert(mvid >= 0);
        model_getvardims(m, mvid, &ni, &nj, NULL);
        for (i = 0; i < ni * nj; ++i)
            src0[i] -= offset0[i];
    }

    interpolate_2d_obs(m, allobs, nobs, obsids, fname, src, dst);
}

/**
 */
void H_sla_bran(dasystem* das, int nobs, int obsids[], char fname[], int mem, int t, char varname[], char varname2[], void* psrc, ENSOBSTYPE dst[])
{
    int o;
    ENSOBSTYPE mean;

    assert(varname2 == NULL);
    H_surf_standard(das, nobs, obsids, fname, mem, t, varname, NULL, psrc, dst);
    if (mem <= 0) {             /* only for background */
        mean = 0.0;
        for (o = 0; o < nobs; ++o)
            mean += dst[o];
        mean /= (float) nobs;
        for (o = 0; o < nobs; ++o)
            dst[o] -= mean;
    }
}

/**
 */
void H_surf_biased(dasystem* das, int nobs, int obsids[], char fname[], int mem, int t, char varname[], char varname2[], void* psrc, ENSOBSTYPE dst[])
{
    model* m = das->m;
    observations* allobs = das->obs;
    float** src = (float**) psrc;
    float** offset = NULL;
    float** bias = NULL;
    int mvid = model_getvarid(m, varname);
    int k;
    char fname2[MAXSTRLEN];
    int ni, nj;
    int i, j;

    assert(mvid >= 0);
    k = grid_gettoplayerid(model_getvargrid(m, mvid));

    model_getvardims(m, mvid, &ni, &nj, NULL);
    bias = alloc2d(nj, ni, sizeof(float));
    if (das->mode == MODE_ENKF)
        model_getmemberfname(m, das->ensdir, varname2, mem, fname2);
    else if (das->mode == MODE_ENOI)
        model_getbgfname(m, das->bgdir, varname2, fname2);
    model_readfield(m, fname2, mem, MAXINT, varname2, k, bias[0]);

    model_readfield(m, fname, mem, t, varname, k, src[0]);

    offset = model_getmodeldata(m, allobs->obstypes[allobs->data[obsids[0]].type].name);
    if (offset != NULL) {
        for (j = 0; j < nj; ++j)
            for (i = 0; i < ni; ++i)
                src[j][i] -= offset[j][i] + bias[j][i];
    } else {
        for (j = 0; j < nj; ++j)
            for (i = 0; i < ni; ++i)
                src[j][i] -= bias[j][i];
    }

    interpolate_2d_obs(m, allobs, nobs, obsids, fname, src, dst);

    free2d(bias);
}

/**
 */
void H_subsurf_standard(dasystem* das, int nobs, int obsids[], char fname[], int mem, int t, char varname[], char varname2[], void* psrc, ENSOBSTYPE dst[])
{
    model* m = das->m;
    observations* allobs = das->obs;
    float*** src = (float***) psrc;
    float*** offset = NULL;

    assert(varname2 == NULL);
    model_read3dfield(m, fname, mem, t, varname, src[0][0]);

    offset = model_getmodeldata(m, allobs->obstypes[allobs->data[obsids[0]].type].name);
    if (offset != NULL) {
        int mvid = model_getvarid(m, varname);
        int ni, nj, nk;
        float* src0 = src[0][0];
        float* offset0 = offset[0][0];
        int i;

        assert(mvid >= 0);
        model_getvardims(m, mvid, &ni, &nj, &nk);
        for (i = 0; i < ni * nj * nk; ++i)
            src0[i] -= offset0[i];
    }

    interpolate_3d_obs(m, allobs, nobs, obsids, fname, src, dst);
}
