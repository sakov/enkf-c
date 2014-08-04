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
#include "ncw.h"
#include "definitions.h"
#include "utils.h"
#include "enkfprm.h"
#include "observations.h"
#include "model.h"
#include "model2obs.h"

/**
 */
static void interpolate_2d_obs(model* m, observations* allobs, int nobs, int obsids[], char fname[], float** v, ENSOBSTYPE out[])
{
    int** mask = model_getnumlevels(m);
    int ni, nj, nk;
    int i;

    model_getdims(m, &ni, &nj, &nk);

    for (i = 0; i < nobs; ++i) {
        int ii = obsids[i];
        observation* o = &allobs->data[ii];

        assert(out[ii] == 0.0);
        out[ii] = interpolate2d(o->fi, o->fj, ni, nj, v, mask);
        if (!isfinite(out[ii])) {
            /*
             * the location is on land due to the round-up error after writing
             * and reading from observations.nc
             */
            o->status = STATUS_ROUNDUP;
            o->value = 0.0;
            o->std = STD_BIG;
            continue;
        }
        if (fabs(out[ii]) > STATE_BIGNUM)
            enkf_quit("obs # %d: forecast = %.3g for \"%s\"; no point to continue", ii, out[ii], fname);
    }
}

/**
 */
static void interpolate_3d_obs(model* m, observations* allobs, int nobs, int obsids[], char fname[], float*** v, ENSOBSTYPE out[])
{
    int** nlevels = model_getnumlevels(m);
    int ni, nj, nk;
    int i;

    model_getdims(m, &ni, &nj, &nk);

    for (i = 0; i < nobs; ++i) {
        int ii = obsids[i];
        observation* o = &allobs->data[ii];

        assert(out[ii] == 0.0);
        out[ii] = interpolate3d(o->fi, o->fj, o->fk, ni, nj, nk, v, nlevels);
        if (!isfinite(out[ii])) {
            /*
             * the location is on land due to the round-up error after writing
             * and reading from observations.nc
             */
            o->status = STATUS_ROUNDUP;
            o->value = 0.0;
            o->std = STD_BIG;
            continue;
        }
        if (fabs(out[ii]) > STATE_BIGNUM)
            enkf_quit("obs # %d: forecast = %.3g for \"%s\"; no point to continue", ii, out[ii], fname);
    }
}

/**
 */
/**
 */
void H_surf_standard(dasystem* das, int nobs, int obsids[], char fname[], int mem, int t, char varname[], void* psrc, ENSOBSTYPE dst[])
{
    model* m = das->m;
    observations* allobs = das->obs;
    float** src = (float**) psrc;

    model_readfield(m, fname, mem, t, varname, 0, src[0]);
    interpolate_2d_obs(m, allobs, nobs, obsids, fname, src, dst);
}

/**
 */
void H_sla_standard(dasystem* das, int nobs, int obsids[], char fname[], int mem, int t, char varname[], void* psrc, ENSOBSTYPE dst[])
{
    model* m = das->m;
    observations* allobs = das->obs;
    float** msl = (float**) model_getmodeldata(m, "MSL");
    float** src = (float**) psrc;
    int ni, nj, nk;
    int i, j;

    model_getdims(m, &ni, &nj, &nk);

    model_readfield(m, fname, mem, t, varname, 0, src[0]);
    for (j = 0; j < nj; ++j)
        for (i = 0; i < ni; ++i)
            src[j][i] -= msl[j][i];
    interpolate_2d_obs(m, allobs, nobs, obsids, fname, src, dst);
}

/**
 */
void H_sla_bran(dasystem* das, int nobs, int obsids[], char fname[], int mem, int t, char varname[], void* psrc, ENSOBSTYPE dst[])
{
    int o;
    ENSOBSTYPE mean;

    H_sla_standard(das, nobs, obsids, fname, mem, t, varname, psrc, dst);
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
void H_subsurf_standard(dasystem* das, int nobs, int obsids[], char fname[], int mem, int t, char varname[], void* psrc, ENSOBSTYPE dst[])
{
    model* m = das->m;
    observations* allobs = das->obs;
    float*** src = (float***) psrc;

    model_read3dfield(m, fname, mem, t, varname, src[0][0]);
    interpolate_3d_obs(m, allobs, nobs, obsids, fname, src, dst);
}
