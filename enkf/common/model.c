/******************************************************************************
 *
 * File:        model.c        
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
#include "nan.h"
#include "definitions.h"
#include "utils.h"
#include "grid.h"
#include "enkfprm.h"
#include "model.h"
#include "allmodels.h"

struct model {
    char* tag;

    int ni;
    int nj;
    int nk;

    void* grid;
    float** depth;
    int** numlevels;

    char* ensdir;
    model_getmemberfname_fn getmemberfname;
    model_getmemberfnameasync_fn getmemberfname_async;
    model_getbgfname_fn getbgfname;
    model_getbgfnameasync_fn getbgfname_async;
    model_readfield_fn readfield;
    model_read3dfield_fn read3dfield;
    model_writefield_fn writefield;
};

/**
 */
model* model_create(enkfprm* prm)
{
    model* m = malloc(sizeof(model));
    modelsetup_fn setupfn = get_modelsetupfn(prm->modeltype);

    m->tag = strdup(prm->modeltype);
    m->grid = grid_create();

    m->numlevels = NULL;
    m->depth = NULL;
    m->ensdir = strdup(prm->ensdir);
    m->getmemberfname = NULL;
    m->getmemberfname_async = NULL;
    m->getbgfname = NULL;
    m->getbgfname_async = NULL;
    m->readfield = NULL;
    m->read3dfield = NULL;
    m->writefield = NULL;

    setupfn(m, prm->gridspec);

    grid_getdimensions(model_getgrid(m), &m->ni, &m->nj, &m->nk);

    return m;
}

/**
 */
void model_destroy(model* m)
{
    grid_destroy(m->grid);

    free(m->tag);
    if (m->depth != NULL)
        free2d(m->depth);
    if (m->numlevels != NULL)
        free2d(m->numlevels);
    free(m->ensdir);
    free(m);
}

/**
 */
void model_setdepth(model* m, float** depth)
{
    m->depth = depth;
}

/**
 */
void model_setnumlevels(model* m, int** numlevels)
{
    m->numlevels = numlevels;
}

/**
 */
void model_setgetmemberfname_fn(model* m, model_getmemberfname_fn fn)
{
    m->getmemberfname = fn;
}

/**
 */
void model_setgetmemberfnameasync_fn(model* m, model_getmemberfnameasync_fn fn)
{
    m->getmemberfname_async = fn;
}

/**
 */
void model_setbgfname_fn(model* m, model_getbgfname_fn fn)
{
    m->getbgfname = fn;
}

/**
 */
void model_setbgfnameasync_fn(model* m, model_getbgfnameasync_fn fn)
{
    m->getbgfname_async = fn;
}

/**
 */
void model_setreadfield_fn(model* m, model_readfield_fn fn)
{
    m->readfield = fn;
}

/**
 */
void model_setread3dfield_fn(model* m, model_read3dfield_fn fn)
{
    m->read3dfield = fn;
}

/**
 */
void model_setwritefield_fn(model* m, model_writefield_fn fn)
{
    m->writefield = fn;
}

/**
 */
void model_getdims(model* m, int* ni, int* nj, int* nk)
{
    *ni = m->ni;
    *nj = m->nj;
    *nk = m->nk;
}

/**
 */
void* model_getgrid(model* m)
{
    return m->grid;
}

/**
 */
int model_getlontype(model* m)
{
    return grid_getlontype(m->grid);
}

/**
 */
float** model_getdepth(model* m)
{
    return m->depth;
}

/**
 */
int** model_getnumlevels(model* m)
{
    return m->numlevels;
}

/**
 */
void model_getmemberfname(model* m, char ensdir[], char varname[], int mem, char fname[])
{
    m->getmemberfname(m, ensdir, varname, mem, fname);
}

/**
 */
int model_getmemberfname_async(model* m, char ensdir[], char varname[], char otname[], int mem, int time, char fname[])
{
    return m->getmemberfname_async(m, ensdir, varname, otname, mem, time, fname);
}

/**
 */
void model_getbgfname(model* m, char ensdir[], char varname[], char fname[])
{
    m->getbgfname(m, ensdir, varname, fname);
}

/**
 */
int model_getbgfname_async(model* m, char ensdir[], char varname[], char otname[], int time, char fname[])
{
    return m->getbgfname_async(m, ensdir, varname, otname, time, fname);
}

/**
 */
int model_ll2fij(model* m, double x, double y, double* fi, double* fj)
{
    int i1, i2, j1, j2;

    grid_getll2fijfn(m->grid) (m->grid, x, y, fi, fj);

    if (isnan(*fi + *fj))
        return STATUS_OUTSIDE;

    i1 = floor(*fi);
    i2 = ceil(*fi);
    j1 = floor(*fj);
    j2 = ceil(*fj);
    if (m->numlevels[j1][i1] == 0 && m->numlevels[j1][i2] == 0 && m->numlevels[j2][i1] == 0 && m->numlevels[j2][i2] == 0) {
        *fi = NaN;
        *fj = NaN;
        return STATUS_LAND;
    }
    return STATUS_OK;
}

/**
 */
int model_fij2ll(model* m, double fi, double fj, double* lon, double* lat)
{
    grid_getfij2llfn(m->grid) (m->grid, fi, fj, lon, lat);

    if (isnan(*lon + *lat))
        return STATUS_OUTSIDE;
    return STATUS_OK;
}

/**
 */
int model_z2fk(model* m, double fi, double fj, double z, double* fk)
{
    int i1, i2, j1, j2, k2;

    if (isnan(fi + fj)) {
        *fk = NaN;
        return STATUS_OUTSIDE;
    }

    grid_getz2fkfn(m->grid) (m->grid, fi, fj, z, fk);

    if (isnan(*fk))
        return STATUS_OUTSIDE;

    i1 = floor(fi);
    i2 = ceil(fi);
    j1 = floor(fj);
    j2 = ceil(fj);
    k2 = ceil(*fk);
    if (m->numlevels[j1][i1] <= k2 && m->numlevels[j1][i2] <= k2 && m->numlevels[j2][i1] <= k2 && m->numlevels[j2][i2] <= k2) {
        *fk = NaN;
        return STATUS_LAND;
    } else if (m->numlevels[j1][i1] <= k2 || m->numlevels[j1][i2] <= k2 || m->numlevels[j2][i1] <= k2 || m->numlevels[j2][i2] <= k2) {
        float** depth = model_getdepth(m);
        int** mask = model_getnumlevels(m);
        int ni, nj, nk;
        double v;

        model_getdims(m, &ni, &nj, &nk);
        v = interpolate2d(fi, fj, ni, nj, depth, mask);
        if (z > v)
            return STATUS_LAND;
    }

    return STATUS_OK;
}

/**
 */
void model_readfield(model* m, char fname[], int mem, int time, char varname[], int k, float* v)
{
    m->readfield(m, fname, mem, time, varname, k, v);
}

/**
 */
void model_read3dfield(model* m, char fname[], int mem, int time, char varname[], float* v)
{
    m->read3dfield(m, fname, mem, time, varname, v);
}

/**
 */
void model_writefield(model* m, char fname[], int time, char varname[], int k, float* v)
{
    m->writefield(m, fname, time, varname, k, v);
}
