/******************************************************************************
 *
 * File:        dasystem.c        
 *
 * Created:     12/2012
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description: code for the `das' (Data Assimilation System) structure. Some
 *              relevant chunks of code has been moved to
 *                ensobs.c
 *                obsstats.c
 *                transforms.c
 *                update.c
 *              , just to reduce the size of dasystem.c.
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
#include <unistd.h>
#include "nan.h"
#include "ncw.h"
#include "hash.h"
#include "version.h"
#include "definitions.h"
#include "utils.h"
#include "distribute.h"
#include "enkfprm.h"
#include "observations.h"
#include "dasystem.h"

#define NPLOGS_INC 10
#define NFIELDS_INC 100
#define MPIIDOFFSET 10000

#if defined(ENKF_CALC)
/**
 */
void das_setobstypes(dasystem* das)
{
    int n = das->obs->nobstypes;
    obstype* types = das->obs->obstypes;
    model* m = das->m;
    int i;

    for (i = 0; i < n; ++i) {
        obstype* type = &types[i];
        int vid = model_getvarid(m, types[i].varname, 1);

        type->vid = vid;
        type->gridid = model_getvargridid(m, vid);

        if (type->offset_fname != NULL) {
            if (type->issurface) {
                float** v = NULL;
                int nx, ny;

                model_getvardims(m, vid, &nx, &ny, NULL);
                v = alloc2d(ny, nx, sizeof(float));
                readfield(type->offset_fname, type->offset_varname, 0, v[0]);
                model_addmodeldata(m, type->name, ALLOCTYPE_2D, v);
            } else {
                float*** v = NULL;
                int nx, ny, nz;

                model_getvardims(m, vid, &nx, &ny, &nz);
                v = alloc3d(nz, ny, nx, sizeof(float));
                readfield(type->offset_fname, type->offset_varname, 0, v[0][0]);
                model_addmodeldata(m, type->name, ALLOCTYPE_3D, v);
            }
        }
    }
}
#endif

/**
 */
dasystem* das_create(enkfprm* prm)
{
    dasystem* das = calloc(1, sizeof(dasystem));
    int i;

    das->prmfname = strdup(prm->fname);
    das->mode = prm->mode;
    das->scheme = prm->scheme;
    if (das->mode == MODE_ENKF || !enkf_fstatsonly)
        das->ensdir = strdup(prm->ensdir);
    if (prm->bgdir != NULL)
        das->bgdir = strdup(prm->bgdir);
    das->nmem = prm->enssize;
#if defined(ENKF_CALC)
    das->obs = obs_create_fromprm(prm);
#endif

    das->m = model_create(prm);
#if defined(ENKF_CALC)
    das_setobstypes(das);
#endif

    das->S = NULL;
    das->s_f = NULL;
    das->std_f = NULL;
    das->s_a = NULL;
    das->std_a = NULL;
    das->s_mode = S_MODE_NONE;
#if defined(ENKF_CALC)
    if (!enkf_fstatsonly) {
        das->kfactor = prm->kfactor;
        das->locrad = prm->locrad;
    } else {
        das->kfactor = NaN;
        das->locrad = NaN;
    }
#elif defined(ENKF_UPDATE)
    das->kfactor = NaN;
    das->locrad = NaN;
#endif
    das->stride = prm->stride;
    if (!enkf_fstatsonly)
        das->fieldbufsize = prm->fieldbufsize;

#if defined(ENKF_CALC)
    das->nregions = prm->nregions;
    if (das->nregions > 0)
        das->regions = malloc(das->nregions * sizeof(region));
    else
        das->regions = NULL;
    for (i = 0; i < das->nregions; ++i) {
        region* dst = &das->regions[i];
        region* src = &prm->regions[i];
        int j;

        dst->name = strdup(src->name);
        dst->x1 = src->x1;
        dst->x2 = src->x2;
        dst->y1 = src->y1;
        dst->y2 = src->y2;
        dst->nzints = src->nzints;
        dst->zints = malloc(dst->nzints * sizeof(zint));
        for (j = 0; j < src->nzints; ++j) {
            dst->zints[j].z1 = src->zints[j].z1;
            dst->zints[j].z2 = src->zints[j].z2;
        }
    }
#endif

    if (prm->nplogs > 0)
        das->ht_plogs = ht_create_i2(prm->nplogs);
    for (i = 0; i < prm->nplogs; ++i) {
        /*
         * ( the grid coordinates of pointlogs are defined on grid #0 )
         */
        void* grid = model_getgridbyid(das->m, 0);
        double lon, lat;
        int key[2];

        grid_fij2xy(grid, (double) prm->plogs[i].i, (double) prm->plogs[i].j, &lon, &lat);

        if (isnan(lon + lat)) {
            enkf_printf("  WARNING: %s: POINTLOG %d %d: point outside the grid\n", das->prmfname, prm->plogs[i].i, prm->plogs[i].j);
            continue;
        }

        if (das->nplogs % NPLOGS_INC == 0)
            das->plogs = realloc(das->plogs, (das->nplogs + NPLOGS_INC) * sizeof(pointlog));

        das->plogs[das->nplogs].i = prm->plogs[i].i;
        das->plogs[das->nplogs].j = prm->plogs[i].j;
        key[0] = das->plogs[das->nplogs].i;
        key[1] = das->plogs[das->nplogs].j;
        ht_insert(das->ht_plogs, key, NULL);
        das->nplogs++;
    }

#if defined(ENKF_CALC)
    das->nbadbatchspecs = prm->nbadbatchspecs;
    if (das->nbadbatchspecs > 0) {
        das->badbatchspecs = malloc(das->nbadbatchspecs * sizeof(badbatchspec));
        for (i = 0; i < das->nbadbatchspecs; ++i) {
            badbatchspec* src = &prm->badbatchspecs[i];
            badbatchspec* dst = &das->badbatchspecs[i];

            dst->obstype = strdup(src->obstype);

            dst->maxbias = src->maxbias;
            dst->maxmad = src->maxmad;
            dst->minnobs = src->minnobs;
        }
    }
#endif

#if defined(ENKF_UPDATE)
    das->updatespec = UPDATE_DEFAULT;
#endif

    return das;
}

/**
 */
void das_destroy(dasystem* das)
{
    int i;

    free(das->prmfname);
    free(das->ensdir);
    if (das->bgdir != NULL)
        free(das->bgdir);
#if defined(ENKF_CALC)
    obs_destroy(das->obs);
#endif
    model_destroy(das->m);
    if (das->S != NULL)
        free2d(das->S);
    if (das->s_f != NULL) {
        free(das->s_f);
        free(das->std_f);
    }
    if (das->s_a != NULL) {
        free(das->s_a);
        free(das->std_a);
    }
    if (das->nregions > 0) {
        for (i = 0; i < das->nregions; ++i) {
            free(das->regions[i].name);
            free(das->regions[i].zints);
        }
        free(das->regions);
    }
    if (das->nplogs > 0) {
        ht_destroy(das->ht_plogs);
        free(das->plogs);
    }
    if (das->nbadbatchspecs > 0) {
        for (i = 0; i < das->nbadbatchspecs; ++i)
            free(das->badbatchspecs[i].obstype);

        free(das->badbatchspecs);
    }
    free(das);
    distribute_free();
}

/**
 */
void das_getnmem(dasystem* das)
{
    model* m = das->m;
    int nvar = model_getnvar(m);

    if (das->mode == MODE_NONE)
        enkf_quit("programming error");

    if (das->mode == MODE_ENOI && enkf_fstatsonly) {
        das->nmem = 1;
        return;
    }

    if (das->ensdir == NULL)
        enkf_quit("programming error");

    das->nmem = 0;
    while (1) {
        char fname[MAXSTRLEN];
        int i;

        for (i = 0; i < nvar; ++i) {
            model_getmemberfname(m, das->ensdir, model_getvarname(m, i), das->nmem + 1, fname);
            if (!file_exists(fname))
                break;
        }
        if (i == nvar)
            das->nmem++;
        else
            break;
    }
}

/** Looks for all horizontal fields of the model to be updated.
 */
void das_getfields(dasystem* das, int gridid, int* nfields, field** fields)
{
    model* m = das->m;
    int nvar = model_getnvar(m);
    int vid;

    assert(*nfields == 0);
    assert(*fields == NULL);

    for (vid = 0; vid < nvar; ++vid) {
        char fname[MAXSTRLEN];
        char* varname = model_getvarname(m, vid);
        int nk, k;

        if (gridid >= 0 && model_getvargridid(m, vid) != gridid)
            continue;

        model_getmemberfname(m, das->ensdir, varname, 1, fname);
        nk = getnlevels(fname, varname);
        for (k = 0; k < nk; ++k) {
            field* f;

            if (*nfields % NFIELDS_INC == 0)
                *fields = realloc(*fields, (*nfields + NFIELDS_INC) * sizeof(field));
            f = &(*fields)[*nfields];
            f->id = *nfields;
            f->varid = vid;
            strcpy(f->varname, varname);
            f->level = k;
            (*nfields)++;
        }
    }
}

/**
 */
void das_getfname_X5(dasystem* das, void* grid, char fname[])
{
    if (model_getngrid(das->m) == 1)
        snprintf(fname, MAXSTRLEN, "%s.nc", FNAMEPREFIX_X5);
    else
        snprintf(fname, MAXSTRLEN, "%s-%d.nc", FNAMEPREFIX_X5, grid_getid(grid));
}

/**
 */
void das_getfname_w(dasystem* das, void* grid, char fname[])
{
    if (model_getngrid(das->m) == 1)
        snprintf(fname, MAXSTRLEN, "%s.nc", FNAMEPREFIX_W);
    else
        snprintf(fname, MAXSTRLEN, "%s-%d.nc", FNAMEPREFIX_W, grid_getid(grid));
}

/**
 */
void das_getfname_stats(dasystem* das, void* grid, char fname[])
{
    if (model_getngrid(das->m) == 1)
        snprintf(fname, MAXSTRLEN, "%s.nc", FNAMEPREFIX_STATS);
    else
        snprintf(fname, MAXSTRLEN, "%s-%d.nc", FNAMEPREFIX_STATS, grid_getid(grid));
}
