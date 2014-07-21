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
dasystem* das_create(enkfprm* prm)
{
    dasystem* das = malloc(sizeof(dasystem));
    int i;

    das->prmfname = strdup(prm->fname);
    das->mode = prm->mode;
    das->scheme = prm->scheme;
    das->target = prm->target;
    das->ensdir = strdup(prm->ensdir);
    if (prm->bgdir != NULL)
        das->bgdir = strdup(prm->bgdir);
    else
        das->bgdir = NULL;
    das->nmem = prm->enssize;
    das->obs = obs_create_fromprm(prm);

    das->m = model_create(prm);

    das->S = NULL;
    das->s_mode = S_MODE_NONE;
    das->Hx = NULL;
    das->s_f = NULL;
    das->std_f = NULL;
    das->s_a = NULL;
    das->std_a = NULL;
    das->kfactor = prm->kfactor;
    das->locrad = prm->locrad;
    das->stride = prm->stride;
    das->nfields = 0;
    das->fields = NULL;
    das->fieldbufsize = prm->fieldbufsize;

    das->nregions = prm->nregions;
    if (das->nregions > 0)
        das->regions = malloc(das->nregions * sizeof(region));
    else
        das->regions = NULL;
    for (i = 0; i < das->nregions; ++i) {
        region* rin = &das->regions[i];
        region* rout = &prm->regions[i];

        rin->name = strdup(rout->name);
        rin->lon1 = rout->lon1;
        rin->lon2 = rout->lon2;
        rin->lat1 = rout->lat1;
        rin->lat2 = rout->lat2;
    }

    das->nplogs = 0;
    das->plogs = NULL;
    for (i = 0; i < prm->nplogs; ++i) {
        double lon, lat;

        model_fij2ll(das->m, (double) prm->plogs[i].i, (double) prm->plogs[i].j, &lon, &lat);
        if (isnan(lon + lat)) {
            enkf_printf("  WARNING: %s: POINTLOG %d %d: point outside the grid\n", das->prmfname, prm->plogs[i].i, prm->plogs[i].j);
            continue;
        }

        if (das->nplogs % NPLOGS_INC == 0)
            das->plogs = realloc(das->plogs, (das->nplogs + NPLOGS_INC) * sizeof(pointlog));

        das->plogs[das->nplogs].i = prm->plogs[i].i;
        das->plogs[das->nplogs].j = prm->plogs[i].j;
        das->nplogs++;
    }

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

    return das;
}
#endif

#if defined(ENKF_POST)
/**
 */
dasystem* das_create(enkfprm* prm)
{
    dasystem* das = malloc(sizeof(dasystem));

    das->prmfname = strdup(prm->fname);
    das->mode = prm->mode;
    das->scheme = prm->scheme;
    das->target = prm->target;
    das->ensdir = strdup(prm->ensdir);
    if (prm->bgdir != NULL)
        das->bgdir = strdup(prm->bgdir);
    else
        das->bgdir = NULL;
    das->nmem = prm->enssize;
    das->obs = NULL;

    das->m = model_create(prm);

    das->S = NULL;
    das->s_mode = S_MODE_NONE;
    das->Hx = NULL;
    das->s_f = NULL;
    das->std_f = NULL;
    das->s_a = NULL;
    das->std_a = NULL;
    das->kfactor = NaN;
    das->locrad = NaN;
    das->stride = prm->stride;
    das->nfields = 0;
    das->fields = NULL;
    das->fieldbufsize = prm->fieldbufsize;

    das->nregions = 0;
    das->regions = NULL;

    das->nbadbatchspecs = 0;
    das->badbatchspecs = NULL;

    das->nplogs = 0;
    das->plogs = NULL;

    return das;
}
#endif

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
    if (das->Hx != NULL)
        free(das->Hx);
    if (das->s_f != NULL) {
        free(das->s_f);
        free(das->std_f);
    }
    if (das->s_a != NULL) {
        free(das->s_a);
        free(das->std_a);
    }
    if (das->fields != NULL)
        free(das->fields);
    if (das->nregions > 0) {
        for (i = 0; i < das->nregions; ++i)
            free(das->regions[i].name);
        free(das->regions);
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
void das_getfields(dasystem* das)
{
    model* m = das->m;
    int nvar = model_getnvar(m);
    int vid;

    assert(das->nfields == 0);
    assert(das->fields == NULL);

    for (vid = 0; vid < nvar; ++vid) {
        char fname[MAXSTRLEN];
        char* varname = model_getvarname(m, vid);
        int nk, k;

        model_getmemberfname(m, das->ensdir, varname, 1, fname);
        nk = getnlevels(fname, varname);
        for (k = 0; k < nk; ++k) {
            field* f;

            if (das->nfields % NFIELDS_INC == 0)
                das->fields = realloc(das->fields, (das->nfields + NFIELDS_INC) * sizeof(field));
            f = &das->fields[das->nfields];
            f->id = das->nfields;
            f->varid = vid;
            strcpy(f->varname, varname);
            f->level = k;
            das->nfields++;
        }
    }
}
