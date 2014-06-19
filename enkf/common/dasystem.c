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

    das->nvar = prm->nvar;
    if (das->nvar > 0) {
        das->varnames = malloc(das->nvar * sizeof(char*));
        das->inflations = malloc(das->nvar * sizeof(float));
        for (i = 0; i < das->nvar; ++i) {
            das->varnames[i] = strdup(prm->varnames[i]);
            das->inflations[i] = prm->inflations[i];
        }
    } else
        das->varnames = NULL;
    das->m = model_create(prm);
    if (prm->msl_fname != NULL) {
        int ni, nj, nk;

        model_getdims(das->m, &ni, &nj, &nk);
        das->msl = alloc2d(nj, ni, sizeof(float));
        readfield(prm->msl_fname, 0, prm->msl_varname, das->msl[0]);
    } else
        das->msl = NULL;
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
    das->regions = malloc(das->nregions * sizeof(region));
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
    obs_destroy(das->obs);
    model_destroy(das->m);
    if (das->nvar > 0) {
        for (i = 0; i < das->nvar; ++i)
            free(das->varnames[i]);
        free(das->varnames);
        free(das->inflations);
    }
    if (das->msl != NULL)
        free2d(das->msl);
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
    free(das);
    distribute_free();
}

/**
 */
void das_getnmem(dasystem* das)
{
    das->nmem = 0;
    while (1) {
        char fname[MAXSTRLEN];
        int i;

        for (i = 0; i < das->nvar; ++i) {
            model_getmemberfname(das->m, das->ensdir, das->varnames[i], das->nmem + 1, fname);
            if (!file_exists(fname))
                break;
        }
        if (i == das->nvar)
            das->nmem++;
        else
            break;
    }
}

/** Looks for all horizontal fields of the model to be updated.
 */
void das_getfields(dasystem* das)
{
    int vid;

    assert(das->nfields == 0);
    assert(das->fields == NULL);

    for (vid = 0; vid < das->nvar; ++vid) {
        char fname[MAXSTRLEN];
        int nk, k;

        model_getmemberfname(das->m, das->ensdir, das->varnames[vid], 1, fname);
        nk = getnlevels(fname, das->varnames[vid]);
        for (k = 0; k < nk; ++k) {
            field* f;

            if (das->nfields % NFIELDS_INC == 0)
                das->fields = realloc(das->fields, (das->nfields + NFIELDS_INC) * sizeof(field));
            f = &das->fields[das->nfields];
            f->id = das->nfields;
            f->varid = vid;
            strcpy(f->varname, das->varnames[vid]);
            f->level = k;
            das->nfields++;
        }
    }
}
