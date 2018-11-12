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
 *              relevant chunks of code have been moved to
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
#include <math.h>
#include <unistd.h>
#include "ncw.h"
#include "hash.h"
#include "definitions.h"
#include "utils.h"
#include "distribute.h"
#include "enkfprm.h"
#include "observations.h"
#include "dasystem.h"

#define NPLOGS_INC 10
#define NFIELDS_INC 100
#define MPIIDOFFSET 10000

/** Determines ensemble size based on existence of forecast files for each
 ** variable. If das->nmem <= 0 then sets the ensemble size to the maximum
 ** available; otherwhile checks if the declared ensemble size is actually
 ** available.
 */
static void das_setnmem(dasystem* das)
{
    model* m = das->m;
    int nvar = model_getnvar(m);
    char fname[MAXSTRLEN];
    int nmem;

    if (das->mode == MODE_NONE)
        enkf_quit("programming error");

    if (das->mode == MODE_ENOI && enkf_fstatsonly) {
        das->nmem = 1;
        return;
    }

    if (das->ensdir == NULL)
        enkf_quit("programming error");

    nmem = 0;
    while (1) {
        int i;

        for (i = 0; i < nvar; ++i) {
            model_getmemberfname(m, das->ensdir, model_getvarname(m, i), nmem + 1, fname);
            if (!file_exists(fname))
                break;
        }
        if (i == nvar)
            nmem++;
        else
            break;
    }
    if (nmem == 0)
        enkf_quit("das_setnmem(): could not find \"%s\"", fname);
    if (das->nmem > 0) {
        if (nmem < das->nmem)
            enkf_quit("das_setnmem(): could not find \"%s\"", fname);
        else if (nmem > das->nmem)
            enkf_printf("      %d members found on disk; ignoring excess to specified ensemble size\n", nmem);
    }
    das->nmem = nmem;
}

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
    das->alpha = prm->alpha;
    das->nmem = prm->enssize;
#if defined(ENKF_CALC)
    das->obs = obs_create_fromprm(prm);
#endif

    das->m = model_create(prm);
#if defined(ENKF_CALC)
    obstypes_set(das->obs->nobstypes, das->obs->obstypes, das->m);

    das->S = NULL;
    das->s_f = NULL;
    das->std_f = NULL;
    das->s_a = NULL;
    das->std_a = NULL;
    das->s_mode = S_MODE_NONE;
#endif

    enkf_printf("  setting the ensemble size:\n");
    das_setnmem(das);
    enkf_printf("    %d member%s\n", das->nmem, das->nmem == 1 ? "" : "s");

#if defined(HE_VIASHMEM)
    {
        int ierror;
        int* recvcounts = NULL;
        int* displs = NULL;

        ierror = MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &das->sm_comm);
        assert(ierror == MPI_SUCCESS);
        ierror = MPI_Comm_rank(das->sm_comm, &das->sm_comm_rank);
        assert(ierror == MPI_SUCCESS);
        das->sm_comm_ranks = malloc(nprocesses * sizeof(int));
        /*
         * build map of local ranks
         */
        das->sm_comm_ranks[rank] = das->sm_comm_rank;
        recvcounts = malloc(nprocesses * sizeof(int));
        displs = malloc(nprocesses * sizeof(int));
        for (i = 0; i < nprocesses; ++i) {
            recvcounts[i] = 1;
            displs[i] = i;
        }
        ierror = MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, das->sm_comm_ranks, recvcounts, displs, MPI_INT, MPI_COMM_WORLD);
        assert(ierror == MPI_SUCCESS);
        das->sm_comm_win = MPI_WIN_NULL;

        /*
         * Create communicator based on sm_comm_rank
         */
        ierror = MPI_Comm_split(MPI_COMM_WORLD, das->sm_comm_rank, rank, &das->node_comm);
        assert(ierror == MPI_SUCCESS);
        ierror = MPI_Comm_rank(das->node_comm, &das->node_comm_rank);
        assert(ierror == MPI_SUCCESS);
        ierror = MPI_Comm_size(das->node_comm, &das->node_comm_size);
        assert(ierror == MPI_SUCCESS);
        if (das->sm_comm_rank != 0) {
            MPI_Comm_free(&das->node_comm);
            das->node_comm_rank = -1;
        }
        das->node_comm_ranks = malloc(nprocesses * sizeof(int));
        das->node_comm_ranks[rank] = das->node_comm_rank;
        ierror = MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, das->node_comm_ranks, recvcounts, displs, MPI_INT, MPI_COMM_WORLD);
        assert(ierror == MPI_SUCCESS);

        free(recvcounts);
        free(displs);
    }
#endif

#if defined(ENKF_CALC)
    if (!enkf_fstatsonly) {
        das->kfactor = prm->kfactor;
    } else {
        das->kfactor = NAN;
    }
#endif
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

        dst->name = strdup(src->name);
        dst->x1 = src->x1;
        dst->x2 = src->x2;
        dst->y1 = src->y1;
        dst->y2 = src->y2;
    }
#endif

    if (prm->nplogs > 0)
        das->ht_plogs = ht_create_s4(prm->nplogs * 2);
    das->plogs = malloc(sizeof(pointlog) * prm->nplogs);
    das->nplogs = prm->nplogs;
    for (i = 0; i < prm->nplogs; ++i) {
        pointlog* src = &prm->plogs[i];
        pointlog* dst = &das->plogs[i];
        void* grid = NULL;
        unsigned short key[4] = { src->i, src->j, 0, 0 };

        dst->id = src->id;
        dst->i = src->i;
        dst->j = src->j;
        if (src->gridname == NULL) {
            dst->gridid = 0;
            grid = model_getgridbyid(das->m, 0);
            dst->gridname = strdup(grid_getname(grid));
        } else {
            dst->gridname = strdup(src->gridname);
            grid = model_getgridbyname(das->m, src->gridname);
            dst->gridid = grid_getid(grid);
        }
        key[3] = dst->gridid;

        grid_ij2xy(grid, dst->i, dst->j, &dst->lon, &dst->lat);

        if (isnan(dst->lon + dst->lat)) {
            enkf_printf("  WARNING: %s: POINTLOG %d %d: point outside the grid \"%s\"\n", das->prmfname, dst->i, dst->j, dst->gridname);
            continue;
        }

        ht_insert(das->ht_plogs, key, dst);
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

    das->ncformat = prm->ncformat;
    das->nccompression = prm->nccompression;

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
#if defined(ENKF_CALC)
    if (das->S != NULL)
        free(das->S);
    if (das->s_f != NULL) {
        free(das->s_f);
        free(das->std_f);
    }
    if (das->s_a != NULL) {
        free(das->s_a);
        free(das->std_a);
    }
#endif
#if defined (HE_VIASHMEM)
    if (das->sm_comm_win != MPI_WIN_NULL)
        MPI_Win_free(&das->sm_comm_win);
    if (das->sm_comm != MPI_COMM_NULL)
        MPI_Comm_free(&das->sm_comm);
    if (das->sm_comm_ranks != NULL)
        free(das->sm_comm_ranks);
    if (das->node_comm != MPI_COMM_NULL)
        MPI_Comm_free(&das->node_comm);
    if (das->node_comm_ranks != NULL)
        free(das->node_comm_ranks);
    if (das->St != NULL)
        free(das->St);
#endif
    if (das->nregions > 0) {
        for (i = 0; i < das->nregions; ++i)
            free(das->regions[i].name);
        free(das->regions);
    }
    if (das->nplogs > 0) {
        ht_destroy(das->ht_plogs);
        for (i = 0; i < das->nplogs; ++i)
            free(das->plogs[i].gridname);
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
            if (nk > 1)
                f->level = k;
            else
                f->level = grid_getsurflayerid(model_getvargrid(m, vid));
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
        snprintf(fname, MAXSTRLEN, "%s.nc", FNAMEPREFIX_DIAG);
    else
        snprintf(fname, MAXSTRLEN, "%s-%d.nc", FNAMEPREFIX_DIAG, grid_getid(grid));
}

/**
 */
void das_getfname_plog(dasystem* das, pointlog* plog, char fname[])
{
    if (model_getngrid(das->m) == 1)
        snprintf(fname, MAXSTRLEN, "%s_%d,%d.nc", FNAMEPREFIX_PLOG, plog->i, plog->j);
    else
        snprintf(fname, MAXSTRLEN, "%s_%d,%d-%d.nc", FNAMEPREFIX_PLOG, plog->i, plog->j, plog->gridid);
}

/** Calculates the mixed layer depth for consistent application of the SST bias
 ** correction.
 */
void das_calcmld(dasystem* das, obstype* ot, float*** src, float** dst)
{
    model* m = das->m;
    int mvid = model_getvarid(m, ot->varnames[0], 1);
    int** nlevels = model_getnumlevels(m, mvid);
    double threshold = ot->mld_threshold;
    int ni, nj, nk, ksurf;
    int i, j, k, kk;

    model_getvardims(m, mvid, &ni, &nj, &nk);
    ksurf = grid_getsurflayerid(model_getvargrid(m, mvid));

    for (j = 0; j < nj; ++j) {
        for (i = 0; i < ni; ++i) {
            double vsurf, vprev, vnow;
            int kprev;
            double fk = NAN, z;

            if (nlevels[j][i] == 0) {
                dst[j][i] = NAN;
                continue;
            }

            vsurf = src[ksurf][j][i];
            vprev = vsurf;
            kprev = ksurf;
            k = ksurf;
            for (kk = 1; kk < nlevels[j][i]; ++kk) {
                k = (ksurf == 0) ? kk : ksurf - kk;
                vnow = src[k][j][i];
                if (fabs(vnow - vsurf) >= threshold) {
                    fk = (ksurf == 0) ? (double) kprev + (threshold - fabs(vprev - vsurf)) / fabs(vnow - vprev) : (double) kprev - (threshold - fabs(vprev - vsurf)) / fabs(vnow - vprev);
                    break;
                }
                kprev = k;
                vprev = vnow;
            }
            if (kk == nlevels[j][i])
                fk = (double) k;
            (void) model_fk2z(m, mvid, i, j, fk, &z);
            dst[j][i] = (float) z;
        }
    }
}
