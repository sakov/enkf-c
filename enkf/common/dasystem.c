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
 * Revisions:   06032020 PS: moved code for setting MPI communicators etc.
 *                from das_create() to enkf_init()
 *
 *****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <unistd.h>
#include <limits.h>
#include "ncw.h"
#include "hash.h"
#include "definitions.h"
#include "utils.h"
#include "ncutils.h"
#include "distribute.h"
#include "enkfprm.h"
#include "observations.h"
#include "dasystem.h"
#include "pointlog.h"

#define NPLOGS_INC 10
#define NFIELDS_INC 100
#define MPIIDOFFSET 10000
#define TEPS 1.0e-3

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

    assert(das->mode != MODE_NONE);

    if (das->mode == MODE_ENOI && enkf_fstatsonly) {
        das->nmem = 1;
        return;
    }

    assert(das->ensdir != NULL);

    nmem = 0;
    while (1) {
        int i;

        for (i = 0; i < nvar; ++i) {
            das_getmemberfname(das, das->ensdir, model_getvarname(m, i), nmem + 1, fname);
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
        else if (nmem > das->nmem) {
            enkf_printf("      %d members found on disk; ignoring excess to specified ensemble size = %d\n", nmem, das->nmem);
            return;
        }
    }
    das->nmem = nmem;
}

#if defined(ENKF_UPDATE)/**
 */
static void get_gridstr(dasystem* das, int gid, char str[])
{
    if (model_getngrid(das->m) == 1)
        strcpy(str, "");
    else
        sprintf(str, "-%d", gid);
}
#endif

/**
 */
dasystem* das_create(enkfprm* prm)
{
    dasystem* das = calloc(1, sizeof(dasystem));
    int ngrid;
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

#if defined(ENKF_CALC)
    if (!enkf_fstatsonly) {
        das->kfactor = prm->kfactor;
    } else {
        das->kfactor = NAN;
    }
#endif
#if defined(USE_SHMEM)
    das->sm_comm_win_S = MPI_WIN_NULL;
    das->St = NULL;
#endif
    if (!enkf_fstatsonly)
        das->fieldbufsize = prm->fieldbufsize;

    /*
     * initialise regions
     */
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

    /*
     * initialise pointlogs
     */
    enkf_printf("  initialising pointlogs:\n");
    das->plogs = malloc(sizeof(pointlog) * prm->nplog);
    das->nplog = prm->nplog;
    ngrid = model_getngrid(das->m);
    for (i = 0; i < prm->nplog; ++i) {
        pointlog* src = &prm->plogs[i];
        pointlog* dst = &das->plogs[i];
        int gid;

        enkf_printf("    pointlog (%.3f, %.3f):\n", src->lon, src->lat);

        dst->id = src->id;
        dst->lon = src->lon;
        dst->lat = src->lat;
        if (src->gridname != NULL) {
            dst->gridname = strdup(src->gridname);
            dst->gridid = grid_getid(model_getgridbyname(das->m, dst->gridname));
        } else {
            dst->gridname = NULL;
            dst->gridid = -1;
        }

        dst->fi = malloc(ngrid * sizeof(double));
        dst->fj = malloc(ngrid * sizeof(double));
        for (gid = 0; gid < model_getngrid(das->m); ++gid) {
            grid* g = model_getgridbyid(das->m, gid);

            if (dst->gridid >= 0 && dst->gridid != gid) {
                dst->fi[gid] = NAN;
                dst->fj[gid] = NAN;
                continue;
            }
#if defined(ENKF_CALC)
            if (grid_xy2fij(g, src->lon, src->lat, &dst->fi[gid], &dst->fj[gid]) != STATUS_OK && gid == dst->gridid)
                enkf_printf("  WARNING: %s: POINTLOG %f %f: point outside the grid \"%s\"\n", das->prmfname, dst->lon, dst->lat, dst->gridname);
#elif defined(ENKF_UPDATE)
            {
                char fname[MAXSTRLEN];
                int ncid, varid;
                char gridstr[SMALLSTRLEN];
                char varname[NC_MAX_NAME];

                das_getfname_plog(das, dst, fname);
                assert(file_exists(fname));

                ncw_open(fname, NC_NOWRITE, &ncid);
                get_gridstr(das, gid, gridstr);
                snprintf(varname, NC_MAX_NAME, "grid%s", gridstr);
                ncw_inq_varid(ncid, varname, &varid);
                ncw_get_att_double(ncid, varid, "fi", &dst->fi[gid]);
                ncw_get_att_double(ncid, varid, "fj", &dst->fj[gid]);
                ncw_close(ncid);
                if (isnan(dst->fi[gid] + dst->fj[gid]) && gid == dst->gridid)
                    enkf_printf("  WARNING: %s: POINTLOG %f %f: point outside the grid \"%s\"\n", das->prmfname, dst->lon, dst->lat, dst->gridname);
            }
#else
            enkf_quit("programming error");
#endif
            enkf_printf("      %s: (i, j) = (%.3f, %.3f)\n", grid_getname(g), dst->fi[gid], dst->fj[gid]);
        }
        for (gid = 0; gid < model_getngrid(das->m); ++gid)
            if (!isnan(dst->fi[gid] + dst->fj[gid]))
                break;
        if (gid == model_getngrid(das->m))
            enkf_quit("%s: POINTLOG %f %f: point outside all grids\n", das->prmfname, dst->lon, dst->lat);
    }

    /*
     * initialise badbatches
     */
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
    model_destroy(das->m);
#if defined(ENKF_CALC)
    obs_destroy(das->obs);
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
#if defined (USE_SHMEM)
    if (das->sm_comm_win_S != MPI_WIN_NULL) {
        MPI_Win_free(&das->sm_comm_win_S);
        assert(das->sm_comm_win_S == MPI_WIN_NULL);
        MPI_Win_free(&das->sm_comm_win_St);
        assert(das->sm_comm_win_St == MPI_WIN_NULL);
    }
    if (das->St != NULL)
        free(das->St);
#endif
    if (das->nregions > 0) {
        for (i = 0; i < das->nregions; ++i)
            free(das->regions[i].name);
        free(das->regions);
    }
    plogs_destroy(das->nplog, das->plogs);
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

        das_getmemberfname(das, das->ensdir, varname, 1, fname);
        nk = ncu_getnlevels(fname, varname);
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
void getfieldfname(char* dir, char* prefix, char* varname, int level, char* fname)
{
    snprintf(fname, MAXSTRLEN, "%s/%s_%s-%03d.nc", dir, prefix, varname, level);
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
    snprintf(fname, MAXSTRLEN, "%s_%.3f,%.3f.nc", FNAMEPREFIX_PLOG, plog->lon, plog->lat);
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

    model_getvargridsize(m, mvid, &ni, &nj, &nk);
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

/**
 */
void das_getmemberfname(dasystem* das, char ensdir[], char varname[], int mem, char fname[])
{
    snprintf(fname, MAXSTRLEN, "%s/mem%03d_%s.nc", ensdir, mem, varname);
}

/**
 */
void das_getbgfname(dasystem* das, char ensdir[], char varname[], char fname[])
{
    snprintf(fname, MAXSTRLEN, "%s/bg_%s.nc", ensdir, varname);
}

#if defined(ENKF_CALC)
/**
 */
int das_getmemberfname_async(dasystem* das, char ensdir[], obstype* ot, int mem, int t, char fname[])
{
    char* alias = ot->alias;
    char* varname = ot->varnames[0];

    snprintf(fname, MAXSTRLEN, "%s/mem%03d_%s_%d.nc", ensdir, mem, alias, t);
    if (!file_exists(fname)) {
        snprintf(fname, MAXSTRLEN, "%s/mem%03d_%s.nc", ensdir, mem, varname);
        return 0;
    }
    /*
     * if the time variable name has been specified for the obs. type -- verify
     * that the time is right for the interval t
     */
    if (ot->async_tname != NULL) {
        int ncid, vid;
        size_t vsize;
        double time, correcttime;

        ncw_open(fname, NC_NOWRITE, &ncid);
        if (!ncw_var_exists(ncid, ot->async_tname))
            enkf_quit("%s: found no time variable \"%s\" specified for observation type \"%s\"", fname, ot->async_tname, ot->name);
        ncw_inq_varid(ncid, ot->async_tname, &vid);
        ncw_inq_varsize(ncid, vid, &vsize);
        if (vsize != 1)
            enkf_quit("%s: dimension of the time variable \"%s\" must be one", fname, ot->async_tname);
        {
            char tunits[MAXSTRLEN];
            double tunits_multiple, tunits_offset;

            ncw_get_var_double(ncid, vid, &time);
            ncw_get_att_text(ncid, vid, "units", tunits);
            tunits_convert(tunits, &tunits_multiple, &tunits_offset);
            time = time * tunits_multiple + tunits_offset;
            correcttime = das->obs->da_time + t * ot->async_tstep;
            if (!ot->async_centred)
                correcttime += 0.5 * ot->async_tstep;
            if (fabs(time - correcttime) > TEPS)
                enkf_quit("%s: \"s\" = %f; expected %f\n", fname, ot->async_tname, time, correcttime);
        }
    }
    return 1;
}

/**
 */
int das_getbgfname_async(dasystem* das, char bgdir[], obstype* ot, int t, char fname[])
{
    char* alias = ot->alias;
    char* varname = ot->varnames[0];

    snprintf(fname, MAXSTRLEN, "%s/bg_%s_%d.nc", bgdir, alias, t);
    if (!file_exists(fname)) {
        snprintf(fname, MAXSTRLEN, "%s/bg_%s.nc", bgdir, varname);
        return 0;
    }
    /*
     * if the time variable name has been specified for the obs. type -- verify
     * that the time is right for the interval t
     */
    if (ot->async_tname != NULL) {
        int ncid, vid;
        size_t vsize;
        double time, correcttime;

        ncw_open(fname, NC_NOWRITE, &ncid);
        if (!ncw_var_exists(ncid, ot->async_tname))
            enkf_quit("%s: found no time variable \"%s\" specified for observation type \"%s\"", fname, ot->async_tname, ot->name);
        ncw_inq_varid(ncid, ot->async_tname, &vid);
        ncw_inq_varsize(ncid, vid, &vsize);
        if (vsize != 1)
            enkf_quit("%s: dimension of the time variable \"%s\" must be one", fname, ot->async_tname);
        {
            char tunits[MAXSTRLEN];
            double tunits_multiple, tunits_offset;

            ncw_get_var_double(ncid, vid, &time);
            ncw_get_att_text(ncid, vid, "units", tunits);
            tunits_convert(tunits, &tunits_multiple, &tunits_offset);
            time = time * tunits_multiple + tunits_offset;
            correcttime = das->obs->da_time + t * ot->async_tstep;
            if (!ot->async_centred)
                correcttime += 0.5 * ot->async_tstep;
            if (fabs(time - correcttime) > TEPS)
                enkf_quit("%s: \"s\" = %f; expected %f\n", fname, ot->async_tname, time, correcttime);
        }
    }
    return 1;
}
#endif
