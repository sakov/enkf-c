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
    int nmem_dynamic;
    int nmem;

    assert(das->mode != MODE_NONE);

    if (das->mode == MODE_ENOI && enkf_fstatsonly) {
        das->nmem = 1;
        return;
    }

    assert(das->ensdir != NULL);
    if (das->mode == MODE_HYBRID)
        assert(das->ensdir2 != NULL);

    /*
     * set das->nmem_dynamic to -1 to enable scanning of das->ensdir
     */
    nmem_dynamic = das->nmem_dynamic;
    das->nmem_dynamic = -1;

    nmem = 0;
    while (1) {
        int i;

        for (i = 0; i < nvar; ++i) {
            das_getmemberfname(das, model_getvarname(m, i), nmem + 1, fname);
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

    if (das->mode == MODE_ENOI || das->mode == MODE_ENKF) {
        if (das->nmem > 0) {
            if (nmem < das->nmem)
                enkf_quit("das_setnmem(): could not find \"%s\"", fname);
            if (nmem > das->nmem)
                enkf_printf("    %d members found on disk; ignoring excess to specified ensemble size = %d\n", nmem, das->nmem);
        } else
            das->nmem = nmem;
        if (nmem == 1)
            enkf_quit("only 1 member found; need at least 2 members to continue");
        if (das->mode == MODE_ENKF) {
            das->nmem_dynamic = das->nmem;
            das->nmem_static = 0;
        }
        if (das->mode == MODE_ENOI) {
            das->nmem_static = das->nmem;
            das->nmem_dynamic = 0;
        }
        return;
    }

    /*
     * das->mode = MODE_HYBRID
     */
    if (nmem_dynamic >= 0 && nmem > nmem_dynamic)
        enkf_printf("    %d dynamic members found on disk; ignoring excess to specified dynamic ensemble size = %d\n", nmem, nmem_dynamic);
    if (nmem_dynamic >= 0 && nmem < nmem_dynamic)
        enkf_quit("das_setnmem(): could not find \"%s\"", fname);
    das->nmem_dynamic = (nmem_dynamic < 0) ? nmem : nmem_dynamic;
    nmem = das->nmem_dynamic;

    while (1) {
        int i;

        for (i = 0; i < nvar; ++i) {
            das_getmemberfname(das, model_getvarname(m, i), nmem + 1, fname);
            if (!file_exists(fname))
                break;
        }
        if (i == nvar)
            nmem++;
        else
            break;
    }
    if (nmem == das->nmem_dynamic)
        enkf_quit("das_setnmem(): could not find \"%s\"", fname);
    if (das->nmem_static >= 0 && nmem > das->nmem_dynamic + das->nmem_static)
        enkf_printf("    %d static members found on disk; ignoring excess to specified static ensemble size = %d\n", nmem - das->nmem_dynamic, das->nmem_static);
    if (das->nmem_static >= 0 && nmem < das->nmem_dynamic + das->nmem_static)
        enkf_quit("das_setnmem(): could not find \"%s\"", fname);
    if (das->nmem_static < 0)
        das->nmem_static = nmem - das->nmem_dynamic;
    das->nmem = das->nmem_dynamic + das->nmem_static;

    if (das->nmem == 1)
        enkf_quit("only 1 member found; need at least 2 members to continue");
    if (das->nmem_dynamic == 1)
        enkf_printf("    only 1 dynamic member found; effectively the system will run in EnOI mode");
    if (das->nmem_static == 1)
        enkf_quit("only 1 static member found; need at least 2 to continue");
}

#if defined(ENKF_UPDATE)
/**
 */
static void get_gridstr(dasystem* das, int gid, char str[])
{
    if (model_getngrid(das->m) == 1)
        str[0] = 0;
    else
        sprintf(str, "-%d", gid);
}
#endif

/**
 */
#if defined(ENKF_UPDATE)
dasystem* das_create(enkfprm* prm, int updatespec)
#else
dasystem* das_create(enkfprm* prm)
#endif
{
    dasystem* das = calloc(1, sizeof(dasystem));

#if defined(ENKF_UPDATE)
    das->updatespec = updatespec;
#endif
    das->prmfname = strdup(prm->fname);
    das->mode = prm->mode;
    das->scheme = prm->scheme;
    if (!(das->mode == MODE_ENOI && enkf_fstatsonly))
        das->ensdir = strdup(prm->ensdir);
    if (prm->ensdir2 != NULL)
        das->ensdir2 = strdup(prm->ensdir2);
    if (prm->bgdir != NULL)
        das->bgdir = strdup(prm->bgdir);
    das->alpha = prm->alpha;
    das->nmem = prm->enssize;
    if (das->mode == MODE_HYBRID) {
        das->nmem_dynamic = prm->enssize_dynamic;
        das->nmem_static = prm->enssize_static;
        das->gamma = prm->gamma;
    }
#if defined(ENKF_CALC)
    das->obs = obs_create_fromprm(prm);
    if (das->strict_time_matching) {
        int otid;

        for (otid = 0; otid < das->obs->nobstypes; ++otid) {
            obstype* ot = &das->obs->obstypes[otid];

            if (ot->isasync && ot->async_tname == NULL)
                enkf_quit("%s: %s: time variable name must be defined for asynchronously assimilated observation types with \"--strict-time-matching\"; see description of entry ASYNC by \"enkf_calc --describe-prm-format obstypes\"", prm->obstypeprm, ot->name);
        }
    }
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
    if (das->mode != MODE_HYBRID)
        enkf_printf("    %d member%s\n", das->nmem, das->nmem == 1 ? "" : "s");
    else {
        enkf_printf("    %d dynamic member%s\n", das->nmem_dynamic, das->nmem == 1 ? "" : "s");
        enkf_printf("    %d static member%s\n", das->nmem_static, das->nmem_static == 1 ? "" : "s");
    }

#if defined(ENKF_CALC)
    if (das->nmem == 1 && !(das->mode == MODE_ENOI && enkf_fstatsonly))
        enkf_quit("CALC is not supposed to be run with 1-member ensemble");
#endif

#if defined(ENKF_CALC)
    if (!enkf_fstatsonly) {
        das->kfactor = prm->kfactor;
    } else {
        das->kfactor = NAN;
    }
#endif
#if defined(USE_SHMEM)
    das->sm_comm_win_S = MPI_WIN_NULL;
    das->sm_comm_win_St = MPI_WIN_NULL;
    das->S = NULL;
    das->St = NULL;
#endif
    if (!enkf_fstatsonly)
        das->fieldbufsize = prm->fieldbufsize;

    /*
     * initialise regions
     */
#if defined(ENKF_CALC)
    das->nregions = prm->nregions;
    if (das->nregions > 0) {
        int i;

        das->regions = malloc(das->nregions * sizeof(region));
        for (i = 0; i < das->nregions; ++i) {
            region* dst = &das->regions[i];
            region* src = &prm->regions[i];

            dst->name = strdup(src->name);
            dst->x1 = src->x1;
            dst->x2 = src->x2;
            dst->y1 = src->y1;
            dst->y2 = src->y2;
        }
    } else
        das->regions = NULL;
#endif

    /*
     * initialise pointlogs
     */
#if defined(ENKF_CALC) || defined(ENKF_UPDATE)
    if (enkf_doplogs && prm->nplog > 0) {
        int ngrid, i;

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

            dst->fij = alloc2d(ngrid, 3, sizeof(double));
            for (gid = 0; gid < model_getngrid(das->m); ++gid) {
                grid* g = model_getgridbyid(das->m, gid);

                if (grid_isempty(g) || (dst->gridid >= 0 && dst->gridid != gid)) {
                    dst->fij[gid][0] = NAN;
                    dst->fij[gid][1] = NAN;
                    dst->fij[gid][2] = NAN;
                    continue;
                }
#if defined(ENKF_CALC)
                {
                    double fij[3] = { NAN, NAN, NAN };

                    if (grid_xy2fij(g, src->lon, src->lat, fij) != STATUS_OK && gid == dst->gridid)
                        enkf_printf("  WARNING: %s: POINTLOG %f %f: point outside the grid \"%s\"\n", das->prmfname, dst->lon, dst->lat, dst->gridname);
                    dst->fij[gid][0] = fij[0];
                    dst->fij[gid][1] = fij[1];
                    dst->fij[gid][2] = fij[2];
                }
#elif defined(ENKF_UPDATE)
                {
                    char fname[MAXSTRLEN];
                    int ncid, varid;
                    char gridstr[SHORTSTRLEN];
                    char varname[NC_MAX_NAME];

                    das_getfname_plog(das, dst, fname);
                    if (!file_exists(fname))
                        enkf_quit("file \"%s\" not found", fname);

                    ncw_open(fname, NC_NOWRITE, &ncid);
                    get_gridstr(das, gid, gridstr);
                    snprintf(varname, NC_MAX_NAME, "grid%s", gridstr);
                    ncw_inq_varid(ncid, varname, &varid);
                    if (ncw_att_exists(ncid, varid, "fi")) {
                        ncw_get_att_double(ncid, varid, "fi", &dst->fij[gid][0]);
                        ncw_get_att_double(ncid, varid, "fj", &dst->fij[gid][1]);
                    } else {
                        ncw_get_att_double(ncid, varid, "fi0", &dst->fij[gid][0]);
                        ncw_get_att_double(ncid, varid, "fi1", &dst->fij[gid][1]);
                    }
                    ncw_close(ncid);
                    if (isnan(dst->fij[gid][0] + dst->fij[gid][1]) && gid == dst->gridid)
                        enkf_printf("  WARNING: %s: POINTLOG %f %f: point outside the grid \"%s\"\n", das->prmfname, dst->lon, dst->lat, dst->gridname);
                }
#else
                enkf_quit("programming error");
#endif
                if (grid_isstructured(g))
                    enkf_printf("      %s: (i, j) = (%.3f, %.3f)\n", grid_getname(g), dst->fij[gid][0], dst->fij[gid][1]);
                else
                    enkf_printf("      %s: (i0, i1, i2) = (%.3f, %.3f, %.3f)\n", grid_getname(g), dst->fij[gid][0], dst->fij[gid][1], dst->fij[gid][2]);
            }
            for (gid = 0; gid < model_getngrid(das->m); ++gid)
                if (!isnan(dst->fij[gid][0] + dst->fij[gid][1]))
                    break;
            if (gid == model_getngrid(das->m))
                enkf_quit("%s: POINTLOG %f %f: point outside all grids\n", das->prmfname, dst->lon, dst->lat);
        }
    }
#endif

    /*
     * initialise badbatches
     */
#if defined(ENKF_CALC)
    das->nbadbatchspecs = prm->nbadbatchspecs;
    if (das->nbadbatchspecs > 0) {
        int i;

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

    das->ncformat = prm->ncformat;
    das->nccompression = prm->nccompression;

    return das;
}

/**
 */
void das_destroy(dasystem* das)
{
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
    }
    if (das->sm_comm_win_St != MPI_WIN_NULL) {
        MPI_Win_free(&das->sm_comm_win_St);
        assert(das->sm_comm_win_St == MPI_WIN_NULL);
    }
    if (das->St != NULL)
        free(das->St);
#endif
    if (das->nregions > 0) {
        int i;

        for (i = 0; i < das->nregions; ++i)
            free(das->regions[i].name);
        free(das->regions);
    }
#if defined(ENKF_CALC) || defined(ENKF_UPDATE)
    plogs_destroy(das->nplog, das->plogs);
#endif
    if (das->nbadbatchspecs > 0) {
        int i;

        for (i = 0; i < das->nbadbatchspecs; ++i)
            free(das->badbatchspecs[i].obstype);

        free(das->badbatchspecs);
    }
    free(das);
    distribute_free();
    if (rank == 0)
        dir_rmifexists(DIRNAME_TMP);
}

/** Looks for all horizontal fields of the model on a specified grid to be
 ** updated.
 * @param das - das structure
 * @param gridid - grid ID (< 0 = all grids)
 * @param nfields - (out) number of fields
 * @param fields - (out) ields
 */
void das_getfields(dasystem* das, int gridid, int* nfields, field** fields)
{
    model* m = das->m;
    int nvar = model_getnvar(m);
    int vid;

    assert(*nfields == 0);
    assert(*fields == NULL);

    if (gridid >= 0 && grid_getaliasid(model_getgridbyid(m, gridid)) >= 0)
        return;

    for (vid = 0; vid < nvar; ++vid) {
        int vargridid = model_getvargridid(m, vid);
        grid* g = model_getgridbyid(m, vargridid);
        int varaliasid = grid_getaliasid(g);
        char* varname = model_getvarname(m, vid);
        char fname[MAXSTRLEN];
        int nj;
        int nk, nkgrid, k;

        if (gridid >= 0 && vargridid != gridid && varaliasid != gridid)
            continue;

        grid_getsize(g, NULL, &nj, &nkgrid);

        das_getmemberfname(das, varname, 1, fname);
        nk = ncu_getnlevels(fname, varname, nj > 0);
        if (nk > 1 && nk != nkgrid)
            enkf_quit("%s: %s: variable nk = %d, grid nk = %d: the number of levels for a variable is supposed to be either that of the grid or 1", fname, grid_getname(g), nk, nkgrid);
        for (k = 0; k < nk; ++k) {
            field* f;

            if (*nfields % NFIELDS_INC == 0)
                *fields = realloc(*fields, (*nfields + NFIELDS_INC) * sizeof(field));
            f = &(*fields)[*nfields];
            f->id = *nfields;
            f->varid = vid;
            f->issurfacevar = (nk == 1);
            strcpy(f->varname, varname);
            f->level = k;
            f->structured = nj > 0;
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
void das_getfname_transforms(dasystem* das, int gridid, char fname[])
{
    if (model_getngrid(das->m) == 1)
        snprintf(fname, MAXSTRLEN, "%s.nc", FNAMEPREFIX_TRANSFORMS);
    else {
        int aliasid = grid_getaliasid(model_getgridbyid(das->m, gridid));

        snprintf(fname, MAXSTRLEN, "%s-%d.nc", FNAMEPREFIX_TRANSFORMS, (aliasid >= 0) ? aliasid : gridid);
    }
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
#if defined(ENKF_PREP) || defined(ENKF_CALC)
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
            {
                int ij[2] = { i, j };

                z = model_fk2z(m, mvid, ij, fk);
            }
            dst[j][i] = (float) z;
        }
    }
}
#endif

/**
 */
void das_getmemberfname(dasystem* das, char varname[], int mem, char fname[])
{
    if (das->mode == MODE_HYBRID && das->nmem_dynamic >= 0 && mem > das->nmem_dynamic)
        snprintf(fname, MAXSTRLEN, "%s/mem%03d_%s.nc", das->ensdir2, mem - das->nmem_dynamic, varname);
    else
        snprintf(fname, MAXSTRLEN, "%s/mem%03d_%s.nc", das->ensdir, mem, varname);
}

/**
 */
void das_getbgfname(dasystem* das, char varname[], char fname[])
{
    snprintf(fname, MAXSTRLEN, "%s/bg_%s.nc", das->bgdir, varname);
}

#if defined(ENKF_CALC)
/**
 */
int das_getmemberfname_async(dasystem* das, obstype* ot, int mem, int t, char fname[], int* r)
{
    char* alias = ot->alias;
    char* varname = ot->varnames[0];
    char* ensdir = das->ensdir;

    *r = -1;

    if (das->mode == MODE_HYBRID && das->nmem_dynamic >= 0 && mem > das->nmem_dynamic) {
        snprintf(fname, MAXSTRLEN, "%s/mem%03d_%s.nc", das->ensdir2, mem - das->nmem_dynamic, varname);
        return 0;
    }

    snprintf(fname, MAXSTRLEN, "%s/mem%03d_%s_%d.nc", ensdir, mem, alias, t);
    if (!file_exists(fname)) {
        char fname2[MAXSTRLEN];

        snprintf(fname2, MAXSTRLEN, "%s/mem%03d_%s_#.nc", ensdir, mem, alias);
        if (!file_exists(fname2)) {
            if (das->strict_time_matching)
                enkf_quit("could not find file \"%s\" or \"%s\", which is necessary to proceed because (1) asynchronous DA is set on for \"%s\" and (2) \"--strict-time-matching\" is used\n", fname, fname2, ot->name);
            snprintf(fname, MAXSTRLEN, "%s/mem%03d_%s.nc", ensdir, mem, varname);
            return 0;
        } else
            strncpy(fname, fname2, MAXSTRLEN);
    }

    /*
     * verify time
     */
    if (ot->async_tname != NULL) {
        int ncid, vid;
        size_t vsize;
        double* time;
        char tunits[MAXSTRLEN];
        size_t attlen;
        double tunits_multiple, tunits_offset;
        double correcttime;
        int i;

        ncw_open(fname, NC_NOWRITE, &ncid);
        if (!ncw_var_exists(ncid, ot->async_tname))
            enkf_quit("%s: found no time variable \"%s\" specified for observation type \"%s\"", fname, ot->async_tname, ot->name);
        ncw_inq_varid(ncid, ot->async_tname, &vid);
        ncw_check_varndims(ncid, vid, 1);
        ncw_inq_varsize(ncid, vid, &vsize);

        time = malloc(vsize * sizeof(double));
        ncw_get_var_double(ncid, vid, time);

        ncw_inq_attlen(ncid, vid, "units", &attlen);
        assert(attlen <= MAXSTRLEN);
        ncw_get_att_text(ncid, vid, "units", tunits);
        tunits_convert(tunits, &tunits_multiple, &tunits_offset);

        correcttime = das->obs->da_time + t * ot->async_tstep;
        if (!ot->async_centred)
            correcttime += 0.5 * ot->async_tstep;

        for (i = 0; i < vsize; ++i) {
            time[i] = time[i] * tunits_multiple + tunits_offset;
            if (fabs(time[i] - correcttime) < TEPS)
                break;
        }
        if (i >= vsize) {
            if (vsize == 1)
                enkf_quit("%s: \"s\" = %f; expected %f\n", fname, ot->async_tname, time, correcttime);
            else
                enkf_quit("%s: time variable \"%s\" has no matching value for asynchronous interval %d (needed time = %f)\n", fname, ot->async_tname, t, correcttime);
        }
        free(time);
        *r = i;
    }
    return 1;
}

/**
 */
int das_getbgfname_async(dasystem* das, obstype* ot, int t, char fname[], int* r)
{
    char* alias = ot->alias;
    char* varname = ot->varnames[0];
    char* bgdir = das->bgdir;

    *r = -1;

    snprintf(fname, MAXSTRLEN, "%s/bg_%s_%d.nc", bgdir, alias, t);
    if (!file_exists(fname)) {
        char fname2[MAXSTRLEN];

        snprintf(fname2, MAXSTRLEN, "%s/bg_%s_#.nc", bgdir, alias);
        if (!file_exists(fname2)) {
            if (das->strict_time_matching)
                enkf_quit("could not find file \"%s\" or \"%s\", which is necessary to proceed because (1) asynchronous DA is set on for \"%s\" and (2) \"--strict-time-matching\" is used\n", fname, fname2, ot->name);
            snprintf(fname, MAXSTRLEN, "%s/bg_%s.nc", bgdir, varname);
            return 0;
        } else
            strncpy(fname, fname2, MAXSTRLEN);
    }
    /*
     * verify time
     */
    if (ot->async_tname != NULL) {
        int ncid, vid;
        size_t vsize;
        double* time;
        char tunits[MAXSTRLEN];
        size_t attlen;
        double tunits_multiple, tunits_offset;
        double correcttime;
        int i;

        ncw_open(fname, NC_NOWRITE, &ncid);
        if (!ncw_var_exists(ncid, ot->async_tname))
            enkf_quit("%s: found no time variable \"%s\" specified for observation type \"%s\"", fname, ot->async_tname, ot->name);
        ncw_inq_varid(ncid, ot->async_tname, &vid);
        ncw_check_varndims(ncid, vid, 1);
        ncw_inq_varsize(ncid, vid, &vsize);

        time = malloc(vsize * sizeof(double));
        ncw_get_var_double(ncid, vid, time);

        ncw_inq_attlen(ncid, vid, "units", &attlen);
        assert(attlen <= MAXSTRLEN);
        ncw_get_att_text(ncid, vid, "units", tunits);
        tunits_convert(tunits, &tunits_multiple, &tunits_offset);

        correcttime = das->obs->da_time + t * ot->async_tstep;
        if (!ot->async_centred)
            correcttime += 0.5 * ot->async_tstep;

        for (i = 0; i < vsize; ++i) {
            time[i] = time[i] * tunits_multiple + tunits_offset;
            if (fabs(time[i] - correcttime) < TEPS)
                break;
        }
        if (i >= vsize) {
            if (vsize == 1)
                enkf_quit("%s: \"s\" = %f; expected %f\n", fname, ot->async_tname, time, correcttime);
            else
                enkf_quit("%s: time variable \"%s\" has no matching value for asynchronous interval %d (needed time = %f)\n", fname, ot->async_tname, t, correcttime);
        }
        free(time);
        *r = i;
    }
    return 1;
}
#endif

/** Calculate dynamic ensemble mean; add it to static ensemble anomalies; scale
 ** static ensemble anomalies.
 *  Note that the allocated size of v should be v[das->nmem + 2][nij] to make it
 *  possible calculating ensemble mean with double precision.
 */
void das_sethybridensemble(dasystem* das, int nij, float** v)
{
    double* vmean = (double*) v[das->nmem];
    int nmem = das->nmem;
    int nmem_d = das->nmem_dynamic;
    int nmem_s = das->nmem_static;
    double k_d = (nmem_d > 1) ? sqrt((double) (nmem - 1) / (double) (nmem_d - 1)) : 0.0;
    double k_s = sqrt(das->gamma * (double) (nmem - 1) / (double) (nmem_s - 1));
    int i, e;

    assert(das->mode == MODE_HYBRID);

    /*
     * calculate dynamic ensemble mean
     */
    memset(vmean, 0, nij * sizeof(double));
    for (e = 0; e < nmem_d; ++e) {
        float* ve = v[e];

        for (i = 0; i < nij; ++i)
            vmean[i] += (double) ve[i];
    }
    for (i = 0; i < nij; ++i)
        vmean[i] /= (double) nmem_d;

    /*
     * set dynamic members
     */
    for (e = 0; e < nmem_d; ++e) {
        float* ve = v[e];

        for (i = 0; i < nij; ++i)
            ve[i] = (ve[i] - vmean[i]) * k_d + vmean[i];
    }

    /*
     * set static members
     */
    for (e = nmem_d; e < nmem; ++e) {
        float* ve = v[e];

        for (i = 0; i < nij; ++i)
            ve[i] = ve[i] * k_s + vmean[i];
    }
}

/** Return 1 if this member is static anomaly, 0 otherwise.
 * @param das - pointer to dasystem
 * @param mem - member number, 1 <= mem <= ensemble size, and mem = -1 is used
 *              for the background
 */
int das_isstatic(dasystem* das, int mem)
{
    return mem > das->nmem_dynamic;
}
