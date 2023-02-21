/******************************************************************************
 *
 * File:        diags.c
 *
 * Created:     15/07/2019
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description: Moved here code from update.c for calculating and writing a
 *              number of diagnostic quantities:
 *                - spread
 *                - inflation
 *                - ensemble correlations with the surface layer
 *
 * Revisions:
 *
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "ncw.h"
#include "distribute.h"
#include "definitions.h"
#include "utils.h"
#include "ncutils.h"
#include "model.h"
#include "dasystem.h"

/** Allocates disk space for ensemble spread.
 */
void das_allocatespread(dasystem* das, char fname[])
{
    model* m = das->m;
    int nvar = model_getnvar(m);
    char fname_src[MAXSTRLEN];
    int ncid, ncid_src;
    int vid;

    if (rank != 0)
        return;

    ncw_create(fname, NC_CLOBBER | das->ncformat, &ncid);
    for (vid = 0; vid < nvar; ++vid) {
        char* varname_src = model_getvarname(m, vid);
        int varid_src;

        das_getmemberfname(das, varname_src, 1, fname_src);
        ncw_open(fname_src, NC_NOWRITE, &ncid_src);
        ncw_inq_varid(ncid_src, varname_src, &varid_src);
        ncw_copy_vardef(ncid_src, varid_src, ncid);
        if ((das->mode == MODE_ENKF || das->mode == MODE_HYBRID) && das->updatespec & UPDATE_DOANALYSISSPREAD) {
            char varname_dst[NC_MAX_NAME];

            strcpy(varname_dst, varname_src);
            if (das->updatespec & UPDATE_OUTPUTINC)
                strncat(varname_dst, "_inc", NC_MAX_NAME - 1);
            else
                strncat(varname_dst, "_an", NC_MAX_NAME - 1);
            ncw_def_var_as(ncid, varname_src, varname_dst);
        }
        ncw_close(ncid_src);
    }
#if defined(DEFLATE_ALL)
    if (das->nccompression > 0)
        ncw_def_deflate(ncid, 0, 1, das->nccompression);
#endif
    ncw_close(ncid);
}

/**
 */
void das_writespread(dasystem* das, int nfields, void** fieldbuffer, field fields[], int isanalysis)
{
    model* m = das->m;
    int nmem = (das->mode == MODE_HYBRID && isanalysis) ? das->nmem_dynamic : das->nmem;
    int ni, nj;
    int fid, e, i, nij;
    double* v1 = NULL;
    double* v2 = NULL;

    model_getvargridsize(m, fields[0].varid, &ni, &nj, NULL);
    nij = ni * nj;
    v1 = malloc(nij * sizeof(double));
    v2 = malloc(nij * sizeof(double));

    for (fid = 0; fid < nfields; ++fid) {
        field* f = &fields[fid];
        float*** v_src = (float***) fieldbuffer[fid];
        char varname[NC_MAX_NAME];

        memset(v1, 0, nij * sizeof(double));
        memset(v2, 0, nij * sizeof(double));

        for (e = 0; e < nmem; ++e) {
            float* v = v_src[e][0];

            for (i = 0; i < nij; ++i) {
                v1[i] += v[i];
                v2[i] += v[i] * v[i];
            }
        }

        for (i = 0; i < nij; ++i) {
            v1[i] /= (double) nmem;
            v2[i] = v2[i] / (double) nmem - v1[i] * v1[i];
            v2[i] = (v2[i] < 0.0) ? 0.0 : sqrt(v2[i]);
            if (fabs(v1[i]) > (double) MAXOBSVAL)
                v2[i] = 0.0;
        }

        strncpy(varname, f->varname, NC_MAX_NAME - 1);
        if (isanalysis) {
            if (das->updatespec & UPDATE_OUTPUTINC)
                strncat(varname, "_inc", NC_MAX_NAME - 1);
            else
                strncat(varname, "_an", NC_MAX_NAME - 1);
        }

        if (!(das->updatespec & UPDATE_DIRECTWRITE)) {
            char fname[MAXSTRLEN];
            int ncid, vid;
            int dimids[2];

            getfieldfname(DIRNAME_TMP, "spread", varname, f->level, fname);

            ncw_create(fname, NC_CLOBBER | das->ncformat, &ncid);
            ncw_def_dim(ncid, "nj", nj, &dimids[0]);
            ncw_def_dim(ncid, "ni", ni, &dimids[1]);
            ncw_def_var(ncid, varname, NC_FLOAT, 2, dimids, &vid);
#if defined(DEFLATE_ALL)
            if (das->nccompression > 0)
                ncw_def_deflate(ncid, 0, 1, das->nccompression);
#endif
            ncw_enddef(ncid);
            ncw_put_var_double(ncid, vid, v2);
            ncw_close(ncid);
        } else {
            float* v = calloc(nij, sizeof(float));

            for (i = 0; i < nij; ++i)
                v[i] = (float) v2[i];

            model_writefieldas(m, FNAME_SPREAD, varname, f->varname, f->level, v, 1);
            free(v);
        }
    }

    free(v1);
    free(v2);
}

/**
 */
void das_assemblespread(dasystem* das)
{
    model* m = das->m;
    int nvar = model_getnvar(m);
    int i;

    if (rank > 0)
        return;

    for (i = 0; i < nvar; ++i) {
        char* varname = model_getvarname(m, i);
        int nlev = ncu_getnlevels(FNAME_SPREAD, varname);
        int ni, nj;
        float* v;
        int k;

        model_getvargridsize(m, i, &ni, &nj, NULL);
        v = malloc(ni * nj * sizeof(float));

        enkf_printf("    %s:", varname);
        for (k = 0; k < nlev; ++k) {
            char fname_src[MAXSTRLEN];
            int ncid_src, vid;

            if (nlev > 1)
                getfieldfname(DIRNAME_TMP, "spread", varname, k, fname_src);
            else
                getfieldfname(DIRNAME_TMP, "spread", varname, grid_getsurflayerid(model_getvargrid(m, i)), fname_src);
            ncw_open(fname_src, NC_NOWRITE, &ncid_src);
            ncw_inq_varid(ncid_src, varname, &vid);
            ncw_get_var_float(ncid_src, vid, v);
            ncw_close(ncid_src);
            file_delete(fname_src);

            model_writefield(m, FNAME_SPREAD, varname, k, v, 1);
            enkf_printf(".");
            enkf_flush();
        }
        enkf_printf("\n");
        enkf_flush();

        if ((das->mode == MODE_ENKF || das->mode == MODE_HYBRID) && das->updatespec & UPDATE_DOANALYSISSPREAD) {
            char varname_an[NC_MAX_NAME];

            strncpy(varname_an, varname, NC_MAX_NAME - 1);
            if (das->updatespec & UPDATE_OUTPUTINC)
                strncat(varname_an, "_inc", NC_MAX_NAME - 1);
            else
                strncat(varname_an, "_an", NC_MAX_NAME - 1);

            enkf_printf("    %s:", varname_an);
            for (k = 0; k < nlev; ++k) {
                char fname_src[MAXSTRLEN];
                int ncid_src, vid;

                if (nlev > 1)
                    getfieldfname(DIRNAME_TMP, "spread", varname_an, k, fname_src);
                else
                    getfieldfname(DIRNAME_TMP, "spread", varname_an, grid_getsurflayerid(model_getvargrid(m, i)), fname_src);
                ncw_open(fname_src, NC_NOWRITE, &ncid_src);
                ncw_inq_varid(ncid_src, varname_an, &vid);
                ncw_get_var_float(ncid_src, vid, v);
                ncw_close(ncid_src);
                file_delete(fname_src);

                model_writefieldas(m, FNAME_SPREAD, varname_an, varname, k, v, 1);
                enkf_printf(".");
                enkf_flush();
            }
        }
        enkf_printf("\n");
        enkf_flush();
        free(v);
    }
}

/** Allocates disk space for inflation magnitudes.
 */
void das_allocateinflation(dasystem* das, char fname[])
{
    model* m = das->m;
    int nvar = model_getnvar(m);
    char fname_src[MAXSTRLEN];
    int ncid, ncid_src;
    int vid;

    if (rank != 0)
        return;

    if (file_exists(fname))
        return;

    ncw_create(fname, NC_CLOBBER | das->ncformat, &ncid);
    for (vid = 0; vid < nvar; ++vid) {
        char* varname_src = model_getvarname(m, vid);
        int varid_src;

        das_getmemberfname(das, varname_src, 1, fname_src);
        ncw_open(fname_src, NC_NOWRITE, &ncid_src);
        ncw_inq_varid(ncid_src, varname_src, &varid_src);
        ncw_copy_vardef(ncid_src, varid_src, ncid);
        ncw_close(ncid_src);
    }
#if defined(DEFLATE_ALL)
    if (das->nccompression > 0)
        ncw_def_deflate(ncid, 0, 1, das->nccompression);
#endif
    ncw_close(ncid);
}

/**
 */
void das_writeinflation(dasystem* das, field* f, int j, float* v)
{
    assert(das->mode == MODE_ENKF || das->mode == MODE_HYBRID);

    if (das->updatespec & UPDATE_DIRECTWRITE)
        ncu_writerow(FNAME_INFLATION, f->varname, f->level, j, v);
    else {
        char fname[MAXSTRLEN];
        int ncid;

        getfieldfname(DIRNAME_TMP, "inflation", f->varname, f->level, fname);

        if (j == 0) {
            int ni, nj;
            int dimids[2];

            model_getvargridsize(das->m, f->varid, &ni, &nj, NULL);

            ncw_create(fname, NC_CLOBBER | das->ncformat, &ncid);
            ncw_def_dim(ncid, "nj", nj, &dimids[0]);
            ncw_def_dim(ncid, "ni", ni, &dimids[1]);
            ncw_def_var(ncid, f->varname, NC_FLOAT, 2, dimids, NULL);
#if defined(DEFLATE_ALL)
            if (das->nccompression > 0)
                ncw_def_deflate(ncid, 0, 1, das->nccompression);
#endif
            ncw_enddef(ncid);
        } else
            ncw_open(fname, NC_WRITE, &ncid);

        ncu_writerow(fname, f->varname, 0, j, v);
        ncw_close(ncid);
    }
}

/**
 */
void das_assembleinflation(dasystem* das)
{
    model* m = das->m;
    int nvar = model_getnvar(m);
    int i;

    assert(das->mode == MODE_ENKF || das->mode == MODE_HYBRID);

    if (rank > 0)
        return;

    for (i = 0; i < nvar; ++i) {
        char* varname = model_getvarname(m, i);
        int nlev, k;
        int ni, nj;
        float* v = NULL;

        enkf_printf("    %s:", varname);
        nlev = ncu_getnlevels(FNAME_INFLATION, varname);

        model_getvargridsize(m, i, &ni, &nj, NULL);
        v = malloc(ni * nj * sizeof(float));

        for (k = 0; k < nlev; ++k) {
            char fname_src[MAXSTRLEN];
            int ncid_src, vid;

            if (nlev > 1)
                getfieldfname(DIRNAME_TMP, "inflation", varname, k, fname_src);
            else
                getfieldfname(DIRNAME_TMP, "inflation", varname, grid_getsurflayerid(model_getvargrid(m, i)), fname_src);
            ncw_open(fname_src, NC_NOWRITE, &ncid_src);

            ncw_inq_varid(ncid_src, varname, &vid);
            ncw_get_var_float(ncid_src, vid, v);
            ncw_close(ncid_src);
            model_writefield(m, FNAME_INFLATION, varname, k, v, 1);
            file_delete(fname_src);

            enkf_printf(".");
        }
        free(v);
        enkf_printf("\n");
    }
}

/** Calculates and writes to disk 3D field of correlation coefficients between
 * surface and other layers of 3D variables.
 */
void das_writevcorrs(dasystem* das)
{
    model* m = das->m;
    int nvar = model_getnvar(m);
    int nfields = 0;
    field* fields = NULL;
    char* varname0 = NULL;
    float** v0 = NULL;
    double* var0 = NULL;
    float** v = NULL;
    double* cor = NULL;
    int ni0 = -1, nj0 = -1, nij0 = -1;
    int fid;

    enkf_printtime("  ");
    enkf_printf("  writing vertical correlations:\n");

    /*
     * allocate disk space
     */
    if (rank == 0 && !file_exists(FNAME_VERTCORR)) {
        int ncid_dst;
        int mvid;

        ncw_create(FNAME_VERTCORR, NC_CLOBBER | das->ncformat, &ncid_dst);
        for (mvid = 0; mvid < nvar; ++mvid) {
            char* varname = model_getvarname(m, mvid);
            char fname_src[MAXSTRLEN];
            int ncid_src, varid_src;

            das_getmemberfname(das, varname, 1, fname_src);
            if (ncu_getnD(fname_src, varname) != 3)
                continue;

            ncw_open(fname_src, NC_NOWRITE, &ncid_src);
            ncw_inq_varid(ncid_src, varname, &varid_src);
            ncw_copy_vardef(ncid_src, varid_src, ncid_dst);
            ncw_close(ncid_src);
        }
#if defined(DEFLATE_ALL)
        if (das->nccompression > 0)
            ncw_def_deflate(ncid_dst, 0, 1, das->nccompression);
#endif
        ncw_close(ncid_dst);
    }

    if (rank == 0)
        dir_createifabsent(DIRNAME_TMP);
#if defined(MPI)
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    das_getfields(das, -1 /* for all grids */ , &nfields, &fields);
    distribute_iterations(0, nfields - 1, nprocesses, rank, "    ");
    enkf_printf("    calculating:");
    enkf_flush();

    for (fid = my_first_iteration; fid <= my_last_iteration; ++fid) {
        field* f = &fields[fid];
        char* varname = f->varname;
        int dimids[2];
        int ksurf;
        int ni, nj, nij = nij0;
        int e, k, i;
        char fname_src[MAXSTRLEN];
        char fname_tile[MAXSTRLEN];
        int ncid_tile, vid_tile;

        {
            char fname[MAXSTRLEN];

            das_getmemberfname(das, varname, 1, fname);
            if (ncu_getnD(fname, varname) != 3)
                continue;
        }

        /*
         * (a new variable)
         */
        if (varname0 == NULL || strcmp(varname0, varname) != 0) {
            /*
             * calculate anomalies and scaled variance at surface
             */
            model_getvargridsize(m, f->varid, &ni, &nj, NULL);
            nij = ni * nj;
            if (ni != ni0 || nj != nj0) {
                if (v != NULL) {
                    free(v0);
                    free(var0);
                    free(v);
                    free(cor);
                }
                /*
                 * for mode = MODE_HYBRID allocate two additional members to
                 * calculate ensemble mean with double precision
                 */
                v0 = alloc2d((das->mode == MODE_HYBRID) ? das->nmem + 2 : das->nmem, nij, sizeof(float));
                var0 = calloc(nij, sizeof(double));
                v = alloc2d((das->mode == MODE_HYBRID) ? das->nmem + 2 : das->nmem, nij, sizeof(float));
                cor = calloc(nij, sizeof(double));
            }

            ksurf = grid_getsurflayerid(model_getvargrid(m, f->varid));
            for (e = 0; e < das->nmem; ++e) {
                int masklog = das_isstatic(das, e + 1);
                char fname[MAXSTRLEN];

                das_getmemberfname(das, varname, e + 1, fname);
                model_readfield(das->m, fname, varname, ksurf, v0[e], masklog);
            }
            if (das->mode == MODE_HYBRID)
                das_sethybridensemble(das, nij, v0);
            for (i = 0; i < nij; ++i) {
                double vmean = 0.0;

                for (e = 0; e < das->nmem; ++e)
                    vmean += (double) v0[e][i];
                vmean /= (double) das->nmem;
                for (e = 0; e < das->nmem; ++e)
                    v0[e][i] -= (float) vmean;
                var0[i] = 0.0;
                for (e = 0; e < das->nmem; ++e)
                    var0[i] += (double) (v0[e][i] * v0[e][i]);
            }
            varname0 = varname;
            ni0 = ni;
            nj0 = nj;
            nij0 = ni * nj;
        }

        /*
         * calculate correlation coefficients at each layer and save to disk
         */
        k = f->level;
        for (e = 0; e < das->nmem; ++e) {
            int masklog = das_isstatic(das, e + 1);

            das_getmemberfname(das, varname, e + 1, fname_src);
            model_readfield(das->m, fname_src, varname, k, v[e], masklog);
        }
        if (das->mode == MODE_HYBRID)
            das_sethybridensemble(das, nij, v);
        for (i = 0; i < nij; ++i) {
            double vmean = 0.0;
            double var = 0.0;

            for (e = 0; e < das->nmem; ++e)
                vmean += (double) v[e][i];
            vmean /= (double) das->nmem;
            for (e = 0; e < das->nmem; ++e)
                v[e][i] -= (float) vmean;
            for (e = 0; e < das->nmem; ++e)
                var += (double) (v[e][i] * v[e][i]);
            cor[i] = 0.0;
            for (e = 0; e < das->nmem; ++e)
                cor[i] += (double) (v[e][i] * v0[e][i]);
            cor[i] /= sqrt(var * var0[i]);
            if (!isnormal(cor[i]))
                cor[i] = 0.0;
        }

        snprintf(fname_tile, MAXSTRLEN, "%s/vcorr_%s-%03d.nc", DIRNAME_TMP, varname, k);
        ncw_create(fname_tile, NC_CLOBBER | das->ncformat, &ncid_tile);
        ncw_def_dim(ncid_tile, "nj", nj, &dimids[0]);
        ncw_def_dim(ncid_tile, "ni", ni, &dimids[1]);
        ncw_def_var(ncid_tile, varname, NC_FLOAT, 2, dimids, &vid_tile);
#if defined(DEFLATE_ALL)
        if (das->nccompression > 0)
            ncw_def_deflate(ncid_tile, 0, 1, das->nccompression);
#endif
        ncw_enddef(ncid_tile);
        ncw_put_var_double(ncid_tile, vid_tile, cor);
        ncw_close(ncid_tile);
        enkf_printf(".");
        enkf_flush();
    }                           /* for fid */
    if (v0 != NULL) {
        free(v0);
        free(var0);
        free(v);
        free(cor);
    }
#if defined(MPI)
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    enkf_printf("\n");

    /*
     * assemble tiles
     */
    enkf_printtime("    ");
    enkf_printf("    assembling:\n");
    if (rank == 0) {
        for (fid = 0; fid < nfields; ++fid) {
            field* f = &fields[fid];
            int ni, nj;
            float* vv = NULL;
            char fname_tile[MAXSTRLEN];
            int ncid_tile, vid_tile;

            {
                char fname[MAXSTRLEN];

                das_getmemberfname(das, f->varname, 1, fname);
                if (ncu_getnD(fname, f->varname) != 3)
                    continue;
            }
            model_getvargridsize(m, f->varid, &ni, &nj, NULL);
            vv = malloc(ni * nj * sizeof(float));

            snprintf(fname_tile, MAXSTRLEN, "%s/vcorr_%s-%03d.nc", DIRNAME_TMP, f->varname, f->level);

            ncw_open(fname_tile, NC_NOWRITE, &ncid_tile);
            ncw_inq_varid(ncid_tile, f->varname, &vid_tile);
            ncw_get_var_float(ncid_tile, vid_tile, vv);
            ncw_close(ncid_tile);
            file_delete(fname_tile);

            model_writefield(m, FNAME_VERTCORR, f->varname, f->level, vv, 1);
            free(vv);
        }
    }
    free(fields);
#if defined(MPI)
    MPI_Barrier(MPI_COMM_WORLD);
#endif
}
