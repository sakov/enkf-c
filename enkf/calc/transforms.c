/******************************************************************************
 *
 * File:        transforms.c        
 *
 * Created:     11/2013
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description: Calculating and saving ensemble transforms.
 *
 * Revisions:
 *
 *****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include "definitions.h"
#include "distribute.h"
#include "utils.h"
#include "ncutils.h"
#if defined(USE_MPIQUEUE)
#include "mpiqueue.h"
#endif
#include "calcs.h"
#include "dasystem.h"
#include "pointlog.h"

#define DFS_MIN 1.0e-6

#if defined(MINIMISE_ALLOC)
#define PLOC_INC 1000
int ploc_allocated1 = 0;
void* storage = NULL;

int ploc_allocated2 = 0;
int* lobs = NULL;
double* lcoeffs = NULL;
#endif

typedef struct {
    long long int nlobs_sum;
    int nlobs_max;
    int n_inv_obs;
    int n_inv_ens;
    int ncell;
} calcstats;

/**
 */
static void nc_createtransforms(dasystem* das, int gridid, size_t nj, size_t ni, int stride)
{
    char fname[MAXSTRLEN];
    int ncid;
    int dimids[4];
    int varid_T, varid_w;

    assert(rank == 0);

    das_getfname_transforms(das, gridid, fname);
    enkf_printf("      creating empty file \"%s\":\n", fname);
    enkf_flush();
    ncw_create(fname, NC_CLOBBER | NC_NOFILL | das->ncformat, &ncid);
    ncw_def_dim(ncid, "nj", nj, &dimids[0]);
    ncw_def_dim(ncid, "ni", ni, &dimids[1]);
    if (das->mode == MODE_ENKF || das->mode == MODE_HYBRID) {
        ncw_def_dim(ncid, "m_dyn", das->nmem_dynamic, &dimids[2]);
        ncw_def_dim(ncid, "m", das->nmem, &dimids[3]);
        ncw_def_var(ncid, "T", NC_FLOAT, 4, dimids, &varid_T);
        {
            size_t chunksize[4] = { 1, ni, das->nmem_dynamic, das->nmem };

            nc_def_var_chunking(ncid, varid_T, NC_CHUNKED, chunksize);
        }
        dimids[2] = dimids[3];
    } else
        ncw_def_dim(ncid, "m", das->nmem, &dimids[2]);
    ncw_def_var(ncid, "w", NC_FLOAT, 3, dimids, &varid_w);
    {
        size_t chunksize[3] = { 1, ni, das->nmem };

        nc_def_var_chunking(ncid, varid_w, NC_CHUNKED, chunksize);
    }
    ncw_put_att_int(ncid, NC_GLOBAL, "stride", 1, &stride);
    ncw_put_att_text(ncid, NC_GLOBAL, "grid_name", grid_getname(model_getgridbyid(das->m, gridid)));
#if defined(DEFLATE_ALL)
    if (das->nccompression > 0)
        ncw_def_deflate(ncid, 0, 1, das->nccompression);
#endif
    ncw_close(ncid);
}

#if !defined(TW_VIAFILE)
/**
 */
static void nc_writetransforms(dasystem* das, int gridid, int j, int ni, float* Tj, float* wj)
{
    char fname[MAXSTRLEN];
    int ncid;
    int varid_T, varid_w;
    size_t start[4], count[4];

    assert(rank == 0);

    das_getfname_transforms(das, gridid, fname);
    ncw_open(fname, NC_WRITE, &ncid);

    start[0] = j;
    start[1] = 0;
    start[2] = 0;
    start[3] = 0;

    count[0] = 1;
    count[1] = ni;
    count[2] = das->nmem_dynamic;
    count[3] = das->nmem;

    if (Tj != NULL) {
        ncw_inq_varid(ncid, "T", &varid_T);
        ncw_put_vara_float(ncid, varid_T, start, count, Tj);
    }

    ncw_inq_varid(ncid, "w", &varid_w);
    count[2] = das->nmem;
    ncw_put_vara_float(ncid, varid_w, start, count, wj);

    ncw_close(ncid);
}
#endif

#if defined(TW_VIAFILE)
/**
 */
static void das_getfname_transformstile(dasystem* das, int gridid, int r, char fname[])
{
    snprintf(fname, MAXSTRLEN, "%s/%s-%d-%03d.nc", DIRNAME_TMP, FNAMEPREFIX_TRANSFORMS, gridid, r);
}

/**
 */
static void nc_createtransformstile(dasystem* das, int gridid, int ni)
{
    char fname[MAXSTRLEN];
    int ncid;
    int varid_T, varid_w;
    int dimids[4];

    das_getfname_transformstile(das, gridid, rank, fname);

    ncw_create(fname, NC_CLOBBER | NC_NOFILL | das->ncformat, &ncid);
    ncw_def_dim(ncid, "nj", NC_UNLIMITED, &dimids[0]);
    ncw_def_dim(ncid, "ni", ni, &dimids[1]);
    ncw_def_var(ncid, "j", NC_INT, 1, &dimids[0], NULL);
    if (das->mode == MODE_ENKF || das->mode == MODE_HYBRID) {
        ncw_def_dim(ncid, "m_dyn", das->nmem_dynamic, &dimids[2]);
        ncw_def_dim(ncid, "m", das->nmem, &dimids[3]);
        ncw_def_var(ncid, "T", NC_FLOAT, 4, dimids, &varid_T);
        dimids[2] = dimids[3];
    } else
        ncw_def_dim(ncid, "m", das->nmem, &dimids[2]);
    ncw_def_var(ncid, "w", NC_FLOAT, 3, dimids, &varid_w);
#if defined(DEFLATE_ALL)
    if (das->nccompression > 0)
        ncw_def_deflate(ncid, 0, 1, das->nccompression);
#endif
    ncw_close(ncid);
}

/**
 */
static void nc_writetransformstile(dasystem* das, int gridid, int iter, int ni, float* Tj, float* wj)
{
    char fname[MAXSTRLEN];
    int ncid;
    int varid_T, varid_w, varid_j;
    size_t start[4], count[4];
    int nr;

    das_getfname_transformstile(das, gridid, rank, fname);
    ncw_open(fname, NC_WRITE, &ncid);

    nr = ncw_inq_nrecords(ncid);
    start[0] = nr;
    start[1] = 0;
    start[2] = 0;
    start[3] = 0;

    count[0] = 1;
    count[1] = ni;
    count[2] = das->nmem_dynamic;
    count[3] = das->nmem;

    if (Tj != NULL) {
        ncw_inq_varid(ncid, "T", &varid_T);
        ncw_put_vara_float(ncid, varid_T, start, count, Tj);
    }

    ncw_inq_varid(ncid, "w", &varid_w);
    count[2] = das->nmem;
    ncw_put_vara_float(ncid, varid_w, start, count, wj);

    ncw_inq_varid(ncid, "j", &varid_j);
    ncw_put_vara_int(ncid, varid_j, start, count, &iter);
    ncw_close(ncid);
}

/**
 */
static void nc_assembletransforms(dasystem* das, int gridid, size_t nj, size_t ni, int stride)
{
    char fname[MAXSTRLEN];
    int ncid, varid_T, varid_w;
    int p;
    float* v_T = NULL;
    float* v_w = NULL;
    int* v_j = NULL;
    int doT = (das->mode == MODE_ENKF || das->mode == MODE_HYBRID);
    int nr_max = 0;

    assert(rank == 0);

    das_getfname_transforms(das, gridid, fname);
    nc_createtransforms(das, gridid, nj, ni, stride);
    ncw_open(fname, NC_WRITE, &ncid);
    if (doT)
        ncw_inq_varid(ncid, "T", &varid_T);
    ncw_inq_varid(ncid, "w", &varid_w);

    enkf_printf("      assembling \"%s\":", fname);

    for (p = 0; p < nprocesses; ++p) {
        char fname_tile[MAXSTRLEN];
        int ncid_tile;
        int varid_T_tile, varid_w_tile, varid_j_tile;
        size_t start[4], count[4];
        int nr, r;

        das_getfname_transformstile(das, gridid, p, fname_tile);
        if (!file_exists(fname_tile))
            continue;

        ncw_open(fname_tile, NC_NOWRITE, &ncid_tile);
        nr = ncw_inq_nrecords(ncid_tile);
        if (nr == 0) {
            ncw_close(ncid_tile);
            file_delete(fname_tile);
            continue;
        }
        if (nr > nr_max) {
            if (doT)
                v_T = realloc(v_T, nr * ni * das->nmem_dynamic * das->nmem * sizeof(float));
            v_w = realloc(v_w, nr * ni * das->nmem * sizeof(float));
            v_j = realloc(v_j, nr * sizeof(int));
            nr_max = nr;
        }
        if (doT) {
            ncw_inq_varid(ncid_tile, "T", &varid_T_tile);
            ncu_readvarfloat(ncid_tile, varid_T_tile, nr * ni * das->nmem * das->nmem_dynamic, v_T);
        }
        ncw_inq_varid(ncid_tile, "w", &varid_w_tile);
        ncu_readvarfloat(ncid_tile, varid_w_tile, nr * ni * das->nmem, v_w);
        ncw_inq_varid(ncid_tile, "j", &varid_j_tile);
        ncw_get_var_int(ncid_tile, varid_j_tile, v_j);
        ncw_close(ncid_tile);
        file_delete(fname_tile);

        for (r = 0; r < nr; ++r) {
            start[0] = v_j[r];
            start[1] = 0;
            start[2] = 0;
            start[3] = 0;

            count[0] = 1;
            count[1] = ni;
            count[2] = das->nmem_dynamic;
            count[3] = das->nmem;

            if (doT)
                ncw_put_vara_float(ncid, varid_T, start, count, &v_T[r * ni * das->nmem * das->nmem_dynamic]);
            count[2] = das->nmem;
            ncw_put_vara_float(ncid, varid_w, start, count, &v_w[r * ni * das->nmem]);
        }
        enkf_printf(".");
        enkf_flush();
    }
    ncw_close(ncid);
    enkf_printf("\n");
    enkf_flush();

    if (doT)
        free(v_T);
    free(v_w);
    free(v_j);
}
#endif                          /* TW_VIAFILE */

#if defined(USE_MPIQUEUE)
/**
 */
static void das_getfname_diagtile(dasystem* das, int gridid, int r, char fname[])
{
    snprintf(fname, MAXSTRLEN, "%s/%s-%d-%03d.nc", DIRNAME_TMP, FNAMEPREFIX_DIAG, gridid, r);
}

/**
 */
static void nc_creatediagtile(dasystem* das, int gridid, int ni)
{
    char fname[MAXSTRLEN];
    int ncid;
    int dimids[2];

    das_getfname_diagtile(das, gridid, rank, fname);
    ncw_create(fname, NC_CLOBBER | das->ncformat, &ncid);
    ncw_def_dim(ncid, "nj", NC_UNLIMITED, &dimids[0]);
    ncw_def_dim(ncid, "ndata", ni * (das->obs->nobstypes + 1), &dimids[1]);
    ncw_def_var(ncid, "j", NC_INT, 1, &dimids[0], NULL);
    ncw_def_var(ncid, "pnlobs", NC_INT, 2, &dimids[0], NULL);
    ncw_def_var(ncid, "pdfs", NC_FLOAT, 2, &dimids[0], NULL);
    ncw_def_var(ncid, "psrf", NC_FLOAT, 2, &dimids[0], NULL);
    ncw_close(ncid);
}

/**
 */
static void nc_writediagtile(dasystem* das, int gridid, int ni, int j, int*** pnlobs, float*** pdfs, float*** psrf)
{
    char fname[MAXSTRLEN];
    int ncid;
    size_t start[2], count[2];
    int varid, nr;

    das_getfname_diagtile(das, gridid, rank, fname);
    ncw_open(fname, NC_WRITE | das->ncformat, &ncid);

    nr = ncw_inq_nrecords(ncid);
    start[0] = nr;
    start[1] = 0;

    count[0] = 1;
    count[1] = ni * (das->obs->nobstypes + 1);

    ncw_inq_varid(ncid, "j", &varid);
    ncw_put_vara_int(ncid, varid, start, count, &j);

    ncw_inq_varid(ncid, "pnlobs", &varid);
    ncw_put_vara_int(ncid, varid, start, count, pnlobs[0][0]);

    ncw_inq_varid(ncid, "pdfs", &varid);
    ncw_put_vara_float(ncid, varid, start, count, pdfs[0][0]);

    ncw_inq_varid(ncid, "psrf", &varid);
    ncw_put_vara_float(ncid, varid, start, count, psrf[0][0]);

    ncw_close(ncid);
}
#endif

/**
 */
static void das_getfname_diag(dasystem* das, int gridid, char fname[])
{
    if (model_getngrid(das->m) == 1)
        snprintf(fname, MAXSTRLEN, "%s.nc", FNAMEPREFIX_DIAG);
    else
        snprintf(fname, MAXSTRLEN, "%s-%d.nc", FNAMEPREFIX_DIAG, gridid);
}

#if defined(USE_MPIQUEUE)
/**
 */
void nc_assemblediag(dasystem* das, int gridid, int nj, int ni, int stride)
{
    int nobstypes = das->obs->nobstypes;
    char fname[MAXSTRLEN];
    int ncid;
    int dimids[3];
    int varid_nlobs, varid_dfs, varid_srf, varid_pnlobs, varid_pdfs, varid_psrf;
    size_t start[3], count[3];
    int p, otid;
    int** pnlobs;
    float** pdfs;
    float** psrf;

    assert(rank == 0);

    das_getfname_diag(das, gridid, fname);
    ncw_create(fname, NC_CLOBBER | das->ncformat, &ncid);

    ncw_def_dim(ncid, "nobstypes", nobstypes, &dimids[0]);
    ncw_def_dim(ncid, "nj", nj, &dimids[1]);
    ncw_def_dim(ncid, "ni", ni, &dimids[2]);
    ncw_def_var(ncid, "nlobs", NC_INT, 2, &dimids[1], &varid_nlobs);
    ncw_def_var(ncid, "dfs", NC_FLOAT, 2, &dimids[1], &varid_dfs);
    ncw_def_var(ncid, "srf", NC_FLOAT, 2, &dimids[1], &varid_srf);
    ncw_def_var(ncid, "pnlobs", NC_INT, 3, dimids, &varid_pnlobs);
    ncw_def_var(ncid, "pdfs", NC_FLOAT, 3, dimids, &varid_pdfs);
    ncw_def_var(ncid, "psrf", NC_FLOAT, 3, dimids, &varid_psrf);
    ncw_put_att_int(ncid, NC_GLOBAL, "stride", 1, &stride);
    for (otid = 0; otid < das->obs->nobstypes; ++otid)
        ncw_put_att_int(ncid, NC_GLOBAL, das->obs->obstypes[otid].name, 1, &otid);
    ncw_enddef(ncid);

    pnlobs = alloc2d(nobstypes + 1, ni, sizeof(int));
    pdfs = alloc2d(nobstypes + 1, ni, sizeof(float));
    psrf = alloc2d(nobstypes + 1, ni, sizeof(float));

    enkf_printf("      assembling \"%s\":", fname);

    for (p = 1; p < nprocesses; ++p) {
        char fname_tile[MAXSTRLEN];
        int ncid_tile;
        int nr, r, j;
        int varid_j_tile, varid_pnlobs_tile, varid_pdfs_tile, varid_psrf_tile;

        das_getfname_diagtile(das, gridid, p, fname_tile);
        if (!file_exists(fname_tile))
            continue;

        ncw_open(fname_tile, NC_NOWRITE, &ncid_tile);
        ncw_inq_varid(ncid_tile, "j", &varid_j_tile);
        ncw_inq_varid(ncid_tile, "pnlobs", &varid_pnlobs_tile);
        ncw_inq_varid(ncid_tile, "pdfs", &varid_pdfs_tile);
        ncw_inq_varid(ncid_tile, "psrf", &varid_psrf_tile);
        nr = ncw_inq_nrecords(ncid_tile);
        if (nr == 0) {
            ncw_close(ncid_tile);
            file_delete(fname_tile);
            continue;
        }
        for (r = 0; r < nr; ++r) {
            start[0] = r;
            start[1] = 0;
            count[0] = 1;
            count[1] = (nobstypes + 1) * ni;

            ncw_get_vara_int(ncid_tile, varid_j_tile, start, count, &j);
            ncw_get_vara_int(ncid_tile, varid_pnlobs_tile, start, count, pnlobs[0]);
            ncw_get_vara_float(ncid_tile, varid_pdfs_tile, start, count, pdfs[0]);
            ncw_get_vara_float(ncid_tile, varid_psrf_tile, start, count, psrf[0]);

            start[0] = 0;
            start[1] = j;
            start[2] = 0;
            count[0] = nobstypes;
            count[1] = 1;
            count[2] = ni;
            ncw_put_vara_int(ncid, varid_nlobs, &start[1], &count[1], pnlobs[nobstypes]);
            ncw_put_vara_float(ncid, varid_dfs, &start[1], &count[1], pdfs[nobstypes]);
            ncw_put_vara_float(ncid, varid_srf, &start[1], &count[1], psrf[nobstypes]);

            ncw_put_vara_int(ncid, varid_pnlobs, start, count, pnlobs[0]);
            ncw_put_vara_float(ncid, varid_pdfs, start, count, pdfs[0]);
            ncw_put_vara_float(ncid, varid_psrf, start, count, psrf[0]);
        }
        ncw_close(ncid_tile);
        file_delete(fname_tile);

        enkf_printf(".");
        enkf_flush();
    }
    ncw_close(ncid);
    enkf_printf("\n");
    enkf_flush();

    free(pnlobs);
    free(pdfs);
    free(psrf);
}
#else                           /* ! defined(USE_MPIQUEUE) */
/**
 */
static void nc_writediag(dasystem* das, char fname[], int nobstypes, int nj, int ni, int stride, int** nlobs, float** dfs, float** srf, int*** pnlobs, float*** pdfs, float*** psrf)
{
    int ncid;
    int dimids[3];
    int varid_nlobs, varid_dfs, varid_srf, varid_pnlobs, varid_pdfs, varid_psrf;
    int i;

    assert(rank == 0);

    enkf_printf("    writing stats to \"%s\":\n", fname);
    ncw_create(fname, NC_CLOBBER | das->ncformat, &ncid);
    ncw_def_dim(ncid, "nobstypes", nobstypes, &dimids[0]);
    ncw_def_dim(ncid, "nj", nj, &dimids[1]);
    ncw_def_dim(ncid, "ni", ni, &dimids[2]);
    ncw_def_var(ncid, "nlobs", NC_INT, 2, &dimids[1], &varid_nlobs);
    ncw_def_var(ncid, "dfs", NC_FLOAT, 2, &dimids[1], &varid_dfs);
    ncw_def_var(ncid, "srf", NC_FLOAT, 2, &dimids[1], &varid_srf);
    ncw_def_var(ncid, "pnlobs", NC_INT, 3, dimids, &varid_pnlobs);
    ncw_def_var(ncid, "pdfs", NC_FLOAT, 3, dimids, &varid_pdfs);
    ncw_def_var(ncid, "psrf", NC_FLOAT, 3, dimids, &varid_psrf);
    ncw_put_att_int(ncid, NC_GLOBAL, "stride", 1, &stride);
    for (i = 0; i < das->obs->nobstypes; ++i)
        ncw_put_att_int(ncid, NC_GLOBAL, das->obs->obstypes[i].name, 1, &i);
    if (das->nccompression > 0)
        ncw_def_deflate(ncid, 0, 1, das->nccompression);
    ncw_enddef(ncid);

    ncw_put_var_int(ncid, varid_nlobs, nlobs[0]);
    ncw_put_var_float(ncid, varid_dfs, dfs[0]);
    ncw_put_var_float(ncid, varid_srf, srf[0]);
    ncw_put_var_int(ncid, varid_pnlobs, pnlobs[0][0]);
    ncw_put_var_float(ncid, varid_pdfs, pdfs[0][0]);
    ncw_put_var_float(ncid, varid_psrf, psrf[0][0]);

    ncw_close(ncid);
}
#endif

#if !defined(TW_VIAFILE)
#if !defined(SHUFFLE_ROWS)
/** Distributes iterations in such a way that for any given i iterations # i 
 * for each process are grouped together.
 */
static void group_iterations(int n, int ids[])
{
    int iter, r, i;

    assert(first_iteration[0] == 0);
    assert(last_iteration[nprocesses - 1] == n - 1);
    for (iter = 0, i = 0; iter <= last_iteration[0]; ++iter)
        for (r = 0; r < nprocesses && i < n; ++r, ++i)
            ids[first_iteration[r] + iter] = i;
}
#endif
#endif

#if defined(MINIMISE_ALLOC)
/**
 */
static void prepare_calcs(size_t ploc, size_t nmem, double** sloc, int** plobs, double*** Sloc, double*** G, double*** M)
{
    if (ploc > ploc_allocated1) {
        size_t size;

        ploc_allocated1 = ploc + PLOC_INC;
        /*
         * sloc [ploc]
         */
        size = ploc_allocated1 * sizeof(double);
        /*
         * plobs [ploc]
         */
        size += ploc_allocated1 * sizeof(int);
        /*
         * Sloc [nmem][ploc]
         */
        size += ploc_allocated1 * nmem * sizeof(double) + nmem * sizeof(void*);
        /*
         * G [ploc][nmem]
         */
        size += ploc_allocated1 * nmem * sizeof(double) + ploc_allocated1 * sizeof(void*);
        /*
         * M [nmem * 2 + 11][nmem] (the surplus is used for work arrays and
         * matrices in calc_wT_etkf()
         */
        size += nmem * (2 * nmem + 11) * sizeof(double) + (2 * nmem + 11) * sizeof(void*);

        storage = realloc(storage, size);
        if (storage == NULL) {
            int errno_saved = errno;

            enkf_quit("prepare_calcs(): realloc(): %s", strerror(errno_saved));
        }
    }

    *sloc = storage;
    *plobs = (void*) &(*sloc)[ploc];
    *Sloc = cast2d(&(*plobs)[ploc], nmem, ploc, sizeof(double));
    *G = cast2d(&(*Sloc)[0][nmem * ploc], ploc, nmem, sizeof(double));
    *M = (void*) &(*G)[0][nmem * ploc];
}
#endif

#if defined(USE_MPIQUEUE)
/*
 * Generally, this is a preferred way of calculating the updates because it
 * reduces idle periods occuring when some of the rows are processed much faster
 * than others.
 */

/** The central DA procedure, where the actual transforms are calculated.
 */
void das_calctransforms(dasystem* das)
{
    model* m = das->m;
    int nmem = das->nmem;
    int nmem_dynamic = das->nmem_dynamic;
    double* w = malloc(nmem * sizeof(double));
    observations* obs = das->obs;
    int ngrid = model_getngrid(m);
    int gid;

    assert(das->s_mode == S_MODE_HA_f);
    das_standardise(das);

    for (gid = 0; gid < ngrid; ++gid) {
        void* g = model_getgridbyid(m, gid);
        char* gridname = grid_getname(g);
        int stride = grid_getstride(g);

        int mni, mnj;
        int nj, ni;
        int* jiter = NULL;
        int* iiter = NULL;

        mpiqueue* queue = NULL;

        /*
         * transforms 
         */
        float*** Tj = NULL;     /* T for one grid row [ni][m x m] */
        double** T = NULL;      /* T for one grid cell [m][m] */

        /*
         * coeffs 
         */
        float** wj = NULL;

        /*
         * stats 
         */
        calcstats stats = { 0, 0, 0, 0, 0 };
        int*** pnlobs = NULL;
        float*** pdfs = NULL;
        float*** psrf = NULL;

        int i, j, ii, jj, ot;

        if (grid_getaliasid(g) >= 0)
            continue;

        /*
         * skip this grid if there are no model variables associated with it
         */
        for (i = 0; i < model_getnvar(m); ++i)
            if (model_getvargridid(m, i) == gid)
                break;
        if (i == model_getnvar(m))
            continue;

        enkf_printf("    calculating transforms for %s:\n", gridname);

        grid_getsize(g, &mni, &mnj, NULL);
        /*
         * a treatment for unstructured grids
         */
        if (mnj <= 0) {
            mnj = mni;
            mni = 1;
        }

        /*
         * work out how to cycle j 
         */
        nj = (mnj - 1) / stride + 1;
        jiter = malloc(nj * sizeof(int));
        for (j = 0, i = 0; j < nj; ++j, i += stride)
            jiter[j] = i;
        ni = (mni - 1) / stride + 1;
        iiter = malloc(ni * sizeof(int));
        for (i = 0, j = 0; i < ni; ++i, j += stride)
            iiter[i] = j;

        enkf_printf("      main cycle for %s (%d x %d local analyses):\n", gridname, nj, ni);
        enkf_flush();
        if (rank == 0)
            dir_createifabsent(DIRNAME_TMP);
        MPI_Barrier(MPI_COMM_WORLD);

        if (nprocesses == 1)
            enkf_quit("\"mpiqueue\" can not be used on a single CPU; run on more than one CPU or recompile without -DUSE_MPIQUEUE flag");
        queue = mpiqueue_create(MPI_COMM_WORLD, nj);

        if (rank == 0)
            mpiqueue_manage(queue);
        else {
            if (das->mode == MODE_ENKF || das->mode == MODE_HYBRID) {
                Tj = alloc3d(ni, nmem_dynamic, nmem, sizeof(float));
                T = alloc2d(nmem, nmem, sizeof(double));
            }
            wj = alloc2d(ni, nmem, sizeof(float));

            pnlobs = alloc3d(obs->nobstypes + 1, 1, ni, sizeof(int));
            pdfs = alloc3d(obs->nobstypes + 1, 1, ni, sizeof(float));
            psrf = alloc3d(obs->nobstypes + 1, 1, ni, sizeof(float));

            nc_createtransformstile(das, gid, ni);
            nc_creatediagtile(das, gid, ni);

            while ((jj = mpiqueue_getjob(queue)) >= 0) {
                /*
                 * j is "j" index of the actual (physical) grid (0 to mnj);
                 * jj is "j" index of the (strided) transforms subrid (0 to nj)
                 */
                j = jiter[jj];
                if (enkf_verbose)
                    printf("        j = %d (%d: %d: %.1f%%)\n", j, rank, jj, 100.0 * (double) jj / (double) nj);
                fflush(stdout);

                for (ii = 0; ii < ni; ++ii) {
                    int ploc = 0;       /* `nlobs' already engaged, using
                                         * another name */

                    double* sloc = NULL;
                    int* plobs = NULL;
                    double** Sloc = NULL;
                    double** G = NULL;
                    double** M = NULL;

#if !defined(MINIMISE_ALLOC)
                    int* lobs = NULL;
                    double* lcoeffs = NULL;
#endif
                    int e, o;

                    i = iiter[ii];

                    {
                        int ij[2] = { i, j };
                        double lon, lat;

                        /*
                         * for unstructured grids -- restore the natural order
                         */
                        if (mni == 1) {
                            ij[0] = j;
                            ij[1] = 0;
                        }

                        grid_ij2xy(g, ij, &lon, &lat);

#if defined(MINIMISE_ALLOC)
                        obs_findlocal(obs, lon, lat, grid_isgeographic(g), grid_getdomainname(g), &ploc, &lobs, &lcoeffs, &ploc_allocated2);
#else
                        obs_findlocal(obs, lon, lat, grid_isgeographic(g), grid_getdomainname(g), &ploc, &lobs, &lcoeffs, NULL);
#endif
                        assert(ploc >= 0 && ploc <= obs->nobs);
                    }
                    if (ploc > stats.nlobs_max)
                        stats.nlobs_max = ploc;
                    stats.nlobs_sum += ploc;
                    stats.ncell++;

                    if (ploc == 0) {
                        if (T != NULL) {
                            memset(Tj[ii][0], 0, nmem_dynamic * nmem * sizeof(float));
                            for (e = 0; e < nmem_dynamic; ++e)
                                Tj[ii][e][e] = (float) 1.0;
                        }
                        memset(wj[ii], 0, nmem * sizeof(float));

                        for (ot = 0; ot <= obs->nobstypes; ot++) {
                            pnlobs[ot][0][ii] = 0;
                            pdfs[ot][0][ii] = 0.0;
                            psrf[ot][0][ii] = 0.0;
                        }
                        continue;
                    }
#if defined(MINIMISE_ALLOC)
                    /*
                     * prepare_calcs() sets Sloc, G and M matrices while trying
                     * to minimise the number of dynamic allocations
                     */
                    prepare_calcs(ploc, nmem, &sloc, &plobs, &Sloc, &G, &M);
#else
                    Sloc = alloc2d(nmem, ploc, sizeof(double));
                    G = alloc2d(ploc, nmem, sizeof(double));
                    sloc = malloc(ploc * sizeof(double));
                    plobs = malloc(ploc * sizeof(int));
#endif

                    for (e = 0; e < nmem; ++e) {
                        float* Se = das->S[e];
                        double* Sloce = Sloc[e];

                        for (o = 0; o < ploc; ++o)
                            Sloce[o] = (double) Se[lobs[o]] * lcoeffs[o];
                    }
                    for (o = 0; o < ploc; ++o)
                        sloc[o] = das->s_f[lobs[o]] * lcoeffs[o];

                    if (das->mode == MODE_ENOI) {
                        calc_G(nmem, ploc, M, Sloc, G);
                        calc_w(nmem, ploc, G, sloc, w);
                    } else if (das->scheme == SCHEME_DENKF) {
                        calc_wT_denkf(nmem, nmem_dynamic, ploc, sloc, Sloc, Sloc, M, G, w, T);
                    } else if (das->scheme == SCHEME_ETKF) {
                        calc_wT_etkf(nmem, nmem_dynamic, ploc, sloc, Sloc, Sloc, M, G, w, T);
                    } else
                        enkf_quit("programming error");

                    if (das->mode == MODE_ENKF || das->mode == MODE_HYBRID) {
                        int e1, e2;

                        /*
                         * incorporate "relaxation to prior"
                         * T = I + alpha * (T - I)
                         */
                        for (e1 = 0; e1 < nmem_dynamic; ++e1) {
                            double* Ti = T[e1];

                            Ti[e1] -= 1.0;
                            for (e2 = 0; e2 < nmem; ++e2)
                                Ti[e2] *= das->alpha;
                            Ti[e1] += 1.0;
                        }

                        /*
                         * convert T to float and store in Tj 
                         */
                        for (e1 = 0; e1 < nmem_dynamic; ++e1)
                            for (e2 = 0; e2 < nmem; ++e2)
                                Tj[ii][e1][e2] = (float) T[e1][e2];
                    }
                    /*
                     * convert w to float and store in wj 
                     */
                    for (e = 0; e < nmem; ++e)
                        wj[ii][e] = (float) w[e];

                    pnlobs[obs->nobstypes][0][ii] = ploc;       /* (ploc > 0
                                                                 * here) */
                    pdfs[obs->nobstypes][0][ii] = traceprod(0, 0, ploc, nmem, G, Sloc, 1);
                    psrf[obs->nobstypes][0][ii] = sqrt(traceprod(0, 1, nmem, ploc, Sloc, Sloc, 1) / pdfs[obs->nobstypes][0][ii]) - 1.0;
                    for (ot = 0; ot < obs->nobstypes; ++ot) {
                        int p = 0;

                        for (o = 0; o < ploc; ++o) {
                            if (obs->data[lobs[o]].type == ot) {
                                plobs[p] = o;
                                p++;
                            }
                        }
                        pnlobs[ot][0][ii] = p;
                        if (p == 0) {
                            pdfs[ot][0][ii] = 0.0;
                            psrf[ot][0][ii] = 0.0;
                        } else {
                            /*
                             * it is used below that local observations in array
                             * plobs are continuous by observation type
                             */
#if 0
                            double** pS = alloc2d(nmem, p, sizeof(double));
                            double** pG = &G[plobs[0]];

                            for (e = 0; e < nmem; ++e)
                                for (o = 0; o < p; ++o)
                                    pS[e][o] = Sloc[e][plobs[o]];
                            pdfs[ot][0][ii] = traceprod(0, 0, p, nmem, pG, pS, 1);
                            if (pdfs[ot][0][ii] > DFS_MIN)
                                psrf[ot][0][ii] = sqrt(traceprod(0, 1, nmem, p, pS, pS, 1) / pdfs[ot][0][ii]) - 1.0;
                            else
                                psrf[ot][0][ii] = 0.0;
#else
                            double** pS = malloc(nmem * sizeof(double*));
                            double** pG = &G[plobs[0]];

                            for (e = 0; e < nmem; ++e)
                                pS[e] = &Sloc[e][plobs[0]];
                            pdfs[ot][0][ii] = traceprod(0, 0, p, nmem, pG, pS, 0);
                            if (pdfs[ot][0][ii] > DFS_MIN)
                                psrf[ot][0][ii] = sqrt(traceprod(0, 1, nmem, p, pS, pS, 0) / pdfs[ot][0][ii]) - 1.0;
                            else
                                psrf[ot][0][ii] = 0.0;
#endif
                            free(pS);
                        }
                    }

                    if (ploc >= nmem)
                        stats.n_inv_ens++;
                    else
                        stats.n_inv_obs++;

#if !defined(MINIMISE_ALLOC)
                    free(G);
                    free(Sloc);
                    free(sloc);
                    free(plobs);
                    free(lobs);
                    free(lcoeffs);
#endif
                }               /* for ii */

                if (das->mode == MODE_ENKF || das->mode == MODE_HYBRID)
                    nc_writetransformstile(das, gid, jj, ni, Tj[0][0], wj[0]);
                else if (das->mode == MODE_ENOI)
                    nc_writetransformstile(das, gid, jj, ni, NULL, wj[0]);
                nc_writediagtile(das, gid, ni, jj, pnlobs, pdfs, psrf);

                mpiqueue_reportjob(queue, jj);
            }                   /* while (jj >= 0) */
            free(pnlobs);
            free(pdfs);
            free(psrf);
        }                       /* rank > 0 */
        mpiqueue_destroy(queue);

        MPI_Barrier(MPI_COMM_WORLD);
        if (rank == 0) {
            nc_assembletransforms(das, gid, nj, ni, stride);
            nc_assemblediag(das, gid, nj, ni, stride);
            dir_rmifexists(DIRNAME_TMP);
        }
        enkf_printf("    finished calculating transforms for %s\n", gridname);
        enkf_flush();
        MPI_Barrier(MPI_COMM_WORLD);

        if (das->mode == MODE_ENKF || das->mode == MODE_HYBRID) {
            free(Tj);
            free(T);
        }
        free(wj);

#if defined(MPI)
        /*
         * merge stats for the report
         */
        if (rank > 0) {
            if (my_number_of_iterations > 0) {
                int ierror = MPI_Send(&stats, sizeof(stats) / sizeof(int), MPI_INT, 0, 99, MPI_COMM_WORLD);

                assert(ierror == MPI_SUCCESS);
            }
        } else {
            int r;

            for (r = 1; r < nprocesses; ++r) {
                calcstats morestats;
                int ierror;

                if (number_of_iterations[r] == 0)
                    continue;

                ierror = MPI_Recv(&morestats, sizeof(stats) / sizeof(int), MPI_INT, r, 99, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                assert(ierror == MPI_SUCCESS);

                stats.nlobs_sum += morestats.nlobs_sum;
                if (morestats.nlobs_max > stats.nlobs_max)
                    stats.nlobs_max = morestats.nlobs_max;
                stats.n_inv_obs += morestats.n_inv_obs;
                stats.n_inv_ens += morestats.n_inv_ens;
                stats.ncell += morestats.ncell;
            }
        }
#endif

        enkf_printf("    summary stats on %s:\n", gridname);

        enkf_printf("      # of local analyses = %d\n", stats.ncell);
        enkf_printf("      average # of local obs = %.1f\n", (double) stats.nlobs_sum / (double) stats.ncell);
        enkf_printf("      # of inversions in obs space = %d\n", stats.n_inv_obs);
        enkf_printf("      # of inversions in ens space = %d\n", stats.n_inv_ens);
        enkf_printf("");

        free(jiter);
        free(iiter);
    }                           /* for gid */

    /*
     * (this block can be put outside the grid loop because the storage sizes
     * do not depend on the grid size)
     */
#if defined(MINIMISE_ALLOC)
    if (storage != NULL) {
        free(storage);
        storage = NULL;
        ploc_allocated1 = 0;
    }
    if (lobs != NULL) {
        free(lobs);
        free(lcoeffs);
        lobs = NULL;
        lcoeffs = NULL;
        ploc_allocated2 = 0;
    }
#endif
    free(w);
}

#else                           /* !defined(USE_MPIQUEUE) */

/** The central DA procedure, where the actual transforms are calculated.
 */
void das_calctransforms(dasystem* das)
{
    model* m = das->m;
    int nmem = das->nmem;
    int nmem_dynamic = das->nmem_dynamic;
    double* w = malloc(nmem * sizeof(double));
    observations* obs = das->obs;
    int ngrid = model_getngrid(m);
    int gid;

    assert(das->s_mode == S_MODE_HA_f);
    das_standardise(das);

    for (gid = 0; gid < ngrid; ++gid) {
        void* g = model_getgridbyid(m, gid);
        char* gridname = grid_getname(g);
        int stride = grid_getstride(g);

        int mni, mnj;
        int nj, ni;
        int* jiter = NULL;
        int* iiter = NULL;

        /*
         * transforms 
         */
        float*** Tj = NULL;     /* T for one grid row [ni][m x m] */
        double** T = NULL;      /* T for one grid cell [m][m] */

        /*
         * coeffs 
         */
        float** wj = NULL;

        /*
         * stats 
         */
        calcstats stats = { 0, 0, 0, 0, 0 };
        int** nlobs = NULL;
        int*** pnlobs = NULL;
        float** dfs = NULL;
        float*** pdfs = NULL;
        float** srf = NULL;
        float*** psrf = NULL;

        int* jpool = NULL;
        int i, j, ii, jj, ot, jjj;

        if (grid_getaliasid(g) >= 0)
            continue;

        /*
         * skip this grid if there are no model variables associated with it
         */
        for (i = 0; i < model_getnvar(m); ++i)
            if (model_getvargridid(m, i) == gid)
                break;
        if (i == model_getnvar(m))
            continue;

        enkf_printf("    calculating transforms for %s:\n", gridname);

        grid_getsize(g, &mni, &mnj, NULL);
        /*
         * a treatment for unstructured grids
         */
        if (mnj <= 0) {
            mnj = mni;
            mni = 1;
        }

        /*
         * work out how to cycle j 
         */
        nj = (mnj - 1) / stride + 1;
        jiter = malloc(nj * sizeof(int));
        for (j = 0, i = 0; j < nj; ++j, i += stride)
            jiter[j] = i;
        ni = (mni - 1) / stride + 1;
        iiter = malloc(ni * sizeof(int));
        for (i = 0, j = 0; i < ni; ++i, j += stride)
            iiter[i] = j;

        distribute_iterations(0, nj - 1, nprocesses, rank, "      ");

        if (das->mode == MODE_ENKF || das->mode == MODE_HYBRID) {
#if !defined(TW_VIAFILE)
            if (rank == 0)
                nc_createtransforms(das, gid, nj, ni, stride);
#else
            if (rank == 0)
                dir_createifabsent(DIRNAME_TMP);
#if defined(MPI)
            MPI_Barrier(MPI_COMM_WORLD);
#endif
            nc_createtransformstile(das, gid, ni);
#endif
            Tj = alloc3d(ni, nmem_dynamic, nmem, sizeof(float));
            T = alloc2d(nmem, nmem, sizeof(double));
        } else if (das->mode == MODE_ENOI) {
#if !defined(TW_VIAFILE)
            if (rank == 0)
                nc_createtransforms(das, gid, nj, ni, stride);
#else
            if (rank == 0)
                dir_createifabsent(DIRNAME_TMP);
#if defined(MPI)
            MPI_Barrier(MPI_COMM_WORLD);
#endif
            nc_createtransformstile(das, gid, ni);
#endif
        } else
            enkf_quit("programming error");
        wj = alloc2d(ni, nmem, sizeof(float));

        if (rank == 0) {
            nlobs = alloc2d(nj, ni, sizeof(int));
            dfs = alloc2d(nj, ni, sizeof(float));
            srf = alloc2d(nj, ni, sizeof(float));
            pnlobs = alloc3d(obs->nobstypes, nj, ni, sizeof(int));
            pdfs = alloc3d(obs->nobstypes, nj, ni, sizeof(float));
            psrf = alloc3d(obs->nobstypes, nj, ni, sizeof(float));
        } else if (my_number_of_iterations > 0) {
            nlobs = alloc2d(my_number_of_iterations, ni, sizeof(int));
            dfs = alloc2d(my_number_of_iterations, ni, sizeof(float));
            srf = alloc2d(my_number_of_iterations, ni, sizeof(float));
            pnlobs = alloc3d(obs->nobstypes, my_number_of_iterations, ni, sizeof(int));
            pdfs = alloc3d(obs->nobstypes, my_number_of_iterations, ni, sizeof(float));
            psrf = alloc3d(obs->nobstypes, my_number_of_iterations, ni, sizeof(float));
        }

        jpool = malloc(nj * sizeof(int));
        if (rank == 0) {
            for (j = 0; j < nj; ++j)
                jpool[j] = j;
            if (nprocesses > 1) {
                /*
                 * If writing of T is organised so that each process writes
                 * independently, then it may have sense to distribute the
                 * iterations randomly, so that the total time taken by each
                 * process is approximately equal. Currently, only the master
                 * writes to T*.nc, therefore shuffling rows makes no
                 * difference.
                 */
#if defined(SHUFFLE_ROWS)
                shuffle(nj, jpool);
#else
#if !defined(TW_VIAFILE)
                group_iterations(nj, jpool);
#endif
#endif
            }
        }
#if defined(MPI)
        {
            int ierror = MPI_Bcast(jpool, nj, MPI_INT, 0, MPI_COMM_WORLD);

            assert(ierror == MPI_SUCCESS);
        }
#endif
        /*
         * main cycle
         */
        enkf_printf("      main cycle for %s (%d x %d local analyses):\n", gridname, nj, ni);
        enkf_flush();
#if defined(MPI)
        MPI_Barrier(MPI_COMM_WORLD);    /* (to sync the logs) */
#endif
        for (jj = my_first_iteration; jj <= my_last_iteration; ++jj) {
            /*
             * j is "j" index of the actual (physical) grid (0 to mnj);
             * jpool[jj] is "j" index of the (strided) transforms subrid (0 to
             * nj)
             */
            j = jiter[jpool[jj]];
            if (enkf_verbose)
                printf("        j = %d (%d: %d: %.1f%%)\n", j, rank, jj, 100.0 * (double) (jj - my_first_iteration + 1) / (double) my_number_of_iterations);
            fflush(stdout);

            if (rank == 0)
                jjj = jpool[jj];
            else
                jjj = jj - my_first_iteration;

            for (ii = 0; ii < ni; ++ii) {
                int ploc = 0;   /* `nlobs' already engaged, using another
                                 * name */

                double* sloc = NULL;
                int* plobs = NULL;
                double** Sloc = NULL;
                double** G = NULL;
                double** M = NULL;

#if !defined(MINIMISE_ALLOC)
                int* lobs = NULL;
                double* lcoeffs = NULL;
#endif
                int e, o;

                i = iiter[ii];

                {
                    int ij[2] = { i, j };
                    double lon, lat;

                    /*
                     * for unstructured grids -- restore the natural order
                     */
                    if (mni == 1) {
                        ij[0] = j;
                        ij[1] = 0;
                    }

                    grid_ij2xy(g, ij, &lon, &lat);

#if defined(MINIMISE_ALLOC)
                    obs_findlocal(obs, lon, lat, grid_isgeographic(g), grid_getdomainname(g), &ploc, &lobs, &lcoeffs, &ploc_allocated2);
#else
                    obs_findlocal(obs, lon, lat, grid_isgeographic(g), grid_getdomainname(g), &ploc, &lobs, &lcoeffs, NULL);
#endif
                    assert(ploc >= 0 && ploc <= obs->nobs);
                }
                if (ploc > stats.nlobs_max)
                    stats.nlobs_max = ploc;
                stats.nlobs_sum += ploc;
                stats.ncell++;

                if (ploc == 0) {
                    if (T != NULL) {
                        memset(Tj[ii][0], 0, nmem_dynamic * nmem * sizeof(float));
                        for (e = 0; e < nmem_dynamic; ++e)
                            Tj[ii][e][e] = (float) 1.0;
                    }
                    memset(wj[ii], 0, nmem * sizeof(float));

                    nlobs[jjj][ii] = 0;
                    dfs[jjj][ii] = 0.0;
                    srf[jjj][ii] = 0.0;
                    for (ot = 0; ot < obs->nobstypes; ot++) {
                        pnlobs[ot][jjj][ii] = 0;
                        pdfs[ot][jjj][ii] = 0.0;
                        psrf[ot][jjj][ii] = 0.0;
                    }
                    continue;
                }
#if defined(MINIMISE_ALLOC)
                /*
                 * prepare_calcs() sets Sloc, G and M matrices while trying
                 * to minimise the number of dynamic allocations
                 */
                prepare_calcs(ploc, nmem, &sloc, &plobs, &Sloc, &G, &M);
#else
                Sloc = alloc2d(nmem, ploc, sizeof(double));
                G = alloc2d(ploc, nmem, sizeof(double));
                sloc = malloc(ploc * sizeof(double));
                plobs = malloc(ploc * sizeof(int));
#endif

                for (e = 0; e < nmem; ++e) {
                    float* Se = das->S[e];
                    double* Sloce = Sloc[e];

                    for (o = 0; o < ploc; ++o)
                        Sloce[o] = (double) Se[lobs[o]] * lcoeffs[o];
                }
                for (o = 0; o < ploc; ++o)
                    sloc[o] = das->s_f[lobs[o]] * lcoeffs[o];

                if (das->mode == MODE_ENOI) {
                    calc_G(nmem, ploc, M, Sloc, G);
                    calc_w(nmem, ploc, G, sloc, w);
                } else if (das->scheme == SCHEME_DENKF) {
                    calc_wT_denkf(nmem, nmem_dynamic, ploc, sloc, Sloc, Sloc, M, G, w, T);
                } else if (das->scheme == SCHEME_ETKF) {
                    calc_wT_etkf(nmem, nmem_dynamic, ploc, sloc, Sloc, Sloc, M, G, w, T);
                } else
                    enkf_quit("programming error");

                if (das->mode == MODE_ENKF || das->mode == MODE_HYBRID) {
                    int e1, e2;

                    /*
                     * incorporate "relaxation to prior" T = I + alpha * (T - I)
                     */
                    for (e1 = 0; e1 < nmem_dynamic; ++e1) {
                        double* Ti = T[e1];

                        Ti[e1] -= 1.0;
                        for (e2 = 0; e2 < nmem; ++e2)
                            Ti[e2] *= das->alpha;
                        Ti[e1] += 1.0;
                    }

                    /*
                     * convert T to float and store in Tj 
                     */
                    for (e1 = 0; e1 < nmem_dynamic; ++e1)
                        for (e2 = 0; e2 < nmem; ++e2)
                            Tj[ii][e1][e2] = (float) T[e1][e2];
                }
                /*
                 * convert w to float and store in wj 
                 */
                for (e = 0; e < nmem; ++e)
                    wj[ii][e] = (float) w[e];

                nlobs[jjj][ii] = ploc;  /* (ploc > 0 here) */
                dfs[jjj][ii] = traceprod(0, 0, ploc, nmem, G, Sloc, 1);
                srf[jjj][ii] = sqrt(traceprod(0, 1, nmem, ploc, Sloc, Sloc, 1) / dfs[jjj][ii]) - 1.0;
                for (ot = 0; ot < obs->nobstypes; ++ot) {
                    int p = 0;

                    for (o = 0; o < ploc; ++o) {
                        if (obs->data[lobs[o]].type == ot) {
                            plobs[p] = o;
                            p++;
                        }
                    }
                    pnlobs[ot][jjj][ii] = p;
                    if (p == 0) {
                        pdfs[ot][jjj][ii] = 0.0;
                        psrf[ot][jjj][ii] = 0.0;
                    } else {
                        /*
                         * it is used below that local observations in array
                         * plobs are continuous by observation type
                         */
#if 0
                        double** pS = alloc2d(nmem, p, sizeof(double));
                        double** pG = &G[plobs[0]];

                        for (e = 0; e < nmem; ++e)
                            for (o = 0; o < p; ++o)
                                pS[e][o] = Sloc[e][plobs[o]];
                        pdfs[ot][jjj][ii] = traceprod(0, 0, p, nmem, pG, pS, 1);
                        if (pdfs[ot][jjj][ii] > DFS_MIN)
                            psrf[ot][jjj][ii] = sqrt(traceprod(0, 1, nmem, p, pS, pS, 1) / pdfs[ot][jjj][ii]) - 1.0;
                        else
                            psrf[ot][jjj][ii] = 0.0;
#else
                        double** pS = malloc(nmem * sizeof(double*));
                        double** pG = &G[plobs[0]];

                        for (e = 0; e < nmem; ++e)
                            pS[e] = &Sloc[e][plobs[0]];
                        pdfs[ot][jjj][ii] = traceprod(0, 0, p, nmem, pG, pS, 0);
                        if (pdfs[ot][jjj][ii] > DFS_MIN)
                            psrf[ot][jjj][ii] = sqrt(traceprod(0, 1, nmem, p, pS, pS, 0) / pdfs[ot][jjj][ii]) - 1.0;
                        else
                            psrf[ot][jjj][ii] = 0.0;
#endif
                        free(pS);
                    }
                }

                if (ploc >= nmem)
                    stats.n_inv_ens++;
                else
                    stats.n_inv_obs++;

#if !defined(MINIMISE_ALLOC)
                free(G);
                free(Sloc);
                free(sloc);
                free(plobs);
                free(lobs);
                free(lcoeffs);
#endif
            }                   /* for ii */

#if defined(MPI)
            if (das->mode == MODE_ENKF || das->mode == MODE_HYBRID) {
#if !defined(TW_VIAFILE)
                if (rank > 0) {
                    if (my_number_of_iterations > 0) {
                        int ierror;

                        ierror = MPI_Send(Tj[0][0], ni * nmem_dynamic * nmem, MPI_FLOAT, 0, jj, MPI_COMM_WORLD);
                        assert(ierror == MPI_SUCCESS);
                        ierror = MPI_Send(wj[0], ni * nmem, MPI_FLOAT, 0, jj, MPI_COMM_WORLD);
                        assert(ierror == MPI_SUCCESS);
                    }
                } else {
                    int r, ierror;

                    /*
                     * write own results 
                     */
                    nc_writetransforms(das, gid, jpool[jj], ni, Tj[0][0], wj[0]);
                    /*
                     * collect and write results from slaves 
                     */
                    for (r = 1; r < nprocesses; ++r) {
                        if (jj > last_iteration[r] - first_iteration[r])
                            continue;
                        /*
                         * (recall that my_number_of_iterations <=
                         * number_of_iterations[0]) 
                         */
                        ierror = MPI_Recv(Tj[0][0], ni * nmem_dynamic * nmem, MPI_FLOAT, r, first_iteration[r] + jj, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        assert(ierror == MPI_SUCCESS);
                        ierror = MPI_Recv(wj[0], ni * nmem, MPI_FLOAT, r, first_iteration[r] + jj, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        assert(ierror == MPI_SUCCESS);
                        nc_writetransforms(das, gid, jpool[first_iteration[r] + jj], ni, Tj[0][0], wj[0]);
                    }
                }
#else                           /* TW_VIAFILE */
                nc_writetransformstile(das, gid, jpool[jj], ni, Tj[0][0], wj[0]);
#endif
            } else if (das->mode == MODE_ENOI) {
#if !defined(TW_VIAFILE)
                if (rank > 0) {
                    if (my_number_of_iterations > 0) {
                        int ierror = MPI_Send(wj[0], ni * nmem, MPI_FLOAT, 0, jj, MPI_COMM_WORLD);

                        assert(ierror == MPI_SUCCESS);
                    }
                } else {
                    int r, ierror;

                    /*
                     * write own results 
                     */
                    nc_writetransforms(das, gid, jpool[jj], ni, NULL, wj[0]);
                    /*
                     * collect and write results from slaves 
                     */
                    for (r = 1; r < nprocesses; ++r) {
                        if (jj > last_iteration[r] - first_iteration[r])
                            continue;
                        /*
                         * (recall that my_number_of_iterations <=
                         * number_of_iterations[0]) 
                         */
                        ierror = MPI_Recv(wj[0], ni * nmem, MPI_FLOAT, r, first_iteration[r] + jj, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                        assert(ierror == MPI_SUCCESS);
                        nc_writetransforms(das, gid, jpool[first_iteration[r] + jj], ni, NULL, wj[0]);
                    }
                }
#else                           /* TW_VIAFILE */
                nc_writetransformstile(das, gid, jpool[jj], ni, NULL, wj[0]);
#endif
            }
#else                           /* no MPI */
            if (das->mode == MODE_ENKF || das->mode == MODE_HYBRID)
                nc_writetransforms(das, gid, jpool[jj], ni, Tj[0][0], wj[0]);
            else if (das->mode == MODE_ENOI)
                nc_writetransforms(das, gid, jpool[jj], ni, NULL, wj[0]);
#endif                          /* if defined(MPI) */
        }                       /* for jj */

#if defined(TW_VIAFILE)
#if defined(MPI)
        MPI_Barrier(MPI_COMM_WORLD);
#endif
        if (rank == 0) {
            nc_assembletransforms(das, gid, nj, ni, stride);
            dir_rmifexists(DIRNAME_TMP);
        }
#endif
#if defined(MPI)
        MPI_Barrier(MPI_COMM_WORLD);
#endif
        enkf_printf("    finished calculating transforms for %s\n", gridname);

        enkf_flush();

#if defined(MPI)
        /*
         * collect stats on master 
         */
        if (rank > 0) {
            if (my_number_of_iterations > 0) {
                int ierror;

                ierror = MPI_Send(nlobs[0], my_number_of_iterations * ni, MPI_INT, 0, 1, MPI_COMM_WORLD);
                assert(ierror == MPI_SUCCESS);
                ierror = MPI_Send(dfs[0], my_number_of_iterations * ni, MPI_FLOAT, 0, 2, MPI_COMM_WORLD);
                assert(ierror == MPI_SUCCESS);
                ierror = MPI_Send(srf[0], my_number_of_iterations * ni, MPI_FLOAT, 0, 3, MPI_COMM_WORLD);
                assert(ierror == MPI_SUCCESS);
                ierror = MPI_Send(pnlobs[0][0], obs->nobstypes * my_number_of_iterations * ni, MPI_INT, 0, 4, MPI_COMM_WORLD);
                assert(ierror == MPI_SUCCESS);
                ierror = MPI_Send(pdfs[0][0], obs->nobstypes * my_number_of_iterations * ni, MPI_FLOAT, 0, 5, MPI_COMM_WORLD);
                assert(ierror == MPI_SUCCESS);
                ierror = MPI_Send(psrf[0][0], obs->nobstypes * my_number_of_iterations * ni, MPI_FLOAT, 0, 6, MPI_COMM_WORLD);
                assert(ierror == MPI_SUCCESS);
            }
        } else {
            int ierror;
            int** buffer_nlobs;
            float** buffer_dfs; /* also used for srf */
            int*** buffer_pnlobs;
            float*** buffer_pdfs;       /* also used for psrf */
            int r;

            for (r = 1; r < nprocesses; ++r) {
                if (number_of_iterations[r] == 0)
                    continue;

                buffer_nlobs = alloc2d(number_of_iterations[r], ni, sizeof(int));
                ierror = MPI_Recv(buffer_nlobs[0], number_of_iterations[r] * ni, MPI_INT, r, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                assert(ierror == MPI_SUCCESS);
                for (jj = first_iteration[r], jjj = 0; jj <= last_iteration[r]; ++jj, ++jjj) {
                    j = jpool[jj];
                    memcpy(nlobs[j], buffer_nlobs[jjj], ni * sizeof(int));
                }
                free(buffer_nlobs);

                buffer_dfs = alloc2d(number_of_iterations[r], ni, sizeof(float));
                ierror = MPI_Recv(buffer_dfs[0], number_of_iterations[r] * ni, MPI_FLOAT, r, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                assert(ierror == MPI_SUCCESS);
                for (jj = first_iteration[r], jjj = 0; jj <= last_iteration[r]; ++jj, ++jjj) {
                    j = jpool[jj];
                    memcpy(dfs[j], buffer_dfs[jjj], ni * sizeof(float));
                }

                ierror = MPI_Recv(buffer_dfs[0], number_of_iterations[r] * ni, MPI_FLOAT, r, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                assert(ierror == MPI_SUCCESS);
                for (jj = first_iteration[r], jjj = 0; jj <= last_iteration[r]; ++jj, ++jjj) {
                    j = jpool[jj];
                    memcpy(srf[j], buffer_dfs[jjj], ni * sizeof(float));
                }
                free(buffer_dfs);

                buffer_pnlobs = alloc3d(obs->nobstypes, number_of_iterations[r], ni, sizeof(int));
                ierror = MPI_Recv(buffer_pnlobs[0][0], obs->nobstypes * number_of_iterations[r] * ni, MPI_INT, r, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                assert(ierror == MPI_SUCCESS);
                for (ot = 0; ot < obs->nobstypes; ot++) {
                    for (jj = first_iteration[r], jjj = 0; jj <= last_iteration[r]; ++jj, ++jjj) {
                        j = jpool[jj];
                        memcpy(pnlobs[ot][j], buffer_pnlobs[ot][jjj], ni * sizeof(int));
                    }
                }
                free(buffer_pnlobs);

                buffer_pdfs = alloc3d(obs->nobstypes, number_of_iterations[r], ni, sizeof(float));
                ierror = MPI_Recv(buffer_pdfs[0][0], obs->nobstypes * number_of_iterations[r] * ni, MPI_FLOAT, r, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                assert(ierror == MPI_SUCCESS);
                for (ot = 0; ot < obs->nobstypes; ot++) {
                    for (jj = first_iteration[r], jjj = 0; jj <= last_iteration[r]; ++jj, ++jjj) {
                        j = jpool[jj];
                        memcpy(pdfs[ot][j], buffer_pdfs[ot][jjj], ni * sizeof(float));
                    }
                }

                ierror = MPI_Recv(buffer_pdfs[0][0], obs->nobstypes * number_of_iterations[r] * ni, MPI_FLOAT, r, 6, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                assert(ierror == MPI_SUCCESS);
                for (ot = 0; ot < obs->nobstypes; ot++) {
                    for (jj = first_iteration[r], jjj = 0; jj <= last_iteration[r]; ++jj, ++jjj) {
                        j = jpool[jj];
                        memcpy(psrf[ot][j], buffer_pdfs[ot][jjj], ni * sizeof(float));
                    }
                }
                free(buffer_pdfs);
            }
        }                       /* rank == 0 */
#endif                          /* MPI */
        if (rank == 0) {
            char fname_diag[MAXSTRLEN];

            das_getfname_diag(das, gid, fname_diag);

            nc_writediag(das, fname_diag, obs->nobstypes, nj, ni, stride, nlobs, dfs, srf, pnlobs, pdfs, psrf);
        }

        if (das->mode == MODE_ENKF || das->mode == MODE_HYBRID) {
            free(Tj);
            free(T);
        }
        free(wj);

#if defined(MPI)
        /*
         * merge stats for the report
         */
        if (rank > 0) {
            if (my_number_of_iterations > 0) {
                int ierror = MPI_Send(&stats, sizeof(stats) / sizeof(int), MPI_INT, 0, 99, MPI_COMM_WORLD);

                assert(ierror == MPI_SUCCESS);
            }
        } else {
            int r;

            for (r = 1; r < nprocesses; ++r) {
                calcstats morestats;
                int ierror;

                if (number_of_iterations[r] == 0)
                    continue;

                ierror = MPI_Recv(&morestats, sizeof(stats) / sizeof(int), MPI_INT, r, 99, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                assert(ierror == MPI_SUCCESS);

                stats.nlobs_sum += morestats.nlobs_sum;
                if (morestats.nlobs_max > stats.nlobs_max)
                    stats.nlobs_max = morestats.nlobs_max;
                stats.n_inv_obs += morestats.n_inv_obs;
                stats.n_inv_ens += morestats.n_inv_ens;
                stats.ncell += morestats.ncell;
            }
        }
#endif

        enkf_printf("    summary stats on %s:\n", gridname);

        enkf_printf("      # of local analyses = %d\n", stats.ncell);
        enkf_printf("      average # of local obs = %.1f\n", (double) stats.nlobs_sum / (double) stats.ncell);
        enkf_printf("      # of inversions in obs space = %d\n", stats.n_inv_obs);
        enkf_printf("      # of inversions in ens space = %d\n", stats.n_inv_ens);
        enkf_printf("");

        if (my_number_of_iterations > 0) {
            free(nlobs);
            free(dfs);
            free(srf);
            free(pnlobs);
            free(pdfs);
            free(psrf);
        }
        free(jiter);
        free(iiter);
        free(jpool);
    }                           /* for gid */

    /*
     * (this block can be put outside the grid loop because the storage sizes
     * do not depend on the grid size)
     */
#if defined(MINIMISE_ALLOC)
    if (storage != NULL) {
        free(storage);
        storage = NULL;
        ploc_allocated1 = 0;
    }
    if (lobs != NULL) {
        free(lobs);
        free(lcoeffs);
        lobs = NULL;
        lcoeffs = NULL;
        ploc_allocated2 = 0;
    }
#endif
    free(w);
}

#endif

/** Calculates transforms for pointlogs. This is done separately from
 * das_calctransfroms() to avoid interpolation related problems. At the moment
 * the procedure cycles serially through the pointlogs, but it would be trivial
 * to parallelise if needed.
 */
void das_calcpointlogtransforms(dasystem* das)
{
    model* m = das->m;
    int nmem = das->nmem;
    int nmem_dynamic = das->nmem_dynamic;
    observations* obs = das->obs;
    double* w = malloc(nmem * sizeof(double));
    double** T = NULL;
    int geographic = grid_isgeographic(model_getgridbyid(m, 0));
    int plogid, e, o, gid;

    if (das->s_mode == S_MODE_HA_f)
        das_standardise(das);

    if (das->mode == MODE_ENKF || das->mode == MODE_HYBRID)
        T = alloc2d(nmem, nmem, sizeof(double));

    for (plogid = 0; plogid < das->nplog; ++plogid) {
        pointlog* plog = &das->plogs[plogid];

#if !defined(MINIMISE_ALLOC)
        int* lobs = NULL;
        double* lcoeffs = NULL;
#endif
        int ploc = 0;

        enkf_printf("    calculating transforms for log point (%.3f,%.3f):", plog->lon, plog->lat);

        /*
         * find all (for all domains) local obs
         */
        obs_findlocal(obs, plog->lon, plog->lat, geographic, NULL, &ploc, &lobs, &lcoeffs, NULL);
        assert(ploc >= 0 && ploc <= obs->nobs);
        enkf_printf(" %d obs\n", ploc);
        /*
         * create a pointlog output file and write all local obs to it
         */
        plog_create(das, plogid, ploc, lobs, lcoeffs);
        if (ploc > 0) {
            free(lobs);
            free(lcoeffs);
            lobs = NULL;
            lcoeffs = NULL;
        }

        /*
         * calculate and write transforms for each grid
         */
        for (gid = 0; gid < model_getngrid(m); ++gid) {
            void* g = model_getgridbyid(m, gid);
            double* sloc = NULL;
            double** Sloc = NULL;
            double** G = NULL;

            if (plog->gridid >= 0 && plog->gridid != gid)
                continue;
            if (grid_getaliasid(g) >= 0)
                continue;

            ploc = 0;
            obs_findlocal(obs, plog->lon, plog->lat, geographic, grid_getdomainname(g), &ploc, &lobs, &lcoeffs, NULL);

            if (ploc == 0) {
                memset(w, 0, nmem * sizeof(double));
                if (T != NULL) {
                    memset(T[0], 0, nmem * nmem * sizeof(double));
                    for (e = 0; e < nmem; ++e)
                        T[e][e] = 1.0;
                }
            } else {
                Sloc = alloc2d(nmem, ploc, sizeof(double));
                G = alloc2d(ploc, nmem, sizeof(double));
                sloc = malloc(ploc * sizeof(double));

                for (e = 0; e < nmem; ++e) {
                    float* Se = das->S[e];
                    double* Sloce = Sloc[e];

                    for (o = 0; o < ploc; ++o)
                        Sloce[o] = (double) Se[lobs[o]] * lcoeffs[o];
                }
                for (o = 0; o < ploc; ++o)
                    sloc[o] = das->s_f[lobs[o]] * lcoeffs[o];

                if (das->mode == MODE_ENOI) {
                    calc_G(nmem, ploc, NULL, Sloc, G);
                    calc_w(nmem, ploc, G, sloc, w);
                } else if (das->scheme == SCHEME_DENKF) {
                    calc_wT_denkf(nmem, nmem_dynamic, ploc, sloc, Sloc, Sloc, NULL, G, w, T);
                } else if (das->scheme == SCHEME_ETKF) {
                    calc_wT_etkf(nmem, nmem_dynamic, ploc, sloc, Sloc, Sloc, NULL, G, w, T);
                } else
                    enkf_quit("programming error");

                if (das->mode == MODE_ENKF || das->mode == MODE_HYBRID) {
                    int e1, e2;

                    for (e1 = 0; e1 < nmem_dynamic; ++e1) {
                        double* Ti = T[e1];

                        Ti[e1] -= 1.0;
                        for (e2 = 0; e2 < nmem; ++e2)
                            Ti[e2] *= das->alpha;
                        Ti[e1] += 1.0;
                    }
                }
            }

            enkf_printf("    writing transforms for point (%.3f,%.3f) on grid \"%s\":", plog->lon, plog->lat, grid_getname(g));
            plog_writetransform(das, plogid, gid, ploc, sloc, (ploc == 0) ? NULL : Sloc[0], w, (T != NULL) ? T[0] : NULL);
            enkf_printf("\n");

            if (ploc > 0) {
                free(G);
                free(sloc);
                free(Sloc);
                free(lobs);
                free(lcoeffs);
                lobs = NULL;
                lcoeffs = NULL;
            }
        }
    }

    if (T != NULL)
        free(T);
    free(w);
}
