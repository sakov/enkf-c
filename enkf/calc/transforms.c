/******************************************************************************
 *
 * File:        transforms.c        
 *
 * Created:     11/2013
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description: Contains procedures for calculating ensemble transforms.
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

/**
 */
static void nc_createtransforms(dasystem* das, int gridid, size_t nj, size_t ni, int stride, int* ncid, int* varid_T, int* varid_w)
{
    char fname[MAXSTRLEN];
    int dimids[4];

    assert(rank == 0);

    das_getfname_transforms(das, gridid, fname);
    enkf_printf("      creating empty file \"%s\":\n", fname);
    enkf_flush();
    ncw_create(fname, NC_CLOBBER | NC_NOFILL | das->ncformat, ncid);
    ncw_def_dim(*ncid, "nj", nj, &dimids[0]);
    ncw_def_dim(*ncid, "ni", ni, &dimids[1]);
    if (das->mode == MODE_ENKF || das->mode == MODE_HYBRID) {
        ncw_def_dim(*ncid, "m_dyn", das->nmem_dynamic, &dimids[2]);
        ncw_def_dim(*ncid, "m", das->nmem, &dimids[3]);
        ncw_def_var(*ncid, "T", NC_FLOAT, 4, dimids, varid_T);
        dimids[2] = dimids[3];
    } else
        ncw_def_dim(*ncid, "m", das->nmem, &dimids[2]);
    ncw_def_var(*ncid, "w", NC_FLOAT, 3, dimids, varid_w);
    ncw_put_att_int(*ncid, NC_GLOBAL, "stride", 1, &stride);
    ncw_put_att_text(*ncid, NC_GLOBAL, "grid_name", grid_getname(model_getgridbyid(das->m, gridid)));
#if defined(DEFLATE_ALL)
    if (das->nccompression > 0)
        ncw_def_deflate(*ncid, 0, 1, das->nccompression);
#endif
    ncw_enddef(*ncid);
}

/**
 */
static void nc_writetransforms(dasystem* das, int ncid, int j, int ni, int varid_T, float* Tj, int varid_w, float* wj)
{
    size_t start[4], count[4];

    assert(rank == 0);

    start[0] = j;
    start[1] = 0;
    start[2] = 0;
    start[3] = 0;

    count[0] = 1;
    count[1] = ni;
    count[2] = das->nmem_dynamic;
    count[3] = das->nmem;

    if (Tj != NULL)
        ncw_put_vara_float(ncid, varid_T, start, count, Tj);

    count[2] = das->nmem;
    ncw_put_vara_float(ncid, varid_w, start, count, wj);
}

#if defined(TW_VIAFILE)
/**
 */
static void das_getfname_transformstile(dasystem* das, int gridid, int r, char fname[])
{
    if (model_getngrid(das->m) == 1)
        snprintf(fname, MAXSTRLEN, "%s/%s-%03d.nc", DIRNAME_TMP, FNAMEPREFIX_TRANSFORMS, r);
    else
        snprintf(fname, MAXSTRLEN, "%s/%s-%d-%03d.nc", DIRNAME_TMP, FNAMEPREFIX_TRANSFORMS, gridid, r);
}

/**
 */
static void nc_createtransformstile(dasystem* das, int gridid, int ni, int* ncid, int* varid_T, int* varid_w)
{
    char fname[MAXSTRLEN];
    int dimids[4];

    das_getfname_transformstile(das, gridid, rank, fname);

    ncw_create(fname, NC_CLOBBER | NC_NOFILL | das->ncformat, ncid);
    ncw_def_dim(*ncid, "nj", my_number_of_iterations, &dimids[0]);
    ncw_def_dim(*ncid, "ni", ni, &dimids[1]);
    if (das->mode == MODE_ENKF || das->mode == MODE_HYBRID) {
        ncw_def_dim(*ncid, "m_dyn", das->nmem_dynamic, &dimids[2]);
        ncw_def_dim(*ncid, "m", das->nmem, &dimids[3]);
        ncw_def_var(*ncid, "T", NC_FLOAT, 4, dimids, varid_T);
        dimids[2] = dimids[3];
    } else
        ncw_def_dim(*ncid, "m", das->nmem, &dimids[2]);
    ncw_def_var(*ncid, "w", NC_FLOAT, 3, dimids, varid_w);
    ncw_put_att_int(*ncid, NC_GLOBAL, "j1", 1, &my_first_iteration);
#if defined(DEFLATE_ALL)
    if (das->nccompression > 0)
        ncw_def_deflate(*ncid, 0, 1, das->nccompression);
#endif
    ncw_enddef(*ncid);
}

/**
 */
static void nc_writetransformstile(dasystem* das, int iter, int ni, float* Tj, float* wj, int ncid, int varid_T, int varid_w)
{
    size_t start[4], count[4];

    start[0] = iter;
    start[1] = 0;
    start[2] = 0;
    start[3] = 0;

    count[0] = 1;
    count[1] = ni;
    count[2] = das->nmem_dynamic;
    count[3] = das->nmem;

    ncw_put_vara_float(ncid, varid_T, start, count, Tj);

    count[2] = das->nmem;
    ncw_put_vara_float(ncid, varid_w, start, count, wj);
}

static void nc_assembletransforms(dasystem* das, int gridid, size_t nj, size_t ni, int stride)
{
    char fname[MAXSTRLEN];
    int ncid, varid_T, varid_w;
    int r;
    float* v_T = NULL;
    float* v_w = NULL;

    assert(rank == 0);

    das_getfname_transforms(das, gridid, fname);
    nc_createtransforms(das, gridid, nj, ni, stride, &ncid, &varid_T, &varid_w);

    enkf_printf("      assembling \"%s\":", fname);

    v_T = malloc(ni * number_of_iterations[0] * das->nmem_dynamic * das->nmem * sizeof(float));
    v_w = malloc(ni * number_of_iterations[0] * das->nmem * sizeof(float));

    for (r = 0; r < nprocesses; ++r) {
        char fname_tile[MAXSTRLEN];
        int ncid_tile;
        int varid_T_tile, varid_w_tile;
        size_t start[4], count[4];

        das_getfname_transformstile(das, gridid, r, fname_tile);

        ncw_open(fname_tile, NC_NOWRITE, &ncid_tile);
        ncw_inq_varid(ncid_tile, "T", &varid_T_tile);
        ncu_readvarfloat(ncid_tile, varid_T_tile, ni * number_of_iterations[r] * das->nmem * das->nmem_dynamic, v_T);
        ncw_inq_varid(ncid_tile, "w", &varid_w_tile);
        ncu_readvarfloat(ncid_tile, varid_w_tile, ni * number_of_iterations[r] * das->nmem, v_w);
        ncw_close(ncid_tile);
        file_delete(fname_tile);

        start[0] = first_iteration[r];
        start[1] = 0;
        start[2] = 0;
        start[3] = 0;

        count[0] = number_of_iterations[r];
        count[1] = ni;
        count[2] = das->nmem_dynamic;
        count[3] = das->nmem;

        ncw_put_vara_float(ncid, varid_T, start, count, v_T);
        count[2] = das->nmem;
        ncw_put_vara_float(ncid, varid_w, start, count, v_w);
        enkf_printf(".");
        enkf_flush();
    }
    ncw_close(ncid);
    enkf_printf("\n");
    enkf_flush();

    free(v_T);
    free(v_w);
}
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
static void nc_writediag(dasystem* das, char fname[], int nobstypes, int nj, int ni, int stride, int** nlobs, float** dfs, float** srf, int*** pnlobs, float*** pdfs, float*** psrf)
{
    int ncid;
    int dimids[3];
    int varid_nlobs, varid_dfs, varid_srf, varid_pnlobs, varid_pdfs, varid_psrf;

    assert(rank == 0);

    enkf_printf("    writing stats to \"%s\":\n", fname);
    ncw_create(fname, NC_CLOBBER | das->ncformat, &ncid);
    ncw_def_dim(ncid, "nobstypes", nobstypes, &dimids[0]);
    ncw_def_dim(ncid, "nj", nj, &dimids[1]);
    ncw_def_dim(ncid, "ni", ni, &dimids[2]);
    ncw_put_att_int(ncid, NC_GLOBAL, "stride", 1, &stride);
    ncw_def_var(ncid, "nlobs", NC_INT, 2, &dimids[1], &varid_nlobs);
    ncw_def_var(ncid, "dfs", NC_FLOAT, 2, &dimids[1], &varid_dfs);
    ncw_def_var(ncid, "srf", NC_FLOAT, 2, &dimids[1], &varid_srf);
    ncw_def_var(ncid, "pnlobs", NC_INT, 3, dimids, &varid_pnlobs);
    ncw_def_var(ncid, "pdfs", NC_FLOAT, 3, dimids, &varid_pdfs);
    ncw_def_var(ncid, "psrf", NC_FLOAT, 3, dimids, &varid_psrf);
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
         * matrices in invsqrtm2() for scheme = ETKF, and just wasted for scheme
         * = DEnKF)
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
        double sfactor = grid_getsfactor(g);
        int stride = grid_getstride(g);

        int mni, mnj;
        int nj, ni;
        int* jiter = NULL;
        int* iiter = NULL;

        /*
         * ids of/in transform file
         */
        int ncid = -1;
        int varid_T = -1;
        int varid_w = -1;

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
                nc_createtransforms(das, gid, nj, ni, stride, &ncid, &varid_T, &varid_w);
#else
            if (rank == 0)
                dir_createifabsent(DIRNAME_TMP);
            MPI_Barrier(MPI_COMM_WORLD);
            nc_createtransformstile(das, gid, ni, &ncid, &varid_T, &varid_w);
#endif
            Tj = alloc3d(ni, nmem_dynamic, nmem, sizeof(float));
            T = alloc2d(nmem, nmem, sizeof(double));
        } else if (das->mode == MODE_ENOI) {
            if (rank == 0)
                nc_createtransforms(das, gid, nj, ni, stride, &ncid, NULL, &varid_w);
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
            int ierror = MPI_Bcast(jpool, nj, MPI_INTEGER, 0, MPI_COMM_WORLD);

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
                    double lon, lat;

                    grid_ij2xy(g, i, j, &lon, &lat);

#if defined(MINIMISE_ALLOC)
                    obs_findlocal(obs, lon, lat, grid_getdomainname(g), &ploc, &lobs, &lcoeffs, &ploc_allocated2);
#else
                    obs_findlocal(obs, lon, lat, grid_getdomainname(g), &ploc, &lobs, &lcoeffs, NULL);
#endif
                    assert(ploc >= 0 && ploc <= obs->nobs);
                }
                if (ploc > stats.nlobs_max)
                    stats.nlobs_max = ploc;
                stats.nlobs_sum += ploc;
                stats.ncell++;

                if (ploc == 0) {
                    if (das->mode == MODE_ENKF || das->mode == MODE_HYBRID) {
                        /*
                         * set T = I 
                         */
                        memset(Tj[ii][0], 0, nmem_dynamic * nmem * sizeof(float));
                        for (e = 0; e < nmem_dynamic; ++e)
                            Tj[ii][e][e] = (float) 1.0;
                    } else if (das->mode == MODE_ENOI)
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
                        Sloce[o] = (double) Se[lobs[o]] * lcoeffs[o] * sfactor;
                }
                for (o = 0; o < ploc; ++o)
                    sloc[o] = das->s_f[lobs[o]] * lcoeffs[o];

                if (das->mode == MODE_ENOI || das->scheme == SCHEME_DENKF) {
                    calc_G(nmem, ploc, M, Sloc, i, j, G);
                    if (das->mode == MODE_ENKF || das->mode == MODE_HYBRID)
                        calc_T_denkf(nmem, ploc, G, Sloc, T);
                } else if (das->scheme == SCHEME_ETKF)
                    calc_GT_etkf(nmem, ploc, M, Sloc, i, j, G, T);
                else
                    enkf_quit("programming error");

                calc_w(nmem, ploc, G, sloc, w);

                if (das->mode == MODE_ENKF || das->mode == MODE_HYBRID) {
                    int e1, e2;

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

                if (ploc > nmem)
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
                    nc_writetransforms(das, ncid, jpool[jj], ni, varid_T, Tj[0][0], varid_w, wj[0]);
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
                        nc_writetransforms(das, ncid, jpool[first_iteration[r] + jj], ni, varid_T, Tj[0][0], varid_w, wj[0]);
                    }
                }
#else                           /* TW_VIAFILE */
                nc_writetransformstile(das, jj - my_first_iteration, ni, Tj[0][0], wj[0], ncid, varid_T, varid_w);
#endif
            } else if (das->mode == MODE_ENOI) {
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
                    nc_writetransforms(das, ncid, jpool[jj], ni, -1, NULL, varid_w, wj[0]);
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
                        nc_writetransforms(das, ncid, jpool[first_iteration[r] + jj], ni, -1, NULL, varid_w, wj[0]);
                    }
                }
            }
#else                           /* no MPI */
            if (das->mode == MODE_ENKF || das->mode == MODE_HYBRID)
                nc_writetransforms(das, ncid, jpool[jj], ni, varid_T, Tj[0][0], varid_w, wj[0]);
            else if (das->mode == MODE_ENOI)
                nc_writetransforms(das, ncid, jpool[jj], ni, -1, NULL, varid_w, wj[0]);
#endif                          /* if defined(MPI) */
        }                       /* for jj */

#if !defined(TW_VIAFILE)
        if (rank == 0)
            ncw_close(ncid);
#else
        if (das->mode == MODE_ENKF || das->mode == MODE_HYBRID || (das->mode == MODE_ENOI && rank == 0))
            ncw_close(ncid);
        if (das->mode == MODE_ENKF || das->mode == MODE_HYBRID) {
            MPI_Barrier(MPI_COMM_WORLD);
            if (rank == 0) {
                nc_assembletransforms(das, gid, nj, ni, stride);

                dir_rmifexists(DIRNAME_TMP);
            }
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
            char fname_stats[MAXSTRLEN];

            das_getfname_stats(das, g, fname_stats);

            nc_writediag(das, fname_stats, obs->nobstypes, nj, ni, stride, nlobs, dfs, srf, pnlobs, pdfs, psrf);
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

/** Calculates transforms for pointlogs. This is done separately from
 * das_calctransfroms() to avoid interpolation related problems. At the moment
 * the procedure cycles serially through the pointlogs, but it would be trivial
 * to parallelise if needed.
 */
void das_calcpointlogtransforms(dasystem* das)
{
    model* m = das->m;
    int nmem = das->nmem;
    double* w = malloc(nmem * sizeof(double));
    observations* obs = das->obs;
    double** X5 = NULL;
    int plogid, e, o, gid;

    if (das->s_mode == S_MODE_HA_f)
        das_standardise(das);

    if (das->mode == MODE_ENKF || das->mode == MODE_HYBRID)
        X5 = alloc2d(nmem, nmem, sizeof(double));

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
#if defined(MINIMISE_ALLOC)
        obs_findlocal(obs, plog->lon, plog->lat, NULL, &ploc, &lobs, &lcoeffs, NULL);
#else
        obs_findlocal(obs, plog->lon, plog->lat, NULL, &ploc, &lobs, &lcoeffs);
#endif
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
            int i = (int) (plog->fi[gid] + 0.5);
            int j = (int) (plog->fj[gid] + 0.5);
            double sfactor = grid_getsfactor(g);
            double* sloc = NULL;
            double** Sloc = NULL;
            double** G = NULL;

            if (plog->gridid >= 0 && plog->gridid != gid)
                continue;

            ploc = 0;
#if defined(MINIMISE_ALLOC)
            obs_findlocal(obs, plog->lon, plog->lat, grid_getdomainname(g), &ploc, &lobs, &lcoeffs, NULL);
#else
            obs_findlocal(obs, plog->lon, plog->lat, grid_getdomainname(g), &ploc, &lobs, &lcoeffs);
#endif

            if (X5 != NULL)
                memset(X5[0], 0, nmem * nmem * sizeof(double));

            if (ploc == 0) {
                if (das->mode == MODE_ENKF || das->mode == MODE_HYBRID)
                    for (e = 0; e < nmem; ++e)
                        X5[e][e] = 1.0;
            } else {
                Sloc = alloc2d(nmem, ploc, sizeof(double));
                sloc = malloc(ploc * sizeof(double));
                G = alloc2d(ploc, nmem, sizeof(double));

                for (e = 0; e < nmem; ++e) {
                    float* Se = das->S[e];
                    double* Sloce = Sloc[e];

                    for (o = 0; o < ploc; ++o)
                        Sloce[o] = (double) Se[lobs[o]] * lcoeffs[o] * sfactor;
                }
                for (o = 0; o < ploc; ++o)
                    sloc[o] = das->s_f[lobs[o]] * lcoeffs[o];

                if (das->mode == MODE_ENOI || das->scheme == SCHEME_DENKF) {
                    calc_G(nmem, ploc, NULL, Sloc, i, j, G);
                    if (das->mode == MODE_ENKF || das->mode == MODE_HYBRID)
                        calc_T_denkf(nmem, ploc, G, Sloc, X5);
                } else if (das->scheme == SCHEME_ETKF)
                    calc_GT_etkf(nmem, ploc, NULL, Sloc, i, j, G, X5);
                else
                    enkf_quit("programming error");

                calc_w(nmem, ploc, G, sloc, w);
                if (das->mode == MODE_ENKF || das->mode == MODE_HYBRID)
                    calc_X5(nmem, das->alpha, w, X5);
            }

            enkf_printf("    writing log for point (%.3f,%.3f) on grid \"%s\":", plog->lon, plog->lat, grid_getname(g));
            plog_writetransform(das, plogid, gid, ploc, sloc, (ploc == 0) ? NULL : Sloc[0], (das->mode == MODE_ENKF || das->mode == MODE_HYBRID) ? X5[0] : w);
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

    if (das->mode == MODE_ENKF || das->mode == MODE_HYBRID)
        free(X5);
    free(w);
}
