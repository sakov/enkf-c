/******************************************************************************
 *
 * File:        ensobs.c        
 *
 * Created:     11/2013
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description: This file contains code for calculating forecast and analysis
 *              ensemble observations HE. The forecast ensemble observations are
 *              calculated using observation functions H from model2obs.c.
 *              They are then used for calculating ensemble transforms X5.
 *              The analysis ensemble observations are calculated by applying
 *              the ensemble transforms X5 to the forecast HE. Both forecast and
 *              analysis ensemble observations are used for calculating
 *              observation statistics in obsstats.c.
 *
 *              For the EnOI the ensemble observations are calculated as
 *              H(E) = H(x) 1' + H(A). This implies that the observation
 *              functions H are linear. For nonlinear H one needs to use
 *              H(E) = H(x 1' + A) instead, but this has little sense as the
 *              EnOI itself implies linearity
 *
 * Revisions:
 *
 *****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <math.h>
#include "distribute.h"
#include "lapack.h"
#include "definitions.h"
#include "utils.h"
#include "grid.h"
#include "model2obs.h"
#include "allhs.h"
#include "dasystem.h"

#define EPSF 1.0e-6f
#if defined(HE_VIASHMEM)
#define HE_NPROCMAX 20
#endif

/**
 */
void das_getHE(dasystem* das)
{
    observations* obs = das->obs;
    model* m = das->m;
    ENSOBSTYPE* Hx = NULL;
    int nobs = obs->nobs;
    int nmem = -1;
    int i, e;

#if defined (HE_VIASHMEM)
    ENSOBSTYPE* SS = NULL;
    MPI_Aint size;
    int ierror;
#endif

    das->s_mode = S_MODE_HE_f;
    if (nobs == 0)
        return;

    if (das->mode == MODE_ENKF || !enkf_fstatsonly) {
        das_setnmem(das);
        assert(das->nmem > 0);
        enkf_printf("    ensemble size = %d\n", das->nmem);
    } else
        das->nmem = 1;
    nmem = das->nmem;

    /*
     * ensemble observation array to be filled 
     */
    assert(das->S == NULL);
#if !defined(HE_VIASHMEM)
    das->S = alloc2d(nmem, nobs, sizeof(ENSOBSTYPE));
#else
    das->S = calloc(nmem, sizeof(void*));
    das->St = calloc(nobs, sizeof(void*));
    size = (das->sm_rank == 0) ? nmem * nobs * sizeof(ENSOBSTYPE) * 2 : 0;
    ierror = MPI_Win_allocate_shared(size, sizeof(ENSOBSTYPE), MPI_INFO_NULL, das->sm_comm, &SS, &das->sm_win);
    assert(ierror == MPI_SUCCESS);

    if (das->sm_rank == 0) {
        memset(SS, 0, size);
    } else {
        int disp_unit;
        MPI_Aint my_size;

        ierror = MPI_Win_shared_query(das->sm_win, 0, &my_size, &disp_unit, &SS);
        assert(ierror == MPI_SUCCESS);
    }
    for (i = 0; i < nmem; ++i)
        das->S[i] = &SS[i * nobs];
    for (i = 0; i < nobs; ++i)
        das->St[i] = &SS[(i + nobs) * nmem];
#endif

    distribute_iterations(0, nmem - 1, nprocesses, rank, "    ");

    if (das->mode == MODE_ENOI)
        Hx = calloc(nobs, sizeof(ENSOBSTYPE));

    for (i = 0; i < obs->nobstypes; ++i) {
        obstype* ot = &obs->obstypes[i];
        H_fn H = NULL;
        int nobs;
        int* obsids;
        char fname[MAXSTRLEN];

        enkf_printf("    %s ", ot->name);
        fflush(stdout);

        /*
         * set H
         */
        H = getH(ot->issurface, ot->hfunction);

        if (ot->isasync) {
            int t1 = get_tshift(ot->date_min, ot->async_tstep, ot->async_centred);
            int t2 = get_tshift(ot->date_max, ot->async_tstep, ot->async_centred);
            int t;

            for (t = t1; t <= t2; ++t) {
                enkf_printf("|");
                obs_find_bytypeandtime(obs, i, t, &nobs, &obsids);
                if (nobs == 0)
                    continue;

                /*
                 * for EnOI it is essential sometimes (e.g. in some bias
                 * correction cases) that the background is interpolated first
                 */
                if (das->mode == MODE_ENOI) {
                    if (enkf_obstype == OBSTYPE_VALUE) {
                        int success = model_getbgfname_async(m, das->bgdir, ot->varnames[0], ot->name, t, fname);

                        H(das, nobs, obsids, fname, -1, t, Hx);
                        enkf_printf((success) ? "A" : "S");
                        fflush(stdout);
                    } else if (enkf_obstype == OBSTYPE_INNOVATION) {
                        Hx[0] = 0;
                        enkf_printf("-");
                        fflush(stdout);
                    }
                }

                if (das->mode == MODE_ENKF || !enkf_fstatsonly) {
                    for (e = my_first_iteration; e <= my_last_iteration; ++e) {
                        int success = model_getmemberfname_async(m, das->ensdir, ot->varnames[0], ot->name, e + 1, t, fname);

                        H(das, nobs, obsids, fname, e + 1, t, das->S[e]);
                        enkf_printf((success) ? "a" : "s");
                        fflush(stdout);
                    }
                }

                free(obsids);
            }
        } else {
            obs_find_bytype(obs, i, &nobs, &obsids);
            if (nobs == 0)
                goto next;

            /*
             * for EnOI it is essential sometimes (e.g. in some bias correction
             * cases) that the background is interpolated first
             */
            if (das->mode == MODE_ENOI) {
                if (enkf_obstype == OBSTYPE_VALUE) {
                    model_getbgfname(m, das->bgdir, ot->varnames[0], fname);
                    H(das, nobs, obsids, fname, -1, INT_MAX, Hx);
                    enkf_printf("+");
                    fflush(stdout);
                } else if (enkf_obstype == OBSTYPE_INNOVATION) {
                    Hx[0] = 0;
                    enkf_printf("-");
                    fflush(stdout);
                }
            }

            if (das->mode == MODE_ENKF || !enkf_fstatsonly) {
                for (e = my_first_iteration; e <= my_last_iteration; ++e) {
                    model_getmemberfname(m, das->ensdir, ot->varnames[0], e + 1, fname);
                    H(das, nobs, obsids, fname, e + 1, INT_MAX, das->S[e]);
                    enkf_printf(".");
                    fflush(stdout);
                }
            }

            free(obsids);
        }

      next:

        enkf_printf("\n");
    }                           /* for i (over obstypes) */

#if defined(MPI)
    if (das->mode == MODE_ENKF || !enkf_fstatsonly) {
#if !defined(HE_VIAFILE)
        /*
         * communicate HE via MPI
         */
        int* recvcounts = calloc(nprocesses, sizeof(int));
        int* displs = calloc(nprocesses, sizeof(int));
        MPI_Datatype mpitype_vec_nobs;
        int ierror;

#if defined(HE_VIASHMEM)
        int ii;
#endif

        ierror = MPI_Type_contiguous(nobs, MPIENSOBSTYPE, &mpitype_vec_nobs);
        assert(ierror == MPI_SUCCESS);

#if !defined(HE_VIASHMEM)
        for (i = 0; i < nprocesses; ++i) {
            recvcounts[i] = number_of_iterations[i];
            displs[i] = first_iteration[i];
        }
        MPI_Barrier(MPI_COMM_WORLD);
        ierror = MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, das->S[0], recvcounts, displs, mpitype_vec_nobs, MPI_COMM_WORLD);
        assert(ierror == MPI_SUCCESS);
#else
        MPI_Barrier(MPI_COMM_WORLD);
        if (das->node_rank >= 0 && das->node_size > 1) {
            for (i = 0, ii = -1; i < nprocesses; ++i) {
                if (das->node_ranks[i] >= 0) {
                    ii = das->node_ranks[i];
                    displs[ii] = first_iteration[i];
                }
                assert(ii >= 0);
                recvcounts[ii] += number_of_iterations[i];
            }
            ierror = MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, das->S[0], recvcounts, displs, mpitype_vec_nobs, das->node_comm);
            assert(ierror == MPI_SUCCESS);
        }
        MPI_Barrier(MPI_COMM_WORLD);
#endif

        free(recvcounts);
        free(displs);
#else                           /* HE_VIAFILE */
        /*
         * communicate HE via file
         */
        {
            int ncid;
            int varid;
            size_t start[2], count[2];

            if (rank == 0) {
                int dimids[2];

                ncw_create(FNAME_HE, NC_CLOBBER | das->ncformat, &ncid);
                ncw_def_dim(FNAME_HE, ncid, "m", nmem, &dimids[0]);
                ncw_def_dim(FNAME_HE, ncid, "p", nobs, &dimids[1]);
                ncw_def_var(FNAME_HE, ncid, "HE", NC_FLOAT, 2, dimids, &varid);
                if (das->nccompression > 0)
                    ncw_def_deflate(ncid, 0, 1, das->nccompression);
                ncw_close(FNAME_HE, ncid);
            }
            MPI_Barrier(MPI_COMM_WORLD);

            ncw_open(FNAME_HE, NC_WRITE, &ncid);
            ncw_inq_varid(FNAME_HE, ncid, "HE", &varid);
            start[0] = my_first_iteration;
            start[1] = 0;
            count[0] = my_last_iteration - my_first_iteration + 1;
            count[1] = nobs;
            ncw_put_vara_float(FNAME_HE, ncid, varid, start, count, das->S[my_first_iteration]);
            ncw_close(FNAME_HE, ncid);
            MPI_Barrier(MPI_COMM_WORLD);

            ncw_open(FNAME_HE, NC_NOWRITE, &ncid);
            ncw_inq_varid(FNAME_HE, ncid, "HE", &varid);
            ncw_get_var_float(FNAME_HE, ncid, varid, das->S[0]);
            ncw_close(FNAME_HE, ncid);
        }
#endif
    }
#endif

    if (das->mode == MODE_ENOI) {
        /*
         * subtract ensemble mean; add background
         */
        if (!enkf_fstatsonly) {
            ENSOBSTYPE* ensmean = calloc(nobs, sizeof(ENSOBSTYPE));

            for (e = 0; e < nmem; ++e) {
                ENSOBSTYPE* Se = das->S[e];

                for (i = 0; i < nobs; ++i)
                    ensmean[i] += Se[i];
            }
            for (i = 0; i < nobs; ++i)
                ensmean[i] /= (ENSOBSTYPE) nmem;

#if defined (HE_VIASHMEM)
            if (das->sm_rank == 0)
#endif
                for (e = 0; e < nmem; ++e) {
                    ENSOBSTYPE* Se = das->S[e];

                    for (i = 0; i < nobs; ++i)
                        Se[i] += Hx[i] - ensmean[i];
                }
#if defined (HE_VIASHMEM)
            MPI_Barrier(das->sm_comm);
#endif
            free(ensmean);
        } else {
#if defined (HE_VIASHMEM)
            if (das->sm_rank == 0)
#endif
                for (e = 0; e < nmem; ++e) {
                    ENSOBSTYPE* Se = das->S[e];

                    for (i = 0; i < nobs; ++i)
                        Se[i] = Hx[i];
                }
#if defined (HE_VIASHMEM)
            MPI_Barrier(das->sm_comm);
#endif
        }
    }

    if (das->mode == MODE_ENOI)
        free(Hx);
}

/**
 */
void das_calcinnandspread(dasystem* das)
{
    observations* obs = das->obs;
    int nobs = obs->nobs;
    int nmem = das->nmem;
    int e, o;

    if (nobs == 0)
        goto finish;

    if (das->s_mode == S_MODE_HE_f) {
        if (das->s_f == NULL) {
            das->s_f = calloc(nobs, sizeof(double));
            assert(das->std_f == NULL);
            das->std_f = calloc(nobs, sizeof(double));
        } else {
            memset(das->s_f, 0, nobs * sizeof(double));
            memset(das->std_f, 0, nobs * sizeof(double));
        }

        /*
         * calculate ensemble mean observations 
         */
        for (e = 0; e < nmem; ++e) {
            ENSOBSTYPE* Se = das->S[e];

            for (o = 0; o < nobs; ++o)
                das->s_f[o] += (double) Se[o];
        }
        for (o = 0; o < nobs; ++o)
            das->s_f[o] /= (double) nmem;

        /*
         * calculate ensemble spread and innovation 
         */
#if defined(HE_VIASHMEM)
        if (das->sm_rank == 0)
#endif
            for (e = 0; e < nmem; ++e) {
                ENSOBSTYPE* Se = das->S[e];

                for (o = 0; o < nobs; ++o)
                    Se[o] -= (ENSOBSTYPE) das->s_f[o];
            }
#if defined(HE_VIASHMEM)
        MPI_Barrier(das->sm_comm);
#endif
        for (e = 0; e < nmem; ++e) {
            ENSOBSTYPE* Se = das->S[e];

            for (o = 0; o < nobs; ++o)
                das->std_f[o] += (double) (Se[o] * Se[o]);
        }

        for (o = 0; o < nobs; ++o) {
            observation* m = &obs->data[o];

            if (m->status != STATUS_OK)
                continue;
            das->std_f[o] = sqrt(das->std_f[o] / (double) (nmem - 1));
            das->s_f[o] = m->value - das->s_f[o];
            if (!isfinite(das->s_f[o]) || fabs(das->s_f[o]) > STATE_BIGNUM) {
                enkf_flush();
                enkf_printf("\n  obs # %d: ", o);
                obs_printob(obs, o);
                enkf_quit("obs # %d: y - Hx_f = %.3g, no point to continue", o, das->s_f[o]);
            }
        }

        das->s_mode = S_MODE_HA_f;
    } else if (das->s_mode == S_MODE_HE_a) {
        if (das->s_a == NULL) {
            das->s_a = calloc(nobs, sizeof(double));
            assert(das->std_a == NULL);
            das->std_a = calloc(nobs, sizeof(double));
        } else {
            memset(das->s_a, 0, nobs * sizeof(double));
            memset(das->std_a, 0, nobs * sizeof(double));
        }

        /*
         * calculate ensemble mean observations 
         */
        for (e = 0; e < nmem; ++e) {
            ENSOBSTYPE* Se = das->S[e];

            for (o = 0; o < nobs; ++o)
                das->s_a[o] += (double) Se[o];
        }
        for (o = 0; o < nobs; ++o)
            das->s_a[o] /= (double) nmem;

        /*
         * calculate ensemble spread and innovation 
         */
#if defined(HE_VIASHMEM)
        if (das->sm_rank == 0)
#endif
            for (e = 0; e < nmem; ++e) {
                ENSOBSTYPE* Se = das->S[e];

                for (o = 0; o < nobs; ++o)
                    Se[o] -= (ENSOBSTYPE) das->s_a[o];
            }
#if defined(HE_VIASHMEM)
        MPI_Barrier(das->sm_comm);
#endif
        for (e = 0; e < nmem; ++e) {
            ENSOBSTYPE* Se = das->S[e];

            for (o = 0; o < nobs; ++o)
                das->std_a[o] += (double) (Se[o] * Se[o]);
        }

        for (o = 0; o < nobs; ++o) {
            observation* m = &obs->data[o];

            if (m->status != STATUS_OK)
                continue;
            das->std_a[o] = sqrt(das->std_a[o] / (double) (nmem - 1));
            das->s_a[o] = m->value - das->s_a[o];
            if (!isfinite(das->s_a[o]) || fabs(das->s_a[o]) > STATE_BIGNUM) {
                enkf_flush();
                enkf_printf("\n  obs # %d: ", o);
                obs_printob(obs, o);
                enkf_quit("obs # %d: y - Hx_a = %d, no point to continue", o);
            }
        }

        das->s_mode = S_MODE_HA_a;
    } else
        enkf_quit("programming error");

  finish:
    if (das->s_mode == S_MODE_HE_f)
        das->s_mode = S_MODE_HA_f;
    else if (das->s_mode == S_MODE_HE_a)
        das->s_mode = S_MODE_HA_a;
}

/** Adds forecast observations and forecast ensemble spread to the observation
 * file.
 */
void das_addforecast(dasystem* das, char fname[])
{
    int ncid;
    int dimid_nobs[1];
    size_t nobs;
    int varid_Hx, varid_spread;
    double* Hx;
    int o;

    if (das->obs->nobs == 0)
        return;
    if (rank != 0)
        return;

    assert(das->s_mode == S_MODE_HA_f);

    ncw_open(fname, NC_WRITE, &ncid);
    if (ncw_var_exists(ncid, "Hx_f")) {
        enkf_printf("    Hx_f already added to \"%s\" (skipping)\n", fname);
        goto finish;
    }

    ncw_inq_dimid(ncid, "nobs", dimid_nobs);
    ncw_inq_dimlen(ncid, dimid_nobs[0], &nobs);
    assert(nobs == das->obs->nobs);
    ncw_redef(ncid);
    ncw_def_var(ncid, "Hx_f", NC_FLOAT, 1, dimid_nobs, &varid_Hx);
    ncw_def_var(ncid, "std_f", NC_FLOAT, 1, dimid_nobs, &varid_spread);
    ncw_enddef(ncid);

    Hx = calloc(nobs, sizeof(double));
    for (o = 0; o < (int) nobs; ++o)
        Hx[o] = das->obs->data[o].value - das->s_f[o];

    ncw_put_var_double(ncid, varid_Hx, Hx);
    ncw_put_var_double(ncid, varid_spread, das->std_f);

    free(Hx);

  finish:
    ncw_close(ncid);
}

/** Modifies observation error so that the increment for this observation would
 * not exceed KFACTOR * <ensemble spread> (all in observation space) after
 * assimilating this observation only. This can be viewed as an adaptive QC.
 */
void das_moderateobs(dasystem* das)
{
    observations* obs = das->obs;
    double kfactor = das->kfactor;
    double* std_new;
    int i;

    if (obs->nobs == 0)
        return;
    if (!isfinite(kfactor))
        return;

    assert(das->s_mode == S_MODE_HA_f);

    std_new = malloc(obs->nobs * sizeof(double));

    for (i = 0; i < obs->nobs; ++i) {
        observation* o = &obs->data[i];
        double svar = das->std_f[i] * das->std_f[i];
        double ovar = o->std * o->std;
        double inn = das->s_f[i];

        if (o->status != STATUS_OK) {
            std_new[i] = o->std;
            continue;
        }

        std_new[i] = sqrt(sqrt((svar + ovar) * (svar + ovar) + svar * inn * inn / kfactor / kfactor) - svar);
        if (svar > 0.0 && std_new[i] * std_new[i] / ovar > 2.0) {
            obs->nmodified++;
            obs->obstypes[o->type].nmodified++;
        }
    }

    for (i = 0; i < obs->nobs; ++i)
        obs->data[i].std = std_new[i];

    enkf_printf("    observations substantially modified:\n");
    for (i = 0; i < obs->nobstypes; ++i)
        enkf_printf("      %s    %7d (%.1f%%)\n", obs->obstypes[i].name, obs->obstypes[i].nmodified, 100.0 * (double) obs->obstypes[i].nmodified / (double) obs->obstypes[i].nobs);
    enkf_printf("      total  %7d (%.1f%%)\n", obs->nmodified, 100.0 * (double) obs->nmodified / (double) obs->nobs);

    free(std_new);
}

/** Replaces observation errors in the observation file with the modified
 * values. The original values are stored as "std_orig".
 */
void das_addmodifiederrors(dasystem* das, char fname[])
{
    int ncid;
    int dimid_nobs[1];
    size_t nobs;
    int varid_std;
    double* std;
    int i;
    double da_julday = NAN;

    if (rank != 0)
        return;

    ncw_open(fname, NC_WRITE, &ncid);
    assert(!ncw_var_exists(ncid, "std_orig"));
    ncw_inq_dimid(ncid, "nobs", dimid_nobs);
    ncw_inq_dimlen(ncid, dimid_nobs[0], &nobs);

    ncw_get_att_double(ncid, NC_GLOBAL, "DA_JULDAY", &da_julday);
    if (!enkf_noobsdatecheck && (isnan(da_julday) || fabs(das->obs->da_date - da_julday) > 1e-6))
        enkf_quit("\"observations.nc\" from a different cycle");

    ncw_redef(ncid);
    ncw_rename_var(ncid, "std", "std_orig");
    ncw_def_var(ncid, "std", NC_FLOAT, 1, dimid_nobs, &varid_std);
    ncw_enddef(ncid);

    std = malloc(nobs * sizeof(double));

    for (i = 0; i < (int) nobs; ++i)
        std[i] = das->obs->data[i].std;
    ncw_put_var_double(ncid, varid_std, std);

    free(std);

    ncw_close(ncid);
}

/**
 */
void das_standardise(dasystem* das)
{
    observations* obs = das->obs;
    double mult = sqrt((double) das->nmem - 1);
    int e, i;

    if (das->s_mode == S_MODE_S_f || das->s_mode == S_MODE_S_a)
        return;

    if (obs->nobs == 0)
        goto finish;

#if defined (HE_VIASHMEM)
    if (das->sm_rank == 0)
#endif
        for (e = 0; e < das->nmem; ++e) {
            ENSOBSTYPE* Se = das->S[e];

            for (i = 0; i < obs->nobs; ++i) {
                observation* o = &obs->data[i];

                Se[i] /= o->std * sqrt(obs->obstypes[o->type].rfactor) * mult;
            }
        }
#if defined (HE_VIASHMEM)
    MPI_Barrier(das->sm_comm);
#endif
    if (das->s_f != NULL) {
        for (i = 0; i < obs->nobs; ++i) {
            observation* o = &obs->data[i];

            das->s_f[i] /= o->std * sqrt(obs->obstypes[o->type].rfactor) * mult;
        }
    }
    if (das->s_a != NULL) {
        for (i = 0; i < obs->nobs; ++i) {
            observation* o = &obs->data[i];

            das->s_a[i] /= o->std * sqrt(obs->obstypes[o->type].rfactor) * mult;
        }
    }

  finish:
    if (das->s_mode == S_MODE_HA_f)
        das->s_mode = S_MODE_S_f;
    else if (das->s_mode == S_MODE_HA_a)
        das->s_mode = S_MODE_S_a;
    else
        enkf_quit("programming error");

}

/**
 */
static int cmp_obs_byij(const void* p1, const void* p2, void* p)
{
    observation* o1 = (observation*) p1;
    observation* o2 = (observation*) p2;
    dasystem* das = (dasystem*) p;
    obstype* obstypes = das->obs->obstypes;
    int gid1 = obstypes[o1->type].gridid;
    int gid2 = obstypes[o2->type].gridid;
    int i1, i2;

    if (gid1 > gid2)
        return 1;
    if (gid1 < gid2)
        return -1;

    i1 = (int) floor(o1->fj + 0.5);
    i2 = (int) floor(o2->fj + 0.5);
    if (i1 > i2)
        return 1;
    if (i1 < i2)
        return -1;

    i1 = (int) floor(o1->fi + 0.5);
    i2 = (int) floor(o2->fi + 0.5);
    if (i1 > i2)
        return 1;
    if (i1 < i2)
        return -1;

    return 0;
}

/**
 */
void das_destandardise(dasystem* das)
{
    observations* obs = das->obs;
    double mult = sqrt((double) das->nmem - 1);
    int e, i;

    if (das->s_mode == S_MODE_HA_f || das->s_mode == S_MODE_HA_a)
        return;

    if (obs->nobs == 0)
        goto finish;

#if defined(HE_VIASHMEM)
    if (das->sm_rank == 0)
#endif
        for (e = 0; e < das->nmem; ++e) {
            ENSOBSTYPE* Se = das->S[e];

            for (i = 0; i < obs->nobs; ++i) {
                observation* o = &obs->data[i];

                Se[i] *= o->std * sqrt(obs->obstypes[o->type].rfactor) * mult;
            }
        }
#if defined (HE_VIASHMEM)
    MPI_Barrier(das->sm_comm);
#endif
    if (das->s_f != NULL) {
        for (i = 0; i < obs->nobs; ++i) {
            observation* o = &obs->data[i];

            das->s_f[i] *= o->std * sqrt(obs->obstypes[o->type].rfactor) * mult;
        }
    }
    if (das->s_a != NULL) {
        for (i = 0; i < obs->nobs; ++i) {
            observation* o = &obs->data[i];

            das->s_a[i] *= o->std * sqrt(obs->obstypes[o->type].rfactor) * mult;
        }
    }

  finish:
    if (das->s_mode == S_MODE_S_f)
        das->s_mode = S_MODE_HA_f;
    else if (das->s_mode == S_MODE_S_a)
        das->s_mode = S_MODE_HA_a;
    else
        enkf_quit("programming error");

}

/** Sorts all observations by j,i so that they can be processed in a single
 * cycle over horizontal grid cells.
 */
static void das_sortobs_byij(dasystem* das)
{
    observations* obs = das->obs;
    int o, e;

    if (das->sort_mode == OBS_SORTMODE_IJ)
        return;
    if (obs->nobs == 0)
        return;

    enkf_printf("    sorting obs by ij:\n");
    qsort_r(obs->data, obs->nobs, sizeof(observation), cmp_obs_byij, das);

    if (das->s_f != NULL) {
        double* s = calloc(obs->nobs, sizeof(double));

        for (o = 0; o < obs->nobs; ++o)
            s[o] = das->s_f[obs->data[o].id];
        memcpy(das->s_f, s, obs->nobs * sizeof(double));

        for (o = 0; o < obs->nobs; ++o)
            s[o] = das->std_f[obs->data[o].id];
        memcpy(das->std_f, s, obs->nobs * sizeof(double));

        free(s);
    }

    if (das->s_a != NULL) {
        double* s = calloc(obs->nobs, sizeof(double));

        for (o = 0; o < obs->nobs; ++o)
            s[o] = das->s_a[obs->data[o].id];
        memcpy(das->s_a, s, obs->nobs * sizeof(double));

        for (o = 0; o < obs->nobs; ++o)
            s[o] = das->std_a[obs->data[o].id];
        memcpy(das->std_a, s, obs->nobs * sizeof(double));

        free(s);
    }
#if defined(HE_VIASHMEM)
    if (das->sm_rank == 0)
#endif
    {
        ENSOBSTYPE* S = calloc(obs->nobs, sizeof(ENSOBSTYPE));

        for (e = 0; e < das->nmem; ++e) {
            ENSOBSTYPE* Se = das->S[e];

            for (o = 0; o < obs->nobs; ++o)
                S[o] = Se[obs->data[o].id];
            memcpy(Se, S, obs->nobs * sizeof(ENSOBSTYPE));
        }
        free(S);
    }
#if defined(HE_VIASHMEM)
    MPI_Barrier(das->sm_comm);
#endif

    das->sort_mode = OBS_SORTMODE_IJ;
}

/**
 */
static void das_changeSmode(dasystem* das, int mode_from, int mode_to)
{
    assert(das->s_mode == mode_from);

    if (das->obs->nobs == 0)
        goto finish;

    if (mode_from == S_MODE_HA_f && mode_to == S_MODE_HE_f) {
        observations* obs = das->obs;
        int e, o;

#if defined(HE_VIASHMEM)
        if (das->sm_rank == 0)
#endif
            for (e = 0; e < das->nmem; ++e) {
                ENSOBSTYPE* Se = das->S[e];

                for (o = 0; o < obs->nobs; ++o)
                    /*
                     * das->s_f is innovation = obs - forecast; hence
                     * forecast = obs - innovation 
                     */
                    Se[o] += obs->data[o].value - das->s_f[o];
            }
#if defined(HE_VIASHMEM)
        MPI_Barrier(das->sm_comm);
#endif
    } else
        enkf_quit("das_changesmode(): transition from mode %d to mode %d is not handled yet\n", mode_from, mode_to);

  finish:
    das->s_mode = mode_to;
}

/** Rolls back sorting of ensemble observations S by i,j, so that S is again
 * consistent with other variables.
 */
static void das_sortobs_byid(dasystem* das)
{
    observations* obs = das->obs;
    int o, e;

    if (das->sort_mode == OBS_SORTMODE_ID)
        return;

    if (das->s_f != NULL) {
        double* s = calloc(obs->nobs, sizeof(double));

        for (o = 0; o < obs->nobs; ++o)
            s[obs->data[o].id] = das->s_f[o];
        memcpy(das->s_f, s, obs->nobs * sizeof(double));

        for (o = 0; o < obs->nobs; ++o)
            s[obs->data[o].id] = das->std_f[o];
        memcpy(das->std_f, s, obs->nobs * sizeof(double));

        free(s);
    }

    if (das->s_a != NULL) {
        double* s = calloc(obs->nobs, sizeof(double));

        for (o = 0; o < obs->nobs; ++o)
            s[obs->data[o].id] = das->s_a[o];
        memcpy(das->s_a, s, obs->nobs * sizeof(double));

        for (o = 0; o < obs->nobs; ++o)
            s[obs->data[o].id] = das->std_a[o];
        memcpy(das->std_a, s, obs->nobs * sizeof(double));

        free(s);
    }
#if defined(HE_VIASHMEM)
    if (das->sm_rank == 0)
#endif
    {
        ENSOBSTYPE* S = calloc(obs->nobs, sizeof(ENSOBSTYPE));

        for (e = 0; e < das->nmem; ++e) {
            ENSOBSTYPE* Se = das->S[e];

            for (o = 0; o < obs->nobs; ++o)
                S[obs->data[o].id] = Se[o];
            memcpy(Se, S, obs->nobs * sizeof(ENSOBSTYPE));
        }
        free(S);
    }
#if defined(HE_VIASHMEM)
    MPI_Barrier(das->sm_comm);
#endif

    /*
     * order obs back by id
     */
    obs_inorder(obs);

    das->sort_mode = OBS_SORTMODE_ID;
}

/**
 */
#if defined(HE_VIASHMEM)
static void gather_St(dasystem* das)
{
    int nmem = das->nmem;
    MPI_Datatype mpitype_vec_nmem;
    int ierror;
    int* recvcounts = calloc(nprocesses, sizeof(int));
    int* displs = calloc(nprocesses, sizeof(int));
    int i, ii;

    ierror = MPI_Type_contiguous(nmem, MPIENSOBSTYPE, &mpitype_vec_nmem);
    assert(ierror == MPI_SUCCESS);

    MPI_Barrier(MPI_COMM_WORLD);
    if (das->node_rank >= 0 && das->node_size > 1) {
        for (i = 0, ii = -1; i < nprocesses; ++i) {
            if (das->node_ranks[i] >= 0) {
                ii = das->node_ranks[i];
                displs[ii] = first_iteration[i];
            }
            assert(ii >= 0);
            recvcounts[ii] += number_of_iterations[i];
        }
        ierror = MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, das->St[0], recvcounts, displs, mpitype_vec_nmem, das->node_comm);
        assert(ierror == MPI_SUCCESS);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    free(displs);
    free(recvcounts);
}
#endif

/** Updates ensemble observations by applying X5.
 *
 * Both update_HE() and update_Hx() use transforms calculated for (int) fj,
 * (int) fi node of the grid. This might be improved in future to yield
 * probably a bit more precise (and better) analysis stats.
 */
static void update_HE(dasystem* das)
{
    model* m = das->m;
    int ngrid = model_getngrid(m);
    int nmem = das->nmem;
    int gid;
    observations* obs = das->obs;
    int nobs = obs->nobs;
    int e, o;
    ENSOBSTYPE* HEi_f;
    ENSOBSTYPE* HEi_a;
    char do_T = 'T';
    float alpha = 1.0f;
    float beta = 0.0f;
    int inc = 1;

    enkf_printf("    updating HE:\n");
    assert(das->s_mode == S_MODE_HE_f);

#if defined(HE_VIASHMEM)
    distribute_iterations(0, nobs - 1, (nprocesses > HE_NPROCMAX) ? HE_NPROCMAX : nprocesses, rank, "    ");
    for (e = 0; e < nmem; ++e)
        for (o = my_first_iteration; o <= my_last_iteration; ++o)
            das->St[o][e] = das->S[e][o];
#else
    my_first_iteration = 0;
    my_last_iteration = nobs - 1;
#endif

#if !defined(HE_VIASHMEM)
    HEi_f = malloc(nmem * sizeof(ENSOBSTYPE));
#endif
    HEi_a = malloc(nmem * sizeof(ENSOBSTYPE));

    /*
     * the following code for interpolation of X5 essentially coincides with
     * that in das_updatefields() 
     */

    for (gid = 0, o = my_first_iteration; gid < ngrid && o <= my_last_iteration; ++gid) {
        void* grid = model_getgridbyid(m, gid);
        int periodic_i = grid_isperiodic_i(grid);
        int stride = grid_getstride(grid);

        char fname_X5[MAXSTRLEN];
        int ncid;
        int varid;
        int dimids[3];
        size_t dimlens[3];
        size_t start[3], count[3];
        float** X5j = NULL;
        float** X5jj = NULL;
        float** X5jj1 = NULL;
        float** X5jj2 = NULL;

        int mni, mnj;
        int* iiter;
        int* jiter;
        int i, j, ni, nj;
        int jj, stepj, ii, stepi;

        if (gid < obs->obstypes[obs->data[o].type].gridid)
            continue;

        das_getfname_X5(das, grid, fname_X5);

        ncw_open(fname_X5, NC_NOWRITE, &ncid);
        ncw_inq_varid(ncid, "X5", &varid);
        ncw_inq_vardimid(ncid, varid, dimids);
        for (i = 0; i < 3; ++i)
            ncw_inq_dimlen(ncid, dimids[i], &dimlens[i]);
        ni = dimlens[1];
        nj = dimlens[0];

        assert((int) dimlens[2] == nmem * nmem);

        jiter = malloc((nj + 1) * sizeof(int)); /* "+ 1" to handle periodic
                                                 * grids */
        iiter = malloc((ni + 1) * sizeof(int));
        for (j = 0, i = 0; j < nj; ++j, i += stride)
            jiter[j] = i;
        for (i = 0, j = 0; i < ni; ++i, j += stride)
            iiter[i] = j;
        if (periodic_i)
            iiter[ni] = iiter[ni - 1] + stride;

        grid_getdims(grid, &mni, &mnj, NULL);

        start[0] = 0;
        start[1] = 0;
        start[2] = 0;
        count[0] = 1;
        count[1] = ni;
        count[2] = nmem * nmem;
        X5j = alloc2d(mni, nmem * nmem, sizeof(float));
        if (stride > 1) {
            X5jj = alloc2d(ni, nmem * nmem, sizeof(float));
            X5jj1 = alloc2d(ni, nmem * nmem, sizeof(float));
            X5jj2 = alloc2d(ni, nmem * nmem, sizeof(float));
            ncw_get_vara_float(ncid, varid, start, count, X5jj2[0]);
        }

        /*
         * jj, ii are the indices of the subsampled grid; i, j are the indices
         * of the model grid 
         */
        for (jj = 0, j = 0; jj < nj; ++jj) {
            for (stepj = 0; stepj < stride && j < mnj; ++stepj, ++j) {

                if ((int) (obs->data[o].fj + 0.5) / stride > jj + 1)
                    continue;

                if (stride == 1) {
                    /*
                     * no interpolation necessary; simply read the ETMs for the
                     * j-th row from disk 
                     */
                    start[0] = j;
                    ncw_get_vara_float(ncid, varid, start, count, X5j[0]);
                } else {
                    /*
                     * the following code interpolates the ETM back to the
                     * original grid, first by j, and then by i 
                     */
                    if (stepj == 0) {
                        memcpy(X5jj[0], X5jj2[0], ni * nmem * nmem * sizeof(float));
                        memcpy(X5jj1[0], X5jj2[0], ni * nmem * nmem * sizeof(float));
                        if (jj < nj - 1) {
                            start[0] = (jj + 1) % nj;
                            ncw_get_vara_float(ncid, varid, start, count, X5jj2[0]);
                        }
                    } else {
                        float weight2 = (float) stepj / stride;
                        float weight1 = (float) 1.0 - weight2;

                        for (ii = 0; ii < ni; ++ii) {
                            float* X5jjii = X5jj[ii];
                            float* X5jj1ii = X5jj1[ii];
                            float* X5jj2ii = X5jj2[ii];

                            for (e = 0; e < nmem * nmem; ++e)
                                X5jjii[e] = X5jj1ii[e] * weight1 + X5jj2ii[e] * weight2;
                        }
                    }

                    for (ii = 0, i = 0; ii < ni; ++ii) {
                        for (stepi = 0; stepi < stride && i < mni; ++stepi, ++i) {
                            if (stepi == 0)
                                memcpy(X5j[i], X5jj[ii], nmem * nmem * sizeof(float));
                            else {
                                float weight2 = (float) stepi / stride;
                                float weight1 = (float) 1.0 - weight2;
                                float* X5jjii1 = X5jj[ii];
                                float* X5ji = X5j[i];
                                float* X5jjii2;

                                if (ii < ni - 1)
                                    X5jjii2 = X5jj[ii + 1];
                                else
                                    X5jjii2 = X5jj[(periodic_i) ? (ii + 1) % ni : ii];

                                for (e = 0; e < nmem * nmem; ++e)
                                    X5ji[e] = X5jjii1[e] * weight1 + X5jjii2[e] * weight2;
                            }
                        }
                    }
                }               /* stride != 1 */

                /*
                 * (at this stage X5j should contain the array of X5 matrices
                 * for the j-th row of the grid) 
                 */

                if (o > my_last_iteration)
                    break;

                for (; o <= my_last_iteration && (int) (obs->data[o].fj + 0.5) == j; ++o) {
                    float inflation0 = NAN;
                    double inf_ratio = NAN;
                    float inflation = NAN;
                    double v1_a = 0.0;

                    model_getvarinflation(m, obs->obstypes[obs->data[o].type].vid, &inflation0, &inf_ratio);

                    /*
                     * HE(i, :) = HE(i, :) * X5 
                     */
                    i = (int) (obs->data[o].fi + 0.5);
                    if (i == mni)
                        i--;
#if defined(HE_VIASHMEM)
                    HEi_f = das->St[o];
#else
                    for (e = 0; e < nmem; ++e)
                        HEi_f[e] = das->S[e][o];
#endif
                    sgemv_(&do_T, &nmem, &nmem, &alpha, X5j[i], &nmem, HEi_f, &inc, &beta, HEi_a, &inc);

                    for (e = 0; e < nmem; ++e)
                        v1_a += HEi_a[e];
                    v1_a /= (double) nmem;

                    if (!isnan(inf_ratio)) {
                        double v1_f = 0.0;
                        double v2_f = 0.0;
                        double v2_a = 0.0;
                        double var_a, var_f;

                        for (e = 0; e < nmem; ++e) {
                            double ve = (double) HEi_f[e];

                            v1_f += ve;
                            v2_f += ve * ve;
                        }
                        v1_f /= (double) nmem;
                        var_f = v2_f / (double) nmem - v1_f * v1_f;

                        for (e = 0; e < nmem; ++e) {
                            double ve = (double) HEi_a[e];

                            v2_a += ve * ve;
                        }
                        var_a = v2_a / (double) nmem - v1_a * v1_a;

                        if (var_a > 0) {
                            /*
                             * (Normal case.) Limit inflation by inf_ratio of
                             * the magnitude of spread reduction.
                             */
                            inflation = (float) (sqrt(var_f / var_a) * inf_ratio + 1.0 - inf_ratio);
                            if (inflation >= inflation0)
                                inflation = inflation0;
                        }
                    } else
                        inflation = inflation0;

                    /*
                     * applying inflation:
                     */
                    if (fabsf(inflation - 1.0f) > EPSF)
                        for (e = 0; e < nmem; ++e)
                            HEi_a[e] = (HEi_a[e] - (float) v1_a) * inflation + v1_a;

#if defined(HE_VIASHMEM)
                    for (e = 0; e < nmem; ++e)
                        das->St[o][e] = HEi_a[e];
#else
                    for (e = 0; e < nmem; ++e)
                        das->S[e][o] = HEi_a[e];
#endif
                }
            }                   /* for stepj */
        }                       /* for jj */

        ncw_close(ncid);

        free(iiter);
        free(jiter);
        free(X5j);
        if (stride > 1) {
            free(X5jj);
            free(X5jj1);
            free(X5jj2);
        }
    }                           /* for gid */

#if defined(HE_VIASHMEM)
    gather_St(das);

    if (rank == 0)
        for (e = 0; e < nmem; ++e)
            for (o = 0; o < nobs; ++o)
                das->S[e][o] = das->St[o][e];
#endif

    free(HEi_a);
#if !defined(HE_VIASHMEM)
    free(HEi_f);
#endif
    das->s_mode = S_MODE_HE_a;
}                               /* update_HE() */

/**
 */
static void update_Hx(dasystem* das)
{
    model* m = das->m;
    int ngrid = model_getngrid(m);
    int nmem = das->nmem;
    int gid;
    observations* obs = das->obs;
    int nobs = obs->nobs;
    int e, o;

    enkf_printf("    updating Hx:\n");
    assert(das->s_mode == S_MODE_HE_f);

#if defined(HE_VIASHMEM)
    distribute_iterations(0, nobs - 1, (nprocesses > HE_NPROCMAX) ? HE_NPROCMAX : nprocesses, rank, "    ");
    for (e = 0; e < nmem; ++e)
        for (o = my_first_iteration; o <= my_last_iteration; ++o)
            das->St[o][e] = das->S[e][o];
#else
    my_first_iteration = 0;
    my_last_iteration = nobs - 1;
#endif

    /*
     * the following code for interpolation of X5 essentially coincides with
     * that in das_updatefields() 
     */

    for (gid = 0, o = my_first_iteration; gid < ngrid && o <= my_last_iteration; ++gid) {
        void* grid = model_getgridbyid(m, gid);
        int periodic_i = grid_isperiodic_i(grid);
        int stride = grid_getstride(grid);

        char fname_w[MAXSTRLEN];
        int ncid;
        int varid;
        int dimids[3];
        size_t dimlens[3];
        size_t start[3], count[3];
        float** wj = NULL;
        float** wjj = NULL;
        float** wjj1 = NULL;
        float** wjj2 = NULL;

        int mni, mnj;
        int* iiter;
        int* jiter;
        int i, j, ni, nj;
        int jj, stepj, ii, stepi;

        if (gid < obs->obstypes[obs->data[o].type].gridid)
            continue;

        das_getfname_w(das, grid, fname_w);

        ncw_open(fname_w, NC_NOWRITE, &ncid);
        ncw_inq_varid(ncid, "w", &varid);
        ncw_inq_vardimid(ncid, varid, dimids);
        for (i = 0; i < 3; ++i)
            ncw_inq_dimlen(ncid, dimids[i], &dimlens[i]);
        ni = dimlens[1];
        nj = dimlens[0];

        assert((int) dimlens[2] == nmem);

        jiter = malloc((nj + 1) * sizeof(int)); /* "+ 1" to handle periodic
                                                 * grids */
        iiter = malloc((ni + 1) * sizeof(int));
        for (j = 0, i = 0; j < nj; ++j, i += stride)
            jiter[j] = i;
        for (i = 0, j = 0; i < ni; ++i, j += stride)
            iiter[i] = j;
        if (periodic_i)
            iiter[ni] = iiter[ni - 1] + stride;

        grid_getdims(grid, &mni, &mnj, NULL);

        start[0] = 0;
        start[1] = 0;
        start[2] = 0;
        count[0] = 1;
        count[1] = ni;
        count[2] = nmem;
        wj = alloc2d(mni, nmem, sizeof(float));
        if (stride > 1) {
            wjj = alloc2d(ni, nmem, sizeof(float));
            wjj1 = alloc2d(ni, nmem, sizeof(float));
            wjj2 = alloc2d(ni, nmem, sizeof(float));
            ncw_get_vara_float(ncid, varid, start, count, wjj2[0]);
        }

        /*
         * jj, ii are the indices of the subsampled grid; i, j are the indices
         * of the model grid 
         */
        for (jj = 0, j = 0; jj < nj; ++jj) {
            for (stepj = 0; stepj < stride && j < mnj; ++stepj, ++j) {

                if ((int) (obs->data[o].fj + 0.5) / stride > jj + 1)
                    continue;

                if (stride == 1) {
                    /*
                     * no interpolation necessary; simply read the ETMs for the
                     * j-th row from disk 
                     */
                    start[0] = j;
                    ncw_get_vara_float(ncid, varid, start, count, wj[0]);
                } else {
                    /*
                     * the following code interpolates the ETM back to the
                     * original grid, first by j, and then by i 
                     */
                    if (stepj == 0) {
                        memcpy(wjj[0], wjj2[0], ni * nmem * sizeof(float));
                        memcpy(wjj1[0], wjj2[0], ni * nmem * sizeof(float));
                        if (jj < nj - 1) {
                            start[0] = (jj + 1) % nj;
                            ncw_get_vara_float(ncid, varid, start, count, wjj2[0]);
                        }
                    } else {
                        float weight2 = (float) stepj / stride;
                        float weight1 = (float) 1.0 - weight2;

                        for (ii = 0; ii < ni; ++ii) {
                            float* wjjii = wjj[ii];
                            float* wjj1ii = wjj1[ii];
                            float* wjj2ii = wjj2[ii];

                            for (e = 0; e < nmem; ++e)
                                wjjii[e] = wjj1ii[e] * weight1 + wjj2ii[e] * weight2;
                        }
                    }

                    for (ii = 0, i = 0; ii < ni; ++ii) {
                        for (stepi = 0; stepi < stride && i < mni; ++stepi, ++i) {
                            if (stepi == 0)
                                memcpy(wj[i], wjj[ii], nmem * sizeof(float));
                            else {
                                float weight2 = (float) stepi / stride;
                                float weight1 = (float) 1.0 - weight2;
                                float* wjjii1 = wjj[ii];
                                float* wji = wj[i];
                                float* wjjii2;

                                if (ii < ni - 1)
                                    wjjii2 = wjj[ii + 1];
                                else
                                    wjjii2 = wjj[(periodic_i) ? (ii + 1) % ni : ii];

                                for (e = 0; e < nmem; ++e)
                                    wji[e] = wjjii1[e] * weight1 + wjjii2[e] * weight2;
                            }
                        }
                    }
                }               /* stride != 1 */

                /*
                 * (at this stage wj should contain the array of b vectors for
                 * the j-th row of the grid) 
                 */

                if (o > my_last_iteration)
                    break;

                for (; o <= my_last_iteration && (int) (obs->data[o].fj + 0.5) == j; ++o) {
                    double dHx = 0.0;
                    double Hx = 0.0;

#if defined(HE_VIASHMEM)
                    for (e = 0; e < nmem; ++e)
                        Hx += das->St[o][e];
#else
                    for (e = 0; e < nmem; ++e)
                        Hx += das->S[e][o];
#endif
                    Hx /= (double) nmem;

                    i = (int) (obs->data[o].fi + 0.5);
                    if (i == mni)
                        i--;
                    /*
                     * HE(i, :) += HA(i, :) * b * 1' 
                     */
#if defined(HE_VIASHMEM)
                    for (e = 0; e < nmem; ++e)
                        dHx += (das->St[o][e] - Hx) * wj[i][e];
                    for (e = 0; e < nmem; ++e)
                        das->St[o][e] += dHx;
#else
                    for (e = 0; e < nmem; ++e)
                        dHx += (das->S[e][o] - Hx) * wj[i][e];
                    for (e = 0; e < nmem; ++e)
                        das->S[e][o] += dHx;
#endif
                }
            }                   /* for stepj */
        }                       /* for jj */

        ncw_close(ncid);

        free(iiter);
        free(jiter);
        free(wj);
        if (stride > 1) {
            free(wjj);
            free(wjj1);
            free(wjj2);
        }
    }                           /* for gid */

#if defined(HE_VIASHMEM)
    gather_St(das);

    if (rank == 0)
        for (e = 0; e < nmem; ++e)
            for (o = 0; o < nobs; ++o)
                das->S[e][o] = das->St[o][e];
#endif
    das->s_mode = S_MODE_HE_a;
}                               /* update_Hx() */

/**
 */
void das_updateHE(dasystem* das)
{
    das_destandardise(das);
    das_sortobs_byij(das);
    das_changeSmode(das, S_MODE_HA_f, S_MODE_HE_f);
    enkf_printtime("    ");
    if (das->mode == MODE_ENKF)
        update_HE(das);
    else
        update_Hx(das);
    enkf_printtime("    ");
    das_calcinnandspread(das);
    das_sortobs_byid(das);
}

/** Add analysed observation estimates and ensemble spread to FNAME_SOBS.
 */
void das_addanalysis(dasystem* das, char fname[])
{
    int ncid;
    int dimid_nobs[1];
    size_t nobs;
    int varid_Hx, varid_spread;
    observation* data = das->obs->data;
    double* s;
    int i;

    if (rank != 0)
        return;

    assert(das->s_mode == S_MODE_HA_a);

    ncw_open(fname, NC_WRITE, &ncid);
    ncw_inq_dimid(ncid, "nobs", dimid_nobs);
    ncw_inq_dimlen(ncid, dimid_nobs[0], &nobs);
    assert(nobs == das->obs->nobs);

    if (ncw_var_exists(ncid, "Hx_a")) {
        enkf_printf("  Hx_a already added to \"%s\", overwriting\n", fname);
        ncw_inq_varid(ncid, "Hx_a", &varid_Hx);
        ncw_inq_varid(ncid, "std_a", &varid_spread);
    } else {
        ncw_redef(ncid);
        ncw_def_var(ncid, "Hx_a", NC_FLOAT, 1, dimid_nobs, &varid_Hx);
        ncw_def_var(ncid, "std_a", NC_FLOAT, 1, dimid_nobs, &varid_spread);
        ncw_enddef(ncid);
    }

    ncw_put_var_double(ncid, varid_spread, das->std_a);
    s = malloc(nobs * sizeof(double));
    /*
     * the obs are still sorted by ij 
     */
    for (i = 0; i < (int) nobs; ++i)
        s[data[i].id] = data[i].value;
    for (i = 0; i < (int) nobs; ++i)
        s[i] -= das->s_a[i];
    ncw_put_var_double(ncid, varid_Hx, s);
    free(s);

    ncw_close(ncid);
}
