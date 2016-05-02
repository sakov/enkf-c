/******************************************************************************
 *
 * File:        ensobs.c        
 *
 * Created:     11/2013
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description: The code in this file structurally belongs to dasystem.c,
 *              and is put in a separate file just to break dasystem.c in
 *              smaller parts.
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
#include "nan.h"
#include "distribute.h"
#include "lapack.h"
#include "definitions.h"
#include "utils.h"
#include "grid.h"
#include "model2obs.h"
#include "allhs.h"
#include "dasystem.h"

#define EPSF 1.0e-6f

/**
 */
void das_getHE(dasystem* das)
{
    observations* obs = das->obs;
    model* m = das->m;
    ENSOBSTYPE* Hx = NULL;
    int i, e;

    das->s_mode = S_MODE_HE_f;
    if (obs->nobs == 0)
        return;

    if (das->nmem <= 0)
        das_getnmem(das);
    enkf_printf("    ensemble size = %d\n", das->nmem);
    assert(das->nmem > 0);

    distribute_iterations(0, das->nmem - 1, nprocesses, rank, "    ");

    /*
     * ensemble observation array to be filled 
     */
    assert(das->S == NULL);
    das->S = alloc2d(das->nmem, obs->nobs, sizeof(ENSOBSTYPE));
    if (das->mode == MODE_ENOI)
        Hx = calloc(obs->nobs, sizeof(ENSOBSTYPE));

    for (i = 0; i < obs->nobstypes; ++i) {
        obstype* ot = &obs->obstypes[i];
        float*** vvv = NULL;
        float** vv = NULL;
        H_fn H = NULL;
        int mvid;
        int ni, nj, nk;
        int nobs;
        int* obsids;
        char fname[MAXSTRLEN];

        enkf_printf("    %s ", ot->name);
        fflush(stdout);

        mvid = model_getvarid(m, obs->obstypes[i].varname, 1);
        if (ot->issurface) {
            model_getvardims(m, mvid, &ni, &nj, NULL);
            vv = alloc2d(nj, ni, sizeof(float));
        } else {
            model_getvardims(m, mvid, &ni, &nj, &nk);
            vvv = alloc3d(nk, nj, ni, sizeof(float));
        }

        /*
         * set H
         */
        H = getH(ot->name, ot->hfunction);

        if (ot->isasync) {
            int t1 = get_tshift(ot->date_min, ot->async_tstep);
            int t2 = get_tshift(ot->date_max, ot->async_tstep);
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
                        int success = model_getbgfname_async(m, das->bgdir, ot->varname, ot->name, t, fname);

                        H(das, nobs, obsids, fname, -1, t, ot->varname, ot->varname2, (ot->issurface) ? (void*) vv : (void*) vvv, Hx);
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
                        int success = model_getmemberfname_async(m, das->ensdir, ot->varname, ot->name, e + 1, t, fname);

                        H(das, nobs, obsids, fname, e + 1, t, ot->varname, ot->varname2, (ot->issurface) ? (void*) vv : (void*) vvv, das->S[e]);
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
                    model_getbgfname(m, das->bgdir, ot->varname, fname);
                    H(das, nobs, obsids, fname, -1, MAXINT, ot->varname, ot->varname2, (ot->issurface) ? (void*) vv : (void*) vvv, Hx);
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
                    model_getmemberfname(m, das->ensdir, ot->varname, e + 1, fname);
                    H(das, nobs, obsids, fname, e + 1, MAXINT, ot->varname, ot->varname2, (ot->issurface) ? (void*) vv : (void*) vvv, das->S[e]);
                    enkf_printf(".");
                    fflush(stdout);
                }
            }

            free(obsids);
        }

      next:

        if (ot->issurface)
            free2d(vv);
        else
            free3d(vvv);
        enkf_printf("\n");
    }                           /* for i (over obstypes) */

#if defined(MPI)
    if (das->mode == MODE_ENKF || !enkf_fstatsonly) {
#if !defined(HE_VIAFILE)
        /*
         * communicate HE via MPI
         */
        int ierror, sendcount, *recvcounts, *displs;

        recvcounts = malloc(nprocesses * sizeof(int));
        displs = malloc(nprocesses * sizeof(int));

        sendcount = my_number_of_iterations * obs->nobs;
        for (i = 0; i < nprocesses; ++i) {
            recvcounts[i] = number_of_iterations[i] * obs->nobs;
            displs[i] = first_iteration[i] * obs->nobs;
        }

        ierror = MPI_Allgatherv(das->S[my_first_iteration], sendcount, MPIENSOBSTYPE, das->S[0], recvcounts, displs, MPIENSOBSTYPE, MPI_COMM_WORLD);
        assert(ierror == MPI_SUCCESS);

        free(recvcounts);
        free(displs);
#else
        /*
         * communicate HE via file
         */
        {
            int ncid;
            int varid;
            size_t start[2], count[2];

            if (rank == 0) {
                int dimids[2];

                ncw_create(FNAME_HE, NC_CLOBBER | NETCDF_FORMAT, &ncid);
                ncw_def_dim(FNAME_HE, ncid, "m", das->nmem, &dimids[0]);
                ncw_def_dim(FNAME_HE, ncid, "p", obs->nobs, &dimids[1]);
                ncw_def_var(FNAME_HE, ncid, "HE", NC_FLOAT, 2, dimids, &varid);
                ncw_close(FNAME_HE, ncid);
            }
            MPI_Barrier(MPI_COMM_WORLD);

            ncw_open(FNAME_HE, NC_WRITE, &ncid);
            ncw_inq_varid(FNAME_HE, ncid, "HE", &varid);
            start[0] = my_first_iteration;
            start[1] = 0;
            count[0] = my_last_iteration - my_first_iteration + 1;
            count[1] = obs->nobs;
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
            double* ensmean = calloc(obs->nobs, sizeof(double));

            for (e = 0; e < das->nmem; ++e) {
                ENSOBSTYPE* Se = das->S[e];

                for (i = 0; i < obs->nobs; ++i)
                    ensmean[i] += Se[i];
            }
            for (i = 0; i < obs->nobs; ++i)
                ensmean[i] /= (double) das->nmem;

            for (e = 0; e < das->nmem; ++e) {
                ENSOBSTYPE* Se = das->S[e];

                for (i = 0; i < obs->nobs; ++i)
                    Se[i] += Hx[i] - ensmean[i];
            }

            free(ensmean);
        } else {
            for (e = 0; e < das->nmem; ++e) {
                ENSOBSTYPE* Se = das->S[e];

                for (i = 0; i < obs->nobs; ++i)
                    Se[i] = Hx[i];
            }
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
    int e, o;

    if (obs->nobs == 0)
        goto finish;

    if (das->s_mode == S_MODE_HE_f) {
        if (das->s_f == NULL) {
            das->s_f = calloc(obs->nobs, sizeof(double));
            assert(das->std_f == NULL);
            das->std_f = calloc(obs->nobs, sizeof(double));
        } else {
            memset(das->s_f, 0, obs->nobs * sizeof(double));
            memset(das->std_f, 0, obs->nobs * sizeof(double));
        }

        /*
         * calculate ensemble mean observations 
         */
        for (e = 0; e < das->nmem; ++e) {
            ENSOBSTYPE* Se = das->S[e];

            for (o = 0; o < obs->nobs; ++o)
                das->s_f[o] += (double) Se[o];
        }
        for (o = 0; o < obs->nobs; ++o)
            das->s_f[o] /= (double) das->nmem;

        /*
         * calculate ensemble spread and innovation 
         */
        for (e = 0; e < das->nmem; ++e) {
            ENSOBSTYPE* Se = das->S[e];

            for (o = 0; o < obs->nobs; ++o) {
                Se[o] -= (ENSOBSTYPE) das->s_f[o];
                das->std_f[o] += (double) (Se[o] * Se[o]);
            }
        }
        for (o = 0; o < obs->nobs; ++o) {
            observation* m = &obs->data[o];

            if (m->status != STATUS_OK)
                continue;
            das->std_f[o] = sqrt(das->std_f[o] / (double) (das->nmem - 1));
            das->s_f[o] = m->value - das->s_f[o];
            if (!isfinite(das->s_f[o]) || fabs(das->s_f[o]) > STATE_BIGNUM)
                enkf_quit("obs # %d: Hx = %d, no point to continue", o);
        }

        das->s_mode = S_MODE_HA_f;
    } else if (das->s_mode == S_MODE_HE_a) {
        if (das->s_a == NULL) {
            das->s_a = calloc(obs->nobs, sizeof(double));
            assert(das->std_a == NULL);
            das->std_a = calloc(obs->nobs, sizeof(double));
        } else {
            memset(das->s_a, 0, obs->nobs * sizeof(double));
            memset(das->std_a, 0, obs->nobs * sizeof(double));
        }

        /*
         * calculate ensemble mean observations 
         */
        for (e = 0; e < das->nmem; ++e) {
            ENSOBSTYPE* Se = das->S[e];

            for (o = 0; o < obs->nobs; ++o)
                das->s_a[o] += (double) Se[o];
        }
        for (o = 0; o < obs->nobs; ++o)
            das->s_a[o] /= (double) das->nmem;

        /*
         * calculate ensemble spread and innovation 
         */
        for (e = 0; e < das->nmem; ++e) {
            ENSOBSTYPE* Se = das->S[e];

            for (o = 0; o < obs->nobs; ++o) {
                Se[o] -= (ENSOBSTYPE) das->s_a[o];
                das->std_a[o] += (double) (Se[o] * Se[o]);
            }
        }
        for (o = 0; o < obs->nobs; ++o) {
            observation* m = &obs->data[o];

            if (m->status != STATUS_OK)
                continue;
            das->std_a[o] = sqrt(das->std_a[o] / (double) (das->nmem - 1));
            das->s_a[o] = m->value - das->s_a[o];
            if (!isfinite(das->s_a[o]) || fabs(das->s_a[o]) > STATE_BIGNUM)
                enkf_quit("obs # %d: Hx = %d, no point to continue", o);
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
        enkf_printf("  Hx already added to \"%s\" (skipping)\n", fname);
        goto finish;
    }

    ncw_inq_dimid(fname, ncid, "nobs", dimid_nobs);
    ncw_inq_dimlen(fname, ncid, dimid_nobs[0], &nobs);
    ncw_redef(fname, ncid);
    ncw_def_var(fname, ncid, "Hx_f", NC_FLOAT, 1, dimid_nobs, &varid_Hx);
    ncw_def_var(fname, ncid, "std_f", NC_FLOAT, 1, dimid_nobs, &varid_spread);
    ncw_enddef(fname, ncid);

    Hx = calloc(nobs, sizeof(double));
    for (o = 0; o < (int) nobs; ++o)
        Hx[o] = das->obs->data[o].value - das->s_f[o];

    ncw_put_var_double(fname, ncid, varid_Hx, Hx);
    ncw_put_var_double(fname, ncid, varid_spread, das->std_f);

    free(Hx);

  finish:
    ncw_close(fname, ncid);
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
    double da_julday = NaN;

    if (rank != 0)
        return;

    ncw_open(fname, NC_WRITE, &ncid);
    ncw_inq_dimid(fname, ncid, "nobs", dimid_nobs);
    ncw_inq_dimlen(fname, ncid, dimid_nobs[0], &nobs);

    ncw_get_att_double(fname, ncid, NC_GLOBAL, "DA_JULDAY", &da_julday);
    if (!enkf_noobsdatecheck && (isnan(da_julday) || fabs(das->obs->da_date - da_julday) > 1e-6))
        enkf_quit("\"observations.nc\" from a different cycle");

    if (ncw_var_exists(ncid, "std_orig")) {
        enkf_printf("    nothing to do\n");
        ncw_close(fname, ncid);
        return;
    }

    ncw_redef(fname, ncid);
    ncw_rename_var(fname, ncid, "std", "std_orig");
    ncw_def_var(fname, ncid, "std", NC_FLOAT, 1, dimid_nobs, &varid_std);
    ncw_enddef(fname, ncid);

    std = malloc(nobs * sizeof(double));

    for (i = 0; i < (int) nobs; ++i)
        std[i] = das->obs->data[i].std;
    ncw_put_var_double(fname, ncid, varid_std, std);

    free(std);

    ncw_close(fname, ncid);
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

    for (e = 0; e < das->nmem; ++e) {
        ENSOBSTYPE* Se = das->S[e];

        for (i = 0; i < obs->nobs; ++i) {
            observation* o = &obs->data[i];

            Se[i] /= o->std * sqrt(obs->obstypes[o->type].rfactor) * mult;
        }
    }
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

    i1 = (int) floor(o1->fj);
    i2 = (int) floor(o2->fj);
    if (i1 > i2)
        return 1;
    if (i1 < i2)
        return -1;

    i1 = (int) floor(o1->fi);
    i2 = (int) floor(o2->fi);
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

    for (e = 0; e < das->nmem; ++e) {
        ENSOBSTYPE* Se = das->S[e];

        for (i = 0; i < obs->nobs; ++i) {
            observation* o = &obs->data[i];

            Se[i] *= o->std * sqrt(obs->obstypes[o->type].rfactor) * mult;
        }
    }
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

        for (e = 0; e < das->nmem; ++e) {
            ENSOBSTYPE* Se = das->S[e];

            for (o = 0; o < obs->nobs; ++o)
                /*
                 * das->s_f is innovation = obs - forecast; hence forecast = obs
                 * - innovation 
                 */
                Se[o] += obs->data[o].value - das->s_f[o];
        }
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

    /*
     * order obs back by id
     */
    obs_inorder(obs);

    das->sort_mode = OBS_SORTMODE_ID;
}

/** Updates ensemble observations by applying X5
 */
static void update_HE(dasystem* das)
{
    model* m = das->m;
    int ngrid = model_getngrid(m);
    int nmem = das->nmem;
    int gid;
    observations* obs = das->obs;
    int e, o;
    float* HEi_f;
    float* HEi_a;
    char do_T = 'T';
    float alpha = 1.0f;
    float beta = 0.0f;
    int inc = 1;

    enkf_printf("    updating HE:\n");
    assert(das->s_mode == S_MODE_HE_f);

    HEi_f = malloc(nmem * sizeof(ENSOBSTYPE));
    HEi_a = malloc(nmem * sizeof(ENSOBSTYPE));

    /*
     * the following code for interpolation of X5 essentially coincides with
     * that in das_updatefields() 
     */

    for (gid = 0, o = 0; gid < ngrid && o < obs->nobs; ++gid) {
        void* grid = model_getgridbyid(m, gid);
        int periodic_i = grid_isperiodic_x(grid);
        int periodic_j = grid_isperiodic_y(grid);

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

        assert(obs->obstypes[obs->data[o].type].gridid == gid);

        das_getfname_X5(das, grid, fname_X5);

        ncw_open(fname_X5, NC_NOWRITE, &ncid);
        ncw_inq_varid(fname_X5, ncid, "X5", &varid);
        ncw_inq_vardimid(fname_X5, ncid, varid, dimids);
        for (i = 0; i < 3; ++i)
            ncw_inq_dimlen(fname_X5, ncid, dimids[i], &dimlens[i]);
        ni = dimlens[1];
        nj = dimlens[0];

        assert((int) dimlens[2] == nmem * nmem);

        jiter = malloc((nj + 1) * sizeof(int)); /* "+ 1" to handle periodic
                                                 * grids */
        iiter = malloc((ni + 1) * sizeof(int));
        for (j = 0, i = 0; j < nj; ++j, i += das->stride)
            jiter[j] = i;
        if (periodic_j)
            jiter[nj] = jiter[nj - 1] + das->stride;
        for (i = 0, j = 0; i < ni; ++i, j += das->stride)
            iiter[i] = j;
        if (periodic_i)
            iiter[ni] = iiter[ni - 1] + das->stride;

        grid_getdims(grid, &mni, &mnj, NULL);

        start[0] = 0;
        start[1] = 0;
        start[2] = 0;
        count[0] = 1;
        count[1] = ni;
        count[2] = nmem * nmem;
        X5j = alloc2d(mni, nmem * nmem, sizeof(float));
        if (das->stride > 1) {
            X5jj = alloc2d(ni, nmem * nmem, sizeof(float));
            X5jj1 = alloc2d(ni, nmem * nmem, sizeof(float));
            X5jj2 = alloc2d(ni, nmem * nmem, sizeof(float));
            ncw_get_vara_float(fname_X5, ncid, varid, start, count, X5jj2[0]);
        }

        /*
         * jj, ii are the indices of the subsampled grid; i, j are the indices
         * of the model grid 
         */
        for (jj = 0, j = 0; jj < nj; ++jj) {
            for (stepj = 0; stepj < das->stride && j < mnj; ++stepj, ++j) {
                if (das->stride == 1) {
                    /*
                     * no interpolation necessary; simply read the ETMs for the
                     * j-th row from disk 
                     */
                    start[0] = j;
                    ncw_get_vara_float(fname_X5, ncid, varid, start, count, X5j[0]);
                } else {
                    /*
                     * the following code interpolates the ETM back to the
                     * original grid, first by j, and then by i 
                     */
                    if (stepj == 0) {
                        memcpy(X5jj[0], X5jj2[0], ni * nmem * nmem * sizeof(float));
                        memcpy(X5jj1[0], X5jj2[0], ni * nmem * nmem * sizeof(float));
                        if (jj < nj - 1 || periodic_j) {
                            start[0] = (jj + 1) % nj;
                            ncw_get_vara_float(fname_X5, ncid, varid, start, count, X5jj2[0]);
                        }
                    } else {
                        float weight2 = (float) stepj / das->stride;
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
                        for (stepi = 0; stepi < das->stride && i < mni; ++stepi, ++i) {
                            if (stepi == 0)
                                memcpy(X5j[i], X5jj[ii], nmem * nmem * sizeof(float));
                            else {
                                float weight2 = (float) stepi / das->stride;
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

                if (o >= obs->nobs)
                    break;
                if ((int) (obs->data[o].fj) > j)
                    continue;

                for (; o < obs->nobs && (int) (obs->data[o].fj) == j; ++o) {
                    float inflation0 = NaN;
                    double inf_ratio = NaN;
                    float inflation = NaN;
                    double v1_a = 0.0;

                    model_getvarinflation(m, obs->obstypes[obs->data[o].type].vid, &inflation0, &inf_ratio);

                    /*
                     * HE(i, :) = HE(i, :) * X5 
                     */
                    i = (int) (obs->data[o].fi);
                    for (e = 0; e < nmem; ++e)
                        HEi_f[e] = das->S[e][o];
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

                    for (e = 0; e < nmem; ++e)
                        das->S[e][o] = HEi_a[e];
                }

            }                   /* for stepj */
        }                       /* for jj */

        ncw_close(fname_X5, ncid);

        free(iiter);
        free(jiter);
        free2d(X5j);
        if (das->stride > 1) {
            free2d(X5jj);
            free2d(X5jj1);
            free2d(X5jj2);
        }
    }                           /* for gid */

    free(HEi_a);
    free(HEi_f);
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
    int e, o;

    enkf_printf("    updating Hx:\n");
    assert(das->s_mode == S_MODE_HE_f);

    /*
     * the following code for interpolation of X5 essentially coincides with
     * that in das_updatefields() 
     */

    for (gid = 0, o = 0; gid < ngrid && o < obs->nobs; ++gid) {
        void* grid = model_getgridbyid(m, gid);
        int periodic_i = grid_isperiodic_x(grid);
        int periodic_j = grid_isperiodic_y(grid);

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

        assert(obs->obstypes[obs->data[o].type].gridid == gid);

        das_getfname_w(das, grid, fname_w);

        ncw_open(fname_w, NC_NOWRITE, &ncid);
        ncw_inq_varid(fname_w, ncid, "w", &varid);
        ncw_inq_vardimid(fname_w, ncid, varid, dimids);
        for (i = 0; i < 3; ++i)
            ncw_inq_dimlen(fname_w, ncid, dimids[i], &dimlens[i]);
        ni = dimlens[1];
        nj = dimlens[0];

        assert((int) dimlens[2] == nmem);

        jiter = malloc((nj + 1) * sizeof(int)); /* "+ 1" to handle periodic
                                                 * grids */
        iiter = malloc((ni + 1) * sizeof(int));
        for (j = 0, i = 0; j < nj; ++j, i += das->stride)
            jiter[j] = i;
        if (periodic_j)
            jiter[nj] = jiter[nj - 1] + das->stride;
        for (i = 0, j = 0; i < ni; ++i, j += das->stride)
            iiter[i] = j;
        if (periodic_i)
            iiter[ni] = iiter[ni - 1] + das->stride;

        grid_getdims(grid, &mni, &mnj, NULL);

        start[0] = 0;
        start[1] = 0;
        start[2] = 0;
        count[0] = 1;
        count[1] = ni;
        count[2] = nmem;
        wj = alloc2d(mni, nmem, sizeof(float));
        if (das->stride > 1) {
            wjj = alloc2d(ni, nmem, sizeof(float));
            wjj1 = alloc2d(ni, nmem, sizeof(float));
            wjj2 = alloc2d(ni, nmem, sizeof(float));
            ncw_get_vara_float(fname_w, ncid, varid, start, count, wjj2[0]);
        }

        /*
         * jj, ii are the indices of the subsampled grid; i, j are the indices
         * of the model grid 
         */
        for (jj = 0, j = 0; jj < nj; ++jj) {
            for (stepj = 0; stepj < das->stride && j < mnj; ++stepj, ++j) {
                if (das->stride == 1) {
                    /*
                     * no interpolation necessary; simply read the ETMs for the
                     * j-th row from disk 
                     */
                    start[0] = j;
                    ncw_get_vara_float(fname_w, ncid, varid, start, count, wj[0]);
                } else {
                    /*
                     * the following code interpolates the ETM back to the
                     * original grid, first by j, and then by i 
                     */
                    if (stepj == 0) {
                        memcpy(wjj[0], wjj2[0], ni * nmem * sizeof(float));
                        memcpy(wjj1[0], wjj2[0], ni * nmem * sizeof(float));
                        if (jj < nj - 1 || periodic_j) {
                            start[0] = (jj + 1) % nj;
                            ncw_get_vara_float(fname_w, ncid, varid, start, count, wjj2[0]);
                        }
                    } else {
                        float weight2 = (float) stepj / das->stride;
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
                        for (stepi = 0; stepi < das->stride && i < mni; ++stepi, ++i) {
                            if (stepi == 0)
                                memcpy(wj[i], wjj[ii], nmem * sizeof(float));
                            else {
                                float weight2 = (float) stepi / das->stride;
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

                if (o >= obs->nobs)
                    break;
                if ((int) (obs->data[o].fj) > j)
                    continue;

                for (; o < obs->nobs && (int) (obs->data[o].fj) == j; ++o) {
                    double dHx = 0.0;
                    double Hx = 0.0;

                    for (e = 0; e < nmem; ++e)
                        Hx += das->S[e][o];
                    Hx /= (double) nmem;

                    i = (int) (obs->data[o].fi);
                    /*
                     * HE(i, :) += HA(i, :) * b * 1' 
                     */
                    for (e = 0; e < nmem; ++e)
                        dHx += (das->S[e][o] - Hx) * wj[i][e];
                    for (e = 0; e < nmem; ++e)
                        das->S[e][o] += dHx;
                }
            }                   /* for stepj */
        }                       /* for jj */

        ncw_close(fname_w, ncid);

        free(iiter);
        free(jiter);
        free2d(wj);
        if (das->stride > 1) {
            free2d(wjj);
            free2d(wjj1);
            free2d(wjj2);
        }
    }                           /* for gid */

    das->s_mode = S_MODE_HE_a;
}                               /* update_Hx() */

/**
 */
void das_updateHE(dasystem* das)
{
    if (rank != 0)
        return;

    das_destandardise(das);
    das_sortobs_byij(das);
    das_changeSmode(das, S_MODE_HA_f, S_MODE_HE_f);
    if (das->mode == MODE_ENKF)
        update_HE(das);
    else
        update_Hx(das);
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
    ncw_inq_dimid(fname, ncid, "nobs", dimid_nobs);
    ncw_inq_dimlen(fname, ncid, dimid_nobs[0], &nobs);

    ncw_redef(fname, ncid);
    ncw_def_var(fname, ncid, "Hx_a", NC_FLOAT, 1, dimid_nobs, &varid_Hx);
    ncw_def_var(fname, ncid, "std_a", NC_FLOAT, 1, dimid_nobs, &varid_spread);
    ncw_enddef(fname, ncid);

    ncw_put_var_double(fname, ncid, varid_spread, das->std_a);
    s = malloc(nobs * sizeof(double));
    /*
     * the obs are still sorted by ij 
     */
    for (i = 0; i < (int) nobs; ++i)
        s[data[i].id] = data[i].value;
    for (i = 0; i < (int) nobs; ++i)
        s[i] -= das->s_a[i];
    ncw_put_var_double(fname, ncid, varid_Hx, s);
    free(s);

    ncw_close(fname, ncid);
}
