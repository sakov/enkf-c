/******************************************************************************
 *
 * File:        obsstats.c        
 *
 * Created:     11/2013
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description: Calculates observation statistics for regions.
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
#include "definitions.h"
#include "utils.h"
#include "hash.h"
#include "grid.h"
#include "dasystem.h"

/* for the report */
#define HT_SIZE 10000

static inline double fsquare(double x)
{
    return x * x;
}

typedef struct {
    int nobs;
    double inn_f;
    double inn_f_abs;
    double inn_a;
    double inn_a_abs;
    double std_f;
    double std_a;
} stats;

/** Prints forecast and analysis observation statistics to stdout.
 */
void das_printobsstats(dasystem* das, int use_rmsd)
{
    observations* obs = das->obs;
    int i, j, otid;

    int inst, ni;
    double* inn_f_inst;
    double* inn_f_abs_inst;
    double* std_f_inst;
    double* inn_a_inst;
    double* inn_a_abs_inst;
    double* std_a_inst;
    int* nobs_inst;
    double (*func) (double) = (use_rmsd) ? fsquare : fabs;

    if (rank != 0)
        return;

    das_destandardise(das);

    ni = st_getsize(obs->instruments);
    inn_f_inst = malloc((ni + 1) * sizeof(double));
    inn_f_abs_inst = malloc((ni + 1) * sizeof(double));
    std_f_inst = malloc((ni + 1) * sizeof(double));
    inn_a_inst = malloc((ni + 1) * sizeof(double));
    inn_a_abs_inst = malloc((ni + 1) * sizeof(double));
    std_a_inst = malloc((ni + 1) * sizeof(double));
    nobs_inst = malloc((ni + 1) * sizeof(int));

    enkf_printf("    region obs.type   # obs.  %cfor.inn.%c %can.inn.%c   for.inn.   an.inn.  for.spread  an.spread\n", (use_rmsd) ? '[' : '|', (use_rmsd) ? ']' : '|', (use_rmsd) ? '[' : '|', (use_rmsd) ? ']' : '|');
    enkf_printf("    ------------------------------------------------------------------------------------------\n");

    for (i = 0; i < das->nregions; ++i) {
        region* r = &das->regions[i];

        enkf_printf("    %s\n", r->name);
        for (otid = 0; otid < obs->nobstypes; ++otid) {
            obstype* ot = &obs->obstypes[otid];
            stats rstats;
            stats* rzstats = NULL;
            int nzints = 0;
            zint* zints = NULL;

            int t1 = -INT_MAX;
            int t2 = -INT_MAX;
            int nt = 0;
            int t;
            double* inn_f_as = NULL;
            double* inn_f_abs_as = NULL;
            double* inn_a_as = NULL;
            double* inn_a_abs_as = NULL;
            double* std_f_as = NULL;
            double* std_a_as = NULL;
            int* nobs_as = NULL;

            if (!ot->issurface) {
                int mvid = model_getvarid(das->m, obs->obstypes[otid].varnames[0], 1);
                void* g = model_getvargrid(das->m, mvid);

                grid_getzints(g, &nzints, &zints);
                if (nzints > 0)
                    rzstats = calloc(nzints, sizeof(stats));
            }

            if (ot->isasync) {
                t1 = get_tshift(ot->time_min, ot->async_tstep, ot->async_centred);
                t2 = get_tshift(ot->time_max, ot->async_tstep, ot->async_centred);
                nt = t2 - t1 + 1;
                if (nt > 1) {
                    inn_f_as = calloc(nt, sizeof(double));
                    inn_f_abs_as = calloc(nt, sizeof(double));
                    inn_a_as = calloc(nt, sizeof(double));
                    inn_a_abs_as = calloc(nt, sizeof(double));
                    std_f_as = calloc(nt, sizeof(double));
                    std_a_as = calloc(nt, sizeof(double));
                    nobs_as = calloc(nt, sizeof(int));
                }
            }

            memset(&rstats, 0, sizeof(stats));
            memset(inn_f_inst, 0, (ni + 1) * sizeof(double));
            memset(inn_f_abs_inst, 0, (ni + 1) * sizeof(double));
            memset(std_f_inst, 0, (ni + 1) * sizeof(double));
            memset(inn_a_inst, 0, (ni + 1) * sizeof(double));
            memset(inn_a_abs_inst, 0, (ni + 1) * sizeof(double));
            memset(std_a_inst, 0, (ni + 1) * sizeof(double));
            memset(nobs_inst, 0, (ni + 1) * sizeof(int));

            for (j = 0; j < obs->nobs; ++j) {
                observation* o = &obs->data[j];

                if (o->status != STATUS_OK)
                    continue;

                if (o->type == otid && o->lat >= r->y1 && o->lat <= r->y2 && inloninterval(o->lon, r->x1, r->x2)) {
                    rstats.inn_f += das->s_f[j];
                    rstats.inn_f_abs += func(das->s_f[j]);
                    rstats.inn_a += das->s_a[j];
                    rstats.inn_a_abs += func(das->s_a[j]);
                    rstats.std_f += das->std_f[j];
                    rstats.std_a += das->std_a[j];
                    rstats.nobs++;
                    if (!ot->issurface) {
                        int ii;

                        for (ii = 0; ii < nzints; ++ii) {
                            if (o->depth >= zints[ii].z1 && o->depth < zints[ii].z2) {
                                rzstats[ii].inn_f += das->s_f[j];
                                rzstats[ii].inn_f_abs += func(das->s_f[j]);
                                rzstats[ii].inn_a += das->s_a[j];
                                rzstats[ii].inn_a_abs += func(das->s_a[j]);
                                rzstats[ii].std_f += das->std_f[j];
                                rzstats[ii].std_a += das->std_a[j];
                                rzstats[ii].nobs++;
                            }
                        }
                    }

                    if (ot->isasync && nt > 1) {
                        t = get_tshift(o->time, ot->async_tstep, ot->async_centred) - t1;
                        if (t < 0)
                            t = 0;
                        else if (t >= nt)
                            t = nt - 1;
                        inn_f_as[t] += das->s_f[j];
                        inn_f_abs_as[t] += func(das->s_f[j]);
                        inn_a_as[t] += das->s_a[j];
                        inn_a_abs_as[t] += func(das->s_a[j]);
                        std_f_as[t] += das->std_f[j];
                        std_a_as[t] += das->std_a[j];
                        nobs_as[t]++;
                    }

                    if (o->instrument >= 0) {
                        inn_f_inst[o->instrument] += das->s_f[j];
                        inn_f_abs_inst[o->instrument] += func(das->s_f[j]);
                        std_f_inst[o->instrument] += das->std_f[j];
                        inn_a_inst[o->instrument] += das->s_a[j];
                        inn_a_abs_inst[o->instrument] += func(das->s_a[j]);
                        std_a_inst[o->instrument] += das->std_a[j];
                        nobs_inst[o->instrument]++;
                    } else {
                        inn_f_inst[ni] += das->s_f[j];
                        inn_f_abs_inst[ni] += func(das->s_f[j]);
                        std_f_inst[ni] += das->std_f[j];
                        inn_a_inst[ni] += das->s_a[j];
                        inn_a_abs_inst[ni] += func(das->s_a[j]);
                        std_a_inst[ni] += das->std_a[j];
                        nobs_inst[ni]++;
                    }
                }
            }
            rstats.inn_f /= (double) rstats.nobs;
            rstats.inn_f_abs /= (double) rstats.nobs;
            rstats.inn_a /= (double) rstats.nobs;
            rstats.inn_a_abs /= (double) rstats.nobs;
            rstats.std_f /= (double) rstats.nobs;
            rstats.std_a /= (double) rstats.nobs;
            if (use_rmsd) {
                rstats.inn_f_abs = sqrt(rstats.inn_f_abs);
                rstats.inn_a_abs = sqrt(rstats.inn_a_abs);
            }

            enkf_printf("           %s      %8d%9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  \n", ot->name, rstats.nobs, rstats.inn_f_abs, rstats.inn_a_abs, rstats.inn_f, rstats.inn_a, rstats.std_f, rstats.std_a);

            /*
             * for asynchronous obs -- print stats for each time interval
             */
            if (ot->isasync && nt > 1 && rstats.nobs > 0) {
                for (t = 0; t < nt; ++t) {
                    inn_f_as[t] /= (double) nobs_as[t];
                    inn_f_abs_as[t] /= (double) nobs_as[t];
                    inn_a_as[t] /= (double) nobs_as[t];
                    inn_a_abs_as[t] /= (double) nobs_as[t];
                    std_f_as[t] /= (double) nobs_as[t];
                    std_a_as[t] /= (double) nobs_as[t];
                    if (use_rmsd) {
                        inn_f_abs_as[t] = sqrt(inn_f_abs_as[t]);
                        inn_a_abs_as[t] = sqrt(inn_a_abs_as[t]);
                    }
                    enkf_printf("           %3d      %8d%9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  \n", t1 + t, nobs_as[t], inn_f_abs_as[t], inn_a_abs_as[t], inn_f_as[t], inn_a_as[t], std_f_as[t], std_a_as[t]);
                }
            }

            /*
             * print stats for each intrument
             */
            for (inst = 0; inst <= ni; ++inst) {
                if (nobs_inst[inst] == 0)
                    continue;
                if (inst == ni && nobs_inst[inst] == rstats.nobs)
                    continue;
                inn_f_inst[inst] /= (double) nobs_inst[inst];
                inn_f_abs_inst[inst] /= (double) nobs_inst[inst];
                std_f_inst[inst] /= (double) nobs_inst[inst];
                inn_a_inst[inst] /= (double) nobs_inst[inst];
                inn_a_abs_inst[inst] /= (double) nobs_inst[inst];
                std_a_inst[inst] /= (double) nobs_inst[inst];
                if (use_rmsd) {
                    inn_f_abs_inst[inst] = sqrt(inn_f_abs_inst[inst]);
                    inn_a_abs_inst[inst] = sqrt(inn_a_abs_inst[inst]);
                }
                enkf_printf("             %-7s%8d%9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  \n", (inst < ni) ? st_findstringbyindex(obs->instruments, inst) : "N/A", nobs_inst[inst], inn_f_abs_inst[inst], inn_a_abs_inst[inst], inn_f_inst[inst], inn_a_inst[inst], std_f_inst[inst], std_a_inst[inst]);
            }

            /*
             * for non-surface obs -- print stats for each vertical interval
             */
            if (!ot->issurface && rstats.nobs > 0) {
                char tag[MAXSTRLEN];
                int ii;

                for (ii = 0; ii < nzints; ++ii) {
                    stats* s = &rzstats[ii];

                    s->inn_f /= (double) s->nobs;
                    s->inn_f_abs /= (double) s->nobs;
                    s->inn_a /= (double) s->nobs;
                    s->inn_a_abs /= (double) s->nobs;
                    s->std_f /= (double) s->nobs;
                    s->std_a /= (double) s->nobs;
                    if (use_rmsd) {
                        s->inn_f_abs = sqrt(s->inn_f_abs);
                        s->inn_a_abs = sqrt(s->inn_a_abs);
                    }
                    snprintf(tag, MAXSTRLEN, "%.0f-%.0fm", zints[ii].z1, zints[ii].z2);
                    enkf_printf("             %-9s%6d%9.3g  %9.3g  %9.3g  %9.3g  %9.3g  %9.3g  \n", tag, s->nobs, s->inn_f_abs, s->inn_a_abs, s->inn_f, s->inn_a, s->std_f, s->std_a);
                }
            }

            if (ot->isasync && nt > 1) {
                free(inn_f_as);
                free(inn_f_abs_as);
                free(inn_a_as);
                free(inn_a_abs_as);
                free(std_f_as);
                free(std_a_as);
                free(nobs_as);
            }
            if (nzints > 0)
                free(rzstats);
        }
    }

    free(inn_f_inst);
    free(inn_f_abs_inst);
    free(std_f_inst);
    free(inn_a_inst);
    free(inn_a_abs_inst);
    free(std_a_inst);
    free(nobs_inst);
}

typedef struct {
    int nobs;
    double inn_f;
    double inn_f_abs;
    double std_f;
} fstats;

/** Prints forecast observation statistics to the stdout.
 */
void das_printfobsstats(dasystem* das, int use_rmsd)
{
    observations* obs = das->obs;
    int i, j, otid;

    int inst, ni;
    double* inn_f_inst;
    double* inn_f_abs_inst;
    double* std_f_inst;
    int* nobs_inst;
    double (*func) (double) = (use_rmsd) ? fsquare : fabs;

    das_destandardise(das);

    if (rank != 0)
        return;

    ni = st_getsize(obs->instruments);
    inn_f_inst = malloc((ni + 1) * sizeof(double));
    inn_f_abs_inst = malloc((ni + 1) * sizeof(double));
    std_f_inst = malloc((ni + 1) * sizeof(double));
    nobs_inst = malloc((ni + 1) * sizeof(int));

    if (das->mode == MODE_ENKF || das->mode == MODE_HYBRID) {
        enkf_printf("    region obs.type   # obs.  %cfor.inn.%c for.inn.   for.spread\n", (use_rmsd) ? '[' : '|', (use_rmsd) ? ']' : '|');
        enkf_printf("    ----------------------------------------------------------\n");
    } else {
        enkf_printf("    region obs.type   # obs.  %cfor.inn.%c for.inn.\n", (use_rmsd) ? '[' : '|', (use_rmsd) ? ']' : '|');
        enkf_printf("    ---------------------------------------------\n");
    }

    for (i = 0; i < das->nregions; ++i) {
        region* r = &das->regions[i];

        enkf_printf("    %s\n", r->name);
        for (otid = 0; otid < obs->nobstypes; ++otid) {
            obstype* ot = &obs->obstypes[otid];
            fstats rstats;
            fstats* rzstats = NULL;
            int nzints = 0;
            zint* zints = NULL;

            int t1 = -INT_MAX;
            int t2 = -INT_MAX;
            int nt = 0;
            int t;
            double* inn_f_as = NULL;
            double* inn_f_abs_as = NULL;
            double* std_f_as = NULL;
            int* nobs_as = NULL;

            if (!ot->issurface) {
                int mvid = model_getvarid(das->m, obs->obstypes[otid].varnames[0], 1);
                void* g = model_getvargrid(das->m, mvid);

                grid_getzints(g, &nzints, &zints);
                if (nzints > 0)
                    rzstats = calloc(nzints, sizeof(fstats));
            }

            if (ot->isasync) {
                t1 = get_tshift(ot->time_min, ot->async_tstep, ot->async_centred);
                t2 = get_tshift(ot->time_max, ot->async_tstep, ot->async_centred);
                nt = t2 - t1 + 1;
                if (nt > 1) {
                    inn_f_as = calloc(nt, sizeof(double));
                    inn_f_abs_as = calloc(nt, sizeof(double));
                    if (das->mode == MODE_ENKF || das->mode == MODE_HYBRID)
                        std_f_as = calloc(nt, sizeof(double));
                    nobs_as = calloc(nt, sizeof(int));
                }
            }

            memset(&rstats, 0, sizeof(fstats));
            memset(inn_f_inst, 0, (ni + 1) * sizeof(double));
            memset(inn_f_abs_inst, 0, (ni + 1) * sizeof(double));
            memset(std_f_inst, 0, (ni + 1) * sizeof(double));
            memset(nobs_inst, 0, (ni + 1) * sizeof(int));

            for (j = 0; j < obs->nobs; ++j) {
                observation* o = &obs->data[j];

                if (o->status != STATUS_OK)
                    continue;

                if (o->type == otid && o->lat >= r->y1 && o->lat <= r->y2 && inloninterval(o->lon, r->x1, r->x2)) {
                    rstats.inn_f += das->s_f[j];
                    rstats.inn_f_abs += func(das->s_f[j]);
                    if (das->mode == MODE_ENKF || das->mode == MODE_HYBRID)
                        rstats.std_f += das->std_f[j];
                    rstats.nobs++;
                    if (!ot->issurface) {
                        int ii;

                        for (ii = 0; ii < nzints; ++ii) {
                            if (o->depth >= zints[ii].z1 && o->depth < zints[ii].z2) {
                                rzstats[ii].inn_f += das->s_f[j];
                                rzstats[ii].inn_f_abs += func(das->s_f[j]);
                                rzstats[ii].std_f += das->std_f[j];
                                rzstats[ii].nobs++;
                            }
                        }
                    }

                    if (ot->isasync && nt > 1) {
                        t = get_tshift(o->time, ot->async_tstep, ot->async_centred) - t1;
                        if (t < 0)
                            t = 0;
                        else if (t >= nt)
                            t = nt - 1;
                        inn_f_as[t] += das->s_f[j];
                        inn_f_abs_as[t] += func(das->s_f[j]);
                        if (das->mode == MODE_ENKF || das->mode == MODE_HYBRID)
                            std_f_as[t] += das->std_f[j];
                        nobs_as[t]++;
                    }

                    if (o->instrument >= 0) {
                        inn_f_inst[o->instrument] += das->s_f[j];
                        inn_f_abs_inst[o->instrument] += func(das->s_f[j]);
                        if (das->mode == MODE_ENKF || das->mode == MODE_HYBRID)
                            std_f_inst[o->instrument] += das->std_f[j];
                        nobs_inst[o->instrument]++;
                    } else {
                        inn_f_inst[ni] += das->s_f[j];
                        inn_f_abs_inst[ni] += func(das->s_f[j]);
                        if (das->mode == MODE_ENKF || das->mode == MODE_HYBRID)
                            std_f_inst[ni] += das->std_f[j];
                        nobs_inst[ni]++;
                    }
                }
            }
            rstats.inn_f /= (double) rstats.nobs;
            rstats.inn_f_abs /= (double) rstats.nobs;
            if (use_rmsd)
                rstats.inn_f_abs = sqrt(rstats.inn_f_abs);
            if (das->mode == MODE_ENKF || das->mode == MODE_HYBRID)
                rstats.std_f /= (double) rstats.nobs;

            if (das->mode == MODE_ENKF || das->mode == MODE_HYBRID)
                enkf_printf("           %s      %8d%9.3g  %9.3g  %9.3g  \n", ot->name, rstats.nobs, rstats.inn_f_abs, rstats.inn_f, rstats.std_f);
            else
                enkf_printf("           %s      %8d%9.3g  %9.3g  \n", ot->name, rstats.nobs, rstats.inn_f_abs, rstats.inn_f);

            /*
             * for asynchronous obs -- print stats for each time interval
             */
            if (ot->isasync && nt > 1 && rstats.nobs > 0) {
                for (t = 0; t < nt; ++t) {
                    inn_f_as[t] /= (double) nobs_as[t];
                    inn_f_abs_as[t] /= (double) nobs_as[t];
                    if (use_rmsd)
                        inn_f_abs_as[t] = sqrt(inn_f_abs_as[t]);
                    if (das->mode == MODE_ENKF || das->mode == MODE_HYBRID)
                        std_f_as[t] /= (double) nobs_as[t];
                    if (das->mode == MODE_ENKF || das->mode == MODE_HYBRID)
                        enkf_printf("           %3d      %8d%9.3g  %9.3g  %9.3g  \n", t1 + t, nobs_as[t], inn_f_abs_as[t], inn_f_as[t], std_f_as[t]);
                    else
                        enkf_printf("           %3d      %8d%9.3g  %9.3g  \n", t1 + t, nobs_as[t], inn_f_abs_as[t], inn_f_as[t]);
                }
            }

            /*
             * print stats for each intrument
             */
            for (inst = 0; inst <= ni; ++inst) {
                if (nobs_inst[inst] == 0)
                    continue;
                if (inst == ni && nobs_inst[inst] == rstats.nobs)
                    continue;
                inn_f_inst[inst] /= (double) nobs_inst[inst];
                inn_f_abs_inst[inst] /= (double) nobs_inst[inst];
                if (use_rmsd)
                    inn_f_abs_inst[inst] = sqrt(inn_f_abs_inst[inst]);
                if (das->mode == MODE_ENKF || das->mode == MODE_HYBRID)
                    std_f_inst[inst] /= (double) nobs_inst[inst];
                if (das->mode == MODE_ENKF || das->mode == MODE_HYBRID)
                    enkf_printf("             %-7s%8d%9.3g  %9.3g  %9.3g  \n", (inst < ni) ? st_findstringbyindex(obs->instruments, inst) : "N/A", nobs_inst[inst], inn_f_abs_inst[inst], inn_f_inst[inst], std_f_inst[inst]);
                else
                    enkf_printf("             %-7s%8d%9.3g  %9.3g  \n", (inst < ni) ? st_findstringbyindex(obs->instruments, inst) : "N/A", nobs_inst[inst], inn_f_abs_inst[inst], inn_f_inst[inst]);
            }

            /*
             * for non-surface obs -- print stats for each vertical interval
             */
            if (!ot->issurface && rstats.nobs > 0) {
                char tag[MAXSTRLEN];
                int ii;

                for (ii = 0; ii < nzints; ++ii) {
                    fstats* s = &rzstats[ii];

                    s->inn_f /= (double) s->nobs;
                    s->inn_f_abs /= (double) s->nobs;
                    if (use_rmsd)
                        s->inn_f_abs = sqrt(s->inn_f_abs);
                    if (das->mode == MODE_ENKF || das->mode == MODE_HYBRID)
                        s->std_f /= (double) s->nobs;

                    snprintf(tag, MAXSTRLEN, "%.0f-%.0fm", zints[ii].z1, zints[ii].z2);
                    if (das->mode == MODE_ENKF || das->mode == MODE_HYBRID)
                        enkf_printf("             %-9s%6d%9.3g  %9.3g  %9.3g  \n", tag, s->nobs, s->inn_f_abs, s->inn_f, s->std_f);
                    else
                        enkf_printf("             %-9s%6d%9.3g  %9.3g  \n", tag, s->nobs, s->inn_f_abs, s->inn_f);
                }
            }

            if (ot->isasync && nt > 1) {
                free(inn_f_as);
                free(inn_f_abs_as);
                if (das->mode == MODE_ENKF || das->mode == MODE_HYBRID)
                    free(std_f_as);
                free(nobs_as);
            }
            if (nzints > 0)
                free(rzstats);
        }
    }

    free(inn_f_inst);
    free(inn_f_abs_inst);
    free(std_f_inst);
    free(nobs_inst);
}

/**
 */
hashtable* das_getbatches(dasystem* das)
{
    observations* obs = das->obs;
    hashtable* batches = ht_create_i1s2(HT_SIZE);
    int key[2] = { -1, -1 };
    short* keys = (short*) key;
    int* nbatch = calloc(obs->nobstypes, sizeof(int));
    int nbatch_total, i;

    enkf_printf("  identifying observation batches:\n");
    nbatch_total = 0;
    for (i = 0; i < obs->nobs; ++i) {
        observation* o = &obs->data[i];

        if (o->status != STATUS_OK)
            continue;

        if (o->fid < 0 || o->batch < 0) /* this superob is not "clean" */
            continue;

        if (o->batch == key[0] && o->type == keys[2] && o->fid == keys[3])
            /*
             * this superob is from the same batch as previous 
             */
            continue;

        key[0] = o->batch;
        keys[2] = o->type;
        keys[3] = o->fid;

        if (ht_find(batches, key) != NULL)
            continue;

        ht_insert(batches, key, o);

        nbatch[o->type]++;
        nbatch_total++;
    }

    for (i = 0; i < obs->nobstypes; ++i)
        enkf_printf("    %-6s  %-4d\n", obs->obstypes[i].name, nbatch[i]);
    enkf_printf("    total:  %-4d\n", nbatch_total);

    free(nbatch);

    return batches;
}

/** (1) Calculate mean innovation and mean absolute innovation for each batch
 ** (2) Write those to FNAME_BATCHES
 ** (3) Detect bad batches
 ** (4) Write those to FNAME_BADBATCHES
 * @param das - dasystem 
 * @param batches - hash table with batches returned by das_getbatches()
 * @return - hash table with bad batches
 */
hashtable* das_processbatches(dasystem* das, hashtable* batches)
{
    observations* obs = das->obs;
    int nbatch = ht_getnentries(batches);
    hashtable* badbatches = NULL;
    int i;

    /*
     * batch stats
     */
    double* inn_f_abs;
    double* inn_f;
    observation** oo;
    int* nobs;

    if (rank != 0 || nbatch == 0)
        return badbatches;

    inn_f_abs = calloc(nbatch, sizeof(double));
    inn_f = calloc(nbatch, sizeof(double));
    oo = calloc(nbatch, sizeof(observation*));
    nobs = calloc(nbatch, sizeof(int));

    /*
     * get batch stats
     */
    for (i = 0; i < obs->nobs; ++i) {
        observation* o = &obs->data[i];
        int key_int[2];
        short* key_short = (short*) key_int;
        int id;

        if (o->status != STATUS_OK)
            continue;

        if (o->fid < 0 || o->batch < 0)
            continue;

        key_int[0] = o->batch;
        key_short[2] = o->type;
        key_short[3] = o->fid;

        id = ht_findid(batches, key_int);
        assert(id >= 0);

        inn_f_abs[id] += fabs(das->s_f[i]);
        inn_f[id] += das->s_f[i];
        nobs[id]++;

        if (oo[id] == NULL)
            oo[id] = ht_find(batches, key_int);
    }
    for (i = 0; i < nbatch; ++i) {
        inn_f_abs[i] /= (double) nobs[i];
        inn_f[i] /= (double) nobs[i];
    }

    /*
     * write batch stats to FNAME_BATCHES
     */
    if (rank == 0) {
        FILE* f = enkf_fopen(FNAME_BATCHES, "w");

        fprintf(f, "     id  obs.type  fid  batch  #obs.   |for.inn.|  for.inn. fname\n");
        for (i = 0; i < nbatch; ++i) {
            observation* o = oo[i];
            char* fname = st_findstringbyindex(obs->datafiles, o->fid);

            fprintf(f, "%7d    %-7s %-4d  %-4d   %-6d %7.3f  %9.3f   %s\n", i, obs->obstypes[o->type].name, o->fid, o->batch, nobs[i], inn_f_abs[i], inn_f[i], fname);
        }
        fclose(f);
    }

    /*
     * detect bad batches and write their stats to FNAME_BADBATCHES
     */
    if (das->nbadbatchspecs > 0) {
        FILE* f = NULL;
        int* nbadbatch = calloc(obs->nobstypes, sizeof(int));

        enkf_printf("  detecting bad batches:\n");
        badbatches = ht_create_i1s2(HT_SIZE / 10);
        if (rank == 0) {
            f = enkf_fopen(FNAME_BADBATCHES, "w");
            fprintf(f, "     id  obs.type  fid  batch  #obs.   |for.inn.|  for.inn. fname\n");
        }
        for (i = 0; i < nbatch; ++i) {
            observation* o = oo[i];
            int j;

            for (j = 0; j < das->nbadbatchspecs; ++j) {
                badbatchspec* bb = &das->badbatchspecs[j];

                if (strcmp(bb->obstype, obs->obstypes[o->type].name) == 0 && nobs[i] >= bb->minnobs && (fabs(inn_f[i]) >= bb->maxbias || inn_f_abs[i] >= bb->maxmad)) {
                    int key_int[2];
                    short* key_short = (short*) key_int;

                    key_int[0] = o->batch;
                    key_short[2] = o->type;
                    key_short[3] = o->fid;
                    ht_insert(badbatches, key_int, o);

                    nbadbatch[o->type]++;

                    if (rank == 0) {
                        char* fname = st_findstringbyindex(obs->datafiles, o->fid);

                        fprintf(f, "%7d    %-7s %-4d  %-4d   %-6d %7.3f  %9.3f   %s\n", i, obs->obstypes[o->type].name, o->fid, o->batch, nobs[i], inn_f_abs[i], inn_f[i], fname);
                    }
                }
            }
        }
        if (rank == 0)
            fclose(f);

        for (i = 0; i < obs->nobstypes; ++i)
            enkf_printf("    %-6s  %-4d\n", obs->obstypes[i].name, nbadbatch[i]);
        enkf_printf("    total:  %-4d\n", ht_getnentries(badbatches));

        free(nbadbatch);
    }

    free(inn_f_abs);
    free(inn_f);
    free(oo);
    free(nobs);

    return badbatches;
}
