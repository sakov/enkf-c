/******************************************************************************
 *
 * File:        obsstats.c        
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
#include "definitions.h"
#include "utils.h"
#include "hash.h"
#include "dasystem.h"

/* for the report */
#define DEPTH_SHALLOW 50.0
#define DEPTH_DEEP 500.0
#define HT_SIZE 5000

/** Prints forecast and analysis observation statistics to stdout.
 */
void das_printobsstats(dasystem* das)
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

    if (rank != 0)
        return;

    das_destandardise(das);

    ni = obs->instruments->n;
    inn_f_inst = malloc((ni + 1) * sizeof(double));
    inn_f_abs_inst = malloc((ni + 1) * sizeof(double));
    std_f_inst = malloc((ni + 1) * sizeof(double));
    inn_a_inst = malloc((ni + 1) * sizeof(double));
    inn_a_abs_inst = malloc((ni + 1) * sizeof(double));
    std_a_inst = malloc((ni + 1) * sizeof(double));
    nobs_inst = malloc((ni + 1) * sizeof(int));

    enkf_printf("    region obs.type   # obs.  |for.inn.| |an.inn.|   for.inn.   an.inn.  for.spread  an.spread\n");
    enkf_printf("    ------------------------------------------------------------------------------------------\n");

    for (i = 0; i < das->nregions; ++i) {
        region* r = &das->regions[i];

        enkf_printf("    %s\n", r->name);
        for (otid = 0; otid < obs->nobstypes; ++otid) {
            obstype* ot = &obs->obstypes[otid];
            double inn_f = 0.0;
            double inn_f_abs = 0.0;
            double inn_a = 0.0;
            double inn_a_abs = 0.0;
            double std_f = 0.0;
            double std_a = 0.0;
            int nobs = 0;
            double inn_f_surf = 0.0;
            double inn_f_abs_surf = 0.0;
            double inn_a_surf = 0.0;
            double inn_a_abs_surf = 0.0;
            double std_f_surf = 0.0;
            double std_a_surf = 0.0;
            int nobs_surf = 0;
            double inn_f_deep = 0.0;
            double inn_f_abs_deep = 0.0;
            double inn_a_deep = 0.0;
            double inn_a_abs_deep = 0.0;
            double std_f_deep = 0.0;
            double std_a_deep = 0.0;
            int nobs_deep = 0;

            int t1 = -MAXINT;
            int t2 = -MAXINT;
            int nt = 0;
            int t;
            double* inn_f_as = NULL;
            double* inn_f_abs_as = NULL;
            double* inn_a_as = NULL;
            double* inn_a_abs_as = NULL;
            double* std_f_as = NULL;
            double* std_a_as = NULL;
            int* nobs_as = NULL;

            if (ot->isasync) {
                t1 = get_tshift(ot->date_min, ot->async_tstep);
                t2 = get_tshift(ot->date_max, ot->async_tstep);
                nt = t2 - t1 + 1;
                inn_f_as = calloc(nt, sizeof(double));
                inn_f_abs_as = calloc(nt, sizeof(double));
                inn_a_as = calloc(nt, sizeof(double));
                inn_a_abs_as = calloc(nt, sizeof(double));
                std_f_as = calloc(nt, sizeof(double));
                std_a_as = calloc(nt, sizeof(double));
                nobs_as = calloc(nt, sizeof(int));
            }

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

                if (o->type == otid && o->lon >= r->lon1 && o->lon <= r->lon2 && o->lat >= r->lat1 && o->lat <= r->lat2) {
                    inn_f += das->s_f[j];
                    inn_f_abs += fabs(das->s_f[j]);
                    inn_a += das->s_a[j];
                    inn_a_abs += fabs(das->s_a[j]);
                    std_f += das->std_f[j];
                    std_a += das->std_a[j];
                    nobs++;
                    if (!ot->issurface) {
                        if (o->depth < DEPTH_SHALLOW) {
                            inn_f_surf += das->s_f[j];
                            inn_f_abs_surf += fabs(das->s_f[j]);
                            inn_a_surf += das->s_a[j];
                            inn_a_abs_surf += fabs(das->s_a[j]);
                            std_f_surf += das->std_f[j];
                            std_a_surf += das->std_a[j];
                            nobs_surf++;
                        } else if (o->depth > DEPTH_DEEP) {
                            inn_f_deep += das->s_f[j];
                            inn_f_abs_deep += fabs(das->s_f[j]);
                            inn_a_deep += das->s_a[j];
                            inn_a_abs_deep += fabs(das->s_a[j]);
                            std_f_deep += das->std_f[j];
                            std_a_deep += das->std_a[j];
                            nobs_deep++;
                        }
                    }

                    if (ot->isasync) {
                        t = get_tshift(o->date, ot->async_tstep) - t1;
                        assert(t >= 0 && t < nt);
                        inn_f_as[t] += das->s_f[j];
                        inn_f_abs_as[t] += fabs(das->s_f[j]);
                        inn_a_as[t] += das->s_a[j];
                        inn_a_abs_as[t] += fabs(das->s_a[j]);
                        std_f_as[t] += das->std_f[j];
                        std_a_as[t] += das->std_a[j];
                        nobs_as[t]++;
                    }

                    if (o->instrument >= 0) {
                        inn_f_inst[o->instrument] += das->s_f[j];
                        inn_f_abs_inst[o->instrument] += fabs(das->s_f[j]);
                        std_f_inst[o->instrument] += das->std_f[j];
                        inn_a_inst[o->instrument] += das->s_a[j];
                        inn_a_abs_inst[o->instrument] += fabs(das->s_a[j]);
                        std_a_inst[o->instrument] += das->std_a[j];
                        nobs_inst[o->instrument]++;
                    } else {
                        inn_f_inst[ni] += das->s_f[j];
                        inn_f_abs_inst[ni] += fabs(das->s_f[j]);
                        std_f_inst[ni] += das->std_f[j];
                        inn_a_inst[ni] += das->s_a[j];
                        inn_a_abs_inst[ni] += fabs(das->s_a[j]);
                        std_a_inst[ni] += das->std_a[j];
                        nobs_inst[ni]++;
                    }
                }
            }
            inn_f /= (double) nobs;
            inn_f_abs /= (double) nobs;
            inn_a /= (double) nobs;
            inn_a_abs /= (double) nobs;
            std_f /= (double) nobs;
            std_a /= (double) nobs;

            if (nobs > 0)
                enkf_printf("           %s      %8d%9.3f  %9.3f  %9.3f  %9.3f  %9.3f  %9.3f  \n", ot->name, nobs, inn_f_abs, inn_a_abs, inn_f, inn_a, std_f, std_a);

            if (ot->isasync && nobs > 0) {
                for (t = 0; t < nt; ++t) {
                    inn_f_as[t] /= (double) nobs_as[t];
                    inn_f_abs_as[t] /= (double) nobs_as[t];
                    inn_a_as[t] /= (double) nobs_as[t];
                    inn_a_abs_as[t] /= (double) nobs_as[t];
                    std_f_as[t] /= (double) nobs_as[t];
                    std_a_as[t] /= (double) nobs_as[t];
                    enkf_printf("           %3d      %8d%9.3f  %9.3f  %9.3f  %9.3f  %9.3f  %9.3f  \n", t1 + t, nobs_as[t], inn_f_abs_as[t], inn_a_abs_as[t], inn_f_as[t], inn_a_as[t], std_f_as[t], std_a_as[t]);
                }
            }

            for (inst = 0; inst <= ni; ++inst) {
                if (nobs_inst[inst] == 0)
                    continue;
                if (inst == ni && nobs_inst[inst] == nobs)
                    continue;
                inn_f_inst[inst] /= (double) nobs_inst[inst];
                inn_f_abs_inst[inst] /= (double) nobs_inst[inst];
                std_f_inst[inst] /= (double) nobs_inst[inst];
                inn_a_inst[inst] /= (double) nobs_inst[inst];
                inn_a_abs_inst[inst] /= (double) nobs_inst[inst];
                std_a_inst[inst] /= (double) nobs_inst[inst];
                enkf_printf("             %-7s%8d%9.3f  %9.3f  %9.3f  %9.3f  %9.3f  %9.3f  \n", (inst < ni) ? st_findstringbyindex(obs->instruments, inst) : "N/A", nobs_inst[inst], inn_f_abs_inst[inst], inn_a_abs_inst[inst], inn_f_inst[inst], inn_a_inst[inst], std_f_inst[inst], std_a_inst[inst]);
            }

            if (!ot->issurface && nobs > 0) {
                inn_f_surf /= (double) nobs_surf;
                inn_f_abs_surf /= (double) nobs_surf;
                inn_a_surf /= (double) nobs_surf;
                inn_a_abs_surf /= (double) nobs_surf;
                std_f_surf /= (double) nobs_surf;
                std_a_surf /= (double) nobs_surf;
                inn_f_deep /= (double) nobs_deep;
                inn_f_abs_deep /= (double) nobs_deep;
                inn_a_deep /= (double) nobs_deep;
                inn_a_abs_deep /= (double) nobs_deep;
                std_f_deep /= (double) nobs_deep;
                std_a_deep /= (double) nobs_deep;

                enkf_printf("             0-%.0fm  %8d%9.3f  %9.3f  %9.3f  %9.3f  %9.3f  %9.3f  \n", DEPTH_SHALLOW, nobs_surf, inn_f_abs_surf, inn_a_abs_surf, inn_f_surf, inn_a_surf, std_f_surf, std_a_surf);
                enkf_printf("             >%.0fm  %8d%9.3f  %9.3f  %9.3f  %9.3f  %9.3f  %9.3f  \n", DEPTH_DEEP, nobs_deep, inn_f_abs_deep, inn_a_abs_deep, inn_f_deep, inn_a_deep, std_f_deep, std_a_deep);
            }

            if (ot->isasync) {
                free(inn_f_as);
                free(inn_f_abs_as);
                free(inn_a_as);
                free(inn_a_abs_as);
                free(std_f_as);
                free(std_a_as);
                free(nobs_as);
            }
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

/** Prints forecast observation statistics to the stdout.
 */
void das_printfobsstats(dasystem* das)
{
    observations* obs = das->obs;
    int i, j, otid;

    int inst, ni;
    double* inn_f_inst;
    double* inn_f_abs_inst;
    double* std_f_inst;
    int* nobs_inst;

    if (rank != 0)
        return;

    das_destandardise(das);

    ni = obs->instruments->n;
    inn_f_inst = malloc((ni + 1) * sizeof(double));
    inn_f_abs_inst = malloc((ni + 1) * sizeof(double));
    std_f_inst = malloc((ni + 1) * sizeof(double));
    nobs_inst = malloc((ni + 1) * sizeof(int));

    if (das->mode == MODE_ENKF) {
        enkf_printf("    region obs.type   # obs.  |for.inn.| for.inn.   for.spread\n");
        enkf_printf("    ----------------------------------------------------------\n");
    } else {
        enkf_printf("    region obs.type   # obs.  |for.inn.| for.inn.\n");
        enkf_printf("    ---------------------------------------------\n");
    }

    for (i = 0; i < das->nregions; ++i) {
        region* r = &das->regions[i];

        enkf_printf("    %s\n", r->name);
        for (otid = 0; otid < obs->nobstypes; ++otid) {
            obstype* ot = &obs->obstypes[otid];
            double inn_f = 0.0;
            double inn_f_abs = 0.0;
            double std_f = 0.0;
            int nobs = 0;
            double inn_f_surf = 0.0;
            double inn_f_abs_surf = 0.0;
            double std_f_surf = 0.0;
            int nobs_surf = 0;
            double inn_f_deep = 0.0;
            double inn_f_abs_deep = 0.0;
            double std_f_deep = 0.0;
            int nobs_deep = 0;

            int t1 = -MAXINT;
            int t2 = -MAXINT;
            int nt = 0;
            int t;
            double* inn_f_as = NULL;
            double* inn_f_abs_as = NULL;
            double* std_f_as = NULL;
            int* nobs_as = NULL;

            if (ot->isasync) {
                t1 = get_tshift(ot->date_min, ot->async_tstep);
                t2 = get_tshift(ot->date_max, ot->async_tstep);
                nt = t2 - t1 + 1;
                inn_f_as = calloc(nt, sizeof(double));
                inn_f_abs_as = calloc(nt, sizeof(double));
                if (das->mode == MODE_ENKF)
                    std_f_as = calloc(nt, sizeof(double));
                nobs_as = calloc(nt, sizeof(int));
            }

            memset(inn_f_inst, 0, (ni + 1) * sizeof(double));
            memset(inn_f_abs_inst, 0, (ni + 1) * sizeof(double));
            memset(std_f_inst, 0, (ni + 1) * sizeof(double));
            memset(nobs_inst, 0, (ni + 1) * sizeof(int));

            for (j = 0; j < obs->nobs; ++j) {
                observation* o = &obs->data[j];

                if (o->status != STATUS_OK)
                    continue;

                if (o->type == otid && o->lon >= r->lon1 && o->lon <= r->lon2 && o->lat >= r->lat1 && o->lat <= r->lat2) {
                    inn_f += das->s_f[j];
                    inn_f_abs += fabs(das->s_f[j]);
                    if (das->mode == MODE_ENKF)
                        std_f += das->std_f[j];
                    nobs++;
                    if (!ot->issurface) {
                        if (o->depth < DEPTH_SHALLOW) {
                            inn_f_surf += das->s_f[j];
                            inn_f_abs_surf += fabs(das->s_f[j]);
                            if (das->mode == MODE_ENKF)
                                std_f_surf += das->std_f[j];
                            nobs_surf++;
                        } else if (o->depth > DEPTH_DEEP) {
                            inn_f_deep += das->s_f[j];
                            inn_f_abs_deep += fabs(das->s_f[j]);
                            if (das->mode == MODE_ENKF)
                                std_f_deep += das->std_f[j];
                            nobs_deep++;
                        }
                    }

                    if (ot->isasync) {
                        t = get_tshift(o->date, ot->async_tstep) - t1;
                        assert(t >= 0 && t < nt);
                        inn_f_as[t] += das->s_f[j];
                        inn_f_abs_as[t] += fabs(das->s_f[j]);
                        if (das->mode == MODE_ENKF)
                            std_f_as[t] += das->std_f[j];
                        nobs_as[t]++;
                    }

                    if (o->instrument >= 0) {
                        inn_f_inst[o->instrument] += das->s_f[j];
                        inn_f_abs_inst[o->instrument] += fabs(das->s_f[j]);
                        if (das->mode == MODE_ENKF)
                            std_f_inst[o->instrument] += das->std_f[j];
                        nobs_inst[o->instrument]++;
                    } else {
                        inn_f_inst[ni] += das->s_f[j];
                        inn_f_abs_inst[ni] += fabs(das->s_f[j]);
                        if (das->mode == MODE_ENKF)
                            std_f_inst[ni] += das->std_f[j];
                        nobs_inst[ni]++;
                    }
                }
            }
            inn_f /= (double) nobs;
            inn_f_abs /= (double) nobs;
            if (das->mode == MODE_ENKF)
                std_f /= (double) nobs;

            if (nobs > 0) {
                if (das->mode == MODE_ENKF)
                    enkf_printf("           %s      %8d%9.3f  %9.3f  %9.3f  \n", ot->name, nobs, inn_f_abs, inn_f, std_f);
                else
                    enkf_printf("           %s      %8d%9.3f  %9.3f  \n", ot->name, nobs, inn_f_abs, inn_f);
            }

            if (ot->isasync && nobs > 0) {
                for (t = 0; t < nt; ++t) {
                    inn_f_as[t] /= (double) nobs_as[t];
                    inn_f_abs_as[t] /= (double) nobs_as[t];
                    if (das->mode == MODE_ENKF)
                        std_f_as[t] /= (double) nobs_as[t];
                    if (das->mode == MODE_ENKF)
                        enkf_printf("           %3d      %8d%9.3f  %9.3f  %9.3f  \n", t1 + t, nobs_as[t], inn_f_abs_as[t], inn_f_as[t], std_f_as[t]);
                    else
                        enkf_printf("           %3d      %8d%9.3f  %9.3f  \n", t1 + t, nobs_as[t], inn_f_abs_as[t], inn_f_as[t]);
                }
            }

            for (inst = 0; inst <= ni; ++inst) {
                if (nobs_inst[inst] == 0)
                    continue;
                if (inst == ni && nobs_inst[inst] == nobs)
                    continue;
                inn_f_inst[inst] /= (double) nobs_inst[inst];
                inn_f_abs_inst[inst] /= (double) nobs_inst[inst];
                if (das->mode == MODE_ENKF)
                    std_f_inst[inst] /= (double) nobs_inst[inst];
                if (das->mode == MODE_ENKF)
                    enkf_printf("             %-7s%8d%9.3f  %9.3f  %9.3f  \n", (inst < ni) ? st_findstringbyindex(obs->instruments, inst) : "N/A", nobs_inst[inst], inn_f_abs_inst[inst], inn_f_inst[inst], std_f_inst[inst]);
                else
                    enkf_printf("             %-7s%8d%9.3f  %9.3f  \n", (inst < ni) ? st_findstringbyindex(obs->instruments, inst) : "N/A", nobs_inst[inst], inn_f_abs_inst[inst], inn_f_inst[inst]);
            }

            if (!ot->issurface && nobs > 0) {
                inn_f_surf /= (double) nobs_surf;
                inn_f_abs_surf /= (double) nobs_surf;
                if (das->mode == MODE_ENKF)
                    std_f_surf /= (double) nobs_surf;
                inn_f_deep /= (double) nobs_deep;
                inn_f_abs_deep /= (double) nobs_deep;
                if (das->mode == MODE_ENKF)
                    std_f_deep /= (double) nobs_deep;

                if (das->mode == MODE_ENKF) {
                    enkf_printf("             0-%.0fm  %8d%9.3f  %9.3f  %9.3f  \n", DEPTH_SHALLOW, nobs_surf, inn_f_abs_surf, inn_f_surf, std_f_surf);
                    enkf_printf("             >%.0fm  %8d%9.3f  %9.3f  %9.3f  \n", DEPTH_DEEP, nobs_deep, inn_f_abs_deep, inn_f_deep, std_f_deep);
                } else {
                    enkf_printf("             0-%.0fm  %8d%9.3f  %9.3f  \n", DEPTH_SHALLOW, nobs_surf, inn_f_abs_surf, inn_f_surf);
                    enkf_printf("             >%.0fm  %8d%9.3f  %9.3f  \n", DEPTH_DEEP, nobs_deep, inn_f_abs_deep, inn_f_deep);
                }
            }

            if (ot->isasync) {
                free(inn_f_as);
                free(inn_f_abs_as);
                if (das->mode == MODE_ENKF)
                    free(std_f_as);
                free(nobs_as);
            }
        }
    }

    free(inn_f_inst);
    free(inn_f_abs_inst);
    free(std_f_inst);
    free(nobs_inst);
}

/** Calculate and print global biases for each batch of observations.
 */
void das_calcbatchstats(dasystem* das, int doprint)
{
    observations* obs = das->obs;
    hashtable* ht = ht_create_i1s2(HT_SIZE);
    int key[2] = { -1, -1 };
    short* keys = (short*) key;
    int* nbatches = calloc(obs->nobstypes, sizeof(int));
    int n, i;

    /*
     * per batch data
     */
    double* inn_f_abs;
    double* inn_f;
    observation** oo;
    int* nobs;

    das_destandardise(das);

    /*
     * calculate the number of batches
     */
    n = 0;
    for (i = 0; i < obs->nobs; ++i) {
        observation* o = &obs->data[i];

        if (o->status != STATUS_OK)
            continue;

        if (o->fid < 0 || o->batch < 0) /* this superob is not "clean" */
            continue;

        if (o->fid == key[0] && o->batch == key[1])     /* same as previous */
            continue;

        key[0] = o->batch;
        keys[2] = o->type;
        keys[3] = o->fid;

        if (ht_find(ht, key) != NULL)
            continue;

        ht_insert(ht, key, o);
        nbatches[o->type]++;
    }

    n = ht_getnentries(ht);
    enkf_printf("  number of batches:\n");
    for (i = 0; i < obs->nobstypes; ++i)
        enkf_printf("    %6s  %4d\n", obs->obstypes[i].name, nbatches[i]);
    enkf_printf("    total:  %4d\n", n);

    if (n > 0) {
        inn_f_abs = calloc(n, sizeof(double));
        inn_f = calloc(n, sizeof(double));
        oo = calloc(n, sizeof(observation*));
        nobs = calloc(n, sizeof(int));

        for (i = 0; i < obs->nobs; ++i) {
            observation* o = &obs->data[i];
            int id;

            if (o->status != STATUS_OK)
                continue;

            if (o->fid < 0 || o->batch < 0)
                continue;

            key[0] = o->batch;
            keys[2] = o->type;
            keys[3] = o->fid;

            id = ht_findid(ht, key);
            assert(id >= 0);

            inn_f_abs[id] += fabs(das->s_f[i]);
            inn_f[id] += das->s_f[i];
            nobs[id]++;

            if (oo[id] == NULL)
                oo[id] = ht_find(ht, key);
        }

        for (i = 0; i < n; ++i) {
            inn_f_abs[i] /= (double) nobs[i];
            inn_f[i] /= (double) nobs[i];
        }

        /*
         * print batch stats
         */
        if (doprint) {
            enkf_printf("  batch statistics:\n");
            enkf_printf("     id  obs.type  fid  batch  # obs.  |for.inn.|  for.inn.\n");
            for (i = 0; i < n; ++i) {
                observation* o = oo[i];

                enkf_printf("%7d    %-7s %-4d  %-5d   %-5d %8.3f  %9.3f\n", i, obs->obstypes[o->type].name, o->fid, o->batch, nobs[i], inn_f_abs[i], inn_f[i]);
            }
        }

        /*
         * identify and report bad batches
         */
        if (das->nbadbatchspecs > 0) {
            FILE* f = enkf_fopen(FNAME_BADBATCHES, "w");

            enkf_printf("  bad batches:\n");
            for (i = 0; i < n; ++i) {
                observation* o = oo[i];
                int j;

                for (j = 0; j < das->nbadbatchspecs; ++j) {
                    badbatchspec* bb = &das->badbatchspecs[j];

                    if (strcmp(bb->obstype, obs->obstypes[o->type].name) == 0 && nobs[i] >= bb->minnobs && (fabs(inn_f[i]) >= bb->maxbias || inn_f_abs[i] >= bb->maxmad)) {
                        char* fname = st_findstringbyindex(obs->datafiles, o->fid);

                        enkf_printf("    %s %s %d %d %.3f %.3f\n", bb->obstype, fname, o->fid, o->batch, inn_f[i], inn_f_abs[i]);
                        fprintf(f, "%s %s %d %d %.3f %.3f\n", bb->obstype, fname, o->fid, o->batch, inn_f[i], inn_f_abs[i]);
                    }
                }
            }
            fclose(f);
        }

        free(inn_f_abs);
        free(inn_f);
        free(oo);
        free(nobs);
    }

    ht_destroy(ht);
    free(nbatches);
}
