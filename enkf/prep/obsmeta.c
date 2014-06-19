/******************************************************************************
 *
 * File:        obsmeta.c        
 *
 * Created:     12/2012
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description:
 *
 * Revisions:
 *
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "definitions.h"
#include "utils.h"
#include "enkfprm.h"
#include "obsmeta.h"

#define OBSMETA_NFILES_INC 10
#define NOBSTYPES_INC 10

/**
 */
static void obsmeta_init(obsmeta* meta)
{
    memset(meta, 0, sizeof(obsmeta));
}

/**
 */
static void obsmeta_addfname(obsmeta* meta, char fname[])
{
    if (meta->nfiles % OBSMETA_NFILES_INC == 0)
        meta->fnames = realloc(meta->fnames, (meta->nfiles + OBSMETA_NFILES_INC) * sizeof(char*));
    meta->fnames[meta->nfiles] = strdup(fname);
    meta->nfiles++;
}

/**
 */
void read_obsmeta(enkfprm* prm, int* nmeta, obsmeta** meta)
{
    char* fname = prm->obsspec;
    FILE* f = NULL;
    char buf[MAXSTRLEN];
    int line;
    obsmeta* m = NULL;
    int i, j;

    *nmeta = 0;
    *meta = NULL;

    f = enkf_fopen(fname, "r");

    line = 0;
    while (fgets(buf, MAXSTRLEN, f) != NULL) {
        char seps[] = " =\t\n";
        char* token;

        line++;
        if (buf[0] == '#')
            continue;
        if ((token = strtok(buf, seps)) == NULL)
            continue;
        if (strcasecmp(token, "PRODUCT") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: PRODUCT entered but not specified", fname, line);
            else {
                *meta = realloc(*meta, (*nmeta + 1) * sizeof(obsmeta));
                m = &(*meta)[*nmeta];
                obsmeta_init(m);
                m->product = strdup(token);
                (*nmeta)++;
            }
        } else if (strcasecmp(token, "READER") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: READER entered but not specified", fname, line);
            else
                m->reader = strdup(token);
        } else if (strcasecmp(token, "TYPE") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: TYPE entered but not specified", fname, line);
            else {
                m->type = strdup(token);
            }
        } else if (strcasecmp(token, "FILE") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: FILE entered but not specified", fname, line);
            else
                obsmeta_addfname(m, token);
        } else if (strcasecmp(token, "ERROR_STD") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: STD entered but not specified", fname, line);
            else {
                double std;
                m->stdtypes = realloc(m->stdtypes, (m->nstds + 1) * sizeof(int));
                m->stds = realloc(m->stds, (m->nstds + 1) * sizeof(void*));
                m->varnames = realloc(m->varnames, (m->nstds + 1) * sizeof(char*));
                if (str2double(token, &std)) {
                    if (!isfinite(std) || std < 0.0)
                        enkf_quit("%s, l.%d: invalid STD value", fname, line);
                    m->stdtypes[m->nstds] = STDTYPE_VALUE;
                    m->stds[m->nstds] = malloc(sizeof(double));
                    *((double*) m->stds[m->nstds]) = std;
                    m->varnames[m->nstds] = NULL;
                } else {
                    m->stdtypes[m->nstds] = STDTYPE_FILE;

                    m->stds[m->nstds] = strdup(token);
                    if ((token = strtok(NULL, seps)) == NULL)
                        enkf_quit("%s, l.%d: variable name not specified", fname, line);
                    m->varnames[m->nstds] = strdup(token);
                }
                m->nstds++;
            }
        } else
            enkf_quit("%s, l.%d: unknown token \"%s\"", fname, line, token);
    }

    fclose(f);

    /*
     * print summary 
     */
    for (i = 0; i < (*nmeta); ++i) {
        obsmeta* m = &(*meta)[i];

        enkf_printf("    PRODUCT = %s\n", m->product);
        if (m->reader == NULL) {
            m->reader = strdup("standard");
            enkf_printf("      (assumed) READER = %s\n", m->reader);
        } else
            enkf_printf("      READER = %s\n", m->reader);
        if (m->type == NULL)
            enkf_quit("%s: obsservation type not specified for product \"%s\"", fname, m->product);
        enkf_printf("      TYPE = %s\n", m->type);
        for (j = 0; j < m->nfiles; ++j)
            enkf_printf("      File:        obsmeta.c = %s\n", m->fnames[j]);
        for (j = 0; j < m->nstds; ++j) {
            if (m->stdtypes[j] == STDTYPE_VALUE)
                enkf_printf("      ERROR_STD = %.3g\n", *((double*) m->stds[j]));
            else if (m->stdtypes[j] == STDTYPE_FILE)
                enkf_printf("      ERROR_STD = %s %s\n", (char*) m->stds[j], m->varnames[j]);
        }
    }

    /*
     * check asynchronous obs types 
     */
    do {
        for (j = 0; j < prm->nasync; ++j) {
            for (i = 0; i < (*nmeta); ++i) {
                obsmeta* m = &(*meta)[i];

                if (strcmp(m->type, prm->async_types[j]) == 0)
                    break;
            }
            if (i == (*nmeta)) {
                int k;

                enkf_printf("    WARNING: %s: asynchronous type \"%s\" not in obs\n", prm->fname, prm->async_types[j]);
                /*
                 * eliminate redundant entry 
                 */
                for (k = j; k < prm->nasync - 1; ++k)
                    prm->async_types[k] = prm->async_types[k + 1];
                prm->nasync--;
            }
        }
    } while (j != prm->nasync);
}

/**
 */
void clean_obsmeta(int n, obsmeta meta[])
{
    int i, j;

    for (i = 0; i < n; ++i) {
        obsmeta* m = &meta[i];

        free(m->product);
        free(m->type);
        for (j = 0; j < m->nfiles; ++j)
            free(m->fnames[j]);
        free(m->fnames);
        for (j = 0; j < m->nstds; ++j) {
            free(m->stds[j]);
            if (m->varnames[j] != NULL)
                free(m->varnames[j]);
        }
        free(m->stdtypes);
        free(m->stds);
        free(m->varnames);
    }
    free(meta);
}
