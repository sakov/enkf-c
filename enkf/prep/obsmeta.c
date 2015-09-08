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
void obsmeta_read(enkfprm* prm, int* nmeta, obsmeta** meta)
{
    char* fname = prm->obsprm;
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
                enkf_quit("%s, l.%d: PRODUCT not specified", fname, line);

            *meta = realloc(*meta, (*nmeta + 1) * sizeof(obsmeta));
            m = &(*meta)[*nmeta];
            obsmeta_init(m);
            m->product = strdup(token);
            (*nmeta)++;
        } else if (strcasecmp(token, "READER") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: READER not specified", fname, line);
            m->reader = strdup(token);
        } else if (strcasecmp(token, "TYPE") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: TYPE not specified", fname, line);
            m->type = strdup(token);
        } else if (strcasecmp(token, "FILE") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: FILE not specified", fname, line);
            obsmeta_addfname(m, token);
        } else if (strcasecmp(token, "ERROR_STD") == 0) {
            metastd* now = NULL;
            double std;

            m->stds = realloc(m->stds, (m->nstds + 1) * sizeof(metastd));
            now = &m->stds[m->nstds];

            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: STD not specified", fname, line);

            if (str2double(token, &std)) {
                if (!isfinite(std) || std < 0.0)
                    enkf_quit("%s, l.%d: invalid STD value", fname, line);
                now->type = STDTYPE_VALUE;
                now->data = malloc(sizeof(double));
                ((double*) now->data)[0] = std;
                now->varname = NULL;
            } else {
                now->type = STDTYPE_FILE;
                now->data = strdup(token);
                if ((token = strtok(NULL, seps)) == NULL)
                    enkf_quit("%s, l.%d: variable name not specified", fname, line);
                now->varname = strdup(token);
            }
            if ((token = strtok(NULL, seps)) == NULL)
                now->op = ARITHMETIC_EQ;
            else if (strncasecmp(token, "EQ", 2) == 0)
                now->op = ARITHMETIC_EQ;
            else if (strncasecmp(token, "PL", 2) == 0)
                now->op = ARITHMETIC_PLUS;
            else if (strncasecmp(token, "MU", 2) == 0)
                now->op = ARITHMETIC_MULT;
            else if (strncasecmp(token, "MI", 2) == 0)
                now->op = ARITHMETIC_MIN;
            else if (strncasecmp(token, "MA", 2) == 0)
                now->op = ARITHMETIC_MAX;
            else
                enkf_quit("%s, l.%d:, unknown operation", fname, line);
            m->nstds++;
        } else if (strcasecmp(token, "PARAMETER") == 0) {
            metapar* now = NULL;

            m->pars = realloc(m->pars, (m->npars + 1) * sizeof(metapar));
            now = &m->pars[m->npars];

            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: parameter name not specified (expected: PARAMETER <name> = <value>)", fname, line);
            now->name = strdup(token);
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: parameter value not specified (expected: PARAMETER <name> = <value>)", fname, line);
            now->value = strdup(token);
            m->npars++;
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
            char operstr[MAXSTRLEN] = "";
            metastd* std = &m->stds[j];

            if (std->op == ARITHMETIC_EQ)
                strcpy(operstr, "EQUAL");
            else if (std->op == ARITHMETIC_PLUS)
                strcpy(operstr, "PLUS");
            else if (std->op == ARITHMETIC_MULT)
                strcpy(operstr, "MULT");
            else if (std->op == ARITHMETIC_MIN)
                strcpy(operstr, "MIN");
            else if (std->op == ARITHMETIC_MAX)
                strcpy(operstr, "MAX");

            if (std->type == STDTYPE_VALUE)
                enkf_printf("      ERROR_STD = %.3g, operation = %s\n", ((double*) std->data)[0], operstr);
            else if (std->type == STDTYPE_FILE)
                enkf_printf("      ERROR_STD = %s %s, operation = %s\n", (char*) std->data, std->varname, operstr);
        }
        for (j = 0; j < m->npars; ++j)
            enkf_printf("      PARAMETER %s = %s\n", m->pars[j].name, m->pars[j].value);
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
void obsmeta_destroy(int n, obsmeta meta[])
{
    int i, j;

    for (i = 0; i < n; ++i) {
        obsmeta* m = &meta[i];

        free(m->product);
        free(m->type);
        if (m->nfiles > 0) {
            for (j = 0; j < m->nfiles; ++j)
                free(m->fnames[j]);
            free(m->fnames);
        }
        if (m->nstds > 0) {
            for (j = 0; j < m->nstds; ++j) {
                metastd* std = &m->stds[j];

                if (std->data != NULL)
                    free(std->data);
                if (std->varname != NULL)
                    free(std->varname);
            }
            free(m->stds);
        }
        if (m->npars > 0) {
            for (j = 0; j < m->npars; ++j) {
                free(m->pars[j].name);
                free(m->pars[j].value);
            }
            free(m->pars);
        }
    }
    free(meta);
}

/**
 */
void obsmeta_describeprm(void)
{
    enkf_printf("\n");
    enkf_printf("  Observation data parameter file format:\n");
    enkf_printf("\n");
    enkf_printf("    PRODUCT   = <product>\n");
    enkf_printf("    READER    = <reader>\n");
    enkf_printf("    TYPE      = <observation type>\n");
    enkf_printf("    FILE      = <data file wildcard> \n");
    enkf_printf("    ...\n");
    enkf_printf("  [ ERROR_STD = { <value> | <data file> } [ EQ* | PL | MU | MI | MA ] ]\n");
    enkf_printf("    ...\n");
    enkf_printf("  [ PARAMETER <name> = <value> ]\n");
    enkf_printf("    ...\n");
    enkf_printf("\n");
    enkf_printf("  [ <more of the above blocks> ]\n");
    enkf_printf("\n");
    enkf_printf("  Notes:\n");
    enkf_printf("    1. { ... | ... | ... } denotes the list of possible choices\n");
    enkf_printf("    2. [ ... ] denotes an optional input\n");
    enkf_printf("    3. * denotes the default value\n");
    enkf_printf("    4. < ... > denotes a description of an entry\n");
    enkf_printf("    5. ... denotes repeating the previous item an arbitrary number of times\n");
    enkf_printf("\n");
}
