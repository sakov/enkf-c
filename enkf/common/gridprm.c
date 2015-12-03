/******************************************************************************
 *
 * File:        gridprm.c
 *
 * Created:     17/12/2014
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description:
 *
 * Revisions:
 *
 *****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "definitions.h"
#include "grid.h"
#include "gridprm.h"
#include "utils.h"

typedef struct {
    char* vtype_tag;
    int vtype;
    char* levelvarnameentry;
} gridvtype_entry;

gridvtype_entry allgridvtypeentries[] = {
    {"Z", GRIDVTYPE_Z, "NUMLEVELSVARNAME"},
    {"SIGMA", GRIDVTYPE_SIGMA, "MASKVARNAME"}
};
int ngridvtypeentries = sizeof(allgridvtypeentries) / sizeof(gridvtype_entry);

/**
 */
void gridprm_create(char* fname, int* ngrid, gridprm** prm)
{
    gridprm* now = NULL;
    FILE* f = NULL;
    char buf[MAXSTRLEN];
    int line;
    int i;

    assert(*ngrid == 0 && *prm == NULL);

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
        if (strcasecmp(token, "NAME") == 0) {
            (*ngrid)++;
            *prm = realloc(*prm, sizeof(gridprm) * (*ngrid));
            now = &(*prm)[*ngrid - 1];
            memset(now, 0, sizeof(gridprm));
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: NAME not specified", fname, line);
            else
                now->name = strdup(token);
            continue;
        } else if (now == NULL)
            enkf_quit("%s, l.%d: entry NAME is required to start grid initialisation", fname, line);

        if (strcasecmp(token, "VTYPE") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: VTYPE not specified", fname, line);
            else if (now->vtype != NULL)
                enkf_quit("%s, l.%d: VTYPE specified twice", fname, line);
            else
                now->vtype = strdup(token);

            for (i = 0; i < ngridvtypeentries; ++i)
                if (strcasecmp(allgridvtypeentries[i].vtype_tag, token) == 0)
                    break;
            if (i == ngridvtypeentries)
                enkf_quit("%s, l %d: VTYPE \"%s\" is unknown", fname, line);
            if (allgridvtypeentries[i].levelvarnameentry != NULL)
                now->levelvarnameentry = strdup(allgridvtypeentries[i].levelvarnameentry);
        } else if (strcasecmp(token, "DATA") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: DATA not specified", fname, line);
            else if (now->fname != NULL)
                enkf_quit("%s, l.%d: DATA specified twice", fname, line);
            else
                now->fname = strdup(token);
        } else if (strcasecmp(token, "XDIMNAME") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: XDIMNAME not specified", fname, line);
            else if (now->xdimname != NULL)
                enkf_quit("%s, l.%d: XDIMNAME specified twice", fname, line);
            else
                now->xdimname = strdup(token);
        } else if (strcasecmp(token, "YDIMNAME") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: YDIMNAME not specified", fname, line);
            else if (now->ydimname != NULL)
                enkf_quit("%s, l.%d: YDIMNAME specified twice", fname, line);
            else
                now->ydimname = strdup(token);
        } else if (strcasecmp(token, "ZDIMNAME") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: ZDIMNAME not specified", fname, line);
            else if (now->zdimname != NULL)
                enkf_quit("%s, l.%d: ZDIMNAME specified twice", fname, line);
            else
                now->zdimname = strdup(token);
        } else if (strcasecmp(token, "XVARNAME") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: XVARNAME not specified", fname, line);
            else if (now->xvarname != NULL)
                enkf_quit("%s, l.%d: XVARNAME specified twice", fname, line);
            else
                now->xvarname = strdup(token);
        } else if (strcasecmp(token, "YVARNAME") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: YVARNAME not specified", fname, line);
            else if (now->yvarname != NULL)
                enkf_quit("%s, l.%d: YVARNAME specified twice", fname, line);
            else
                now->yvarname = strdup(token);
        } else if (strcasecmp(token, "ZVARNAME") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: ZVARNAME not specified", fname, line);
            else if (now->zvarname != NULL)
                enkf_quit("%s, l.%d: ZVARNAME specified twice", fname, line);
            else
                now->zvarname = strdup(token);
        } else if (strcasecmp(token, "DEPTHVARNAME") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: DEPTHVARNAME not specified", fname, line);
            else if (now->depthvarname != NULL)
                enkf_quit("%s, l.%d: DEPTHVARNAME specified twice", fname, line);
            else
                now->depthvarname = strdup(token);
        } else if (now->levelvarnameentry != NULL && strcasecmp(token, now->levelvarnameentry) == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_printf("%s, l.%d: \"%s\" not specified", fname, line, now->levelvarnameentry);
            else if (now->levelvarname != NULL)
                enkf_quit("%s, l.%d: \"%s\" specified twice", fname, line, now->levelvarnameentry);
            else
                now->levelvarname = strdup(token);
        }
    }

    fclose(f);

    for (i = 0; i < *ngrid; ++i) {
        gridprm* now = &(*prm)[i];

        if (now->vtype == NULL)
            enkf_quit("%s: VTYPE not specified for grid \"%s\"", fname, now->name);
        if (now->fname == NULL)
            enkf_quit("%s: DATA not specified for grid \"%s\"", fname, now->name);
        if (now->xdimname == NULL)
            enkf_quit("%s: XDIMNAME not specified for grid \"%s\"", fname, now->name);
        if (now->ydimname == NULL)
            enkf_quit("%s: YDIMNAME not specified for grid \"%s\"", fname, now->name);
        if (now->zdimname == NULL)
            enkf_quit("%s: ZDIMNAME not specified for grid \"%s\"", fname, now->name);
        if (now->xvarname == NULL)
            enkf_quit("%s: XVARNAME not specified for grid \"%s\"", fname, now->name);
        if (now->yvarname == NULL)
            enkf_quit("%s: YVARNAME not specified for grid \"%s\"", fname, now->name);
        if (now->zvarname == NULL)
            enkf_quit("%s: ZVARNAME not specified for grid \"%s\"", fname, now->name);
    }
}

/**
 */
void gridprm_destroy(int ngrid, gridprm prm[])
{
    int i;

    for (i = 0; i < ngrid; ++i) {
        gridprm* now = &prm[i];

        free(now->name);
        free(now->fname);
        free(now->xdimname);
        free(now->ydimname);
        free(now->zdimname);
        free(now->xvarname);
        free(now->yvarname);
        free(now->zvarname);
        free(now->depthvarname);
        if (now->levelvarnameentry != NULL) {
            free(now->levelvarname);
            free(now->levelvarnameentry);
        }
    }
    free(prm);
}

/**
 */
void gridprm_print(gridprm* prm, char offset[])
{
    enkf_printf("%sgrid prm info:\n", offset);
    enkf_printf("%s  NAME = \"%s\"\n", offset, prm->name);
    enkf_printf("%s  FILE = \"%s\"\n", offset, prm->fname);
    enkf_printf("%s  XDIMNAME = \"%s\"\n", offset, prm->xdimname);
    enkf_printf("%s  YDIMNAME = \"%s\"\n", offset, prm->ydimname);
    enkf_printf("%s  ZDIMNAME = \"%s\"\n", offset, prm->zdimname);
    enkf_printf("%s  XVARNAME = \"%s\"\n", offset, prm->xvarname);
    enkf_printf("%s  YVARNAME = \"%s\"\n", offset, prm->yvarname);
    enkf_printf("%s  ZVARNAME = \"%s\"\n", offset, prm->zvarname);
    enkf_printf("%s  DEPTHVARNAME = \"%s\"\n", offset, prm->depthvarname);
    if (prm->levelvarnameentry != NULL && prm->levelvarname != NULL)
        enkf_printf("%s  %s = \"%s\"\n", offset, prm->levelvarnameentry, prm->levelvarname);
}

/**
 */
int gridprm_getvtype(gridprm* prm)
{
    int i;

    for (i = 0; i < ngridvtypeentries; ++i) {
        if (strcasecmp(prm->vtype, allgridvtypeentries[i].vtype_tag) == 0)
            return allgridvtypeentries[i].vtype;
    }

    return GRIDVTYPE_NONE;
}
