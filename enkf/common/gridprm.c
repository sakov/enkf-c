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
#include "gridprm.h"
#include "utils.h"

/**
 */
void gridprm_create(char* fname, int* ngrid, gridprm** prm, char* levelvarnameentry)
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
            now->levelvarnameentry = strdup(levelvarnameentry);
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: NAME not specified", fname, line);
            else
                now->name = strdup(token);
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
        } else if (strcasecmp(token, levelvarnameentry) == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: \"%s\" not specified", fname, line, levelvarnameentry);
            else if (now->levelvarname != NULL)
                enkf_quit("%s, l.%d: \"%s\" specified twice", fname, line, levelvarnameentry);
            else
                now->levelvarname = strdup(token);
        }
    }

    fclose(f);

    for (i = 0; i < *ngrid; ++i) {
        gridprm* now = &(*prm)[i];

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
        if (now->depthvarname == NULL)
            if (now->levelvarname == NULL)
                enkf_quit("%s: \"%s\" not specified for grid \"%s\"", fname, now->levelvarnameentry, now->name);
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
        free(now->levelvarname);
        free(now->levelvarnameentry);
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
    enkf_printf("%s  %s = \"%s\"\n", offset, prm->levelvarnameentry, prm->levelvarname);
}
