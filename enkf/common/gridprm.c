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
#include <math.h>
#include "definitions.h"
#include "grid.h"
#include "gridprm.h"
#include "utils.h"

#define NINC 10

typedef struct {
    char* vtype_tag;
    int vtype;
    char* levelvarnameentry;
} gridvtype_entry;

gridvtype_entry allgridvtypeentries[] = {
    {"NONE", GRIDVTYPE_NONE, "MASKVARNAME"},
    {"Z", GRIDVTYPE_Z, "NUMLEVELSVARNAME"},
    {"SIGMA", GRIDVTYPE_SIGMA, "MASKVARNAME"},
    {"HYBRID", GRIDVTYPE_HYBRID, "MASKVARNAME"}
};

int ngridvtypeentries = sizeof(allgridvtypeentries) / sizeof(gridvtype_entry);

/**
 */
void gridprm_create(char* fname, int* ngrid, gridprm** prm)
{
    gridprm* now = NULL;
    FILE* f = NULL;
    int ignore = 0;
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
            now->prmfname = strdup(fname);
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: NAME not specified", fname, line);
            else
                now->name = strdup(token);
            ignore = 0;
            if ((token = strtok(NULL, seps)) != NULL) {
                if (strcasecmp("PREP", token) != 0 && strcasecmp("CALC", token) != 0)
                    enkf_quit("%s, l.%d: grid qualifier = \"%s\"; must be either \"PREP\" or \"CALC\"\n", fname, line, token);
#if defined(ENKF_PREP)
                if (strcasecmp("PREP", token) != 0) {
                    (*ngrid)--;
                    now = NULL;
                    ignore = 1;
                }
#elif defined(ENKF_CALC) || defined(ENKF_UPDATE)
                if (strcasecmp("CALC", token) != 0) {
                    (*ngrid)--;
                    now = NULL;
                    ignore = 1;
                }
#endif
            }
            continue;
        } else if (now == NULL) {
            if (ignore)
                continue;
            else
                enkf_quit("%s, l.%d: entry NAME is required to start grid initialisation", fname, line);
        }

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
                enkf_quit("%s, l %d: VTYPE \"%s\" is unknown", fname, line, token);
            if (allgridvtypeentries[i].levelvarnameentry != NULL)
                now->levelvarnameentry = strdup(allgridvtypeentries[i].levelvarnameentry);
        } else if (strcasecmp(token, "DATA") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: DATA not specified", fname, line);
            else if (now->fname != NULL)
                enkf_quit("%s, l.%d: DATA specified twice", fname, line);
            else
                now->fname = strdup(token);
        } else if (strcasecmp(token, "HGRIDFROM") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: HGRIDFROM not specified", fname, line);
            else if (now->aliasname != NULL)
                enkf_quit("%s, l.%d: HGRIDFROM specified twice", fname, line);
            else
                now->aliasname = strdup(token);
        } else if (strcasecmp(token, "GEOGRAPHIC") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: GEOGRAPHIC not specified", fname, line);
            if (!str2bool(token, &now->geographic))
                enkf_quit("%s, l.%d: could not convert \"%s\" to boolean", fname, line, token);
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
            if (now->vtype == NULL)
                enkf_quit("%s, l.%d: VTYPE must be set first", fname, line);
            if (strcasecmp(now->vtype, "Z") != 0)
                enkf_quit("%s, l.%d: VTYPE must be set to \"Z\" for entry ZVARNAME (currently \"%s\")", fname, line, now->vtype);
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: ZVARNAME not specified", fname, line);
            else if (now->zvarname != NULL)
                enkf_quit("%s, l.%d: ZVARNAME specified twice", fname, line);
            else
                now->zvarname = strdup(token);
        } else if (strcasecmp(token, "ZCVARNAME") == 0) {
            if (now->vtype == NULL)
                enkf_quit("%s, l.%d: VTYPE must be set first", fname, line);
            if (strcasecmp(now->vtype, "Z") != 0)
                enkf_quit("%s, l.%d: VTYPE must be set to \"Z\" for entry ZVARNAME (currently \"%s\")", fname, line, now->vtype);
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: ZCVARNAME not specified", fname, line);
            else if (now->zcvarname != NULL)
                enkf_quit("%s, l.%d: ZCVARNAME specified twice", fname, line);
            else
                now->zcvarname = strdup(token);
        } else if (strcasecmp(token, "CVARNAME") == 0) {
            if (now->vtype == NULL)
                enkf_quit("%s, l.%d: VTYPE must be set first", fname, line);
            if (strcasecmp(now->vtype, "SIGMA") != 0)
                enkf_quit("%s, l.%d: VTYPE must be set to \"SIGMA\" for entry CVARNAME (currently \"%s\")", fname, line, now->vtype);
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: CVARNAME not specified", fname, line);
            else if (now->cvarname != NULL)
                enkf_quit("%s, l.%d: CVARNAME specified twice", fname, line);
            else
                now->cvarname = strdup(token);
        } else if (strcasecmp(token, "CCVARNAME") == 0) {
            if (now->vtype == NULL)
                enkf_quit("%s, l.%d: VTYPE must be set first", fname, line);
            if (strcasecmp(now->vtype, "SIGMA") != 0)
                enkf_quit("%s, l.%d: VTYPE must be set to \"SIGMA\" for entry CCVARNAME (currently \"%s\")", fname, line, now->vtype);
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: CCVARNAME not specified", fname, line);
            else if (now->ccvarname != NULL)
                enkf_quit("%s, l.%d: CCVARNAME specified twice", fname, line);
            else
                now->ccvarname = strdup(token);
        } else if (strcasecmp(token, "SVARNAME") == 0) {
            if (now->vtype == NULL)
                enkf_quit("%s, l.%d: VTYPE must be set first", fname, line);
            if (strcasecmp(now->vtype, "SIGMA") != 0)
                enkf_quit("%s, l.%d: VTYPE must be set to \"SIGMA\" for entry SVARNAME (currently \"%s\")", fname, line, now->vtype);
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: SVARNAME not specified", fname, line);
            else if (now->svarname != NULL)
                enkf_quit("%s, l.%d: SVARNAME specified twice", fname, line);
            else
                now->svarname = strdup(token);
        } else if (strcasecmp(token, "SCVARNAME") == 0) {
            if (now->vtype == NULL)
                enkf_quit("%s, l.%d: VTYPE must be set first", fname, line);
            if (strcasecmp(now->vtype, "SIGMA") != 0)
                enkf_quit("%s, l.%d: VTYPE must be set to \"SIGMA\" for entry SCVARNAME (currently \"%s\")", fname, line, now->vtype);
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: SCVARNAME not specified", fname, line);
            else if (now->scvarname != NULL)
                enkf_quit("%s, l.%d: SCVARNAME specified twice", fname, line);
            else
                now->scvarname = strdup(token);
        } else if (strcasecmp(token, "HCVARNAME") == 0) {
            if (now->vtype == NULL)
                enkf_quit("%s, l.%d: VTYPE must be set first", fname, line);
            if (strcasecmp(now->vtype, "SIGMA") != 0)
                enkf_quit("%s, l.%d: VTYPE must be set to \"SIGMA\" (currently \"%s\")", fname, line, now->vtype);
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: HCVARNAME not specified", fname, line);
            else if (now->hcvarname != NULL)
                enkf_quit("%s, l.%d: HCVARNAME specified twice", fname, line);
            else
                now->hcvarname = strdup(token);
        } else if (strcasecmp(token, "DEPTHVARNAME") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: DEPTHVARNAME not specified", fname, line);
            else if (now->depthvarname != NULL)
                enkf_quit("%s, l.%d: DEPTHVARNAME specified twice", fname, line);
            else
                now->depthvarname = strdup(token);
        } else if (strcasecmp(token, "STRIDE") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: STRIDE not specified", fname, line);
            if (now->stride != 0)
                enkf_quit("%s, l.%d: STRIDE specified twice", fname, line);
            if (!str2int(token, &now->stride))
                enkf_quit("%s, l.%d: could not convert \"%s\" to int", fname, line, token);
        } else if (strcasecmp(token, "SOBSTRIDE") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: SOBSTRIDE not specified", fname, line);
            if (now->sob_stride != 0)
                enkf_quit("%s, l.%d: SOBSTRIDE specified twice", fname, line);
            if (!str2int(token, &now->sob_stride))
                enkf_quit("%s, l.%d: could not convert SOBSTRIDE value", fname, line);
        } else if (now->levelvarnameentry != NULL && strcasecmp(token, now->levelvarnameentry) == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_printf("%s, l.%d: \"%s\" not specified", fname, line, now->levelvarnameentry);
            else if (now->levelvarname != NULL)
                enkf_quit("%s, l.%d: \"%s\" specified twice", fname, line, now->levelvarnameentry);
            else
                now->levelvarname = strdup(token);
        } else if (strncasecmp(token, "VDIR", 4) == 0) {
            if (now->vdirection != NULL)
                enkf_quit("%s, l.%d: VDIR specified twice", fname, line);
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: VDIR not specified", fname, line);
            if (strncasecmp(token, "TOSURF", 6) != 0 && strncasecmp(token, "FROMSURF", 8) != 0)
                enkf_quit("%s, l.%d: unknown entry for VDIR, fname, line");
            now->vdirection = strdup(token);
        } else if (strcasecmp(token, "AVARNAME") == 0) {
            if (now->vtype == NULL)
                enkf_quit("%s, l.%d: VTYPE must be set first", fname, line);
            if (strcasecmp(now->vtype, "HYBRID") != 0)
                enkf_quit("%s, l.%d: VTYPE must be set to \"HYBRID\" (currently \"%s\")", fname, line, now->vtype);
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: AVARNAME not specified", fname, line);
            else if (now->avarname != NULL)
                enkf_quit("%s, l.%d: AVARNAME specified twice", fname, line);
            else
                now->avarname = strdup(token);
        } else if (strcasecmp(token, "BVARNAME") == 0) {
            if (now->vtype == NULL)
                enkf_quit("%s, l.%d: VTYPE must be set first", fname, line);
            if (strcasecmp(now->vtype, "HYBRID") != 0)
                enkf_quit("%s, l.%d: VTYPE must be set to \"HYBRID\" (currently \"%s\")", fname, line, now->vtype);
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: BVARNAME not specified", fname, line);
            else if (now->bvarname != NULL)
                enkf_quit("%s, l.%d: BVARNAME specified twice", fname, line);
            else
                now->bvarname = strdup(token);
        } else if (strcasecmp(token, "ACVARNAME") == 0) {
            if (now->vtype == NULL)
                enkf_quit("%s, l.%d: VTYPE must be set first", fname, line);
            if (strcasecmp(now->vtype, "HYBRID") != 0)
                enkf_quit("%s, l.%d: VTYPE must be set to \"HYBRID\" (currently \"%s\")", fname, line, now->vtype);
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: ACVARNAME not specified", fname, line);
            else if (now->acvarname != NULL)
                enkf_quit("%s, l.%d: ACVARNAME specified twice", fname, line);
            else
                now->acvarname = strdup(token);
        } else if (strcasecmp(token, "BCVARNAME") == 0) {
            if (now->vtype == NULL)
                enkf_quit("%s, l.%d: VTYPE must be set first", fname, line);
            if (strcasecmp(now->vtype, "HYBRID") != 0)
                enkf_quit("%s, l.%d: VTYPE must be set to \"HYBRID\" (currently \"%s\")", fname, line, now->vtype);
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: BCVARNAME not specified", fname, line);
            else if (now->bcvarname != NULL)
                enkf_quit("%s, l.%d: BCVARNAME specified twice", fname, line);
            else
                now->bcvarname = strdup(token);
        } else if (strcasecmp(token, "P1VARNAME") == 0) {
            if (now->vtype == NULL)
                enkf_quit("%s, l.%d: VTYPE must be set first", fname, line);
            if (strcasecmp(now->vtype, "HYBRID") != 0)
                enkf_quit("%s, l.%d: VTYPE must be set to \"HYBRID\" (currently \"%s\")", fname, line, now->vtype);
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: P1VARNAME not specified", fname, line);
            else if (now->p1varname != NULL)
                enkf_quit("%s, l.%d: P1VARNAME specified twice", fname, line);
            else
                now->p1varname = strdup(token);
        } else if (strcasecmp(token, "P2VARNAME") == 0) {
            if (now->vtype == NULL)
                enkf_quit("%s, l.%d: VTYPE must be set first", fname, line);
            if (strcasecmp(now->vtype, "HYBRID") != 0)
                enkf_quit("%s, l.%d: VTYPE must be set to \"HYBRID\" (currently \"%s\")", fname, line, now->vtype);
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: P2VARNAME not specified", fname, line);
            else if (now->p2varname != NULL)
                enkf_quit("%s, l.%d: P2VARNAME specified twice", fname, line);
            else
                now->p2varname = strdup(token);
        } else if (strcasecmp(token, "ZSTATINTS") == 0) {
            char zseps[] = " =\t\n[](){}";

            while ((token = strtok(NULL, zseps)) != NULL) {
                if (now->nzints % NINC == 0)
                    now->zints = realloc(now->zints, (now->nzints + NINC) * sizeof(zint));
                if (!str2double(token, &now->zints[now->nzints].z1))
                    enkf_quit("%s, l.%d: could not convert \"%s\" to double", fname, line, token);
                if ((token = strtok(NULL, zseps)) == NULL)
                    enkf_quit("%s, l.%d: maximal depth/height for an interval not specified", fname, line);
                if (!str2double(token, &now->zints[now->nzints].z2))
                    enkf_quit("%s, l.%d: could not convert \"%s\" to double", fname, line, token);
                now->nzints++;
            }
            /*
             * a temporal setting, to indicate that an empty range has been
             * entered (as opposed to no entry)
             */
            if (now->nzints == 0)
                now->nzints = -1;
        } else if (strcasecmp(token, "DOMAIN") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: DOMAIN not specified", fname, line);
            else if (now->p1varname != NULL)
                enkf_quit("%s, l.%d: DOMAIN specified twice", fname, line);
            else
                now->domainname = strdup(token);
        } else
            enkf_quit("%s, l.%d: unexpected token \"%s\"", fname, line, token);
    }
    fclose(f);

    for (i = 0; i < *ngrid; ++i) {
        gridprm* gprm = &(*prm)[i];

        if (gprm->vtype == NULL)
            enkf_quit("%s: VTYPE not specified for grid \"%s\"", fname, gprm->name);
        if (gprm->fname == NULL)
            enkf_quit("%s: DATA not specified for grid \"%s\"", fname, gprm->name);
        if (gprm->aliasname != NULL && (gprm->xvarname != NULL || gprm->yvarname != NULL))
            enkf_quit("%s: %s: either HGRIDFROM or XVARNAME and YVARNAME should be specified", fname, gprm->name);
        if (gprm->domainname == NULL)
            gprm->domainname = strdup("Default");
        if (strcasecmp(gprm->vtype, "NONE") == 0);
        else if (strcasecmp(gprm->vtype, "Z") == 0) {
            if (gprm->zvarname == NULL)
                enkf_quit("%s: %s: ZVARNAME must be specified for Z grids", fname, gprm->name);
        } else if (strcasecmp(gprm->vtype, "SIGMA") == 0) {
            if (gprm->cvarname == NULL)
                enkf_quit("%s: %s: CVARNAME must be specified for sigma grids", fname, gprm->name);
        } else if (strcasecmp(gprm->vtype, "HYBRID") == 0) {
            if (gprm->avarname == NULL)
                enkf_quit("%s: %s: AVARNAME must be specified for hybrid grids", fname, gprm->name);
            if (gprm->bvarname == NULL)
                enkf_quit("%s: %s: BVARNAME must be specified for hybrid grids", fname, gprm->name);
            if (gprm->p1varname == NULL)
                enkf_quit("%s: %s: P1VARNAME must be specified for hybrid grids", fname, gprm->name);
            if (gprm->p2varname == NULL)
                enkf_quit("%s: %s: P2VARNAME must be specified for hybrid grids", fname, gprm->name);
        } else
            enkf_quit("vertical type \"%s\" specified for grid \"%s\" is unknown", gprm->vtype, gprm->name);
        if (gprm->aliasname == NULL && gprm->xvarname == NULL)
            enkf_quit("%s: XVARNAME not specified for grid \"%s\"", fname, gprm->name);
        if (gprm->aliasname == NULL && gprm->yvarname == NULL)
            enkf_quit("%s: YVARNAME not specified for grid \"%s\"", fname, gprm->name);
        if (strcasecmp(gprm->vtype, "NONE") != 0) {
            if (gprm->vdirection == NULL)
                gprm->vdirection = strdup("FROMSURF");
            if (gprm->nzints == 0) {
                gprm->nzints = 3;
                gprm->zints = malloc(gprm->nzints * sizeof(zint));
                gprm->zints[0].z1 = 0.0;
                gprm->zints[0].z2 = DEPTH_SHALLOW;
                gprm->zints[1].z1 = DEPTH_SHALLOW;
                gprm->zints[1].z2 = DEPTH_DEEP;
                gprm->zints[2].z1 = DEPTH_DEEP;
                gprm->zints[2].z2 = DEPTH_MAX;
            } else if (gprm->nzints < 0)
                gprm->nzints = 0;
        } else
            gprm->zints = 0;
    }
}

/**
 */
void gridprm_destroy(int ngrid, gridprm prm[])
{
    int i;

    for (i = 0; i < ngrid; ++i) {
        gridprm* now = &prm[i];

        free(now->prmfname);
        free(now->name);
        free(now->fname);
        if (now->aliasname != NULL)
            free(now->aliasname);
        else {
            free(now->xvarname);
            free(now->yvarname);
        }
        if (now->vtype != NULL)
            free(now->vtype);
        if (now->zvarname != NULL)
            free(now->zvarname);
        if (now->zcvarname != NULL)
            free(now->zcvarname);
        if (now->hcvarname != NULL)
            free(now->hcvarname);
        if (now->depthvarname != NULL)
            free(now->depthvarname);
        if (now->vdirection != NULL)
            free(now->vdirection);
        if (now->avarname != NULL)
            free(now->avarname);
        if (now->bvarname != NULL)
            free(now->bvarname);
        if (now->p1varname != NULL)
            free(now->p1varname);
        if (now->p2varname != NULL)
            free(now->p2varname);
        if (now->levelvarnameentry != NULL) {
            free(now->levelvarname);
            free(now->levelvarnameentry);
        }
        if (now->nzints > 0)
            free(now->zints);
        free(now->domainname);
    }
    free(prm);
}

/**
 */
void gridprm_print(gridprm* prm, char offset[])
{
    int i;

    enkf_printf("%sgrid prm info:\n", offset);
    enkf_printf("%s  NAME = \"%s\"\n", offset, prm->name);
    enkf_printf("%s  DOMAIN = %s\n", offset, prm->domainname);
    enkf_printf("%s  DATA = \"%s\"\n", offset, prm->fname);
    enkf_printf("%s  VTYPE = \"%s\"\n", offset, prm->vtype);
    enkf_printf("%s  VDIR = %s\n", offset, prm->vdirection);
    enkf_printf("%s  GEOGRAPHIC = %s\n", offset, (prm->geographic) ? "yes" : "no");
    if (strcasecmp(prm->vtype, "NONE") == 0);
    else if (strcasecmp(prm->vtype, "Z") == 0) {
        enkf_printf("%s  ZVARNAME = \"%s\"\n", offset, prm->zvarname);
        if (prm->zcvarname != NULL)
            enkf_printf("%s  ZCVARNAME = \"%s\"\n", offset, prm->zcvarname);
        else
            enkf_printf("%s  ZCVARNAME = <none>\n", offset);
    } else if (strcasecmp(prm->vtype, "SIGMA") == 0) {
        enkf_printf("%s  CVARNAME = \"%s\"\n", offset, prm->cvarname);
        if (prm->ccvarname != NULL)
            enkf_printf("%s  CCVARNAME = \"%s\"\n", offset, prm->ccvarname);
        else
            enkf_printf("%s  CCVARNAME = <none>\n", offset);
        if (prm->svarname != NULL)
            enkf_printf("%s  SVARNAME = \"%s\"\n", offset, prm->svarname);
        else
            enkf_printf("%s  SVARNAME = <none>\n", offset);
        if (prm->scvarname != NULL)
            enkf_printf("%s  SCVARNAME = \"%s\"\n", offset, prm->scvarname);
        else
            enkf_printf("%s  SCVARNAME = <none>\n", offset);
        if (prm->hcvarname != 0)
            enkf_printf("%s  HCVARNAME = \"%s\"\n", offset, prm->hcvarname);
        else
            enkf_printf("%s  HCVARNAME = <none> (meaning hc = 0)\n", offset);
    } else if (strcasecmp(prm->vtype, "HYBRID") == 0) {
        enkf_printf("%s    AVARNAME = \"%s\"\n", offset, prm->avarname);
        enkf_printf("%s    BVARNAME = \"%s\"\n", offset, prm->bvarname);
        enkf_printf("%s    P1VARNAME = \"%s\"\n", offset, prm->p1varname);
        enkf_printf("%s    P2VARNAME = \"%s\"\n", offset, prm->p2varname);
    } else
        enkf_quit("programming error");
    if (prm->depthvarname != NULL)
        enkf_printf("%s  DEPTHVARNAME = \"%s\"\n", offset, prm->depthvarname);
    if (prm->levelvarnameentry != NULL && prm->levelvarname != NULL)
        enkf_printf("%s  %s = \"%s\"\n", offset, prm->levelvarnameentry, prm->levelvarname);
    if (prm->vdirection != NULL)
        enkf_printf("%s  VDIR = \"%s\"\n", offset, prm->vdirection);
    if (prm->xvarname != NULL) {
        enkf_printf("%s  XVARNAME = \"%s\"\n", offset, prm->xvarname);
        enkf_printf("%s  YVARNAME = \"%s\"\n", offset, prm->yvarname);
    } else if (prm->aliasname != NULL)
        enkf_printf("%s  HGRIDFROM = \"%s\"\n", offset, prm->aliasname);
    if (prm->stride != 0)
        enkf_printf("%s  STRIDE = %d\n", offset, prm->stride);
    if (prm->sob_stride != 0)
        enkf_printf("%s  SOBSTRIDE = %d\n", offset, prm->sob_stride);
    if (prm->nzints != 0) {
        enkf_printf("%s  ZSTATINTS = ", offset);
        for (i = 0; i < prm->nzints; ++i)
            enkf_printf("[%.0f %.0f] ", prm->zints[i].z1, prm->zints[i].z2);
        enkf_printf("\n");
    }
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

    return GRIDVTYPE_UNDEFINED;
}
