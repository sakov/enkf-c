/******************************************************************************
 *
 * File:        obstypes.c        
 *
 * Created:     18/08/2014
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description: Handles "obstype" structure, which provides interface between
 *              observations and model.
 *
 * Revisions:
 *
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include "kdtree.h"
#include "definitions.h"
#include "utils.h"
#include "grid.h"
#include "model.h"
#include "enkfprm.h"
#include "obstypes.h"

#define NINC 10

/**
 */
static void obstype_new(obstype* type, int i, char* name)
{
    type->id = i;
    type->issurface = -1;
    type->name = strdup(name);
    type->nvar = 0;
    type->varnames = NULL;
    type->issurface = -1;
    type->offset_fname = NULL;
    type->offset_varname = NULL;
    type->mld_varname = NULL;
    type->mld_threshold = NAN;
    type->hfunction = NULL;
    type->allowed_min = -DBL_MAX;
    type->allowed_max = DBL_MAX;
    type->isasync = 0;
    type->async_tstep = NAN;
    type->async_centred = 1;    /* interval 0 is centred at assimilation
                                 * time, suitable for instantaneous fields */
    type->nlocrad = 0;
    type->locrad = NULL;
    type->locweight = NULL;
    type->rfactor = 1.0;
    type->nlobsmax = -1;
    type->vid = -1;
    type->gridid = -1;
    type->xmin = -DBL_MAX;
    type->xmax = DBL_MAX;
    type->ymin = -DBL_MAX;
    type->ymax = DBL_MAX;
    type->zmin = -DBL_MAX;
    type->zmax = DBL_MAX;
    type->nobs = 0;
    type->ngood = 0;
    type->noutside_grid = 0;
    type->noutside_obsdomain = 0;
    type->noutside_obswindow = 0;
    type->nland = 0;
    type->nshallow = 0;
    type->nbadbatch = 0;
    type->nrange = 0;
    type->nthinned = 0;
    type->nsubgrid = 0;
    type->nmodified = 0;
    type->windowmin = NAN;
    type->windowmax = NAN;
    type->date_min = DBL_MAX;
    type->date_max = -DBL_MAX;
    type->ndomains = 0;
    type->domainnames = NULL;
}

/**
 */
static void obstype_check(obstype* type)
{
    assert(type->name != NULL);
    if (type->issurface == -1)
        enkf_quit("\"%s\": ISSURFACE not specified\n", type->name);
    if (type->nvar == 0)
        enkf_quit("\"%s\": VAR not specified\n", type->name);
    if (type->rfactor <= 0)
        enkf_quit("\"%s\": RFACTOR = %f\n", type->name, type->rfactor);
    if (type->nlobsmax < 0)
        enkf_quit("\"%s\": NLOBSMAX = %d\n", type->name, type->nlobsmax);
    if (type->hfunction == NULL)
        enkf_quit("\"%s\": HFUNCTION not specified\n", type->name);
}

/**
 */
static void obstype_print(obstype* type)
{
    int i;

    enkf_printf("    NAME = %s\n", type->name);
    if (type->ndomains > 0) {
        enkf_printf("    DOMAINS =");
        for (i = 0; i < type->ndomains; ++i)
            enkf_printf(" %s\n", type->domainnames[i]);
        enkf_printf("\n");
    }
    enkf_printf("    ISSURFACE = %d\n", (type->issurface) ? 1 : 0);
    enkf_printf("      VAR =");
    for (i = 0; i < type->nvar; ++i)
        enkf_printf(" %s", type->varnames[i]);
    enkf_printf("\n");
    enkf_printf("      ID = %d\n", type->id);
    if (type->offset_fname != NULL)
        enkf_printf("      OFFSET = %s %s\n", type->offset_fname, type->offset_varname);
    enkf_printf("      HFUNCTION = %s\n", type->hfunction);
    enkf_printf("      ALLOWED MIN = %.3g\n", type->allowed_min);
    enkf_printf("      ALLOWED MAX = %.3g\n", type->allowed_max);
    enkf_printf("      ASYNCHRONOUS = %s", (type->isasync) ? "yes" : "no");
    if (type->isasync)
        enkf_printf(", DT = %.3f (%s)\n", type->async_tstep, (type->async_centred) ? "centre-aligned" : "endpoint-aligned");
    else
        enkf_printf("\n");
    enkf_printf("      LOCRAD  =");
    for (i = 0; i < type->nlocrad; ++i)
        enkf_printf(" %.3g", type->locrad[i]);
    enkf_printf("\n");
    enkf_printf("      LOCWEIGHT = ");
    for (i = 0; i < type->nlocrad; ++i)
        enkf_printf(" %.3g", type->locweight[i]);
    enkf_printf("\n");
    enkf_printf("      RFACTOR = %.3g\n", type->rfactor);
    if (type->nlobsmax != INT_MAX)
        enkf_printf("      NLOBSMAX = %d\n", type->nlobsmax);
    if (type->xmin > -DBL_MAX || type->xmax < DBL_MAX || type->ymin > -DBL_MAX || type->ymax < DBL_MAX || type->zmin > -DBL_MAX || type->zmax < DBL_MAX)
        enkf_printf("      SPATIAL DOMAIN = %.3g %.3g %.3g %.3g %.3g %.3g\n", type->xmin, type->xmax, type->ymin, type->ymax, type->zmin, type->zmax);
    if (isfinite(type->windowmin)) {
        enkf_printf("      WINDOWMIN = %.3f", type->windowmin);
        enkf_printf("      WINDOWMAX = %.3f", type->windowmax);
    }
}

/**
 */
void obstypes_read(enkfprm* prm, char fname[], int* n, obstype** types)
{
    FILE* f = NULL;
    char buf[MAXSTRLEN];
    int line;
    obstype* now = NULL;
    int i;

    assert(*n == 0);

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
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: NAME not specified", fname, line);
            if (obstype_getid(*n, *types, token, 0) >= 0)
                enkf_quit("%s: l.%d: type \"%s\" already specified", fname, line, token);

            *types = realloc(*types, (*n + 1) * sizeof(obstype));
            now = &(*types)[*n];
            obstype_new(now, *n, token);
            (*n)++;
            continue;
        }

        if (now == NULL)
            enkf_quit("%s, l.%d: NAME not specified", fname, line);

        if (strcasecmp(token, "ISSURFACE") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: ISSURFACE not specified", fname, line);
            if (now->issurface >= 0)
                enkf_quit("%s, l.%d: ISSURFACE already specified", fname, line);
            if (!str2bool(token, &now->issurface))
                enkf_quit("%s, l.%d: could not convert \"%s\" to boolean", fname, line, token);
        } else if (strcasecmp(token, "VAR") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: VAR not specified", fname, line);
            if (now->varnames != NULL)
                enkf_quit("%s, l.%d: VAR already specified", fname, line);
            while (token != NULL) {
                if (now->nvar % NINC == 0)
                    now->varnames = realloc(now->varnames, (now->nvar + NINC) * sizeof(char*));
                now->varnames[now->nvar] = strdup(token);
                now->nvar++;
                token = strtok(NULL, seps);
            }
        } else if (strcasecmp(token, "OFFSET") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: OFFSET file name not specified", fname, line);
            if (now->offset_fname != NULL)
                enkf_quit("%s, l.%d: OFFSET entry already specified", fname, line);
            now->offset_fname = strdup(token);
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: OFFSET variable name not specified", fname, line);
            now->offset_varname = strdup(token);
        } else if (strcasecmp(token, "MLD_VARNAME") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: MLD_VARNAME not specified", fname, line);
            if (now->mld_varname != NULL)
                enkf_quit("%s, l.%d: MLD_VARNAME entry already specified", fname, line);
            now->mld_varname = strdup(token);
        } else if (strcasecmp(token, "MLD_THRESH") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: MLD_THRESH not specified", fname, line);
            if (!isnan(now->mld_threshold))
                enkf_quit("%s, l.%d: MLD_THRESH already specified", fname, line);
            if (!str2double(token, &now->mld_threshold))
                enkf_quit("%s, l.%d: could not convert \"%s\" to double", fname, line, token);
        } else if (strcasecmp(token, "HFUNCTION") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: HFUNCTION not specified", fname, line);
            if (now->hfunction != NULL)
                enkf_quit("%s, l.%d: HFUNCTION already specified", fname, line);
            now->hfunction = strdup(token);
        } else if (strcasecmp(token, "MINVALUE") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: MINVALUE not specified", fname, line);
            if (!str2double(token, &now->allowed_min))
                enkf_quit("%s, l.%d: could not convert \"%s\" to double", fname, line, token);
        } else if (strcasecmp(token, "MAXVALUE") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: MAXVALUE not specified", fname, line);
            if (!str2double(token, &now->allowed_max))
                enkf_quit("%s, l.%d: could not convert \"%s\" to double", fname, line, token);
        } else if (strncasecmp(token, "ASYNC", 5) == 0) {
            now->isasync = 1;
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: ASYNC time interval not specified", fname, line);
            if (!str2double(token, &now->async_tstep))
                enkf_quit("%s, l.%d: could not convert \"%s\" to double", fname, line, token);
            if (now->async_tstep == 0.0)
                now->isasync = 0;
            else if (now->async_tstep < 0.0)
                enkf_quit("%s, l.%d: negative length of asynchronous time interval", fname, line);
            if ((token = strtok(NULL, seps)) != NULL) {
                if (token[0] == 'c' || token[0] == 'C')
                    now->async_centred = 1;
                else if (token[0] == 'e' || token[0] == 'E')
                    now->async_centred = 0;
                else
                    enkf_quit("%s, l.%d: the asynchronous intervals can be either \"c\" (centre-aligned) or \"e\" (endpoint-aligned)", fname, line);
            }
        } else if (strcasecmp(token, "LOCRAD") == 0) {
            int sid = 0;

            if (now->nlocrad > 0)
                now->locrad = malloc(sizeof(double) * now->nlocrad);
            while ((token = strtok(NULL, seps)) != NULL) {
                if (now->nlocrad == sid) {
                    now->locrad = realloc(now->locrad, sizeof(double) * (now->nlocrad + 1));
                    now->nlocrad++;
                }
                if (!str2double(token, &now->locrad[sid]))
                    enkf_quit("%s, l.%d: could not convert \"%s\" to double", fname, line, token);
                sid++;
            }
            if (now->nlocrad > sid)
                enkf_quit("%s, l.%d: LOCRAD not specified or its dimension does not match that of RFACTOR", fname, line);
        } else if (strcasecmp(token, "LOCWEIGHT") == 0) {
            int sid = 0;

            if (now->nlocrad > 0) {
                if (now->locweight == NULL)
                    now->locweight = malloc(sizeof(double) * now->nlocrad);
            } else
                enkf_quit("%s, l.%d: LOCRAD must be entered before LOCWEIGHT", fname, line);
            while ((token = strtok(NULL, seps)) != NULL) {
                if (now->nlocrad == sid) {
                    now->locweight = realloc(now->locweight, sizeof(double) * (now->nlocrad + 1));
                    now->nlocrad++;
                }
                if (!str2double(token, &now->locweight[sid]))
                    enkf_quit("%s, l.%d: could not convert \"%s\" to double", fname, line, token);
                sid++;
            }
            if (now->nlocrad > sid)
                enkf_quit("%s, l.%d: LOCWEIGHT entered but not specified or its dimension does not match that of LOCRAD", fname, line);
        } else if (strcasecmp(token, "RFACTOR") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: RFACTOR not specified", fname, line);
            if (!str2double(token, &now->rfactor))
                enkf_quit("%s, l.%d: could not convert \"%s\" to double", fname, line, token);
        } else if (strcasecmp(token, "NLOBSMAX") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: NLOBSMAX not specified", fname, line);
            if (!str2int(token, &now->nlobsmax))
                enkf_quit("%s, l.%d: could not convert \"%s\" to int", fname, line, token);
        } else if (strcasecmp(token, "XMIN") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: XMIN not specified", fname, line);
            if (!str2double(token, &now->xmin))
                enkf_quit("%s, l.%d: could not convert \"%s\" to double", fname, line, token);
        } else if (strcasecmp(token, "XMAX") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: XMAX not specified", fname, line);
            if (!str2double(token, &now->xmax))
                enkf_quit("%s, l.%d: could not convert \"%s\" to double", fname, line, token);
        } else if (strcasecmp(token, "YMIN") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: YMIN not specified", fname, line);
            if (!str2double(token, &now->ymin))
                enkf_quit("%s, l.%d: could not convert \"%s\" to double", fname, line, token);
        } else if (strcasecmp(token, "YMAX") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: YMAX not specified", fname, line);
            if (!str2double(token, &now->ymax))
                enkf_quit("%s, l.%d: could not convert \"%s\" to double", fname, line, token);
        } else if (strcasecmp(token, "ZMIN") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: ZMIN not specified", fname, line);
            if (!str2double(token, &now->zmin))
                enkf_quit("%s, l.%d: could not convert \"%s\" to double", fname, line, token);
        } else if (strcasecmp(token, "ZMAX") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: ZMAX not specified", fname, line);
            if (!str2double(token, &now->zmax))
                enkf_quit("%s, l.%d: could not convert \"%s\" to double", fname, line, token);
        } else if (strcasecmp(token, "WINDOWMIN") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: WINDOWMIN not specified", fname, line);
            else if (!str2double(token, &now->windowmin))
                enkf_quit("%s, l.%d: could convert WINDOWMIN entry", fname, line);
        } else if (strcasecmp(token, "WINDOWMAX") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: WINDOWMAX not specified", fname, line);
            else if (!str2double(token, &now->windowmax))
                enkf_quit("%s, l.%d: could convert WINDOWMAX entry", fname, line);
        } else if (strcasecmp(token, "DOMAINS") == 0) {
            if (now->ndomains != 0)
                enkf_quit("%s, l.%d: DOMAINS specified twice", fname, line);
            while ((token = strtok(NULL, seps)) != NULL) {
                if (now->ndomains % NINC == 0)
                    now->domainnames = realloc(now->domainnames, (now->ndomains + NINC) * sizeof(char*));
                now->domainnames[now->ndomains] = strdup(token);
                now->ndomains++;
            }
            if (now->ndomains == 0)
                enkf_quit("%s, l.%d: DOMAINS not specified", fname, line);
        } else
            enkf_quit("%s, l.%d: unknown token \"%s\"", fname, line, token);
    }
    fclose(f);

    for (i = 0; i < *n; ++i) {
        obstype* type = &(*types)[i];

        if (type->nlocrad == 0) {
            int j;

#if defined(ENKF_CALC)
            if (prm->nlocrad == 0 && !enkf_fstatsonly)
                enkf_quit("LOCRAD specified neither in %s or %s for obstype %s", prm->fname, fname, type->name);
#endif
            type->nlocrad = prm->nlocrad;
            type->locrad = malloc(type->nlocrad * sizeof(double));
            type->locweight = malloc(type->nlocrad * sizeof(double));
            for (j = 0; j < type->nlocrad; ++j) {
                type->locrad[j] = prm->locrad[j];
                type->locweight[j] = prm->locweight[j];
            }
        } else if (type->locweight == NULL) {
            if (type->nlocrad == 1) {
                type->locweight = malloc(sizeof(double));
                type->locweight[0] = 1.0;
            } else
                enkf_quit("%s: LOCWEIGHT not specified for multi-scale observation type \"%s\"", fname, type->name);
        } else {
            double sum = 0.0;
            int j;

            for (j = 0; j < type->nlocrad; ++j)
                sum += type->locweight[j];
            assert(sum > 0.0);
            for (j = 0; j < type->nlocrad; ++j)
                type->locweight[j] /= sum;
        }
        type->rfactor *= prm->rfactor_base;

        if (type->nlobsmax < 0)
            type->nlobsmax = prm->nlobsmax;
    }

    for (i = 0; i < *n; ++i) {
        obstype_check(&(*types)[i]);
        obstype_print(&(*types)[i]);
    }
}

/**
 */
void obstypes_describeprm(void)
{
    enkf_printf("\n");
    enkf_printf("  Observation types parameter file format:\n");
    enkf_printf("\n");
    enkf_printf("    NAME        = <name>\n");
    enkf_printf("  [ DOMAINS     = <domain name> ... ]\n");
    enkf_printf("    ISSURFACE   = {0 | 1}\n");
    enkf_printf("    VAR         = <model variable name> ...\n");
    enkf_printf("  [ OFFSET      = <file name> <variable name> ]    (none*)\n");
    enkf_printf("  [ MLD_VARNAME = <model varname> ]                (none*)\n");
    enkf_printf("  [ MLD_THRESH  = <threshold> ]                    (NaN*)\n");
    enkf_printf("    HFUNCTION   = <H function name>\n");
    enkf_printf("  [ ASYNC       = <time interval> [c*|e]]          (synchronous*)\n");
    enkf_printf("  [ LOCRAD      = <loc. radius in km> ... ]\n");
    enkf_printf("  [ LOCWEIGHT   = <loc. weight> ... ]              (# LOCRAD > 1)\n");
    enkf_printf("  [ RFACTOR     = <rfactor> ]                      (1*)\n");
    enkf_printf("  [ NLOBSMAX    = <max. number of local obs. ]\n");
    enkf_printf("  [ MINVALUE    = <minimal allowed value> ]        (-inf*)\n");
    enkf_printf("  [ MAXVALUE    = <maximal allowed value> ]        (+inf*)\n");
    enkf_printf("  [ XMIN        = <minimal allowed X coordinate> ] (-inf*)\n");
    enkf_printf("  [ XMAX        = <maximal allowed X coordinate> ] (+inf*)\n");
    enkf_printf("  [ YMIN        = <minimal allowed Y coordinate> ] (-inf*)\n");
    enkf_printf("  [ YMAX        = <maximal allowed Y coordinate> ] (+inf*)\n");
    enkf_printf("  [ ZMIN        = <minimal allowed Z coordinate> ] (-inf*)\n");
    enkf_printf("  [ ZMAX        = <maximal allowed Z coordinate> ] (+inf*)\n");
    enkf_printf("  [ WINDOWMIN   = <start of obs window in days from analysis> ] (-inf*)\n");
    enkf_printf("  [ WINDOWMAX   = <end of obs window in days from analysis> ]   (+inf*)\n");
    enkf_printf("\n");
    enkf_printf("  [ <more of the above blocks> ]\n");
    enkf_printf("\n");
    enkf_printf("  Notes:\n");
    enkf_printf("    1. { ... | ... | ... } denotes the list of possible choices\n");
    enkf_printf("    2. [ ... ] denotes an optional input\n");
    enkf_printf("    3. ( ... ) is a note\n");
    enkf_printf("    4. * denotes the default value\n");
    enkf_printf("    5. < ... > denotes a description of an entry\n");
    enkf_printf("\n");
}

/**
 */
void obstypes_set(int n, obstype* types, model* m)
{
    int i;

    for (i = 0; i < n; ++i) {
        obstype* ot = &types[i];
        int vid = model_getvarid(m, types[i].varnames[0], 1);
        int j;

        ot->vid = vid;
        ot->gridid = model_getvargridid(m, vid);
        ot->sob_stride = grid_getsobstride(model_getgridbyid(m, ot->gridid));
        if (ot->ndomains > 0)
            for (j = 0; j < ot->ndomains; ++j)
                if (model_getdomainid(m, ot->domainnames[j]) < 0)
                    enkf_quit("OBSTYPE = %s: no grid is associated with domain \"%s\"\n", ot->name, ot->domainnames[j]);

#if defined(ENKF_CALC)
        char tag[MAXSTRLEN];

        snprintf(tag, MAXSTRLEN, "%s:OFFSET", ot->name);
        if (ot->offset_fname != NULL) {
            if (ot->issurface || !is3d(ot->offset_fname, ot->offset_varname)) {
                float** v = NULL;
                int nx, ny;

                model_getvardims(m, vid, &nx, &ny, NULL);
                v = alloc2d(ny, nx, sizeof(float));
                readfield(ot->offset_fname, ot->offset_varname, 0, nx, ny, 1, v[0]);
                model_adddata(m, tag, vid, ALLOCTYPE_2D, v);
            } else {
                float*** v = NULL;
                int nx, ny, nz;

                model_getvardims(m, vid, &nx, &ny, &nz);
                v = alloc3d(nz, ny, nx, sizeof(float));
                read3dfield(ot->offset_fname, ot->offset_varname, nx, ny, nz, v[0][0]);
                model_adddata(m, tag, vid, ALLOCTYPE_3D, v);
            }
        }
#endif
    }
}

/**
 */
void obstypes_destroy(int n, obstype* types)
{
    int i, j;

    if (n == 0)
        return;

    for (i = 0; i < n; ++i) {
        obstype* ot = &types[i];

        free(ot->name);
        for (j = 0; j < ot->nvar; ++j)
            free(ot->varnames[j]);
        free(ot->varnames);
        free(ot->hfunction);
        if (ot->offset_fname != NULL) {
            free(ot->offset_fname);
            free(ot->offset_varname);
        }
        if (ot->mld_varname != NULL)
            free(ot->mld_varname);
        free(ot->locrad);
        free(ot->locweight);
        if (ot->ndomains > 0) {
            for (j = 0; j < ot->ndomains; ++j)
                free(ot->domainnames[j]);
            free(ot->domainnames);
        }
    }

    free(types);
}

/**
 */
int obstype_getid(int n, obstype types[], char* name, int hastosucceed)
{
    int i;

    for (i = 0; i < n; ++i)
        if (strcmp(types[i].name, name) == 0)
            return i;

    if (hastosucceed)
        enkf_quit("failed to identify observation type \"%s\"", name);

    return -1;
}

/**
 */
double obstype_calclcoeff(obstype* type, double dist)
{
    double sum = 0.0;
    int i;

    for (i = 0; i < type->nlocrad; ++i)
        if (dist <= type->locrad[i])
            sum += type->locweight[i] * taper_gc(dist / type->locrad[i]);
    return sum;
}

/**
 */
double obstype_getmaxlocrad(obstype* type)
{
    double max = 0.0;
    int i;

    for (i = 0; i < type->nlocrad; ++i)
        if (type->locrad[i] > max)
            max = type->locrad[i];
    return max;
}
