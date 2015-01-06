/******************************************************************************
 *
 * File:        obstypes.c        
 *
 * Created:     18/08/2014
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
#include <values.h>
#include <assert.h>
#include "string.h"
#include "nan.h"
#include "definitions.h"
#include "utils.h"
#include "obstypes.h"

/**
 */
int obstype_getid(int n, obstype types[], char* name)
{
    int i;

    for (i = 0; i < n; ++i)
        if (strcmp(types[i].name, name) == 0)
            return i;

    return -1;
}

/**
 */
static void obstype_new(obstype* type, int i, char* name)
{
    type->id = i;
    type->name = strdup(name);
    type->varname = NULL;
    type->varname2 = NULL;
    type->issurface = -1;
    type->offset_fname = NULL;
    type->offset_varname = NULL;
    type->hfunction = NULL;
    type->allowed_min = -DBL_MAX;
    type->allowed_max = DBL_MAX;
    type->isasync = 0;
    type->async_tstep = NaN;
    type->rfactor = 1.0;
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
    type->noutside = 0;
    type->nland = 0;
    type->nshallow = 0;
    type->nbadbatch = 0;
    type->nroundup = 0;
    type->nrange = 0;
    type->nmodified = 0;
    type->date_min = DBL_MAX;
    type->date_max = -DBL_MAX;
}

/**
 */
static void obstype_check(obstype* type)
{
    assert(type->name != NULL);
    if (type->issurface < 0)
        enkf_quit("\"%s\": ISSURFACE not specified\n");
    if (type->varname == NULL)
        enkf_quit("\"%s\": VAR not specified\n");
    if (type->hfunction == NULL)
        enkf_quit("\"%s\": HFUNCTION not specified\n");
}

/**
 */
static void obstype_print(obstype* type)
{
    enkf_printf("    NAME = %s\n", type->name);
    enkf_printf("      VAR = %s\n", type->varname);
    if (type->varname2 != NULL)
        enkf_printf("      SECONDARY VAR = %s\n", type->varname2);
    enkf_printf("      ID = %d\n", type->id);
    enkf_printf("      ISSURFACE = %s\n", (type->issurface) ? "yes" : "no");
    if (type->offset_fname != NULL)
        enkf_printf("      OFFSET = %s %s\n", type->offset_fname, type->offset_varname);
    enkf_printf("      HFUNCTION = %s\n", type->hfunction);
    enkf_printf("      ALLOWED MIN = %.3g\n", type->allowed_min);
    enkf_printf("      ALLOWED MAX = %.3g\n", type->allowed_max);
    enkf_printf("      ASYNCHRONOUS = %s", (type->isasync) ? "yes" : "no");
    if (type->isasync)
        enkf_printf(", DT = %.3f\n", type->async_tstep);
    else
        enkf_printf("\n");
    enkf_printf("      RFACTOR = %.3g\n", type->rfactor);
    if (type->xmin > -DBL_MAX || type->xmax < DBL_MAX || type->ymin > -DBL_MAX || type->ymax < DBL_MAX || type->zmin > -DBL_MAX || type->zmax < DBL_MAX)
        enkf_printf("      DOMAIN = %.3g %.3g %.3g %.3g %.3g %.3g\n", type->xmin, type->xmax, type->ymin, type->ymax, type->zmin, type->zmax);
}

/**
 */
void obstypes_read(char fname[], int* n, obstype** types, double rfactor_base)
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
            if (obstype_getid(*n, *types, token) >= 0)
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
            now->issurface = read_bool(token);
            if (now->issurface < 0)
                enkf_quit("%s, l.%d: could not convert \"%s\" to boolean", fname, line, token);
        } else if (strcasecmp(token, "VAR") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: VAR not specified", fname, line);
            if (now->varname != NULL)
                enkf_quit("%s, l.%d: VAR already specified", fname, line);
            now->varname = strdup(token);
        } else if (strcasecmp(token, "VAR2") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: VAR2 not specified", fname, line);
            if (now->varname2 != NULL)
                enkf_quit("%s, l.%d: VAR2 already specified", fname, line);
            now->varname2 = strdup(token);
        } else if (strcasecmp(token, "OFFSET") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: OFFSET file name not specified", fname, line);
            if (now->offset_fname != NULL)
                enkf_quit("%s, l.%d: OFFSET entry already specified", fname, line);
            now->offset_fname = strdup(token);
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: OFFSET variable name not specified", fname, line);
            now->offset_varname = strdup(token);
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
        } else if (strcasecmp(token, "RFACTOR") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: RFACTOR not specified", fname, line);
            if (!str2double(token, &now->rfactor))
                enkf_quit("%s, l.%d: could not convert \"%s\" to double", fname, line, token);
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
        } else
            enkf_quit("%s, l.%d: unknown token \"%s\"", fname, line, token);
    }
    fclose(f);

    for (i = 0; i < *n; ++i)
        (*types)[i].rfactor *= rfactor_base;

    for (i = 0; i < *n; ++i) {
        obstype_check(&(*types)[i]);
        obstype_print(&(*types)[i]);
    }
}

/**
 */
void obstypes_destroy(int n, obstype* types)
{
    int i;

    for (i = 0; i < n; ++i) {
        obstype* type = &types[i];

        free(type->name);
        free(type->varname);
        if (type->varname2 != NULL)
            free(type->varname2);
        free(type->hfunction);
        if (type->offset_fname != NULL) {
            free(type->offset_fname);
            free(type->offset_varname);
        }
    }

    free(types);
}

/**
 */
void obstypes_describeprm(void)
{
    enkf_printf("\n");
    enkf_printf("  Observation types parameter file format:\n");
    enkf_printf("\n");
    enkf_printf("    NAME      = <name>\n");
    enkf_printf("    VAR       = <model variable name>\n");
    enkf_printf("  [ VAR2       = <model variable name> ]\n");
    enkf_printf("    ISSURFACE = { yes | no }\n");
    enkf_printf("  [ OFFSET    = <file name> <variable name> ]    (none*)\n");
    enkf_printf("    HFUNCTION = <H function name>\n");
    enkf_printf("  [ ASYNC     = <time interval> ]                (synchronous*)\n");
    enkf_printf("  [ RFACTOR   = <rfactor> ]                      (1*)\n");
    enkf_printf("  [ MINVALUE  = <minimal allowed value> ]        (-inf*)\n");
    enkf_printf("  [ MAXVALUE  = <maximal allowed value> ]        (+inf*)\n");
    enkf_printf("  [ XMIN      = <minimal allowed X coordinate> ] (-inf*)\n");
    enkf_printf("  [ XMAX      = <maximal allowed X coordinate> ] (+inf*)\n");
    enkf_printf("  [ YMIN      = <minimal allowed Y coordinate> ] (-inf*)\n");
    enkf_printf("  [ YMAX      = <maximal allowed Y coordinate> ] (+inf*)\n");
    enkf_printf("  [ ZMIN      = <minimal allowed Z coordinate> ] (-inf*)\n");
    enkf_printf("  [ ZMAX      = <maximal allowed Z coordinate> ] (+inf*)\n");
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
