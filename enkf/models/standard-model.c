/******************************************************************************
 *
 * File:        standard-model.c        
 *
 * Created:     09/12/2012
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description: Contains standard model procedures
 *
 * Revisions:
 *
 *****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#if !defined(NO_GRIDUTILS)
#include <gridnodes.h>
#endif
#include "ncw.h"
#include "definitions.h"
#include "utils.h"
#include "grid.h"
#include "enkfprm.h"
#include "model.h"
#include "standard-model.h"

/**
 */
void standardmodel_getmemberfname(model* m, char ensdir[], char varname[], int mem, char fname[])
{
    snprintf(fname, MAXSTRLEN, "%s/mem%03d_%s.nc", ensdir, mem, varname);
}

/**
 */
int standardmodel_getmemberfname_async(model* m, char ensdir[], char varname[], char otname[], int mem, int t, char fname[])
{
    snprintf(fname, MAXSTRLEN, "%s/mem%03d_%s_%d.nc", ensdir, mem, varname, t);
    if (!file_exists(fname)) {
        snprintf(fname, MAXSTRLEN, "%s/mem%03d_%s.nc", ensdir, mem, varname);
        return 0;
    }
    return 1;
}

/**
 */
void standardmodel_getbgfname(model* m, char ensdir[], char varname[], char fname[])
{
    snprintf(fname, MAXSTRLEN, "%s/bg_%s.nc", ensdir, varname);
}

/**
 */
int standardmodel_getbgfname_async(model* m, char bgdir[], char varname[], char otname[], int t, char fname[])
{
    snprintf(fname, MAXSTRLEN, "%s/bg_%s_%d.nc", bgdir, varname, t);
    if (!file_exists(fname)) {
        snprintf(fname, MAXSTRLEN, "%s/bg_%s.nc", bgdir, varname);
        return 0;
    }
    return 1;
}

/**
 */
void standardmodel_adddata_2D(model* m, char* token, char* fname, int line)
{
    char seps[] = " =\t\n";
    char tag[MAXSTRLEN];
    char vfname[MAXSTRLEN];
    float** v = NULL;
    int mvid;
    int nx, ny;

    strncpy(tag, token, MAXSTRLEN);

    enkf_printf("    adding custom data \"%s\":\n", token);
    if ((token = strtok(NULL, seps)) == NULL)
        enkf_quit("%s, l.%d: file for \"%s\" not specified", fname, line, tag);
    enkf_printf("      file = %s\n", token);
    strcpy(vfname, token);
    if ((token = strtok(NULL, seps)) == NULL)
        enkf_quit("%s, l.%d: variable for \"%s\" not specified", fname, line, tag);
    enkf_printf("      variable = %s\n", token);

    mvid = model_getvarid(m, token);
    assert(mvid >= 0);
    model_getvardims(m, mvid, &nx, &ny, NULL);
    v = alloc2d(ny, nx, sizeof(float));
    readfield(vfname, 0, token, v[0]);
    model_addmodeldata(m, tag, ALLOCTYPE_2D, v);
}

/**
 */
void standardmodel_readfield(model* m, char fname[], int mem, int time, char varname[], int k, float* v)
{
    readfield(fname, k, varname, v);
}

/**
 */
void standardmodel_read3dfield(model* m, char fname[], int mem, int time, char varname[], float* v)
{
    read3dfield(fname, varname, v);
}

/**
 */
void standardmodel_writefield(model* m, char fname[], int time, char varname[], int k, float* v)
{
    writefield(fname, k, varname, v);
}
