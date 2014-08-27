/******************************************************************************
 *
 * File:        model.c        
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "nan.h"
#include "definitions.h"
#include "utils.h"
#include "grid.h"
#include "enkfprm.h"
#include "model.h"
#include "allmodels.h"

#define NMODELDATA_INC 10
#define NVAR_INC 10

typedef struct {
    char* tag;
    void* data;
} modeldata;

typedef struct {
    int id;
    char* name;
    double inflation;
} variable;

struct model {
    char* name;
    char* type;

    int nvar;
    variable* vars;

    void* grid;

    int ndata;
    modeldata* data;

    model_getmemberfname_fn getmemberfname;
    model_getmemberfnameasync_fn getmemberfname_async;
    model_getbgfname_fn getbgfname;
    model_getbgfnameasync_fn getbgfname_async;
    model_readfield_fn readfield;
    model_read3dfield_fn read3dfield;
    model_writefield_fn writefield;
};

/**
 */
static void variable_new(variable * v, int id, char* name)
{
    v->id = id;
    v->name = strdup(name);
    v->inflation = 1.0;
}

/**
 */
static void variables_destroy(int n, variable * vars)
{
    int i;

    for (i = 0; i < n; ++i) {
        variable* v = &vars[i];

        free(v->name);
    }

    free(vars);
}

/**
 */
model* model_create(enkfprm* prm)
{
    model* m = calloc(1, sizeof(model));
    char* modelprm = prm->modelprm;
    char* gridprm = prm->gridprm;

    /*
     * read model parameter file
     */
    {
        FILE* f = NULL;
        char buf[MAXSTRLEN];
        int line;
        variable* now = NULL;

        /*
         * get model tag and type
         */
        f = enkf_fopen(modelprm, "r");
        line = 0;
        while (fgets(buf, MAXSTRLEN, f) != NULL) {
            char seps[] = " =\t\n";
            char* token = NULL;

            line++;
            if (buf[0] == '#')
                continue;
            if ((token = strtok(buf, seps)) == NULL)
                continue;
            if (strcasecmp(token, "NAME") == 0) {
                if ((token = strtok(NULL, seps)) == NULL)
                    enkf_quit("%s, l.%d: NAME not specified", modelprm, line);
                else if (m->name != NULL)
                    enkf_quit("%s, l.%d: NAME specified twice", modelprm, line);
                else
                    m->name = strdup(token);
            } else if (strcasecmp(token, "TYPE") == 0) {
                if ((token = strtok(NULL, seps)) == NULL)
                    enkf_quit("%s, l.%d: TYPE not specified", modelprm, line);
                else if (m->type != NULL)
                    enkf_quit("%s, l.%d: TYPE specified twice", modelprm, line);
                else
                    m->type = strdup(token);
            } else if (strncasecmp(token, "VAR", 3) == 0) {
                if ((token = strtok(NULL, seps)) == NULL)
                    enkf_quit("%s, l.%d: VAR not specified", modelprm, line);
                if (m->nvar % NVAR_INC == 0)
                    m->vars = realloc(m->vars, (m->nvar + NVAR_INC) * sizeof(variable));
                now = &m->vars[m->nvar];
                variable_new(now, m->nvar, token);
                m->nvar++;
            } else if (strcasecmp(token, "INFLATION") == 0) {
                if (now == NULL)
                    enkf_quit("%s, l.%d: VAR not specified", modelprm, line);
                if (now->inflation != 1.0)
                    enkf_quit("%s, l.%d: INFLATION already specified for \"%s\"", modelprm, line, now->name);
                if ((token = strtok(NULL, seps)) == NULL)
                    enkf_quit("%s, l.%d: INFLATION not specified", modelprm, line);
                if (!str2double(token, &now->inflation))
                    enkf_quit("%s, l.%d: could not convert \"%s\" to double", modelprm, line, token);
            } else if (strcasecmp(token, "DATA") == 0)
                continue;
            else
                enkf_quit("%s, l.%d: unknown token \"%s\"", modelprm, line, token);
        }                       /* while reading modelprm */
        fclose(f);
        assert(m->name != NULL);
        assert(m->type != NULL);
        assert(m->nvar > 0);
    }

    /*
     * set inflations
     */
    {
        int i;

        for (i = 0; i < m->nvar; ++i)
            m->vars[i].inflation *= prm->inflation_base;
        prm->inflation_base = NaN;
    }

    /*
     * set the grid
     */
    get_modelsetgridfn(m->type) (m, gridprm);
    assert(m->grid !=NULL);

    /*
     * finish the model setup
     */
    get_modelsetupfn(m->type) (m, modelprm);
    assert(m->getmemberfname != NULL);
    assert(m->getmemberfname_async != NULL);
    assert(m->getbgfname != NULL);
    assert(m->getbgfname_async != NULL);
    assert(m->readfield != NULL);
    assert(m->read3dfield != NULL);
    assert(m->writefield != NULL);

    return m;
}

/**
 */
static void model_freemodeldata(model* m)
{
    int i;

    for (i = 0; i < m->ndata; ++i) {
        modeldata* data = &m->data[i];

        if (strcasecmp(data->tag, "MSL") == 0)
            free2d(data->data);
        else
            free(data->data);
        free(data->tag);
    }
    free(m->data);
    m->ndata = 0;
}

/**
 */
void model_destroy(model* m)
{
    free(m->name);
    free(m->type);
    grid_destroy(m->grid);

    variables_destroy(m->nvar, m->vars);
    model_freemodeldata(m);
    free(m);
}

/**
 */
void model_print(model* m, char offset[])
{
    int i;

    enkf_printf("%smodel info:\n", offset);
    enkf_printf("%s  name = %s\n", offset, m->name);
    enkf_printf("%s  type = %s\n", offset, m->type);
    enkf_printf("%s  %d variables:\n", offset, m->nvar);
    for (i = 0; i < m->nvar; ++i) {
        variable* v = &m->vars[i];

        enkf_printf("%s    %s:\n", offset, v->name);
        enkf_printf("%s      inflation = %.3f\n", offset, v->inflation);
    }
    enkf_printf("%s  %d modeldata:\n", offset, m->ndata);
    for (i = 0; i < m->ndata; ++i)
        enkf_printf("%s    %s\n", offset, m->data[i].tag);
}

/**
 */
void model_describeprm(void)
{
    enkf_printf("\n");
    enkf_printf("  Model parameter file format:\n");
    enkf_printf("\n");
    enkf_printf("    NAME        = <name>\n");
    enkf_printf("    TYPE        = <type>\n");
    enkf_printf("  [ DATA <tag>  = <data> ]\n");
    enkf_printf("    ...\n");
    enkf_printf("\n");
    enkf_printf("    VAR         = <name>\n");
    enkf_printf("  [ INFLATION   = <value> ]\n");
    enkf_printf("\n");
    enkf_printf("  [ <more of the above blocks> ]\n");
    enkf_printf("\n");
    enkf_printf("  Notes:\n");
    enkf_printf("    1. [ ... ] denotes an optional input\n");
    enkf_printf("    2. < ... > denotes a description of an entry\n");
    enkf_printf("    3. ... denotes repeating the previous item an arbitrary number of times\n");
    enkf_printf("\n");
}

/**
 */
void model_setgetmemberfname_fn(model* m, model_getmemberfname_fn fn)
{
    m->getmemberfname = fn;
}

/**
 */
void model_setgetmemberfnameasync_fn(model* m, model_getmemberfnameasync_fn fn)
{
    m->getmemberfname_async = fn;
}

/**
 */
void model_setbgfname_fn(model* m, model_getbgfname_fn fn)
{
    m->getbgfname = fn;
}

/**
 */
void model_setbgfnameasync_fn(model* m, model_getbgfnameasync_fn fn)
{
    m->getbgfname_async = fn;
}

/**
 */
void model_setreadfield_fn(model* m, model_readfield_fn fn)
{
    m->readfield = fn;
}

/**
 */
void model_setread3dfield_fn(model* m, model_read3dfield_fn fn)
{
    m->read3dfield = fn;
}

/**
 */
void model_setwritefield_fn(model* m, model_writefield_fn fn)
{
    m->writefield = fn;
}

void model_setgrid(model* m, void* g)
{
    m->grid = g;
}

/**
 */
void model_addmodeldata(model* m, char tag[], void* data)
{
    if (m->ndata % NMODELDATA_INC == 0)
        m->data = realloc(m->data, (m->ndata + NMODELDATA_INC) * sizeof(modeldata));
    m->data[m->ndata].tag = strdup(tag);
    m->data[m->ndata].data = data;
    m->ndata++;
}

/**
 */
void* model_getmodeldata(model* m, char tag[])
{
    int i;

    for (i = 0; i < m->ndata; ++i) {
        modeldata* data = &m->data[i];

        if (strcasecmp(data->tag, tag) == 0)
            return data->data;
    }

    enkf_quit("  getmodeldata(): could not find data \"%s\" for model \"%s\"", tag, m->name);
    return NULL;
}

/**
 */
int model_getnvar(model* m)
{
    return m->nvar;
}

/**
 */
char* model_getvarname(model* m, int varid)
{
    return m->vars[varid].name;
}

/**
 */
float model_getvarinflation(model* m, int varid)
{
    return m->vars[varid].inflation;
}

/**
 */
void model_getdims(model* m, int* ni, int* nj, int* nk)
{
    grid_getdims(m->grid, ni, nj, nk);
}

/**
 */
void* model_getgrid(model* m)
{
    return m->grid;
}

/**
 */
int model_getlontype(model* m)
{
    return grid_getlontype(m->grid);
}

/**
 */
float** model_getdepth(model* m)
{
    return grid_getdepth(m->grid);
}

/**
 */
int** model_getnumlevels(model* m)
{
    return grid_getnumlevels(m->grid);
}

/**
 */
void model_getmemberfname(model* m, char ensdir[], char varname[], int mem, char fname[])
{
    m->getmemberfname(m, ensdir, varname, mem, fname);
}

/**
 */
int model_getmemberfname_async(model* m, char ensdir[], char varname[], char otname[], int mem, int time, char fname[])
{
    return m->getmemberfname_async(m, ensdir, varname, otname, mem, time, fname);
}

/**
 */
void model_getbgfname(model* m, char ensdir[], char varname[], char fname[])
{
    m->getbgfname(m, ensdir, varname, fname);
}

/**
 */
int model_getbgfname_async(model* m, char ensdir[], char varname[], char otname[], int time, char fname[])
{
    return m->getbgfname_async(m, ensdir, varname, otname, time, fname);
}

/**
 */
int model_xy2fij(model* m, double x, double y, double* fi, double* fj)
{
    int** numlevels = grid_getnumlevels(m->grid);
    int lontype = grid_getlontype(m->grid);
    int i1, i2, j1, j2;

    if (lontype == LONTYPE_180) {
        if (x > 180.0)
            x -= 360.0;
    } else if (lontype == LONTYPE_360) {
        if (x < 0.0)
            x += 360.0;
    }

    grid_getxy2fijfn(m->grid) (m->grid, x, y, fi, fj);

    if (isnan(*fi + *fj))
        return STATUS_OUTSIDEGRID;

    i1 = floor(*fi);
    i2 = ceil(*fi);
    j1 = floor(*fj);
    j2 = ceil(*fj);
    if (numlevels[j1][i1] == 0 && numlevels[j1][i2] == 0 && numlevels[j2][i1] == 0 && numlevels[j2][i2] == 0) {
        *fi = NaN;
        *fj = NaN;
        return STATUS_LAND;
    }
    return STATUS_OK;
}

/**
 */
int model_fij2xy(model* m, double fi, double fj, double* x, double* y)
{
    grid_getfij2xyfn(m->grid) (m->grid, fi, fj, x, y);

    if (isnan(*x + *y))
        return STATUS_OUTSIDEGRID;
    return STATUS_OK;
}

/**
 */
int model_z2fk(model* m, double fi, double fj, double z, double* fk)
{
    int** numlevels = grid_getnumlevels(m->grid);
    int i1, i2, j1, j2, k2;

    if (isnan(fi + fj)) {
        *fk = NaN;
        return STATUS_OUTSIDEGRID;
    }

    grid_getz2fkfn(m->grid) (m->grid, fi, fj, z, fk);

    if (isnan(*fk))
        return STATUS_OUTSIDEGRID;

    i1 = floor(fi);
    i2 = ceil(fi);
    j1 = floor(fj);
    j2 = ceil(fj);
    k2 = ceil(*fk);
    if (numlevels[j1][i1] <= k2 && numlevels[j1][i2] <= k2 && numlevels[j2][i1] <= k2 && numlevels[j2][i2] <= k2) {
        *fk = NaN;
        return STATUS_LAND;
    } else if (numlevels[j1][i1] <= k2 || numlevels[j1][i2] <= k2 || numlevels[j2][i1] <= k2 || numlevels[j2][i2] <= k2) {
        float** depth = model_getdepth(m);
        int ni, nj, nk;
        double v;

        model_getdims(m, &ni, &nj, &nk);
        v = interpolate2d(fi, fj, ni, nj, depth, numlevels);
        if (z > v)
            return STATUS_LAND;
    }

    return STATUS_OK;
}

/**
 */
void model_readfield(model* m, char fname[], int mem, int time, char varname[], int k, float* v)
{
    m->readfield(m, fname, mem, time, varname, k, v);
}

/**
 */
void model_read3dfield(model* m, char fname[], int mem, int time, char varname[], float* v)
{
    m->read3dfield(m, fname, mem, time, varname, v);
}

/**
 */
void model_writefield(model* m, char fname[], int time, char varname[], int k, float* v)
{
    m->writefield(m, fname, time, varname, k, v);
}
