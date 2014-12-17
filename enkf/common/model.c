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
#define GRID_INC 10

typedef struct {
    char* tag;
    int alloctype;
    void* data;
} modeldata;

struct variable {
    int id;
    char* name;
    int gridid;
    double inflation;
};

struct model {
    char* name;
    char* type;

    int nvar;
    variable* vars;

    int ngrid;
    void** grids;

    int ndata;
    modeldata* data;

    model_getmemberfname_fn getmemberfname;
    model_getmemberfnameasync_fn getmemberfname_async;
    model_getbgfname_fn getbgfname;
    model_getbgfnameasync_fn getbgfname_async;
    model_readfield_fn readfield;
    model_read3dfield_fn read3dfield;
    model_writefield_fn writefield;
    model_adddata_fn adddatafn;
};

/**
 */
static void variable_new(variable * v, int id, char* name)
{
    v->id = id;
    v->name = strdup(name);
    v->gridid = -1;
    v->inflation = NaN;
}

/**
 */
static void model_destroyvars(model* m)
{
    int i;

    for (i = 0; i < m->nvar; ++i)
        free(m->vars[i].name);
    free(m->vars);
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
         * get model tag, type and variables
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
                    enkf_quit("%s, l.%d: TYPE already specified", modelprm, line);
                else
                    m->type = strdup(token);
                /*
                 * set the grid
                 */
                get_modelsetgridfn(m->type) (m, gridprm);
                /*
                 * do custom model setup
                 */
                get_modelsetupfn(m->type) (m, modelprm);
            } else if (strncasecmp(token, "VAR", 3) == 0) {
                int i;

                if ((token = strtok(NULL, seps)) == NULL)
                    enkf_quit("%s, l.%d: VAR not specified", modelprm, line);
                for (i = 0; i < m->nvar; ++i)
                    if (strcasecmp(m->vars[i].name, token) == 0)
                        enkf_quit("%s, l.%d: VAR \"%s\" already specified", modelprm, line, token);
                if (m->nvar % NVAR_INC == 0)
                    m->vars = realloc(m->vars, (m->nvar + NVAR_INC) * sizeof(variable));
                now = &m->vars[m->nvar];
                variable_new(now, m->nvar, token);
                m->nvar++;
            } else if (strcasecmp(token, "GRID") == 0) {
                int i;

                if (now == NULL)
                    enkf_quit("%s, l.%d: VAR not specified", modelprm, line);
                if (now->gridid >= 0)
                    enkf_quit("%s, l.%d: GRID already specified for \"%s\"", modelprm, line, now->name);
                if ((token = strtok(NULL, seps)) == NULL)
                    enkf_quit("%s, l.%d: GRID not specified", modelprm, line);
                for (i = 0; i < m->ngrid; ++i)
                    if (strcasecmp(token, grid_getname(m->grids[i])) == 0) {
                        now->gridid = i;
                        break;
                    }
                if (i == m->ngrid)
                    enkf_quit("%s, l.%d: grid \"%s\" not specified", modelprm, line, token);
            } else if (strcasecmp(token, "INFLATION") == 0) {
                if (now == NULL)
                    enkf_quit("%s, l.%d: VAR not specified", modelprm, line);
                if (!isnan(now->inflation))
                    enkf_quit("%s, l.%d: INFLATION already specified for \"%s\"", modelprm, line, now->name);
                if ((token = strtok(NULL, seps)) == NULL)
                    enkf_quit("%s, l.%d: INFLATION not specified", modelprm, line);
                if (!str2double(token, &now->inflation))
                    enkf_quit("%s, l.%d: could not convert \"%s\" to double", modelprm, line, token);
            } else if (strcasecmp(token, "DATA") == 0) {
                ;
            } else
                enkf_quit("%s, l.%d: unknown token \"%s\"", modelprm, line, token);
        }                       /* while reading modelprm */
        fclose(f);
        assert(m->name != NULL);
        assert(m->type != NULL);
        assert(m->nvar > 0);
        {
            int i;

            for (i = 0; i < m->nvar; ++i)
                if (m->vars[i].gridid == -1) {
                    if (m->ngrid == 1)
                        m->vars[i].gridid = 0;
                    else
                        enkf_quit("%s: grid not specified for variable \"%s\"\n", modelprm, m->vars[i].name);
                }
        }
        /*
         * get model data (this needs variables to be already set)
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
            else if (strcasecmp(token, "DATA") == 0) {
                if ((token = strtok(NULL, seps)) == NULL)
                    enkf_quit("%s, l.%d: DATA tag not specified", modelprm, line);
                model_addcustomdata(m, token, modelprm, line);
            }
        }
    }

    /*
     * set inflations
     */
    {
        int i;

        for (i = 0; i < m->nvar; ++i)
            if (isnan(m->vars[i].inflation))
                m->vars[i].inflation = prm->inflation_base;
        prm->inflation_base = NaN;
    }

    model_print(m, "    ");

    assert(m->ngrid > 0);

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

        if (data->alloctype == ALLOCTYPE_1D)
            free(data->data);
        else if (data->alloctype == ALLOCTYPE_2D)
            free2d(data->data);
        else if (data->alloctype == ALLOCTYPE_3D)
            free3d(data->data);
        else
            enkf_quit("programming error");
        free(data->tag);
    }
    free(m->data);
    m->ndata = 0;
}

/**
 */
static void model_destroygrids(model* m)
{
    int i;

    for (i = 0; i < m->ngrid; ++i)
        grid_destroy(m->grids[i]);

    free(m->grids);
}

/**
 */
void model_destroy(model* m)
{
    free(m->name);
    free(m->type);
    model_destroygrids(m);
    model_destroyvars(m);
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
    for (i = 0; i < m->ndata; ++i) {
        enkf_printf("%s    %s:\n", offset, m->data[i].tag);
        if (m->data[i].alloctype == ALLOCTYPE_1D)
            enkf_printf("%s      type = 1D\n", offset);
        else if (m->data[i].alloctype == ALLOCTYPE_2D)
            enkf_printf("%s      type = 2D\n", offset);
        else if (m->data[i].alloctype == ALLOCTYPE_3D)
            enkf_printf("%s      type = 3D\n", offset);
    }
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
    enkf_printf("    GRID        = <name>\n");
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

/**
 */
void model_setgrid(model* m, void* g)
{
    if (m->ngrid % GRID_INC == 0)
        m->grids = realloc(m->grids, (m->ngrid + GRID_INC) * sizeof(void*));
    m->grids[m->ngrid] = g;
    m->ngrid++;
}

/**
 */
void model_setadddata_fn(model* m, model_adddata_fn fn)
{
    m->adddatafn = fn;
}

/**
 */
void model_addmodeldata(model* m, char tag[], int alloctype, void* data)
{
    int i;

    for (i = 0; i < m->ndata; ++i)
        if (strcmp(tag, m->data[i].tag) == 0)
            enkf_quit("model data tag \"%s\" already in use", tag);

    if (m->ndata % NMODELDATA_INC == 0)
        m->data = realloc(m->data, (m->ndata + NMODELDATA_INC) * sizeof(modeldata));
    m->data[m->ndata].tag = strdup(tag);
    m->data[m->ndata].alloctype = alloctype;
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
char* model_gettype(model* m)
{
    return m->type;
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
int model_getvarid(model* m, char* varname)
{
    int i;

    for (i = 0; i < m->nvar; ++i)
        if (strcmp(m->vars[i].name, varname) == 0)
            return i;

    return -1;
}

/**
 */
float model_getvarinflation(model* m, int varid)
{
    return m->vars[varid].inflation;
}

/**
 */
void model_getvardims(model* m, int vid, int* ni, int* nj, int* nk)
{
    grid_getdims(m->grids[m->vars[vid].gridid], ni, nj, nk);
}

/**
 */
void* model_getvargrid(model* m, int vid)
{
    return m->grids[m->vars[vid].gridid];
}

/**
 */
int model_getvargridid(model* m, int vid)
{
    return m->vars[vid].gridid;
}

/**
 */
int model_getngrid(model* m)
{
    return m->ngrid;
}

/**
 */
void* model_getgridbyid(model* m, int gridid)
{
    return m->grids[gridid];
}

/**
 */
void* model_getgridbyname(model* m, char name[])
{
    int i;

    for (i = 0; i < m->ngrid; ++i)
        if (strcmp(grid_getname(m->grids[i]), name) == 0)
            return m->grids[i];
    return NULL;
}

/**
 */
int model_getlontype(model* m, int vid)
{
    return grid_getlontype(m->grids[m->vars[vid].gridid]);
}

/**
 */
float** model_getdepth(model* m, int vid)
{
    return grid_getdepth(m->grids[m->vars[vid].gridid]);
}

/**
 */
int** model_getnumlevels(model* m, int vid)
{
    return grid_getnumlevels(m->grids[m->vars[vid].gridid]);
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
int model_xy2fij(model* m, int vid, double x, double y, double* fi, double* fj)
{
    void* grid = m->grids[m->vars[vid].gridid];
    int** numlevels = grid_getnumlevels(grid);
    int lontype = grid_getlontype(grid);
    int i1, i2, j1, j2;

    if (lontype == LONTYPE_180) {
        if (x > 180.0)
            x -= 360.0;
    } else if (lontype == LONTYPE_360) {
        if (x < 0.0)
            x += 360.0;
    }

    grid_xy2fij(grid, x, y, fi, fj);

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
int model_fij2xy(model* m, int vid, double fi, double fj, double* x, double* y)
{
    void* grid = m->grids[m->vars[vid].gridid];

    grid_fij2xy(grid, fi, fj, x, y);

    if (isnan(*x + *y))
        return STATUS_OUTSIDEGRID;
    return STATUS_OK;
}

/**
 */
int model_z2fk(model* m, int vid, double fi, double fj, double z, double* fk)
{
    void* grid = m->grids[m->vars[vid].gridid];
    int** numlevels = grid_getnumlevels(grid);
    int i1, i2, j1, j2, k2;

    if (isnan(fi + fj)) {
        *fk = NaN;
        return STATUS_OUTSIDEGRID;
    }

    grid_z2fk(grid, fi, fj, z, fk);

    if (isnan(*fk))
        return STATUS_OUTSIDEGRID;

    if (grid_getvtype(grid) == GRIDVTYPE_SIGMA)
         return STATUS_OK;

    /*
     * a depth check for z-grid:
     */
    i1 = floor(fi);
    i2 = ceil(fi);
    j1 = floor(fj);
    j2 = ceil(fj);
    k2 = floor(*fk);
    if (numlevels[j1][i1] <= k2 && numlevels[j1][i2] <= k2 && numlevels[j2][i1] <= k2 && numlevels[j2][i2] <= k2) {
        *fk = NaN;
        return STATUS_LAND;
    } else if (numlevels[j1][i1] <= k2 || numlevels[j1][i2] <= k2 || numlevels[j2][i1] <= k2 || numlevels[j2][i2] <= k2) {
        float** depth = grid_getdepth(grid);
        int ni, nj;
        double v;

        grid_getdims(grid, &ni, &nj, NULL);

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

/**
 */
void model_addcustomdata(model* m, char* token, char* fname, int line)
{
    m->adddatafn(m, token, fname, line);
}
