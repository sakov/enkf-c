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
#include "definitions.h"
#include "utils.h"
#include "grid.h"
#include "gridprm.h"
#include "enkfprm.h"
#include "model.h"

#define NMODELDATA_INC 10
#define NVAR_INC 10
#define GRID_INC 10

typedef struct {
    char* tag;
    int vid;                    /* (to access grid dimensions) */
    int alloctype;
    void* data;
} modeldata;

struct variable {
    int id;
    char* name;
    int gridid;
    double inflation;
    double inf_ratio;
    /*
     * if not NaN then the code will "propagate" the variable as follows:
     *   v <- deflation * v + (1 - deflation^2)^1/2 * s, 
     * where s ~ N(0, sigma) is a field-wide value (different across the
     * ensemble)
     */
    double deflation;
    double sigma;
};

struct model {
    char* name;

    int nvar;
    variable* vars;

    int ngrid;
    void** grids;

    int ndata;
    modeldata* data;
};

/**
 */
static void variable_new(variable * v, int id, char* name)
{
    v->id = id;
    v->name = strdup(name);
    v->gridid = -1;
    v->inflation = NAN;
    v->inf_ratio = NAN;
    v->deflation = NAN;
    v->sigma = NAN;
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
static void model_setgrids(model* m, char gfname[])
{
    int ngrid = 0;
    gridprm* prm = NULL;
    int i;

    gridprm_create(gfname, &ngrid, &prm);
    assert(ngrid > 0);

    for (i = 0; i < ngrid; ++i) {
        grid* g = NULL;

        g = grid_create(&prm[i], i);
        grid_settocartesian_fn(g, ll2xyz);
        model_setgrid(m, g);
    }

    gridprm_destroy(ngrid, prm);
}

static void model_checkvars(model* m, char* modelprm)
{
    int i;

    for (i = 0; i < m->nvar; ++i) {
        variable* v = &m->vars[i];

        if (!isfinite(v->inflation) || v->inflation <= 0)
            enkf_quit("\"%s\": \"%s\": inflation = %.3g\n", modelprm, v->name, v->inflation);
        if (v->deflation <= 0)
            enkf_quit("\"%s\": \"%s\": deflation = %.3g\n", modelprm, v->name, v->deflation);
    }
}

/**
 */
model* model_create(enkfprm* prm)
{
    model* m = calloc(1, sizeof(model));
    char* modelprm = prm->modelprm;
    char* gridprm = prm->gridprm;
    int i;

    model_setgrids(m, gridprm);

    for (i = 0; i < m->ngrid; ++i) {
        grid* g = m->grids[i];

        if (grid_getstride(g) == 0)
            grid_setstride(g, prm->stride);
    }

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
            } else if (strncasecmp(token, "VAR", 3) == 0) {
                if ((token = strtok(NULL, seps)) == NULL)
                    enkf_quit("%s, l.%d: VAR not specified", modelprm, line);
                for (i = 0; i < m->nvar; ++i)
                    if (strcmp(m->vars[i].name, token) == 0)
                        enkf_quit("%s, l.%d: VAR \"%s\" already specified", modelprm, line, token);
                if (m->nvar % NVAR_INC == 0)
                    m->vars = realloc(m->vars, (m->nvar + NVAR_INC) * sizeof(variable));
                now = &m->vars[m->nvar];
                variable_new(now, m->nvar, token);
                m->nvar++;
            } else if (strcasecmp(token, "GRID") == 0) {
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
                if ((token = strtok(NULL, seps)) != NULL) {
                    if (!str2double(token, &now->inf_ratio))
                        enkf_quit("%s, l.%d: could not convert \"%s\" to double", modelprm, line, token);
                }
            } else if (strcasecmp(token, "RANDOMISE") == 0 || strcasecmp(token, "RANDOMIZE") == 0) {
                if (!isnan(now->deflation))
                    enkf_quit("%s, l.%d: randomisation multiple already specified for \"%s\"", modelprm, line, now->name);
                if ((token = strtok(NULL, seps)) == NULL)
                    enkf_quit("%s, l.%d: RANDOMISE parameters not specified", modelprm, line);
                if (!str2double(token, &now->deflation))
                    enkf_quit("%s, l.%d: could not convert \"%s\" to double", modelprm, line, token);
                if ((token = strtok(NULL, seps)) == NULL)
                    enkf_quit("%s, l.%d: RANDOMISE STD not specified", modelprm, line);
                if (!str2double(token, &now->sigma))
                    enkf_quit("%s, l.%d: could not convert \"%s\" to double", modelprm, line, token);
            } else
                enkf_quit("%s, l.%d: unknown token \"%s\"", modelprm, line, token);
        }                       /* while reading modelprm */
        fclose(f);
        assert(m->name != NULL);
        assert(m->nvar > 0);
        {
            for (i = 0; i < m->nvar; ++i)
                if (m->vars[i].gridid == -1) {
                    if (m->ngrid == 1)
                        m->vars[i].gridid = 0;
                    else
                        enkf_quit("%s: grid not specified for variable \"%s\"\n", modelprm, m->vars[i].name);
                }
        }
    }

    /*
     * set inflations
     */
    {
        for (i = 0; i < m->nvar; ++i)
            if (isnan(m->vars[i].inflation)) {
                m->vars[i].inflation = prm->inflation_base;
                m->vars[i].inf_ratio = prm->inf_ratio;
            }
        prm->inflation_base = NAN;
        prm->inf_ratio = NAN;
    }
    /*
     * check randomisation parameters
     */
    {
        for (i = 0; i < m->nvar; ++i) {
            if (isnan(m->vars[i].deflation))
                continue;
            if (m->vars[i].deflation < 0.0 || m->vars[i].deflation > 1.0)
                enkf_quit("deflation = %.3g for variable \"%s\" (should be between 0 and 1)", m->vars[i].deflation, m->vars[i].name);
        }
    }
    model_print(m, "    ");

    model_checkvars(m, modelprm);

    assert(m->ngrid > 0);

    return m;
}

/**
 */
static void model_freemodeldata(model* m)
{
    int i;

    if (m->ndata == 0)
        return;

    for (i = 0; i < m->ndata; ++i) {
        modeldata* data = &m->data[i];

        if (data->alloctype == ALLOCTYPE_1D)
            free(data->data);
        else if (data->alloctype == ALLOCTYPE_2D)
            free(data->data);
        else if (data->alloctype == ALLOCTYPE_3D)
            free(data->data);
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
    enkf_printf("%s  %d variables:\n", offset, m->nvar);
    for (i = 0; i < m->nvar; ++i) {
        variable* v = &m->vars[i];

        enkf_printf("%s    %s:\n", offset, v->name);
        enkf_printf("%s      grid = \"%s\"\n", offset, grid_getname(model_getgridbyid(m, v->gridid)));
        if (isnan(v->inf_ratio))
            enkf_printf("%s      inflation = %.3f PLAIN\n", offset, v->inflation);
        else
            enkf_printf("%s      inflation = %.3f %.2f\n", offset, v->inflation, v->inf_ratio);
        if (!isnan(v->deflation))
            enkf_printf("%s      randomise: deflation = %.3f, sigma = %.3f\n", offset, v->deflation, v->sigma);
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
    enkf_printf("    NAME      = <name>\n");
    enkf_printf("    VAR       = <name>\n");
    enkf_printf("  [ GRID      = <name> ]                    (# grids > 1)\n");
    enkf_printf("  [ INFLATION = <value> [<value> | PLAIN] ]\n");
    enkf_printf("  [ RANDOMISE <deflation> <sigma> ]\n");
    enkf_printf("\n");
    enkf_printf("  [ <more of the above blocks> ]\n");
    enkf_printf("\n");
    enkf_printf("  Notes:\n");
    enkf_printf("    1. [ ... ] denotes an optional input\n");
    enkf_printf("    2. < ... > denotes a description of an entry\n");
    enkf_printf("    3. ( ... ) is a note\n");
    enkf_printf("    4. ... denotes repeating the previous item an arbitrary number of times\n");
    enkf_printf("\n");
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
void model_adddata(model* m, char tag[], int vid, int alloctype, void* data)
{
    modeldata* mdata;
    int i;

    for (i = 0; i < m->ndata; ++i)
        if (strcmp(tag, m->data[i].tag) == 0)
            enkf_quit("model data tag \"%s\" already in use", tag);

    if (m->ndata % NMODELDATA_INC == 0)
        m->data = realloc(m->data, (m->ndata + NMODELDATA_INC) * sizeof(modeldata));

    mdata = &m->data[m->ndata];
    mdata->tag = strdup(tag);
    mdata->vid = vid;
    mdata->alloctype = alloctype;
    mdata->data = data;
    m->ndata++;
}

/**
 */
void model_addorreplacedata(model* m, char tag[], int vid, int alloctype, void* data)
{
    modeldata* mdata;
    int i;
    int ni, nj, nk;

    for (i = 0; i < m->ndata; ++i)
        if (strcmp(tag, m->data[i].tag) == 0)
            break;

    mdata = &m->data[i];
    if (i == m->ndata) {
        if (m->ndata % NMODELDATA_INC == 0)
            m->data = realloc(m->data, (m->ndata + NMODELDATA_INC) * sizeof(modeldata));
        mdata->tag = strdup(tag);
        mdata->vid = vid;
        mdata->alloctype = alloctype;
        m->ndata++;
    } else
        assert(mdata->alloctype == alloctype);

    model_getvardims(m, mdata->vid, &ni, &nj, &nk);
    if (mdata->alloctype == ALLOCTYPE_1D)
        memcpy(mdata->data, data, nk * sizeof(float));
    else if (mdata->alloctype == ALLOCTYPE_2D)
        mdata->data = copy2d(data, nj, ni, sizeof(float));
    else if (mdata->alloctype == ALLOCTYPE_3D)
        mdata->data = copy3d(data, nk, nj, ni, sizeof(float));
    else
        enkf_quit("programming error");
}

/**
 */
void* model_getdata(model* m, char tag[])
{
    int i;

    for (i = 0; i < m->ndata; ++i) {
        modeldata* data = &m->data[i];

        if (strcasecmp(data->tag, tag) == 0)
            return data->data;
    }

    return NULL;
}

/**
 */
int model_getdataalloctype(model* m, char tag[])
{
    int i;

    for (i = 0; i < m->ndata; ++i) {
        modeldata* data = &m->data[i];

        if (strcasecmp(data->tag, tag) == 0)
            return data->alloctype;
    }

    return ALLOCTYPE_NONE;
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
int model_getvarid(model* m, char* varname, int hastosucceed)
{
    int i;

    if (varname == NULL)
        return -1;

    for (i = 0; i < m->nvar; ++i)
        if (strcmp(m->vars[i].name, varname) == 0)
            return i;

    if (hastosucceed)
        enkf_quit("model_getvarid(): can not find variable \"%s\" in the model", varname);

    return -1;
}

/**
 */
void model_getvarinflation(model* m, int varid, float* inflation, double* inf_ratio)
{
    *inflation = (float) m->vars[varid].inflation;
    *inf_ratio = m->vars[varid].inf_ratio;
}

/**
 */
double model_getvardeflation(model* m, int varid)
{
    return m->vars[varid].deflation;
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

    enkf_quit("model_getgridbyname(): found no grid named \"%s\"", name);
    return NULL;
}

/**
 */
double model_getlonbase(model* m, int vid)
{
    return grid_getlonbase(m->grids[m->vars[vid].gridid]);
}

/**
 */
float** model_getdepth(model* m, int vid, int musthave)
{
    void* grid = m->grids[m->vars[vid].gridid];
    float** depth = grid_getdepth(grid);

    if (musthave && depth == NULL)
        enkf_quit("DEPTHVARNAME not specified for grid \"%s\"", grid_getname(grid));

    return depth;
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
    snprintf(fname, MAXSTRLEN, "%s/mem%03d_%s.nc", ensdir, mem, varname);
}

/**
 */
int model_getmemberfname_async(model* m, char ensdir[], char varname[], char otname[], int mem, int t, char fname[])
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
void model_getbgfname(model* m, char ensdir[], char varname[], char fname[])
{
    snprintf(fname, MAXSTRLEN, "%s/bg_%s.nc", ensdir, varname);
}

/**
 */
int model_getbgfname_async(model* m, char bgdir[], char varname[], char otname[], int t, char fname[])
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
int model_xy2fij(model* m, int vid, double x, double y, double* fi, double* fj)
{
    void* grid = m->grids[m->vars[vid].gridid];
    int isperiodic_x = grid_isperiodic_x(grid);
    int** numlevels = grid_getnumlevels(grid);
    double lonbase = grid_getlonbase(grid);
    int ni, nj;
    int i1, i2, j1, j2;

    if (!isnan(lonbase)) {
        if (x < lonbase)
            x += 360.0;
        else if (x >= lonbase + 360.0)
            x -= 360.0;
    }

    grid_xy2fij(grid, x, y, fi, fj);

    if (isnan(*fi + *fj))
        return STATUS_OUTSIDEGRID;

    /*
     * Note that this section should be consistent with similar sections in
     * interpolate2d() and interpolate3d().
     */
    i1 = floor(*fi);
    i2 = ceil(*fi);
    j1 = floor(*fj);
    j2 = ceil(*fj);

    model_getvardims(m, vid, &ni, &nj, NULL);
    if (i1 == -1)
        i1 = (isperiodic_x) ? ni - 1 : i2;
    if (i2 == ni)
        i2 = (isperiodic_x) ? 0 : i1;
    if (j1 == -1)
        j1 = j2;
    if (j2 == nj)
        j2 = j1;

    if (numlevels[j1][i1] == 0 && numlevels[j1][i2] == 0 && numlevels[j2][i1] == 0 && numlevels[j2][i2] == 0) {
        *fi = NAN;
        *fj = NAN;
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
int model_ij2xy(model* m, int vid, int i, int j, double* x, double* y)
{
    void* grid = m->grids[m->vars[vid].gridid];

    grid_ij2xy(grid, i, j, x, y);

    if (isnan(*x + *y))
        return STATUS_OUTSIDEGRID;
    return STATUS_OK;
}

/**
 */
int model_z2fk(model* m, int vid, double fi, double fj, double z, double* fk)
{
    void* grid = m->grids[m->vars[vid].gridid];
    int isperiodic_x = grid_isperiodic_x(grid);
    int** numlevels = grid_getnumlevels(grid);
    int ni, nj;
    int i1, i2, j1, j2, k2;

    if (isnan(fi + fj)) {
        *fk = NAN;
        return STATUS_OUTSIDEGRID;
    }

    grid_z2fk(grid, fi, fj, z, fk);

    if (isnan(*fk))
        return STATUS_OUTSIDEGRID;

    if (grid_getvtype(grid) == GRIDVTYPE_SIGMA || grid_getdepth(grid) == NULL)
         return STATUS_OK;

    /*
     * a depth check for z-grid:
     */
    model_getvardims(m, vid, &ni, &nj, NULL);
    i1 = floor(fi);
    i2 = ceil(fi);
    if (i1 == -1)
        i1 = (isperiodic_x) ? ni - 1 : i2;
    if (i2 == ni)
        i2 = (isperiodic_x) ? 0 : i1;
    j1 = floor(fj);
    j2 = ceil(fj);
    if (j1 == -1)
        j1 = j2;
    if (j2 == nj)
        j2 = j1;
    k2 = floor(*fk);
    if (numlevels[j1][i1] <= k2 && numlevels[j1][i2] <= k2 && numlevels[j2][i1] <= k2 && numlevels[j2][i2] <= k2) {
        *fk = NAN;
        return STATUS_LAND;
    } else if (numlevels[j1][i1] <= k2 || numlevels[j1][i2] <= k2 || numlevels[j2][i1] <= k2 || numlevels[j2][i2] <= k2) {
        float** depth = grid_getdepth(grid);
        int ni, nj;
        double v;

        grid_getdims(grid, &ni, &nj, NULL);

        v = interpolate2d(fi, fj, ni, nj, depth, numlevels, grid_isperiodic_x(grid));

        if (z > v)
            return STATUS_LAND;
    }

    return STATUS_OK;
}

/**
 */
int model_fk2z(model* m, int vid, int i, int j, double fk, double* z)
{
    grid* g = m->grids[m->vars[vid].gridid];
    float** depth;
    int ni, nj;

    grid_getdims(g, &ni, &nj, NULL);
    if (i < 0 || j < 0 || i >= ni || j >= nj) {
        *z = NAN;
        return STATUS_OUTSIDEGRID;
    }
    grid_fk2z(g, i, j, fk, z);
    depth = grid_getdepth(g);
    if (*z > depth[j][i]) {
        *z = NAN;
        return STATUS_OUTSIDEGRID;
    }

    return STATUS_OK;
}

/**
 */
void model_readfield(model* m, char fname[], int time, char varname[], int k, float* v)
{
    readfield(fname, varname, k, v);
}

/**
 */
void model_read3dfield(model* m, char fname[], int time, char varname[], float* v)
{
    read3dfield(fname, varname, v);
}

/**
 */
void model_writefield(model* m, char fname[], int time, char varname[], int k, float* v)
{
    writefield(fname, varname, k, v);
}

/**
 */
void model_randomisefield(model* m, int varid, float** v)
{
    float deflation = (float) m->vars[varid].deflation;
    float sigma = (float) m->vars[varid].sigma;
    float s;
    int ni, nj;
    int i, j;
    double tmp[2];

    get_normalpair(tmp);
    s = (float) (sqrt(1.0 - deflation * deflation) * tmp[0]) * sigma;

    model_getvardims(m, varid, &ni, &nj, NULL);
    for (j = 0; j < nj; ++j)
        for (i = 0; i < ni; ++i)
            v[j][i] = deflation * v[j][i] + s;
}
