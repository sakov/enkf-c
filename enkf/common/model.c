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
#include "stringtable.h"
#include "definitions.h"
#include "utils.h"
#include "ncutils.h"
#include "grid.h"
#include "gridprm.h"
#include "enkfprm.h"
#include "model.h"

#define NMODELDATA_INC 10
#define NVAR_INC 10

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
    int applylog;
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

    stringtable* domains;
};

/**
 */
static void variable_new(variable* v, int id, char* name)
{
    v->id = id;
    v->name = strdup(name);
    v->gridid = -1;
    v->inflation = NAN;
    v->inf_ratio = NAN;
    v->applylog = 0;
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
    int i;

    grids_create(prm, &m->ngrid, &m->grids);

    assert(m->ngrid > 0);

    m->domains = st_create("domains");
    for (i = 0; i < m->ngrid; ++i)
        st_add_ifabsent(m->domains, grid_getdomainname(m->grids[i]), -1);

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
                    if (strcmp(token, grid_getname(m->grids[i])) == 0) {
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
                    enkf_quit("%s, l.%d: INFLATION value not specified", modelprm, line);
                if (!str2double(token, &now->inflation))
                    enkf_quit("%s, l.%d: could not convert \"%s\" to double", modelprm, line, token);
                if ((token = strtok(NULL, seps)) != NULL) {
                    if (!str2double(token, &now->inf_ratio))
                        enkf_quit("%s, l.%d: could not convert \"%s\" to double", modelprm, line, token);
                }
            } else if (strcasecmp(token, "APPLYLOG") == 0) {
                if (now == NULL)
                    enkf_quit("%s, l.%d: VAR not specified", modelprm, line);
                if ((token = strtok(NULL, seps)) == NULL)
                    enkf_quit("%s, l.%d: APPLYLOG value not specified", modelprm, line);
                now->applylog = istrue(token) ? 1 : 0;
                if (now->applylog && prm->mode != MODE_ENKF && !enkf_allowenoilog)

                    enkf_quit("%s, l.%d: to use APPLYLOG in EnOI or Hybrid modes the static ensemble must be in log space. If you are aware of this and wish to proceed use command line option \"--allow-logspace-with-static-ens\"", modelprm, line);
            } else if (strcasecmp(token, "RANDOMISE") == 0 || strcasecmp(token, "RANDOMIZE") == 0) {
                if (!isnan(now->deflation))
                    enkf_quit("%s, l.%d: randomisation multiple already specified for \"%s\"", modelprm, line, now->name);
                if ((token = strtok(NULL, seps)) == NULL)
                    enkf_quit("%s, l.%d: RANDOMISE parameters not specified", modelprm, line);
                if (!str2double(token, &now->deflation))
                    enkf_quit("%s, l.%d: could not convert \"%s\" to double", modelprm, line, token);
                if ((token = strtok(NULL, seps)) == NULL)
                    enkf_quit("%s, l.%d: RANDOMISE STD value not specified", modelprm, line);
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
                m->vars[i].inflation = prm->inflation;
                m->vars[i].inf_ratio = prm->inf_ratio;
            }
        prm->inflation = NAN;
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

    model_checkvars(m, modelprm);

    model_print(m, "  ");

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
void model_destroygrids(model* m)
{
    grids_destroy(m->ngrid, m->grids);
    m->ngrid = 0;
}

/**
 */
void model_destroy(model* m)
{
    free(m->name);
    model_destroygrids(m);
    model_destroyvars(m);
    model_freemodeldata(m);
    st_destroy(m->domains);
    free(m);
}

#if defined(ENKF_CALC)
int model_destroygxytrees(model* m)
{
    return grids_destroyhtrees(m->ngrid, m->grids);
}
#endif

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
        if (v->applylog)
            enkf_printf("%s      applylog = true\n", offset);
        if (!isnan(v->deflation))
            enkf_printf("%s      randomise: deflation = %.3f, sigma = %.3f\n", offset, v->deflation, v->sigma);
    }
    if (m->ndata > 0) {
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
    enkf_printf("  [ APPLYLOG  = <YES | NO*> ]\n");
    enkf_printf("  [ RANDOMISE <deflation> <sigma> ]\n");
    enkf_printf("\n");
    enkf_printf("  [ <more of the above blocks> ]\n");
    enkf_printf("\n");
    enkf_printf("  Notes:\n");
    enkf_printf("    1. [ ... ] denotes an optional input\n");
    enkf_printf("    2. < ... > denotes a description of an entry\n");
    enkf_printf("    3. ( ... ) is a note\n");
    enkf_printf("    4. * denotes the default value\n");
    enkf_printf("\n");
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

    if (i == m->ndata) {
        if (m->ndata % NMODELDATA_INC == 0)
            m->data = realloc(m->data, (m->ndata + NMODELDATA_INC) * sizeof(modeldata));
        mdata = &m->data[i];
        mdata->tag = strdup(tag);
        mdata->vid = vid;
        mdata->alloctype = alloctype;
        m->ndata++;
    } else {
        mdata = &m->data[i];
        assert(m->data[i].alloctype == alloctype);
    }

    model_getvargridsize(m, mdata->vid, &ni, &nj, &nk);
    if (mdata->alloctype == ALLOCTYPE_1D)
        memcpy(mdata->data, data, nk * sizeof(float));
    else if (mdata->alloctype == ALLOCTYPE_2D)
        mdata->data = copy2d(data, nj, ni, sizeof(float));
    else if (mdata->alloctype == ALLOCTYPE_3D)
        mdata->data = copy3d(data, nk, nj, ni, sizeof(float));
    else
        enkf_quit("programming error");
}

/** Get model data.
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
double model_getvarsigma(model* m, int varid)
{
    return m->vars[varid].sigma;
}

/**
 */
void model_getvargridsize(model* m, int vid, int* ni, int* nj, int* nk)
{
    grid_getsize(m->grids[m->vars[vid].gridid], ni, nj, nk);
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
int model_getvarislog(model* m, int vid)
{
    return m->vars[vid].applylog;
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
    void* g = m->grids[m->vars[vid].gridid];
    float** depth = grid_getdepth(g);

    if (musthave && depth == NULL)
        enkf_quit("DEPTHVARNAME not specified for grid \"%s\"", grid_getname(g));

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
int model_getdomainid(model* m, char* domainname)
{
    return st_findindexbystring(m->domains, domainname);
}

/**
 */
int model_hasunstructured(model* m)
{
    int i;

    for (i = 0; i < m->ngrid; ++i)
        if (!grid_isstructured(m->grids[i]))
            return 1;

    return 0;
}

/**
 */
#if defined(ENKF_PREP) || defined(ENKF_CALC)
int model_xy2fij(model* m, int vid, double x, double y, double* fij)
{
    return grid_xy2fij(m->grids[m->vars[vid].gridid], x, y, fij);
}
#endif

/**
 */
#if defined(ENKF_PREP) || defined(ENKF_CALC)
int model_ij2xy(model* m, int vid, int* ij, double* x, double* y)
{
    void* g = m->grids[m->vars[vid].gridid];

    grid_ij2xy(g, ij, x, y);

    if (isnan(*x + *y))
        return STATUS_OUTSIDEGRID;
    return STATUS_OK;
}
#endif

/**
 */
#if defined(ENKF_PREP) || defined(ENKF_CALC)
int model_z2fk(model* m, int vid, double* fij, double z, double* fk)
{
    return grid_z2fk(m->grids[m->vars[vid].gridid], fij, z, fk);
}
#endif

/**
 */
#if defined(ENKF_PREP) || defined(ENKF_CALC)
int model_z2fk_f(model* m, int vid, double* fij, double z, float* fk)
{
    return grid_z2fk_f(m->grids[m->vars[vid].gridid], fij, z, fk);
}
#endif

/**
 */
#if defined(ENKF_PREP) || defined(ENKF_CALC)
double model_fk2z(model* m, int vid, int* ij, double fk)
{
    double z;

    (void) grid_fk2z(m->grids[m->vars[vid].gridid], ij, fk, &z);

    return z;
}
#endif

/**
 */
void model_readfield(model* m, char fname[], char varname[], int k, float* v, int ignorelog)
{
    int ni, nj, nk;
    int mvid = model_getvarid(m, varname, 1);

    model_getvargridsize(m, mvid, &ni, &nj, &nk);
    assert(k < nk);
    ncu_readfield(fname, varname, k, ni, nj, nk, v);

    if (m->vars[mvid].applylog && !ignorelog) {
        size_t nij = ni * nj;
        size_t i;

        for (i = 0; i < nij; ++i)
            v[i] = log10(v[i]);
    }
}

/**
 */
void model_read3dfield(model* m, char fname[], char varname[], float* v, int ignorelog)
{
    int ni, nj, nk;
    int mvid = model_getvarid(m, varname, 1);

    model_getvargridsize(m, mvid, &ni, &nj, &nk);
    ncu_read3dfield(fname, varname, ni, nj, nk, v);

    if (m->vars[mvid].applylog && !ignorelog) {
        size_t nijk = (size_t) ((nj > 0) ? ni * nj * nk : ni * nk);
        size_t i;

        for (i = 0; i < nijk; ++i)
            v[i] = log10(v[i]);
    }
}

/**
 */
void model_writefield(model* m, char fname[], char varname[], int k, float* v, int ignorelog)
{
    int ni, nj, nk;
    int mvid = model_getvarid(m, varname, 1);

    model_getvargridsize(m, mvid, &ni, &nj, &nk);
    assert(k < nk);

    if (m->vars[mvid].applylog && !ignorelog) {
        size_t nij = (size_t) ((nj > 0) ? ni * nj : ni);
        size_t i;

        for (i = 0; i < nij; ++i)
            v[i] = exp10(v[i]);
    }

    ncu_writefield(fname, varname, (k >= 0) ? k : 0, ni, nj, (k >= 0) ? nk : 1, v);
}

/**
 */
void model_writefieldas(model* m, char fname[], char varname[], char varnameas[], int k, float* v, int ignorelog)
{
    int ni, nj, nk;
    int mvid = model_getvarid(m, varnameas, 1);

    model_getvargridsize(m, mvid, &ni, &nj, &nk);
    assert(k < nk);

    if (m->vars[mvid].applylog && !ignorelog) {
        size_t nij = (size_t) ((nj > 0) ? ni * nj : ni);
        size_t i;

        for (i = 0; i < nij; ++i)
            v[i] = exp10(v[i]);
    }

    ncu_writefield(fname, varname, (k >= 0) ? k : 0, ni, nj, (k >= 0) ? nk : 1, v);
}
