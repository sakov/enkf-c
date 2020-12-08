/******************************************************************************
 *
 * File:        model.h        
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

#if !defined(_MODEL_H)

struct model;
typedef struct model model;

struct variable;
typedef struct variable variable;

model* model_create(enkfprm* prm);
void model_destroy(model* m);
void model_destroygrids(model* m);

#if defined(ENKF_CALC)
int model_destroygxytrees(model* m);
#endif
void model_print(model* m, char offset[]);
void model_describeprm(void);

void model_adddata(model* m, char tag[], int vid, int alloctype, void* data);
void model_addorreplacedata(model* m, char tag[], int vid, int alloctype, void* data);
void* model_getdata(model* m, char tag[]);
int model_getdataalloctype(model* m, char tag[]);

int model_getnvar(model* m);
char* model_getvarname(model* m, int varid);
int model_getvarid(model* m, char* varname, int hastosucceed);
void model_getvarinflation(model* m, int varid, float* inflation, double* inf_ratio);
double model_getvardeflation(model* m, int varid);
double model_getvarsigma(model* m, int varid);
void model_getvargridsize(model* m, int vid, int* ni, int* nj, int* nk);
void* model_getvargrid(model* m, int vid);
int model_getvargridid(model* m, int vid);
int model_getvarislog(model* m, int vid);
int model_getngrid(model* m);
void* model_getgridbyid(model* m, int gridid);
void* model_getgridbyname(model* m, char name[]);
double model_getlonbase(model* m, int vid);
float** model_getdepth(model* m, int vid, int musthave);
int** model_getnumlevels(model* m, int vid);
int model_getdomainid(model* m, char* domainname);

int model_xy2fij(model* m, int vid, double x, double y, double* fi, double* fj);
int model_fij2xy(model* m, int vid, double fi, double fj, double* x, double* y);
int model_ij2xy(model* m, int vid, int i, int j, double* x, double* y);
int model_z2fk(model* m, int vid, double fi, double fj, double z, double* fk);
int model_fk2z(model* m, int vid, int i, int j, double fk, double* z);
void model_readfield(model* m, char fname[], char varname[], int k, float* v);
void model_read3dfield(model* m, char fname[], char varname[], float* v);
void model_writefield(model* m, char fname[], char varname[], int k, float* v, int ignorelog);
void model_writefieldas(model* m, char fname[], char varname[], char varnameas[], int k, float* v, int ignorelog);
void model_addcustomdata(model* m, char* token, char* fname, int line);
void model_addcustomtaper(model* m, variable* var, char* token, char* fname, int line);

#define _MODEL_H
#endif
