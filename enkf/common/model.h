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

#include "enkfprm.h"

struct model;
typedef struct model model;

struct variable;
typedef struct variable variable;

typedef void (*model_getmemberfname_fn) (model* m, char dir[], char varname[], int mem, char fname[]);
typedef int (*model_getmemberfnameasync_fn) (model* m, char dir[], char varname[], char otname[], int mem, int t, char fname[]);
typedef void (*model_getbgfname_fn) (model* m, char dir[], char varname[], char fname[]);
typedef int (*model_getbgfnameasync_fn) (model* m, char dir[], char varname[], char otnamename[], int t, char fname[]);
typedef void (*model_readfield_fn) (model* m, char fname[], int mem, int time, char varname[], int k, float* v);
typedef void (*model_read3dfield_fn) (model* m, char fname[], int mem, int time, char varname[], float* v);
typedef void (*model_writefield_fn) (model* m, char fname[], int time, char varname[], int k, float* v);
typedef void (*model_adddata_fn) (model* m, char* token, char* fname, int line);

model* model_create(enkfprm* prm);
void model_destroy(model* m);
void model_print(model* m, char offset[]);
void model_describeprm(void);

void model_setgetmemberfname_fn(model* m, model_getmemberfname_fn fn);
void model_setgetmemberfnameasync_fn(model* m, model_getmemberfnameasync_fn fn);
void model_setbgfname_fn(model* m, model_getbgfname_fn fn);
void model_setbgfnameasync_fn(model* m, model_getbgfnameasync_fn fn);
void model_setreadfield_fn(model* m, model_readfield_fn fn);
void model_setread3dfield_fn(model* m, model_read3dfield_fn fn);
void model_setwritefield_fn(model* m, model_writefield_fn fn);
void model_setadddata_fn(model* m, model_adddata_fn fn);

void model_setgrid(model* m, void* g);

void model_addmodeldata(model* m, char tag[], int alloctype, void* data);
void* model_getmodeldata(model* m, char tag[]);

int model_getnvar(model* m);
char* model_getvarname(model* m, int varid);
int model_getvarid(model* m, char* varname);
float model_getvarinflation(model* m, int varid);
void model_getvardims(model* m, int vid, int* ni, int* nj, int* nk);
void* model_getvargrid(model* m, int vid);
int model_getvargridid(model* m, int vid);
int model_getngrid(model* m);
void* model_getgridbyid(model* m, int gridid);
void* model_getgridbyname(model* m, char name[]);
int model_getlontype(model* m, int vid);
float** model_getdepth(model* m, int vid);
int** model_getnumlevels(model* m, int vid);
void model_getmemberfname(model* m, char dir[], char varname[], int mem, char fname[]);
int model_getmemberfname_async(model* m, char dir[], char varname[], char otname[], int mem, int time, char fname[]);
void model_getbgfname(model* m, char dir[], char varname[], char fname[]);
void model_getspreadfname(model* m, char varname[], char fname[]);
int model_getbgfname_async(model* m, char dir[], char varname[], char otname[], int time, char fname[]);

int model_xy2fij(model* m, int vid, double x, double y, double* fi, double* fj);
int model_fij2xy(model* m, int vid, double fi, double fj, double* x, double* y);
int model_z2fk(model* m, int vid, double fi, double fj, double z, double* fk);
void model_readfield(model* m, char fname[], int mem, int time, char varname[], int k, float* v);
void model_read3dfield(model* m, char fname[], int mem, int time, char varname[], float* v);
void model_writefield(model* m, char fname[], int time, char varname[], int k, float* v);
void model_addcustomdata(model* m, char* token, char* fname, int line);
void model_addcustomtaper(model* m, variable * var, char* token, char* fname, int line);

#define _MODEL_H
#endif
