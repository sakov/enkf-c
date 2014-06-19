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

typedef void (*model_getmemberfname_fn) (model* m, char dir[], char varname[], int mem, char fname[]);
typedef int (*model_getmemberfnameasync_fn) (model* m, char dir[], char varname[], char otname[], int mem, int t, char fname[]);
typedef void (*model_getbgfname_fn) (model* m, char dir[], char varname[], char fname[]);
typedef int (*model_getbgfnameasync_fn) (model* m, char dir[], char varname[], char otnamename[], int t, char fname[]);
typedef void (*model_readfield_fn) (model* m, char fname[], int mem, int time, char varname[], int k, float* v);
typedef void (*model_read3dfield_fn) (model* m, char fname[], int mem, int time, char varname[], float* v);
typedef void (*model_writefield_fn) (model* m, char fname[], int time, char varname[], int k, float* v);

model* model_create(enkfprm* prm);
void model_destroy(model* m);

void model_setdepth(model* m, float** depth);
void model_setnumlevels(model* m, int** numlevels);
void model_setgetmemberfname_fn(model* m, model_getmemberfname_fn fn);
void model_setgetmemberfnameasync_fn(model* m, model_getmemberfnameasync_fn fn);
void model_setbgfname_fn(model* m, model_getbgfname_fn fn);
void model_setbgfnameasync_fn(model* m, model_getbgfnameasync_fn fn);
void model_setreadfield_fn(model* m, model_readfield_fn fn);
void model_setread3dfield_fn(model* m, model_read3dfield_fn fn);
void model_setwritefield_fn(model* m, model_writefield_fn fn);

void model_getdims(model* m, int* ni, int* nj, int* nk);
void* model_getgrid(model* m);
int model_getlontype(model* m);
float** model_getdepth(model* m);
int** model_getnumlevels(model* m);
void model_getmemberfname(model* m, char dir[], char varname[], int mem, char fname[]);
int model_getmemberfname_async(model* m, char dir[], char varname[], char otname[], int mem, int time, char fname[]);
void model_getbgfname(model* m, char dir[], char varname[], char fname[]);
void model_getspreadfname(model* m, char varname[], char fname[]);
int model_getbgfname_async(model* m, char dir[], char varname[], char otname[], int time, char fname[]);

int model_ll2fij(model* m, double lon, double lat, double* fi, double* fj);
int model_fij2ll(model* m, double fi, double fj, double* lon, double* lat);
int model_z2fk(model* m, double fi, double fj, double z, double* fk);
void model_readfield(model* m, char fname[], int mem, int time, char varname[], int k, float* v);
void model_read3dfield(model* m, char fname[], int mem, int time, char varname[], float* v);
void model_writefield(model* m, char fname[], int time, char varname[], int k, float* v);

#define _MODEL_H
#endif
