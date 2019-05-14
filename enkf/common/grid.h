/******************************************************************************
 *
 * File:        grid.h        
 *
 * Created:     12/2012
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description:
 *
 * Revisions:   10/12/2014 PS: moved here gridprm stuff from z-model and 
 *                sigma-model.c; some reshuffling of mapping functions.
 *
 *****************************************************************************/

#if !defined(_GRID_H)
#include "kdtree.h"

#define GRIDHTYPE_NONE 0
#define GRIDHTYPE_LATLON 1
#if !defined(NO_GRIDUTILS)
#define GRIDHTYPE_CURVILINEAR 2
#endif

#define GRIDVTYPE_NONE 0
#define GRIDVTYPE_Z 1
#define GRIDVTYPE_SIGMA 2
#define GRIDVTYPE_HYBRID 3

typedef struct {
    double z1, z2;
} zint;

struct grid;
typedef struct grid grid;

grid* grid_create(void* prm, int id);
void grid_destroy(grid* g);
void grid_print(grid* g, char offset[]);
void grid_describeprm(void);

void grid_getsize(grid* g, int* ni, int* nj, int* nk);
int grid_getsurflayerid(grid* g);
char* grid_getname(grid* g);
int grid_getid(grid* g);
int grid_gethtype(grid* g);
int grid_getvtype(grid* g);
float** grid_getdepth(grid* g);
int** grid_getnumlevels(grid* g);
double grid_getlonbase(grid* g);
int grid_getstride(grid* g);
void grid_setstride(grid* g, int stride);
double grid_getsfactor(grid* g);
void grid_getzints(grid* g, int* nzints, zint* zints[]);
char* grid_getdomainname(grid* g);

int grid_xy2fij(grid* g, double x, double y, double* fi, double* fj);
int grid_z2fk(grid* g, double fi, double fj, double z, double* fk);
void grid_fij2xy(grid* g, double fi, double fj, double* x, double* y);
void grid_ij2xy(grid* g, int i, int j, double* x, double* y);
int grid_fk2z(grid* g, int i, int j, double fk, double* z);

#if defined(ENKF_CALC)
kdtree* grid_gettree(grid* g);
#endif
int grid_isperiodic_i(grid* g);

/*
 * stuff to handle an array of grids
 */
void grids_create(char gprmfname[], int stride, int* ngrid, void*** grids);
void grids_destroy(int ngrid, void** grids);

#define _GRID_H
#endif
