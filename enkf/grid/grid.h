/******************************************************************************
 *
 * File:        grid.h        
 *
 * Created:     12/2012
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Purpose:     Header for `grid' object (generic grid).
 *
 * Description:
 *
 * Revisions:
 *
 *****************************************************************************/

#if !defined(_GRID_H)
#include "kdtree.h"

#define GRIDHTYPE_UNDEFINED 0
/*
 * (1) GRIDHTYPE_1D and GRIDHTYPE_2D are generic types used in UPDATE to save
 *     resources when creating actual grid objects is not required.
 * (2) GRIDHTYPE_1D and GRIDHTYPE_UNSTRUCTURED are used for unstructured 2D
 *     grids; (at the moment) for these 2 types there must be
 *     grid->hgrid->nj = 0, which may be used as a test for a grid being
 *     structured.
 */
#define GRIDHTYPE_1D 1
#define GRIDHTYPE_2D 2
#define GRIDHTYPE_RECTANGULAR 3
#define GRIDHTYPE_CURVILINEAR 4
#define GRIDHTYPE_UNSTRUCTURED 5

#define GRIDVTYPE_UNDEFINED 0
/*
 * GRIDVTYPE_NONE is supposed to be used for 2D variables on grids without 
 * vertical structure (to avoid creating Z dimension of size 1). It needs
 * to be tested.
 */
#define GRIDVTYPE_NONE 1
#define GRIDVTYPE_Z 2
#define GRIDVTYPE_SIGMA 3
#define GRIDVTYPE_HYBRID 4
#define GRIDVTYPE_NUMERIC 5

struct grid;
typedef struct grid grid;

grid* grid_create(void* prm, int id, void** grids);
void grid_destroy(grid* g);
grid* setdepth(grid* g, void* depth);

void* grid_getvgrid(grid* g);
int grid_gethtype(grid* g);
int grid_getvtype(grid* g);
void* grid_getdepth(grid* g);

void grid_print(grid* g, char offset[]);
void grid_describeprm(void);

void grid_getsize(grid* g, int* ni, int* nj, int* nk);
int grid_getsurflayerid(grid* g);
char* grid_getname(grid* g);
int grid_getid(grid* g);
int grid_gethtype(grid* g);
int grid_getvtype(grid* g);
void* grid_getnumlevels(grid* g);
double grid_getlonbase(grid* g);
int grid_getstride(grid* g);
void grid_getzints(grid* g, int* nzints, zint* zints[]);
char* grid_getdomainname(grid* g);
int grid_getaliasid(grid* g);

#if defined(ENKF_PREP) || defined (ENKF_CALC)
int grid_xy2fij(grid* g, double x, double y, double* fij);
void grid_ij2xy(grid* g, int* ij, double* x, double* y);
int grid_fk2z(grid* g, int* ij, double fk, double* z);
int grid_island(grid* g, double* fij, double fk);
#endif
int grid_z2fk(grid* g, double* fij, double z, double* fk);
int grid_z2fk_f(grid* g, double* fij, double z, float* fk);
int grid_isperiodic_i(grid* g);
int grid_isstructured(grid* g);
int grid_isgeographic(grid* g);
float grid_interpolate2d(grid* g, double* fij, void* v);
float grid_interpolate3d(grid* g, double* fij, float fk, void* v);

#if defined(ENKF_CALC)
kdtree* grid_gettreeXYZ(grid* g, int createifnull);
void grid_destroytreeXYZ(grid* g);
#endif

/*
 * stuff to handle an array of grids
 */
void grids_create(enkfprm* prm, int* ngrid, void*** grids);
void grids_destroy(int ngrid, void** grids);

#if defined(ENKF_CALC)
int grids_destroyhtrees(int ngrid, void** grids);
#endif

#define _GRID_H
#endif
