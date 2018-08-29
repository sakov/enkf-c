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

typedef void (*grid_tocartesian_fn) (double in[2], double out[3]);

grid* grid_create(void* prm, int id);
void grid_destroy(grid* g);
void grid_print(grid* g, char offset[]);
void grid_describeprm(void);
void grid_settocartesian_fn(grid* g, grid_tocartesian_fn fn);

void grid_getdims(grid* g, int* ni, int* nj, int* nk);
int grid_getsurflayerid(grid* g);
char* grid_getname(grid* g);
int grid_getid(grid* g);
int grid_gethtype(grid* g);
int grid_getvtype(grid* g);
float** grid_getdepth(grid* g);
float** grid_getangle(grid* g);
int** grid_getnumlevels(grid* g);
double grid_getlonbase(grid* g);
int grid_getstride(grid* g);
void grid_setstride(grid* g, int stride);
int grid_getsobstride(grid* g);
void grid_setsobstride(grid* g, int sobstride);
double grid_getsfactor(grid* g);
void grid_getzints(grid* g, int* nzints, zint* zints[]);

int grid_xy2fij(grid* g, double x, double y, double* fi, double* fj);
int grid_z2fk(grid* g, double fi, double fj, double z, double* fk);
void grid_fij2xy(grid* g, double fi, double fj, double* x, double* y);
void grid_ij2xy(grid* g, int i, int j, double* x, double* y);
int grid_fk2z(grid* g, int i, int j, double fk, double* z);
void grid_tocartesian(grid* g, double in[2], double out[3]);

int grid_isperiodic_i(grid* g);

/*
 * stuff to handle an array of grids
 */
void grids_create(char gprmfname[], int stride, int sob_stride, int* ngrid, void*** grids);
void grids_destroy(int ngrid, void** grids);

#define _GRID_H
#endif
