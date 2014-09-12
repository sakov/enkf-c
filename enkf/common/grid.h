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
 * Revisions:
 *
 *****************************************************************************/

#if !defined(_GRID_H)

#define GRIDHTYPE_NONE 0
#define GRIDHTYPE_LATLON_REGULAR 1
#define GRIDHTYPE_LATLON_IRREGULAR 2
#define GRIDHTYPE_CURVILINEAR 3

typedef void (*grid_xy2fij_fn) (void* grid, double x, double y, double* fi, double* fj);
typedef void (*grid_z2fk_fn) (void* grid, double fi, double fj, double z, double* fk);
typedef void (*grid_fij2xy_fn) (void* grid, double fi, double fj, double* x, double* y);
typedef void (*grid_tocartesian_fn) (double* in, double* out);

struct grid;
typedef struct grid grid;

grid* grid_create(char name[]);
void grid_destroy(grid* g);
void grid_print(grid* g, char offset[]);
void grid_describeprm(void);

void grid_setcoords(grid* g, int type, int periodic_x, int periodic_y, int nx, int ny, int nz, void* x, void* y, double* z);
void grid_setdepth(grid* g, float** depth);
void grid_setnumlevels(grid* g, int** numlevels);
void grid_settocartesian_fn(grid* g, grid_tocartesian_fn fn);

void grid_getdims(grid* g, int* ni, int* nj, int* nk);
char* grid_getname(grid* g);
int grid_gethtype(grid* g);
float** grid_getdepth(grid* g);
int** grid_getnumlevels(grid* g);
int grid_getlontype(grid* g);
grid_xy2fij_fn grid_getxy2fijfn(grid* g);
grid_z2fk_fn grid_getz2fkfn(grid* g);
grid_fij2xy_fn grid_getfij2xyfn(grid* g);
int grid_isperiodic_x(grid* g);
int grid_isperiodic_y(grid* g);
void grid_tocartesian(grid* g, double* in, double* out);

#define _GRID_H
#endif
