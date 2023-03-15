/******************************************************************************
 *
 * File:           triangulation.h
 *
 * Created:        25/01/2022
 *
 * Author:         Pavel Sakov
 *                 Bureau of Meteorology
 *
 * Purpose:        Header for triangulation related code.
 *
 * Description:    The triangulation code in EnKF-C provides (1) triangle
 *                 search by location (xy2i) and (2) calculation of barycentric
 *                 coordinates of a point (xy2bc).
 *
 * Revisions:      
 *
 *****************************************************************************/

#if !defined(_TRIANGULATION_H)
#define _TRIANGULATION_H

#if defined(USE_SHMEM)
#include <mpi.h>
extern MPI_Comm sm_comm;
extern int sm_comm_rank;
#endif

typedef struct {
    double x;
    double y;
} point;

typedef struct {
    int vids[3];
} triangle;

typedef struct {
    int tids[3];
} triangle_neighbours;

typedef struct {
    int npoints;
    point* points;
    double xmin;
    double xmax;
    double ymin;
    double ymax;
    int seed;

    int ntriangles;
    triangle* triangles;
    triangle_neighbours* neighbours;    /* for triangulation_xy2i() */

#if defined(USE_SHMEM)
    MPI_Win sm_win_tridata;
#endif
} triangulation;

/* It is possible to set the quit procedure. By default, the internal procedure
 * is used.
 */
typedef void (*triangulation_quit_fn) (char* format, ...);
void triangulation_set_quitfn(triangulation_quit_fn quit_fn);

triangulation* triangulation_read(char* fname, char* name_x, char* name_y, char* name_tri, char* name_neigh);
void triangulation_destroy(triangulation* d);
void triangulation_getpoints(triangulation* d, int* npoint, point** points);
point* triangulation_getpoint(triangulation* d, int i);
int triangulation_xy2i(triangulation* d, point* p);
int triangulation_xy2bc(triangulation* d, point* p, int* tid, double* bc);
void triangulation_getminmax(triangulation* d, double* xmin, double* xmax, double* ymin, double* ymax);

#endif
