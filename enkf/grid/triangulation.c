/******************************************************************************
 *
 * File:           triangulation.c
 *
 * Created:        25/01/2022
 *
 * Author:         Pavel Sakov
 *                 Bureau of Meteorology
 *
 * Purpose:        Triangulation handling code needed for unstructured grids.
 *
 * Description:    Note that there is no code for triangulation as such, which
 *                 must be done outside of EnKF-C. The triangulation code in
 *                 EnKF-C provides (1) triangle search by location (xy2i) and
 *                 (2) calculation of barycentric coordinates of a point
 *                 (xy2bc).
 *
 * Revisions:      
 *
 *****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#include <assert.h>
#include "ncw.h"
#include "triangulation.h"

/**
 */
static void quit(char* format, ...)
{
    va_list args;

    fflush(stdout);

    fprintf(stderr, "\n\n  ERROR: triangulation: ");
    va_start(args, format);
    vfprintf(stderr, format, args);
    va_end(args);
    fprintf(stderr, "\n\n");
    fflush(NULL);
    exit(1);
}

/**
 */
static triangulation* triangulation_create()
{
    triangulation* d = calloc(1, sizeof(triangulation));

    d->xmin = DBL_MAX;
    d->xmax = -DBL_MAX;
    d->ymin = DBL_MAX;
    d->ymax = -DBL_MAX;

    return d;
}

/**
 */
triangulation* triangulation_read(char* fname, char* name_x, char* name_y, char* name_tri, char* name_nei)
{
    triangulation* d;
    int ncid;
    int vid_x, vid_y, vid_tri, vid_nei, dimid_np, dimid_tri, tmp;
    size_t np, ntri, nnei;
    int i;

    void* storage = NULL;
    size_t storagesize = 0;

    ncw_open(fname, NC_NOWRITE, &ncid);

    ncw_inq_varid(ncid, name_x, &vid_x);
    ncw_check_varndims(ncid, vid_x, 1);
    ncw_inq_varid(ncid, name_y, &vid_y);
    ncw_check_varndims(ncid, vid_y, 1);
    ncw_inq_vardimid(ncid, vid_x, &dimid_np);
    ncw_inq_vardimid(ncid, vid_y, &tmp);
    assert(dimid_np == tmp);
    ncw_inq_dimlen(ncid, dimid_np, &np);

    ntri = 0;
    nnei = 0;
    if (name_tri != NULL) {
        ncw_inq_varid(ncid, name_tri, &vid_tri);
        ncw_check_varndims(ncid, vid_tri, 2);
        ncw_inq_vardimid(ncid, vid_tri, &dimid_tri);
        ncw_inq_dimlen(ncid, dimid_tri, &ntri);
    }
    if (name_nei != NULL)
        nnei = ntri;

    storagesize = np * sizeof(point) + ntri * sizeof(triangle) + nnei * sizeof(triangle_neighbours);

    d = triangulation_create();
    d->npoints = np;
#if defined(USE_SHMEM)
    {
        int ierror;

        assert(sizeof(MPI_Aint) == sizeof(size_t));
        if (sm_comm_rank == 0) {
            ierror = MPI_Win_allocate_shared(storagesize, sizeof(double), MPI_INFO_NULL, sm_comm, &storage, &d->sm_win_tridata);
            assert(ierror == MPI_SUCCESS);
        } else {
            MPI_Aint my_size;
            int disp_unit;

            ierror = MPI_Win_allocate_shared(0, sizeof(double), MPI_INFO_NULL, sm_comm, &storage, &d->sm_win_tridata);
            assert(ierror == MPI_SUCCESS);
            ierror = MPI_Win_shared_query(d->sm_win_tridata, 0, &my_size, &disp_unit, &storage);
            assert(ierror == MPI_SUCCESS);
            assert(my_size = storagesize);
        }
    }
#else
    storage = malloc(storagesize);
#endif
    d->points = storage;

#if defined(USE_SHMEM)
    if (sm_comm_rank == 0) {
#else
    {
#endif
        double* x = malloc(2 * np * sizeof(double));
        double* y = &x[np];

        ncw_get_var_double(ncid, vid_x, x);
        ncw_get_var_double(ncid, vid_y, y);
        for (i = 0; i < np; ++i) {
            d->points[i].x = x[i];
            d->points[i].y = y[i];
        }
        free(x);

        if (name_tri != NULL) {
            d->ntriangles = ntri;
            d->triangles = (triangle*) &d->points[np];
            ncw_get_var(ncid, vid_tri, d->triangles);

            if (name_nei != NULL) {
                ncw_inq_varid(ncid, name_nei, &vid_nei);
                d->neighbours = (triangle_neighbours*) &d->triangles[ntri];
                ncw_get_var(ncid, vid_nei, d->neighbours);
            }
        }
        ncw_close(ncid);
    }
#if defined(USE_SHMEM)
    MPI_Win_fence(0, d->sm_win_tridata);
    MPI_Barrier(sm_comm);
#endif

    d->xmin = DBL_MAX;
    d->xmax = -DBL_MAX;
    d->ymin = DBL_MAX;
    d->ymax = -DBL_MAX;
    for (i = 0; i < np; ++i) {
        point* p = &d->points[i];

        if (p->x < d->xmin)
            d->xmin = p->x;
        if (p->x > d->xmax)
            d->xmax = p->x;
        if (p->y < d->ymin)
            d->ymin = p->y;
        if (p->y > d->ymax)
            d->ymax = p->y;
    }

    return d;
}

/**
 */
void triangulation_destroy(triangulation* d)
{
    if (d == NULL)
        return;

#if !defined(USE_SHMEM)
    if (d->triangles != NULL)
        free(d->points);
#else
    MPI_Win_free(&d->sm_win_tridata);
#endif
    free(d);
}

/**
 */
void triangulation_getpoints(triangulation* d, int* npoints, point** points)
{
    *npoints = d->npoints;
    *points = d->points;
}

/**
 */
point* triangulation_getpoint(triangulation* d, int i)
{
    if (i < 0 || i >= d->npoints)
        return NULL;
    return &d->points[i];
}

/** Returns whether the point p is on the right side of the vector (p0, p1).
 */
static int onrightside(point* p, point* p0, point* p1)
{
    return (p1->x - p->x) * (p0->y - p->y) > (p0->x - p->x) * (p1->y - p->y);
}

/* Finds triangle specified point belongs to (if any).
 *
 * @param d Triangulation triangulation
 * @param p Point to be mapped
 * @param seed Triangle index to start with
 * @return Triangle id if successful, -1 otherwhile
 */
int triangulation_xy2i(triangulation* d, point* p)
{
    int id = d->seed;
    triangle* t;
    int i;

    if (d->triangles == NULL || d->neighbours == NULL)
        quit("triangulation_xy2i(): incomplete structure");

    if (p->x < d->xmin || p->x > d->xmax || p->y < d->ymin || p->y > d->ymax)
        return -1;

    t = &d->triangles[id];
    do {
        for (i = 0; i < 3; ++i) {
            int i1 = (i + 1) % 3;

            if (onrightside(p, &d->points[t->vids[i]], &d->points[t->vids[i1]])) {
                id = d->neighbours[id].tids[(i + 2) % 3];
                if (id < 0)
                    return id;
                t = &d->triangles[id];
                break;
            }
        }
    } while (i < 3);

    d->seed = id;

    return id;
}

/** Calculate barycentric coordinates of a point.
 ** ids[3] and bcs[2] should be pre-allocated.
 * @param d - triangulation structure
 * @param p - point
 * @param ids - (output) [3] indices of the vertices of the triangle the point
 *              belongs to
 * @param bcs - (output) [2] barycentric coordinates for the first two vertices
 * @return 1 on success, 0 on failure
 */
int triangulation_xy2bc(triangulation* d, point* p, int* tid, double* bcs)
{
    if (d->triangles == NULL || d->neighbours == NULL)
        quit("triangulation_xy2bc(): incomplete structure");

    *tid = triangulation_xy2i(d, p);
    if (*tid < 0) {
        bcs[0] = NAN;
        bcs[1] = NAN;

        return 0;
    }

    {
        triangle* t = &d->triangles[*tid];
        point* p0 = &d->points[t->vids[0]];
        point* p1 = &d->points[t->vids[1]];
        point* p2 = &d->points[t->vids[2]];
        double det = (p0->x - p2->x) * (p1->y - p2->y) - (p1->x - p2->x) * (p0->y - p2->y);

        bcs[0] = ((p->x - p2->x) * (p1->y - p2->y) - (p->y - p2->y) * (p1->x - p2->x)) / det;
        bcs[1] = ((p->x - p2->x) * (p2->y - p0->y) - (p->y - p2->y) * (p2->x - p0->x)) / det;
        bcs[2] = 1.0 - bcs[0] - bcs[1];
    }

    return 1;
}

/**
 */
void triangulation_getminmax(triangulation* d, double* xmin, double* xmax, double* ymin, double* ymax)
{
    if (xmin != NULL)
        *xmin = d->xmin;
    if (xmax != NULL)
        *xmax = d->xmax;
    if (ymin != NULL)
        *ymin = d->ymin;
    if (ymax != NULL)
        *ymax = d->ymax;
}
