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

#if defined(__GNUC__)
static void quit_def(char* format, ...) __attribute__ ((format(printf, 1, 2)));
#else
static void quit_def(char* format, ...);
#endif
static triangulation_quit_fn quit = quit_def;

/**
 */
static void quit_def(char* format, ...)
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
void triangulation_set_quitfn(triangulation_quit_fn quit_fn)
{
    quit = quit_fn;
}

/**
 */
static triangulation* triangulation_create(void)
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
static void check_inttype(int ncid, int varid, char* fname, char* varname)
{
    nc_type type;

    ncw_inq_vartype(ncid, varid, &type);
    if (type != NC_BYTE && type != NC_SHORT && type != NC_INT && type != NC_UBYTE && type != NC_USHORT && type != NC_UINT && type != NC_INT64 && type != NC_UINT64)
        quit("%s: \"%s\": must be of integer type", fname, varname);
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
    if (dimid_np != tmp)
        quit("%s: \"%s\" and \"%s\" must have the same dimension", fname, name_x, name_y);
    ncw_inq_dimlen(ncid, dimid_np, &np);

    ntri = 0;
    nnei = 0;
    if (name_tri != NULL) {
        size_t dimlen[2];

        ncw_inq_varid(ncid, name_tri, &vid_tri);
        ncw_check_varndims(ncid, vid_tri, 2);
        ncw_inq_vardimid(ncid, vid_tri, &dimid_tri);
        ncw_inq_dimlen(ncid, dimid_tri, &ntri);
        dimlen[0] = ntri;
        dimlen[1] = 3;
        ncw_check_vardims(ncid, vid_tri, 2, dimlen);
        check_inttype(ncid, vid_tri, fname, name_tri);
        if (name_nei != NULL) {
            ncw_inq_varid(ncid, name_nei, &vid_nei);
            ncw_check_vardims(ncid, vid_nei, 2, dimlen);
            check_inttype(ncid, vid_nei, fname, name_nei);
            nnei = ntri;
        }
    }

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
            assert(my_size == storagesize);
        }
    }
    MPI_Win_fence(0, d->sm_win_tridata);
#else
    storage = malloc(storagesize);
#endif

    d->points = storage;
#if defined(USE_SHMEM)
    if (sm_comm_rank == 0)
#endif
    {
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
            d->triangles = (triangle*) &d->points[np];
            ncw_get_var_int(ncid, vid_tri, (int*) d->triangles);

            if (name_nei != NULL) {
                d->neighbours = (triangle_neighbours*) &d->triangles[ntri];
                ncw_get_var_int(ncid, vid_nei, (int*) d->neighbours);
            }
        }
    }
#if defined(USE_SHMEM)
    else if (name_tri != NULL) {
        d->triangles = (triangle*) &d->points[np];
        if (name_nei != NULL)
            d->neighbours = (triangle_neighbours*) &d->triangles[ntri];
    }
#endif
    if (name_tri != NULL)
        d->ntriangles = ntri;
    ncw_close(ncid);
#if defined(USE_SHMEM)
    MPI_Win_fence(0, d->sm_win_tridata);
#endif

    /*
     * check index ranges (indices must be 0-based; negative neighbour ids
     * indicate boundary edges)
     */
    if (d->triangles != NULL) {
        int vmin = INT_MAX, vmax = INT_MIN;

        for (i = 0; i < (int) ntri; ++i) {
            int k;

            for (k = 0; k < 3; ++k) {
                int v = d->triangles[i].vids[k];

                if (v < vmin)
                    vmin = v;
                if (v > vmax)
                    vmax = v;
            }
        }
        if (ntri > 0 && (vmin < 0 || vmax >= (int) np))
            quit("%s: \"%s\": vertex ids must be within [0, %zu] (0-based), found [%d, %d]%s", fname, name_tri, np - 1, vmin, vmax, (vmin == 1 && vmax == (int) np) ? " (1-based?)" : "");
        if (d->neighbours != NULL) {
            vmax = INT_MIN;
            for (i = 0; i < (int) ntri; ++i) {
                int k;

                for (k = 0; k < 3; ++k)
                    if (d->neighbours[i].tids[k] > vmax)
                        vmax = d->neighbours[i].tids[k];
            }
            if (ntri > 0 && vmax >= (int) ntri)
                quit("%s: \"%s\": triangle ids must be < %zu (0-based, negative for no neighbour), found %d", fname, name_nei, ntri, vmax);
        }
    }

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
    free(d->points);
#else
    MPI_Win_free(&d->sm_win_tridata);
#endif
    free(d->bstart);
    free(d->btris);
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

/** Returns whether the point p is inside (or on the boundary of) triangle t
 ** (assuming the same vertex orientation as the walk in triangulation_xy2i()).
 */
static int intriangle(triangulation* d, point* p, triangle* t)
{
    int i;

    for (i = 0; i < 3; ++i)
        if (onrightside(p, &d->points[t->vids[i]], &d->points[t->vids[(i + 1) % 3]]))
            return 0;
    return 1;
}

/** Gets index range of buckets for interval [x1, x2].
 */
static void getbrange(double x1, double x2, double xmin, double dx, int n, int* i1, int* i2)
{
    double f1 = floor((x1 - xmin) / dx);
    double f2 = floor((x2 - xmin) / dx);

    *i1 = (f1 < 0.0) ? 0 : (f1 > (double) (n - 1)) ? n - 1 : (int) f1;
    *i2 = (f2 < 0.0) ? 0 : (f2 > (double) (n - 1)) ? n - 1 : (int) f2;
}

/** Builds a bucket grid with lists of triangles whose bounding boxes overlap
 ** each bucket.
 */
static void build_buckets(triangulation* d)
{
    int nt = d->ntriangles;
    double dx = d->xmax - d->xmin;
    double dy = d->ymax - d->ymin;
    size_t nb, b, *pos;
    int t;

    if (!(dx > 0.0))
        dx = 1.0;
    if (!(dy > 0.0))
        dy = 1.0;
    d->nbx = (int) ceil(sqrt((double) (nt > 0 ? nt : 1) * dx / dy));
    if (d->nbx < 1)
        d->nbx = 1;
    else if (d->nbx > nt && nt > 0)
        d->nbx = nt;
    d->nby = (int) ceil((double) (nt > 0 ? nt : 1) / (double) d->nbx);
    if (d->nby < 1)
        d->nby = 1;
    d->bdx = dx / (double) d->nbx;
    d->bdy = dy / (double) d->nby;
    nb = (size_t) d->nbx * (size_t) d->nby;

    d->bstart = calloc(nb + 1, sizeof(size_t));
    for (t = 0; t < nt; ++t) {
        triangle* tr = &d->triangles[t];
        double x1 = DBL_MAX, x2 = -DBL_MAX, y1 = DBL_MAX, y2 = -DBL_MAX;
        int i1, i2, j1, j2, i, j, k;

        for (k = 0; k < 3; ++k) {
            point* p = &d->points[tr->vids[k]];

            if (!(p->x >= x1))
                x1 = p->x;
            if (!(p->x <= x2))
                x2 = p->x;
            if (!(p->y >= y1))
                y1 = p->y;
            if (!(p->y <= y2))
                y2 = p->y;
        }
        if (!isfinite(x1) || !isfinite(x2) || !isfinite(y1) || !isfinite(y2))
            continue;
        getbrange(x1, x2, d->xmin, d->bdx, d->nbx, &i1, &i2);
        getbrange(y1, y2, d->ymin, d->bdy, d->nby, &j1, &j2);
        for (j = j1; j <= j2; ++j)
            for (i = i1; i <= i2; ++i)
                d->bstart[(size_t) j * d->nbx + i + 1]++;
    }
    for (b = 0; b < nb; ++b)
        d->bstart[b + 1] += d->bstart[b];
    d->btris = malloc((d->bstart[nb] > 0 ? d->bstart[nb] : 1) * sizeof(int));
    pos = malloc(nb * sizeof(size_t));
    for (b = 0; b < nb; ++b)
        pos[b] = d->bstart[b];
    for (t = 0; t < nt; ++t) {
        triangle* tr = &d->triangles[t];
        double x1 = DBL_MAX, x2 = -DBL_MAX, y1 = DBL_MAX, y2 = -DBL_MAX;
        int i1, i2, j1, j2, i, j, k;

        for (k = 0; k < 3; ++k) {
            point* p = &d->points[tr->vids[k]];

            if (!(p->x >= x1))
                x1 = p->x;
            if (!(p->x <= x2))
                x2 = p->x;
            if (!(p->y >= y1))
                y1 = p->y;
            if (!(p->y <= y2))
                y2 = p->y;
        }
        if (!isfinite(x1) || !isfinite(x2) || !isfinite(y1) || !isfinite(y2))
            continue;
        getbrange(x1, x2, d->xmin, d->bdx, d->nbx, &i1, &i2);
        getbrange(y1, y2, d->ymin, d->bdy, d->nby, &j1, &j2);
        for (j = j1; j <= j2; ++j)
            for (i = i1; i <= i2; ++i)
                d->btris[pos[(size_t) j * d->nbx + i]++] = t;
    }
    free(pos);
}

/** Finds triangle containing point p using the bucket grid.
 */
static int xy2i_buckets(triangulation* d, point* p)
{
    int i, j;
    size_t b, ii;

    if (d->bstart == NULL)
        build_buckets(d);

    getbrange(p->x, p->x, d->xmin, d->bdx, d->nbx, &i, &i);
    getbrange(p->y, p->y, d->ymin, d->bdy, d->nby, &j, &j);
    b = (size_t) j* d->nbx + i;

    for (ii = d->bstart[b]; ii < d->bstart[b + 1]; ++ii) {
        int id = d->btris[ii];

        if (intriangle(d, p, &d->triangles[id]))
            return id;
    }

    return -1;
}

/* Finds triangle specified point belongs to (if any).
 *
 * First tries a straight-line walk from the last found triangle; if the walk
 * hits the boundary (e.g. islands, bays) or does not converge, falls back to
 * a search with a bucket grid.
 *
 * @param d Triangulation triangulation
 * @param p Point to be mapped
 * @return Triangle id if successful, -1 otherwhile
 */
int triangulation_xy2i(triangulation* d, point* p)
{
    int id = d->seed;
    int maxniter;
    int niter;
    triangle* t;
    int i;

    if (d->triangles == NULL || d->neighbours == NULL)
        quit("triangulation_xy2i(): incomplete structure");

    if (!(p->x >= d->xmin && p->x <= d->xmax && p->y >= d->ymin && p->y <= d->ymax))
        return -1;
    if (d->ntriangles <= 0)
        return -1;

    if (id < 0 || id >= d->ntriangles)
        id = 0;
    maxniter = 100 + 10 * (int) sqrt((double) d->ntriangles);
    t = &d->triangles[id];
    niter = 0;
    do {
        for (i = 0; i < 3; ++i) {
            int i1 = (i + 1) % 3;

            if (onrightside(p, &d->points[t->vids[i]], &d->points[t->vids[i1]])) {
                id = d->neighbours[id].tids[(i + 2) % 3];
                break;
            }
        }
        if (i < 3) {
            if (id < 0 || ++niter > maxniter) {
                id = xy2i_buckets(d, p);
                break;
            }
            t = &d->triangles[id];
        }
    } while (i < 3);

    if (id >= 0)
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
