/******************************************************************************
 *
 * File:        gxy_curv.c
 *
 * Created:     27/11/2019
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Purpose:     Code for 'gxy_curv' object (handling of structured curvilinear
 *              grids).
 *
 * Description: The code has been extracted from and replaces that from
 *              gridutils library.
 *
 * Revisions:   
 *
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "kdtree.h"
#include "definitions.h"
#include "utils.h"
#include "grid.h"
#include "hgrid.h"
#include "gxy_curv.h"

#define EPS 1.0e-8
#define EPS_ZERO 1.0e-5

struct gxy_curv {
    hgrid* parent;
    int ni;
    int nj;
    double** x;
    double** y;
    int sign;
    kdtree* nodetreeXY;
#if defined(USE_SHMEM)
    MPI_Win sm_comm_win;
#endif
};

/**
 */
gxy_curv* gxy_curv_create(hgrid* hg, int ni, int nj, double** x, double** y, int** mask)
{
    gxy_curv* gxy = malloc(sizeof(gxy_curv));
    char kdname[MAXSTRLEN];
    double* nodecoords[3];

    assert(x != NULL && y != NULL);

    hg->parent = hg;
    gxy->ni = ni;
    gxy->nj = nj;
    gxy->x = x;
    gxy->y = y;
    gxy->sign = 0;

    snprintf(kdname, MAXSTRLEN - 4, "%s_XY", grid_getname(hg->parent));
    if (!hg->geographic) {
        nodecoords[0] = x[0];
        nodecoords[1] = y[0];
        gxy->nodetreeXY = kd_create(kdname, 2);
    } else {
        int i, c;

        for (c = 0; c < 3; ++c)
            nodecoords[c] = malloc(ni * nj * sizeof(double));

        for (i = 0; i < ni * nj; ++i) {
            double ll[2] = { x[0][i], y[0][i] };
            double xyz[3];

            ll2xyz(ll, xyz);
            for (c = 0; c < 3; ++c)
                nodecoords[c][i] = xyz[c];
        }
        gxy->nodetreeXY = kd_create(kdname, 3);
    }
#if defined(USE_SHMEM)
    {
        MPI_Aint size;
        int ierror;
        size_t nnodes, i;

        if (mask == NULL)
            nnodes = ni * nj;
        else
            for (i = 0, nnodes = 0; i < ni * nj; ++i)
                if (mask[i] != 0)
                    nnodes++;

        size = kd_getstoragesize(gxy->nodetreeXY, nnodes);
        if (sm_comm_rank == 0) {
            void* storage = NULL;

            assert(sizeof(MPI_Aint) == sizeof(size_t));
            ierror = MPI_Win_allocate_shared(size, sizeof(double), MPI_INFO_NULL, sm_comm, &storage, &gxy->sm_comm_win);
            assert(ierror == MPI_SUCCESS);
            kd_setstorage(gxy->nodetreeXY, nnodes, storage, 1);
            kd_insertnodes(gxy->nodetreeXY, ni * nj, nodecoords, NULL, (mask != NULL) ? mask[0] : NULL, 1);
        } else {
            MPI_Aint my_size;
            void* storage = NULL;
            int disp_unit;

            ierror = MPI_Win_allocate_shared(0, sizeof(double), MPI_INFO_NULL, sm_comm, &storage, &gxy->sm_comm_win);
            assert(ierror == MPI_SUCCESS);
            ierror = MPI_Win_shared_query(gxy->sm_comm_win, 0, &my_size, &disp_unit, &storage);
            assert(ierror == MPI_SUCCESS);
            assert(my_size = size);
            kd_setstorage(gxy->nodetreeXY, nnodes, storage, 0);
            kd_syncsize(gxy->nodetreeXY);
        }
        MPI_Win_fence(0, gxy->sm_comm_win);
        MPI_Barrier(sm_comm);
    }
#else
    kd_insertnodes(gxy->nodetreeXY, ni * nj, nodecoords, NULL, (mask != NULL) ? mask[0] : NULL, 1);
#endif

    if (hg->geographic) {
        free(nodecoords[0]);
        free(nodecoords[1]);
        free(nodecoords[2]);
    }

    return gxy;
}

int gxy_curv_destroykdtree(gxy_curv* gxy)
{
    if (gxy->nodetreeXY == NULL)
        return 0;
    kd_destroy(gxy->nodetreeXY);
    gxy->nodetreeXY = NULL;
#if defined(USE_SHMEM)
    MPI_Win_free(&gxy->sm_comm_win);
    assert(gxy->sm_comm_win == MPI_WIN_NULL);
#endif

    return 1;
}

/**
 */
void gxy_curv_destroy(gxy_curv* gxy)
{
    if (gxy == NULL)
        return;

    free(gxy->x);
    free(gxy->y);
    gxy_curv_destroykdtree(gxy);
    free(gxy);
}

/**
 */
double** gxy_curv_getx(gxy_curv* gxy)
{
    return gxy->x;
}

/**
 */
double** gxy_curv_gety(gxy_curv* gxy)
{
    return gxy->y;
}

/**
 */
int gxy_curv_getni(gxy_curv* gxy)
{
    return gxy->ni;
}

/**
 */
int gxy_curv_getnj(gxy_curv* gxy)
{
    return gxy->nj;
}

/** Test whether the point (x, y) is inside quadrilateral {(px[0], py[0]),
 ** ..., (px[3], py[3])}. Based on Algorithm 1 from Symmetry 2018, 10, 477;
 ** doi:10.3390/sym10100477.
 */
static int incell(double x, double y, double* px, double* py)
{
    int count, i, i1;
    double v1, v2, u1, u2, f;

    for (i = 0, count = 0; i < 4; ++i) {
        i1 = (i + 1) % 4;
        v1 = py[i] - y;
        v2 = py[i1] - y;
        if ((v1 < 0.0 && v2 < 0.0) || (v1 > 0.0 && v2 > 0.0))
            continue;
        u1 = px[i] - x;
        u2 = px[i1] - x;
        f = u1 * v2 - u2 * v1;
        if (v2 > 0.0 && v1 <= 0.0) {
            if (f > 0.0)
                count++;
            else if (f == 0.0)
                return 1;
        } else if (v1 > 0.0 && v2 <= 0.0) {
            if (f < 0.0)
                count++;
            else if (f == 0.0)
                return 1;
        } else if (v2 == 0.0 && v1 < 0.0) {
            if (f == 0.0)
                return 1;
        } else if (v1 == 0.0 && v2 < 0.0) {
            if (f == 0.0)
                return 1;
        } else if (v1 == 0.0 && v2 == 0.0) {
            if (u2 <= 0.0 && u1 >= 0.0)
                return 1;
            else if (u1 <= 0.0 && u2 >= 0.0)
                return 1;
        }
    }

    return count % 2;
}

/**
 */
static int gxy_curv_xy2ij(gxy_curv* gxy, double x, double y, int* iout, int* jout)
{
    double* minmax;
    double pos[3];
    size_t nearest;
    size_t id;
    int i, j, i1, i2, j1, j2;
    double px[4], py[4];

    minmax = kd_getminmax(gxy->nodetreeXY);
    if (!gxy->parent->geographic) {
        pos[0] = x;
        pos[1] = y;
        if (x < minmax[0] || y < minmax[1] || x > minmax[2] || y > minmax[3])
            return 0;
    } else {
        double ll[2] = { x, y };

        ll2xyz(ll, pos);
        if (pos[0] < minmax[0] || pos[1] < minmax[1] || pos[0] > minmax[3] || pos[1] > minmax[4])
            return 0;
    }

    /*
     * this is a rather expensive call, O(log N)
     */
    nearest = kd_findnearestnode(gxy->nodetreeXY, pos);
    id = kd_getnodedata(gxy->nodetreeXY, nearest);

    j = id / (gxy->ni);
    i = id % (gxy->ni);

    /*
     * the code below could be somewhat optimised, I guess, but it would not
     * matter much because of the rather expensive search above
     */
    i1 = (i > 0) ? i - 1 : i;
    i2 = (i < gxy->ni - 1) ? i + 1 : i;
    j1 = (j > 0) ? j - 1 : j;
    j2 = (j < gxy->nj - 1) ? j + 1 : j;

    for (j = j1; j <= j2 - 1; ++j) {
        for (i = i1; i <= i2 - 1; ++i) {
            if (!isfinite(gxy->x[j][i]) || !isfinite(gxy->x[j][i + 1]) || !isfinite(gxy->x[j + 1][i + 1]) || !isfinite(gxy->x[j + 1][i]))
                continue;

            px[0] = gxy->x[j][i];
            py[0] = gxy->y[j][i];
            px[1] = gxy->x[j][i + 1];
            py[1] = gxy->y[j][i + 1];
            px[2] = gxy->x[j + 1][i + 1];
            py[2] = gxy->y[j + 1][i + 1];
            px[3] = gxy->x[j + 1][i];
            py[3] = gxy->y[j + 1][i];

            if (incell(x, y, px, py)) {
                *iout = i;
                *jout = j;
                return 1;
            }
        }
    }

    return 0;
}

/**
 */
static int calc_branch(gxy_curv* gxy, double x, double y)
{
    double** gx = gxy->x;
    double** gy = gxy->y;
    int sign = 1;
    int i, j;
    double error[2];

    /*
     * try xy2ij() first 
     */
    if (gxy_curv_xy2ij(gxy, x, y, &i, &j) == 0)
        return 0;               /* failed */

    {
        double a = gx[j][i] - gx[j][i + 1] - gx[j + 1][i] + gx[j + 1][i + 1];
        double b = gx[j][i + 1] - gx[j][i];
        double c = gx[j + 1][i] - gx[j][i];
        double d = gx[j][i];
        double e = gy[j][i] - gy[j][i + 1] - gy[j + 1][i] + gy[j + 1][i + 1];
        double f = gy[j][i + 1] - gy[j][i];
        double g = gy[j + 1][i] - gy[j][i];
        double h = gy[j][i];

        double A = a * f - b * e;

        double B, C;
        int k;

        /*
         * normally one checks A before calling calc_branch() 
         */
        if (fabs(A) < EPS_ZERO)
            return 0;           /* failed */

        B = e * x - a * y + a * h - d * e + c * f - b * g;
        C = g * x - c * y + c * h - d * g;

        for (k = 0; k < 2; ++k) {
            double u = (-B + sign * sqrt(B * B - 4.0 * A * C)) / (2.0 * A);
            double v_denom = a * u + c;
            double v = (fabs(v_denom) < EPS_ZERO) ? (y - f * u - h) / (e * u + g) : (x - b * u - d) / v_denom;

            error[k] = 0.0;

            if (u < 0.0)
                error[k] -= u;
            else if (u > 1.0)
                error[k] += (u - 1.0);
            if (v < 0.0)
                error[k] -= v;
            else if (v > 1.0)
                error[k] += (v - 1.0);

            sign = -1;
        }
    }

    if (error[0] < error[1])
        return 1;
    return -1;
}

/** Maps physical coordinates to fractional grid indices.
 */
int gxy_curv_xy2fij(gxy_curv* gxy, double x, double y, double* fij)
{
    double** gx = gxy->x;
    double** gy = gxy->y;
    int i, j;

    fij[0] = NAN;
    fij[1] = NAN;

    if (gxy_curv_xy2ij(gxy, x, y, &i, &j) == 0)
        return 0;               /* failed */

    {
        double a = gx[j][i] - gx[j][i + 1] - gx[j + 1][i] + gx[j + 1][i + 1];
        double b = gx[j][i + 1] - gx[j][i];
        double c = gx[j + 1][i] - gx[j][i];
        double d = gx[j][i];
        double e = gy[j][i] - gy[j][i + 1] - gy[j + 1][i] + gy[j + 1][i + 1];
        double f = gy[j][i + 1] - gy[j][i];
        double g = gy[j + 1][i] - gy[j][i];
        double h = gy[j][i];

        double A = a * f - b * e;
        double B = e * x - a * y + a * h - d * e + c * f - b * g;
        double C = g * x - c * y + c * h - d * g;

        double u, v, d1, d2;

        if (fabs(A) < EPS_ZERO)
            u = -C / B * (1.0 + A * C / B / B);
        else {
            if (gxy->sign == 0) {
                gxy->sign = calc_branch(gxy, x, y);
                if (gxy->sign == 0)
                    return 0;   /* failed */
            }
            u = (-B + gxy->sign * sqrt(B * B - 4.0 * A * C)) / (2.0 * A);
        }
        d1 = a * u + c;
        d2 = e * u + g;
        v = (fabs(d2) > fabs(d1)) ? (y - f * u - h) / d2 : (x - b * u - d) / d1;

        /*
         * this can happen e.g. when there is a singularity (longitude
         * discontinuity) in this cell
         */
        if (hypot(u - 0.5, v - 0.5) > 1)
            return 0;

        if (u < 0.0)
            u = 0.0;
        else if (u >= 1.0)
            u = 1.0 - EPS;
        if (v < 0.0)
            v = 0.0;
        else if (v >= 1.0)
            v = 1.0 - EPS;

        if ((int) ((double) i + u) > i)
            u = u - EPS;
        if ((int) ((double) j + v) > j)
            v = v - EPS;
        assert((int) ((double) i + u) == i);
        assert((int) ((double) j + v) == j);

        fij[0] = i + u;
        fij[1] = j + v;
    }

    return 1;
}

/** Maps fractional grid indices to physical coordinates.
 */
int gxy_curv_fij2xy(gxy_curv* gxy, double fi, double fj, double* x, double* y)
{
    double** gx = gxy->x;
    double** gy = gxy->y;
    int i, j;
    double u, v;
    double a, b, c, d, e, f, g, h;

    if (fi < 0.0 || fi > (double) gxy->ni - 1.0 || fj < 0.0 || fj > (double) gxy->nj - 1.0)
        return 0;

    u = fi - floor(fi);
    v = fj - floor(fj);
    i = (int) fi;
    j = (int) fj;
    if (u == 0.0 && v == 0.0) {
        *x = gx[j][i];
        *y = gy[j][i];
    } else if (u == 0.0) {
        *x = gx[j + 1][i] * v + gx[j][i] * (1.0 - v);
        *y = gy[j + 1][i] * v + gy[j][i] * (1.0 - v);
    } else if (v == 0.0) {
        *x = gx[j][i + 1] * u + gx[j][i] * (1.0 - u);
        *y = gy[j][i + 1] * u + gy[j][i] * (1.0 - u);
    } else {
        a = gx[j][i] - gx[j][i + 1] - gx[j + 1][i] + gx[j + 1][i + 1];
        b = gx[j][i + 1] - gx[j][i];
        c = gx[j + 1][i] - gx[j][i];
        d = gx[j][i];
        e = gy[j][i] - gy[j][i + 1] - gy[j + 1][i] + gy[j + 1][i + 1];
        f = gy[j][i + 1] - gy[j][i];
        g = gy[j + 1][i] - gy[j][i];
        h = gy[j][i];

        *x = a * u * v + b * u + c * v + d;
        *y = e * u * v + f * u + g * v + h;
    }

    return 1;
}

/**
 */
kdtree* gxy_curv_gettree(gxy_curv* gxy)
{
    return gxy->nodetreeXY;
}
