/******************************************************************************
 *
 * File:        gxy_curv2.c
 *
 * Created:     06/05/2024
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Purpose:     Handles geographic curvilinear grids using stereographic
 *              projections. This code is supposed to be indifferent to
 *              grid singularities in lon/lat (such as discontinuity in node
 *              latitudes on ORCA grids).
 *
 * Description: Derived from gxy_curv.c.
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
#include "gxy_curv2.h"

#define EPS 1.0e-8
#define EPS_ZERO 1.0e-5

struct gxy_curv2 {
    char* name;
    int ni;
    int nj;
    double** x_orig;
    double** y_orig;
    /*
     * every field below is duplicated for two projections:
     * [0] is related to sterographic projection from the North Pole;
     * [1] - from the South Pole
     */
    double** x[2];
    double** y[2];
    int sign[2];
    kdtree* nodetreeXY[2];
#if defined(USE_SHMEM)
    MPI_Win sm_comm_win[2];
#endif
};

/**
 */
gxy_curv2* gxy_curv2_create(void* g, int ni, int nj, double** x, double** y, int** mask)
{
    char name[MAXSTRLEN];
    gxy_curv2* gxy = calloc(sizeof(gxy_curv2), 1);
    int proj;

    snprintf(name, MAXSTRLEN - 4, "%s_XY2", grid_getname(g));

    gxy->name = strdup(name);
    gxy->ni = ni;
    gxy->nj = nj;
    gxy->x_orig = x;
    gxy->y_orig = y;

    for (proj = 0; proj < 2; ++proj) {
        double* nodecoords[2];
        int i;

        gxy->x[proj] = alloc2d(nj, ni, sizeof(double));
        gxy->y[proj] = alloc2d(nj, ni, sizeof(double));

        nodecoords[0] = gxy->x[proj][0];
        nodecoords[1] = gxy->y[proj][0];

        for (i = 0; i < ni * nj; ++i) {
            double ll[2] = { x[0][i], (proj == 0) ? y[0][i] : -y[0][i] };
            double xyz[3];

            ll2xyz(ll, xyz);
            nodecoords[0][i] = xyz[0] / (REARTH - xyz[2]);
            nodecoords[1][i] = xyz[1] / (REARTH - xyz[2]);
        }
        gxy->nodetreeXY[proj] = kd_create(name, 2);
#if defined(USE_SHMEM)
        {
            MPI_Aint size;
            int ierror;
            size_t nnodes;

            if (mask == NULL)
                nnodes = ni * nj;
            else
                for (i = 0, nnodes = 0; i < ni * nj; ++i)
                    if (mask[i] != 0)
                        nnodes++;

            size = kd_getstoragesize(gxy->nodetreeXY[proj], nnodes);

            if (sm_comm_rank == 0) {
                void* storage = NULL;

                assert(sizeof(MPI_Aint) == sizeof(size_t));
                ierror = MPI_Win_allocate_shared(size, sizeof(double), MPI_INFO_NULL, sm_comm, &storage, &gxy->sm_comm_win[proj]);
                assert(ierror == MPI_SUCCESS);
                kd_setstorage(gxy->nodetreeXY[proj], nnodes, storage, 1);
                kd_insertnodes(gxy->nodetreeXY[proj], ni * nj, nodecoords, NULL, (mask != NULL) ? mask[0] : NULL, 1);
            } else {
                MPI_Aint my_size;
                void* storage = NULL;
                int disp_unit;

                ierror = MPI_Win_allocate_shared(0, sizeof(double), MPI_INFO_NULL, sm_comm, &storage, &gxy->sm_comm_win[proj]);
                assert(ierror == MPI_SUCCESS);
                ierror = MPI_Win_shared_query(gxy->sm_comm_win[proj], 0, &my_size, &disp_unit, &storage);
                assert(ierror == MPI_SUCCESS);
                assert(my_size = size);
                kd_setstorage(gxy->nodetreeXY[proj], nnodes, storage, 0);
                kd_syncsize(gxy->nodetreeXY[proj]);
            }
            MPI_Win_fence(0, gxy->sm_comm_win[proj]);
            MPI_Barrier(sm_comm);
        }
#else
        kd_insertnodes(gxy->nodetreeXY[proj], ni * nj, nodecoords, NULL, (mask != NULL) ? mask[0] : NULL, 1);
#endif
    }

    return gxy;
}

/**
 */
int gxy_curv2_destroykdtree(gxy_curv2* gxy)
{
    int proj;

    for (proj = 0; proj < 2; ++proj) {
        /*
         * calling this function means that no more mappings are intended,
         * therefore free grid arrays on stereographic projections
         */
        free(gxy->x[proj]);
        gxy->x[proj] = NULL;
        free(gxy->y[proj]);
        gxy->y[proj] = NULL;

        if (gxy->nodetreeXY[proj] == NULL)
            continue;
        kd_destroy(gxy->nodetreeXY[proj]);
        gxy->nodetreeXY[proj] = NULL;
#if defined(USE_SHMEM)
        MPI_Win_free(&gxy->sm_comm_win[proj]);
        assert(gxy->sm_comm_win[proj] == MPI_WIN_NULL);
#endif
    }

    return 1;
}

/**
 */
void gxy_curv2_destroy(gxy_curv2* gxy)
{
    if (gxy == NULL)
        return;

    free(gxy->name);
    free(gxy->x_orig);
    free(gxy->y_orig);
    gxy_curv2_destroykdtree(gxy);
    free(gxy);
}

/**
 */
double** gxy_curv2_getx(gxy_curv2* gxy)
{
    return gxy->x_orig;
}

/**
 */
double** gxy_curv2_gety(gxy_curv2* gxy)
{
    return gxy->y_orig;
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

/** This function is called only from gxy_curv2_xy2fij(); x and y are assumed
 ** to be already in sterographic projection
 */
static int gxy_curv2_xy2ij(gxy_curv2* gxy, int proj, double x, double y, int* iout, int* jout)
{
    double** xx = gxy->x[proj];
    double** yy = gxy->y[proj];
    kdtree* tree = gxy->nodetreeXY[proj];

    int i, j;
    int i1, i2, j1, j2;
    double px[4], py[4];

    {
        double* minmax = kd_getminmax(tree);

        if (x < minmax[0] || y < minmax[1] || x > minmax[2] || y > minmax[3])
            return 0;
    }

    /*
     * this is a rather expensive call, O(log N)
     */
    {
        double pos[2] = { x, y };
        size_t nearest;
        size_t id;

        nearest = kd_findnearestnode(tree, pos);
        id = kd_getnodedata(tree, nearest);

        j = id / (gxy->ni);
        i = id % (gxy->ni);
    }

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
            if (!isfinite(xx[j][i]) || !isfinite(xx[j][i + 1]) || !isfinite(xx[j + 1][i + 1]) || !isfinite(xx[j + 1][i]))
                continue;

            px[0] = xx[j][i];
            py[0] = yy[j][i];
            px[1] = xx[j][i + 1];
            py[1] = yy[j][i + 1];
            px[2] = xx[j + 1][i + 1];
            py[2] = yy[j + 1][i + 1];
            px[3] = xx[j + 1][i];
            py[3] = yy[j + 1][i];

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
static int calc_branch(gxy_curv2* gxy, int proj, double x, double y)
{
    double** gx = gxy->x[proj];
    double** gy = gxy->y[proj];
    int sign = 1;
    int i, j;
    double error[2];

    /*
     * try xy2ij() first 
     */
    if (gxy_curv2_xy2ij(gxy, proj, x, y, &i, &j) == 0)
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
int gxy_curv2_xy2fij(gxy_curv2* gxy, double x, double y, double* fij)
{
    double** gx;
    double** gy;
    int proj;
    double ll[2], pos[3];
    int i, j;

    fij[0] = NAN;
    fij[1] = NAN;

    proj = (y < 0.0) ? 0 : 1;
    ll[0] = x;
    ll[1] = (y < 0.0) ? y : -y;

    ll2xyz(ll, pos);
    x = pos[0] / (REARTH - pos[2]);
    y = pos[1] / (REARTH - pos[2]);

    if (gxy_curv2_xy2ij(gxy, proj, x, y, &i, &j) == 0)
        return 0;               /* failed */

    gx = gxy->x[proj];
    gy = gxy->y[proj];
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
            if (gxy->sign[proj] == 0) {
                gxy->sign[proj] = calc_branch(gxy, proj, x, y);
                if (gxy->sign == 0)
                    return 0;   /* failed */
            }
            u = (-B + gxy->sign[proj] * sqrt(B * B - 4.0 * A * C)) / (2.0 * A);
        }
        d1 = a * u + c;
        d2 = e * u + g;
        v = (fabs(d2) > fabs(d1)) ? (y - f * u - h) / d2 : (x - b * u - d) / d1;

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
