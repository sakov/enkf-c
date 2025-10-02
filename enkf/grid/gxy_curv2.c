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
 *              projections. This code is supposed to provide robust handling
 *              for grids with  singularities in lon/lat (such as ORCA grids).
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
#include "hgrid.h"
#include "gxy_curv2.h"

#define EPS 1.0e-8
#define EPS_ZERO 1.0e-5
#define USE_LL2XYZ 0

/*
 * "gxy_curv2" projects grids to the "North" and "South" stereographic
 * projections (SGPs). (By North SGP we call the projection from the North Pole
 * onto the plane tangent at the South Pole.)
 *
 * A point is mapped using the NSGP for points in the southern hemisphere, and
 * SSGP for points in the northern hemisphere. This design with two projections
 * targets global grids, and is likely to be somewhat excessive for local grids.
 *
 * Compared to "gxy_curv", "gxy_curv2" uses quite a bit more memory. In
 * addition to the original lon/lat grid coordinates (necessary for ij2xy
 * mappings) it also needs to carry two Kd-trees and two arrays of grid
 * cordinates. Speed-wise it should be as good as the original "gxy_curv".
 */

struct gxy_curv2 {
    hgrid* parent;
    int ni;
    int nj;
    double** x_orig;
    double** y_orig;
    /*
     * every field below is duplicated for two projections:
     * [0] is related to the North SGP;
     * [1] - to the South SGP
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
gxy_curv2* gxy_curv2_create(hgrid* hg, int ni, int nj, double** x, double** y, int** mask)
{
    char kdname[MAXSTRLEN];
    gxy_curv2* gxy = calloc(sizeof(gxy_curv2), 1);
    int proj;

    gxy->parent = hg;
    gxy->ni = ni;
    gxy->nj = nj;
    gxy->x_orig = x;
    gxy->y_orig = y;

    snprintf(kdname, MAXSTRLEN - 4, "%s_XY2", grid_getname(hg->parent));
    for (proj = 0; proj < 2; ++proj) {
        double* nodecoords[2];
        int sign, i;

        gxy->x[proj] = alloc2d(nj, ni, sizeof(double));
        gxy->y[proj] = alloc2d(nj, ni, sizeof(double));

        nodecoords[0] = gxy->x[proj][0];
        nodecoords[1] = gxy->y[proj][0];

        sign = (proj == 0) ? 1 : -1;
#if USE_LL2XYZ
        for (i = 0; i < ni * nj; ++i) {
            double ll[2] = { x[0][i], y[0][i] * sign };
            double xyz[3];

            ll2xyz(ll, xyz);
            nodecoords[0][i] = xyz[0] / (REARTH - xyz[2]);
            nodecoords[1][i] = xyz[1] / (REARTH - xyz[2]);
        }
#else
        for (i = 0; i < ni * nj; ++i) {
            double r = 1.0 / tan(M_PI_4 - sign * DEG2RAD * y[0][i] / 2.0);
            double phi = DEG2RAD * x[0][i];

            nodecoords[0][i] = r * cos(phi);
            nodecoords[1][i] = r * sin(phi);
        }
#endif
        gxy->nodetreeXY[proj] = kd_create(kdname, 2);
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
void gxy_curv2_destroykdtree(gxy_curv2* gxy)
{
    int proj;

    if (gxy->x[0] == NULL)
        return;

    for (proj = 0; proj < 2; ++proj) {
        /*
         * calling this function means that no more mappings are intended,
         * therefore free grid arrays on stereographic projections
         */
        free(gxy->x[proj]);
        gxy->x[proj] = NULL;
        free(gxy->y[proj]);
        gxy->y[proj] = NULL;

        kd_destroy(gxy->nodetreeXY[proj]);
        gxy->nodetreeXY[proj] = NULL;
#if defined(USE_SHMEM)
        MPI_Win_free(&gxy->sm_comm_win[proj]);
        assert(gxy->sm_comm_win[proj] == MPI_WIN_NULL);
#endif
    }
}

/**
 */
void gxy_curv2_destroy(gxy_curv2* gxy)
{
    if (gxy == NULL)
        return;

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
    int i, j;

    fij[0] = NAN;
    fij[1] = NAN;

    proj = (y < 0.0) ? 0 : 1;

#if USE_LL2XYZ
    {
        double ll[2] = { x, -fabs(y) };
        double pos[3];

        ll2xyz(ll, pos);

        x = pos[0] / (REARTH - pos[2]);
        y = pos[1] / (REARTH - pos[2]);
    }
#else
    {
        double r = 1.0 / tan(M_PI_4 + DEG2RAD * fabs(y) / 2.0);
        double phi = DEG2RAD * x;

        x = r * cos(phi);
        y = r * sin(phi);
    }
#endif

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
                if (gxy->sign[proj] == 0)
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
