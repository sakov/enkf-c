/******************************************************************************
 *
 * File:        gxy_rect.h        
 *
 * Created:     25/01/2022
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Purpose:     Code for `gxy_rect' object (handling of "rectangular"
 *              horizontal grids).
 *
 * Description: By "rectangular" grid we understand a structured quadrilateral
 *              grid with node rows and columns aligned with coordinates.
 *
 * Revisions:
 *
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "utils.h"
#include "hgrid.h"
#include "gxy_rect.h"

#define EPS_LON 1.0e-3

/*
 * Note: physical coordinates are mapped to fractional indices by linear
 * interpolation between nodes (this is equivalent to using cell edges at
 * midpoints between nodes); therefore cell edges are not stored.
 */
struct gxy_rect {
    hgrid* parent;
    int ni;
    int nj;
    int regular_x;
    int regular_y;

    double* x;
    double* y;
};

/** Checks that coordinates are strictly monotonic.
 */
static void check_monotonic(int n, double* v, char* name)
{
    int ascending = (v[n - 1] > v[0]);
    int i;

    for (i = 0; i < n - 1; ++i)
        if (!(ascending ? (v[i + 1] > v[i]) : (v[i + 1] < v[i])))
            enkf_quit("non-monotonic %s coordinate for a rect grid (%d: %f, %d: %f)", name, i, v[i], i + 1, v[i + 1]);
}

/**
 */
gxy_rect* gxy_rect_create(hgrid* hg, int ni, int nj, double* x, double* y)
{
    gxy_rect* gxy = malloc(sizeof(gxy_rect));
    int i;
    double dx, dy;

    if (ni < 2 || nj < 2)
        enkf_quit("rect grid must have at least 2 nodes in each direction (ni = %d, nj = %d)", ni, nj);

    gxy->parent = hg;

    /*
     * x
     */
    gxy->ni = ni;
    gxy->x = x;
    check_monotonic(ni, x, "X");

    if (hg->geographic) {
        if (fabs(fabs(gxy->x[ni - 1] - gxy->x[0]) - 360.0) < EPS_LON)
            hg->periodic_i = 1; /* closed grid */
        else if (fabs(fabs(2.0 * gxy->x[ni - 1] - gxy->x[ni - 2] - gxy->x[0]) - 360.0) < EPS_LON) {
            hg->periodic_i = 2; /* non-closed grid (used e.g. by MOM) */
            gxy->x = realloc(gxy->x, (ni + 1) * sizeof(double));
            gxy->x[ni] = 2.0 * gxy->x[ni - 1] - gxy->x[ni - 2];
        } else
            hg->periodic_i = 0;
    } else
        hg->periodic_i = 0;

    dx = (gxy->x[ni - 1] - gxy->x[0]) / (double) (ni - 1);
    for (i = 1; i < (int) ni; ++i)
        if (fabs(gxy->x[i] - gxy->x[i - 1] - dx) / fabs(dx) > EPS_LON)
            break;
    gxy->regular_x = (i == ni);

    /*
     * y
     */
    gxy->nj = nj;
    gxy->y = y;
    check_monotonic(nj, y, "Y");

    dy = (y[nj - 1] - y[0]) / (double) (nj - 1);
    for (i = 1; i < (int) nj; ++i)
        if (fabs(y[i] - y[i - 1] - dy) / fabs(dy) > EPS_LON)
            break;
    gxy->regular_y = (i == nj);

    return gxy;
}

/**
 */
double* gxy_rect_getx(gxy_rect* gxy)
{
    return gxy->x;
}

/**
 */
double* gxy_rect_gety(gxy_rect* gxy)
{
    return gxy->y;
}

/**
 */
int gxy_rect_getni(gxy_rect* gxy)
{
    return gxy->ni;
}

/**
 */
int gxy_rect_getnj(gxy_rect* gxy)
{
    return gxy->nj;
}

/**
 */
static double fi2x(int n, double* v, double fi, int periodic)
{
    double ifrac;
    int i;

    if (n < 2)
        return NAN;

    if (periodic == 2)
        n += 1;

    if (periodic) {
        if (fi < 0.0)
            fi += (double) (n - 1);
        else if (fi >= (double) (n - 1))
            fi -= (double) (n - 1);
    } else {
        if (fi < 0.0 || fi > (double) (n - 1))
            return NAN;
    }

    ifrac = fi - floor(fi);
    i = (int) fi;
    if (ifrac == 0.0)           /* (also covers i = n - 1) */
        return v[i];
    else
        return v[i] + ifrac * (v[i + 1] - v[i]);
}

/** Convert fractional indices to physical coordinates.
 */
void gxy_rect_fij2xy(gxy_rect* gxy, double fi, double fj, double* x, double* y)
{
    *x = fi2x(gxy->ni, gxy->x, fi, gxy->parent->periodic_i);
    *y = fi2x(gxy->nj, gxy->y, fj, 0);
}

/** Find fractional index for a coordinate in a regularly spaced 1D grid.
 */
#if defined(ENKF_PREP) || defined(ENKF_CALC)
static double x2fi_reg(int n, double* v, double x, int periodic)
{
    double fi;

    if (n < 2)
        return NAN;

    if (periodic == 2)
        n += 1;

    fi = (x - v[0]) / (v[n - 1] - v[0]) * (double) (n - 1);

    if (periodic) {
        fi = fmod(fi, (double) (n - 1));
        if (fi < 0.0)
            fi += (double) (n - 1);
        if (fi >= (double) (n - 1))     /* (rounding) */
            fi = 0.0;
    } else if (fi < 0.0 || fi > (double) (n - 1))
        return NAN;

    return fi;
}
#endif

/** Gets fractional index of a coordinate for a 1D irregular grid, by linear
 ** interpolation between nodes.
 * @param n Number of grid nodes (vertical layers)
 * @param v Coordinates of the nodes (layer centres)
 * @param x Input coordinate
 * @param periodic Flag for grid periodicity
 * @param lonbase lon range = [lonbase, lonbase + 360)
 * @return Fractional index for `x'
 */
#if defined(ENKF_PREP) || defined(ENKF_CALC)
static double x2fi_irreg(int n, double v[], double x, int periodic, double lonbase)
{
    int ascending, i1, i2, imid;
    double fi;

    if (n < 2)
        return NAN;

    if (periodic == 2)
        n += 1;

    if (!isnan(lonbase)) {
        x = lonbase + fmod(x - lonbase, 360.0);
        if (x < lonbase)
            x += 360.0;
    }

    ascending = (v[n - 1] > v[0]);

    if (periodic) {
        double vmin = (ascending) ? v[0] : v[n - 1];
        double vmax = (ascending) ? v[n - 1] : v[0];

        x = vmin + fmod(x - vmin, 360.0);
        if (x < vmin)
            x += 360.0;
        if (x > vmax)           /* (rounding gap at the seam) */
            x = vmax;
    }

    if (ascending) {
        if (x < v[0] || x > v[n - 1])
            return NAN;
    } else {
        if (x > v[0] || x < v[n - 1])
            return NAN;
    }

    i1 = 0;
    i2 = n - 1;
    while (i2 - i1 > 1) {
        imid = (i1 + i2) / 2;
        if ((ascending) ? (x > v[imid]) : (x < v[imid]))
            i1 = imid;
        else
            i2 = imid;
    }
    fi = (double) i1 + (x - v[i1]) / (v[i1 + 1] - v[i1]);
    if (periodic && fi >= (double) (n - 1))
        fi -= (double) (n - 1);

    return fi;
}
#endif

/** Convert physical coordinates to fractional indices.
 */
#if defined(ENKF_PREP) || defined(ENKF_CALC)
void gxy_rect_xy2fij(gxy_rect* gxy, double x, double y, double* fij)
{
    if (gxy->regular_x)
        fij[0] = x2fi_reg(gxy->ni, gxy->x, x, gxy->parent->periodic_i);
    else
        fij[0] = x2fi_irreg(gxy->ni, gxy->x, x, gxy->parent->periodic_i, gxy->parent->lonbase);
    if (gxy->regular_y)
        fij[1] = x2fi_reg(gxy->nj, gxy->y, y, 0);
    else
        fij[1] = x2fi_irreg(gxy->nj, gxy->y, y, 0, NAN);
}
#endif

/**
 */
void gxy_rect_destroy(gxy_rect* gxy)
{
    free(gxy->x);
    free(gxy->y);
    free(gxy);
}
