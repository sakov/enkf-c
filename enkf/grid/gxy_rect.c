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

struct gxy_rect {
    hgrid* parent;
    int ni;
    int nj;
    int regular_x;
    int regular_y;

    double* x;
    double* y;
    double* xc;
    double* yc;
};

/**
 */
gxy_rect* gxy_rect_create(hgrid* hg, int ni, int nj, double* x, double* y)
{
    gxy_rect* gxy = malloc(sizeof(gxy_rect));
    int i, ascending;
    double dx, dy;

    assert(ni >= 2 && nj >= 2);

    gxy->parent = hg;

    /*
     * x
     */
    gxy->ni = ni;
    gxy->x = x;
    gxy->xc = malloc((ni + 1) * sizeof(double));
    gxy->xc[0] = x[0] * 1.5 - x[1] * 0.5;
    for (i = 1; i < ni + 1; ++i)
        gxy->xc[i] = 2 * x[i - 1] - gxy->xc[i - 1];

    if (hg->geographic) {
        if (fabs(fmod(gxy->x[ni - 1] - gxy->x[0] + EPS_LON / 2.0, 360.0)) < EPS_LON)
            hg->periodic_i = 1; /* closed grid */
        else if (fabs(fmod(2.0 * gxy->x[ni - 1] - gxy->x[ni - 2] - gxy->x[0] + EPS_LON / 2.0, 360.0)) < EPS_LON) {
            hg->periodic_i = 2; /* non-closed grid (used e.g. by MOM) */
            gxy->x = realloc(gxy->x, (ni + 1) * sizeof(double));
            gxy->x[ni] = 2.0 * gxy->x[ni - 1] - gxy->x[ni - 2];
            gxy->xc = realloc(gxy->xc, (ni + 2) * sizeof(double));
            gxy->xc[ni + 1] = 2.0 * gxy->xc[ni] - gxy->xc[ni - 1];
        } else
            hg->periodic_i = 0;
    } else
        hg->periodic_i = 0;

    dx = (gxy->x[ni - 1] - gxy->x[0]) / (double) (ni - 1);
    for (i = 1; i < (int) ni; ++i)
        if (fabs(gxy->x[i] - gxy->x[i - 1] - dx) / fabs(dx) > EPS_LON)
            break;
    gxy->regular_x = (i == ni);

    ascending = (gxy->x[ni - 1] > gxy->x[0]);
    if (ascending) {
        for (i = 0; i < ni - 1; ++i)
            if ((gxy->x[i + 1] - gxy->x[i]) < 0.0)
                enkf_quit("non-monotonic X coordinate for a rect grid\n");
    } else {
        for (i = 0; i < ni - 1; ++i)
            if ((gxy->x[i + 1] - gxy->x[i]) > 0.0)
                enkf_quit("non-monotonic X coordinate for a rect grid\n");
    }

    /*
     * y
     */
    gxy->nj = nj;
    gxy->y = y;
    gxy->yc = malloc((nj + 1) * sizeof(double));
    gxy->yc[0] = y[0] * 1.5 - y[1] * 0.5;
    for (i = 1; i < nj + 1; ++i)
        gxy->yc[i] = 2 * y[i - 1] - gxy->yc[i - 1];

    dy = (y[nj - 1] - y[0]) / (double) (nj - 1);
    for (i = 1; i < (int) nj; ++i)
        if (fabs(y[i] - y[i - 1] - dy) / fabs(dy) > EPS_LON)
            break;
    gxy->regular_y = (i == nj);

    ascending = (y[nj - 1] > y[0]);
    if (ascending) {
        for (i = 0; i < nj - 1; ++i)
            if ((y[i + 1] - y[i]) < 0.0)
                enkf_quit("non-monotonic Y coordinate for a rect grid\n");
    } else {
        for (i = 0; i < nj - 1; ++i)
            if ((y[i + 1] - y[i]) > 0.0)
                enkf_quit("non-monotonic Y coordinate for a rect grid\n");
    }

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

    if (fi < 0.0) {
        if (!periodic)
            return NAN;
        else
            fi += (double) (n - 1);
    } else if (fi >= (double) (n - 1)) {
        if (!periodic)
            return NAN;
        else
            fi -= (double) (n - 1);
    }

    return fi;
}
#endif

/** Gets fractional index of a coordinate for a 1D irregular grid.
 * @param n Number of grid nodes (vertical layers)
 * @param v Coordinates of the nodes (layer centres)
 * @param vb Coordinates of the cell/layer boundaries [n + 1]
 * @param x Input coordinate
 * @param periodic Flag for grid periodicity
 * @param lonbase lon range = [lonbase, lonbase + 360)
 * @return Fractional index for `x'
 */
#if defined(ENKF_PREP) || defined(ENKF_CALC)
static double x2fi_irreg(int n, double v[], double vb[], double x, int periodic, double lonbase)
{
    int ascending, i1, i2, imid;

    if (n < 2)
        return NAN;

    if (periodic == 2)
        n += 1;

    if (!isnan(lonbase)) {
        if (x < lonbase)
            x += 360.0;
        else if (x >= lonbase + 360.0)
            x -= 360.0;
    }

    ascending = (v[n - 1] > v[0]);

    if (ascending) {
        if ((x < vb[0]) || x > vb[n])
            return NAN;
    } else {
        if ((x > vb[0]) || x < vb[n])
            return NAN;
    }

    i1 = 0;
    i2 = n - 1;
    if (ascending) {
        while (1) {
            imid = (i1 + i2) / 2;
            if (imid == i1)
                break;
            if (x > vb[imid])
                i1 = imid;
            else
                i2 = imid;
        }
        if (x < vb[i1 + 1])
            return (double) i1 + (x - v[i1]) / (vb[i1 + 1] - vb[i1]);
        else
            return (double) i1 + 0.5 + (x - vb[i1 + 1]) / (vb[i1 + 2] - vb[i1 + 1]);
    } else {
        while (1) {
            imid = (i1 + i2) / 2;
            if (imid == i1)
                break;
            if (x > vb[imid])
                i2 = imid;
            else
                i1 = imid;
        }
        if (x > vb[i1 + 1])
            return (double) i1 + (x - v[i1]) / (vb[i1 + 1] - vb[i1]);
        else
            return (double) i1 + 0.5 + (x - vb[i1 + 1]) / (vb[i1 + 2] - vb[i1 + 1]);
    }
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
        fij[0] = x2fi_irreg(gxy->ni, gxy->x, gxy->xc, x, gxy->parent->periodic_i, gxy->parent->lonbase);
    if (gxy->regular_y)
        fij[1] = x2fi_reg(gxy->nj, gxy->y, y, 0);
    else
        fij[1] = x2fi_irreg(gxy->nj, gxy->y, gxy->yc, y, 0, NAN);
}
#endif

/**
 */
void gxy_rect_destroy(gxy_rect* gxy)
{
    free(gxy->x);
    free(gxy->y);
    free(gxy->xc);
    free(gxy->yc);
    free(gxy);
}
