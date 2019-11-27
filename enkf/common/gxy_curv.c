/******************************************************************************
 *
 * File:        gxy_curv.c
 *
 * Created:     27/11/2019
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description: Contains code for mappings between indices and physical
 *              coordinates of a curvilinear (or, rather, quadrilateral) grid.
 *              The code has been extracted from and replaces that from
 *              gridutils library.
 *
 * Revisions:   
 *
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "kdtree.h"
#include "gxy_curv.h"

#define EPS 1.0e-8
#define EPS_ZERO 1.0e-5

struct gxy_curv {
    int ni;
    int nj;
    double** x;
    double** y;
    int sign;
    kdtree* nodetree;
};

/**
 */
gxy_curv* gxy_curv_create(int ni, int nj, double** x, double** y)
{
    gxy_curv* gxy = malloc(sizeof(gxy_curv));

    gxy->ni = ni;
    gxy->nj = nj;
    gxy->x = x;
    gxy->y = y;
    gxy->sign = 0;

    assert(gxy->x != NULL && gxy->y != NULL);
    {
        double* data[2];

        data[0] = x[0];
        data[1] = y[0];
        gxy->nodetree = kd_create(2);
        kd_insertnodes(gxy->nodetree, ni * nj, data, 1 /* shuffle */ );
    }

    return gxy;
}

/**
 */
void gxy_curv_destroy(gxy_curv* gxy)
{
    free(gxy->x);
    free(gxy->y);
    kd_destroy(gxy->nodetree);
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

/**
 */
static inline int onleft(double x0, double y0, double x1, double y1, double x2, double y2)
{
    double a = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0);

    if (a > 0.0)
        return 1;
    else if (a < 0.0)
        return -1;
    else
        return 0;
}

/**
 */
static int inquadri(double x, double y, double* px, double* py)
{
    int count, i;

    for (i = 0, count = 0; i < 5; ++i) {
        if (py[i] <= y) {
            if (py[i + 1] > y)
                if (onleft(px[i], py[i], px[i + 1], py[i + 1], x, y) > 0)
                    ++count;
        } else {
            if (py[i + 1] <= y)
                if (onleft(px[i], py[i], px[i + 1], py[i + 1], x, y) < 0)
                    --count;
        }
    }

    return count;
}

/**
 */
static int gxy_curv_xy2ij(gxy_curv* gxy, double x, double y, int* iout, int* jout)
{
    double* minmax;
    double pos[2];
    size_t nearest;
    size_t id;
    int i, j, i1, i2, j1, j2;
    double px[5], py[5];
    int success = 0;

    minmax = kd_getminmax(gxy->nodetree);
    if (x < minmax[0] || y < minmax[1] || x > minmax[2] || y > minmax[3])
        return success;

    pos[0] = x;
    pos[1] = y;
    nearest = kd_findnearestnode(gxy->nodetree, pos);
    id = kd_getnodeorigid(gxy->nodetree, nearest);

    j = id / (gxy->ni);
    i = id % (gxy->ni);

    i1 = (i > 0) ? i - 1 : i;
    i2 = (i < gxy->ni - 1) ? i + 1 : i;
    j1 = (j > 0) ? j - 1 : j;
    j2 = (j < gxy->nj - 1) ? j + 1 : j;

    for (j = j1; j <= j2 - 1; ++j)
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
            px[4] = gxy->x[j][i];
            py[4] = gxy->y[j][i];

            if (inquadri(x, y, px, py)) {
                success = 1;
                *iout = i;
                *jout = j;
                goto finish;
            }
        }

  finish:
    return success;
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
int gxy_curv_xy2fij(gxy_curv* gxy, double x, double y, double* fi, double* fj)
{
    double** gx = gxy->x;
    double** gy = gxy->y;
    int i, j;

    *fi = NAN;
    *fj = NAN;

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

        if (u < 0.0)
            u = 0.0;
        else if (u >= 1.0)
            u = 1.0 - EPS;
        if (v < 0.0)
            v = 0.0;
        else if (v >= 1.0)
            v = 1.0 - EPS;

        *fi = i + u;
        *fj = j + v;
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

    /*
     * Trim I to range 0 to nce1 
     */
    if (fi < 0 || fi >= gxy->ni || fj < 0 || fj >= gxy->nj)
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
