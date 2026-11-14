/******************************************************************************
 *
 * File:        gxy_unstr.h        
 *
 * Created:     25/01/2022
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Purpose:     Code for `gxy_unstr' object (handling of unstructured
 *              horizontal grids).
 *
 * Description:
 *
 * Revisions:
 *
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <float.h>
#include "triangulation.h"
#include "gxy_unstr.h"

struct gxy_unstr {
    triangulation* d;
};

/**
 */
gxy_unstr* gxy_unstr_create(triangulation* d)
{
    gxy_unstr* gxy = malloc(sizeof(gxy_unstr));

    gxy->d = d;

    return gxy;
}

/**
 */
void gxy_unstr_destroy(gxy_unstr* gxy)
{
    triangulation_destroy(gxy->d);
    free(gxy);
}

/**
 */
triangulation* gxy_unstr_gettriangulation(gxy_unstr* gxy)
{
    return gxy->d;
}

/**
 */
point* gxy_unstr_getpoint(gxy_unstr* gxy, int i)
{
    if (i < 0 || i >= gxy->d->npoints)
        return NULL;
    return &gxy->d->points[i];
}

/** Maps physical coordinates to "fractional indices" fij[3]. For each vertex
 ** of the triangle containing the point fij[i] = vertex id + barycentric
 ** coordinate, with the barycentric coordinate clamped to [0, 1 - eps] so that
 ** (int) fij[i] always gives the vertex id.
 */
int gxy_unstr_xy2fij(gxy_unstr* gxy, double x, double y, double* fij)
{
    point p = { x, y };
    triangle* t;
    int tid;
    double bc[3];
    int i;

    if (!triangulation_xy2bc(gxy->d, &p, &tid, bc) || !isfinite(bc[0]) || !isfinite(bc[1]) || !isfinite(bc[2])) {
        fij[0] = NAN;
        fij[1] = NAN;
        fij[2] = NAN;
        return 0;               /* failed */
    }
    t = &gxy->d->triangles[tid];

    for (i = 0; i < 3; ++i) {
        int vid = t->vids[i];

        /*
         * eps is relative so that vid + bc never rounds up to vid + 1
         */
        double eps = 4.0 * DBL_EPSILON * ((double) vid + 1.0);
        double w = bc[i];

        if (w < 0.0)
            w = 0.0;
        else if (w > 1.0 - eps)
            w = 1.0 - eps;
        fij[i] = (double) vid + w;
        assert((int) fij[i] == vid);
    }

    return 1;
}
