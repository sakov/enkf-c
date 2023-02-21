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
#include "triangulation.h"
#include "gxy_unstr.h"

struct gxy_unstr {
    triangulation* d;
};

/**
 */
gxy_unstr* gxy_unstr_create(void* grid, triangulation* d)
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

/**
 */
int gxy_unstr_xy2fij(gxy_unstr* gxy, double x, double y, double* fij)
{
    point p = { x, y };
    triangle* t;
    int tid;
    double bc[3];
    int i;

    if (!triangulation_xy2bc(gxy->d, &p, &tid, bc)) {
        fij[0] = NAN;
        fij[1] = NAN;
        fij[2] = NAN;
        return 0;               /* failed */
    }
    t = &gxy->d->triangles[tid];

    for (i = 0; i < 3; ++i) {
        fij[i] = (double) t->vids[i] + bc[i];
        assert((int) fij[i] == t->vids[i]);
    }

    return 1;
}
