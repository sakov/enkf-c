/******************************************************************************
 *
 * File:        grid.c        
 *
 * Created:     12/2012
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description:
 *
 * Revisions:   15/02/2016 PS x2fi_reg() now returns NaN if estimated fi either
 *                < 0.0 or > (double) (n - 1). (In the previous version indices
 *                -0.5 < fi <= 0.0 were mapped to 0.0, and indices
 *                (double) (n - 1) < fi < (double) n - 0.5 were mapped to
 *                (double) (n - 1).)
 *              11/12/2017 PS Added struct gz_hybrid.
 *              25/01/2018 PS Modified z2fk_basic() to handle the case
 *                zt[i] != 0.5 * (zc[i] + zc[i + 1])
 *              30/01/2018 PS Modified fk2z() to match the above changes.
 *              15/06/2018 PS Added st and sc arrays to struct gz_sigma to hold
 *                the possibly stretched sigma coordinates of layers, and
 *                modified the rest of the code accordingly
 *
 *****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include "ncw.h"
#include "ncutils.h"
#include "definitions.h"
#include "utils.h"
#include "grid.h"
#include "gxy_curv.h"
#include "gridprm.h"

#define EPS_LON 1.0e-3
#define EPS_IJ 1.0e-3
#define EPS_Z 1.0e-3

#define GRIDVDIR_NONE 0
#define GRIDVDIR_FROMSURF 0
#define GRIDVDIR_TOSURF 1

#define GRID_INC 10

typedef struct {
    int ni;
    int nj;
} gxy_none;

typedef struct {
    int ni;
    int nj;
    int periodic_i;
    int regular_i;
    int regular_y;

    double* x;
    double* y;
    double* xc;
    double* yc;
} gxy_simple;

typedef struct {
    int nk;
    int vdirection;
    double* zt;
    double* zc;
} gz_z;

/*
 * This structure manages the "new" ROMS vertical coordinate as described in
 * https://www.myroms.org/wiki/Vertical_S-coordinate, Eq. 2 and by Shchepetkin
 * in https://www.myroms.org/forum/viewtopic.php?f=20&t=2189. With hc = 0 it
 * turns into "pure" sigma. To be consistent with the rest of the code we use
 * ct and cc for the stretching function C at the layer centres and corners
 * (these arrays correspond to variables Cs_r and Cs_w in ROMS restarts).
 *
 * Some recent stretching options in ROMS also require stretching of the sigma
 * coordinate, which previously used to be uniform -1 <= s < = 0. This why there
 * are st and sc arrays added.
 */
typedef struct {
    int ni;
    int nj;
    int nk;
    int vdirection;
    double hc;
    double* ct;
    double* cc;
    double* st;
    double* sc;

    double fi_prev;
    double fj_prev;
    double* zt;
    double* zc;
} gz_sigma;

/*
 * See e.g. Eckerman (2009), Eq. 6 
 * https://journals.ametsoc.org/doi/pdf/10.1175/2008MWR2537.1
 */
typedef struct {
    int ni;
    int nj;
    int nk;
    int vdirection;
    double* at;
    double* ac;
    double* bt;
    double* bc;
    float** p1;
    float** p2;

    double fi_prev;
    double fj_prev;
    double* pt;
    double* pc;
} gz_hybrid;

struct grid {
    char* name;
    int id;
    int htype;                  /* horizontal type */
    int vtype;                  /* vertical type */

    void* gridnodes_xy;         /* (the structure is defined by `htype') */
    double lonbase;             /* (lon range = [lonbase, lonbase + 360)] */

    void* gridnodes_z;

    /*
     * `numlevels' can hold either the number of levels (z-model) or the
     * land mask (sigma-model)
     */
    int** numlevels;
    float** depth;

    /*
     * `stride' for calculating ensemble transforms. "0" means to use the
     * common value defined in the top prm file. 
     */
    int stride;

    /*
     * "Spread factor", normally set to 1. Introduced to adjust relative spread
     * of the ocean and atmospheric parts in climate models. Applies to all
     * variables of the grid.
     */
    double sfactor;

    /*
     * Vertical intervals for observation statistics
     */
    int nzints;
    zint* zints;

    char* domainname;

    /*
     * used for calculating forecast obs with finite footprint
     */
    kdtree* nodetreeXYZ;
};

/**
 */
static gxy_none* gxy_none_create(int ni, int nj)
{
    gxy_none* gxy = malloc(sizeof(gxy_none));

    gxy->ni = ni;
    gxy->nj = nj;

    return gxy;
}

/**
 */
static gxy_simple* gxy_simple_create(int ni, int nj, double* x, double* y)
{
    gxy_simple* gxy = malloc(sizeof(gxy_simple));
    int i, ascending;
    double dx, dy;

    assert(ni >= 2 && nj >= 2);

    /*
     * x
     */
    gxy->ni = ni;
    gxy->x = x;
    gxy->xc = malloc((ni + 1) * sizeof(double));
    gxy->xc[0] = x[0] * 1.5 - x[1] * 0.5;
    for (i = 1; i < ni + 1; ++i)
        gxy->xc[i] = 2 * x[i - 1] - gxy->xc[i - 1];

    if (fabs(fmod(gxy->x[ni - 1] - gxy->x[0] + EPS_LON / 2.0, 360.0)) < EPS_LON)
        gxy->periodic_i = 1;    /* closed grid */
    else if (fabs(fmod(2.0 * gxy->x[ni - 1] - gxy->x[ni - 2] - gxy->x[0] + EPS_LON / 2.0, 360.0)) < EPS_LON) {
        gxy->periodic_i = 2;    /* non-closed grid (used e.g. by MOM) */
        gxy->x = realloc(gxy->x, (ni + 1) * sizeof(double));
        gxy->x[ni] = 2.0 * gxy->x[ni - 1] - gxy->x[ni - 2];
        gxy->xc = realloc(gxy->xc, (ni + 2) * sizeof(double));
        gxy->xc[ni + 1] = 2.0 * gxy->xc[ni] - gxy->xc[ni - 1];
    } else
        gxy->periodic_i = 0;

    dx = (gxy->x[ni - 1] - gxy->x[0]) / (double) (ni - 1);
    for (i = 1; i < (int) ni; ++i)
        if (fabs(gxy->x[i] - gxy->x[i - 1] - dx) / fabs(dx) > EPS_LON)
            break;
    gxy->regular_i = (i == ni);

    ascending = (gxy->x[ni - 1] > gxy->x[0]);
    if (ascending) {
        for (i = 0; i < ni - 1; ++i)
            if ((gxy->x[i + 1] - gxy->x[i]) < 0.0)
                enkf_quit("non-monotonic X coordinate for a simple grid\n");
    } else {
        for (i = 0; i < ni - 1; ++i)
            if ((gxy->x[i + 1] - gxy->x[i]) > 0.0)
                enkf_quit("non-monotonic X coordinate for a simple grid\n");
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
                enkf_quit("non-monotonic Y coordinate for a simple grid\n");
    } else {
        for (i = 0; i < nj - 1; ++i)
            if ((y[i + 1] - y[i]) > 0.0)
                enkf_quit("non-monotonic Y coordinate for a simple grid\n");
    }

    return gxy;
}

/**
 */
void gxy_simple_destroy(gxy_simple* gxy)
{
    free(gxy->x);
    free(gxy->y);
    free(gxy->xc);
    free(gxy->yc);
    free(gxy);
}

#define EPSZ 0.001

/**
 */
static gz_z* gz_z_create(grid* g, int nk, double* z, int nkc, double* zc, char* vdirection)
{
    gz_z* gz = malloc(sizeof(gz_z));
    int i;

    if (strncasecmp(vdirection, "FROMSURF", 8) == 0)
        gz->vdirection = GRIDVDIR_FROMSURF;
    else if (strncasecmp(vdirection, "TOSURF", 6) == 0)
        gz->vdirection = GRIDVDIR_TOSURF;
    else
        enkf_quit("programming error");

    gz->zt = z;
    gz->nk = nk;
    /*
     * check monotonicity
     */
    for (i = 0; i < nk - 2; ++i)
        if ((z[i + 1] - z[i]) * (z[i + 2] - z[i + 1]) < 0.0)
            enkf_quit("%s: non-monotonic Z grid", g->name);
    /*
     * change sign if negative
     */
    {
        double sum = 0.0;

        for (i = 0; i < nk; ++i)
            sum += z[i];
        if (sum < 0.0)
            for (i = 0; i < nk; ++i)
                z[i] = -z[i];
    }

    gz->zc = malloc((nk + 1) * sizeof(double));
    if (zc == NULL) {
        if (z[nk - 1] >= z[0]) {
            gz->zc[0] = 0.0;
            for (i = 1; i <= nk; ++i)
                gz->zc[i] = 2.0 * z[i - 1] - gz->zc[i - 1];
        } else {
            gz->zc[nk] = 0.0;
            for (i = nk - 1; i >= 0; --i)
                gz->zc[i] = 2.0 * z[i] - gz->zc[i + 1];
        }
    } else {
        if (nkc == nk) {        /* e.g. in MOM */
            memcpy(gz->zc, zc, nk * sizeof(double));
            gz->zc[nk] = 2.0 * z[nk - 1] - zc[nk - 1];
        } else if (nkc == nk + 1)
            memcpy(gz->zc, zc, nkc * sizeof(double));
        free(zc);
    }

    if (gz->vdirection == GRIDVDIR_FROMSURF)
        assert(gz->zt[0] <= gz->zt[nk - 1] && gz->zc[0] < gz->zc[nk]);
    else if (gz->vdirection == GRIDVDIR_TOSURF)
        assert(gz->zt[0] >= gz->zt[nk - 1] && gz->zc[0] > gz->zc[nk]);
    for (i = 0; i < nk; ++i)
        assert((gz->zt[i] - gz->zc[i]) * (gz->zc[i + 1] - gz->zt[i]) > 0.0);

    return gz;
}

/**
 */
static gz_sigma* gz_sigma_create(grid* g, int ni, int nj, int nk, double* st, double* sc, double* ct, double* cc, double hc, char* vdirection)
{
    gz_sigma* gz = malloc(sizeof(gz_sigma));
    int i;

    gz->ni = ni;
    gz->nj = nj;
    gz->nk = nk;
    gz->hc = hc;
    if (strncasecmp(vdirection, "FROMSURF", 8) == 0)
        gz->vdirection = GRIDVDIR_FROMSURF;
    else if (strncasecmp(vdirection, "TOSURF", 6) == 0)
        gz->vdirection = GRIDVDIR_TOSURF;
    else
        enkf_quit("programming error");

    assert(ct != NULL);
    /*
     * check monotonicity
     */
    for (i = 0; i < nk - 2; ++i)
        if ((ct[i + 1] - ct[i]) * (ct[i + 2] - ct[i + 1]) < 0.0)
            enkf_quit("%s: non-monotonic Cs_t", g->name);

    gz->ct = ct;

    if (cc != NULL)
        gz->cc = cc;
    else
        gz->cc = malloc((nk + 1) * sizeof(double));
    /*
     * change sign if negative
     */
    {
        double sum = 0.0;

        for (i = 0; i < nk; ++i)
            sum += ct[i];
        if (sum < 0.0) {
            for (i = 0; i < nk; ++i)
                ct[i] = -ct[i];
            if (cc != NULL)
                for (i = 0; i <= nk; ++i)
                    cc[i] = -cc[i];
        }
    }
    /*
     * build Cs_w if necessary
     */
    if (cc == NULL) {
        if (ct[0] > ct[nk - 1]) {
            gz->cc[0] = 1.0;
            gz->cc[nk - 1] = 0.0;
        } else {
            gz->cc[0] = 0.0;
            gz->cc[nk - 1] = 1.0;
        }
        for (i = 1; i < nk; ++i)
            gz->cc[i] = (gz->ct[i - 1] + gz->ct[i]) / 2.0;
    }

    /*
     * read sigma arrays if necessary (assume uniform stretching by default)
     */
    if (st == NULL) {
        gz->st = malloc(nk * sizeof(double));
        if (gz->vdirection == GRIDVDIR_FROMSURF)
            for (i = 0; i < nk; ++i)
                gz->st[i] = ((double) i + 0.5) / (double) gz->nk;
        else
            for (i = 0; i < nk; ++i)
                gz->st[i] = 1.0 - ((double) i + 0.5) / (double) gz->nk;
        assert(sc == NULL);
        gz->sc = malloc((nk + 1) * sizeof(double));
        if (gz->vdirection == GRIDVDIR_FROMSURF)
            for (i = 0; i <= nk; ++i)
                gz->sc[i] = (double) i / (double) gz->nk;
        else
            for (i = 0; i <= nk; ++i)
                gz->sc[i] = 1.0 - (double) i / (double) gz->nk;
    } else {
        double sum = 0.0;

        gz->st = st;
        assert(sc != NULL);
        gz->sc = sc;
        /*
         * check monotonicity
         */
        for (i = 0; i < nk - 2; ++i)
            if ((st[i + 1] - st[i]) * (st[i + 2] - st[i + 1]) < 0.0)
                enkf_quit("%s: non-monotonic s_rho coordinate (entry SVARNAME in the grid parameter file)", g->name);
        /*
         * change sign if negative
         */
        for (i = 0; i < nk; ++i)
            sum += st[i];
        if (sum < 0.0)
            for (i = 0; i < nk; ++i)
                st[i] = -st[i];
    }
    /*
     * build s_w if necessary
     */
    if (sc == NULL) {
        if (st[0] > st[nk - 1]) {
            gz->sc[0] = 1.0;
            gz->sc[nk - 1] = 0.0;
        } else {
            gz->sc[0] = 0.0;
            gz->sc[nk - 1] = 1.0;
        }
        for (i = 1; i < nk; ++i)
            gz->sc[i] = (gz->st[i - 1] + gz->st[i]) / 2.0;
    } else {
        double sum = 0.0;

        /*
         * check monotonicity
         */
        for (i = 0; i < nk - 1; ++i)
            if ((sc[i + 1] - sc[i]) * (sc[i + 2] - sc[i + 1]) < 0.0)
                enkf_quit("%s: non-monotonic s_w coordinate (entry SCVARNAME in the grid parameter file)", g->name);
        /*
         * change sign if negative
         */
        for (i = 0; i < nk; ++i)
            sum += sc[i];
        if (sum < 0.0)
            for (i = 0; i <= nk; ++i)
                sc[i] = -sc[i];
    }

    gz->zt = malloc(nk * sizeof(double));
    gz->zc = malloc((nk + 1) * sizeof(double));
    gz->fi_prev = NAN;
    gz->fj_prev = NAN;

    if (gz->vdirection == GRIDVDIR_FROMSURF) {
        assert(gz->ct[0] <= gz->ct[nk - 1] && gz->cc[0] < gz->cc[nk]);
        assert(gz->st[0] <= gz->st[nk - 1] && gz->sc[0] < gz->sc[nk]);
    } else if (gz->vdirection == GRIDVDIR_TOSURF) {
        assert(gz->ct[0] >= gz->ct[nk - 1] && gz->cc[0] > gz->cc[nk]);
        assert(gz->st[0] >= gz->st[nk - 1] && gz->sc[0] > gz->sc[nk]);
    }
    for (i = 0; i < nk; ++i)
        assert((gz->ct[i] - gz->cc[i]) * (gz->cc[i + 1] - gz->ct[i]) > 0.0);
    for (i = 0; i < nk; ++i)
        assert((gz->st[i] - gz->sc[i]) * (gz->sc[i + 1] - gz->st[i]) > 0.0);

    return gz;
}

/**
 */
static gz_hybrid* gz_hybrid_create(int ni, int nj, int nk, double* a, double* b, double* ac, double* bc, float** p1, float** p2, char* vdirection)
{
    gz_hybrid* gz = malloc(sizeof(gz_hybrid));
    int i;

    if (strncasecmp(vdirection, "FROMSURF", 8) == 0)
        gz->vdirection = GRIDVDIR_FROMSURF;
    else if (strncasecmp(vdirection, "TOSURF", 6) == 0)
        gz->vdirection = GRIDVDIR_TOSURF;
    else
        enkf_quit("programming error");

    gz->ni = ni;
    gz->nj = nj;
    gz->nk = nk;
    gz->at = a;
    gz->bt = b;
    gz->p1 = p1;
    gz->p2 = p2;

    if (ac != NULL) {
        assert(bc != NULL);

        gz->ac = a;
        gz->bc = b;
    } else {
        assert(bc == NULL);

        gz->ac = malloc((nk + 1) * sizeof(double));
        gz->bc = malloc((nk + 1) * sizeof(double));

        gz->ac[0] = 1.5 * gz->at[0] - 0.5 * gz->at[1];
        if (gz->ac[0] < 0.0)
            gz->ac[0] = 0.0;
        gz->ac[nk] = 1.5 * gz->at[nk - 1] - 0.5 * gz->at[nk - 2];
        if (gz->ac[nk] < 0.0)
            gz->ac[nk] = 0.0;
        for (i = 1; i < nk; ++i)
            gz->ac[i] = (gz->at[i - 1] + gz->at[i]) / 2.0;
        gz->bc[0] = 1.5 * gz->bt[0] - 0.5 * gz->bt[1];
        if (gz->bc[0] < 0.0)
            gz->bc[0] = 0.0;
        gz->bc[nk] = 1.5 * gz->bt[nk - 1] - 0.5 * gz->bt[nk - 2];
        if (gz->bc[nk] < 0.0)
            gz->bc[nk] = 0.0;
        for (i = 1; i < nk; ++i)
            gz->bc[i] = (gz->bt[i - 1] + gz->bt[i]) / 2.0;
    }

    gz->pt = malloc(nk * sizeof(double));
    gz->pc = malloc((nk + 1) * sizeof(double));
    gz->fi_prev = NAN;
    gz->fj_prev = NAN;

    return gz;
}

/**
 */
static void gz_z_destroy(gz_z * gz)
{
    free(gz->zt);
    free(gz->zc);
    free(gz);
}

/**
 */
static void gz_sigma_destroy(gz_sigma* gz)
{
    free(gz->ct);
    free(gz->cc);
    free(gz->zt);
    free(gz->zc);
    free(gz->st);
    free(gz->sc);
    free(gz);
}

/**
 */
static void gz_hybrid_destroy(gz_hybrid* gz)
{
    free(gz->at);
    free(gz->bt);
    free(gz->ac);
    free(gz->bc);
    free(gz->p1);
    free(gz->p2);
    free(gz->pt);
    free(gz->pc);
    free(gz);
}

/**
 */
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

/** Converts fractional indices to physical coordinates for "simple" horizontal
 ** grid.
 */
static void gs_fij2xy(void* p, double fi, double fj, double* x, double* y)
{
    gxy_simple* gxy = (gxy_simple*) ((grid*) p)->gridnodes_xy;

    *x = fi2x(gxy->ni, gxy->x, fi, gxy->periodic_i);
    *y = fi2x(gxy->nj, gxy->y, fj, 0);
}

/** Gets fractional index of a coordinate for a 1D irregular grid.
 * @param n Number of grid nodes (vertical layers)
 * @param v Coordinates of the nodes (layer centres)
 * @param vb Coordinates of the cell/layer boundaries [n + 1]
 * @param x Input coordinate
 * @param periodic Flag for grid periodicity
 * @param lonbase lon range = [lonbase, lonbase + 360)
 * @return Fractional index for `x'
 */
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

/** Converts physical coordinates to fractional indices for "simple" horizontal
 ** grid.
 */
static void gs_xy2fij(void* p, double x, double y, double* fi, double* fj)
{
    gxy_simple* nodes = (gxy_simple*) ((grid*) p)->gridnodes_xy;
    double lonbase = ((grid*) p)->lonbase;

    if (nodes->regular_i)
        *fi = x2fi_reg(nodes->ni, nodes->x, x, nodes->periodic_i);
    else
        *fi = x2fi_irreg(nodes->ni, nodes->x, nodes->xc, x, nodes->periodic_i, lonbase);
    if (nodes->regular_y)
        *fj = x2fi_reg(nodes->nj, nodes->y, y, 0);
    else
        *fj = x2fi_irreg(nodes->nj, nodes->y, nodes->yc, y, 0, NAN);
}

/** Maps vertical coordinate to fractional cell index.
 */
static double z2fk_basic(int n, double* zt, double* zc, double z)
{
    int ascending, i1, i2, imid;

    ascending = (zt[n - 1] > zt[0]) ? 1 : 0;

    if (ascending) {
        if (z < zc[0])
            return 0.0;
        if (z > zc[n])
            return NAN;
    } else {
        if (z < zc[n])
            return 0.0;
        if (z > zc[0])
            return NAN;
    }

    i1 = 0;
    i2 = n - 1;
    if (ascending) {
        while (1) {
            imid = (i1 + i2) / 2;
            if (imid == i1)
                break;
            if (z > zt[imid])
                i1 = imid;
            else
                i2 = imid;
        }
        if (z < zc[i1 + 1])
            return (double) i1 + 0.5 * (z - zt[i1]) / (zc[i1 + 1] - zt[i1]);
        else
            return (double) i1 + 0.5 + 0.5 * (z - zc[i1 + 1]) / (zt[i1 + 1] - zc[i1 + 1]);
    } else {
        while (1) {
            imid = (i1 + i2) / 2;
            if (imid == i1)
                break;
            if (z > zt[imid])
                i2 = imid;
            else
                i1 = imid;
        }
        if (z > zc[i1 + 1])
            return (double) i1 + 0.5 * (z - zt[i1]) / (zc[i1 + 1] - zt[i1]);
        else
            return (double) i1 + 0.5 + 0.5 * (z - zc[i1 + 1]) / (zt[i1 + 1] - zc[i1 + 1]);
    }
}

/**
 */
static void gz_z_z2fk(void* p, double fi, double fj, double z, double* fk)
{
    grid* g = (grid*) p;
    gz_z* gz = g->gridnodes_z;

    *fk = z2fk_basic(gz->nk, gz->zt, gz->zc, z);
}

/**
 */
static void gz_sigma_z2fk(void* p, double fi, double fj, double z, double* fk)
{
    grid* g = (grid*) p;
    gz_sigma* gz = g->gridnodes_z;

    if (g->depth == NULL)
        enkf_quit("%s: DEPTHVARNAME must be entered for SIGMA grids to assimilate subsurface obs", g->name);

    if (isnan(gz->fi_prev) || fabs(fi - gz->fi_prev) > EPS_IJ || fabs(fj - gz->fj_prev) > EPS_IJ) {
        double h = (double) interpolate2d(fi, fj, gz->ni, gz->nj, g->depth, g->numlevels, grid_isperiodic_i(g));
        int i;

        if (gz->hc != 0.0) {
            for (i = 0; i < gz->nk; ++i)
                gz->zt[i] = h * (gz->hc * gz->st[i] + h * gz->ct[i]) / (gz->hc + h);
            for (i = 0; i <= gz->nk; ++i)
                gz->zc[i] = h * (gz->hc * gz->sc[i] + h * gz->cc[i]) / (gz->hc + h);
        } else {
            for (i = 0; i < gz->nk; ++i)
                gz->zt[i] = h * gz->ct[i];
            for (i = 0; i <= gz->nk; ++i)
                gz->zc[i] = h * gz->cc[i];
        }
    }

    *fk = z2fk_basic(gz->nk, gz->zt, gz->zc, z);
}

/**
 */
static void gz_hybrid_z2fk(void* p, double fi, double fj, double z, double* fk)
{
    grid* g = (grid*) p;
    gz_hybrid* gz = (gz_hybrid*) g->gridnodes_z;

    if (isnan(gz->fi_prev) || fabs(fi - gz->fi_prev) > EPS_IJ || fabs(fj - gz->fj_prev) > EPS_IJ) {

        double p1 = interpolate2d(fi, fj, gz->ni, gz->nj, gz->p1, g->numlevels, grid_isperiodic_i(g));
        double p2 = interpolate2d(fi, fj, gz->ni, gz->nj, gz->p2, g->numlevels, grid_isperiodic_i(g));
        int nk = gz->nk;
        int i;

        for (i = 0; i < nk; ++i)
            gz->pt[i] = gz->at[i] + gz->bt[i] * (p1 - p2);
        for (i = 0; i <= nk; ++i)
            gz->pc[i] = gz->ac[i] + gz->bc[i] * (p1 - p2);
        gz->fi_prev = fi;
        gz->fj_prev = fj;
    }
    *fk = z2fk_basic(gz->nk, gz->pt, gz->pc, z);
}

/**
 */
#if !defined(ENKF_UPDATE)
static void grid_setlonbase(grid* g)
{
    double xmin = DBL_MAX;
    double xmax = -DBL_MAX;

    if (g->htype == GRIDHTYPE_LATLON) {
        gxy_simple* gxy = (gxy_simple*) g->gridnodes_xy;
        double* x = gxy->x;
        int ni = gxy->ni;

        if (gxy->periodic_i == 2)       /* periodic non-closed grid */
            ni++;

        if (xmin > x[0])
            xmin = x[0];
        if (xmin > x[ni - 1])
            xmin = x[ni - 1];
        if (xmax < x[0])
            xmax = x[0];
        if (xmax < x[ni - 1])
            xmax = x[ni - 1];
    } else if (g->htype == GRIDHTYPE_CURVILINEAR) {
        gxy_curv* gxy = (gxy_curv*) g->gridnodes_xy;
        double** x = gxy_curv_getx(gxy);
        int ni = gxy_curv_getni(gxy);
        int nj = gxy_curv_getnj(gxy);
        int i, j;

        for (j = 0; j < nj; ++j) {
            for (i = 0; i < ni; ++i) {
                if (xmin > x[j][i])
                    xmin = x[j][i];
                if (xmax < x[j][i])
                    xmax = x[j][i];
            }
        }
    }
    if (xmin >= 0.0 && xmax <= 360.0)
        g->lonbase = 0.0;
    else if (xmin >= -180.0 && xmax <= 180.0)
        g->lonbase = -180.0;
    else
        g->lonbase = xmin;
}
#endif

/**
 */
static void grid_sethgrid(grid* g, int htype, int ni, int nj, void* x, void* y)
{
    g->htype = htype;
    if (htype == GRIDHTYPE_NONE)
        g->gridnodes_xy = gxy_none_create(ni, nj);
    else if (htype == GRIDHTYPE_LATLON)
        g->gridnodes_xy = gxy_simple_create(ni, nj, x, y);
#if defined(ENKF_PREP) || defined(ENKF_CALC)
    else if (htype == GRIDHTYPE_CURVILINEAR)
        g->gridnodes_xy = gxy_curv_create(ni, nj, x, y, g->numlevels);
#endif
    else
        enkf_quit("programming error");
#if !defined(ENKF_UPDATE)
    grid_setlonbase(g);
#endif
}

/**
 */
grid* grid_create(void* p, int id)
{
    gridprm* prm = (gridprm*) p;
    grid* g = calloc(1, sizeof(grid));
    char* fname = prm->fname;
    int ncid;
    int varid_x, varid_y;
    int ndims_x, ndims_y;
    size_t ni, nj, nk;
    int varid_depth, varid_numlevels;

    g->name = strdup(prm->name);
    g->id = id;
    g->domainname = strdup(prm->domainname);    /* ("Default" by default) */
    g->vtype = gridprm_getvtype(prm);
    if (prm->stride != 0)
        g->stride = prm->stride;
    g->sfactor = prm->sfactor;
    g->nzints = prm->nzints;
    if (prm->nzints > 0) {
        g->zints = malloc(g->nzints * sizeof(zint));
        memcpy(g->zints, prm->zints, g->nzints * sizeof(zint));
    }

    ncw_open(fname, NC_NOWRITE, &ncid);
    ncw_inq_varid(ncid, prm->xvarname, &varid_x);
    ncw_inq_varndims(ncid, varid_x, &ndims_x);
    ncw_inq_varid(ncid, prm->yvarname, &varid_y);
    ncw_inq_varndims(ncid, varid_y, &ndims_y);

    /*
     * set horizontal grid
     */
    if (ndims_x == 1 && ndims_y == 1) {
        int dummy;
        double* x;
        double* y;

        ncw_inq_vardims(ncid, varid_x, 1, &dummy, &ni);
        ncw_inq_vardims(ncid, varid_y, 1, &dummy, &nj);

        x = malloc(ni * sizeof(double));
        y = malloc(nj * sizeof(double));

        ncw_get_var_double(ncid, varid_x, x);
        ncw_get_var_double(ncid, varid_y, y);

        grid_sethgrid(g, GRIDHTYPE_LATLON, ni, nj, x, y);
    } else if (ndims_x == 2 && ndims_y == 2) {
#if defined(ENKF_UPDATE)
        int dummy;
        size_t dimlen[2];

        ncw_inq_vardims(ncid, varid_x, 2, &dummy, dimlen);
        ncw_check_vardims(ncid, varid_y, 2, dimlen);
        nj = dimlen[0];
        ni = dimlen[1];

        grid_sethgrid(g, GRIDHTYPE_NONE, ni, nj, NULL, NULL);
#else
        int dummy;
        double** x;
        double** y;
        size_t dimlen[2];

        ncw_inq_vardims(ncid, varid_x, 2, &dummy, dimlen);
        ncw_check_vardims(ncid, varid_y, 2, dimlen);
        nj = dimlen[0];
        ni = dimlen[1];

        x = alloc2d(nj, ni, sizeof(double));
        y = alloc2d(nj, ni, sizeof(double));

        ncu_readvardouble(ncid, varid_x, ni * nj, x[0]);
        ncu_readvardouble(ncid, varid_y, ni * nj, y[0]);

        grid_sethgrid(g, GRIDHTYPE_CURVILINEAR, ni, nj, x, y);
#endif
    } else
        enkf_quit("%s: could not determine the horizontal grid type", fname);

    /*
     * set vertical grid
     */
    if (prm->depthvarname != NULL) {
        size_t dimlen[2] = { nj, ni };

        g->depth = alloc2d(nj, ni, sizeof(float));
        ncw_inq_varid(ncid, prm->depthvarname, &varid_depth);
        ncw_check_vardims(ncid, varid_depth, 2, dimlen);
        ncw_get_var_float(ncid, varid_depth, g->depth[0]);
    }
    if (g->vtype == GRIDVTYPE_NONE);
    else if (g->vtype == GRIDVTYPE_Z) {
        int varid;
        int dummy;
        double* z = NULL;
        double* zc = NULL;
        size_t nkc = 0;

        ncw_inq_varid(ncid, prm->zvarname, &varid);
        ncw_inq_vardims(ncid, varid, 1, &dummy, &nk);
        z = malloc(nk * sizeof(double));
        ncw_get_var_double(ncid, varid, z);
        if (prm->zcvarname != NULL) {
            ncw_inq_varid(ncid, prm->zcvarname, &varid);
            ncw_inq_vardims(ncid, varid, 1, &dummy, &nkc);
            /*
             * (nkc = nk in MOM)
             */
            assert(nkc == nk || nkc == nk + 1);
            zc = malloc(nkc * sizeof(double));
            ncw_get_var_double(ncid, varid, zc);
        }
        g->gridnodes_z = gz_z_create(g, nk, z, nkc, zc, prm->vdirection);
        if (zc != NULL)
            free(zc);
    } else if (g->vtype == GRIDVTYPE_SIGMA) {
        int varid;
        int dummy;
        double* ct = NULL;
        double* cc = NULL;
        double* st = NULL;
        double* sc = NULL;
        double hc = 0.0;

        ncw_inq_varid(ncid, prm->cvarname, &varid);
        ncw_inq_vardims(ncid, varid, 1, &dummy, &nk);
        ct = malloc(nk * sizeof(double));
        ncw_get_var_double(ncid, varid, ct);
        if (prm->ccvarname != NULL) {
            size_t nkc = nk + 1;

            ncw_inq_varid(ncid, prm->ccvarname, &varid);
            ncw_check_vardims(ncid, varid, 1, &nkc);
            cc = malloc(nkc * sizeof(double));
            ncw_get_var_double(ncid, varid, cc);
        }
        if (prm->hcvarname != NULL) {
            ncw_inq_varid(ncid, prm->hcvarname, &varid);
            ncw_check_varndims(ncid, varid, 0);
            ncw_get_var_double(ncid, varid, &hc);
        }
        if (prm->svarname != NULL) {
            st = malloc(nk * sizeof(double));
            ncw_inq_varid(ncid, prm->svarname, &varid);
            ncw_check_vardims(ncid, varid, 1, &nk);
            ncw_get_var_double(ncid, varid, st);
        }
        if (prm->scvarname != NULL) {
            size_t nkc = nk + 1;

            sc = malloc((nk + 1) * sizeof(double));
            ncw_inq_varid(ncid, prm->scvarname, &varid);
            ncw_check_vardims(ncid, varid, 1, &nkc);
            ncw_get_var_double(ncid, varid, sc);
        }

        g->gridnodes_z = gz_sigma_create(g, ni, nj, nk, st, sc, ct, cc, hc, prm->vdirection);
    } else if (g->vtype == GRIDVTYPE_HYBRID) {
        double* a = NULL;
        double* b = NULL;
        double* ac = NULL;
        double* bc = NULL;
        float** p1 = alloc2d(nj, ni, sizeof(float));
        float** p2 = alloc2d(nj, ni, sizeof(float));
        int varid;
        int dummy;
        size_t dimlen[2];

        ncw_inq_varid(ncid, prm->avarname, &varid);
        ncw_inq_vardims(ncid, varid, 1, &dummy, &nk);
        a = malloc(nk * sizeof(double));
        ncw_get_var_double(ncid, varid, a);

        ncw_inq_varid(ncid, prm->bvarname, &varid);
        ncw_check_vardims(ncid, varid, 1, &nk);
        b = malloc(nk * sizeof(double));
        ncw_get_var_double(ncid, varid, b);

        if (prm->acvarname != NULL) {
            size_t nkc = nk + 1;

            assert(prm->bcvarname != NULL);

            ncw_inq_varid(ncid, prm->acvarname, &varid);
            ncw_check_vardims(ncid, varid, 1, &nkc);
            ac = malloc(nkc * sizeof(double));
            ncw_get_var_double(ncid, varid, ac);

            ncw_inq_varid(ncid, prm->bcvarname, &varid);
            ncw_check_vardims(ncid, varid, 1, &nkc);
            bc = malloc(nk * sizeof(double));
            ncw_get_var_double(ncid, varid, bc);
        }

        ncw_inq_varid(ncid, prm->p1varname, &varid);
        dimlen[0] = nj;
        dimlen[1] = ni;
        ncw_check_vardims(ncid, varid, 2, dimlen);
        ncw_get_var_float(ncid, varid, p1[0]);

        ncw_inq_varid(ncid, prm->p2varname, &varid);
        ncw_check_vardims(ncid, varid, 2, dimlen);
        ncw_get_var_float(ncid, varid, p2[0]);

        g->gridnodes_z = gz_hybrid_create(ni, nj, nk, a, b, ac, bc, p1, p2, prm->vdirection);
    } else
        enkf_quit("not implemented");

    if (prm->levelvarname != NULL) {
        size_t dimlen[2] = { nj, ni };

        g->numlevels = alloc2d(nj, ni, sizeof(int));
        ncw_inq_varid(ncid, prm->levelvarname, &varid_numlevels);
        ncw_check_varndims(ncid, varid_numlevels, 2);
        ncw_check_vardims(ncid, varid_numlevels, 2, dimlen);
        ncw_get_var_int(ncid, varid_numlevels, g->numlevels[0]);
        if (g->vtype == GRIDVTYPE_SIGMA || g->vtype == GRIDVTYPE_HYBRID) {
            int* numlevels = g->numlevels[0];
            int i;

            for (i = 0; i < ni * nj; ++i) {
                if (numlevels[i] == 0)
                    continue;
                else if (numlevels[i] == 1)
                    numlevels[i] = nk;
                else
                    enkf_quit("%s: the grid mask (variable MASKVARNAME = \"%s\") contains a value other than 0 or 1", g->name, prm->levelvarname);
            }
        }
    }
    ncw_close(ncid);

    if (g->numlevels == NULL) {
        g->numlevels = alloc2d(nj, ni, sizeof(int));
        if (g->vtype == GRIDVTYPE_SIGMA || g->vtype == GRIDVTYPE_HYBRID) {
            int i, j;

            for (j = 0; j < nj; ++j)
                for (i = 0; i < ni; ++i)
                    if (g->depth == NULL || g->depth[j][i] > 0.0)
                        g->numlevels[j][i] = nk;
        } else {
            int i, j;

            if (g->depth != NULL) {
                for (j = 0; j < nj; ++j) {
                    for (i = 0; i < ni; ++i) {
                        double depth = g->depth[j][i];
                        double fk = NAN;

                        if (depth > 0.0) {
                            gz_z_z2fk(g, j, i, depth, &fk);
                            g->numlevels[j][i] = ceil(fk + 0.5);
                        }
                    }
                }
            } else {
                for (j = 0; j < nj; ++j)
                    for (i = 0; i < ni; ++i)
                        g->numlevels[j][i] = nk;
            }
        }
    }

    gridprm_print(prm, "    ");
    grid_print(g, "    ");

    return g;
}

/**
 */
void grid_destroy(grid* g)
{
    free(g->name);
    free(g->domainname);
    if (g->htype == GRIDHTYPE_NONE)
        free(g->gridnodes_xy);
    else if (g->htype == GRIDHTYPE_LATLON)
        gxy_simple_destroy(g->gridnodes_xy);
#if defined(ENKF_PREP) || defined(ENKF_CALC)
    else if (g->htype == GRIDHTYPE_CURVILINEAR)
        gxy_curv_destroy(g->gridnodes_xy);
#endif
    else
        enkf_quit("programming_error");
    if (g->gridnodes_z != NULL) {
        if (g->vtype == GRIDVTYPE_Z)
            gz_z_destroy(g->gridnodes_z);
        else if (g->vtype == GRIDVTYPE_SIGMA)
            gz_sigma_destroy(g->gridnodes_z);
        else if (g->vtype == GRIDVTYPE_HYBRID)
            gz_hybrid_destroy(g->gridnodes_z);
        else
            enkf_quit("not implemented");
    }
    if (g->numlevels != NULL)
        free(g->numlevels);
    if (g->depth != NULL)
        free(g->depth);
    if (g->nzints > 0)
        free(g->zints);
#if defined(ENKF_CALC)
    if (g->nodetreeXYZ != NULL)
        kd_destroy(g->nodetreeXYZ);
#endif

    free(g);
}

/**
 */
void grid_print(grid* g, char offset[])
{
    int ni, nj, nk;

    enkf_printf("%sgrid info:\n", offset);
    switch (g->htype) {
    case GRIDHTYPE_LATLON:
        enkf_printf("%s  hor type = LATLON\n", offset);
        enkf_printf("%s  periodic by X = %s\n", offset, grid_isperiodic_i(g) ? "yes" : "no");
        break;
    case GRIDHTYPE_CURVILINEAR:
        enkf_printf("%s  hor type = CURVILINEAR\n", offset);
        break;
    default:
        enkf_printf("%s  h type = NONE\n", offset);
    }
    grid_getsize(g, &ni, &nj, &nk);
    enkf_printf("%s  dims = %d x %d x %d\n", offset, ni, nj, nk);
    if (!isnan(g->lonbase))
        enkf_printf("%s  longitude range = [%.3f, %.3f]\n", offset, g->lonbase, g->lonbase + 360.0);
    else
        enkf_printf("%s  longitude range = any\n", offset);
    switch (g->vtype) {
    case GRIDVTYPE_NONE:
        enkf_printf("%s  v type = NONE\n", offset);
        break;
    case GRIDVTYPE_Z:
        enkf_printf("%s  v type = Z\n", offset);
        break;
    case GRIDVTYPE_SIGMA:
        enkf_printf("%s  v type = SIGMA\n", offset);
        break;
    case GRIDVTYPE_HYBRID:
        enkf_printf("%s  v type = HYBRID\n", offset);
        break;
    default:
        enkf_printf("%s  v type = UNDEFINED\n", offset);
    }
    if (g->vtype == GRIDVTYPE_NONE);
    else if (g->vtype == GRIDVTYPE_Z) {
        gz_z* nodes = g->gridnodes_z;

        enkf_printf("%s  v dir = %s\n", offset, (nodes->vdirection == GRIDVDIR_FROMSURF) ? "FROMSURF" : "TOSURF");
    } else if (g->vtype == GRIDVTYPE_SIGMA) {
        gz_sigma* nodes = g->gridnodes_z;

        enkf_printf("%s  hc = %f\n", offset, nodes->hc);
        enkf_printf("%s  v dir = %s\n", offset, (nodes->vdirection == GRIDVDIR_FROMSURF) ? "FROMSURF" : "TOSURF");
    } else if (g->vtype == GRIDVTYPE_HYBRID) {
        gz_hybrid* nodes = g->gridnodes_z;

        enkf_printf("%s  v dir = %s\n", offset, (nodes->vdirection == GRIDVDIR_FROMSURF) ? "FROMSURF" : "TOSURF");
    } else
        enkf_quit("not implemented");
    if (g->stride != 1)
        enkf_printf("%s  STRIDE = %d\n", offset, g->stride);
    if (g->sfactor != 1.0)
        enkf_printf("%s  SFACTOR = %f\n", offset, g->sfactor);
}

/**
 */
void grid_describeprm(void)
{
    enkf_printf("\n");
    enkf_printf("  Grid parameter file format:\n");
    enkf_printf("\n");
    enkf_printf("    NAME             = <name> [ PREP | CALC ]\n");
    enkf_printf("  [ DOMAIN           = <domain name> ]\n");
    enkf_printf("    VTYPE            = { z | sigma | hybrid }\n");
    enkf_printf("  [ VDIR             = { fromsurf* | tosurf } ]\n");
    enkf_printf("    DATA             = <data file name>\n");
    enkf_printf("    XVARNAME         = <X variable name>\n");
    enkf_printf("    YVARNAME         = <Y variable name>\n");
    enkf_printf("    ZVARNAME         = <Z variable name>                  (z)\n");
    enkf_printf("  [ ZCVARNAME        = <ZC variable name> ]               (z)\n");
    enkf_printf("    CVARNAME         = <Cs_rho variable name>             (sigma)\n");
    enkf_printf("  [ CCVARNAME        = <Cs_w variable name> ]             (sigma)\n");
    enkf_printf("  [ SVARNAME         = <s_rho variable name> ]            (uniform*) (sigma)\n");
    enkf_printf("  [ SCVARNAME        = <s_w variable name> ]              (uniform*) (sigma)\n");
    enkf_printf("  [ HCVARNAME        = <hc variable name> ]               (0.0*) (sigma)\n");
    enkf_printf("  [ DEPTHVARNAME     = <depth variable name> ]            (z | sigma)\n");
    enkf_printf("  [ NUMLEVELSVARNAME = <# of levels variable name> ]      (z)\n");
    enkf_printf("  [ MASKVARNAME      = <land mask variable name> ]        (sigma | hybrid)\n");
    enkf_printf("    AVARNAME         = <A variable name>                  (hybrid)\n");
    enkf_printf("    BVARNAME         = <B variable name>                  (hybrid)\n");
    enkf_printf("  [ ACVARNAME        = <AC variable name> ]               (hybrid)\n");
    enkf_printf("  [ BCVARNAME        = <BC variable name> ]               (hybrid)\n");
    enkf_printf("    P1VARNAME        = <P1 variable name>                 (hybrid)\n");
    enkf_printf("    P2VARNAME        = <P2 variable name>                 (hybrid)\n");
    enkf_printf("  [ STRIDE           = <stride for ensemble transforms> ] (1*)\n");
    enkf_printf("  [ SOBSTRIDE        = <stride for superobing> ]          (1*)\n");
    enkf_printf("  [ SFACTOR          = <spread factor> ]                  (1.0*)\n");
    enkf_printf("  [ ZSTATINTS        = [<z1> <z2>] ... ]\n");
    enkf_printf("\n");
    enkf_printf("  [ <more of the above blocks> ]\n");
    enkf_printf("\n");
    enkf_printf("  Notes:\n");
    enkf_printf("    1. < ... > denotes a description of an entry\n");
    enkf_printf("    2. [ ... ] denotes an optional input\n");
    enkf_printf("    3. (...) is a note\n");
    enkf_printf("    4. * denotes the default value\n");
    enkf_printf("\n");
}

/**
 */
void grid_getsize(grid* g, int* ni, int* nj, int* nk)
{
    if (ni != NULL) {
        if (g->htype == GRIDHTYPE_LATLON) {
            gxy_simple* nodes = (gxy_simple*) g->gridnodes_xy;

            *ni = nodes->ni;
            *nj = nodes->nj;
#if defined(ENKF_PREP) || defined(ENKF_CALC)
        } else if (g->htype == GRIDHTYPE_CURVILINEAR) {
            gxy_curv* nodes = (gxy_curv*) g->gridnodes_xy;

            *ni = gxy_curv_getni(nodes);
            *nj = gxy_curv_getnj(nodes);
#endif
        } else if (g->htype == GRIDHTYPE_NONE) {
            gxy_none* gxy = (gxy_none *) g->gridnodes_xy;

            *ni = gxy->ni;
            *nj = gxy->nj;
        } else
            enkf_quit("programming error");
    }

    if (nk != NULL) {
        if (g->vtype == GRIDVTYPE_Z)
            *nk = ((gz_z *) g->gridnodes_z)->nk;
        else if (g->vtype == GRIDVTYPE_SIGMA)
            *nk = ((gz_sigma*) g->gridnodes_z)->nk;
        else if (g->vtype == GRIDVTYPE_HYBRID)
            *nk = ((gz_hybrid*) g->gridnodes_z)->nk;
        else if (g->vtype == GRIDVTYPE_NONE)
            *nk = 1;
        else
            enkf_quit("programming error");
    }
}

/**
 */
int grid_getsurflayerid(grid* g)
{
    if (g->vtype == GRIDVTYPE_NONE)
        return 0;
    else if (g->vtype == GRIDVTYPE_Z) {
        gz_z* gz = g->gridnodes_z;

        if (gz->vdirection == GRIDVDIR_FROMSURF)
            return 0;
        else if (gz->vdirection == GRIDVDIR_TOSURF)
            return gz->nk - 1;
        else
            enkf_quit("programming error");
    } else if (g->vtype == GRIDVTYPE_SIGMA) {
        gz_sigma* gz = g->gridnodes_z;

        if (gz->vdirection == GRIDVDIR_FROMSURF)
            return 0;
        else if (gz->vdirection == GRIDVDIR_TOSURF)
            return gz->nk - 1;
        else
            enkf_quit("programming error");
    } else if (g->vtype == GRIDVTYPE_HYBRID) {
        gz_hybrid* gz = g->gridnodes_z;

        if (gz->vdirection == GRIDVDIR_FROMSURF)
            return 0;
        else if (gz->vdirection == GRIDVDIR_TOSURF)
            return gz->nk - 1;
        else
            enkf_quit("programming error");
    } else
        enkf_quit("not implemented");

    return -1;
}

/**
 */
char* grid_getname(grid* g)
{
    return g->name;
}

/**
 */
int grid_getid(grid* g)
{
    return g->id;
}

/**
 */
int grid_gethtype(grid* g)
{
    return g->htype;
}

/**
 */
int grid_getvtype(grid* g)
{
    return g->vtype;
}

/**
 */
float** grid_getdepth(grid* g)
{
    return g->depth;
}

/**
 */
int** grid_getnumlevels(grid* g)
{
    return g->numlevels;
}

/**
 */
double grid_getlonbase(grid* g)
{
    return g->lonbase;
}

/**
 */
int grid_getstride(grid* g)
{
    return g->stride;
}

/**
 */
void grid_setstride(grid* g, int stride)
{
    g->stride = stride;
}

/**
 */
double grid_getsfactor(grid* g)
{
    return g->sfactor;
}

/**
 */
void grid_getzints(grid* g, int* nzints, zint* zints[])
{
    *nzints = g->nzints;
    *zints = g->zints;
}

/**
 */
int grid_xy2fij(grid* g, double x, double y, double* fi, double* fj)
{
    int isperiodic_i = grid_isperiodic_i(g);
    int ni = -1, nj = -1;
    int i1, i2, j1, j2;

    if (!isnan(g->lonbase)) {
        if (x < g->lonbase)
            x += 360.0;
        else if (x >= g->lonbase + 360.0)
            x -= 360.0;
    }

    if (g->htype == GRIDHTYPE_LATLON)
        gs_xy2fij(g, x, y, fi, fj);
#if defined(ENKF_PREP) || defined(ENKF_CALC)
    else if (g->htype == GRIDHTYPE_CURVILINEAR)
        (void) gxy_curv_xy2fij(g->gridnodes_xy, x, y, fi, fj);
#endif
    else
        enkf_quit("programming error");
    if (isnan(*fi + *fj))
        return STATUS_OUTSIDEGRID;

    /*
     * This procedure is used to map observations within a grid. The observation
     * will then be possibly collated (superobed) and stored in observations.nc.
     * In particular, the fractional grid coordinates of the observation will be
     * stored as NC_FLOAT. This onse in a lifetime causes a mismatch between i1,
     * i2, j1, j2 values below obtained from double precision coordinates and
     * those obtained from single precision ccordinates and then used for
     * getting the forecast observation values by interpolating the
     * corresponding model field. Therefore this bit:
     */
    *fi = (float) (*fi);
    *fj = (float) (*fj);

    /*
     * Note that this section should be consistent with similar sections in
     * interpolate2d() and interpolate3d().
     */
    i1 = floor(*fi);
    i2 = ceil(*fi);
    j1 = floor(*fj);
    j2 = ceil(*fj);

    grid_getsize(g, &ni, &nj, NULL);
    if (i1 == -1)
        i1 = (isperiodic_i) ? ni - 1 : i2;
    if (i2 == ni)
        i2 = (isperiodic_i) ? 0 : i1;
    if (j1 == -1)
        j1 = j2;
    if (j2 == nj)
        j2 = j1;

    if (g->numlevels[j1][i1] == 0 && g->numlevels[j1][i2] == 0 && g->numlevels[j2][i1] == 0 && g->numlevels[j2][i2] == 0) {
        *fi = NAN;
        *fj = NAN;
        return STATUS_LAND;
    }
    return STATUS_OK;
}

/**
 */
int grid_z2fk(grid* g, double fi, double fj, double z, double* fk)
{
    int isperiodic_i = grid_isperiodic_i(g);
    int ksurf = grid_getsurflayerid(g);
    int ni = -1, nj = -1;
    int i1, i2, j1, j2, k;

    if (isnan(fi + fj)) {
        *fk = NAN;
        return STATUS_OUTSIDEGRID;
    }

    if (g->vtype == GRIDVTYPE_Z)
        gz_z_z2fk(g, fi, fj, z, fk);
    else if (g->vtype == GRIDVTYPE_SIGMA)
        gz_sigma_z2fk(g, fi, fj, z, fk);
    else if (g->vtype == GRIDVTYPE_HYBRID)
        gz_hybrid_z2fk(g, fi, fj, z, fk);
    else
        enkf_quit("not implemented");

    if (isnan(*fk))
        return STATUS_OUTSIDEGRID;

    if (g->vtype != GRIDVTYPE_Z || g->depth == NULL)
        return STATUS_OK;

    /*
     * see the note on the similar bit in grid_xy2fij()
     */
    *fk = (float) (*fk);

    /*
     * a depth check for z-grid:
     */
    grid_getsize(g, &ni, &nj, NULL);
    i1 = floor(fi);
    i2 = ceil(fi);
    if (i1 == -1)
        i1 = (isperiodic_i) ? ni - 1 : i2;
    if (i2 == ni)
        i2 = (isperiodic_i) ? 0 : i1;
    j1 = floor(fj);
    j2 = ceil(fj);
    if (j1 == -1)
        j1 = j2;
    if (j2 == nj)
        j2 = j1;
    k = (ksurf == 0) ? floor(*fk) : ksurf - ceil(*fk);
    if (g->numlevels[j1][i1] <= k && g->numlevels[j1][i2] <= k && g->numlevels[j2][i1] <= k && g->numlevels[j2][i2] <= k) {
        *fk = NAN;
        return STATUS_LAND;
    } else if (g->numlevels[j1][i1] <= k || g->numlevels[j1][i2] <= k || g->numlevels[j2][i1] <= k || g->numlevels[j2][i2] <= k) {
        double v;

        v = interpolate2d(fi, fj, ni, nj, g->depth, g->numlevels, isperiodic_i);

        if (z > v)
            return STATUS_LAND;
    }

    return STATUS_OK;
}

/**
 */
void grid_fij2xy(grid* g, double fi, double fj, double* x, double* y)
{
    if (g->htype == GRIDHTYPE_LATLON)
        gs_fij2xy(g, fi, fj, x, y);
#if defined(ENKF_PREP) || defined(ENKF_CALC)
    else if (g->htype == GRIDHTYPE_CURVILINEAR)
        (void) gxy_curv_fij2xy(g->gridnodes_xy, fi, fj, x, y);
#endif
    else
        enkf_quit("programming error");
}

/**
 */
void grid_ij2xy(grid* g, int i, int j, double* x, double* y)
{
    if (g->htype == GRIDHTYPE_LATLON) {
        gxy_simple* gs = (gxy_simple*) g->gridnodes_xy;

        if (i < 0 || j < 0 || i >= gs->ni || j >= gs->nj) {
            *x = NAN;
            *y = NAN;
        } else {
            *x = gs->x[i];
            *y = gs->y[j];
        }
#if defined(ENKF_PREP) || defined(ENKF_CALC)
    } else if (g->htype == GRIDHTYPE_CURVILINEAR) {
        gxy_curv* gxy = g->gridnodes_xy;

        if (i < 0 || j < 0 || i >= gxy_curv_getni(gxy) || j >= gxy_curv_getnj(gxy)) {
            *x = NAN;
            *y = NAN;
        } else {
            *x = gxy_curv_getx(gxy)[j][i];
            *y = gxy_curv_gety(gxy)[j][i];
        }
#endif
    } else
        enkf_quit("programming error");
}

/**
 */
int grid_fk2z(grid* g, int i, int j, double fk, double* z)
{
    int ni = -1, nj = -1;

    grid_getsize(g, &ni, &nj, NULL);
    if (i < 0 || j < 0 || i >= ni || j >= nj) {
        *z = NAN;
        return STATUS_OUTSIDEGRID;
    }

    fk += 0.5;

    if (g->vtype == GRIDVTYPE_Z) {
        gz_z* gz = g->gridnodes_z;
        double* zc = gz->zc;
        int nt = gz->nk;

        if (fk <= 0.0)
            *z = zc[0];
        else if (fk >= nt)
            *z = zc[nt];
        else {
            int k = (int) floor(fk);
            double dk = fk - (double) k;

            if (dk < 0.5)
                *z = gz->zc[k] + dk * (gz->zt[k] - gz->zc[k]) / 0.5;
            else
                *z = gz->zt[k] + (dk - 0.5) * (gz->zc[k + 1] - gz->zt[k]) / 0.5;
        }
    } else if (g->vtype == GRIDVTYPE_SIGMA) {
        gz_sigma* gz = (gz_sigma*) g->gridnodes_z;
        double h;

        if (g->depth == NULL)
            enkf_quit("%s: DEPTHVARNAME must be entered for SIGMA grids to assimilate subsurface obs", g->name);
        h = g->depth[j][i];

        if (isnan(gz->fi_prev) || fabs((double) i - gz->fi_prev) > EPS_IJ || fabs((double) j - gz->fj_prev) > EPS_IJ) {
            int k;

            if (gz->hc != 0.0) {
                for (k = 0; k < gz->nk; ++k)
                    gz->zt[k] = h * (gz->hc * (1.0 - gz->st[k]) + h * gz->ct[k]) / (gz->hc + h);
                for (k = 0; k <= gz->nk; ++k)
                    gz->zc[k] = h * (gz->hc * (1.0 - gz->sc[k]) + h * gz->cc[k]) / (gz->hc + h);
            } else {
                for (k = 0; k < gz->nk; ++k)
                    gz->zt[k] = h * gz->ct[k];
                for (k = 0; k <= gz->nk; ++k)
                    gz->zc[k] = h * gz->cc[k];
            }
        }

        if (fk <= 0.0)
            *z = gz->zc[0];
        else if (fk >= gz->nk)
            *z = gz->zc[gz->nk];
        else {
            int k = (int) floor(fk);
            double dk = fk - (double) k;

            if (dk < 0.5)
                *z = gz->zc[k] + dk * (gz->zt[k] - gz->zc[k]) / 0.5;
            else
                *z = gz->zt[k] + (dk - 0.5) * (gz->zc[k + 1] - gz->zt[k]) / 0.5;
        }
    } else if (g->vtype == GRIDVTYPE_HYBRID) {
        gz_hybrid* gz = (gz_hybrid*) g->gridnodes_z;
        int nk = gz->nk;

        if (isnan(gz->fi_prev) || fabs((double) i - gz->fi_prev) > EPS_IJ || fabs((double) j - gz->fj_prev) > EPS_IJ) {
            double p1 = gz->p1[j][i];
            double p2 = gz->p2[j][i];
            int k;

            for (k = 0; k < nk; ++k)
                gz->pt[k] = gz->at[k] + gz->bt[k] * (p1 - p2);
            for (k = 0; k <= nk; ++k)
                gz->pc[k] = gz->ac[k] + gz->bc[k] * (p1 - p2);
            gz->fi_prev = (double) i;
            gz->fj_prev = (double) j;
        }

        if (fk <= 0.0)
            *z = gz->pc[0];
        else if (fk >= nk)
            *z = gz->pc[nk];
        else {
            int k = (int) floor(fk);
            double dk = fk - (double) k;

            if (dk < 0.5)
                *z = gz->pc[k] + dk * (gz->pt[k] - gz->pc[k]) / 0.5;
            else
                *z = gz->pt[k] + (dk - 0.5) * (gz->pc[k + 1] - gz->pt[k]) / 0.5;
        }
    } else
        enkf_quit("programming error");

    if (g->depth != NULL && *z > g->depth[j][i]) {
        *z = NAN;
        return STATUS_OUTSIDEGRID;
    }

    return STATUS_OK;

}

/**
 */
int grid_isperiodic_i(grid* g)
{
    if (g->htype == GRIDHTYPE_LATLON)
        return ((gxy_simple*) ((grid*) g)->gridnodes_xy)->periodic_i;

    return 0;
}

/**
 */
static void grids_addgrid(int* ngrid, void*** grids, void* g)
{
    if (*ngrid % GRID_INC == 0)
        (*grids) = realloc((*grids), (*ngrid + GRID_INC) * sizeof(void*));
    (*grids)[*ngrid] = g;
    (*ngrid)++;
}

/**
 */
void grids_create(char gprmfname[], int stride, int* ngrid, void*** grids)
{
    int n = 0;
    gridprm* gprm = NULL;
    int i;

    gridprm_create(gprmfname, &n, &gprm);
    assert(n > 0);

    for (i = 0; i < n; ++i) {
        grid* g = NULL;

        g = grid_create(&gprm[i], i);
        grids_addgrid(ngrid, grids, g);

        if (grid_getstride(g) == 0)
            grid_setstride(g, stride);
    }
    gridprm_destroy(n, gprm);
}

/**
 */
void grids_destroy(int ngrid, void** grids)
{
    int i;

    for (i = 0; i < ngrid; ++i)
        grid_destroy(grids[i]);
    free(grids);
}

/**
 */
char* grid_getdomainname(grid* g)
{
    return g->domainname;
}

#if defined(ENKF_CALC)
/**
 */
kdtree* grid_gettreeXYZ(grid* g)
{
    kdtree* tree;
    int ni, nj;
    size_t* ids;

    if (g->nodetreeXYZ != NULL)
        return g->nodetreeXYZ;

    tree = kd_create(3);
    grid_getsize(g, &ni, &nj, NULL);
    ids = malloc(ni * nj * sizeof(size_t));
    if (g->htype == GRIDHTYPE_LATLON) {
        gxy_simple* gxy = (gxy_simple*) g->gridnodes_xy;
        size_t ii, n;

        for (ii = 0, n = 0; ii < ni * nj; ++ii) {
            if (g->numlevels[0][ii] == 0)
                continue;
            ids[n] = ii;
            n++;
        }
        shuffle(n, ids);
        for (ii = 0; ii < n; ++ii) {
            int id = ids[ii];
            double ll[2], xyz[3];

            ll[0] = gxy->x[id % ni];
            ll[1] = gxy->y[id / ni];
            ll2xyz(ll, xyz);
            kd_insertnode(tree, xyz, ids[ii]);
        }
    } else if (g->htype == GRIDHTYPE_CURVILINEAR) {
        gxy_curv* gxy = (gxy_curv*) g->gridnodes_xy;
        double** x = gxy_curv_getx(gxy);
        double** y = gxy_curv_gety(gxy);
        size_t ii, n;

        for (ii = 0, n = 0; ii < ni * nj; ++ii) {
            if (g->numlevels[0][ii] == 0 || isnan(x[0][ii]))
                continue;
            ids[n] = ii;
            n++;
        }
        shuffle(n, ids);
        for (ii = 0; ii < n; ++ii) {
            double ll[2], xyz[3];

            ll[0] = x[0][ids[ii]];
            ll[1] = y[0][ids[ii]];
            ll2xyz(ll, xyz);
            kd_insertnode(tree, xyz, ids[ii]);
        }
    }
    free(ids);
    g->nodetreeXYZ = tree;

    return g->nodetreeXYZ;
}
#endif
