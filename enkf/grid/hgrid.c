/******************************************************************************
 *
 * File:        hgrid.c
 *
 * Created:     25/01/2022
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Purpose:     Code for `hgrid' object (generic horizontal grid).
 *
 * Description: Specific grid types are handled via field `gxy' of `hgrid'.
 *
 * Revisions:
 *
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <float.h>
#include "definitions.h"
#include "utils.h"
#include "ncw.h"
#include "ncutils.h"
#include "grid.h"
#include "gridprm.h"
#include "gxy_rect.h"
#include "gxy_curv.h"
#include "gxy_unstr.h"
#include "hgrid.h"

/*
 * for DIAG the grid stuff is the same as for UPDATE
 */
#if defined(ENKF_DIAG)
#define ENKF_UPDATE
#endif

struct gxy_1d {
    int n;
};

struct gxy_2d {
    int ni;
    int nj;
};

/**
 */
static gxy_2d* gxy_2d_create(int ni, int nj)
{
    gxy_2d* gxy = malloc(sizeof(gxy_2d));

    gxy->ni = ni;
    gxy->nj = nj;

    return gxy;
}

/**
 */
#if defined(ENKF_UPDATE)
static gxy_1d* gxy_1d_create(int n)
{
    gxy_1d* gxy = malloc(sizeof(gxy_1d));

    gxy->n = n;

    return gxy;
}
#endif

#if !defined(ENKF_UPDATE)
static void hgrid_setlonbase(hgrid* hg)
{
    double xmin = DBL_MAX;
    double xmax = -DBL_MAX;

    if (hg->type == GRIDHTYPE_RECTANGULAR) {
        gxy_rect* gxy = hg->gxy;
        double* x = gxy_rect_getx(gxy);
        int ni = gxy_rect_getni(gxy);

        if (gxy_rect_getperiodic_i(gxy) == 2)   /* periodic non-closed grid */
            ni++;

        if (xmin > x[0])
            xmin = x[0];
        if (xmin > x[ni - 1])
            xmin = x[ni - 1];
        if (xmax < x[0])
            xmax = x[0];
        if (xmax < x[ni - 1])
            xmax = x[ni - 1];
    } else if (hg->type == GRIDHTYPE_CURVILINEAR) {
        gxy_curv* gxy = (gxy_curv*) hg->gxy;

        if (!gxy_curv_isgeographic(gxy)) {
            double* minmax = kd_getminmax(gxy_curv_gettree(gxy));

            xmin = minmax[0];
            xmax = minmax[2];
        } else {
            int i;
            int nij = gxy_curv_getni(gxy) * gxy_curv_getnj(gxy);
            double* x = gxy_curv_getx(gxy)[0];

            for (i = 0; i < nij; ++i) {
                if (x[i] < xmin)
                    xmin = x[i];
                if (x[i] > xmax)
                    xmax = x[i];
            }
        }
    } else if (hg->type == GRIDHTYPE_UNSTRUCTURED)
        triangulation_getminmax(gxy_unstr_gettriangulation(hg->gxy), &xmin, &xmax, NULL, NULL);

    if (xmin >= 0.0 && xmax <= 360.0)
        hg->lonbase = 0.0;
    else if (xmin >= -180.0 && xmax <= 180.0)
        hg->lonbase = -180.0;
    else
        hg->lonbase = xmin;
}
#endif

hgrid* hgrid_create(void* p, void* g)
{
    gridprm* prm = p;
    hgrid* hg;
    int ncid;
    int varid_x, varid_y;
    int ndims_x, ndims_y;
    size_t ni = 0, nj = 0;

    hg = calloc(1, sizeof(hgrid));
    hg->type = gridprm_gethtype(prm);
    hg->parent = g;
    hg->lonbase = NAN;

    if (hg->type == GRIDHTYPE_NONE) {
        hg->ni = 1;
        hg->nj = 1;
        hg->gxy = NULL;

        return hg;
    }

    ncw_open(prm->gdatafname, NC_NOWRITE, &ncid);

    ncw_inq_varid(ncid, prm->xvarname, &varid_x);
    ncw_inq_varndims(ncid, varid_x, &ndims_x);
    ncw_inq_varid(ncid, prm->yvarname, &varid_y);
    ncw_inq_varndims(ncid, varid_y, &ndims_y);

    if (ndims_x == 1 && ndims_y == 1) {
        int dimid_x, dimid_y;

        ncw_inq_vardimid(ncid, varid_x, &dimid_x);
        ncw_inq_vardimid(ncid, varid_y, &dimid_y);

        if (dimid_x != dimid_y) {
            double* x;
            double* y;

            if (hg->type == GRIDHTYPE_UNDEFINED)
                hg->type = GRIDHTYPE_RECTANGULAR;
            else if (hg->type != GRIDHTYPE_RECTANGULAR)
                enkf_quit("%s: grid \"%s\": 1-dimensional coordinates of different size are only possible for HTYPE = RECT", prm->prmfname, prm->name);

            ncw_inq_vardims(ncid, varid_x, 1, NULL, &ni);
            ncw_inq_vardims(ncid, varid_y, 1, NULL, &nj);

            x = malloc(ni * sizeof(double));
            y = malloc(nj * sizeof(double));

            ncu_readvardouble(ncid, varid_x, ni, x);
            ncu_readvardouble(ncid, varid_y, nj, y);

            hg->gxy = gxy_rect_create(g, ni, nj, x, y);
            hg->periodic_i = gxy_rect_getperiodic_i(hg->gxy);
        } else {
            if (hg->type == GRIDHTYPE_UNDEFINED)
                hg->type = GRIDHTYPE_UNSTRUCTURED;
            else if (hg->type != GRIDHTYPE_UNSTRUCTURED)
                enkf_quit("%s: grid \"%s\": 1-dimensional coordinates of same size are only possible for HTYPE = UNSTRUCT", prm->prmfname, prm->name);
#if defined(ENKF_UPDATE)
            ncw_inq_vardims(ncid, varid_x, 1, NULL, &ni);

            hg->type = GRIDHTYPE_1D;
            hg->gxy = gxy_1d_create(ni);
#else
            triangulation* d = triangulation_read(prm->gdatafname, prm->xvarname, prm->yvarname, prm->trivarname, prm->neivarname);

            hg->gxy = gxy_unstr_create(g, d);
            ni = d->npoints;
#endif
        }
    } else if (ndims_x == 2 && ndims_y == 2) {
        if (hg->type == GRIDHTYPE_UNDEFINED) {
#if defined(ENKF_UPDATE)
            hg->type = GRIDHTYPE_2D;
#else
            hg->type = GRIDHTYPE_CURVILINEAR;
#endif
        } else if (hg->type != GRIDHTYPE_RECTANGULAR && hg->type != GRIDHTYPE_CURVILINEAR)
            enkf_quit("%s: grid \"%s\": 2-dimensional coordinates are only possible for HTYPE = CURV or HTYPE = RECT", prm->prmfname, prm->name);

        if (hg->type == GRIDHTYPE_2D) {
            size_t dimlen[2];

            ncw_inq_vardims(ncid, varid_x, 2, NULL, dimlen);
            ncw_check_vardims(ncid, varid_y, 2, dimlen);
            nj = dimlen[0];
            ni = dimlen[1];

            hg->type = GRIDHTYPE_2D;
            hg->gxy = gxy_2d_create(ni, nj);
        } else {
            double** x;
            double** y;
            double* xx;
            double* yy;
            size_t dimlen[2];

            ncw_inq_vardims(ncid, varid_x, 2, NULL, dimlen);
            ncw_check_vardims(ncid, varid_y, 2, dimlen);
            nj = dimlen[0];
            ni = dimlen[1];

            x = alloc2d(nj, ni, sizeof(double));
            y = alloc2d(nj, ni, sizeof(double));

            ncu_readvardouble(ncid, varid_x, ni * nj, x[0]);
            ncu_readvardouble(ncid, varid_y, ni * nj, y[0]);

            if (hg->type == GRIDHTYPE_RECTANGULAR) {
                int i, j;

                for (j = 1; j < nj; ++j)
                    for (i = 0; i < ni; ++i)
                        if (x[j][i] != x[0][i])
                            enkf_quit("%s: grid \"%s\": for HTYPE = RECT with 2-dimensional coordinates X coordinates must be identical for each grid row", prm->prmfname, prm->name);
                for (j = 0; j < nj; ++j)
                    for (i = 1; i < ni; ++i)
                        if (y[j][0] != y[j][i])
                            enkf_quit("%s: grid \"%s\": for HTYPE = RECT with 2-dimensional coordinates Y coordinates must be identical for each grid column", prm->prmfname, prm->name);

                xx = malloc(ni * sizeof(double));
                yy = malloc(nj * sizeof(double));
                for (i = 0; i < ni; ++i)
                    xx[i] = x[0][i];
                for (j = 0; j < nj; ++j)
                    yy[j] = y[j][0];
                free(x);
                free(y);
                hg->gxy = gxy_rect_create(g, ni, nj, xx, yy);
                hg->periodic_i = gxy_rect_getperiodic_i(hg->gxy);
            } else {
#if defined(ENKF_PREP) || defined(ENKF_CALC)
                hg->gxy = gxy_curv_create(g, ni, nj, x, y, grid_getnumlevels(g), prm->geographic);
#else
                enkf_quit("programming error");
#endif
            }
        }
    }
    ncw_close(ncid);

    hg->ni = ni;
    hg->nj = nj;

#if defined(ENKF_PREP) || defined(ENKF_CALC)
    hgrid_setlonbase(hg);
#endif

    return hg;
}

/**
 */
void hgrid_destroy(hgrid* hg)
{
    if (hg->type == GRIDHTYPE_RECTANGULAR)
        gxy_rect_destroy(hg->gxy);
#if defined(ENKF_PREP) || defined(ENKF_CALC)
    else if (hg->type == GRIDHTYPE_CURVILINEAR)
        gxy_curv_destroy(hg->gxy);
    else if (hg->type == GRIDHTYPE_UNSTRUCTURED)
        gxy_unstr_destroy(hg->gxy);
#else
    else if (hg->type == GRIDHTYPE_CURVILINEAR || hg->type == GRIDHTYPE_UNSTRUCTURED)
        enkf_quit("programming error");
#endif
#if defined(ENKF_CALC)
    if (hg->nodetreeXYZ != NULL)
        kd_destroy(hg->nodetreeXYZ);
#endif
    free(hg);
}

/**
 */
#if defined(ENKF_PREP) || defined(ENKF_CALC)
int hgrid_destroynodetree(hgrid* hg)
{
    if (hg->type == GRIDHTYPE_CURVILINEAR) {
        gxy_curv_destroykdtree(hg->gxy);
        return 1;
    }
    return 0;
}
#endif

/**
 */
void hgrid_describe(hgrid* hg, char* offset)
{
    switch (hg->type) {
    case GRIDHTYPE_1D:
        enkf_printf("%s  h type = 1D (no coords)\n", offset);
        break;
    case GRIDHTYPE_2D:
        enkf_printf("%s  h type = 2D (no coords)\n", offset);
        break;
    case GRIDHTYPE_RECTANGULAR:
        enkf_printf("%s  h type = RECTANGULAR\n", offset);
        enkf_printf("%s  periodic by X = %s\n", offset, (gxy_rect_getperiodic_i(hg->gxy)) ? "yes" : "no");
        break;
    case GRIDHTYPE_CURVILINEAR:
        enkf_printf("%s  h type = CURVILINEAR\n", offset);
#if defined(ENKF_PREP) || defined(ENKF_CALC)
        enkf_printf("%s  geographic = %s\n", offset, (gxy_curv_isgeographic(hg->gxy)) ? "yes" : "no");
#endif
        break;
    case GRIDHTYPE_UNSTRUCTURED:
        enkf_printf("%s  h type = UNSTRUCTURED\n", offset);
        break;
    default:
        enkf_printf("%s  h type = NONE\n", offset);
    }
    if (!isnan(hg->lonbase))
        enkf_printf("%s  longitude range = [%.3f, %.3f]\n", offset, hg->lonbase, hg->lonbase + 360.0);
    else
        enkf_printf("%s  longitude range = any\n", offset);
#if defined(ENKF_PREP) || defined(ENKF_CALC)
    if (hg->type == GRIDHTYPE_CURVILINEAR) {
        char offset2[SHORTSTRLEN];

        sprintf(offset2, "%s%s", offset, "  ");
        kd_printinfo(gxy_curv_gettree(hg->gxy), offset2);
    }
#endif
}

/**
 */
#if defined(ENKF_PREP) || defined(ENKF_CALC)
int hgrid_xy2fij(hgrid* hg, void* mask, double x, double y, double* fij)
{
    if (!isnan(hg->lonbase)) {
        if (x < hg->lonbase)
            x += 360.0;
        else if (x >= hg->lonbase + 360.0)
            x -= 360.0;
    }

    if (hg->type == GRIDHTYPE_RECTANGULAR || hg->type == GRIDHTYPE_CURVILINEAR) {
        int i1, i2, j1, j2;

        /*
         * fij[2] is redundant for structured grids. isnan(fij[2]) is used as
         * a test for structured/unstrctured grids in cmp_obs() in
         * enkf_prep.c:main().
         */
        fij[2] = NAN;

        if (hg->type == GRIDHTYPE_RECTANGULAR) {
            gxy_rect_xy2fij(hg->gxy, x, y, fij);
        } else if (hg->type == GRIDHTYPE_CURVILINEAR)
            (void) gxy_curv_xy2fij(hg->gxy, x, y, fij);
        if (isnan(fij[0] + fij[1]))
            return STATUS_OUTSIDEGRID;

        /*
         * Checking that the observation is in watered point.
         *
         * Note that this section should be consistent with similar sections in
         * interpolate2d() and interpolate3d().
         */
        i1 = floor(fij[0]);
        i2 = ceil(fij[0]);
        j1 = floor(fij[1]);
        j2 = ceil(fij[1]);

        if (i1 == -1)
            i1 = (hg->periodic_i) ? hg->ni - 1 : i2;
        if (i2 == hg->ni)
            i2 = (hg->periodic_i) ? 0 : i1;
        if (j1 == -1)
            j1 = j2;
        if (j2 == hg->nj)
            j2 = j1;

        {
            int** numlevels = mask;

            if (numlevels[j1][i1] == 0 && numlevels[j1][i2] == 0 && numlevels[j2][i1] == 0 && numlevels[j2][i2] == 0) {
                fij[0] = NAN;
                fij[1] = NAN;
                return STATUS_LAND;
            }
        }
    } else if (hg->type == GRIDHTYPE_UNSTRUCTURED) {
        int* numlevels = mask;

        if (!gxy_unstr_xy2fij(hg->gxy, x, y, fij))
            return STATUS_OUTSIDEGRID;

        if (numlevels[(int) fij[0]] == 0 && numlevels[(int) fij[1]] == 0 && numlevels[(int) fij[2]] == 0) {
            fij[0] = NAN;
            fij[1] = NAN;
            fij[2] = NAN;
            return STATUS_LAND;
        }
    } else
        enkf_quit("programming error");

    return STATUS_OK;
}
#endif

/**
 */
#if defined(ENKF_PREP) || defined(ENKF_CALC)
void hgrid_ij2xy(hgrid* hg, int* ij, double* x, double* y)
{
    if (hg->type == GRIDHTYPE_NONE) {
        *x = NAN;
        *y = NAN;
    } else if (hg->type == GRIDHTYPE_RECTANGULAR) {
        int i = ij[0];
        int j = ij[1];

        if (i < 0 || j < 0 || i >= hg->ni || j >= hg->nj) {
            *x = NAN;
            *y = NAN;
        } else {
            *x = gxy_rect_getx(hg->gxy)[i];
            *y = gxy_rect_gety(hg->gxy)[j];
        }
    } else if (hg->type == GRIDHTYPE_CURVILINEAR) {
        int i = ij[0];
        int j = ij[1];

        if (i < 0 || j < 0 || i >= hg->ni || j >= hg->nj) {
            *x = NAN;
            *y = NAN;
        } else {
            *x = gxy_curv_getx(hg->gxy)[j][i];
            *y = gxy_curv_gety(hg->gxy)[j][i];
        }
    } else if (hg->type == GRIDHTYPE_UNSTRUCTURED) {
        point* p = gxy_unstr_getpoint(hg->gxy, ij[0]);

        if (p == NULL) {
            *x = NAN;
            *y = NAN;
        } else {
            *x = p->x;
            *y = p->y;
        }
    } else
        enkf_quit("programming error");
}
#endif

/**
 */
#if defined(ENKF_CALC)
kdtree* hgrid_gettreeXYZ(hgrid* hg, int createifnull)
{
    char name[MAXSTRLEN];
    kdtree* tree;

    if (!createifnull || hg->nodetreeXYZ != NULL)
        return hg->nodetreeXYZ;

    snprintf(name, MAXSTRLEN - 4, "%s_XYZ", grid_getname(hg->parent));
    tree = kd_create(name, 3);
    if (hg->type == GRIDHTYPE_RECTANGULAR) {
        double* x = gxy_rect_getx(hg->gxy);
        double* y = gxy_rect_gety(hg->gxy);
        int** numlevels = grid_getnumlevels(hg->parent);
        size_t* ids;
        size_t ii, n;

        ids = malloc(hg->ni * hg->nj * sizeof(size_t));

        for (ii = 0, n = 0; ii < hg->ni * hg->nj; ++ii) {
            if (numlevels[0][ii] == 0)
                continue;
            ids[n] = ii;
            n++;
        }
        shuffle(n, ids);
        for (ii = 0; ii < n; ++ii) {
            int id = ids[ii];
            double ll[2], xyz[3];

            ll[0] = x[id % hg->ni];
            ll[1] = y[id / hg->ni];
            ll2xyz(ll, xyz);
            kd_insertnode(tree, xyz, ids[ii]);
        }
        free(ids);
    } else if (hg->type == GRIDHTYPE_CURVILINEAR) {
        double** x = gxy_curv_getx(hg->gxy);
        double** y = gxy_curv_gety(hg->gxy);
        int** numlevels = grid_getnumlevels(hg->parent);
        size_t* ids;
        size_t ii, n;

        ids = malloc(hg->ni * hg->nj * sizeof(size_t));

        for (ii = 0, n = 0; ii < hg->ni * hg->nj; ++ii) {
            if (numlevels[0][ii] == 0 || isnan(x[0][ii]))
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
        free(ids);
    } else if (hg->type == GRIDHTYPE_UNSTRUCTURED) {
        int* numlevels = grid_getnumlevels(hg->parent);
        size_t* ids;
        size_t ii, n;

        ids = malloc(hg->ni * sizeof(size_t));

        for (ii = 0, n = 0; ii < hg->ni; ++ii) {
            if (numlevels[ii] == 0)
                continue;
            ids[n] = ii;
            n++;
        }
        shuffle(n, ids);
        for (ii = 0; ii < n; ++ii) {
            point* p = gxy_unstr_getpoint(hg->gxy, ids[ii]);
            double ll[2] = { p->x, p->y };
            double xyz[3];

            ll2xyz(ll, xyz);
            kd_insertnode(tree, xyz, ids[ii]);
        }
        free(ids);
    } else
        enkf_quit("programming error");
    kd_finalise(tree);
    hg->nodetreeXYZ = tree;

    return hg->nodetreeXYZ;
}
#endif

/**
 */
#if defined(ENKF_CALC)
void hgrid_destroytreeXYZ(hgrid* hg)
{
    kd_destroy(hg->nodetreeXYZ);
    hg->nodetreeXYZ = NULL;
}
#endif
