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
 * Revisions:
 *
 *****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <values.h>
#include <string.h>
#include "ncw.h"
#include "definitions.h"
#include "utils.h"
#include "grid.h"
#include "gridprm.h"
#if !defined(NO_GRIDUTILS)
#include <gridnodes.h>
#endif
#include "gridmap.h"
#include "nan.h"

#define EPSLON 1.0e-3

typedef struct {
    int nx;
    int ny;
    int periodic_x;
    int periodic_y;
    int regular;

    double* x;
    double* y;
    double* xc;
    double* yc;
} gnxy_simple;

#if !defined(NO_GRIDUTILS)
typedef struct {
    gridnodes* gn;
    gridmap* gm;
} gnxy_curv;
#else
#define NT_NONE 0
#endif

typedef struct {
    int nz;
    double* zt;
    double* zc;
} gnz;

typedef void (*grid_xy2fij_fn) (void* grid, double x, double y, double* fi, double* fj);
typedef void (*grid_z2fk_fn) (void* grid, double fi, double fj, double z, double* fk);
typedef void (*grid_fij2xy_fn) (void* grid, double fi, double fj, double* x, double* y);

struct grid {
    char* name;
    int id;
    int htype;                  /* horizontal type */
    int vtype;                  /* vertical type */

    grid_xy2fij_fn xy2fij_fn;
    grid_z2fk_fn z2fk_fn;
    grid_fij2xy_fn fij2xy_fn;
    grid_tocartesian_fn tocartesian_fn;

    void* gridnodes_xy;         /* (the structure is defined by `htype') */
    int lontype;                /* (range: any, [-180,180) or [0,360)) */

    gnz* gridnodes_z;

    /*
     * `numlevels' can hold either the number of levels (z-model) or the
     * land mask (sigma-model)
     */
    int** numlevels;
    float** depth;
};

/**
 */
static gnxy_simple* gnxy_simple_create(int nx, int ny, double* x, double* y, int periodic_x, int periodic_y, int regular)
{
    gnxy_simple* nodes = malloc(sizeof(gnxy_simple));
    int i;

    nodes->nx = nx;
    nodes->ny = ny;
    nodes->x = x;
    nodes->y = y;
    nodes->periodic_x = periodic_x;
    nodes->periodic_y = periodic_y;
    nodes->regular = regular;

    if (regular) {
        nodes->xc = NULL;
        nodes->yc = NULL;
    } else {
        nodes->xc = malloc((nx + 1) * sizeof(double));
        nodes->xc[0] = x[0] * 1.5 - x[1] * 0.5;
        for (i = 1; i < nx + 1; ++i)
            nodes->xc[i] = 2 * x[i - 1] - nodes->xc[i - 1];

        nodes->yc = malloc((ny + 1) * sizeof(double));
        nodes->yc[0] = y[0] * 1.5 - y[1] * 0.5;
        for (i = 1; i < ny + 1; ++i)
            nodes->yc[i] = 2 * y[i - 1] - nodes->yc[i - 1];
    }

    return nodes;
}

/**
 */
void gnxy_simple_destroy(gnxy_simple* nodes)
{
    free(nodes->x);
    free(nodes->y);
    if (!nodes->regular) {
        free(nodes->xc);
        free(nodes->yc);
    }
    free(nodes);
}

#if !defined(NO_GRIDUTILS)
/**
 */
static gnxy_curv* gnxy_curv_create(int nodetype, int nx, int ny, double** x, double** y)
{
    gnxy_curv* nodes = malloc(sizeof(gnxy_curv));

    if (nodetype == NT_CEN) {
        gridnodes* gn_new;

        nodes->gn = gridnodes_create2(nx, ny, NT_CEN, x, y);
        gn_new = gridnodes_transform(nodes->gn, NT_COR);
        gridnodes_destroy(nodes->gn);
        nodes->gn = gn_new;
    } else if (nodetype == NT_COR)
        nodes->gn = gridnodes_create2(nx, ny, NT_COR, x, y);
    else
        enkf_quit("unknown node type for horizontal curvilinear grid");
    gridnodes_validate(nodes->gn);
    nodes->gm = gridmap_build(gridnodes_getnce1(nodes->gn), gridnodes_getnce2(nodes->gn), gridnodes_getx(nodes->gn), gridnodes_gety(nodes->gn));

    return nodes;
}

/**
 */
void gnxy_curv_destroy(gnxy_curv* nodes)
{
    gridnodes_destroy(nodes->gn);
    gridmap_destroy(nodes->gm);
    free(nodes);
}
#endif

#define EPSZ 0.001

/**
 */
static gnz* gnz_create(int nz, double* z)
{
    gnz* nodes = malloc(sizeof(gnz));
    int i;

    /*
     * reverse z if it is negative
     */
    for (i = 0; i < nz; ++i)
        if (z[i] < -EPSZ)
            break;
    if (i < nz)
        for (i = 0; i < nz; ++i)
            z[i] = -z[i];
    for (i = 0; i < nz; ++i)
        if (z[i] < -EPSZ)
            enkf_quit("layer centre coordinates should be either positive or negative only");

    nodes->zt = z;
    nodes->nz = nz;

    /*
     * this code is supposed to work both for z and sigma grids
     */
    nodes->zc = malloc((nz + 1) * sizeof(double));
    if (z[nz - 1] >= z[0]) {
        /*
         * layer 0 at surface
         */
        nodes->zc[0] = 0.0;
        for (i = 1; i <= nz; ++i) {
            nodes->zc[i] = 2.0 * z[i - 1] - nodes->zc[i - 1];
            /*
             * layer boundary should be above the next layer centre
             */
            if (i < nz)
                assert(nodes->zc[i] < z[i]);
        }
    } else {
        /*
         * layer 0 the deepest
         */
        nodes->zc[nz] = 0.0;
        for (i = nz - 1; i >= 0; --i) {
            nodes->zc[i] = 2.0 * z[i] - nodes->zc[i + 1];
            /*
             * layer boundary should be above the next layer centre
             */
            if (i > 0)
                assert(nodes->zc[i] < z[i - 1]);
        }
    }

    return nodes;
}

/**
 */
static void gnz_destroy(gnz* nodes)
{
    free(nodes->zt);
    free(nodes->zc);
    free(nodes);
}

/**
 */
static double x2fi_reg(int n, double* v, double x, int periodic)
{
    double fi;

    if (n < 2)
        return NaN;

    if (periodic) {
        if (x < 0.0)
            x = x + 360.0;
        else if (x >= 360.0)
            x = x - 360.0;
    }

    fi = (x - v[0]) / (v[n - 1] - v[0]) * (double) (n - 1);

    if (fi < -0.5 || fi > (double) n - 0.5)
        return NaN;

    if (!periodic) {
        if (fi < 0.0)
            fi = 0.0;
        else if (fi > (double) (n - 1))
            fi = (double) (n - 1);
    }

    return fi;
}

/**
 */
static void g1_xy2fij(void* p, double x, double y, double* fi, double* fj)
{
    gnxy_simple* nodes = (gnxy_simple*) ((grid*) p)->gridnodes_xy;

    *fi = x2fi_reg(nodes->nx, nodes->x, x, nodes->periodic_x);
    *fj = x2fi_reg(nodes->ny, nodes->y, y, nodes->periodic_y);
}

/**
 */
static double fi2x(int n, double* v, double fi, int periodic)
{
    double ifrac;
    int i;

    if (n < 2)
        return NaN;

    if (fi < -1.0 || fi > (double) n)
        return NaN;

    ifrac = fi - floor(fi);

    if (fi < 0.0 || fi > (double) (n - 1)) {
        double v1, v2, x;

        if (!periodic)
            return NaN;

        v1 = v[n - 1];
        v2 = v[0];
        if (v1 - v2 > 180.0)
            v2 += 360.0;
        x = v1 + ifrac * (v2 - v1);

        return (x < 360.0) ? x : x - 360.0;
    }

    i = (int) fi;
    if (ifrac == 0.0)
        return v[i];
    else
        return v[i] + ifrac * (v[i + 1] - v[i]);
}

/**
 */
static void g12_fij2xy(void* p, double fi, double fj, double* x, double* y)
{
    gnxy_simple* nodes = (gnxy_simple*) ((grid*) p)->gridnodes_xy;

    *x = fi2x(nodes->nx, nodes->x, fi, nodes->periodic_x);
    *y = fi2x(nodes->ny, nodes->y, fj, nodes->periodic_y);
}

/** Gets fractional index of a coordinate for a 1D irregular grid.
 * @param n Number of grid nodes (vertical layers)
 * @param v Coordinates of the nodes (layer centres)
 * @param vb Coordinates of the cell/layer boundaries [n + 1]
 * @param x Input coordinate
 * @param periodic Flag for grid periodicity
 * @param lontype Range for longitude: 0 - any, 1 - [-180, 180), 2 - [0, 360)
 * @return Fractional index for `x'
 */
static double x2fi_irreg(int n, double v[], double vb[], double x, int periodic, int lontype)
{
    int ascending, i1, i2, imid;

    if (n < 2)
        return NaN;

    if (lontype == LONTYPE_180) {
        if (x < -180)
            x = x + 360.0;
        else if (x >= 180.0)
            x = x - 360.0;
    } else if (lontype == LONTYPE_360) {
        if (x < 0.0)
            x = x + 360.0;
        else if (x >= 360.0)
            x = x - 360.0;
    }

    ascending = (v[n - 1] > v[0]) ? 1 : 0;

    if ((ascending && x < v[0]) || (!ascending && x > v[0])) {
        double fi = (x - v[0]) / (v[1] - v[0]);

        if (periodic) {
            assert(fi >= -0.5);
            return fi;
        }
        return (fi >= -0.5) ? 0.0 : NaN;
    }
    if ((ascending && x > v[n - 1]) || (!ascending && x < v[n - 1])) {
        double fi = (double) (n - 1) + (x - v[n - 1]) / (v[n - 1] - v[n - 2]);

        if (periodic) {
            assert(fi < (double) n - 0.5);
            return fi;
        }
        return (fi <= (double) n - 0.5) ? (double) (n - 1) : NaN;
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

/**
 */
static void g2_xy2fij(void* p, double x, double y, double* fi, double* fj)
{
    gnxy_simple* nodes = (gnxy_simple*) ((grid*) p)->gridnodes_xy;
    int lontype = ((grid*) p)->lontype;

    *fi = x2fi_irreg(nodes->nx, nodes->x, nodes->xc, x, nodes->periodic_x, lontype);
    *fj = x2fi_irreg(nodes->ny, nodes->y, nodes->yc, y, nodes->periodic_y, LONTYPE_NONE);
}

#if !defined(NO_GRIDUTILS)
/**
 */
static void gc_xy2fij(void* p, double x, double y, double* fi, double* fj)
{
    gnxy_curv* nodes = (gnxy_curv*) ((grid*) p)->gridnodes_xy;

    if (gridmap_xy2fij(nodes->gm, x, y, fi, fj) == 1)
        return;                 /* success */

    *fi = NaN;
    *fj = NaN;
}

/**
 */
static void gc_fij2xy(void* p, double fi, double fj, double* x, double* y)
{
    gnxy_curv* nodes = (gnxy_curv*) ((grid*) p)->gridnodes_xy;

    gridmap_fij2xy(nodes->gm, fi, fj, x, y);
}
#endif

/**
 */
static double z2fk_basic(int n, double* zt, double* zc, double z)
{
    int ascending, i1, i2, imid;

    ascending = (zt[n - 1] > zt[0]) ? 1 : 0;

    if (ascending) {
        if (z < zc[0])
            return 0.0;
        if (z > zc[n])
            return NaN;
    } else {
        if (z < zc[n])
            return 0.0;
        if (z > zc[0])
            return NaN;
    }

    i1 = 0;
    i2 = n - 1;
    if (ascending) {
        while (1) {
            imid = (i1 + i2) / 2;
            if (imid == i1)
                break;
            if (z > zc[imid])
                i1 = imid;
            else
                i2 = imid;
        }
    } else {
        while (1) {
            imid = (i1 + i2) / 2;
            if (imid == i1)
                break;
            if (z > zc[imid])
                i2 = imid;
            else
                i1 = imid;
        }
    }

    if (z < zc[i1 + 1])
        return (double) i1 + (z - zt[i1]) / (zc[i1 + 1] - zc[i1]);
    else
        return (double) i1 + 0.5 + (z - zc[i1 + 1]) / (zc[i1 + 2] - zc[i1 + 1]);
}

/**
 */
static void z2fk(void* p, double fi, double fj, double z, double* fk)
{
    grid* g = (grid*) p;
    gnz* nodes = g->gridnodes_z;

    /*
     * for sigma coordinates convert `z' to sigma
     */
    if (g->vtype == GRIDVTYPE_SIGMA) {
        int ni = 0, nj = 0;
        double depth;

        if (g->htype == GRIDHTYPE_LATLON_REGULAR || g->htype == GRIDHTYPE_LATLON_IRREGULAR) {
            gnxy_simple* nodes = (gnxy_simple*) g->gridnodes_xy;

            ni = nodes->nx;
            nj = nodes->ny;
#if !defined(NO_GRIDUTILS)
        } else if (g->htype == GRIDHTYPE_CURVILINEAR) {
            gnxy_curv* nodes = (gnxy_curv*) g->gridnodes_xy;

            ni = gridnodes_getnx(nodes->gn);
            nj = gridnodes_getny(nodes->gn);
#endif
        } else
            enkf_quit("programming error");

        depth = (double) interpolate2d(fi, fj, ni, nj, g->depth, g->numlevels, grid_isperiodic_x(g), grid_isperiodic_y(g));
        z /= depth;
    }

    *fk = z2fk_basic(nodes->nz, nodes->zt, nodes->zc, z);
}

/**
 */
static void grid_setlontype(grid* g)
{
    double xmin = DBL_MAX;
    double xmax = -DBL_MAX;

    if (g->htype == GRIDHTYPE_LATLON_REGULAR || g->htype == GRIDHTYPE_LATLON_IRREGULAR) {
        double* x = ((gnxy_simple*) g->gridnodes_xy)->x;
        int nx = ((gnxy_simple*) g->gridnodes_xy)->nx;

        if (xmin < x[0])
            xmin = x[0];
        if (xmin < x[nx - 1])
            xmin = x[nx - 1];
        if (xmax > x[0])
            xmax = x[0];
        if (xmax > x[nx - 1])
            xmax = x[nx - 1];
#if !defined(NO_GRIDUTILS)
    } else if (g->htype == GRIDHTYPE_CURVILINEAR) {
        double** x = gridnodes_getx(((gnxy_curv*) g->gridnodes_xy)->gn);
        int nx = gridnodes_getnce1(((gnxy_curv*) g->gridnodes_xy)->gn);
        int ny = gridnodes_getnce2(((gnxy_curv*) g->gridnodes_xy)->gn);
        int i, j;

        for (j = 0; j < ny; ++j) {
            for (i = 0; i < nx; ++i) {
                if (xmin < x[j][i])
                    xmin = x[j][i];
                if (xmax > x[j][i])
                    xmax = x[j][i];
            }
        }
#endif
    }
    if (xmin < 0.0 && xmax <= 180.0)
        g->lontype = LONTYPE_180;
    else if (xmin >= 0 && xmax <= 360.0)
        g->lontype = LONTYPE_360;
}

/**
 */
static void grid_setcoords(grid* g, int htype, int hnodetype, int periodic_x, int periodic_y, int nx, int ny, int nz, void* x, void* y, double* z)
{
    g->htype = htype;
    if (htype == GRIDHTYPE_LATLON_REGULAR) {
        g->xy2fij_fn = g1_xy2fij;
        g->fij2xy_fn = g12_fij2xy;
        g->gridnodes_xy = gnxy_simple_create(nx, ny, x, y, periodic_x, periodic_y, 1);
    } else if (htype == GRIDHTYPE_LATLON_IRREGULAR) {
        g->xy2fij_fn = g2_xy2fij;
        g->fij2xy_fn = g12_fij2xy;
        g->gridnodes_xy = gnxy_simple_create(nx, ny, x, y, periodic_x, periodic_y, 0);
#if !defined(NO_GRIDUTILS)
    } else if (htype == GRIDHTYPE_CURVILINEAR) {
        g->xy2fij_fn = gc_xy2fij;
        g->fij2xy_fn = gc_fij2xy;
        g->gridnodes_xy = gnxy_curv_create(hnodetype, nx, ny, x, y);
#if defined(GRIDNODES_WRITE)
        {
            char fname[MAXSTRLEN];

            snprintf(fname, MAXSTRLEN, "gridnodes-%d.txt", grid_getid(g));
            if (!file_exists(fname))
                gridnodes_write(((gnxy_curv*) g->gridnodes_xy)->gn, fname, CT_XY);
        }
#endif
#endif
    } else
        enkf_quit("programming error");

    grid_setlontype(g);
    g->z2fk_fn = z2fk;

    g->gridnodes_z = gnz_create(nz, z);
}

/**
 */
grid* grid_create(void* p, int id)
{
    gridprm* prm = (gridprm*) p;
    grid* g = calloc(1, sizeof(grid));
    char* fname = prm->fname;
    int ncid;
    int dimid_x, dimid_y, dimid_z;
    int varid_x, varid_y, varid_z;
    int ndims_x, ndims_y, ndims_z;
    size_t nx, ny, nz;
    int varid_depth, varid_numlevels;

    g->name = strdup(prm->name);
    g->id = id;
    g->vtype = gridprm_getvtype(prm);

    ncw_open(fname, NC_NOWRITE, &ncid);
    ncw_inq_dimid(fname, ncid, prm->xdimname, &dimid_x);
    ncw_inq_dimid(fname, ncid, prm->ydimname, &dimid_y);
    ncw_inq_dimid(fname, ncid, prm->zdimname, &dimid_z);
    ncw_inq_dimlen(fname, ncid, dimid_x, &nx);
    ncw_inq_dimlen(fname, ncid, dimid_y, &ny);
    ncw_inq_dimlen(fname, ncid, dimid_z, &nz);

    ncw_inq_varid(fname, ncid, prm->xvarname, &varid_x);
    ncw_inq_varid(fname, ncid, prm->yvarname, &varid_y);
    ncw_inq_varid(fname, ncid, prm->zvarname, &varid_z);

    ncw_inq_varndims(fname, ncid, varid_x, &ndims_x);
    ncw_inq_varndims(fname, ncid, varid_y, &ndims_y);
    ncw_inq_varndims(fname, ncid, varid_z, &ndims_z);

    if (ndims_x == 1 && ndims_y == 1) {
        double* x;
        double* y;
        double* z;
        int i;
        double dx, dy;
        int periodic_x;

        x = malloc(nx * sizeof(double));
        y = malloc(ny * sizeof(double));
        z = malloc(nz * sizeof(double));

        ncw_get_var_double(fname, ncid, varid_x, x);
        ncw_get_var_double(fname, ncid, varid_y, y);
        ncw_get_var_double(fname, ncid, varid_z, z);

        periodic_x = fabs(fmod(2.0 * x[nx - 1] - x[nx - 2], 360.0) - x[0]) < EPSLON;

        dx = (x[nx - 1] - x[0]) / (double) (nx - 1);
        for (i = 1; i < (int) nx; ++i)
            if (fabs(x[i] - x[i - 1] - dx) / fabs(dx) > EPSLON)
                break;
        if (i != (int) nx)
            grid_setcoords(g, GRIDHTYPE_LATLON_IRREGULAR, NT_NONE, periodic_x, 0, nx, ny, nz, x, y, z);
        else {
            dy = (y[ny - 1] - y[0]) / (double) (ny - 1);
            for (i = 1; i < (int) ny; ++i)
                if (fabs(y[i] - y[i - 1] - dy) / fabs(dy) > EPSLON)
                    break;
            if (i != (int) ny)
                grid_setcoords(g, GRIDHTYPE_LATLON_IRREGULAR, NT_NONE, periodic_x, 0, nx, ny, nz, x, y, z);
            else
                grid_setcoords(g, GRIDHTYPE_LATLON_REGULAR, NT_NONE, periodic_x, 0, nx, ny, nz, x, y, z);
        }
    }
#if !defined(NO_GRIDUTILS)
    else if (ndims_x == 2 && ndims_y == 2) {
        double** x;
        double** y;
        double* z;

        x = alloc2d(ny, nx, sizeof(double));
        y = alloc2d(ny, nx, sizeof(double));
        z = malloc(nz * sizeof(double));

        ncw_get_var_double(fname, ncid, varid_x, x[0]);
        ncw_get_var_double(fname, ncid, varid_y, y[0]);
        ncw_get_var_double(fname, ncid, varid_z, z);

        grid_setcoords(g, GRIDHTYPE_CURVILINEAR, NT_COR, 0, 0, nx, ny, nz, x, y, z);
    }
#endif
    else
        enkf_quit("%s: could not determine the grid type", fname);

    if (prm->depthvarname != NULL) {
        float** depth = alloc2d(ny, nx, sizeof(float));

        ncw_inq_varid(fname, ncid, prm->depthvarname, &varid_depth);
        ncw_get_var_float(fname, ncid, varid_depth, depth[0]);
        g->depth = depth;
    }

    if (prm->levelvarname != NULL) {
        g->numlevels = alloc2d(ny, nx, sizeof(int));
        ncw_inq_varid(fname, ncid, prm->levelvarname, &varid_numlevels);
        ncw_get_var_int(fname, ncid, varid_numlevels, g->numlevels[0]);
        if (g->vtype == GRIDVTYPE_SIGMA) {
            int i, j;

            for (j = 0; j < ny; ++j)
                for (i = 0; i < nx; ++i)
                    g->numlevels[j][i] *= nz;
        }
    }
    ncw_close(fname, ncid);

    gridprm_print(prm, "    ");
    grid_print(g, "    ");

    return g;
}

/**
 */
void grid_destroy(grid* g)
{
    free(g->name);
    if (g->htype == GRIDHTYPE_LATLON_REGULAR || g->htype == GRIDHTYPE_LATLON_IRREGULAR)
        gnxy_simple_destroy(g->gridnodes_xy);
#if !defined(NO_GRIDUTILS)
    else if (g->htype == GRIDHTYPE_CURVILINEAR)
        gnxy_curv_destroy(g->gridnodes_xy);
#endif
    else
        enkf_quit("programming_error");
    if (g->gridnodes_z != NULL)
        gnz_destroy(g->gridnodes_z);
    if (g->numlevels != NULL)
        free2d(g->numlevels);
    if (g->depth != NULL)
        free2d(g->depth);

    free(g);
}

/**
 */
void grid_print(grid* g, char offset[])
{
    int nx, ny, nz;

    enkf_printf("%sgrid info:\n", offset);
    switch (g->htype) {
    case GRIDHTYPE_LATLON_REGULAR:
        enkf_printf("%s  hor type = LATLON_REGULAR\n", offset);
        break;
    case GRIDHTYPE_LATLON_IRREGULAR:
        enkf_printf("%s  hor type = LATLON_IRREGULAR\n", offset);
        break;
#if !defined(NO_GRIDUTILS)
    case GRIDHTYPE_CURVILINEAR:
        enkf_printf("%s  hor type = CURVILINEAR\n", offset);
        break;
#endif
    default:
        enkf_printf("%s  h type = NONE\n", offset);
    }
    enkf_printf("%s  periodic by X = %s\n", offset, grid_isperiodic_x(g) ? "yes" : "no");
    enkf_printf("%s  periodic by Y = %s\n", offset, grid_isperiodic_y(g) ? "yes" : "no");
    grid_getdims(g, &nx, &ny, &nz);
    enkf_printf("%s  dims = %d x %d x %d\n", offset, nx, ny, nz);
    if (g->lontype == LONTYPE_180)
        enkf_printf("%s  longitude range = [-180, 180]\n", offset);
    else if (g->lontype == LONTYPE_360)
        enkf_printf("%s  longitude range = [0, 360]\n", offset);
    else if (g->lontype == LONTYPE_NONE)
        enkf_printf("%s  longitude range = any\n", offset);
    switch (g->vtype) {
    case GRIDVTYPE_Z:
        enkf_printf("%s  vert type = Z\n", offset);
        break;
    case GRIDVTYPE_SIGMA:
        enkf_printf("%s  vert type = SIGMA\n", offset);
        break;
    default:
        enkf_printf("%s  vert type = NONE\n", offset);
    }
}

/**
 */
void grid_describeprm(void)
{
    enkf_printf("\n");
    enkf_printf("  Grid parameter file format:\n");
    enkf_printf("\n");
    enkf_printf("    NAME             = <name>\n");
    enkf_printf("    VTYPE            = { z | sigma }\n");
    enkf_printf("    DATA             = <data file name>\n");
    enkf_printf("    XDIMNAME         = <x dimension name>\n");
    enkf_printf("    YDIMNAME         = <y dimension name>\n");
    enkf_printf("    ZDIMNAME         = <z dimension name>\n");
    enkf_printf("    XVARNAME         = <x variable name>\n");
    enkf_printf("    YVARNAME         = <y variable name>\n");
    enkf_printf("    ZVARNAME         = <z variable name>\n");
    enkf_printf("    DEPTHVARNAME     = <depth variable name>\n");
    enkf_printf("    NUMLEVELSVARNAME = <# of levels variable name> (z)\n");
    enkf_printf("    MASKVARNAME      = <land mask variable name> (sigma)\n");
    enkf_printf("\n");
    enkf_printf("  [ <more of the above blocks> ]\n");
    enkf_printf("\n");
    enkf_printf("  Notes:\n");
    enkf_printf("    1. < ... > denotes a description of an entry\n");
    enkf_printf("\n");
}

/**
 */
void grid_getdims(grid* g, int* ni, int* nj, int* nk)
{
    if (ni != NULL) {
        if (g->htype == GRIDHTYPE_LATLON_REGULAR || g->htype == GRIDHTYPE_LATLON_IRREGULAR) {
            gnxy_simple* nodes = (gnxy_simple*) g->gridnodes_xy;

            *ni = nodes->nx;
            *nj = nodes->ny;
#if !defined(NO_GRIDUTILS)
        } else if (g->htype == GRIDHTYPE_CURVILINEAR) {
            gnxy_curv* nodes = (gnxy_curv*) g->gridnodes_xy;

            *ni = gridnodes_getnx(nodes->gn);
            *nj = gridnodes_getny(nodes->gn);
#endif
        } else
            enkf_quit("programming error");
    }
    if (nk != NULL)
        *nk = g->gridnodes_z->nz;
}

/**
 */
int grid_gettoplayerid(grid* g)
{
    gnz* nodes = g->gridnodes_z;
    int kmax = nodes->nz - 1;
    double* z = nodes->zt;

    return (z[0] < z[kmax]) ? 0 : kmax;
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
int grid_getlontype(grid* g)
{
    return g->lontype;
}

/**
 */
void grid_xy2fij(grid* g, double x, double y, double* fi, double* fj)
{
    g->xy2fij_fn(g, x, y, fi, fj);
}

/**
 */
void grid_z2fk(grid* g, double fi, double fj, double z, double* fk)
{
    g->z2fk_fn(g, fi, fj, z, fk);
}

/**
 */
void grid_fij2xy(grid* g, double fi, double fj, double* x, double* y)
{
    g->fij2xy_fn(g, fi, fj, x, y);
}

/**
 */
int grid_isperiodic_x(grid* g)
{
    if (g->htype == GRIDHTYPE_LATLON_REGULAR || g->htype == GRIDHTYPE_LATLON_IRREGULAR) {
        gnxy_simple* nodes = (gnxy_simple*) ((grid*) g)->gridnodes_xy;

        return nodes->periodic_x;
    }

    return 0;
}

/**
 */
int grid_isperiodic_y(grid* g)
{
    if (g->htype == GRIDHTYPE_LATLON_REGULAR || g->htype == GRIDHTYPE_LATLON_IRREGULAR) {
        gnxy_simple* nodes = (gnxy_simple*) ((grid*) g)->gridnodes_xy;

        return nodes->periodic_y;
    }

    return 0;
}

/**
 */
void grid_settocartesian_fn(grid* g, grid_tocartesian_fn fn)
{
    g->tocartesian_fn = fn;
}

/**
 */
void grid_tocartesian(grid* g, double* in, double* out)
{
    g->tocartesian_fn(in, out);
}
