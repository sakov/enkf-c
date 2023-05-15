/******************************************************************************
 *
 * File:        grid.c        
 *
 * Created:     12/2012
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Purpose:     Code for `grid' object (generic grid).
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
#include <float.h>
#include <string.h>
#include <stdint.h>
#include "ncw.h"
#include "ncutils.h"
#include "definitions.h"
#include "utils.h"
#include "grid.h"
#include "vgrid.h"
#include "hgrid.h"
#include "enkfprm.h"
#include "gridprm.h"

#define EPS_IJ 1.0e-3
#define EPS_Z 1.0e-3

#define GRID_INC 10

struct grid {
    void* das;
    char* name;
    int id;
    int aliasid;

    vgrid* vgrid;
    hgrid* hgrid;

    /*
     * `numlevels' can hold either the number of levels ("z" or "numeric"
     * vertical grids) or the land mask ("sigma" vertical grids)
     */
    void* numlevels;            /* int** or int* */
    void* depth;                /* float** or float* */

    /*
     * `stride' for calculating ensemble transforms. "0" means to use the
     * common value defined in the top prm file. 
     */
    int stride;

    /*
     * Vertical intervals for observation statistics
     */
    int nzints;
    zint* zints;

    char* domainname;
};

/**
 */
grid* grid_create(void* p, int id, void** grids)
{
    gridprm* prm = (gridprm*) p;
    grid* g = calloc(1, sizeof(grid));

    g->das = NULL;
    g->name = strdup(prm->name);
    g->id = id;
    g->aliasid = -1;
    g->domainname = strdup(prm->domainname);    /* ("ALL" by default) */
    if (prm->stride != 0) {
        if (g->hgrid->type == GRIDHTYPE_UNSTRUCTURED && prm->stride != 1)
            enkf_quit("grid \"%s\": for unstructured grids STRIDE can be equal to 1 only", g->name);
        g->stride = prm->stride;
    }
    g->nzints = prm->nzints;
    if (prm->nzints > 0) {
        g->zints = malloc(g->nzints * sizeof(zint));
        memcpy(g->zints, prm->zints, g->nzints * sizeof(zint));
    }

    /*
     * set horizontal grid
     */
    if (prm->aliasname == NULL) {
        g->hgrid = hgrid_create(prm, g);
    } else {
        int i;

        for (i = 0; i < id; ++i) {
            grid* othergrid = (grid*) grids[i];

            if (strcmp(othergrid->name, prm->aliasname) == 0) {
                g->aliasid = othergrid->id;
                g->hgrid = othergrid->hgrid;

                break;
            }
        }
        if (i == id)
            enkf_quit("%s: %s: no NetCDF file or grid \"%s\" found\n", prm->prmfname, prm->name, prm->aliasname);
    }
    if (g->hgrid->type == GRIDHTYPE_UNSTRUCTURED || g->hgrid->type == GRIDHTYPE_1D)
        g->stride = 1;

    /*
     * set vertical grid
     */
    g->vgrid = vgrid_create(prm, g);

    /*
     * set depth data
     */
    if (prm->depthvarname != NULL) {
        int htype = grid_gethtype(g);
        int ncid, varid;

        ncw_open(prm->fname, NC_NOWRITE, &ncid);
        ncw_inq_varid(ncid, prm->depthvarname, &varid);
        if (htype == GRIDHTYPE_RECTANGULAR || htype == GRIDHTYPE_CURVILINEAR || htype == GRIDHTYPE_2D) {
            int ni, nj;
            float** depth;

            grid_getsize(g, &ni, &nj, NULL);
            {
                size_t dimlen[2] = { nj, ni };

                ncw_check_vardims(ncid, varid, 2, dimlen);
            }
            depth = alloc2d(nj, ni, sizeof(float));
            ncu_readvarfloat(ncid, varid, ni * nj, depth[0]);
            g->depth = depth;
            if (vgrid_signchanged(g->vgrid)) {
                float* v = depth[0];
                int i;

                for (i = 0; i < ni * nj; ++i)
                    v[i] = -v[i];
            }
        } else if (htype == GRIDHTYPE_UNSTRUCTURED || htype == GRIDHTYPE_1D) {
            int ni;

            grid_getsize(g, &ni, NULL, NULL);
            {
                size_t dimlen = ni;

                ncw_check_vardims(ncid, varid, 1, &dimlen);
            }
            g->depth = calloc(ni, sizeof(float));
            ncu_readvarfloat(ncid, varid, ni, g->depth);
            if (vgrid_signchanged(g->vgrid)) {
                float* v = g->depth;
                int i;

                for (i = 0; i < ni; ++i)
                    v[i] = -v[i];
            }
        } else
            enkf_quit("programming error");
        ncw_close(ncid);
    }

    /*
     * set 3D mask
     */
    if (prm->levelvarname != NULL) {
        int htype = grid_gethtype(g);
        int vtype = grid_getvtype(g);
        int ncid, varid;

        ncw_open(prm->fname, NC_NOWRITE, &ncid);
        ncw_inq_varid(ncid, prm->levelvarname, &varid);
        if (htype == GRIDHTYPE_RECTANGULAR || htype == GRIDHTYPE_CURVILINEAR || htype == GRIDHTYPE_2D) {
            int ni, nj, nk;
            int** numlevels;

            grid_getsize(g, &ni, &nj, &nk);
            {
                size_t dimlen[2] = { nj, ni };

                ncw_check_vardims(ncid, varid, 2, dimlen);
            }
            numlevels = alloc2d(nj, ni, sizeof(int));
            ncw_get_var_int(ncid, varid, numlevels[0]);
            if (vtype == GRIDVTYPE_SIGMA || vtype == GRIDVTYPE_HYBRID) {
                int i;

                for (i = 0; i < ni * nj; ++i) {
                    if (numlevels[0][i] != 0)
                        numlevels[0][i] = nk;
                }
            }
            g->numlevels = numlevels;
        } else if (htype == GRIDHTYPE_1D || htype == GRIDHTYPE_UNSTRUCTURED) {
            int ni, nk;
            int* numlevels;

            grid_getsize(g, &ni, NULL, &nk);
            {
                size_t dimlen = ni;

                ncw_check_vardims(ncid, varid, 1, &dimlen);
            }
            numlevels = calloc(ni, sizeof(int));
            ncw_get_var_int(ncid, varid, numlevels);
            if (vtype == GRIDVTYPE_SIGMA || vtype == GRIDVTYPE_HYBRID) {
                int i;

                for (i = 0; i < ni; ++i) {
                    if (numlevels[i] != 0)
                        numlevels[i] = nk;
                }
            }
            g->numlevels = numlevels;
        } else
            enkf_quit("programming error");
        ncw_close(ncid);
    }

    /*
     * set 3D mask if still necessary
     */
    if (g->numlevels == NULL) {
        int htype = grid_gethtype(g);
        int vtype = grid_getvtype(g);

        if (htype == GRIDHTYPE_RECTANGULAR || htype == GRIDHTYPE_CURVILINEAR || htype == GRIDHTYPE_2D) {
            int** numlevels;
            int ni, nj, nk;
            int i, j;

            grid_getsize(g, &ni, &nj, &nk);
            numlevels = alloc2d(nj, ni, sizeof(int));
            if (vtype == GRIDVTYPE_SIGMA || vtype == GRIDVTYPE_HYBRID) {
                for (j = 0; j < nj; ++j)
                    for (i = 0; i < ni; ++i)
                        if (g->depth == NULL || ((float**) g->depth)[j][i] > 0.0)
                            numlevels[j][i] = nk;
            } else {
                if (g->depth != NULL) {
                    for (j = 0; j < nj; ++j) {
                        for (i = 0; i < ni; ++i) {
                            double depth = ((float**) g->depth)[j][i];
                            double fk = NAN;

                            if (depth > 0.0) {
                                grid_z2fk(g, NULL, depth, &fk);
                                numlevels[j][i] = ceil(fk + 0.5);
                            }
                        }
                    }
                } else {
                    for (j = 0; j < nj; ++j)
                        for (i = 0; i < ni; ++i)
                            numlevels[j][i] = nk;
                }
            }
            g->numlevels = numlevels;
        } else if (htype == GRIDHTYPE_1D || htype == GRIDHTYPE_UNSTRUCTURED) {
            int* numlevels;
            int ni, nk;
            int i;

            grid_getsize(g, &ni, NULL, &nk);
            numlevels = calloc(ni, sizeof(int));
            if (vtype == GRIDVTYPE_SIGMA || vtype == GRIDVTYPE_HYBRID) {
                for (i = 0; i < ni; ++i)
                    if (g->depth == NULL || ((float*) g->depth)[i] > 0.0)
                        numlevels[i] = nk;
            } else {
                if (g->depth != NULL) {
                    for (i = 0; i < ni; ++i) {
                        double depth = ((float*) g->depth)[i];
                        double fk = NAN;

                        if (depth > 0.0) {
                            grid_z2fk(g, NULL, depth, &fk);
                            numlevels[i] = ceil(fk + 0.5);
                        }
                    }
                } else {
                    for (i = 0; i < ni; ++i)
                        numlevels[i] = nk;
                }
            }
            g->numlevels = numlevels;
        } else
            enkf_quit("programming error");
    }

    /*
     * print settings and results to log
     */
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
    vgrid_destroy(g->vgrid);

    if (g->aliasid < 0)
        hgrid_destroy(g->hgrid);

    if (g->numlevels != NULL)
        free(g->numlevels);
    if (g->depth != NULL)
        free(g->depth);
    if (g->nzints > 0)
        free(g->zints);
    free(g);
}

/**
 */
void* grid_getvgrid(grid* g)
{
    return g->vgrid;
}

/**
 */
int grid_gethtype(grid* g)
{
    return g->hgrid->type;
}

/**
 */
int grid_getvtype(grid* g)
{
    return g->vgrid->type;
}

/**
 */
void* grid_getdepth(grid* g)
{
    return g->depth;
}

/**
 */
void grid_print(grid* g, char offset[])
{
    int ni, nj, nk;

    if (rank != 0)
        return;

    enkf_printf("%sgrid info:\n", offset);
    grid_getsize(g, &ni, &nj, &nk);
    enkf_printf("%s  size = %d x %d x %d\n", offset, ni, nj, nk);
    if (g->aliasid < 0)
        hgrid_describe(g->hgrid, offset);

    else
        enkf_printf("%s  horizontal grid -- aliased to grid #%d\n", offset, g->aliasid);
    vgrid_describe(g->vgrid, offset);

    if (g->stride != 1)
        enkf_printf("%s  STRIDE = %d\n", offset, g->stride);
#if defined(ENKF_CALC)
    if (g->hgrid->nodetreeXYZ != NULL)
        kd_printinfo(g->hgrid->nodetreeXYZ, "      ");
#endif
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
    enkf_printf("    DATA             = <data file name>\n");
    enkf_printf("    (either)\n");
    enkf_printf("      HTYPE          = { rect | curv | unstr }\n");
    enkf_printf("      XVARNAME       = <X variable name>\n");
    enkf_printf("      YVARNAME       = <Y variable name>\n");
    enkf_printf("      (if htype = unstr)\n");
    enkf_printf("      [ TRIVARNAME    = <triangle vertice IDS variable name> ]\n");
    enkf_printf("      [ TRINEIVARNAME = <triangle neighbour IDS variable name> ]\n");
    enkf_printf("      (end if)\n");
    enkf_printf("    (or)\n");
    enkf_printf("      HGRIDFROM      = <grid name>\n");
    enkf_printf("    (end either)\n");
    enkf_printf("    VTYPE            = { z | sigma | hybrid | numeric | none }\n");
    enkf_printf("  [ VDIR             = { fromsurf* | tosurf } ]\n");
    enkf_printf("  [ GEOGRAPHIC       = { yes* | no } ]\n");
    enkf_printf("    (if vtype = z)\n");
    enkf_printf("      ZVARNAME       = <Z variable name>\n");
    enkf_printf("    [ ZCVARNAME      = <ZC variable name> ]\n");
    enkf_printf("    [ NUMLEVELSVARNAME = <# of levels variable name> ]\n");
    enkf_printf("    [ DEPTHVARNAME   = <depth variable name> ]\n");
    enkf_printf("    (else if vtype = sigma)\n");
    enkf_printf("      CVARNAME       = <Cs_rho variable name>\n");
    enkf_printf("    [ CCVARNAME      = <Cs_w variable name> ]\n");
    enkf_printf("    [ SVARNAME       = <s_rho variable name> ]            (uniform*)\n");
    enkf_printf("    [ SCVARNAME      = <s_w variable name> ]              (uniform*)\n");
    enkf_printf("    [ HCVARNAME      = <hc variable name> ]               (0.0*)\n");
    enkf_printf("    [ DEPTHVARNAME   = <depth variable name> ]\n");
    enkf_printf("    [ MASKVARNAME    = <land mask variable name> ]\n");
    enkf_printf("    (else if vtype = hybrid)\n");
    enkf_printf("      AVARNAME       = <A variable name>\n");
    enkf_printf("      BVARNAME       = <B variable name>\n");
    enkf_printf("    [ ACVARNAME      = <AC variable name> ]\n");
    enkf_printf("    [ BCVARNAME      = <BC variable name> ]\n");
    enkf_printf("      P1VARNAME      = <P1 variable name>\n");
    enkf_printf("      P2VARNAME      = <P2 variable name>\n");
    enkf_printf("    [ MASKVARNAME    = <land mask variable name> ]\n");
    enkf_printf("    (else if vtype = numeric)\n");
    enkf_printf("      ZVARNAME       = <Z variable name>\n");
    enkf_printf("    [ ZCVARNAME      = <ZC variable name> ]\n");
    enkf_printf("    [ NUMLEVELSVARNAME = <# of levels variable name> ]\n");
    enkf_printf("    [ DEPTHVARNAME   = <depth variable name> ]\n");
    enkf_printf("    (end if)\n");
    enkf_printf("  [ STRIDE           = <stride for ensemble transforms> ] (1*)\n");
    enkf_printf("  [ SOBSTRIDE        = <stride for superobing> ]          (1*)\n");
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
    if (ni != NULL)
        *ni = g->hgrid->ni;

    if (nj != NULL)
        *nj = g->hgrid->nj;

    if (nk != NULL)
        *nk = g->vgrid->nk;
}

/**
 */
int grid_getsurflayerid(grid* g)
{
    return (g->vgrid->direction == GRIDVDIR_FROMSURF) ? 0 : g->vgrid->nk - 1;
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
void* grid_getnumlevels(grid* g)
{
    return g->numlevels;
}

/**
 */
double grid_getlonbase(grid* g)
{
    return g->hgrid->lonbase;
}

/**
 */
int grid_getstride(grid* g)
{
    return g->stride;
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
#if defined(ENKF_PREP) || defined (ENKF_CALC)
int grid_xy2fij(grid* g, double x, double y, double* fij)
{
    return hgrid_xy2fij(g->hgrid, g->numlevels, x, y, fij);
}
#endif

/**
 */
int grid_z2fk(grid* g, double* fij, double z, double* fk)
{
    int vtype = g->vgrid->type;
    int htype = g->hgrid->type;

    if (fij != NULL && isnan(fij[0] + fij[1])) {
        *fk = NAN;
        return STATUS_OUTSIDEGRID;
    }

    vgrid_z2fk(g->vgrid, fij, z, fk);

    if (isnan(*fk))
        return STATUS_OUTSIDEGRID;

    /*
     * Check depth for z-grid.
     */
    if (vtype == GRIDVTYPE_Z && fij != NULL) {
        if (htype == GRIDHTYPE_RECTANGULAR || htype == GRIDHTYPE_CURVILINEAR) {
            int** numlevels = g->numlevels;
            int ksurf = grid_getsurflayerid(g);
            int i1 = floor(fij[0]);
            int i2 = ceil(fij[0]);
            int j1 = floor(fij[1]);
            int j2 = ceil(fij[1]);
            int k;
            int ni, nj;

            grid_getsize(g, &ni, &nj, NULL);
            if (i1 == -1)
                i1 = (g->hgrid->periodic_i) ? ni - 1 : i2;
            if (i2 == ni)
                i2 = (g->hgrid->periodic_i) ? 0 : i1;
            if (j1 == -1)
                j1 = j2;
            if (j2 == nj)
                j2 = j1;

            /*
             * see the note on the similar bit in grid_xy2fij()
             */
            *fk = (float) (*fk);
            k = (ksurf == 0) ? floor(*fk) : ksurf - ceil(*fk);

            if (numlevels[j1][i1] <= k && numlevels[j1][i2] <= k && numlevels[j2][i1] <= k && numlevels[j2][i2] <= k) {
                *fk = NAN;
                return STATUS_LAND;
            } else if (numlevels[j1][i1] <= k || numlevels[j1][i2] <= k || numlevels[j2][i1] <= k || numlevels[j2][i2] <= k) {
                double v = grid_interpolate2d(g, fij, g->depth);

                if (z > v)
                    return STATUS_LAND;
            }
        } else if (htype == GRIDHTYPE_UNSTRUCTURED) {
            int* numlevels = g->numlevels;
            int ksurf = grid_getsurflayerid(g);
            int i0 = floor(fij[0]);
            int i1 = floor(fij[1]);
            int i2 = floor(fij[2]);
            int k;

            *fk = (float) (*fk);
            k = (ksurf == 0) ? floor(*fk) : ksurf - ceil(*fk);

            if (numlevels[i0] <= k && numlevels[i1] <= k && numlevels[i2] <= k) {
                *fk = NAN;
                return STATUS_LAND;
            } else if (numlevels[i0] <= k || numlevels[i1] <= k || numlevels[i2] <= k) {
                double v = grid_interpolate2d(g, fij, g->depth);

                if (z > v)
                    return STATUS_LAND;
            }
        }
    }

    return STATUS_OK;
}

/**
 */
#if defined(ENKF_PREP) || defined(ENKF_CALC)
int grid_z2fk_f(grid* g, double* fij, double z, float* fk)
{
    double fk_d;
    int status = grid_z2fk(g, fij, z, &fk_d);

    *fk = (float) fk_d;

    return status;
}
#endif

/**
 */
#if defined(ENKF_PREP) || defined(ENKF_CALC)
void grid_ij2xy(grid* g, int* ij, double* x, double* y)
{
    hgrid_ij2xy(g->hgrid, ij, x, y);
}
#endif

/**
 */
#if defined(ENKF_PREP) || defined(ENKF_CALC)
/** Check whether a point in grid space can be considered watered for
 ** interpolation purposes. The answer is positive if at least one of
 ** (generally) 8 corner nodes of the 3D grid cell the point belongs to is
 ** watered. For points  exactly on the cell face/edge/node the number of nodes
 ** checked effectively reduces to 4/2/1. 
 */
int grid_island(grid* g, double* fij, double fk)
{
    if (g->hgrid->type == GRIDHTYPE_RECTANGULAR || g->hgrid->type == GRIDHTYPE_CURVILINEAR) {
        int i1 = (int) floor(fij[0]);
        int i2 = (int) ceil(fij[0]);
        int j1 = (int) floor(fij[1]);
        int j2 = (int) ceil(fij[1]);
        int ksurf = grid_getsurflayerid(g);
        int** numlevels = g->numlevels;
        int k;

        if (isnan(fk))
            k = ksurf;
        else {
            if (ceil(fk) != floor(fk))
                k = (ksurf == 0) ? ceil(fk) : ksurf - floor(fk);
            else
                k = (ksurf == 0) ? ceil(fk) + 1 : ksurf - floor(fk) + 1;
        }

        if (i1 == -1)
            i1 = (g->hgrid->periodic_i) ? g->hgrid->ni - 1 : i2;

        if (i2 == g->hgrid->ni)
            i2 = (g->hgrid->periodic_i) ? 0 : i1;
        if (j1 == -1)
            j1 = j2;
        if (j2 == g->hgrid->nj)
            j2 = j1;

        if (k == 0)
            return numlevels[j1][i1] == 0 && numlevels[j1][i2] == 0 && numlevels[j2][i1] == 0 && numlevels[j2][i2] == 0;
        else
            return numlevels[j1][i1] < k && numlevels[j1][i2] < k && numlevels[j2][i1] < k && numlevels[j2][i2] < k;
    } else if (g->hgrid->type == GRIDHTYPE_UNSTRUCTURED) {
        int i1 = (int) floor(fij[0]);
        int i2 = (int) ceil(fij[0]);
        int ksurf = grid_getsurflayerid(g);
        int* numlevels = g->numlevels;
        int k;

        if (ceil(fk) != floor(fk))
            k = (ksurf == 0) ? ceil(fk) : ksurf - floor(fk);
        else
            k = (ksurf == 0) ? ceil(fk) + 1 : ksurf - floor(fk) + 1;

        if (i1 == -1)
            i1 = i2;
        if (i2 == g->hgrid->ni)
            i2 = i1;

        if (k == 0)
            return numlevels[i1] == 0 && numlevels[i2] == 0;
        else
            return numlevels[i1] < k && numlevels[i2] < k;
    }

    enkf_quit("programming error");
    return 0;
}
#endif

/**
 */
float grid_interpolate2d(grid* g, double* fij, void* v)
{
    if (v == NULL)
        return NAN;

    if (g->hgrid->type == GRIDHTYPE_RECTANGULAR || g->hgrid->type == GRIDHTYPE_CURVILINEAR || g->hgrid->type == GRIDHTYPE_2D) {
        int ni, nj;

        grid_getsize(g, &ni, &nj, NULL);
        return interpolate2d_structured(fij, ni, nj, v, g->numlevels, g->hgrid->periodic_i);
    } else if (g->hgrid->type == GRIDHTYPE_UNSTRUCTURED) {
        return interpolate2d_unstructured(fij, v, g->numlevels);
    } else
        enkf_quit("programming error");

    return NAN;
}

/**
 */
#if defined(ENKF_PREP) || defined(ENKF_CALC)
float grid_interpolate3d(grid* g, double* fij, float fk, void* v)
{
    if (v == NULL)
        return NAN;

    if (g->hgrid->type == GRIDHTYPE_RECTANGULAR || g->hgrid->type == GRIDHTYPE_CURVILINEAR)
        return interpolate3d_structured(fij, (double) fk, g->hgrid->ni, g->hgrid->nj, g->vgrid->nk, grid_getsurflayerid(g), v, g->numlevels, g->hgrid->periodic_i);
    else if (g->hgrid->type == GRIDHTYPE_UNSTRUCTURED)
        return interpolate3d_unstructured(fij, (double) fk, g->vgrid->nk, grid_getsurflayerid(g), v, g->numlevels);
    else
        enkf_quit("programming error");

    return NAN;
}
#endif

/**
 */
int grid_isperiodic_i(grid* g)
{
    return g->hgrid->periodic_i;
}

/**
 */
int grid_isstructured(grid* g)
{
    if (g->hgrid->type == GRIDHTYPE_RECTANGULAR || g->hgrid->type == GRIDHTYPE_CURVILINEAR || g->hgrid->type == GRIDHTYPE_2D)
        return 1;
    if (g->hgrid->type == GRIDHTYPE_UNSTRUCTURED || g->hgrid->type == GRIDHTYPE_1D)
        return 0;

    enkf_quit("programming error");
    return 1;
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
void grids_create(enkfprm* prm, int* ngrid, void*** grids)
{
    int n = 0;
    gridprm* gprm = NULL;
    int i;

    gridprm_create(prm, &n, &gprm);
    assert(n > 0);

    for (i = 0; i < n; ++i) {
        grid* g = NULL;

        g = grid_create(&gprm[i], i, *grids);
        grids_addgrid(ngrid, grids, g);

        if (grid_getstride(g) == 0)
            g->stride = prm->stride;
    }
    gridprm_destroy(n, gprm);
}

/**
 */
void grids_destroy(int ngrid, void** grids)
{
    int i;

    if (ngrid == 0)
        return;
    for (i = 0; i < ngrid; ++i)
        grid_destroy(grids[i]);
    free(grids);
}

/**
 */
#if defined(ENKF_CALC)
kdtree* grid_gettreeXYZ(grid* g, int createifnull)
{
    return hgrid_gettreeXYZ(g->hgrid, createifnull);
}
#endif

/**
 */
#if defined(ENKF_CALC)
void grid_destroytreeXYZ(grid* g)
{
    hgrid_destroytreeXYZ(g->hgrid);
}
#endif

/**
 */
#if defined(ENKF_CALC)
int grids_destroyhtrees(int ngrid, void** grids)
{
    int hadtrees, i;

    for (i = 0, hadtrees = 0; i < ngrid; ++i) {
        grid* g = grids[i];

        if (hgrid_destroynodetree(g->hgrid))
             hadtrees = 1;
    }

    return hadtrees;
}
#endif

/**
 */
char* grid_getdomainname(grid* g)
{
    return g->domainname;
}

/**
 */
int grid_getaliasid(grid* g)
{
    return g->aliasid;
}
