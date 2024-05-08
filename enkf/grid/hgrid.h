/******************************************************************************
 *
 * File:        hgrid.h
 *
 * Created:     25/01/2022
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Purpose:     Header for `hgrid' object (generic horizontal grid).
 *
 * Description: Specific grid types are handled via field `gxy'.
 *
 * Revisions:
 *
 *****************************************************************************/

#if !defined(_HGRID_H)

typedef struct {
    void* parent;
    int type;
    int ni;
    int nj;
    /*
     * flag - whether the grid cordinates are lon/lat
     */
    int geographic;
    /*
     * longitude range = [lonbase, lonbase + 360]
     * NAN if not relevant for the grid
     */
    double lonbase;
    /*
     * flag - whether the grid is periodic in I direction
     * at the moment can be != 0 for rectangular geographic grids only
     */
    int periodic_i;
    void* gxy;
#if defined(ENKF_CALC)
    /*
     * used for calculating forecast obs with finite footprint
     */
    kdtree* nodetreeXYZ;
#endif
} hgrid;

struct gxy_1d;
typedef struct gxy_1d gxy_1d;

struct gxy_2d;
typedef struct gxy_2d gxy_2d;

hgrid* hgrid_create(void* prm, void* grid);
void hgrid_destroy(hgrid* hg);
void hgrid_describe(hgrid* hg, char* offset);

#if defined(ENKF_PREP) || defined(ENKF_CALC)
int hgrid_destroynodetree(hgrid* hg);
int hgrid_xy2fij(hgrid* hg, void* numlevels, double x, double y, double* fij);
void hgrid_ij2xy(hgrid* hg, int* ij, double* x, double* y);
#endif
#if defined(ENKF_CALC)
kdtree* hgrid_gettreeXYZ(hgrid* hg, int createifnull);
void hgrid_destroytreeXYZ(hgrid* hg);
#endif

#define _HGRID_H
#endif
