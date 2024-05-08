/******************************************************************************
 *
 * File:        gxy_curv2.h
 *
 * Created:     06/05/2024
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Purpose:     Handles geographic curvilinear grids using stereographic
 *              projections.
 *
 * Description: Derived from gxy_curv.h.
 *
 * Revisions:   
 *
 *****************************************************************************/

#if !defined(_GXY_CURV2_H)
struct gxy_curv2;
typedef struct gxy_curv2 gxy_curv2;

gxy_curv2* gxy_curv2_create(hgrid* hg, int ni, int nj, double** x, double** y, int** mask);
void gxy_curv2_destroy(gxy_curv2* gxy);
int gxy_curv2_destroykdtree(gxy_curv2* gxy);
double** gxy_curv2_getx(gxy_curv2* gxy);
double** gxy_curv2_gety(gxy_curv2* gxy);
int gxy_curv2_xy2fij(gxy_curv2* gxy, double x, double y, double* fij);

#define _GXY_CURV2_H
#endif
