/******************************************************************************
 *
 * File:        gxy_rect.h        
 *
 * Created:     25/01/2022
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Purpose:     Header for `gxy_rect' object (handling of "rectangular"
 *              horizontal grids).
 *
 * Description: By "rectangular" grid we understand a structured quadrilateral
 *              grid with node rows and columns aligned with coordinates.
 *
 * Revisions:
 *
 *****************************************************************************/

#if !defined(_GXY_RECT_H)
struct gxy_rect;
typedef struct gxy_rect gxy_rect;

gxy_rect* gxy_rect_create(hgrid* hg, int ni, int nj, double* x, double* y);
void gxy_rect_destroy(gxy_rect* gxy);
double* gxy_rect_getx(gxy_rect* gxy);
double* gxy_rect_gety(gxy_rect* gxy);
int gxy_rect_getni(gxy_rect* gxy);
int gxy_rect_getnj(gxy_rect* gxy);
void gxy_rect_xy2fij(gxy_rect* gxy, double x, double y, double* fij);
void gxy_rect_fij2xy(gxy_rect* gxy, double fi, double fj, double* x, double* y);

#define _GXY_RECT_H
#endif
