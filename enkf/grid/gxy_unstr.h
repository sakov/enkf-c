/******************************************************************************
 *
 * File:        gxy_unstr.h        
 *
 * Created:     25/01/2022
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Purpose:     Header for `gxy_unstr' object (handling of unstructured
 *              horizontal grids).
 *
 * Description:
 *
 * Revisions:
 *
 *****************************************************************************/

#if !defined(_GXY_UNSTR_H)
#include "triangulation.h"

struct gxy_unstr;
typedef struct gxy_unstr gxy_unstr;

gxy_unstr* gxy_unstr_create(triangulation* d);
void gxy_unstr_destroy(gxy_unstr* gxy);
triangulation* gxy_unstr_gettriangulation(gxy_unstr* gxy);
point* gxy_unstr_getpoint(gxy_unstr* gxy, int i);
int gxy_unstr_xy2fij(gxy_unstr* gxy, double x, double y, double* fij);

#define _GXY_UNSTR_H
#endif
