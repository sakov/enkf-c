/******************************************************************************
 *
 * File:        vgrid.h
 *
 * Created:     25/01/2022
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Purpose:     Header for `vgrid' object (generic vertical grid).
 *
 * Description: Specific grid types are handled via field `gz'.
 *
 * Revisions:
 *
 *****************************************************************************/

#if !defined(_VGRID_H)

#define GRIDVDIR_FROMSURF 0
#define GRIDVDIR_TOSURF 1

typedef struct {
    void* parent;
    int type;
    int nk;
    int direction;
    void* gz;
} vgrid;

vgrid* vgrid_create(void* gridprm, void* grid);
void vgrid_destroy(vgrid* vg);
void vgrid_describe(vgrid* vg, char* offset);
void vgrid_z2fk(vgrid* vg, double* fij, double z, double* fk);
int vgrid_signchanged(vgrid* vg);

#define _VGRID_H
#endif
