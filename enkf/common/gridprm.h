/******************************************************************************
 *
 * File:        gridprm.h
 *
 * Created:     17/12/2014
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description:
 *
 * Revisions:
 *
 *****************************************************************************/

#if !defined(_GRIDPRM_H)

typedef struct {
    char* name;
    char* vtype;
#if !defined(NO_GRIDUTILS)
    char maptype;
#endif
    char* fname;
    char* xvarname;
    char* yvarname;
    char* angle;
    /*
     * Z stuff: zvarname and zcvarname hold the names of the layer centre
     * and layer corner corrdinate variable.
     */
    char* zvarname;
    char* zcvarname;

    /*
     * Sigma stuff:
     * This is generalised sigma as described in
     * https://https://www.myroms.org/wiki/Vertical_S-coordinate
     * For "normal" sigma specify C arrays (Cs_rho and Cs_w in ROMS) only via
     * entries CVARNAME (and, optionally) CCVARNAME along with WDIR ("fromsurf"
     * by default).
     */
    char* svarname;             /* variable name for sigma coordinate of
                                 * layer centres for non-uniform mappings */
    char* scvarname;            /* variable name for sigma coordinate of
                                 * layer corners for non-uniform mappings */
    char* hcvarname;            /* parameter for ROMS sigma coords */
    char* cvarname;             /* Cs_rho */
    char* ccvarname;            /* Cs_w */

    /*
     * Hybrid stuff: p(k) = A(k) + B(k) * (p1 - p2)
     * https://journals.ametsoc.org/doi/pdf/10.1175/2008MWR2537.1
     * See e.g. Eckerman (2009), eq. 6.
     */
    char* avarname;
    char* bvarname;
    char* acvarname;
    char* bcvarname;
    char* p1varname;            /* in atmosphere -- P of the surface layer */
    char* p2varname;            /* in atmosphere -- P of the top layer */

    char* depthvarname;
    char* levelvarnameentry;
    char* levelvarname;
    char* vdirection;

    /*
     * DA related stuff attached to grids
     */
    int stride;
    int sob_stride;
    double sfactor;

    /*
     * Vertical intervals for observation statistics
     */
    int nzints;
    zint* zints;
} gridprm;

void gridprm_create(char* fname, int* ngrid, gridprm** prm);
void gridprm_destroy(int ngrid, gridprm prm[]);
void gridprm_print(gridprm* prm, char offset[]);
int gridprm_getvtype(gridprm* prm);

#define _GRIDPRM_H
#endif
