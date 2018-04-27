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
    char* zvarname;
    char* zcvarname;
    char* hcvarname;            /* for ROMS sigma coords */
    char* depthvarname;
    char* levelvarnameentry;
    char* levelvarname;
    char* vdirection;

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
