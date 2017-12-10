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
    char* xdimname;
    char* ydimname;
    char* zdimname;
    char* xvarname;
    char* yvarname;
    char* zvarname;
    char* depthvarname;
    char* levelvarnameentry;
    char* levelvarname;
    /*
     * Hybrid stuff: p(k) = A(k) + B(k) * (p1 - p2)
     * See e.g. Eckerman (2008), eq. 6.
     */
    char* avarname;
    char* bvarname;
    char* p1varname;            /* in atmosphere -- P of the surface layer */
    char* p2varname;            /* in atmosphere - P of the top layer */
    /*
     * hybrid stuff end 
     */
    int stride;
    double sfactor;
} gridprm;

void gridprm_create(char* fname, int* ngrid, gridprm** prm);
void gridprm_destroy(int ngrid, gridprm prm[]);
void gridprm_print(gridprm* prm, char offset[]);
int gridprm_getvtype(gridprm* prm);

#define _GRIDPRM_H
#endif
