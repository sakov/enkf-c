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
} gridprm;

void gridprm_create(char* fname, int* ngrid, gridprm** prm);
void gridprm_destroy(int ngrid, gridprm prm[]);
void gridprm_print(gridprm* prm, char offset[]);
int gridprm_getvtype(gridprm* prm);

#define _GRIDPRM_H
#endif
