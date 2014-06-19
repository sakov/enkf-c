/******************************************************************************
 *
 * File:        obsmeta.h        
 *
 * Created:     12/2012
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description:
 *
 * Revisions:
 *
 *****************************************************************************/

#if !defined(_OBSMETA_H)

typedef struct {
    char* product;
    char* reader;
    char* type;
    int nfiles;
    char** fnames;
    int nstds;
    int* stdtypes;
    void** stds;
    char** varnames;
} obsmeta;

void read_obsmeta(enkfprm* prm, int* nmeta, obsmeta** meta);
void clean_obsmeta(int n, obsmeta meta[]);

#define _OBSMETA_H
#endif
