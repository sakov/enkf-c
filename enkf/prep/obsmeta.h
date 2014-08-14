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

#define STDTYPE_NONE 0
#define STDTYPE_VALUE 1
#define STDTYPE_FILE 2

#define ARITHMETIC_EQ 0
#define ARITHMETIC_PLUS 1
#define ARITHMETIC_MULT 2
#define ARITHMETIC_MIN 3
#define ARITHMETIC_MAX 4

typedef struct {
    char* product;
    char* reader;
    char* type;
    int nfiles;
    char** fnames;
    int nstds;
    int* stdtypes;
    int* stdops;
    void** stds;
    char** varnames;
} obsmeta;

void read_obsmeta(enkfprm* prm, int* nmeta, obsmeta** meta);
void clean_obsmeta(int n, obsmeta meta[]);

#define _OBSMETA_H
#endif
