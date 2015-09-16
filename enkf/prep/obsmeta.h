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

#include "enkfprm.h"

#define STDTYPE_NONE 0
#define STDTYPE_VALUE 1
#define STDTYPE_FILE 2

#define ARITHMETIC_EQ 0
#define ARITHMETIC_PLUS 1
#define ARITHMETIC_MULT 2
#define ARITHMETIC_MIN 3
#define ARITHMETIC_MAX 4

typedef struct {
    int type;
    int op;
    void* data;                 /* double[1] or char* */
    char* varname;              /* only needed if data = char* */
} metastd;

typedef struct {
    char* name;
    char* value;
} metapar;

typedef struct {
    char* product;
    char* reader;
    char* type;
    int nfiles;
    char** fnames;
    /*
     * variables for specifying observation error
     */
    int nstds;
    metastd* stds;
    /*
     * extra parameters
     */
    int npars;
    metapar* pars;
} obsmeta;

void obsmeta_read(enkfprm* prm, int* nmeta, obsmeta** meta);
void obsmeta_destroy(int n, obsmeta meta[]);
void obsmeta_describeprm(void);

#define _OBSMETA_H
#endif
