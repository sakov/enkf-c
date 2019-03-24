/******************************************************************************
 *
 * File:        obsprm.h        
 *
 * Created:     12/2012
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description:
 *
 * Revisions:   18062018 PS Renamed from obsmeta.h
 *
 *****************************************************************************/

#if !defined(_OBSPRM_H)

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
    char* prmfname;
    char* product;
    char* reader;
    char* type;
    int nfiles;
    char** fnames;
    /*
     * variables for specifying observation error
     */
    int nestds;
    metastd* estds;
    /*
     * extra parameters
     */
    int npars;
    metapar* pars;
} obsmeta;

void obsprm_read(char fname[], int* nmeta, obsmeta** meta);
void obsprm_destroy(int n, obsmeta meta[]);
void obsprm_describeprm(void);

#define _OBSPRM_H
#endif
