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
} std_entry;

typedef struct {
    char* name;
    char* value;
} par_entry;

typedef struct {
    int id;
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
    std_entry* estds;
    /*
     * extra parameters
     */
    int npars;
    par_entry* pars;
} obssection;

void obsprm_read(char fname[], int* nsection, obssection** sections, int* nexclude, obsregion** exclude);
void obsprm_destroy(int nsection, obssection sections[], int nexclude, obsregion* exclude);
void obsprm_describeprm(void);

#define _OBSPRM_H
#endif
