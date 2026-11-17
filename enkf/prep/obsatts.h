/******************************************************************************
 *
 * File:        obsatts.h        
 *
 * Created:     17/11/2026
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description:
 *
 * Revisions:   Observation reader helper
 *
 *****************************************************************************/

#if !defined(_OBSATTS_H)

typedef struct {
    char* varname;
    int sign;                   /* +1 for ADDVAR, -1 for SUBVAR */
    double* v;                  /* flat [nobs] */
} obsaddvar;

typedef struct {
    /*
     * parsed parameter strings 
     */
    char* stdname;
    char* estdname;
    char* batchname;
    char* instattname;
    char* instprefix;
    char instrument[MAXSTRLEN];
    int naddvar;
    obsaddvar* addvars;

    /*
     * loaded observation attributes 
     */
    double var_estd;
    double* std;                /* flat [nobs] or NULL */
    double* estd;               /* flat [nobs] or NULL */
    int* batch;                 /* flat [nobs] or NULL */
    int nqcflagvars;
    char** qcflagvarnames;
    uint32_t* qcflagmasks;
    int32_t** qcflag;           /* [nqcflagvars][nobs] or NULL */
    int instid;
} obsatts;

void obsatts_init(obsatts* oa);
int obsatts_parse_param(obsatts* oa, const char* name, char* value);
void obsatts_read(obsatts* oa, obssection* section, int ncid, const char* varname, size_t nobs, size_t nprof, size_t nz, int check_allmissing, int allow_char_batch_qc, observations* obs);
int obsatts_valid(const obsatts* oa, size_t idx);
void obsatts_apply(const obsatts* oa, size_t idx, int default_batch, observation* o);
void obsatts_free(obsatts* oa);

#define _OBSATTS_H
#endif
