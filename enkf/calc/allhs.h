/******************************************************************************
 *
 * File:        allhs.h        
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

#if !defined(_ALLHS_H)

#include "dasystem.h"

/* I would prefer to avoid dependence of on `dasystem', and rather depend on
 * `model' and `observations'. But some functions need additional data, e.g.
 * mapping of SLA requires MSL, and those with bias estimation require the bias
 * fields. So `dasystem' is supposed to provide a generic access to the
 * required information. */

typedef void (*H_fn) (dasystem* das, int nobs, int obsids[], char fname[], int mem, int t, char varname[], void* psrc, ENSOBSTYPE dst[]);

void describe_hentries(char* obstypename);
H_fn getH(char obstypename[], char mappingname[]);

#define _ALLHS_H
#endif
