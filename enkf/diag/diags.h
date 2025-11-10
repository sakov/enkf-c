/******************************************************************************
 *
 * File:        diags.h        
 *
 * Created:     15/07/2019
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description:
 *
 * Revisions:
 *
 *****************************************************************************/

#if !defined(_DIAGS_H)

#define CALC_COV 0
#define CALC_CORR 1
#define CALC_SENS 2

#if defined(ENS_DIAG)
void das_writespread(dasystem* das);
#endif
void das_writevcorrs(dasystem* das);
void das_writevcorrs_with(dasystem* das, char* varname, int level, int calctype);

#define _DIAGS_H
#endif
