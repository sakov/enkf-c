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

void das_writespread(dasystem* das);
void das_writevcorrs(dasystem* das);
void das_writevcorrs_with(dasystem* das, char* varname, int level, int docorr);

#define _DIAGS_H
#endif
