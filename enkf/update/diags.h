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

void das_allocatespread(dasystem* das, char fname[]);
void das_writespread(dasystem* das, int nfields, void** fieldbuffer, field fields[], int isanalysis);
void das_assemblespread(dasystem* das);
void das_allocateinflation(dasystem* das, char fname[]);
void das_writeinflation(dasystem* das, field* f, int j, float* v);
void das_assembleinflation(dasystem* das);
void das_writevcorrs(dasystem* das);

#define _DIAGS_H
#endif
