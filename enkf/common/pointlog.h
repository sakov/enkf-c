/******************************************************************************
 *
 * File:        pointlog.h       
 *
 * Created:     7/10/2013
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description:
 *
 * Revisions:
 *
 *****************************************************************************/

#if !defined(_POINTLOG_H)

void plog_write(dasystem* das, int i, int j, double lon, double lat, double depth, int p, int* lobs, double* lcoeffs, double* s, double* S, double* transform);
void plog_writeactualtransform(dasystem* das, int i, int j, float* transform);
void plog_definestatevars(dasystem* das);
void plog_writestatevars(dasystem* das, int nfields, void** fieldbuffer, field* fields, int isanalysis);
void plog_assemblestatevars(dasystem* das);

#define _POINTLOG_H
#endif
