/******************************************************************************
 *
 * File:        z-model.h        
 *
 * Created:     12/2012
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description:
 *
 * Revisions:   PS 1.7.2014 -- renamed from mom4.h
 *
 *****************************************************************************/

#if !defined(_ZMODEL_H)

void zmodel_setgrid(model* m, char gridprm[]);
void zmodel_setup(model* m, char modelprm[]);

#define _ZMODEL_H
#endif
