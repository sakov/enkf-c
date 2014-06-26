/******************************************************************************
 *
 * File:        mom4.h        
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

#if !defined(_MOM4_H)

void mom4_setgrid(model* m, char gridprm[]);
void mom4_setup(model* m, char modelprm[]);

#define _MOM4_H
#endif
