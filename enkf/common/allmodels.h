/******************************************************************************
 *
 * File:        allmodels.h        
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

#if !defined(_ALLMODELS_H)

#include "model.h"

typedef void (*modelsetup_fn) (model* m, char fname[]);

modelsetup_fn get_modelsetgridfn(char modeltype[]);
modelsetup_fn get_modelsetupfn(char modeltype[]);

#define _ALLMODELS_H
#endif
