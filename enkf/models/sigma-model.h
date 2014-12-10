/******************************************************************************
 *
 * File:        sigma-model.h        
 *
 * Created:     03/12/2014
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description:
 *
 * Revisions:
 *
 *****************************************************************************/

#if !defined(_SIGMA_MODEL_H)

void smodel_setgrids(model* m, char gridprm[]);
void smodel_setup(model* m, char modelprm[]);

#define _SIGMA_MODEL_H
#endif
