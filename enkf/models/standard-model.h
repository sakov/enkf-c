/******************************************************************************
 *
 * File:        standard-model.h       
 *
 * Created:     09/12/2012
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description: Contains standard model procedures
 *
 * Revisions:
 *
 *****************************************************************************/
#if !defined(_STANDARD_MODEL_H)

void standardmodel_getmemberfname(model* m, char ensdir[], char varname[], int mem, char fname[]);
int standardmodel_getmemberfname_async(model* m, char ensdir[], char varname[], char otname[], int mem, int t, char fname[]);
void standardmodel_getbgfname(model* m, char ensdir[], char varname[], char fname[]);
int standardmodel_getbgfname_async(model* m, char bgdir[], char varname[], char otname[], int t, char fname[]);
void standardmodel_adddata_2D(model* m, char* token, char* fname, int line);
void standardmodel_readfield(model* m, char fname[], int mem, int time, char varname[], int k, float* v);
void standardmodel_read3dfield(model* m, char fname[], int mem, int time, char varname[], float* v);
void standardmodel_writefield(model* m, char fname[], int time, char varname[], int k, float* v);

#define _STANDARD_MODEL_H
#endif
