/******************************************************************************
 *
 * File:        prep.h        
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

#if !defined(_PREP_H)

void obs_add(observations* obs, model* m, obsmeta* meta);
void get_obsfiles(obsmeta* meta, int* nfiles, char*** fnames);
void print_obsstats(observations* obs, observations* sobs);

#define _PREP_H
#endif
