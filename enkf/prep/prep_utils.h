/******************************************************************************
 *
 * File:        prep_utils.h        
 *
 * Created:     12/2012
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description:
 *
 * Revisions:   22.07.2015 PS Renamed from prep.h
 *
 *****************************************************************************/

#if !defined(_PREP_UTILS_H)

void obs_add(observations* obs, model* m, obsmeta* meta);
void get_obsfiles(obsmeta* meta, int* nfiles, char*** fnames);
void print_obsstats(observations* obs, observations* sobs);

#define _PREP_UTILS_H
#endif
