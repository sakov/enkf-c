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
int obs_checkforland(observations* obs, model* m);
void get_obsfiles(obsmeta* meta, int* nfiles, char*** fnames);
void print_obsstats(observations* obs, observations* sobs);
char* get_lonname(int ncid, char* lonname);
char* get_latname(int ncid, char* latname);
char* get_zname(int ncid, char* zname);
char* get_timename(int ncid, char* timename);
void get_qcflags(obsmeta* meta, int* nqcflags, char*** qcflagname, uint32_t** qcflagvals);

#define _PREP_UTILS_H
#endif
