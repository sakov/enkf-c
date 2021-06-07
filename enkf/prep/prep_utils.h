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

void obs_add(observations* obs, model* m, obsmeta* meta, int nexclude, obsregion* exclude);
int obs_checkforland(observations* obs, model* m);
void get_obsfiles(obsmeta* meta, int* nfiles, char*** fnames);
void print_obsstats(observations* obs, observations* sobs);
char* get_lonname(int ncid, char* lonname);
char* get_latname(int ncid, char* latname);
char* get_zname(int ncid, char* zname);
int get_insttag(int ncid, char* varname, char* insttag);
void get_qcflags(obsmeta* meta, int* nqcflags, char*** qcflagname, uint32_t** qcflagvals);
void get_time(obsmeta* meta, int ncid, size_t* size, double** time);
void list_readers(void);
void describe_reader(char readertag[]);

#define _PREP_UTILS_H
#endif
