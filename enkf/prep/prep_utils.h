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

void obs_add(observations* obs, model* m, obssection* section, int nexclude, obsregion* exclude);
int obs_checkforland(observations* obs, model* m);
void get_obsfiles(obssection* section, int* nfiles, char*** fnames);
void print_obsstats(observations* obs, observations* sobs);
char* get_lonname(int ncid, char* lonname);
char* get_latname(int ncid, char* latname);
char* get_zname(int ncid, char* zname);
int get_insttag(int ncid, char* varname, char* insttag);
void get_qcflags(obssection* section, int* nqcflags, char*** qcflagname, uint32_t** qcflagvals);
void get_time(obssection* section, int ncid, size_t* size, double** time);
void list_readers(void);
void describe_reader(char readertag[]);
void describe_commonreaderparams(void);
void describe_commongenericreaderparams(void);

#define _PREP_UTILS_H
#endif
