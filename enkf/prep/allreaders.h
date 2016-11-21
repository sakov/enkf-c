/******************************************************************************
 *
 * File:        allreaders.h        
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

#if !defined(_ALLREADERS_H)

typedef void (*obsread_fn) (char* fname, int fid, obsmeta* meta, model* m, observations* obs);

void describe_readers(void);
obsread_fn get_obsreadfn(obsmeta* m);

void reader_rads_standard(char* fname, int fid, obsmeta* meta, model* m, observations* obs);
void reader_rads_standard2(char* fname, int fid, obsmeta* meta, model* m, observations* obs);
void reader_navo_standard(char* fname, int fid, obsmeta* meta, model* m, observations* obs);
void reader_windsat_standard(char* fname, int fid, obsmeta* meta, model* m, observations* obs);
void reader_cars_standard(char* fname, int fid, obsmeta* meta, model* m, observations* obs);
void reader_mmt_standard(char* fname, int fid, obsmeta* meta, model* m, observations* obs);
void reader_amsr2_standard(char* fname, int fid, obsmeta* meta, model* m, observations* obs);
void reader_amsre_standard(char* fname, int fid, obsmeta* meta, model* m, observations* obs);
void reader_pathfinder_standard(char* fname, int fid, obsmeta* meta, model* m, observations* obs);
void reader_aquarius_standard(char* fname, int fid, obsmeta* meta, model* m, observations* obs);
void reader_smos_standard(char* fname, int fid, obsmeta* meta, model* m, observations* obs);

#define _ALLREADERS_H
#endif
