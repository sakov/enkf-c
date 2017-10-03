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

/*
 * generic readers
 */
void reader_xy_scattered(char* fname, int fid, obsmeta* meta, model* m, observations* obs);
void reader_xy_gridded(char* fname, int fid, obsmeta* meta, model* m, observations* obs);
void reader_xyz_scattered(char* fname, int fid, obsmeta* meta, model* m, observations* obs);
void reader_xyz_gridded(char* fname, int fid, obsmeta* meta, model* m, observations* obs);

/*
 * readers for particular products
 */
void reader_rads_standard(char* fname, int fid, obsmeta* meta, model* m, observations* obs);
void reader_rads_standard2(char* fname, int fid, obsmeta* meta, model* m, observations* obs);
void reader_navo_standard(char* fname, int fid, obsmeta* meta, model* m, observations* obs);
void reader_windsat_standard(char* fname, int fid, obsmeta* meta, model* m, observations* obs);
void reader_cars_standard(char* fname, int fid, obsmeta* meta, model* m, observations* obs);
void reader_mmt_standard(char* fname, int fid, obsmeta* meta, model* m, observations* obs);
void reader_amsr2_standard(char* fname, int fid, obsmeta* meta, model* m, observations* obs);
void reader_amsre_standard(char* fname, int fid, obsmeta* meta, model* m, observations* obs);
void reader_pathfinder_standard(char* fname, int fid, obsmeta* meta, model* m, observations* obs);
void reader_h8_standard(char* fname, int fid, obsmeta* meta, model* m, observations* obs);
void reader_viirs_standard(char* fname, int fid, obsmeta* meta, model* m, observations* obs);

#define _ALLREADERS_H
#endif
