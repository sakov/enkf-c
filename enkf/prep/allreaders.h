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

typedef void (*obsread_fn) (char* fname, int fid, obsmeta* meta, grid* g, observations* obs);

obsread_fn get_obsreadfn(obsmeta* m);

/*
 * generic readers
 */
void reader_scattered(char* fname, int fid, obsmeta* meta, grid* g, observations* obs);
void reader_gridded_xy(char* fname, int fid, obsmeta* meta, grid* g, observations* obs);
void reader_gridded_xyz(char* fname, int fid, obsmeta* meta, grid* g, observations* obs);
void reader_gridded_xyh(char* fname, int fid, obsmeta* meta, grid* g, observations* obs);
void reader_z(char* fname, int fid, obsmeta* meta, grid* g, observations* obs);

/*
 * readers for particular products
 */
void reader_navo(char* fname, int fid, obsmeta* meta, grid* g, observations* obs);
void reader_windsat(char* fname, int fid, obsmeta* meta, grid* g, observations* obs);
void reader_cars(char* fname, int fid, obsmeta* meta, grid* g, observations* obs);
void reader_mmt(char* fname, int fid, obsmeta* meta, grid* g, observations* obs);
void reader_amsr2(char* fname, int fid, obsmeta* meta, grid* g, observations* obs);
void reader_amsre(char* fname, int fid, obsmeta* meta, grid* g, observations* obs);
void reader_cmems(char* fname, int fid, obsmeta* meta, grid* g, observations* obs);

void reader_scattered_describe(void);
void reader_gridded_xy_describe(void);
void reader_gridded_xyz_describe(void);
void reader_gridded_xyh_describe(void);
void reader_z_describe(void);
void reader_navo_describe(void);
void reader_windsat_describe(void);
void reader_cars_describe(void);
void reader_mmt_describe(void);
void reader_amsr2_describe(void);
void reader_amsre_describe(void);
void reader_cmems_describe(void);

#define _ALLREADERS_H
#endif
