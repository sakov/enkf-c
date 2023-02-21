/******************************************************************************
 *
 * File:        ncutils.h        
 *
 * Created:     6/9/2019
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description: Intermediate level NetCDF r/w procedures. Among other things
 *              these procedures are supposed to handle the following
 *              attributes: _FillValue, missing_value, valid_range, valid_min,
 *              valid_max, add_offset, and scale_factor.
 *
 * Revisions:
 *
 *****************************************************************************/

#if !defined(_NCUTILS_H)

typedef void (*ncu_quit_fn) (char* format, ...);

void ncu_set_quitfn(ncu_quit_fn quit_fn);

/*
 * model variable information
 */
int ncu_getnlevels(char fname[], char varname[], int isstructured);
int ncu_getnD(char fname[], char varname[]);

/*
 * generic read procedures
 */
void ncu_readvarfloat(int ncid, int varid, size_t n, float v[]);
void ncu_readvardouble(int ncid, int varid, size_t n, double v[]);

/*
 * model r/w procedures
 */
void ncu_readfield(char fname[], char varname[], int k, int ni, int nj, int nk, float* v);
void ncu_writefield(char fname[], char varname[], int k, int ni, int nj, int nk, float* v);
void ncu_writerow(char fname[], char varname[], int k, int j, int ni, int nj, int nk, float* v);
void ncu_read3dfield(char fname[], char varname[], int ni, int nj, int nk, float* v);

#define _NCUTILS_H
#endif
