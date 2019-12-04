/******************************************************************************
 *
 * File:        ncutils.h        
 *
 * Created:     6/9/2019
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description: Moved NetCDF procedures from utils.h and prep_utils.h to here
 *
 * Revisions:
 *
 *****************************************************************************/

#if !defined(_NCUTILS_H)

typedef void (*ncu_quit_fn) (char* format, ...);

void ncu_set_quitfn(ncu_quit_fn quit_fn);
int ncu_getnlevels(char fname[], char varname[]);
void ncu_readfield(char fname[], char varname[], int k, int ni, int nj, int nk, float* v);
void ncu_writefield(char fname[], char varname[], int k, int ni, int nj, int nk, float* v);
void ncu_writerow(char fname[], char varname[], int k, int j, float* v);
void ncu_read3dfield(char fname[], char varname[], int ni, int nj, int nk, float* v);
int ncu_getnD(char fname[], char varname[]);
void ncu_readvarfloat(int ncid, int varid, int n, float v[]);
void ncu_readvardouble(int ncid, int varid, int n, double v[]);

#define _NCUTILS_H
#endif
