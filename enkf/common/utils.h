/******************************************************************************
 *
 * File:        utils.h        
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

#if !defined(_UTILS_H)

#include <unistd.h>

void enkf_init(int* argc, char*** argv);
void enkf_finish(void);
void enkf_quit(char* format, ...);
int enkf_printf(const char* format, ...);
void enkf_flush(void);
FILE* enkf_fopen(const char fname[], const char* mode);
void enkf_printtime(const char offset[]);
void enkf_printcompileflags(const char offset[]);
void enkf_printflags(const char offset[]);
void enkf_printversion(void);
void find_files(char template[], int* nfiles, char*** fnames);
void tunits_convert(char tunits[], double* tunits_multiple, double* tunits_offset);
void file_delete(char* fname);
int file_exists(char fname[]);
int str2double(char token[], double* value);
int str2int(char* token, int* value);
void* alloc2d(int n2, int n1, size_t unitsize);
void free2d(void* pp);
void* alloc3d(int n3, int n2, int n1, size_t unitsize);
void free3d(void* ppp);
int nc_isunlimdimid(int ncid, int dimid);
int getnlevels(char fname[], char varname[]);
void readfield(char fname[], int k, char varname[], float* v);
void writefield(char fname[], int k, char varname[], float* v);
void read3dfield(char fname[], char varname[], float* v);
double date_str2dbl(char strdate[]);
int get_tshift(double date, double tstep);
void print_matrix_double(int n, int m, double** A, char offset[]);
void print_matrix_float(int n, int m, float** A, char offset[]);
void print_vector_double(int n, double* a, char offset[]);
void print_vector_float(int n, float* a, char offset[]);
ENSOBSTYPE interpolate2d(double fi, double fj, int ni, int nj, float** v, int** mask);
ENSOBSTYPE interpolate3d(double fi, double fj, double fk, int ni, int nj, int nk, float*** v, int** nlevels);

#define _UTILS_H
#endif
