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
 * Revisions:   6/9/2019 PS - moved NetCDF procedures to ncutils.h
 *
 *****************************************************************************/

#if !defined(_UTILS_H)

#include <unistd.h>
#include "kdtree.h"

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
int dir_createifabsent(char dirname[]);
void dir_rmifexists(char dirname[]);
int str2double(char token[], double* value);
int str2float(char* token, float* value);
int str2int(char* token, int* value);
int str2bool(char* token, int* value);
void* alloc2d(size_t nj, size_t ni, size_t unitsize);
void* cast2d(void* p, size_t nj, size_t ni, size_t unitsize);
void* copy2d(void** src, size_t nj, size_t ni, size_t unitsize);
void* alloc3d(size_t n3, size_t nj, size_t ni, size_t unitsize);
void* copy3d(void*** src, size_t nk, size_t nj, size_t ni, size_t unitsize);
double date2day(char strdate[]);
int get_tshift(double date, double tstep, int centred);
void print_matrix_double(int n, int m, double** A, char offset[]);
void print_matrix_float(int n, int m, float** A, char offset[]);
void print_vector_double(int n, double* a, char offset[]);
void print_vector_float(int n, float* a, char offset[]);
ENSOBSTYPE interpolate2d(double fi, double fj, int ni, int nj, float** v, int** mask, int periodic_i);
ENSOBSTYPE interpolate3d(double fi, double fj, double fk, int ni, int nj, int nk, int ktop, float*** v, int** nlevels, int periodic_i);
ENSOBSTYPE average2d(int n, size_t* ids, float** v);
int island(double fi, double fj, double fk, int ni, int nj, int ksurf, int** numlevels, int periodic_i);
double taper_gc(double x);
void ll2xyz(double in[2], double out[3]);
void print_commandinfo(int argc, char* argv[]);
void get_normalpair(double x[]);
int istrue(char str[]);
int inloninterval(double lon, double lon1, double lon2);
void shuffle(size_t n, size_t ids[]);
void kd_printinfo(kdtree* tree, char* offset);

void print_memory_usage();

#if defined(ENKF_PREP) || defined(ENKF_CALC)
#if defined(INTERNAL_QSORT_R)
typedef int (*__compar_d_fn_t) (const void*, const void*, void*);
static inline void qsort_r(void* base, size_t n, size_t width, __compar_d_fn_t cmp, void* p)
{
    int nested_cmp(const void* a, const void* b) {
        return cmp(a, b, p);
    }
    qsort(base, n, width, nested_cmp);
}
#endif
#endif

#define _UTILS_H
#endif
