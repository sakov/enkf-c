/******************************************************************************
 *
 * File:        utils.c        
 *
 * Created:     12/2012
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description: Miscellaneous utilities for EnKF-C package.
 *
 * Revisions:   5/2016 PS - modified alloc2d() and alloc3d(). They now work
 *                with only one internal allocation, and no longer require
 *                special dellocation 
 *
 *****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <errno.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <float.h>
#include <execinfo.h>
#include <dirent.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <stdint.h>
#if !defined(NO_GRIDUTILS)
#include <guquit.h>
#endif
#include "ncw.h"
#include "definitions.h"
#include "version.h"
#include "utils.h"

#define FILE_FIND_INC 10
#define EPS_DOUBLE 4.0e-15
#define EPS_FLOAT  1.0e-7
#define BACKTRACE_SIZE 50

long int seed_rand48 = 0;

/**
 */
void enkf_quit(char* format, ...)
{
    va_list args;

    fflush(stdout);

    fprintf(stderr, "\n\n  ERROR: enkf: CPU #%d: ", rank);
    va_start(args, format);
    vfprintf(stderr, format, args);
    va_end(args);
    fprintf(stderr, "\n\n");
    enkf_printtime("  ");

    if (enkf_exitaction == EXITACTION_BACKTRACE) {
        void* buffer[BACKTRACE_SIZE];
        size_t size;
        char** strings;
        size_t i;

        fprintf(stderr, "\n  I am CPU #%d, now printing the backtrace:\n\n", rank);
        size = backtrace(buffer, BACKTRACE_SIZE);
        strings = backtrace_symbols(buffer, size);
        fprintf(stderr, "  obtained %zd stack frames:\n", size);
        for (i = 0; i < size; i++)
            fprintf(stderr, "%s\n", strings[i]);
        free(strings);
    } else if (enkf_exitaction == EXITACTION_SEGFAULT) {
        fprintf(stderr, "\n  I am CPU #%d, now generating a segfault:\n\n", rank);
    }
    fflush(NULL);               /* flush all streams */
#if defined(MPI)
    MPI_Abort(MPI_COMM_WORLD, 1);       /* kill all MPI jobs */
#else
    abort();                    /* raise SIGABRT for debugging */
#endif

    exit(1);
}

/**
 */
static void randomise_rand48(void)
{
    char fname[] = "/dev/urandom";
    FILE* f;
    size_t status;

    if (seed_rand48 != 0)
        return;

    f = enkf_fopen(fname, "r");

    status = fread(&seed_rand48, sizeof(seed_rand48), 1, f);
    if (status != 1) {
        int errno_saved = errno;

        enkf_quit("randomise_rand48(): could not read from \"%s\": %s", fname, strerror(errno_saved));
    }
    fclose(f);
    srand(seed_rand48);
}

/**
 */
void enkf_init(int* argc, char*** argv)
{
#if defined(MPI)
#include <unistd.h>
    int size;

    MPI_Init(argc, argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocesses);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (*argc > 1 && nprocesses > 1) {
        enkf_printf("  MPI: initialised %d process(es)\n", nprocesses);
        MPI_Barrier(MPI_COMM_WORLD);
        if (enkf_verbose) {
            printf("  MPI: rank = %d, PID = %d\n", rank, getpid());
            fflush(NULL);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Type_size(MPI_INTEGER, &size);
    assert(size == sizeof(int));
    MPI_Type_size(MPIENSOBSTYPE, &size);
    assert(size == sizeof(ENSOBSTYPE));
#endif

    ncw_set_quitfn(enkf_quit);
#if !defined(NO_GRIDUTILS)
    gu_setquitfn(enkf_quit);
#endif

    /*
     * initialise the random number generator to a random state for each cpu
     */
    randomise_rand48();
}

/**
 */
void enkf_finish(void)
{
#if defined(MPI)
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    enkf_printtime("  ");
    enkf_printf("  finished\n");
    enkf_flush();
#if defined(MPI)
    MPI_Finalize();
#endif
}

/**
 */
int enkf_printf(const char* format, ...)
{
    va_list args;
    int status = 0;

    if (enkf_verbose && rank == 0) {
        va_start(args, format);
        status = vprintf(format, args);
        va_end(args);
    }

    return status;
}

/**
 */
void enkf_flush(void)
{
    if (rank == 0)
        fflush(stdout);
}

/**
 */
FILE* enkf_fopen(const char* fname, const char* mode)
{
    FILE* f = NULL;

    if ((f = fopen(fname, mode)) == NULL) {
        int errno_saved = errno;

        enkf_quit("enkf_fopen(): could not open \"%s\": %s", fname, strerror(errno_saved));
    }

    return f;
}

/**
 */
void enkf_printtime(const char offset[])
{
    time_t t;
    struct tm tm;

    if (rank != 0 || !enkf_verbose)
        return;

    t = time(NULL);
    tm = *localtime(&t);

    printf("%s%d-%02d-%02d %02d:%02d:%02d\n", offset, tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
}

/**
 */
void enkf_printcompileflags(const char offset[])
{
    enkf_printf("%scompile flags:\n", offset);
#if defined(ENKF_CALC)
#if defined(CHECK_X5)
    enkf_printf("%s  CHECK_X5         = [+]\n", offset);
#else
    enkf_printf("%s  CHECK_X5         = [-]\n", offset);
#endif
#if defined(CHECK_G)
    enkf_printf("%s  CHECK_G          = [+]\n", offset);
#else
    enkf_printf("%s  CHECK_G          = [-]\n", offset);
#endif
#if defined(SHUFFLE_ROWS)
    enkf_printf("%s  SHUFFLE_ROWS     = [+]\n", offset);
#else
    enkf_printf("%s  SHUFFLE_ROWS     = [-]\n", offset);
#endif
#if defined(HE_VIAFILE)
    enkf_printf("%s  HE_VIAFILE       = [+]\n", offset);
#else
    enkf_printf("%s  HE_VIAFILE       = [-]\n", offset);
#endif
#if defined(HE_VIASHMEM)
    enkf_printf("%s  HE_VIASHMEM      = [+]\n", offset);
#else
    enkf_printf("%s  HE_VIASHMEM      = [-]\n", offset);
#endif
#if defined(MINIMISE_ALLOC)
    enkf_printf("%s  MINIMISE_ALLOC   = [+]\n", offset);
#else
    enkf_printf("%s  MINIMISE_ALLOC   = [-]\n", offset);
#endif
#if defined(OBS_SHUFFLE)
    enkf_printf("%s  OBS_SHUFFLE      = [+]\n", offset);
#else
    enkf_printf("%s  OBS_SHUFFLE      = [-]\n", offset);
#endif
#endif                          /* ENKF_CALC */
#if defined(ENKF_PREP) || defined(ENKF_CALC)
#if defined(GRIDNODES_WRITE)
    enkf_printf("%s  GRIDNODES_WRITE  = [+]\n", offset);
#else
    enkf_printf("%s  GRIDNODES_WRITE  = [-]\n", offset);
#endif
#if defined(INTERNAL_QSORT_R)
    enkf_printf("%s  INTERNAL_QSORT_R = [+]\n", offset);
#else
    enkf_printf("%s  INTERNAL_QSORT_R = [-]\n", offset);
#endif
#if defined(NO_GRIDUTILS)
    enkf_printf("%s  NO_GRIDUTILS     = [+]\n", offset);
#else
    enkf_printf("%s  NO_GRIDUTILS     = [-]\n", offset);
#endif
#endif                          /* ENKF_PREP || ENKF_CALC */
#if defined(ENKF_UPDATE)
#if defined(NCW_SKIPSINGLE)
    enkf_printf("%s  NCW_SKIPSINGLE = [+]\n", offset);
#else
    enkf_printf("%s  NCW_SKIPSINGLE = [-]\n", offset);
#endif
#endif                          /* ENKF_UPDATE */
}

/**
 */
void enkf_printflags(const char offset[])
{
    enkf_printf("%sEnKF flags:\n", offset);
    if (enkf_exitaction == EXITACTION_BACKTRACE)
        enkf_printf("%s  enkf_exitaction  = [BACKTRACE]\n", offset);
    else if (enkf_exitaction == EXITACTION_SEGFAULT)
        enkf_printf("%s  enkf_exitaction  = [SEGFAULT]\n", offset);
    if (enkf_obstype == OBSTYPE_VALUE)
        enkf_printf("%s  enkf_obstype     = [VALUE]\n", offset);
    else if (enkf_obstype == OBSTYPE_INNOVATION)
        enkf_printf("%s  enkf_obstype     = [INNOVATION]\n", offset);
    enkf_printcompileflags(offset);
}

/**
 */
void enkf_printversion(void)
{
    enkf_printf("  EnKF version %s\n", ENKF_VERSION);
    enkf_printcompileflags("  ");
}

#if 1
#include <glob.h>
/** Find files matching the template using glob.
 */
void find_files(char* template, int* nfiles, char*** fnames)
{
    glob_t gl;
    int status;

    status = glob(template, GLOB_PERIOD | GLOB_TILDE_CHECK, NULL, &gl);
    if (status == GLOB_NOSPACE || status == GLOB_ABORTED || status == GLOB_ERR) {
        int errno_saved = errno;

        enkf_quit("failed looking for \"%s\": %s", template, strerror(errno_saved));
    }

    if (gl.gl_pathc > 0) {
        int i, ii;

        if (*nfiles == 0)
            *fnames = NULL;
        *fnames = realloc(*fnames, (*nfiles + gl.gl_pathc) * sizeof(void*));
        for (i = 0, ii = *nfiles; i < gl.gl_pathc; ++i, ++ii)
            (*fnames)[ii] = strdup(gl.gl_pathv[i]);
        *nfiles += gl.gl_pathc;
        globfree(&gl);
    } else if (strncasecmp(template, "http", 4) == 0) {
        if (*nfiles == 0)
            *fnames = NULL;
        *fnames = realloc(*fnames, (*nfiles + 1) * sizeof(void*));
        (*fnames)[*nfiles] = strdup(template);
        *nfiles += 1;
    }
}
#else
/** Find files matching the template using "ls".
 */
void find_files(char* template, int* nfiles, char*** fnames)
{
    char command[MAXSTRLEN];

    FILE* in;
    char buf[MAXSTRLEN];
    char* eol;

    snprintf(command, MAXSTRLEN, "ls -1 %s", template);
    fflush(stdout);
    in = popen(command, "r");
    if (in == NULL)
        return;

    while (fgets(buf, MAXSTRLEN, in) != NULL) {
        if (*nfiles % FILE_FIND_INC == 0)
            *fnames = realloc(*fnames, (*nfiles + FILE_FIND_INC) * sizeof(char*));
        eol = strchr(buf, '\n');
        if (eol != NULL)
            *eol = 0;
        (*fnames)[*nfiles] = strdup(buf);
        (*nfiles)++;
    }

    pclose(in);
}
#endif

/*
** scalar date routines    --    public domain by Ray Gardner
** Numerically, these will work over the range 1/01/01 thru 14699/12/31.
** Practically, these only work from the beginning of the Gregorian 
** calendar thru 14699/12/31.  The Gregorian calendar took effect in
** much of Europe in about 1582, some parts of Germany in about 1700, in
** England and the colonies in about 1752ff, and in Russia in 1918.
*/
static int isleap(unsigned yr)
{
    return yr % 400 == 0 || (yr % 4 == 0 && yr % 100 != 0);
}
static unsigned months_to_days(unsigned month)
{
    return (month * 3057 - 3007) / 100;
}
static long years_to_days(unsigned yr)
{
    return yr * 365L + yr / 4 - yr / 100 + yr / 400;
}
static long ymd_to_scalar(unsigned yr, unsigned mo, unsigned day)
{
    long scalar;

    scalar = day + months_to_days(mo);
    if (mo > 2)                 /* adjust if past February */
        scalar -= isleap(yr) ? 1 : 2;
    yr--;
    scalar += years_to_days(yr);
    return scalar;
}

/**
 */
static long int daydiff(unsigned int y1, unsigned int m1, unsigned int d1, unsigned int y2, unsigned int m2, unsigned int d2)
{
    long int dn1, dn2;

    dn1 = ymd_to_scalar(y1, m1, d1);
    dn2 = ymd_to_scalar(y2, m2, d2);

    return dn1 - dn2;
}

/** Calculates transform from `tunits' to "days from BASEDAY/BASEMONTH/BASYEAR".
 * This is a very basic procedure, do not expect it to work smoothly if things
 * get complicated. The conversion equation is
 *
 * <new values> = <old values> * tunits_multiple + tunits_offset
 *
 * @param tunits - input, e.g. "seconds since 1981-01-01 00:00:00"
 * @param tunits_multiple - output, multiple
 * @param tunits_offset - output, offset
 */
void tunits_convert(char* tunits, double* tunits_multiple, double* tunits_offset)
{
    char* startdate;
    int year, month, day, h, m, s;
    char* token;
    char seps_date[] = " -\n";
    char seps_time[] = " :\n";

    if (strncasecmp(tunits, "fraction of a ", 14) == 0)
        tunits += 14;

    if (strncasecmp(tunits, "sec", 3) == 0)
        *tunits_multiple = 1.0 / 86400.0;
    else if (strncasecmp(tunits, "hou", 3) == 0)
        *tunits_multiple = 1.0 / 24.0;
    else if (strncasecmp(tunits, "day", 3) == 0)
        *tunits_multiple = 1.0;
    else
        enkf_quit("can not interpret time units \"%s\"", tunits);

    startdate = strstr(tunits, "since");
    if (startdate == NULL)
        enkf_quit("can not interpret time units \"%s\"", tunits);
    startdate += strlen("since ");
    if ((token = strtok(startdate, seps_date)) == NULL)
        enkf_quit("can not interpret time units \"%s\"", tunits);
    if (!str2int(token, &year))
        enkf_quit("could not convert \"%s\" to time units", tunits);
    if ((token = strtok(NULL, seps_date)) == NULL)
        enkf_quit("can not interpret time units \"%s\"", tunits);
    if (!str2int(token, &month))
        enkf_quit("could not convert \"%s\" to time units", tunits);
    if ((token = strtok(NULL, seps_date)) == NULL)
        enkf_quit("can not interpret time units \"%s\"", tunits);
    if (!str2int(token, &day))
        enkf_quit("could not convert \"%s\" to time units", tunits);
    h = 0;
    m = 0;
    s = 0;
    if ((token = strtok(NULL, seps_time)) != NULL)
        if (!str2int(token, &h))
            enkf_quit("could not convert \"%s\" to time units", tunits);
    if ((token = strtok(NULL, seps_time)) != NULL)
        if (!str2int(token, &m))
            enkf_quit("could not convert \"%s\" to time units", tunits);
    if ((token = strtok(NULL, seps_time)) != NULL)
        if (!str2int(token, &s))
            enkf_quit("could not convert \"%s\" to time units", tunits);

    *tunits_offset = (double) daydiff(year, month, day, BASEYEAR, BASEMONTH, BASEDAY) - (double) h / 24.0 - (double) m / 1440.0 - (double) s / 86400.0;
}

/**
 */
int file_exists(char* fname)
{
    FILE* f;

    f = fopen(fname, "r");
    if (f == NULL)
        return 0;
    fclose(f);
    return 1;
}

/**
 */
void file_delete(char* fname)
{
    int status = -1;

    status = unlink(fname);
    if (status != 0) {
        int errno_saved = errno;

        enkf_quit("could not delete file \"%s\": %s", fname, strerror(errno_saved));
    }
}

/**
 */
static int dir_exists(char dirname[])
{
    DIR* d;

    d = opendir(dirname);
    if (d == NULL)
        return 0;
    closedir(d);
    return 1;
}

/**
 */
int dir_createifabsent(char dirname[])
{
    int status;

    if (dir_exists(dirname))
        return 1;
    status = mkdir(dirname, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    if (status == 0)
        return 1;
    else {
        int errno_saved = errno;

        enkf_quit("could not create directory \"%s\": %s", dirname, strerror(errno_saved));
    }
    return 0;
}

/**
 */
void dir_rmifexists(char dirname[])
{
    int status;

    if (!dir_exists(dirname))
        return;
    status = rmdir(dirname);
    if (status != 0) {
        int errno_saved = errno;

        enkf_quit("could not remove directory \"%s\": %s", dirname, strerror(errno_saved));
    }
}

/**
 */
int str2double(char* token, double* value)
{
    char* end = NULL;

    if (token == NULL) {
        *value = NAN;
        return 0;
    }

    *value = strtod(token, &end);

    if (end == token) {
        *value = NAN;
        return 0;
    }

    return 1;
}

/**
 */
int str2float(char* token, float* value)
{
    char* end = NULL;

    if (token == NULL) {
        *value = NAN;
        return 0;
    }

    *value = strtod(token, &end);

    if (end == token) {
        *value = NAN;
        return 0;
    }

    return 1;
}

/**
 */
int str2int(char* token, int* value)
{
    long int tmp;
    char* end = NULL;

    if (token == NULL) {
        *value = INT_MAX;
        return 0;
    }

    tmp = strtol(token, &end, 10);

    if (end == token || tmp > INT_MAX || tmp < INT_MIN) {
        *value = INT_MAX;
        return 0;
    }

    *value = (int) tmp;
    return 1;
}

/**
 */
int str2bool(char* token, int* value)
{
    if (token[0] == 'y' || token[0] == 'Y')
        *value = 1;
    else if (token[0] == 'n' || token[0] == 'N')
        *value = 0;
    else if (!str2int(token, value))
        return 0;

    if (*value == 0 || *value == 1)
        return 1;
    return 0;
}

/** Allocates ni x nj matrix of something and fills it with zeros. An element
 ** (i,j) will be accessed as [j][i]. For deallocation use free().
 *
 * @param nj Dimension 2
 * @param ni Dimension 1
 * @param unitsize Size of one matrix element in bytes
 * @return Matrix
 */
void* alloc2d(size_t nj, size_t ni, size_t unitsize)
{
    size_t size;
    void* p;
    void** pp;
    int i;

    if (ni <= 0 || nj <= 0)
        enkf_quit("alloc2d(): invalid size (nj = %d, ni = %d)", nj, ni);

    size = nj * sizeof(void*) + nj * ni * unitsize;
    if ((p = malloc(size)) == NULL) {
        int errno_saved = errno;

        enkf_quit("alloc2d(): %s", strerror(errno_saved));
    }
    memset(p, 0, size);

    pp = p;
    p = &((size_t*) p)[nj];
    for (i = 0; i < nj; ++i)
        pp[i] = &((char*) p)[i * ni * unitsize];

    return pp;
}

/** Copies 2D matrix.
 * @param src Source matrix
 * @param nj Dimension 2
 * @param ni Dimension 1
 * @param unitsize Size of one matrix element in bytes
 * @return Copy matrix
 */
void* copy2d(void** src, size_t nj, size_t ni, size_t unitsize)
{
    size_t size;
    void* p;
    void** pp;
    size_t i;

    if (ni <= 0 || nj <= 0)
        enkf_quit("copy2d(): invalid size (nj = %d, ni = %d)", nj, ni);

    size = nj * ni * unitsize + nj * sizeof(void*);
    if ((p = malloc(size)) == NULL) {
        int errno_saved = errno;

        enkf_quit("copy2d(): %s", strerror(errno_saved));
    }

    pp = p;
    p = &((size_t*) p)[nj];
    for (i = 0; i < nj; ++i)
        pp[i] = &((char*) p)[i * ni * unitsize];
    memcpy(p, src[0], nj * ni * unitsize);

    return pp;
}

/** Allocates ni x nj x nk matrix of something and fills it with zeros. An
 * element (i,j,k) will be accessed as [k][j][i]. For deallocation use free().
 *
 * @param nk Dimension 3
 * @param nj Dimension 2
 * @param ni Dimension 1
 * @param unitsize Size of one matrix element in bytes
 * @return Matrix
 */
void* alloc3d(size_t nk, size_t nj, size_t ni, size_t unitsize)
{
    size_t size;
    void* p;
    void** pp;
    void*** ppp;
    size_t i;

    if (nk <= 0 || nj <= 0 || ni <= 0)
        enkf_quit("alloc3d(): invalid size (nk = %d, nj = %d, ni = %d)", nk, nj, ni);

    size = nk * (nj + 1) * sizeof(void*) + nk * nj * ni * unitsize;
    if ((p = malloc(size)) == NULL) {
        int errno_saved = errno;

        enkf_quit("alloc3d(): %s", strerror(errno_saved));
    }
    memset(p, 0, size);

    ppp = p;
    pp = &((void**) p)[nk];
    p = &((size_t*) p)[nk + nk * nj];
    for (i = 0; i < nk; ++i)
        ppp[i] = &pp[i * nj];
    for (i = 0; i < nk * nj; ++i)
        pp[i] = &((char*) p)[i * ni * unitsize];

    return ppp;
}

/** Copies nk x nj x ni array of something.
 * @param src Source matrix
 * @param nk Dimension 3
 * @param nj Dimension 2
 * @param ni Dimension 1
 * @param unitsize Size of one matrix element in bytes
 * @return Copy matrix
 */
void* copy3d(void*** src, size_t nk, size_t nj, size_t ni, size_t unitsize)
{
    size_t size;
    void* p;
    void** pp;
    void*** ppp;
    int i;

    if (nk <= 0 || nj <= 0 || ni <= 0)
        enkf_quit("copy3d(): invalid size (nk = %d, nj = %d, ni = %d)", nk, nj, ni);

    size = nk * (nj + 1) * sizeof(void*) + nk * nj * ni * unitsize;
    if ((p = malloc(size)) == NULL) {
        int errno_saved = errno;

        enkf_quit("copy3d(): %s", strerror(errno_saved));
    }

    ppp = p;
    pp = &((void**) p)[nk];
    p = &((size_t*) p)[nk + nk * nj];
    for (i = 0; i < nk; i++)
        ppp[i] = &pp[i * nj];
    for (i = 0; i < nk * nj; i++)
        pp[i] = &((char*) p)[i * ni * unitsize];
    memcpy(p, src[0][0], nk * nj * unitsize);

    return ppp;
}

/**
 */
int getnlevels(char fname[], char varname[])
{
    int ncid;
    int varid;
    int ndims;
    size_t dimlen[4];
    int hasrecorddim;

    ncw_open(fname, NC_NOWRITE, &ncid);
    ncw_inq_varid(ncid, varname, &varid);
    ncw_inq_vardims(ncid, varid, 4, &ndims, dimlen);
    hasrecorddim = ncw_var_hasunlimdim(ncid, varid);
    ncw_close(ncid);

    if (ndims > 4)
        enkf_quit("%s: %s: EnKF-C does not know how to handle more than 4-dimensional variables\n", fname, varname);
    if (ndims == 4) {
        if (!hasrecorddim)
            enkf_quit("%s: %s: expect an unlimited dimension to be present for a 4-dimensional variable\n", fname, varname);
        return (int) dimlen[1];
    }
    if (ndims == 3)
        return (hasrecorddim) ? 1 : (int) dimlen[0];
    if (ndims == 2)
        return (hasrecorddim) ? 0 : 1;

    return 0;
}

/** Reads one horizontal field (layer) for a variable from a NetCDF file.
 ** Verifies that the field dimensions are ni x nj.
 */
void readfield(char fname[], char varname[], int k, int ni, int nj, int nk, float* v)
{
    int ncid;
    int varid;
    int ndims;
    size_t dimlen[4];
    size_t start[4], count[4];
    int i, n;
    int hasrecorddim;

    ncw_open(fname, NC_NOWRITE, &ncid);
    ncw_inq_varid(ncid, varname, &varid);
    ncw_inq_vardims(ncid, varid, 4, &ndims, dimlen);
    hasrecorddim = ncw_var_hasunlimdim(ncid, varid);
    if (hasrecorddim)
        if (dimlen[0] == 0)
            enkf_quit("%s: %s: empty record dimension", fname, varname);

    if (ndims == 4) {
        if (!hasrecorddim)
            enkf_quit("%s: %s: expect an unlimited dimension to be present for a 4-dimensional variable\n", fname, varname);
        start[0] = dimlen[0] - 1;
        if (dimlen[1] != nk) {
            if (dimlen[1] != 1)
                enkf_quit("\"%s\": vertical dimension of variable \"%s\" (nk = %d) does not match grid dimension (nk = %d)", fname, varname, dimlen[1], nk);
            else
                k = 0;          /* ignore k */
        }
        start[1] = k;
        start[2] = 0;
        start[3] = 0;
        count[0] = 1;
        count[1] = 1;
        count[2] = dimlen[2];
        count[3] = dimlen[3];
        if (dimlen[3] != ni || dimlen[2] != nj)
            enkf_quit("\"%s\": horizontal dimensions of variable \"%s\" (ni = %d, nj = %d) do not match grid dimensions (ni = %d, nj = %d)", fname, varname, dimlen[3], dimlen[2], ni, nj);
    } else if (ndims == 3) {
        if (!hasrecorddim) {
            if (dimlen[0] != nk && !(dimlen[0] == 1 && (k == 0 || k == nk - 1)))
                enkf_quit("\"%s\": vertical dimension of variable \"%s\" (nk = %d) does not match grid dimension (nk = %d)", fname, varname, dimlen[0], nk);
            start[0] = k;
            start[1] = 0;
            start[2] = 0;
            count[0] = 1;
            count[1] = dimlen[1];
            count[2] = dimlen[2];
        } else {
            /*
             * ignore k in this case
             */
            start[0] = dimlen[0] - 1;
            start[1] = 0;
            start[2] = 0;
            count[0] = 1;
            count[1] = dimlen[1];
            count[2] = dimlen[2];
        }
        if (dimlen[2] != ni || dimlen[1] != nj)
            enkf_quit("\"%s\": horizontal dimensions of variable \"%s\" (ni = %d, nj = %d) do not match grid dimensions (ni = %d, nj = %d)", fname, varname, dimlen[2], dimlen[1], ni, nj);
    } else if (ndims == 2) {
        if (hasrecorddim)
            enkf_quit("%s: can not read a layer from a 1D variable \"%s\"", fname, varname);
        /*
         * ignore k
         */
        start[0] = 0;
        start[1] = 0;
        count[0] = dimlen[0];
        count[1] = dimlen[1];
        if (dimlen[1] != ni || dimlen[0] != nj)
            enkf_quit("\"%s\": horizontal dimensions of variable \"%s\" (ni = %d, nj = %d) do not match grid dimensions (ni = %d, nj = %d)", fname, varname, dimlen[1], dimlen[0], ni, nj);
    } else
        enkf_quit("%s: can not read 2D field for \"%s\": # of dimensions = %d", fname, varname, ndims);

    ncw_get_vara_float(ncid, varid, start, count, v);

    n = 1;
    for (i = 0; i < ndims; ++i)
        n *= count[i];

    if (ncw_att_exists(ncid, varid, "valid_range")) {
        nc_type xtype;
        size_t len;

        ncw_inq_att(ncid, varid, "valid_range", &xtype, &len);
        if (len != 2)
            enkf_quit("%s: %s: \"valid_range\" attribute must have 2 elements", fname, varname);
        if (xtype == NC_SHORT) {
            int valid_range[2];
            int fill_value = SHRT_MIN;
            short* vv = malloc(n * sizeof(short));

            ncw_get_vara_short(ncid, varid, start, count, vv);

            if (ncw_att_exists(ncid, varid, "_FillValue"))
                ncw_get_att_int(ncid, varid, "_FillValue", &fill_value);
            else if (ncw_att_exists(ncid, varid, "missing_value"))
                ncw_get_att_int(ncid, varid, "missing_value", &fill_value);

            ncw_get_att_int(ncid, varid, "valid_range", valid_range);
            for (i = 0; i < n; ++i)
                if (vv[i] == fill_value || vv[i] < valid_range[0] || vv[i] > valid_range[1])
                    v[i] = NAN;
            free(vv);
        } else if (xtype == NC_FLOAT || xtype == NC_DOUBLE) {
            float valid_range[2];
            float fill_value = -FLT_MAX;

            /*
             * (we assume that double numbers outside of [-FLT_MAX, FLT_MAX]
             * interval are converted to either -FLT_MAX or FLT_MAX)
             */

            if (ncw_att_exists(ncid, varid, "_FillValue"))
                ncw_get_att_float(ncid, varid, "_FillValue", &fill_value);
            else if (ncw_att_exists(ncid, varid, "missing_value"))
                ncw_get_att_float(ncid, varid, "missing_value", &fill_value);

            ncw_get_att_float(ncid, varid, "valid_range", valid_range);
            for (i = 0; i < n; ++i)
                if (v[i] == fill_value || v[i] < valid_range[0] || v[i] > valid_range[1])
                    v[i] = NAN;
        } else
            enkf_quit("%s: %s: can not properly read NetCDF types other than NC_SHORT, NC_FLOAT, or NC_DOUBLE yet", fname, varname);
    }

    if (ncw_att_exists(ncid, varid, "scale_factor")) {
        float scale_factor;

        ncw_get_att_float(ncid, varid, "scale_factor", &scale_factor);
        for (i = 0; i < n; ++i)
            v[i] *= scale_factor;
    }

    if (ncw_att_exists(ncid, varid, "add_offset")) {
        float add_offset;

        ncw_get_att_float(ncid, varid, "add_offset", &add_offset);

        for (i = 0; i < n; ++i)
            v[i] += add_offset;
    }

    ncw_close(ncid);
}

/** Writes one horizontal field (layer) for a variable to a NetCDF file.
 */
void writefield(char fname[], char varname[], int k, int ni, int nj, int nk, float* v)
{
    int ncid;
    int varid;
    int ndims;
    size_t dimlen[4];
    size_t start[4], count[4];
    int i, n;
    int hasrecorddim;

    ncw_open(fname, NC_WRITE, &ncid);
    ncw_inq_varid(ncid, varname, &varid);
    ncw_inq_vardims(ncid, varid, 4, &ndims, dimlen);
    hasrecorddim = ncw_var_hasunlimdim(ncid, varid);

    if (ndims == 4) {
        if (!hasrecorddim)
            enkf_quit("%s: %s: expect an unlimited dimension to be present for a 4-dimensional variable\n", fname, varname);
        start[0] = (dimlen[0] == 0) ? 0 : dimlen[0] - 1;
        if (dimlen[1] != nk)
            enkf_quit("\"%s\": vertical dimension of variable \"%s\" (nk = %d) does not match grid dimension (nk = %d)", fname, varname, dimlen[1], nk);
        start[1] = k;
        start[2] = 0;
        start[3] = 0;
        count[0] = 1;
        count[1] = 1;
        count[2] = dimlen[2];
        count[3] = dimlen[3];
        if (dimlen[3] != ni || dimlen[2] != nj)
            enkf_quit("\"%s\": horizontal dimensions of variable \"%s\" (ni = %d, nj = %d) do not match grid dimensions (ni = %d, nj = %d)", fname, varname, dimlen[3], dimlen[2], ni, nj);
    } else if (ndims == 3) {
        if (!hasrecorddim) {
            if (dimlen[0] != nk && !(dimlen[0] == 1 && (k == 0 || k == nk - 1)))
                enkf_quit("\"%s\": vertical dimension of variable \"%s\" (nk = %d) does not match grid dimension (nk = %d)", fname, varname, dimlen[0], nk);
            start[0] = k;
            start[1] = 0;
            start[2] = 0;
            count[0] = 1;
            count[1] = dimlen[1];
            count[2] = dimlen[2];
        } else {
            /*
             * 2D variable, ignore k
             */
            start[0] = (dimlen[0] == 0) ? 0 : dimlen[0] - 1;
            start[1] = 0;
            start[2] = 0;
            count[0] = 1;
            count[1] = dimlen[1];
            count[2] = dimlen[2];
        }
        if (dimlen[2] != ni || dimlen[1] != nj)
            enkf_quit("\"%s\": horizontal dimensions of variable \"%s\" (ni = %d, nj = %d) do not match grid dimensions (ni = %d, nj = %d)", fname, varname, dimlen[2], dimlen[1], ni, nj);
    } else if (ndims == 2) {
        if (hasrecorddim)
            enkf_quit("%s: can not write a layer from a 1D variable \"%s\"", fname, varname);
        /*
         * ignore k
         */
        start[0] = 0;
        start[1] = 0;
        count[0] = dimlen[0];
        count[1] = dimlen[1];
        if (dimlen[1] != ni || dimlen[0] != nj)
            enkf_quit("\"%s\": horizontal dimensions of variable \"%s\" (ni = %d, nj = %d) do not match grid dimensions (ni = %d, nj = %d)", fname, varname, dimlen[1], dimlen[0], ni, nj);
    } else
        enkf_quit("%s: can not write 2D field for \"%s\": # of dimensions = %d", fname, varname, ndims);

    n = 1;
    for (i = 0; i < ndims; ++i)
        n *= count[i];

    if (ncw_att_exists(ncid, varid, "add_offset")) {
        float add_offset;

        ncw_get_att_float(ncid, varid, "add_offset", &add_offset);

        for (i = 0; i < n; ++i)
            v[i] -= add_offset;
    }

    if (ncw_att_exists(ncid, varid, "scale_factor")) {
        float scale_factor;

        ncw_get_att_float(ncid, varid, "scale_factor", &scale_factor);
        for (i = 0; i < n; ++i)
            v[i] /= scale_factor;
    }

    if (ncw_att_exists(ncid, varid, "valid_range")) {
        nc_type xtype;
        size_t len;

        ncw_inq_att(ncid, varid, "valid_range", &xtype, &len);
        if (xtype == NC_SHORT) {
            int valid_range[2];
            int fill_value = SHRT_MIN;

            if (ncw_att_exists(ncid, varid, "_FillValue"))
                ncw_get_att_int(ncid, varid, "_FillValue", &fill_value);
            else if (ncw_att_exists(ncid, varid, "missing_value"))
                ncw_get_att_int(ncid, varid, "missing_value", &fill_value);

            ncw_get_att_int(ncid, varid, "valid_range", valid_range);
            for (i = 0; i < n; ++i) {
                if (v[i] == fill_value)
                    continue;
                if ((int) v[i] < valid_range[0])
                    v[i] = (float) valid_range[0];
                else if ((int) v[i] > valid_range[1])
                    v[i] = (float) valid_range[1];
            }
        } else if (xtype == NC_FLOAT) {
            float valid_range[2];
            float fill_value = -FLT_MAX;

            if (ncw_att_exists(ncid, varid, "_FillValue"))
                ncw_get_att_float(ncid, varid, "_FillValue", &fill_value);
            else if (ncw_att_exists(ncid, varid, "missing_value"))
                ncw_get_att_float(ncid, varid, "missing_value", &fill_value);

            ncw_get_att_float(ncid, varid, "valid_range", valid_range);
            for (i = 0; i < n; ++i) {
                if (v[i] == fill_value)
                    continue;
                if (v[i] < valid_range[0])
                    v[i] = valid_range[0];
                else if (v[i] > valid_range[1])
                    v[i] = valid_range[1];
            }
        } else
            enkf_quit("%s: %s: output for types other than NC_SHORT and NC_FLOAT is not handled", fname, varname);
    }

    ncw_put_vara_float(ncid, varid, start, count, v);
    ncw_close(ncid);
}

/** Writes one row of a horizontal field (layer) for a variable to a NetCDF
 *  file.
 */
void writerow(char fname[], char varname[], int k, int j, float* v)
{
    int ncid;
    int varid;
    int ndims;
    size_t dimlen[4];
    size_t start[4], count[4];
    int i, n;
    int hasrecorddim;

    ncw_open(fname, NC_WRITE, &ncid);
    ncw_inq_varid(ncid, varname, &varid);
    ncw_inq_vardims(ncid, varid, 4, &ndims, dimlen);
    hasrecorddim = ncw_var_hasunlimdim(ncid, varid);

    if (ndims == 4) {
        if (!hasrecorddim)
            enkf_quit("%s: %s: expect an unlimited dimension to be present for a 4-dimensional variable\n", fname, varname);
        if (k >= dimlen[1])
            enkf_quit("%s: %s: the length of dimension 1 (%d) is not sufficient to read layer %d", fname, varname, dimlen[1], k);
        start[0] = (dimlen[0] == 0) ? 0 : dimlen[0] - 1;
        start[1] = k;
        start[2] = j;
        start[3] = 0;
        count[0] = 1;
        count[1] = 1;
        count[2] = 1;
        count[3] = dimlen[3];
    } else if (ndims == 3) {
        if (!hasrecorddim) {
            assert(k < dimlen[0]);
            start[0] = k;
            start[1] = j;
            start[2] = 0;
            count[0] = 1;
            count[1] = 1;
            count[2] = dimlen[2];
        } else {
            /*
             * 2D variable, ignore k
             */
            start[0] = (dimlen[0] == 0) ? 0 : dimlen[0] - 1;
            start[1] = j;
            start[2] = 0;
            count[0] = 1;
            count[1] = 1;
            count[2] = dimlen[2];
        }
    } else if (ndims == 2) {
        if (hasrecorddim)
            enkf_quit("%s: can not write a row %d to a 1D variable \"%s\"", fname, j, varname);
        if (k > 0)
            enkf_quit("%s: can not write layer %d row %d to a 2D variable \"%s\"", fname, k, j, varname);
        start[0] = j;
        start[1] = 0;
        count[0] = 1;
        count[1] = dimlen[1];
    } else
        enkf_quit("%s: can not write a row to 2D field for \"%s\": # of dimensions = %d", fname, varname, ndims);

    n = count[ndims - 1];

    if (ncw_att_exists(ncid, varid, "add_offset")) {
        float add_offset;

        ncw_get_att_float(ncid, varid, "add_offset", &add_offset);

        for (i = 0; i < n; ++i)
            v[i] -= add_offset;
    }

    if (ncw_att_exists(ncid, varid, "scale_factor")) {
        float scale_factor;

        ncw_get_att_float(ncid, varid, "scale_factor", &scale_factor);
        for (i = 0; i < n; ++i)
            v[i] /= scale_factor;
    }

    if (ncw_att_exists(ncid, varid, "valid_range")) {
        nc_type xtype;
        size_t len;

        ncw_inq_att(ncid, varid, "valid_range", &xtype, &len);
        if (xtype == NC_SHORT) {
            int valid_range[2];
            int fill_value = SHRT_MIN;

            if (ncw_att_exists(ncid, varid, "_FillValue"))
                ncw_get_att_int(ncid, varid, "_FillValue", &fill_value);
            else if (ncw_att_exists(ncid, varid, "missing_value"))
                ncw_get_att_int(ncid, varid, "missing_value", &fill_value);

            ncw_get_att_int(ncid, varid, "valid_range", valid_range);
            for (i = 0; i < n; ++i) {
                if (v[i] == fill_value)
                    continue;
                if ((int) v[i] < valid_range[0])
                    v[i] = (float) valid_range[0];
                else if ((int) v[i] > valid_range[1])
                    v[i] = (float) valid_range[1];
            }
        } else if (xtype == NC_FLOAT) {
            float valid_range[2];
            float fill_value = -FLT_MAX;

            if (ncw_att_exists(ncid, varid, "_FillValue"))
                ncw_get_att_float(ncid, varid, "_FillValue", &fill_value);
            else if (ncw_att_exists(ncid, varid, "missing_value"))
                ncw_get_att_float(ncid, varid, "missing_value", &fill_value);

            ncw_get_att_float(ncid, varid, "valid_range", valid_range);
            for (i = 0; i < n; ++i) {
                if (v[i] == fill_value)
                    continue;
                if (v[i] < valid_range[0])
                    v[i] = valid_range[0];
                else if (v[i] > valid_range[1])
                    v[i] = valid_range[1];
            }
        } else
            enkf_quit("%s: %s: output for types other than NC_SHORT and NC_FLOAT are not handled yet", fname, varname);
    }

    ncw_put_vara_float(ncid, varid, start, count, v);
    ncw_close(ncid);
}

/**
 */
void read3dfield(char* fname, char* varname, int ni, int nj, int nk, float* v)
{
    int ncid;
    int varid;
    int ndims;
    size_t dimlen[4];
    size_t start[4], count[4];
    int i, n;
    int hasrecorddim;

    ncw_open(fname, NC_NOWRITE, &ncid);
    ncw_inq_varid(ncid, varname, &varid);
    ncw_inq_vardims(ncid, varid, 4, &ndims, dimlen);
    hasrecorddim = ncw_var_hasunlimdim(ncid, varid);
    if (hasrecorddim)
        if (dimlen[0] == 0)
            enkf_quit("%s: %s: empty record dimension", fname, varname);

    if (ndims == 4) {
        if (!hasrecorddim)
            enkf_quit("%s: %s: expect an unlimited dimension to be present for a 4-dimensional variable\n", fname, varname);
        start[0] = dimlen[0] - 1;
        start[1] = 0;
        start[2] = 0;
        start[3] = 0;
        count[0] = 1;
        count[1] = dimlen[1];
        count[2] = dimlen[2];
        count[3] = dimlen[3];
        if (dimlen[3] != ni || dimlen[2] != nj || dimlen[1] != nk)
            enkf_quit("\"%s\": horizontal dimensions of variable \"%s\" (ni = %d, nj = %d, nk = %d) do not match grid dimensions (ni = %d, nj = %d, nk = %d)", fname, varname, dimlen[3], dimlen[2], dimlen[1], ni, nj, nk);
    } else if (ndims == 3) {
        if (hasrecorddim)
            enkf_quit("%s: %s: can not read 3D field because the variable is only 2D", fname, varname);
        start[0] = 0;
        start[1] = 0;
        start[2] = 0;
        count[0] = dimlen[0];
        count[1] = dimlen[1];
        count[2] = dimlen[2];
        if (dimlen[2] != ni || dimlen[1] != nj || dimlen[0] != nk)
            enkf_quit("\"%s\": horizontal dimensions of variable \"%s\" (ni = %d, nj = %d, nk = %d) do not match grid dimensions (ni = %d, nj = %d, nk = %d)", fname, varname, dimlen[2], dimlen[1], dimlen[0], ni, nj, nk);
    } else
        enkf_quit("%s: can not read 3D field for \"%s\": # of dimensions = %d", fname, varname, ndims);

    ncw_get_vara_float(ncid, varid, start, count, v);

    n = 1;
    for (i = 0; i < ndims; ++i)
        n *= count[i];

    if (ncw_att_exists(ncid, varid, "valid_range")) {
        nc_type xtype;
        size_t len;

        ncw_inq_att(ncid, varid, "valid_range", &xtype, &len);
        if (len != 2)
            enkf_quit("%s: %s: \"valid_range\" attribute must have 2 elements", fname, varname);
        if (xtype == NC_SHORT) {
            int valid_range[2];
            int fill_value = SHRT_MIN;
            short* vv = malloc(n * sizeof(short));

            ncw_get_vara_short(ncid, varid, start, count, vv);

            if (ncw_att_exists(ncid, varid, "_FillValue"))
                ncw_get_att_int(ncid, varid, "_FillValue", &fill_value);
            else if (ncw_att_exists(ncid, varid, "missing_value"))
                ncw_get_att_int(ncid, varid, "missing_value", &fill_value);

            ncw_get_att_int(ncid, varid, "valid_range", valid_range);
            for (i = 0; i < n; ++i)
                if (vv[i] == fill_value || vv[i] < valid_range[0] || vv[i] > valid_range[1])
                    v[i] = NAN;
            free(vv);
        } else if (xtype == NC_FLOAT || xtype == NC_DOUBLE) {
            float valid_range[2];
            float fill_value = -FLT_MAX;

            /*
             * (we assume that double numbers outside of [-FLT_MAX, FLT_MAX]
             * interval are converted to either -FLT_MAX or FLT_MAX)
             */

            if (ncw_att_exists(ncid, varid, "_FillValue"))
                ncw_get_att_float(ncid, varid, "_FillValue", &fill_value);
            else if (ncw_att_exists(ncid, varid, "missing_value"))
                ncw_get_att_float(ncid, varid, "missing_value", &fill_value);

            ncw_get_att_float(ncid, varid, "valid_range", valid_range);
            for (i = 0; i < n; ++i)
                if (v[i] == fill_value || v[i] < valid_range[0] || v[i] > valid_range[1])
                    v[i] = NAN;
        } else
            enkf_quit("%s: %s: can not properly read NetCDF types other than NC_SHORT, NC_FLOAT, or NC_DOUBLE yet", fname, varname);
    }

    if (ncw_att_exists(ncid, varid, "scale_factor")) {
        float scale_factor;

        ncw_get_att_float(ncid, varid, "scale_factor", &scale_factor);
        for (i = 0; i < n; ++i)
            v[i] *= scale_factor;
    }

    if (ncw_att_exists(ncid, varid, "add_offset")) {
        float add_offset;

        ncw_get_att_float(ncid, varid, "add_offset", &add_offset);

        for (i = 0; i < n; ++i)
            v[i] += add_offset;
    }

    ncw_close(ncid);
}

/** Tries to determine whether the variable is 3D or 2D.
 */
int is3d(char fname[], char varname[])
{
    int ncid;
    int varid;
    int ndims;
    int hasrecorddim;

    ncw_open(fname, NC_NOWRITE, &ncid);
    ncw_inq_varid(ncid, varname, &varid);
    ncw_inq_varndims(ncid, varid, &ndims);
    if (ndims > 4)
        enkf_quit("%s: %s: do not know how to read a %d-dimensional variable", fname, varname, ndims);
    hasrecorddim = ncw_var_hasunlimdim(ncid, varid);
    ncw_close(ncid);

    ndims -= hasrecorddim;
    if (ndims < 2 || ndims > 3)
        enkf_quit("%s: %s: a 2D or 3D variable expected", fname, varname);

    return ndims == 3;
}

/**
 */
/**
 */
double date2day(char strdate[])
{
    char buf[MAXSTRLEN];
    char* token;
    char seps[] = " ";
    char seps2[] = "\n";
    double day, offset, multiple;

    strncpy(buf, strdate, MAXSTRLEN);
    if ((token = strtok(buf, seps)) == NULL)
        enkf_quit("could not understand date \"%s\"", strdate);
    if (!str2double(token, &day))
        enkf_quit("could not convert date \"%s\" to double", token);
    if ((token = strtok(NULL, seps2)) == NULL)
        enkf_quit("could not understand date \"%s\"", strdate);
    tunits_convert(token, &multiple, &offset);
    day = day * multiple + offset;

    return day;
}

/**
 */
int get_tshift(double date, double tstep, int centred)
{
    double offset = (centred) ? 0.5 : 0.0;

    return (int) floor(date / tstep + offset);
}

/** For debugging purposes - to be called from GDB.
 */
void print_matrix_double(int n, int m, double** A, char offset[])
{
    int i, j;

    for (i = 0; i < n; ++i) {
        printf("%s", offset);
        for (j = 0; j < m; ++j)
            printf("%10.5g ", fabs(A[j][i]) < EPS_DOUBLE ? 0.0 : A[j][i]);
        printf("\n");
    }
}

/** For debugging purposes - to be called from GDB.
 */
void print_matrix_float(int n, int m, float** A, char offset[])
{
    int i, j;

    for (i = 0; i < n; ++i) {
        printf("%s", offset);
        for (j = 0; j < m; ++j)
            printf("%10.5g ", fabsf(A[j][i]) < EPS_FLOAT ? 0.0 : (double) A[j][i]);
        printf("\n");
    }
}

/** For debugging purposes - to be called from GDB.
 */
void print_vector_double(int n, double* a, char offset[])
{
    int i;

    printf("%s", offset);
    for (i = 0; i < n; ++i)
        printf("%10.5g ", fabs(a[i]) < EPS_DOUBLE ? 0.0 : a[i]);
    printf("\n");
}

/** For debugging purposes - to be called from GDB.
 */
void print_vector_float(int n, float* a, char offset[])
{
    int i;

    printf("%s", offset);
    for (i = 0; i < n; ++i)
        printf("%10.5g ", fabsf(a[i]) < EPS_FLOAT ? 0.0 : (double) a[i]);
    printf("\n");
}

/**
 */
ENSOBSTYPE interpolate2d(double fi, double fj, int ni, int nj, float** v, int** mask, int periodic_i)
{
    int i1 = (int) floor(fi);
    double wi1 = ceil(fi) - fi;
    int i2 = (int) ceil(fi);
    double wi2 = fi - floor(fi);
    int j1 = (int) floor(fj);
    double wj1 = ceil(fj) - fj;
    int j2 = (int) ceil(fj);
    double wj2 = fj - floor(fj);
    double sum = 0.0;
    double w = 0.0;
    double ww;

    if (i1 == i2)
        wi1 = 1.0;
    if (j1 == j2)
        wj1 = 1.0;

    /*
     * Note that this section should be consistent with the similar section in 
     * model_xy2fij().
     */
    if (i1 == -1)
        i1 = (periodic_i) ? ni - 1 : i2;
    if (i2 == ni)
        i2 = (periodic_i) ? 0 : i1;
    if (j1 == -1)
        j1 = j2;
    if (j2 == nj)
        j2 = j1;

    assert(i1 >= 0 && i2 < ni && j1 >= 0 && j2 < nj);

    if (mask[j1][i1]) {
        ww = wj1 * wi1;
        sum += v[j1][i1] * ww;
        w += ww;
    }
    if (mask[j1][i2]) {
        ww = wj1 * wi2;
        sum += v[j1][i2] * ww;
        w += ww;
    }
    if (mask[j2][i1]) {
        ww = wj2 * wi1;
        sum += v[j2][i1] * ww;
        w += ww;
    }
    if (mask[j2][i2]) {
        ww = wj2 * wi2;
        sum += v[j2][i2] * ww;
        w += ww;
    }
    sum = sum / w;

    return (ENSOBSTYPE) sum;
}

ENSOBSTYPE average2d(size_t* ids, int n, float** v)
{
    double sum = 0.0;
    float* v0 = v[0];
    int i;

    for (i = 0; i < n; ++i)
        sum += v0[ids[i]];

    return (ENSOBSTYPE) (sum / (double) n);
}

/** A part of interpolate2d() that looks at mask in adjacent nodes only.
 */
int island(double fi, double fj, double fk, int ni, int nj, int ksurf, int** numlevels, int periodic_i)
{
    int i1 = (int) floor(fi);
    int i2 = (int) ceil(fi);
    int j1 = (int) floor(fj);
    int j2 = (int) ceil(fj);

    if (i1 == -1)
        i1 = (periodic_i) ? ni - 1 : i2;
    if (i2 == ni)
        i2 = (periodic_i) ? 0 : i1;
    if (j1 == -1)
        j1 = j2;
    if (j2 == nj)
        j2 = j1;

    if (ksurf > 0)
        fk = (double) ksurf - fk;
    if (fk == 0.0)
        return !(numlevels[j1][i1] || numlevels[j1][i2] || numlevels[j2][i1] || numlevels[j2][i2]);
    else {
        int k = ceil(fk);

        return !(numlevels[j1][i1] > k || numlevels[j1][i2] > k || numlevels[j2][i1] > k || numlevels[j2][i2] > k);
    }
}

/** Linearly interpolates a 3D field to fractional coordinates in index space.
 *  Assumes that integer k indices correspond to layer centres. E.g. for 
 *  fk = 1.2 the vertical weights are 0.8 of layer 1 and 0.2 of layer 2.
 */
ENSOBSTYPE interpolate3d(double fi, double fj, double fk, int ni, int nj, int nk, int ktop, float*** v, int** nlevels, int periodic_i)
{
    int i1 = (int) floor(fi);
    double wi1 = ceil(fi) - fi;
    int i2 = (int) ceil(fi);
    double wi2 = fi - floor(fi);
    int j1 = (int) floor(fj);
    double wj1 = ceil(fj) - fj;
    int j2 = (int) ceil(fj);
    double wj2 = fj - floor(fj);
    int k1, k2;
    int k1top, k2top;           /* layer number from the top */
    double wk1, wk2;
    double sum = 0.0;
    double w = 0.0;
    double ww;

    assert(ktop == 0 || ktop == nk - 1);

    /*
     * It is assumed that -0.5 <= fk <= nk - 0.5; so, when -0.5 <= fk <= 0 or
     * nk - 1 <= fk <= nk - 0.5 -- do not interpolate, take the layer value.
     */
    if (fk < 0.0)
        fk = 0.0;
    if (fk > (double) (nk - 1))
        fk = (double) (nk - 1);
    k1 = (int) floor(fk);
    k1top = abs(ktop - k1);
    wk1 = ceil(fk) - fk;
    k2 = (int) ceil(fk);
    k2top = abs(ktop - k2);
    wk2 = fk - floor(fk);

    if (i1 == i2)
        wi1 = 1.0;
    if (j1 == j2)
        wj1 = 1.0;
    if (k1 == k2)
        wk1 = 1.0;

    /*
     * Note that this section should be consistent with the similar section in 
     * model_xy2fij().
     */
    if (i1 == -1)
        i1 = (periodic_i) ? ni - 1 : i2;
    if (i2 == ni)
        i2 = (periodic_i) ? 0 : i1;
    if (j1 == -1)
        j1 = j2;
    if (j2 == nj)
        j2 = j1;

    assert(i1 >= 0 && i2 < ni && j1 >= 0 && j2 < nj && k1 >= 0 && k2 < nk);

    if (nlevels[j1][i1] > k1top) {
        ww = wj1 * wi1 * wk1;
        sum += v[k1][j1][i1] * ww;
        w += ww;
    }
    if (nlevels[j1][i2] > k1top) {
        ww = wj1 * wi2 * wk1;
        sum += v[k1][j1][i2] * ww;
        w += ww;
    }
    if (nlevels[j2][i1] > k1top) {
        ww = wj2 * wi1 * wk1;
        sum += v[k1][j2][i1] * ww;
        w += ww;
    }
    if (nlevels[j2][i2] > k1top) {
        ww = wj2 * wi2 * wk1;
        sum += v[k1][j2][i2] * ww;
        w += ww;
    }
    if (nlevels[j1][i1] > k2top) {
        ww = wj1 * wi1 * wk2;
        sum += v[k2][j1][i1] * ww;
        w += ww;
    }
    if (nlevels[j1][i2] > k2top) {
        ww = wj1 * wi2 * wk2;
        sum += v[k2][j1][i2] * ww;
        w += ww;
    }
    if (nlevels[j2][i1] > k2top) {
        ww = wj2 * wi1 * wk2;
        sum += v[k2][j2][i1] * ww;
        w += ww;
    }
    if (nlevels[j2][i2] > k2top) {
        ww = wj2 * wi2 * wk2;
        sum += v[k2][j2][i2] * ww;
        w += ww;
    }
    sum = sum / w;

    return (ENSOBSTYPE) sum;
}

/** Gaspary & Cohn's taper function.
 * @param x Support radius
 * @return Taper coefficient
 */
double taper_gc(double x)
{
    double x2, x3;

    assert(x >= 0 && x <= 1.0 + 1.0e-8);

    if (x >= 1.0)               /* handle possible round-up error */
        return 0.0;

    x *= 2.0;
    x2 = x * x;
    x3 = x2 * x;
    if (x < 1.0)
        return 1.0 + x2 * (-x3 / 4.0 + x2 / 2.0) + x3 * (5.0 / 8.0) - x2 * (5.0 / 3.0);
    return x2 * (x3 / 12.0 - x2 / 2.0) + x3 * (5.0 / 8.0) + x2 * (5.0 / 3.0) - x * 5.0 + 4.0 - (2.0 / 3.0) / x;
}

/** Converts from geographic to cartesian coordinates.
 * @param in Input: {lon, lat}
 * @param out Output: {x, y, z}
 */
void ll2xyz(double in[2], double out[3])
{
    double lon = in[0] * DEG2RAD;
    double lat = in[1] * DEG2RAD;
    double coslat = cos(lat);

    out[0] = REARTH * sin(lon) * coslat;
    out[1] = REARTH * cos(lon) * coslat;
    out[2] = REARTH * sin(lat);
}

/**
 */
void print_commandinfo(int argc, char* argv[])
{
    int i;
    char cwd[MAXSTRLEN];

    enkf_printf("    command = \"%s", argv[0]);
    for (i = 1; i < argc; ++i)
        enkf_printf(" %s", argv[i]);
    enkf_printf("\"\n");
    if (getcwd(cwd, MAXSTRLEN) != NULL)
        enkf_printf("    dir = \"%s\"\n", cwd);
}

/**
 */
void get_normalpair(double x[])
{
    double u = sqrt(-2.0 * log(1.0 - drand48()));
    double v = TWOPI * drand48();

    x[0] = u * cos(v);
    x[1] = u * sin(v);
}

/**
 */
int istrue(char str[])
{
    if (str == NULL)
        return 0;
    if (strcmp(str, "1") == 0)
        return 1;
    if (strcasecmp(str, "y") == 0)
        return 1;
    if (strcasecmp(str, "yes") == 0)
        return 1;
    if (strcasecmp(str, "t") == 0)
        return 1;
    if (strcasecmp(str, "true") == 0)
        return 1;
    return 0;
}

/**
 */
int inloninterval(double lon, double lon1, double lon2)
{
    while (lon2 < lon1)
        lon2 += 360.0;
    while (lon - 360.0 > lon1)
        lon = lon - 360.0;
    while (lon < lon1)
        lon += 360.0;
    return (lon <= lon2);
}

void shuffle(size_t n, size_t ids[])
{
    size_t i;

    for (i = 0; i < n; ++i) {
        size_t ii = (size_t) ((double) n * drand48());
        size_t tmp = ids[i];

        ids[i] = ids[ii];
        ids[ii] = tmp;
    }
}
