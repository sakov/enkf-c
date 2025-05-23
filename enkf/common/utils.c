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
 *              6/9/2019 PS - moved NetCDF procedures to ncutils.c
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
#include <glob.h>
#include <ftw.h>
#if defined(MPI)
#include <mpi.h>
#endif
#include "ncw.h"
#include "ncutils.h"
#include "triangulation.h"
#if defined(MPI) && defined(USE_MPIQUEUE)
#include "mpiqueue.h"
#endif
#include "definitions.h"
#include "version.h"
#include "utils.h"

#define FILE_FIND_INC 10
#define EPS_DOUBLE 4.0e-15
#define EPS_FLOAT  1.0e-7
#define EPS_FIJ 1.0e-4
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
static char* get_command(int argc, char* argv[])
{
    char* cmd = NULL;
    int len = 0;
    int i;

    len = strlen(argv[0]);
    for (i = 1; i < argc; ++i)
        len += strlen(argv[i]) + 1;
    len++;

    cmd = malloc(len);
    strcpy(cmd, argv[0]);
    len = strlen(argv[0]);
    for (i = 1; i < argc; ++i) {
        cmd[len] = ' ';
        len++;
        strcpy(&cmd[len], argv[i]);
        len += strlen(argv[i]);
    }
    cmd[len] = 0;

    return cmd;
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

    if (*argc == 2 && strcmp((*argv)[1], "--version") == 0)
        return;

    if (*argc > 1) {
        enkf_printf("  MPI: initialised %d process(es)\n", nprocesses);
        MPI_Barrier(MPI_COMM_WORLD);
#if defined(USE_SHMEM)
        /*
         * initialise communicators for handling shared memory stuff
         */
        {
            int ierror;
            int* recvcounts = NULL;
            int* displs = NULL;
            int i;

            ierror = MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &sm_comm);
            assert(ierror == MPI_SUCCESS);
            ierror = MPI_Comm_rank(sm_comm, &sm_comm_rank);
            assert(ierror == MPI_SUCCESS);
            ierror = MPI_Comm_size(sm_comm, &sm_comm_size);
            assert(ierror == MPI_SUCCESS);
            sm_comm_ranks = malloc(nprocesses * sizeof(int));
            /*
             * build map of local (i.e. within the node the core belongs to)
             * ranks
             */
            sm_comm_ranks[rank] = sm_comm_rank;
            recvcounts = malloc(nprocesses * sizeof(int));
            displs = malloc(nprocesses * sizeof(int));
            for (i = 0; i < nprocesses; ++i) {
                recvcounts[i] = 1;
                displs[i] = i;
            }
            ierror = MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, sm_comm_ranks, recvcounts, displs, MPI_INT, MPI_COMM_WORLD);
            assert(ierror == MPI_SUCCESS);
            /*
             * create communicators based on local ranks
             */
            ierror = MPI_Comm_split(MPI_COMM_WORLD, sm_comm_rank, rank, &node_comm);
            assert(ierror == MPI_SUCCESS);
            ierror = MPI_Comm_rank(node_comm, &node_comm_rank);
            assert(ierror == MPI_SUCCESS);
            ierror = MPI_Comm_size(node_comm, &node_comm_size);
            assert(ierror == MPI_SUCCESS);
            /*
             * Free communicators for local ranks other than 0. The communicator
             * for local rank 0 will be used for gathering S and St.
             */
            if (sm_comm_rank != 0) {
                MPI_Comm_free(&node_comm);
                node_comm_rank = -1;
            }
            node_comm_ranks = malloc(nprocesses * sizeof(int));
            node_comm_ranks[rank] = node_comm_rank;
            ierror = MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, node_comm_ranks, recvcounts, displs, MPI_INT, MPI_COMM_WORLD);
            assert(ierror == MPI_SUCCESS);

            free(recvcounts);
            free(displs);

            enkf_printf("  Using MPI-3 shared memory:\n");
            enkf_printf("    sm_comm size = %d\n", sm_comm_size);
            enkf_printf("    node_comm size = %d\n", node_comm_size);
        }
        if (enkf_verbose > 1) {
            MPI_Barrier(MPI_COMM_WORLD);
            printf("  MPI: rank = %3d, sm_comm_rank = %2d, node_comm_rank = %2d, PID = %d\n", rank, sm_comm_rank, node_comm_rank, getpid());
            fflush(NULL);
        }
#else
        if (enkf_verbose > 1) {
            MPI_Barrier(MPI_COMM_WORLD);
            printf("  MPI: rank = %d, PID = %d\n", rank, getpid());
            fflush(NULL);
        }
#endif                          /* USE_SHMEM */
    }

    MPI_Type_size(MPI_INT, &size);
    assert(size == sizeof(int));
    MPI_Type_size(MPI_FLOAT, &size);
    assert(size == sizeof(float));
#else
    if (*argc == 2 && strcmp((*argv)[1], "--version") == 0)
        return;
    if (enkf_verbose > 1)
        printf("  rank = %d, PID = %d\n", rank, getpid());
#endif                          /* MPI */

#if (defined(ENKF_CALC) && defined(TW_VIAFILE)) || defined(ENKF_UPDATE)
    if (rank == 0 && file_exists(DIRNAME_TMP))
        enkf_quit("directory \"%s\" already exists", DIRNAME_TMP);
#endif

    ncw_set_quitfn(enkf_quit);
    ncu_set_quitfn(enkf_quit);
#if defined(ENKF_CALC)
    kd_set_quitfn(enkf_quit);
#endif
#if defined(ENKF_PREP) || defined(ENKF_CALC)
    triangulation_set_quitfn(enkf_quit);
#endif
#if defined(MPI) && defined(USE_MPIQUEUE)
    mpiqueue_setquitfn(enkf_quit);
#endif

    /*
     * initialise the random number generator to a random state for each cpu
     */
    randomise_rand48();

    enkf_cmd = get_command(*argc, *argv);
    enkf_cwd = getcwd(NULL, 0);
}

/**
 */
void enkf_finish(void)
{
#if defined(MPI)
    MPI_Barrier(MPI_COMM_WORLD);
#endif
#if defined (USE_SHMEM)
    if (sm_comm != MPI_COMM_NULL)
        MPI_Comm_free(&sm_comm);
    if (sm_comm_ranks != NULL)
        free(sm_comm_ranks);
    if (node_comm != MPI_COMM_NULL)
        MPI_Comm_free(&node_comm);
    if (node_comm_ranks != NULL)
        free(node_comm_ranks);
#endif
    enkf_printtime("  ");
    if (enkf_cmd != NULL)
        free(enkf_cmd);
    if (enkf_cwd != NULL)
        free(enkf_cwd);
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

    if (enkf_verbose < 0 || (enkf_verbose && rank == 0)) {
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
#if defined(ENKF_PREP)
    enkf_printf("%senkf_prep compile flags:\n", offset);
#elif defined(ENKF_CALC)
    enkf_printf("%senkf_calc compile flags:\n", offset);
#elif defined(ENKF_UPDATE)
    enkf_printf("%senkf_update compile flags:\n", offset);
#elif defined(ENKF_DIAG)
    enkf_printf("%senkf_diag compile flags:\n", offset);
#else
    enkf_quit("programming error");
#endif
#if defined(ENKF_CALC)
#if defined(SHUFFLE_ROWS)
    enkf_printf("%s  SHUFFLE_ROWS     = [+]\n", offset);
#else
    enkf_printf("%s  SHUFFLE_ROWS     = [-]\n", offset);
#endif
#if defined(USE_SHMEM)
    enkf_printf("%s  USE_SHMEM        = [+]\n", offset);
#else
    enkf_printf("%s  USE_SHMEM        = [-]\n", offset);
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
#if defined(TW_VIAFILE)
    enkf_printf("%s  TW_VIAFILE       = [+]\n", offset);
#else
    enkf_printf("%s  TW_VIAFILE       = [-]\n", offset);
#endif
#if defined(USE_MPIQUEUE)
    enkf_printf("%s  USE_MPIQUEUE     = [+]\n", offset);
#else
    enkf_printf("%s  USE_MPIQUEUE     = [-]\n", offset);
#endif
#endif                          /* ENKF_CALC */
#if defined(ENKF_PREP) || defined(ENKF_CALC)
#if defined(INTERNAL_QSORT_R)
    enkf_printf("%s  INTERNAL_QSORT_R = [+]\n", offset);
#else
    enkf_printf("%s  INTERNAL_QSORT_R = [-]\n", offset);
#endif
#endif
#if defined(ENKF_UPDATE) || defined(ENKF_DIAG)
#if defined(NCW_SKIPSINGLE)
    enkf_printf("%s  NCW_SKIPSINGLE   = [+]\n", offset);
#else
    enkf_printf("%s  NCW_SKIPSINGLE   = [-]\n", offset);
#endif
#endif
#if defined(ENKF_UPDATE)
#if defined(USE_MPIQUEUE)
    enkf_printf("%s  USE_MPIQUEUE     = [+]\n", offset);
#else
    enkf_printf("%s  USE_MPIQUEUE     = [-]\n", offset);
#endif
#endif                          /* ENKF_UPDATE */
#if defined(ENKF_CALC) || defined(ENKF_UPDATE)
#if defined(DEFLATE_ALL)
    enkf_printf("%s  DEFLATE_ALL      = [+]\n", offset);
#else
    enkf_printf("%s  DEFLATE_ALL      = [-]\n", offset);
#endif
#endif
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
    enkf_printf("  EnKF-C version %s\n", ENKF_VERSION);
    enkf_printcompileflags("  ");
}

/** Find files matching the template using glob.
 */
void find_files(char* template, int* nfiles, char*** fnames)
{
    glob_t gl;
    int status;

    status = glob(template, GLOB_BRACE | GLOB_PERIOD | GLOB_TILDE_CHECK, NULL, &gl);
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

    if (strncasecmp(tunits, "microsec", 8) == 0)
        *tunits_multiple = 1.0 / 86400000000.0;
    else if (strncasecmp(tunits, "millisec", 8) == 0)
        *tunits_multiple = 1.0 / 86400000.0;
    else if (strncasecmp(tunits, "sec", 3) == 0)
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

    *tunits_offset = (double) daydiff(year, month, day, BASEYEAR, BASEMONTH, BASEDAY) + (double) h / 24.0 + (double) m / 1440.0 + (double) s / 86400.0;
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
void file_delete(char fname[])
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
void file_rename(char oldname[], char newname[])
{
    int status = -1;

    status = rename(oldname, newname);
    if (status != 0) {
        int errno_saved = errno;

        enkf_quit("could not rename file \"%s\" to \"%s\": %s", oldname, newname, strerror(errno_saved));
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

        enkf_quit("dir_rmifexists(): \"%s\": %s", dirname, strerror(errno_saved));
    }
}

/**
 */
static int rmentry(const char* entry, const struct stat* sb, int typeflag, struct FTW* ftwbuf)
{
    return remove(entry);
}

/**
 */
void dir_rmallifexists(char dirname[])
{
    int status;

    if (!dir_exists(dirname))
        return;
    status = nftw(dirname, rmentry, 64, FTW_DEPTH | FTW_PHYS);
    if (status != 0) {
        int errno_saved = errno;

        enkf_quit("dir_rmallifexists(): \"%s\": %s", dirname, strerror(errno_saved));
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
    if (token[0] == 'y' || token[0] == 'Y' || token[0] == 't' || token[0] == 'T')
        *value = 1;
    else if (token[0] == 'n' || token[0] == 'N' || token[0] == 'f' || token[0] == 'F')
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

/** Casts an ni x nj matrix onto a pre-allocated storage and fills it with
 ** zeros.
 * @param p Address of allocated bloc
 * @param nj Dimension 2
 * @param ni Dimension 1
 * @param unitsize Size of one matrix element in bytes
 * @return Matrix (= (unit**) p)
 */
void* cast2d(void* p, size_t nj, size_t ni, size_t unitsize)
{
    size_t size;
    void** pp;
    int i;

    if (ni <= 0 || nj <= 0)
        enkf_quit("alloc2d(): invalid size (nj = %d, ni = %d)", nj, ni);

    size = nj * sizeof(void*) + nj * ni * unitsize;
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
double date2day(char* fname, char* strdate)
{
    char buf[MAXSTRLEN];
    char* token;
    char seps[] = " ";
    char seps2[] = "\n";
    double day, offset, multiple;

    strncpy(buf, strdate, MAXSTRLEN - 1);
    if ((token = strtok(buf, seps)) == NULL)
        enkf_quit("%s: date2day(): could not understand date \"%s\"", fname, strdate);
    if (!str2double(token, &day))
        enkf_quit("%s: date2day(): \"%s\": could not convert \"%s\" to double", fname, strdate, token);
    if ((token = strtok(NULL, seps2)) == NULL) {
#if 0
        enkf_quit("%s: %s: date2day(): could not understand date \"%s\"", fname, strdat);
#else
        enkf_geophysical = 0;
#endif
    } else {
        tunits_convert(token, &multiple, &offset);
        day = day * multiple + offset;
    }

    return day;
}

/**
 */
int get_tshift(double reltime, double tstep, int centred)
{
    double offset = (centred) ? 0.5 : 0.0;

    return (int) floor(reltime / tstep + offset);
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
float interpolate2d_structured(double* fij, int ni, int nj, float** v, int** mask, int periodic_i)
{
    int i1, i2, j1, j2;
    double wi1, wi2, wj1, wj2;
    double sum, w, ww;

    double fi = fij[0];
    double fj = fij[1];

    /*
     * (very rarely) superobs need to be fixed because of the round-off errors
     */
    if (fi >= (double) ni)
        fi = (double) ni - EPS_FLOAT;
    if (fj >= (double) nj)
        fj = (double) nj - EPS_FLOAT;

    i1 = (int) floor(fi);
    wi1 = ceil(fi) - fi;
    i2 = (int) ceil(fi);
    wi2 = fi - floor(fi);
    j1 = (int) floor(fj);
    wj1 = ceil(fj) - fj;
    j2 = (int) ceil(fj);
    wj2 = fj - floor(fj);

    if (i1 == i2)
        wi1 = 1.0;
    if (j1 == j2)
        wj1 = 1.0;

    /*
     * Note that this section should be consistent with the similar section in 
     * grid_xy2fij().
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

    sum = 0.0;
    w = 0.0;
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
    sum /= w;

    return (float) sum;
}

/**
 */
float interpolate2d_unstructured(double* fi, float* v, int* mask)
{
    double sum = 0.0, w = 0.0;
    int id;

    for (id = 0; id < 3; ++id) {
        int ii = (int) fi[id];

        if (mask[ii]) {
            double ww = fi[id] - floor(fi[id]);

            sum += ww * v[(int) fi[id]];
            w += ww;
        }
    }
    sum /= w;

    return (float) sum;
}

/**
 */
float average(int n, size_t* ids, float* v)
{
    double sum = 0.0;
    int i;

    for (i = 0; i < n; ++i)
        sum += v[ids[i]];

    return (float) (sum / (double) n);
}

/** Linearly interpolates a 3D field to fractional coordinates in index space.
 *  Assumes that integer k indices correspond to layer centres. E.g. for 
 *  fk = 1.2 the vertical weights are 0.8 of layer 1 and 0.2 of layer 2.
 */
float interpolate3d_structured(double* fij, double fk, int ni, int nj, int nk, int ktop, float*** v, int** nlevels, int periodic_i)
{
    int i1, i2, j1, j2, k1, k2;
    double wi1, wi2, wj1, wj2, wk1, wk2;
    int k1top, k2top;           /* layer number from the top */
    double sum, w, ww;

    double fi = fij[0];
    double fj = fij[1];

    /*
     * (very rarely) superobs need to be fixed because of the round-off errors
     */
    if (fi > (double) ni)
        fi = (double) ni - EPS_FLOAT;
    if (fj >= (double) nj)
        fj = (double) nj - EPS_FLOAT;

    i1 = (int) floor(fi);
    wi1 = ceil(fi) - fi;
    i2 = (int) ceil(fi);
    wi2 = fi - floor(fi);
    j1 = (int) floor(fj);
    wj1 = ceil(fj) - fj;
    j2 = (int) ceil(fj);
    wj2 = fj - floor(fj);

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
     * grid_xy2fij().
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

    sum = 0.0;
    w = 0.0;
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

    return (float) sum;
}

/**
 */
float interpolate3d_unstructured(double* fi, double fk, int nk, int ktop, float** v, int* nlevels)
{
    double sum = 0.0;
    double w = 0.0;
    int k1, k2;
    int k1top, k2top;           /* layer number from the top */
    double wk1, wk2;
    int id;

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

    for (id = 0; id < 3; ++id) {
        int ii = (int) fi[id];

        if (nlevels[ii] > k1top) {
            double wi = fi[id] - floor(fi[id]);
            double ww = wi * wk1;

            sum += ww * v[k1][(int) fi[id]];
            w += ww;
        }
        if (nlevels[ii] > k2top) {
            double wi = fi[id] - floor(fi[id]);
            double ww = wi * wk2;

            sum += ww * v[k2][(int) fi[id]];
            w += ww;
        }
    }

    sum /= w;

    return (float) sum;
}

/** Gaspary & Cohn's taper function.
 * @param x Support radius
 * @return Taper coefficient
 */
double taper_gc(double x)
{
    double x2, x3, f;

    assert(x >= 0 && x <= 1.0 + 1.0e-8);

    if (x >= 1.0)               /* handle possible round-up error */
        return 0.0;

    x *= 2.0;
    x2 = x * x;
    x3 = x2 * x;
    if (x < 1.0)
        f = 1.0 + x2 * (-x3 / 4.0 + x2 / 2.0) + x3 * (5.0 / 8.0) - x2 * (5.0 / 3.0);
    else
        f = x2 * (x3 / 12.0 - x2 / 2.0) + x3 * (5.0 / 8.0) + x2 * (5.0 / 3.0) - x * 5.0 + 4.0 - (2.0 / 3.0) / x;

    return (f >= 0.0) ? f : 0.0;
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

    out[0] = REARTH * cos(lon) * coslat;
    out[1] = REARTH * sin(lon) * coslat;
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

/**
 */
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

#if defined(ENKF_PREP) || defined(ENKF_CALC)
/**
 */
void kd_printinfo(kdtree* tree, char* offset)
{
    size_t nnode;
    size_t nalloc;

    if (rank != 0)
        return;
    if (tree == NULL)
        return;

    nnode = kd_getsize(tree);
    nalloc = kd_getnalloc(tree);

    enkf_printf("%skdtree \"%s\":\n", offset == NULL ? "" : offset, kd_getname(tree));
    enkf_printf("%s  %zu nodes", offset == NULL ? "" : offset, nnode);
    if (nalloc != 0) {
        if (nnode != nalloc)
            enkf_printf("(%zu allocated)\n", nalloc);
        else
            enkf_printf("\n");
        enkf_printf("%s  %zu bytes\n", offset == NULL ? "" : offset, kd_getstoragesize(tree, nalloc));
    } else {
        enkf_printf("\n");
        enkf_printf("%s  %zu bytes (externally allocated)\n", offset == NULL ? "" : offset, kd_getstoragesize(tree, 0));
    }
    enkf_flush();
}
#endif

#include <memory.h>

/**
 */
static int get_memory_usage_kb(size_t* vmrss_kb, size_t* vmsize_kb)
{
    /*
     * Get the the current process' status file from the proc filesystem 
     */
    FILE* procfile = fopen("/proc/self/status", "r");

    size_t to_read = 8192;
    char buffer[to_read];

    /*
     * (dummy if to avoid warning from GCC) 
     */
    if (fread(buffer, sizeof(char), to_read, procfile));
    fclose(procfile);

    int found_vmrss = 0;
    int found_vmsize = 0;
    char* search_result;

    /*
     * Look through proc status contents line by line 
     */
    char delims[] = "\n";
    char* line = strtok(buffer, delims);

    while (line != NULL && (found_vmrss == 0 || found_vmsize == 0)) {
        search_result = strstr(line, "VmRSS:");
        if (search_result != NULL) {
            sscanf(line, "%*s %zu", vmrss_kb);
            found_vmrss = 1;
        }

        search_result = strstr(line, "VmSize:");
        if (search_result != NULL) {
            sscanf(line, "%*s %zu", vmsize_kb);
            found_vmsize = 1;
        }

        line = strtok(NULL, delims);
    }

    return (found_vmrss == 1 && found_vmsize == 1) ? 0 : 1;
}

/**
 */
static int get_cluster_memory_usage_kb(size_t* vmrss_per_process, size_t* vmsize_per_process, int root, int np)
{
    size_t vmrss_kb;
    size_t vmsize_kb;
    int ret_code = get_memory_usage_kb(&vmrss_kb, &vmsize_kb);

    if (ret_code != 0)
        enkf_quit("could not gather memory usage: %s\n", strerror(ret_code));

#if defined(MPI)
    MPI_Gather(&vmrss_kb, 1, MPI_UNSIGNED_LONG, vmrss_per_process, 1, MPI_UNSIGNED_LONG, root, MPI_COMM_WORLD);

    MPI_Gather(&vmsize_kb, 1, MPI_UNSIGNED_LONG, vmsize_per_process, 1, MPI_UNSIGNED_LONG, root, MPI_COMM_WORLD);
#else
    vmrss_per_process[0] = vmrss_kb;
    vmsize_per_process[0] = vmsize_kb;
#endif

    return 0;
}

/**
 */
void print_memory_usage(void)
{
    size_t vmrss_per_process[nprocesses];
    size_t vmsize_per_process[nprocesses];
    size_t vmrss, vmsize;

    get_cluster_memory_usage_kb(vmrss_per_process, vmsize_per_process, 0, nprocesses);

    if (rank == 0) {
        int i;

        enkf_printf("  memory usage:\n");
        if (nprocesses > 1) {
            for (i = 0, vmrss = 0, vmsize = 0; i < nprocesses; i++) {
                enkf_printf("    process %03d: VmRSS = %zu kB, VmSize = %zu kB\n", i, vmrss_per_process[i], vmsize_per_process[i]);
                vmrss += vmrss_per_process[i];
                vmsize += vmsize_per_process[i];
            }
            enkf_printf("    total: VmRSS = %zu kB, VmSize = %zu kB\n", vmrss, vmsize);
        } else
            enkf_printf("    VmRSS = %zu kB, VmSize = %zu kB\n", vmrss_per_process[0], vmsize_per_process[0]);
    }
}

/**
 */
void print_memory_avail(void)
{
    if (rank == 0) {
        size_t npage = sysconf(_SC_PHYS_PAGES);
        size_t page_size = sysconf(_SC_PAGESIZE);

        enkf_printf("  available RAM per node:\n");
        enkf_printf("    no. pages = %zu\n", npage);
        enkf_printf("    page size = %zu\n", page_size);
        enkf_printf("    total = %zu\n", npage * page_size);
#if defined(USE_SHMEM)
        enkf_printf("    (%zu bytes per core)\n", npage * page_size / sm_comm_size);
#endif
    }
}

/**
 */
void enkf_writeinfo(char* fname)
{
    if (rank == 0) {
        int ncid;

        ncw_open(fname, NC_WRITE, &ncid);
        if (!ncw_att_exists(ncid, NC_GLOBAL, "EnKF-C version")) {
            ncw_put_att_text(ncid, NC_GLOBAL, "EnKF-C version", ENKF_VERSION);
            ncw_put_att_text(ncid, NC_GLOBAL, "command", enkf_cmd);
            ncw_put_att_text(ncid, NC_GLOBAL, "wdir", enkf_cwd);
        }
        ncw_close(ncid);
    }
}
