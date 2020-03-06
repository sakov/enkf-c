/******************************************************************************
 *
 * File:        dasystem.h        
 *
 * Created:     12/2012
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description:
 *
 * Revisions:   06032020 PS: moved MPI communicators etc. to global.c
 *
 *****************************************************************************/

#if !defined(_DASYSTEM_H)

#include "ncw.h"
#include "observations.h"
#include "model.h"

#define S_MODE_NONE 0
#define S_MODE_HE_f 1
#define S_MODE_HA_f 2
#define S_MODE_S_f  3
#define S_MODE_HE_a 4
#define S_MODE_HA_a 5
#define S_MODE_S_a  6

typedef struct {
    int id;
    int varid;
    char varname[NC_MAX_NAME];
    int level;
} field;

typedef struct {
    char* prmfname;
    int mode;
    int scheme;
    double alpha;               /* moderating multiple */
    char* ensdir;
    char* bgdir;

    int nmem;

    model* m;

#if defined(ENKF_CALC)
    observations* obs;
    /*
     * Currently there are two observation sorting modes used by the code: by ID
     * and by IJ. Correspondingly, the `sort_mode' flag is set to either
     * OBS_SORTMODE_ID or OBS_SORTMODE_IJ. Note that at any time the
     * correspondence between array of observations obs->data and arrays of
     * ensemble observations/innovations below is assumed and should be
     * maintained: obs->data[i] corresponds to S[:][i], s_f[i], std_f[i],
     * s_a[i] and std_a[i].
     */
    int sort_mode;

    /*
     * The ensemble observations `S' and innovations `s' can be either in
     * standardised or not standardised mode. Further, `S' is too big to be
     * copied, instead it is modified when transferred between states. The
     * states of `s' and 'S' are described by `s_mode':
     *
     *   s_mode      s_f, s_a and S      S contains    S relates to
     *               are standardised    anomalies     forecast
     * S_MODE_HE_f         no                no           yes
     * S_MODE_HA_f         no                yes          yes
     * S_MODE_S_f          yes               yes          yes
     * S_MODE_HE_a         no                no           no
     * S_MODE_HA_a         no                yes          no
     * S_MODE_S_a          yes               yes          no
     */
    ENSOBSTYPE** S;             /* HE or HA or S [mem][obs] */
    double* s_f;                /* innovation */
    double* std_f;              /* ensemble spread */
    double* s_a;                /* innovation */
    double* std_a;              /* ensemble spread */
    int s_mode;

    double kfactor;
#endif
#if defined(HE_VIASHMEM)
    MPI_Win sm_comm_win_S;
    ENSOBSTYPE** St;            /* (S transposed) */
#endif
    int fieldbufsize;

    int nregions;
    region* regions;

    int nplog;
    pointlog* plogs;

    int nbadbatchspecs;
    badbatchspec* badbatchspecs;

    int updatespec;             /* binary flags */

    int ncformat;
    int nccompression;
} dasystem;

dasystem* das_create(enkfprm* prm);
void das_destroy(dasystem* das);

void das_getHE(dasystem* das);
void das_addanalysis(dasystem* das, char fname[]);
void das_addforecast(dasystem* das, char fname[]);
void das_addmodifiederrors(dasystem* das, char fname[]);
void das_writevcorrs(dasystem* das);
void das_calcinnandspread(dasystem* das);
void das_calctransforms(dasystem* das);
void das_calcpointlogtransforms(dasystem* das);
void das_getfields(dasystem* das, int gridid, int* nfield, field** fields);
void getfieldfname(char* dir, char* prefix, char* varname, int level, char* fname);
void das_moderateobs(dasystem* das);
void das_calcbatchstats(dasystem* das, int doprint);
void das_printobsstats(dasystem* das, int use_rmsd);
void das_printfobsstats(dasystem* das, int use_rmsd);
void das_readobs(dasystem* das, char fname[]);
void das_standardise(dasystem* das);
void das_destandardise(dasystem* das);
void das_update(dasystem* das);
void das_updateHE(dasystem* das);

void das_getfname_X5(dasystem* das, void* grid, char fname[]);
void das_getfname_w(dasystem* das, void* grid, char fname[]);
void das_getfname_stats(dasystem* das, void* grid, char fname[]);
void das_getfname_plog(dasystem* das, pointlog* plog, char fname[]);

void das_calcmld(dasystem* das, obstype* ot, float*** src, float** dst);

void plog_create(dasystem* das, int plogid, int ploc, int* lobs, double* lcoeffs);
void plog_writetransform(dasystem* das, int plogid, int gid, int ploc, double* s, double* S, double* transform);
void plog_definestatevars(dasystem* das);
void plog_writestatevars(dasystem* das, int nfields, void** fieldbuffer, field* fields, int isanalysis);
void plog_assemblestatevars(dasystem* das);

void das_getmemberfname(dasystem* das, char ensdir[], char varname[], int mem, char fname[]);
int das_getmemberfname_async(dasystem* das, char ensdir[], obstype* ot, int mem, int t, char fname[]);
void das_getbgfname(dasystem* das, char ensdir[], char varname[], char fname[]);
int das_getbgfname_async(dasystem* das, char bgdir[], obstype* ot, int t, char fname[]);

#define _DASYSTEM_H
#endif
