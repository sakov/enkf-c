
/******************************************************************************
 *
 * File:        observations.h        
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

#if !defined(_OBSERVATIONS_H)

#include "stringtable.h"
#include "hash.h"
#include "kdtree.h"
#include "enkfprm.h"
#include "grid.h"
#include "model.h"
#include "obstypes.h"

/*
 * (we keep the integer types signed for compatibility with netcdf 3)
 */
typedef struct {
    /*
     * observation index
     */
    int id;
    /*
     * for primary observations - the sequential number of observation; for
     * superobs - the original ID of the very first observation collated into
     * this observation
     */
    int id_orig;
    short int type;
    short int product;
    short int instrument;
    short int fid;
    int batch;
    double value;
    double estd;
    double footprint;
    double lon;
    double lat;
    double depth;
    double model_depth;
    double fi;
    double fj;
    double fk;
    double time;
    int status;                 /* 0 = OK */
    /*
     * auxiliary information:
     *  - # of obs for a sob
     *  - sob id for the original ob
     */
    int aux;
} observation;

typedef struct {
    stringtable* products;
    stringtable* instruments;
    stringtable* datafiles;

    int nobstypes;
    obstype* obstypes;
#if defined(ENKF_CALC)
    kdtree** loctrees;
    int** obsids;
#if defined(USE_SHMEM)
    MPI_Win* sm_comm_wins_kd;
    MPI_Win sm_comm_win_data;
#endif
#endif

    double da_time;             /* fractional days since 00:00:00
                                 * BASEDAY-BASEMONTH-BASEYEAR */
    char* datestr;

    int allobs;                 /* flag - whether to keep obs outside model
                                 * grid */
    int nallocated;

    int nobs;
    observation* data;
    int compacted;
    int has_nonpointobs;

    int ngood;
    int noutside_grid;
    int noutside_obsdomain;
    int noutside_obswindow;
    int nland;
    int nshallow;
    int nbadbatch;
    int nrange;
    int nthinned;
    int nexcluded;
    int nmodified;

    hashtable* badbatches;

    int ncformat;
    int nccompression;

#if defined(ENKF_PREP)
    void* model;
#endif
} observations;

observations* obs_create(void);
observations* obs_create_fromprm(enkfprm* prm);
observations* obs_create_fromdata(observations* parentobs, int nobs, observation data[]);
void obs_destroy(observations* obs);
void obs_checkalloc(observations* obs);
void obs_addtype(observations* obs, obstype* src);
void obs_compact(observations* obs);
void obs_inorder(observations* obs);
void obs_calcstats(observations* obs);
void obs_markbadbatches(observations* obs);
void obs_read(observations* obs, char fname[]);
void obs_write(observations* obs, char fname[]);
void obs_writeaux(observations* obs, char fname[]);
int obs_modifiederrors_alreadywritten(observations* obs, char fname[]);

#if defined(ENKF_PREP)
void obs_superob(observations* obs, __compar_d_fn_t cmp_obs, observations** sobs, int sobid, int do_thin);
#endif
void obs_find_bytype(observations* obs, int type, int* nobs, int** obsids);
void obs_find_bytypeandtime(observations* obs, int type, int time, int* nobs, int** obsids);
void obs_printob(observations* obs, int id);

#if defined(ENKF_CALC)
void obs_createkdtrees(observations* obs);
void obs_destroykdtrees(observations* obs);

#if defined(MINIMISE_ALLOC)
void obs_findlocal(observations* obs, double lon, double lat, char* dimainname, int* n, int** ids, double** lcoeffs, int* ploc_allocated);
#else
void obs_findlocal(observations* obs, double lon, double lat, char* domainname, int* n, int** ids, double** lcoeffs);
#endif
#endif

#define _OBSERVATIONS_H
#endif
