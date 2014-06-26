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

#include "definitions.h"
#include "stringtable.h"
#include "kdtree.h"
#include "enkfprm.h"
#include "grid.h"

typedef struct {
    int type;
    int product;
    int instrument;
    int id;
    /*
     * for primary observations - the original ID corresponding to the number
     * of the primary observation during the very first read of data files; for
     * superobs - the original ID of the very first observation of those
     * combined
     */
    int id_orig;
    double value;
    double std;
    double lon;
    double lat;
    double depth;
    double fi;
    double fj;
    double fk;
    double date;
    int status;                 /* 0 = OK */
    /*
     * auxiliary information:
     *  - # of obs for a sob
     *  - sob id for the original ob
     */
    int aux;
} measurement;

typedef struct {
    int id;
    char* name;
    int issurface;
    char* varname;
    char* hfunction;
    double allowed_min;
    double allowed_max;
    int isasync;
    double async_tstep;
    int nobs;
    int ngood;
    int noutside;
    int nland;
    int nshallow;
    int nrange;
    int nmodified;
    double date_min;
    double date_max;
    double rfactor;
} obstype;

typedef struct {
    stringtable* types;
    stringtable* products;
    stringtable* instruments;

    int nobstypes;
    obstype* obstypes;

    double da_date;             /* days since 00:00:00
                                 * BASEDAY-BASEMONTH-BASEYEAR */
    char* datestr;

    int allobs;                 /* flag - whether to keep obs outside model
                                 * grid */
    int nobs;
    measurement* data;
    int compacted;

    int stride;

    int hasstats;
    int ngood;
    int noutside;
    int nland;
    int nshallow;
    int nrange;
    int nmodified;

    kdtree* tree;
} observations;

typedef struct {
    char* typename;
    int issurface;
    double minval;
    double maxval;
} obstypedesc;

extern obstypedesc otdescs[];
extern int notdescs;

observations* obs_create(void);
observations* obs_create_fromprm(enkfprm* prm);
observations* obs_create_fromdata(observations* parentobs, int nobs, measurement data[]);
void obs_destroy(observations* obs);
void obs_addtype(observations* obs, char name[], int issurface, char varname[], char hfunction[], double rfactor, int isasync, double async_tstep);
void obs_checklon(observations* obs);
void obs_compact(observations* obs);
void obs_calcstats(observations* obs);
void obs_write(observations* obs, char fname[]);
void obs_superob(observations* obs, __compar_d_fn_t cmp_obs, observations** sobs, int sobid);
void obs_writeaux(observations* obs, char fname[]);
void obs_find_bytype(observations* obs, int type, int* nobs, int** obsids, int fstatsonly);
void obs_find_bytypeandtime(observations* obs, int type, int time, int* nobs, int** obsids, int fstatsonly);
void obs_printob(observations* obs, int id);
void obs_createkdtree(observations* obs, grid* g);
void obs_findlocal(observations* obs, grid* g, double lon, double lat, double locrad, int* n, int** ids, double** lcoeffs);

#define _OBSERVATIONS_H
#endif
