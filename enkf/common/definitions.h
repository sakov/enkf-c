/******************************************************************************
 *
 * File:        definitions.h        
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

#if !defined(_DEFINITIONS_H)

#ifdef DMALLOC
#include "dmalloc.h"
#endif

#define MAXSTRLEN 2048
/*
 * Memory allocation increment for observations. One might increase it for
 * large systems and possibly reduce for small systems.
 */
#define NOBS_INC 500000

#define BASEYEAR 1990
#define BASEMONTH 1
#define BASEDAY 1

#define FNAME_OBS "observations-orig.nc"
#define FNAME_SOBS "observations.nc"
#define FNAMEPREFIX_X5 "X5"
#define FNAMEPREFIX_W "w"
#define FNAMEPREFIX_DIAG "enkf_diag"
#define FNAMEPREFIX_PLOG "pointlog"
#define FNAME_SPREAD "spread.nc"
#define FNAME_HE "HE.nc"
#define FNAME_BADBATCHES "badbatches.out"
#define FNAME_INFLATION "inflation.nc"

#define STATUS_OK 0             /* do not change */
#define STATUS_OUTSIDEGRID 1
#define STATUS_LAND 2
#define STATUS_SHALLOW 3
#define STATUS_RANGE 4
#define STATUS_BADBATCH 5
#define STATUS_OUTSIDEOBSDOMAIN 6
/*
 * The fractional grid indices can (very rarely) jump to a neighbour cell after
 * being saved to observations.nc and read back. If they end up on land - mark
 * them with this status.
 */
#define STATUS_ROUNDUP 7

#define LONTYPE_NONE 0
#define LONTYPE_180 1
#define LONTYPE_360 2

#define MODE_NONE 0
#define MODE_ENKF 1
#define MODE_ENOI 2

#define ALPHA_DEFAULT 1.0

#define SCHEME_NONE 0
#define SCHEME_DENKF 1
#define SCHEME_ETKF 2
#define SCHEME_DEFAULT SCHEME_DENKF

#define EXITACTION_BACKTRACE 0
#define EXITACTION_SEGFAULT 1
#define EXITACTION_DEFAULT EXITACTION_BACKTRACE

#define OBSTYPE_VALUE 0
#define OBSTYPE_INNOVATION 1

#define INFLATION_PLAIN 0
#define INFLATION_SPREADLIMITED 1
#define INFLATION_DEFAULT INFLATION_SPREADLIMITED
/*
 * Maximal allowed default ratio of (inflation - 1) / (spread_f / spread_a - 1).
 * Use of this limiter can be overriden by using qualifier PLAIN for inflation.
 */
#define INFRATIO_DEFAULT 1.0

#define OBS_SORTMODE_ID 0
#define OBS_SORTMODE_IJ 1

#define ENSOBSTYPE float
#define MPIENSOBSTYPE MPI_FLOAT

#define STATE_BIGNUM 1.0e3
#define STD_BIG 1.0e10

#define ALLOCTYPE_1D 0
#define ALLOCTYPE_2D 1
#define ALLOCTYPE_3D 2

#define UPDATE_DOFIELDS       (1 << 0)
#define UPDATE_DOSPREAD       (1 << 1)
#define UPDATE_DOINFLATION    (1 << 2)
#define UPDATE_DOPLOGS        (1 << 3)
#define UPDATE_LEAVETILES     (1 << 4)
#define UPDATE_OUTPUTINC      (1 << 5)
#define UPDATE_SEPARATEOUTPUT (1 << 6)
#define UPDATE_DIRECTWRITE    (1 << 7)
#define UPDATE_DEFAULT        (UPDATE_DOFIELDS | UPDATE_DOPLOGS | UPDATE_SEPARATEOUTPUT)

/*
 * the default vertical split for 3D fields in obs stats
 */
#define DEPTH_SHALLOW 50.0
#define DEPTH_DEEP 500.0
#define DEPTH_MAX 9999.0

/* 
 *it is assumed that if |value| > MAXOBSVAL, then it is filled with the
 * missing value
 */
#define MAXOBSVAL 900.0

#define REARTH 6371.0
#define DEG2RAD (M_PI / 180.0)
#define TWOPI (M_PI * 2.0)

#if defined(MPI)
#include <mpi.h>
#endif
extern int nprocesses;
extern int rank;

extern int enkf_obstype;
extern int enkf_exitaction;
extern int enkf_verbose;
extern int enkf_separateout;
extern int enkf_directwrite;
extern int enkf_nomeanupdate;
extern int enkf_fstatsonly;
extern int enkf_noobsdatecheck;
extern int enkf_considersubgridvar;

#define _DEFINITIONS_H
#endif
