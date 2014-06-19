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

#define MAXSTRLEN 2048
#define NOBS_INC 200000

#define BASEYEAR 1990
#define BASEMONTH 1
#define BASEDAY 1

#define FNAME_OBS "observations-orig.nc"
#define FNAME_SOBS "observations.nc"
#define FNAME_X5 "X5.nc"
#define FNAME_W "w.nc"
#define FNAME_ENKFSTATS "enkf_diag.nc"
#define FNAME_SPREAD "spread.nc"
#define FNAME_HE "HE.nc"

#define STATUS_OK 0             /* do not change */
#define STATUS_OUTSIDE 1
#define STATUS_LAND 2
#define STATUS_SHALLOW 3
#define STATUS_RANGE 4

#define LONTYPE_NONE 0
#define LONTYPE_180 1
#define LONTYPE_360 2

#define STDTYPE_NONE 0
#define STDTYPE_VALUE 1
#define STDTYPE_FILE 2

#define MODE_NONE 0
#define MODE_ENKF 1
#define MODE_ENOI 2

#define TARGET_NONE 0
#define TARGET_ANALYSIS 1
#define TARGET_INCREMENT 2
#define TARGET_DEFAULT TARGET_ANALYSIS

#define SCHEME_NONE 0
#define SCHEME_DENKF 1
#define SCHEME_ETKF 2
#define SCHEME_DEFAULT SCHEME_DENKF

#define EXITACTION_BACKTRACE 0
#define EXITACTION_SEGFAULT 1
#define EXITACTION_DEFAULT EXITACTION_BACKTRACE
extern int enkf_exitaction;

#define OBSTYPE_VALUE 0
#define OBSTYPE_INNOVATION 1
extern int enkf_obstype;

#define ENSOBSTYPE float
#define MPIENSOBSTYPE MPI_FLOAT

#define STATE_BIGNUM 1.0e3

/* it is assumed that if |value| > MAXOBSVAL, then it is filled with the
 * missing value */
#define MAXOBSVAL 900.0

#define DEG2RAD (M_PI / 180.0)
#define REARTH 6371.0
#define TWOPI (M_PI * 2.0)

#if defined(MPI)
#include <mpi.h>
#endif
extern int nprocesses;
extern int rank;

extern int enkf_separateout;
extern int enkf_directwrite;
extern int enkf_nomeanupdate;

#define _DEFINITIONS_H
#endif
