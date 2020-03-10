/******************************************************************************
 *
 * File:        global.c        
 *
 * Created:     19/12/2013
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description: Stores global variables and flags
 *
 * Revisions:   06032020 PS: moved here communicators etc. from struct dasystem
 *
  *****************************************************************************/

#include "definitions.h"

int nprocesses = 1;
int rank = 0;

#if defined(HE_VIASHMEM)
/*
 * "sm" below stands for "shared memory". The shared memory is allocated
 * on each node to hold HE (S) and (HE)^T (S^T) objects.
 */
MPI_Comm sm_comm;
int sm_comm_rank = -1;
int* sm_comm_ranks = NULL;

/*
 * The node communicator includes the first core on each node. It is created
 * to gather S and S^T.
 */
MPI_Comm node_comm = MPI_COMM_WORLD;
int node_comm_rank = -1;
int node_comm_size = 0;
int* node_comm_ranks = NULL;
#endif

int enkf_obstype = OBSTYPE_VALUE;
int enkf_exitaction = EXITACTION_DEFAULT;
int enkf_verbose = 1;
int enkf_nomeanupdate = 0;
int enkf_fstatsonly = 0;
int enkf_noobsdatecheck = 0;
int enkf_considersubgridvar = 0;
