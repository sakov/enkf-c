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

#include <stdlib.h>
#include "definitions.h"

int nprocesses = 1;
int rank = 0;

#if defined(USE_SHMEM)
/*
 * "sm" below stands for "shared memory". The shared memory is allocated
 * on each node to hold HE (S) and (HE)^T (S^T) objects. From v1.109.0 and
 * v1.110.0 it is also used for holding large kd-trees.
 */
MPI_Comm sm_comm;
int sm_comm_rank = -1;
int sm_comm_size = 0;
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
int enkf_fstatsonly = 0;
int enkf_noobsdatecheck = 0;
int enkf_considersubgridvar = 0;
int enkf_doplogs = 1;
int enkf_allowenoilog = 0;
