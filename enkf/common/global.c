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
 * Revisions:
  *****************************************************************************/

#include "definitions.h"

int nprocesses = 1;
int rank = 0;

int enkf_obstype = OBSTYPE_VALUE;
int enkf_exitaction = EXITACTION_DEFAULT;
int enkf_verbose = 1;
int enkf_separateout = 0;
int enkf_directwrite = 0;
int enkf_nomeanupdate = 0;
int enkf_fstatsonly = 0;
