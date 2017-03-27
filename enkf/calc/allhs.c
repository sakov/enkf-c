/******************************************************************************
 *
 * File:        allhs.c        
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "definitions.h"
#include "utils.h"
#include "model2obs.h"
#include "allhs.h"

typedef struct {
    char* obstypename;
    char* mappingname;
    H_fn H;
} H_entry;

H_entry allhentries[] = {
    {"SST", "standard", H_surf_standard},
    {"SST", "biased", H_surf_biased},
    {"SLA", "standard", H_surf_standard},
    {"SLA", "biased", H_surf_biased},
    {"TEM", "standard", H_subsurf_standard},
    {"TEM", "wsurfbias", H_subsurf_wsurfbias},
    {"SAL", "standard", H_subsurf_standard},
    {"SSS", "standard", H_surf_standard}
};

/** List H-functions for all observation types (obstypename = NULL) or for
 ** observations of particular type specified by `obstypename'.
 */
void describe_hentries(char* obstypename)
{
    int nhentries = sizeof(allhentries) / sizeof(H_entry);
    int i;

    if (obstypename == NULL) {
        enkf_printf("  Available H functions:\n");
        enkf_printf("    obs. type  function\n");
        enkf_printf("    -------------------\n");
        for (i = 0; i < nhentries; ++i)
            enkf_printf("     %5s     %s\n", allhentries[i].obstypename, allhentries[i].mappingname);
    } else {
        enkf_printf("  Available H functions for %s:\n", obstypename);
        for (i = 0; i < nhentries; ++i)
            if (strcmp(obstypename, allhentries[i].obstypename) == 0)
                enkf_printf("    %s\n", allhentries[i].mappingname);
    }
}

/**
 */
H_fn getH(obstype* ot, char mappingname[])
{
    int nhentries = sizeof(allhentries) / sizeof(H_entry);
    int i;

    for (i = 0; i < nhentries; ++i)
        if ((strcmp(allhentries[i].obstypename, ot->name) == 0 && strcmp(allhentries[i].mappingname, mappingname) == 0) || (ot->aliasname != NULL && strcmp(allhentries[i].obstypename, ot->aliasname) == 0 && strcmp(allhentries[i].mappingname, mappingname) == 0))
            return allhentries[i].H;

    enkf_printf("\n\n  ERROR: no H function \"%s\" for observation type \"%s\"\n\n", mappingname, ot->name);
    if (ot->aliasname != NULL)
        enkf_printf("    or for the aliased observation type \"%s\"\n", ot->aliasname);
    describe_hentries(ot->name);
    enkf_quit("getH(): bailing out");
    return NULL;
}
