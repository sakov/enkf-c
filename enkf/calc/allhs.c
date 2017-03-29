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
    char* obstypekind;
    char* mappingname;
    H_fn H;
} H_entry;

H_entry allhentries[] = {
    {"SURFACE", "standard", H_surf_standard},
    {"SURFACE", "biased", H_surf_biased},
    {"SUBSURFACE", "standard", H_subsurf_standard},
    {"SUBSURFACE", "wsurfbias", H_subsurf_wsurfbias}
};

/** List H-functions for all observation types (obstypename = NULL) or for
 ** observations of particular type specified by `obstypename'.
 */
void describe_hentries(char* obstypekind)
{
    int nhentries = sizeof(allhentries) / sizeof(H_entry);
    int i;

    if (obstypekind == NULL) {
        enkf_printf("  Available H functions:\n");
        enkf_printf("    obs. type kind    H function\n");
        enkf_printf("    ----------------------------\n");
        for (i = 0; i < nhentries; ++i)
            enkf_printf("     %10s     %s\n", allhentries[i].obstypekind, allhentries[i].mappingname);
    } else {
        enkf_printf("  Available H functions for observation kind \"%s\":\n", obstypekind);
        for (i = 0; i < nhentries; ++i)
            if (strcmp(obstypekind, allhentries[i].obstypekind) == 0)
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
        if (strcmp(allhentries[i].obstypekind, ot->kind) == 0 && strcmp(allhentries[i].mappingname, mappingname) == 0)
            return allhentries[i].H;

    enkf_printf("\n\n  ERROR: no H function \"%s\" for observation kind \"%s\"\n\n", mappingname, ot->kind);
    describe_hentries(ot->kind);
    enkf_quit("getH(): bailing out");
    return NULL;
}
