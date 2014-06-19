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
    {"SLA", "standard", H_sla_standard},
    {"TEM", "standard", H_subsurf_standard},
    {"SAL", "standard", H_subsurf_standard}
};

/**
 */
void describe_hentries(void)
{
    int nhentries = sizeof(allhentries) / sizeof(H_entry);
    int i;

    enkf_printf("  Available H functions:\n");
    enkf_printf("    obs. type  function\n");
    enkf_printf("    -------------------\n");
    for (i = 0; i < nhentries; ++i)
        enkf_printf("     %5s     %s\n", allhentries[i].obstypename, allhentries[i].mappingname);
}

/**
 */
H_fn getH(char obstypename[], char mappingname[])
{
    int nhentries = sizeof(allhentries) / sizeof(H_entry);
    int i;

    for (i = 0; i < nhentries; ++i)
        if (strcmp(allhentries[i].obstypename, obstypename) == 0 && strcmp(allhentries[i].mappingname, mappingname) == 0)
            return allhentries[i].H;

    enkf_printf("\n\n  ERROR: no H function \"%s\" for observation type \"%s\"\n\n", mappingname, obstypename);
    describe_hentries();
    enkf_quit("bailing out");

    return NULL;
}
