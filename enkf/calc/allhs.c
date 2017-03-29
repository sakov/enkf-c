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
    int issurface;
    char* H_tag;
    H_fn H;
} H_entry;

H_entry allhentries[] = {
    {1, "standard", H_surf_standard},
    {1, "biased", H_surf_biased},
    {0, "standard", H_subsurf_standard},
    {0, "wsurfbias", H_subsurf_wsurfbias}
};

/** List H-functions for all observation types (obstypename = NULL) or for
 ** observations of particular type specified by `obstypename'.
 */
void describe_hentries(int issurface)
{
    int nhentries = sizeof(allhentries) / sizeof(H_entry);
    int i;

    if (issurface == -1) {
        enkf_printf("  Available H functions:\n");
        enkf_printf("    surface flag    H function\n");
        enkf_printf("    ----------------------------\n");
        for (i = 0; i < nhentries; ++i)
            enkf_printf("     %5d     %s\n", allhentries[i].issurface, allhentries[i].H_tag);
    } else {
        enkf_printf("  Available H functions for %s observations:\n", (issurface) ? "surface" : "subsurface");
        for (i = 0; i < nhentries; ++i)
            if (issurface == allhentries[i].issurface)
                enkf_printf("    %s\n", allhentries[i].H_tag);
    }
}

/**
 */
H_fn getH(obstype* ot, char H_tag[])
{
    int nhentries = sizeof(allhentries) / sizeof(H_entry);
    int i;

    for (i = 0; i < nhentries; ++i)
        if (allhentries[i].issurface == ot->issurface && strcmp(allhentries[i].H_tag, H_tag) == 0)
            return allhentries[i].H;

    enkf_printf("\n\n  ERROR: no H function \"%s\" for %s observations\n\n", H_tag, (ot->issurface) ? "surface" : "subsurface");
    describe_hentries(ot->issurface);
    enkf_quit("getH(): bailing out");
    return NULL;
}
