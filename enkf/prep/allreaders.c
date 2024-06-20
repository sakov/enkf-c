/******************************************************************************
 *
 * File:        allreaders.c        
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "stringtable.h"
#include "kdtree.h"
#include "definitions.h"
#include "utils.h"
#include "obsprm.h"
#include "model.h"
#include "grid.h"
#include "observations.h"
#include "prep_utils.h"
#include "allreaders.h"

typedef struct {
    char* product;
    char* reader;
    obsread_fn readfn;
    void (*describefn) (void);
} obsreader_entry;

obsreader_entry allreaders[] = {
    {"ANY", "scattered", reader_scattered, reader_scattered_describe},
    {"ANY", "gridded_xy", reader_gridded_xy, reader_gridded_xy_describe},
    {"ANY", "gridded_xyz", reader_gridded_xyz, reader_gridded_xyz_describe},
    {"ANY", "gridded_xyh", reader_gridded_xyh, reader_gridded_xyh_describe},
    {"ANY", "z", reader_z, reader_z_describe},
    {"NAVO", "navo", reader_navo, reader_navo_describe},
    {"WINDSAT", "windsat", reader_windsat, reader_windsat_describe},
    {"MMT", "mmt", reader_mmt, reader_mmt_describe},
    {"AMSR2", "amsr2", reader_amsr2, reader_amsr2_describe},
    {"AMSRE", "amsre", reader_amsre, reader_amsre_describe},
    {"CMEMS", "cmems", reader_cmems, reader_cmems_describe},
    {"EN4", "en4", reader_en4, reader_en4_describe}
};

/**
 */
void list_readers(void)
{
    int nreaders = sizeof(allreaders) / sizeof(obsreader_entry);
    int i;

    enkf_printf("  generic readers:\n");
    for (i = 0; i < nreaders; ++i)
        if (strcmp(allreaders[i].product, "ANY") == 0)
            enkf_printf("    %s\n", allreaders[i].reader);
    enkf_printf("  custom readers:\n");
    for (i = 0; i < nreaders; ++i)
        if (strcmp(allreaders[i].product, "ANY") != 0)
            enkf_printf("    %s\n", allreaders[i].reader);
}

/**
 */
static void describe_readers(void)
{
    int nreaders = sizeof(allreaders) / sizeof(obsreader_entry);
    int i;
    const char anystr[] = "(ANY)";

    enkf_printf("  Available readers:\n");
    enkf_printf("    product   reader\n");
    enkf_printf("    ----------------\n");
    for (i = 0; i < nreaders; ++i)
        if (strcmp(allreaders[i].product, "ANY") == 0)
            enkf_printf("    %-9s %s\n", anystr, allreaders[i].reader);
        else
            enkf_printf("    %-9s %s\n", allreaders[i].product, allreaders[i].reader);
}

/**
 */
obsread_fn get_obsreadfn(obsmeta* m)
{
    int nreaders = sizeof(allreaders) / sizeof(obsreader_entry);
    int i;

    for (i = 0; i < nreaders; ++i)
        if ((strcmp(allreaders[i].product, "ANY") == 0 || strcasecmp(m->product, allreaders[i].product) == 0) && strcmp(m->reader, allreaders[i].reader) == 0)
            return allreaders[i].readfn;

    enkf_printf("\n\n  ERROR: no observation reader \"%s\" for product \"%s\"\n\n", m->reader, m->product);
    describe_readers();
    enkf_quit("bailing out", m->product);
    return NULL;
}

/**
 */
void describe_reader(char readertag[])
{
    int nreaders = sizeof(allreaders) / sizeof(obsreader_entry);
    int i;

    for (i = 0; i < nreaders; ++i)
        if (strcmp(allreaders[i].reader, readertag) == 0) {
            allreaders[i].describefn();
            return;
        }
    enkf_printf("  describe_reader(): reader \"%s\" not found\n", readertag);
}
