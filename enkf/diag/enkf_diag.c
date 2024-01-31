/******************************************************************************
 *
 * File:        enkf_diag.c        
 *
 * Created:     30/01/2024
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description: Performs various diagnostics of the ensemble.
 *
 * Revisions:
 *
 *****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "version.h"
#include "definitions.h"
#include "utils.h"
#include "enkfprm.h"
#include "dasystem.h"
#include "diags.h"

#define NINC 10

/**
 */
static void usage(void)
{
    enkf_printf("  Usage: enkf_diag <prm file> [<options>]\n");
    enkf_printf("  Options:\n");
    enkf_printf("  --calculate-spread\n");
    enkf_printf("      calculate ensemble spread and write to %s\n", FNAME_SPREAD);
    enkf_printf("  --calculate-vertical-correlations\n");
    enkf_printf("      calculate correlation coefficients between surface and other layers of\n");
    enkf_printf("      3D variables and write to %s\n", FNAME_VERTCORR);
    enkf_printf("  --calculate-vertical-correlations-with <varname1> <layer1>\n");
    enkf_printf("    [<varname2> <layer2>] [...]\n");
    enkf_printf("      calculate correlation coefficients between specified layer of a 3D\n");
    enkf_printf("      variable and all layers of 3D variables on the same horizontal grid\n");
    enkf_printf("      and write to %s-<varname>-<layer>.nc\n", FNAMEPREFIX_VERTCORRWITH);
    enkf_printf("  --calculate-vertical-covariances-with <varname1> <layer1>\n");
    enkf_printf("    [<varname2> <layer2>] [...]\n");
    enkf_printf("      calculate covariances between specified layer of a 3D\n");
    enkf_printf("      variable and all layers of 3D variables on the same horizontal grid\n");
    enkf_printf("      and write to %s-<varname>-<layer>.nc\n", FNAMEPREFIX_VERTCOVWITH);
    enkf_printf("  --version\n");
    enkf_printf("      print version and exit\n");

    exit(0);
}

/**
 */
static void parse_commandline(int argc, char* argv[], char** fname, int* dospread, int* dovcorrs, int* nvcorrwith, char*** vcorrwith, int** kvcorrwith, int* nvcovwith, char*** vcovwith, int** kvcovwith)
{
    int i;
    
    if (argc < 2)
        usage();

    i = 1;
    while (i < argc) {
        if (argv[i][0] != '-') {
            if (*fname == NULL) {
                *fname = argv[i];
                i++;
                continue;
            } else
                usage();
        } else if (strcmp(argv[i], "--calculate-spread") == 0) {
            *dospread = 1;
            i++;
        } else if (strcmp(argv[i], "--calculate-vertical-correlations") == 0) {
            *dovcorrs = 1;
            i++;
        } else if (strcmp(argv[i], "--calculate-vertical-correlations-with") == 0) {
            i++;
            while (i < argc && argv[i][0] != '-') {
                if (*nvcorrwith % NINC == 0) {
                    *vcorrwith = realloc(*vcorrwith, (*nvcorrwith + NINC) * sizeof(void*));
                    *kvcorrwith = realloc(*kvcorrwith, (*nvcorrwith + NINC) * sizeof(int));
                }
                (*vcorrwith)[*nvcorrwith] = argv[i];
                i++;
                if (argv[i][0] == '-')
                    enkf_quit("parse_commandline(): expected layer ID after variable name \"%s\"", (*vcorrwith)[*nvcorrwith]);
                if (!str2int(argv[i], &(*kvcorrwith)[*nvcorrwith]))
                    enkf_quit("parse_commandline(): could no convert \"%s\" to int", argv[i][0]);
                i++;
                (*nvcorrwith)++;
            }
            if (*nvcorrwith == 0)
                enkf_quit("parse_commandline(): no variable name and layer ID found after \"--calculate-vertical-correlations-with\"");
        } else if (strcmp(argv[i], "--calculate-vertical-covariances-with") == 0) {
            i++;
            while (i < argc && argv[i][0] != '-') {
                if (*nvcovwith % NINC == 0) {
                    *vcovwith = realloc(*vcovwith, (*nvcovwith + NINC) * sizeof(void*));
                    *kvcovwith = realloc(*kvcovwith, (*nvcovwith + NINC) * sizeof(int));
                }
                (*vcovwith)[*nvcovwith] = argv[i];
                i++;
                if (argv[i][0] == '-')
                    enkf_quit("parse_commandline(): expected layer ID after variable name \"%s\"", (*vcovwith)[*nvcovwith]);
                if (!str2int(argv[i], &(*kvcovwith)[*nvcovwith]))
                    enkf_quit("parse_commandline(): could no convert \"%s\" to int", argv[i][0]);
                i++;
                (*nvcovwith)++;
            }
            if (*nvcovwith == 0)
                enkf_quit("parse_commandline(): no variable name and layer ID found after \"--calculate-vertical-covariances-with\"");
        } else if (strcmp(argv[i], "--version") == 0) {
            enkf_printversion();
            exit(0);
        } else
            enkf_quit("parse_commandline(): option \"%s\" not recognised", argv[i]);
    }
    
    if (*fname == NULL)
        enkf_quit("parse_commandline(): parameter file not specified");
}

/**
 */
int main(int argc, char* argv[])
{
    char* fname_prm = NULL;
    int dospread = 0;
    int dovcorrs = 0;
    int nvcorrwith = 0;
    char** vcorrwith = NULL;
    int* kvcorrwith = NULL;
    int nvcovwith = 0;
    char** vcovwith = NULL;
    int* kvcovwith = NULL;
    enkfprm* prm = NULL;
    dasystem* das = NULL;
    int i;

    enkf_init(&argc, &argv);

    parse_commandline(argc, argv, &fname_prm, &dospread, &dovcorrs, &nvcorrwith, &vcorrwith, &kvcorrwith, &nvcovwith, &vcovwith, &kvcovwith);

    enkf_printf("  running DIAG for EnKF-C version %s:\n", ENKF_VERSION);
    print_commandinfo(argc, argv);
    enkf_printtime("  ");

    if (dospread == 0 && dovcorrs == 0 && nvcorrwith == 0 && nvcovwith == 0) {
        enkf_printf("  nothing to do\n");
        goto finish;
    }

    enkf_printf("  reading system specs from \"%s\":\n", fname_prm);
    prm = enkfprm_read(fname_prm);
    enkfprm_print(prm, "    ");

    enkf_printf("  initialising the system:\n");
    das = das_create(prm);
    enkfprm_destroy(prm);

    if (dospread)
        das_writespread(das);

    if (dovcorrs)
        das_writevcorrs(das);

    for (i = 0; i < nvcorrwith; ++i)
        das_writevcorrs_with(das, vcorrwith[i], kvcorrwith[i], 1);

    for (i = 0; i < nvcovwith; ++i)
        das_writevcorrs_with(das, vcovwith[i], kvcovwith[i], 0);

    das_destroy(das);
    if (nvcorrwith > 0) {
        free(vcorrwith);
        free(kvcorrwith);
    }
    if (nvcovwith > 0) {
        free(vcovwith);
        free(kvcovwith);
    }

 finish:
    enkf_finish();

    return 0;
}
