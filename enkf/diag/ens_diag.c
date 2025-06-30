/******************************************************************************
 *
 * File:        ens_diag.c        
 *
 * Created:     30/01/2024
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description: Performs various diagnostics of the ensemble.
 *
 * Revisions:   30/06/2024 Renamed from enkf_diag.c to ens_diag.c.
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
    enkf_printf("  Usage: ens_diag <prm file> [<options>]\n");
    enkf_printf("  Options:\n");
    enkf_printf("  --calculate-spread\n");
    enkf_printf("      calculate ensemble spread and write to %s\n", FNAME_SPREAD);
    enkf_printf("  --calculate-vertical-correlations\n");
    enkf_printf("      calculate correlation coefficients between surface and other layers of\n");
    enkf_printf("      3D variables and write to %s\n", FNAME_VERTCORR);
    enkf_printf("  --calculate-vertical-correlations-with <varname1> <layer1>\n");
    enkf_printf("    [<varname2> <layer2>] [...]\n");
    enkf_printf("      calculate correlation coefficients between specified field (a layer of\n");
    enkf_printf("      a variable) and all other fields on the same horizontal grid and\n");
    enkf_printf("      write to %s-<varname>-<layer>.nc\n", FNAMEPREFIX_VERTCORRWITH);
    enkf_printf("  --calculate-vertical-covariances-with <varname1> <layer1>\n");
    enkf_printf("    [<varname2> <layer2>] [...]\n");
    enkf_printf("      calculate covariances between specified field and all other fields\n");
    enkf_printf("      on the same horizontal grid and write to %s-<varname>-<layer>.nc\n", FNAMEPREFIX_VERTCOVWITH);
    enkf_printf("  --calculate-vertical-sensitivities-with <varname1> <layer1>\n");
    enkf_printf("    [<varname2> <layer2>] [...]\n");
    enkf_printf("      calculate sensitivities between specified field and all other fields\n");
    enkf_printf("      on the same horizontal grid and write to %s-<varname>-<layer>.nc\n", FNAMEPREFIX_VERTCOVWITH);
    enkf_printf("  --describe-prm-format [main|model|grid]\n");
    enkf_printf("      describe format of a parameter file and exit\n");
    enkf_printf("  --version\n");
    enkf_printf("      print version and exit\n");

    exit(0);
}

/**
 */
static void parse_commandline(int argc, char* argv[], char** fname, int* dospread, int* dovcorrs, int* nvcorrwith, char*** vcorrwith, int** kvcorrwith, int* nvcovwith, char*** vcovwith, int** kvcovwith, int* nvsenswith, char*** vsenswith, int** kvsenswith)
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
                if (i == argc || argv[i][0] == '-')
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
                if (i == argc || argv[i][0] == '-')
                    enkf_quit("parse_commandline(): expected layer ID after variable name \"%s\"", (*vcovwith)[*nvcovwith]);
                if (!str2int(argv[i], &(*kvcovwith)[*nvcovwith]))
                    enkf_quit("parse_commandline(): could no convert \"%s\" to int", argv[i][0]);
                i++;
                (*nvcovwith)++;
            }
            if (*nvcovwith == 0)
                enkf_quit("parse_commandline(): no variable name and layer ID found after \"--calculate-vertical-covariances-with\"");
        } else if (strcmp(argv[i], "--calculate-vertical-sensitivities-with") == 0) {
            i++;
            while (i < argc && argv[i][0] != '-') {
                if (*nvsenswith % NINC == 0) {
                    *vsenswith = realloc(*vsenswith, (*nvsenswith + NINC) * sizeof(void*));
                    *kvsenswith = realloc(*kvsenswith, (*nvsenswith + NINC) * sizeof(int));
                }
                (*vsenswith)[*nvsenswith] = argv[i];
                i++;
                if (i == argc || argv[i][0] == '-')
                    enkf_quit("parse_commandline(): expected layer ID after variable name \"%s\"", (*vsenswith)[*nvsenswith]);
                if (!str2int(argv[i], &(*kvsenswith)[*nvsenswith]))
                    enkf_quit("parse_commandline(): could no convert \"%s\" to int", argv[i][0]);
                i++;
                (*nvsenswith)++;
            }
            if (*nvsenswith == 0)
                enkf_quit("parse_commandline(): no variable name and layer ID found after \"--calculate-vertical-sensitivities-with\"");
        } else if (strcmp(argv[i], "--version") == 0) {
            enkf_printversion();
            exit(0);
        } else if (strcmp(argv[i], "--describe-prm-format") == 0) {
            if (i < argc - 1) {
                if (strcmp(argv[i + 1], "main") == 0)
                    enkfprm_describeprm_ensdiag();
                else if (strcmp(argv[i + 1], "model") == 0)
                    model_describeprm();
                else if (strcmp(argv[i + 1], "grid") == 0)
                    grid_describeprm();
                else
                    usage();
            } else
                enkfprm_describeprm_ensdiag();
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
    int nvsenswith = 0;
    char** vsenswith = NULL;
    int* kvsenswith = NULL;
    enkfprm* prm = NULL;
    dasystem* das = NULL;
    int i;

    parse_commandline(argc, argv, &fname_prm, &dospread, &dovcorrs, &nvcorrwith, &vcorrwith, &kvcorrwith, &nvcovwith, &vcovwith, &kvcovwith, &nvsenswith, &vsenswith, &kvsenswith);

    enkf_init(&argc, &argv);

    enkf_printf("  running DIAG for EnKF-C version %s:\n", ENKF_VERSION);
    print_commandinfo(argc, argv);
    enkf_printtime("  ");

    if (dospread == 0 && dovcorrs == 0 && nvcorrwith == 0 && nvcovwith == 0 && nvsenswith == 0) {
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
        das_writevcorrs_with(das, vcorrwith[i], kvcorrwith[i], CALC_CORR);

    for (i = 0; i < nvcovwith; ++i)
        das_writevcorrs_with(das, vcovwith[i], kvcovwith[i], CALC_COV);

    for (i = 0; i < nvsenswith; ++i)
        das_writevcorrs_with(das, vsenswith[i], kvsenswith[i], CALC_SENS);

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
