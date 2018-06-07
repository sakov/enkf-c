/******************************************************************************
 *
 * File:        enkf_update.c        
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
#include <math.h>
#include "version.h"
#include "definitions.h"
#include "utils.h"
#include "enkfprm.h"
#include "dasystem.h"

/**
 */
static void usage()
{
    enkf_printf("  Usage: enkf_update <prm file> [<options>]\n");
    enkf_printf("  Options:\n");
    enkf_printf("  --calculate-spread\n");
    enkf_printf("      calculate ensemble spread and write to %s\n", FNAME_SPREAD);
    enkf_printf("  --calculate-forecast-spread\n");
    enkf_printf("      calculate forecast ensemble spread only and write to %s\n", FNAME_SPREAD);
    enkf_printf("  --calculate-vertical-correlations\n");
    enkf_printf("      calculate correlation coefficients between surface and other layers of\n");
    enkf_printf("      3D variables and write to %s\n", FNAME_VERTCORR);
    enkf_printf("  --describe-prm-format [main|model|grid]\n");
    enkf_printf("      describe format of a parameter file and exit\n");
    enkf_printf("  --direct-write\n");
    enkf_printf("      write fields directly to the output file (default: write to tiles first)\n");
    enkf_printf("  --joint-output\n");
    enkf_printf("      append analyses to forecast files (default: write to separate files)\n");
    enkf_printf("  --leave-tiles\n");
    enkf_printf("      do not delete tiles\n");
    enkf_printf("  --no-fields-write\n");
    enkf_printf("      do not write analysis fields (only point logs and/or spread)\n");
    enkf_printf("  --output-increment\n");
    enkf_printf("      output analysis increment (default: output analysis)\n");
    enkf_printf("  --write-inflation\n");
    enkf_printf("      write adaptive inflation magnitudes to %s\n", FNAME_INFLATION);
    enkf_printf("  --version\n");
    enkf_printf("      print version and exit\n");

    exit(0);
}

/**
 */
static void parse_commandline(int argc, char* argv[], char** fname, int* updatespec)
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
            *updatespec |= (UPDATE_DOFORECASTSPREAD | UPDATE_DOANALYSISSPREAD);
            i++;
            continue;
        } else if (strcmp(argv[i], "--calculate-forecast-spread") == 0) {
            *updatespec |= UPDATE_DOFORECASTSPREAD;
            i++;
            continue;
        } else if (strcmp(argv[i], "--calculate-vertical-correlations") == 0) {
            *updatespec |= UPDATE_DOVERTCORRS;
            i++;
        } else if (strcmp(argv[i], "--describe-prm-format") == 0) {
            if (i < argc - 1) {
                if (strcmp(argv[i + 1], "main") == 0)
                    enkfprm_describeprm();
                else if (strcmp(argv[i + 1], "model") == 0)
                    model_describeprm();
                else if (strcmp(argv[i + 1], "grid") == 0)
                    grid_describeprm();
                else
                    usage();
            } else
                enkfprm_describeprm();
            exit(0);
        } else if (strcmp(argv[i], "--direct-write") == 0) {
            *updatespec |= UPDATE_DIRECTWRITE;
            i++;
            continue;
        } else if (strcmp(argv[i], "--joint-output") == 0) {
            *updatespec &= ~UPDATE_SEPARATEOUTPUT;
            i++;
            continue;
        } else if (strcmp(argv[i], "--leave-tiles") == 0) {
            *updatespec |= UPDATE_LEAVETILES;
            i++;
            continue;
        } else if (strcmp(argv[i], "--no-fields-write") == 0) {
            *updatespec &= ~UPDATE_DOFIELDS;
            i++;
            continue;
        } else if (strcmp(argv[i], "--output-increment") == 0) {
            *updatespec |= UPDATE_OUTPUTINC;
            i++;
            continue;
        } else if (strcmp(argv[i], "--separate-output") == 0) {
            /*
             * (historic option, became default, do nothing)
             */
            i++;
            continue;
        } else if (strcmp(argv[i], "--write-inflation") == 0) {
            *updatespec |= UPDATE_DOINFLATION;
            i++;
            continue;
        } else if (strcmp(argv[i], "--version") == 0) {
            enkf_printversion();
            exit(0);
        } else
            enkf_quit("command line: option \"%s\" not recognised", argv[i]);
    }

    if (*fname == NULL)
        enkf_quit("command line: parameter file not specified");
}

/**
 */
static void describe_updatespec(int updatespec)
{
    enkf_printf("  update specs:\n");
    enkf_printf("    do model fields  = %s\n", (updatespec & UPDATE_DOFIELDS) ? "[+]" : "[-]");
    enkf_printf("    do spread        = %s\n", (updatespec & UPDATE_DOSPREAD) ? "[+]" : "[-]");
    enkf_printf("    do pointlogs     = %s\n", (updatespec & UPDATE_DOSPREAD) ? "[+]" : "[-]");
    if (updatespec & UPDATE_DIRECTWRITE)
        enkf_printf("    direct write     = [+]\n");
    if (!(updatespec & UPDATE_DIRECTWRITE) && updatespec & UPDATE_LEAVETILES)
        enkf_printf("    leave tiles      = [+]\n");
    if (updatespec & UPDATE_DOFIELDS) {
        if (updatespec & UPDATE_OUTPUTINC)
            enkf_printf("    output increment = [+]\n");
        enkf_printf("    separate output  = %s\n", (updatespec & UPDATE_SEPARATEOUTPUT) ? "[+]" : "[-]");
    }
}

/**
 */
int main(int argc, char* argv[])
{
    char* fname_prm = NULL;
    int updatespec = UPDATE_DEFAULT;
    enkfprm* prm = NULL;
    dasystem* das = NULL;

    enkf_init(&argc, &argv);

    parse_commandline(argc, argv, &fname_prm, &updatespec);

    enkf_printf("  running UPDATE for EnKF-C version %s:\n", ENKF_VERSION);
    print_commandinfo(argc, argv);
    enkf_printtime("  ");

    enkf_printf("  reading system specs from \"%s\":\n", fname_prm);
    prm = enkfprm_read(fname_prm);
    enkfprm_print(prm, "    ");

    if (nprocesses == 1 && !(updatespec & UPDATE_DIRECTWRITE)) {
        enkf_printf("  nproc = 1 -> using direct write\n");
        updatespec |= UPDATE_DIRECTWRITE;
    }
    if (updatespec & UPDATE_DOINFLATION) {
        if (prm->mode != MODE_ENKF) {
            enkf_printf("  prm: mode = EnOI -> cancelling writing of inflation\n");
            updatespec &= ~UPDATE_DOINFLATION;
        } else if (isnan(prm->inf_ratio)) {
            enkf_printf("  prm: inflation mode = plain -> cancelling writing of inflation\n");
            updatespec &= ~UPDATE_DOINFLATION;
        }
    }
    if (prm->nplogs == 0)
        updatespec &= ~UPDATE_DOPLOGS;

    describe_updatespec(updatespec);
    if ((updatespec & (UPDATE_DOFIELDS | UPDATE_DOSPREAD | UPDATE_DOPLOGS | UPDATE_DOINFLATION | UPDATE_DOVERTCORRS)) == 0)
        enkf_quit("nothing to do");

    enkf_printf("  initialising the system:\n");
    das = das_create(prm);
    enkfprm_destroy(prm);
    das->updatespec = updatespec;

    if (updatespec & (UPDATE_DOFIELDS | UPDATE_DOSPREAD | UPDATE_DOPLOGS | UPDATE_DOINFLATION)) {
        if (das->mode == MODE_ENKF)
            enkf_printf("  updating the ensemble:\n");
        else if (das->mode == MODE_ENOI)
            enkf_printf("  updating the model state:\n");
        das_update(das);
        enkf_flush();
    }

    if (das->updatespec & UPDATE_DOVERTCORRS)
        das_writevcorrs(das);

    das_destroy(das);

    enkf_finish();

    return 0;
}
