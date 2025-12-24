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
#include "pointlog.h"

/**
 */
static void usage(void)
{
    enkf_printf("  Usage: enkf_update <prm file> [<options>]\n");
    enkf_printf("  Options:\n");
    enkf_printf("  --allow-logspace-with-static-ens\n");
    enkf_printf("      confirm that static ensemble is conditioned for using log space\n");
    enkf_printf("  --calculate-spread\n");
    enkf_printf("      calculate forecast ensemble spread and write to %s\n", FNAME_SPREAD);
    enkf_printf("  --describe-prm-format [main|model|grid]\n");
    enkf_printf("      describe format of a parameter file and exit\n");
    enkf_printf("  --direct-write\n");
    enkf_printf("      write fields directly to the output file (default: write to tiles first)\n");
    enkf_printf("  --no-fields-write\n");
    enkf_printf("      do not write analysis fields (only diagnostic data)\n");
    enkf_printf("  --no-update\n");
    enkf_printf("      exclude tasks that require ensemble update\n");
    enkf_printf("  --output-increment\n");
    enkf_printf("      output analysis increment (default: output analysis)\n");
    enkf_printf("  --write-inflation\n");
    enkf_printf("      write capped inflation magnitudes to %s\n", FNAME_INFLATION);
    enkf_printf("  --version\n");
    enkf_printf("      print version and exit\n");

    exit(0);
}

/**
 */
static void parse_commandline(int argc, char* argv[], char** fname, int* updatespec)
{
    int no_update = 0;
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
        } else if (strcmp(argv[i], "--allow-logspace-with-static-ens") == 0) {
            enkf_allowenoilog = 1;
            i++;
            continue;
        } else if (strcmp(argv[i], "--calculate-spread") == 0) {
            *updatespec |= (UPDATE_DOSPREAD);
            i++;
            continue;
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
        } else if (strcmp(argv[i], "--no-fields-write") == 0) {
            *updatespec &= ~UPDATE_DOFIELDS;
            i++;
            continue;
        } else if (strcmp(argv[i], "--no-update") == 0) {
            no_update = 1;
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
            enkf_quit("parse_commandline(): option \"%s\" not recognised", argv[i]);
    }

    if (*fname == NULL)
        enkf_quit("parse_commandline(): parameter file not specified");

    if (no_update) {
        *updatespec &= ~UPDATE_NEEDAN;
        *updatespec &= ~UPDATE_DOPLOGSAN;
        *updatespec &= ~UPDATE_DOINFLATION;
    }
}

/**
 */
static void describe_updatespec(int updatespec)
{
    enkf_printf("  update specs:\n");
    enkf_printf("    do model fields   = %s\n", (updatespec & UPDATE_DOFIELDS) ? "[+]" : "[-]");
    enkf_printf("    do spread         = %s\n", (updatespec & UPDATE_DOSPREAD) ? "[+]" : "[-]");
    enkf_printf("    do pointlogs      = %s\n", (updatespec & UPDATE_DOPLOGS) ? "[+]" : "[-]");
    if (updatespec & UPDATE_DIRECTWRITE)
        enkf_printf("    direct write      = [+]\n");
    if (updatespec & UPDATE_DOFIELDS) {
        if (updatespec & UPDATE_OUTPUTINC)
            enkf_printf("    output increment  = [+]\n");
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

    parse_commandline(argc, argv, &fname_prm, &updatespec);

    enkf_init(&argc, &argv);

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
        if (prm->mode == MODE_ENOI) {
            enkf_printf("  prm: mode = EnOI -> cancelling writing of inflation\n");
            updatespec &= ~UPDATE_DOINFLATION;
        } else if (isnan(prm->inf_ratio)) {
            enkf_printf("  prm: inflation mode = plain -> cancelling writing of inflation\n");
            updatespec &= ~UPDATE_DOINFLATION;
        }
    }
    if (prm->nplog == 0)
        updatespec &= ~UPDATE_DOPLOGS;

    describe_updatespec(updatespec);
    if ((updatespec & (UPDATE_DOFIELDS | UPDATE_DOSPREAD | UPDATE_DOPLOGS | UPDATE_DOINFLATION)) == 0)
        enkf_quit("nothing to do");

    enkf_printf("  initialising the system:\n");
    das = das_create(prm, updatespec);
    enkfprm_destroy(prm);

#if defined(MPI)
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    if (das->updatespec & UPDATE_DOPLOGS && rank == 0) {
        enkf_printf("  defining state variables in point logs:");
        plogs_definestatevars(das);
        enkf_printf("\n");
        enkf_flush();
        dir_createifabsent(DIRNAME_TMP);
    }

    if (das->updatespec & UPDATE_DOPLOGSFC) {
        enkf_printf("  writing forecast variables to point logs:\n");
        enkf_printtime("    ");
        plogs_writestatevars(das, 0);
    }

    if (das->updatespec & UPDATE_DOSPREAD && rank == 0) {
        enkf_printf("  allocating disk space for spread:\n");
        enkf_printtime("    ");
        das_allocatedst(das, ROOTNAME_SPREAD);
        enkf_flush();
        if (rank == 0 && !(das->updatespec & UPDATE_DIRECTWRITE))
            dir_createifabsent(DIRNAME_TMP);
    }

    if (das->updatespec & (UPDATE_DOFIELDS | UPDATE_DOPLOGSAN | UPDATE_DOINFLATION)) {
        if (das->mode == MODE_ENKF) {
            if (updatespec & UPDATE_NEEDAN)
                enkf_printf("  updating ensemble:\n");
            else
                enkf_printf("  processing ensemble:\n");
        } else if (das->mode == MODE_ENOI) {
            if (updatespec & UPDATE_NEEDAN)
                enkf_printf("  updating model state:\n");
            else
                enkf_printf("  processing model state and/or ensemble:\n");
        }
        enkf_printtime("    ");
        das_update(das);
        enkf_flush();
    }

    das_destroy(das);

    enkf_finish();

    return 0;
}
