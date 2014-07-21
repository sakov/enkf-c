/******************************************************************************
 *
 * File:        enkf_prep.c        
 *
 * Created:     12/2012
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description:
 *
 * Revisions:
  *****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "nan.h"
#include "version.h"
#include "definitions.h"
#include "utils.h"
#include "enkfprm.h"
#include "obsmeta.h"
#include "grid.h"
#include "model.h"
#include "observations.h"
#include "prep.h"

/*
 * put all obs int observations-orig.nc (default = only local obs)
 */
int log_all_obs = 0;

/*
 * describe this superobservation and exit
 */
int describe_superob_id = -1;

/*
 * to communicate to sort()
 */
obstype* obstypes_p = NULL;

/**
 */
static void usage()
{
    enkf_printf("  Usage: enkf_prep <prm file> [<options>]\n");
    enkf_printf("  Options:\n");
    enkf_printf("  --log-all-obs\n");
    enkf_printf("      put all obs into observations-orig.nc (default: local obs only)\n");
    enkf_printf("  --describe-prm-format\n");
    enkf_printf("      describe format of the parameter file and exit\n");
    enkf_printf("  --describe-superob <sob #>\n");
    enkf_printf("      print composition of this superobservation and exit\n");
    enkf_printf("  --version\n");
    enkf_printf("      print version and exit\n");

    exit(0);
}

/**
 */
static void parse_commandline(int argc, char* argv[], char** fname)
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
        } else if (strcmp(argv[i], "--describe-prm-format") == 0) {
            enkfprm_describe();
            exit(0);
        } else if (strcmp(argv[i], "--describe-superob") == 0) {
            i++;
            if (i >= argc)
                usage();
            if (!str2int(argv[i], &describe_superob_id))
                enkf_quit("usage: could not convert \"%s\" to integer", argv[i]);
            i++;
            continue;
        } else if (strcmp(argv[i], "--log-all-obs") == 0) {
            log_all_obs = 1;
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
static int cmp_obs(const void* p1, const void* p2, void* p)
{
    observation* m1 = (observation*) p1;
    observation* m2 = (observation*) p2;
    int stride = ((observations*) p)->stride;
    double offset;
    int i1, i2;

    if (m1->type > m2->type)
        return 1;
    if (m1->type < m2->type)
        return -1;

    if (obstypes_p[m1->type].isasync) {
        i1 = get_tshift(m1->date, obstypes_p[m1->type].async_tstep);
        i2 = get_tshift(m2->date, obstypes_p[m2->type].async_tstep);
        if (i1 > i2)
            return 1;
        if (i2 > i1)
            return -1;
    }

    /*
     * the idea behind offset is to center superobservation domains in the case
     * of even stride on grid cells rather than on cell corners
     */
    offset = (double) ((stride - 1) % 2) / 2.0;

    i1 = (int) floor(m1->fi + offset) / stride;
    i2 = (int) floor(m2->fi + offset) / stride;
    if (i1 > i2)
        return 1;
    if (i1 < i2)
        return -1;

    i1 = (int) floor(m1->fj + offset) / stride;
    i2 = (int) floor(m2->fj + offset) / stride;
    if (i1 > i2)
        return 1;
    if (i1 < i2)
        return -1;

    i1 = (int) floor(m1->fk);
    i2 = (int) floor(m2->fk);
    if (i1 > i2)
        return 1;
    if (i1 < i2)
        return -1;

    return 0;
}

/**
 */
int main(int argc, char* argv[])
{
    char* fname_prm = NULL;
    enkfprm* prm = NULL;
    int nmeta = 0;
    obsmeta* meta = NULL;
    observations* obs = NULL;
    observations* sobs = NULL;
    model* m = NULL;
    int i;

    parse_commandline(argc, argv, &fname_prm);

    enkf_init(&argc, &argv);
    enkf_printf("  running PREP for EnKF version %s:\n", ENKF_VERSION);
    enkf_printtime("  ");

    enkf_printf("  reading prep specs from \"%s\":\n", fname_prm);
    prm = enkfprm_read(fname_prm);
    enkfprm_print(prm, "    ");

    enkf_printf("  reading observation specs from \"%s\":\n", prm->obsspec);
    read_obsmeta(prm, &nmeta, &meta);

    enkf_printf("  setting the model grid:\n");
    m = model_create(prm);

    obs = obs_create_fromprm(prm);
    enkfprm_destroy(prm);
    obs->allobs = log_all_obs;

    enkf_printf("  reading observations:\n");
    for (i = 0; i < nmeta; i++)
        obs_add(obs, m, &meta[i]);
    obs_markbadbatches(obs);
    clean_obsmeta(nmeta, meta);
    obs_checklon(obs);
    obs_compact(obs);
    obs_calcstats(obs);

    if (describe_superob_id < 0) {
        enkf_printf("  writing observations to \"%s\":\n", FNAME_OBS);
        obs_write(obs, FNAME_OBS);
    }

    enkf_printf("  superobing:\n");
    obstypes_p = obs->obstypes;
    obs_superob(obs, cmp_obs, &sobs, describe_superob_id);
    /*
     * write superob indices to the file with original observations 
     */
    obs_writeaux(obs, FNAME_OBS);

    if (describe_superob_id < 0) {
        enkf_printf("  writing superobservations to \"%s\":\n", FNAME_SOBS);
        obs_write(sobs, FNAME_SOBS);
    }

    enkf_printf("  printing observation statistics:\n");
    print_obsstats(obs, sobs);

    obs_destroy(obs);
    obs_destroy(sobs);
    model_destroy(m);

    enkf_finish();

    return 0;
}
