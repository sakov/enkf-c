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
#include <assert.h>
#include <math.h>
#include "version.h"
#include "definitions.h"
#include "utils.h"
#include "enkfprm.h"
#include "obsprm.h"
#include "grid.h"
#include "model.h"
#include "observations.h"
#include "prep_utils.h"

/*
 * describe this superobservation and exit (-1 = normal workflow)
 */
int describe_superob_id = -1;

/*
 * superobing can be switched off
 */
int do_superob = 1;

/*
 * superobing across instruments can be switched on
 */
int do_superob_acrossinst = 0;

/*
 * superobing across batches can be switched off
 */
int do_superob_acrossbatches = 1;

/*
 * writing of the original obs can be swithched on
 * 0 = do not write
 * 1 = write original obs within model domain
 * 2 = write all original obs
 */
int write_orig_obs = 0;

/*
 * thinning of obs with identical positions in the same time window can be
 * switched off
 */
int do_thin = 1;

/**
 */
static void usage()
{
    enkf_printf("  Usage: enkf_prep <prm file> [<options>]\n");
    enkf_printf("  Options:\n");
    enkf_printf("  --allow-logspace-with-static-ens\n");
    enkf_printf("      confirm that static ensemble is conditioned for using log space\n");
    enkf_printf("  --consider-subgrid-variability\n");
    enkf_printf("      increase error of superobservations according to subgrid variability\n");
    enkf_printf("  --describe-prm-format [main|model|grid|obstypes|obsdata]\n");
    enkf_printf("      describe format of a parameter file and exit\n");
    enkf_printf("  --describe-reader <reader>\n");
    enkf_printf("      decribe reader and exit\n");
    enkf_printf("  --describe-superob <sob #>\n");
    enkf_printf("      print composition of this superobservation and exit\n");
    enkf_printf("  --list-readers\n");
    enkf_printf("      list available readers and exit\n");
    enkf_printf("  --log-all-obs\n");
    enkf_printf("      write all obs to %s (default: obs within model domain only)\n", FNAME_OBS);
    enkf_printf("  --no-superobing\n");
    enkf_printf("  --no-superobing-across-batches\n");
    enkf_printf("  --no-thinning\n");
    enkf_printf("  --superob-across-instruments\n");
    enkf_printf("  --write-orig-obs\n");
    enkf_printf("      write original obs within model domain to %s\n", FNAME_OBS);
    enkf_printf("  --write-all-orig-obs\n");
    enkf_printf("      write all original obs to %s\n", FNAME_OBS);
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
        } else if (strcmp(argv[i], "--allow-logspace-with-static-ens") == 0) {
            enkf_allowenoilog = 1;
            i++;
            continue;
        } else if (strcmp(argv[i], "--consider-subgrid-variability") == 0) {
            enkf_considersubgridvar = 1;
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
                else if (strcmp(argv[i + 1], "obstypes") == 0)
                    obstypes_describeprm();
                else if (strcmp(argv[i + 1], "obsdata") == 0)
                    obsprm_describeprm();
                else
                    usage();
            } else
                enkfprm_describeprm();
            exit(0);
        } else if (strcmp(argv[i], "--describe-reader") == 0) {
            i++;
            if (i >= argc)
                enkf_quit("usage: reader name not found");
            describe_reader(argv[i]);
            exit(0);
        } else if (strcmp(argv[i], "--describe-superob") == 0) {
            i++;
            if (i >= argc)
                enkf_quit("usage: superob ID not found");
            if (!str2int(argv[i], &describe_superob_id))
                enkf_quit("usage: could not convert \"%s\" to integer", argv[i]);
            i++;
            continue;
        } else if (strcmp(argv[i], "--list-readers") == 0) {
            list_readers();
            exit(0);
        } else if (strcmp(argv[i], "--no-superobing") == 0) {
            do_superob = 0;
            i++;
            continue;
        } else if (strcmp(argv[i], "--no-superobing-across-batches") == 0) {
            do_superob_acrossbatches = 0;
            i++;
            continue;
        } else if (strcmp(argv[i], "--no-thinning") == 0) {
            do_thin = 0;
            i++;
            continue;
        } else if (strcmp(argv[i], "--superob-across-instruments") == 0) {
            do_superob_acrossinst = 1;
            i++;
            continue;
        } else if (strcmp(argv[i], "--write-orig-obs") == 0) {
            write_orig_obs = 1;
            i++;
            continue;
        } else if (strcmp(argv[i], "--write-all-orig-obs") == 0) {
            write_orig_obs = 2;
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
    observation* o1 = (observation*) p1;
    observation* o2 = (observation*) p2;
    observations* obs = (observations*) p;
    obstype* ot;
    int stride;
    double offset;
    int i1, i2;

    if (o1->type > o2->type)
        return 1;
    if (o1->type < o2->type)
        return -1;

    if (!do_superob_acrossinst) {
        if (o1->instrument > o2->instrument)
            return 1;
        if (o1->instrument < o2->instrument)
            return -1;
    }

    if (!do_superob_acrossbatches) {
        if (o1->batch > o2->batch)
            return 1;
        if (o1->batch < o2->batch)
            return -1;
    }

    ot = &obs->obstypes[o1->type];
    if (ot->isasync) {
        i1 = get_tshift(o1->time, ot->async_tstep, ot->async_centred);
        i2 = get_tshift(o2->time, ot->async_tstep, ot->async_centred);
        if (i1 > i2)
            return 1;
        if (i2 > i1)
            return -1;
    }

    if (o1->footprint > o2->footprint)
        return 1;
    else if (o1->footprint < o2->footprint)
        return -1;

    stride = ot->sob_stride;
    if (stride == 0) {          /* no superobing on this grid */
        if (o1->id > o2->id)
            return 1;
        else if (o1->id < o2->id)
            return -1;
        return 0;
    }

    /*
     * The offset below is supposed to align the superobservation cell
     * boundaries with the cell boundaries.
     */
    offset = (double) (stride % 2) / 2.0;

    i1 = (int) floor(o1->fi + offset) / stride;
    i2 = (int) floor(o2->fi + offset) / stride;
    if (i1 > i2)
        return 1;
    if (i1 < i2)
        return -1;

    i1 = (int) floor(o1->fj + offset) / stride;
    i2 = (int) floor(o2->fj + offset) / stride;
    if (i1 > i2)
        return 1;
    if (i1 < i2)
        return -1;

    i1 = (int) floor(o1->fk + 0.5);
    i2 = (int) floor(o2->fk + 0.5);
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
    model* m;
    int nmeta = 0;
    obsmeta* meta = NULL;
    int nexclude = 0;
    obsregion* exclude = NULL;
    observations* obs = NULL;
    observations* sobs = NULL;
    int i;

    enkf_init(&argc, &argv);

    parse_commandline(argc, argv, &fname_prm);

    enkf_printf("  running PREP for EnKF-C version %s:\n", ENKF_VERSION);
    print_commandinfo(argc, argv);
    enkf_printtime("  ");

    enkf_printf("  reading prep specs from \"%s\":\n", fname_prm);
    prm = enkfprm_read(fname_prm);
    enkfprm_print(prm, "    ");

    enkf_printf("  reading observation specs from \"%s\":\n", prm->obsprm);
    obsprm_read(prm->obsprm, &nmeta, &meta, &nexclude, &exclude);

    enkf_printf("  creating model:\n");
    m = model_create(prm);
    enkf_printf("  creating observations:\n");
    obs = obs_create_fromprm(prm);
    obstypes_set(obs->nobstypes, obs->obstypes, m);
    obs->allobs = (write_orig_obs == 2) ? 1 : 0;
    obs->model = m;

    for (i = 0; i < nexclude; ++i) {
        if (strcasecmp(exclude[i].otname, "ALL") == 0)
            exclude[i].otid = -1;
        else
            exclude[i].otid = obstype_getid(obs->nobstypes, obs->obstypes, exclude[i].otname, 1);
    }

    enkf_printf("  reading observations:\n");
    for (i = 0; i < nmeta; i++)
        obs_add(obs, m, &meta[i], nexclude, exclude);
    obs_markbadbatches(obs);
    obsprm_destroy(nmeta, meta, nexclude, exclude);
    enkf_printtime("  ");
    enkf_printf("  compacting obs:\n");
    obs_compact(obs);
    obs_calcstats(obs);

    if (do_superob) {
        enkf_printtime("  ");
        enkf_printf("  sorting:\n");
        enkf_flush();
        qsort_r(obs->data, obs->ngood, sizeof(observation), cmp_obs, obs);
        enkf_printtime("  ");
        enkf_printf("  superobing:\n");
        enkf_flush();
        obs_superob(obs, cmp_obs, &sobs, describe_superob_id, do_thin);

        if (describe_superob_id >= 0)
            goto finalise;

        enkf_printtime("  ");
        enkf_printf("  checking for superobs on land:\n");
        enkf_flush();
        if (obs_checkforland(sobs, m)) {
            obs_compact(sobs);
            for (i = 0; i < sobs->nobs; ++i)
                if (sobs->data[i].status != STATUS_OK)
                    break;
            assert(i != sobs->nobs);
            enkf_printf("    deleted %d observation(s)\n", sobs->nobs - i);
            sobs->nobs = i;
            obs_calcstats(sobs);
        } else
            enkf_printf("    all good\n");

        enkf_printtime("  ");
        enkf_printf("  writing superobservations to \"%s\":\n", FNAME_SOBS);
        obs_write(sobs, FNAME_SOBS);
    } else {
        observation* data = malloc(obs->ngood * sizeof(observation));

        memcpy(data, obs->data, obs->ngood * sizeof(observation));
        sobs = obs_create_fromdata(obs, obs->ngood, data);
        obs_write(sobs, FNAME_SOBS);
        goto finalise;
    }

    if (write_orig_obs && describe_superob_id < 0) {
        enkf_printf("  writing observations to \"%s\":\n", FNAME_OBS);
        obs_write(obs, FNAME_OBS);
    }

    /*
     * write superob indices to the file with original observations 
     */
    if (write_orig_obs && describe_superob_id < 0 && do_superob)
        obs_writeaux(obs, FNAME_OBS);

    if (enkf_considersubgridvar) {
        int firsttime = 1;

        for (i = 0; i < sobs->nobstypes; ++i) {
            obstype* ot = &sobs->obstypes[i];

            if (ot->nsubgrid > 0) {
                if (firsttime) {
                    enkf_printf("  # obs with increased error due to subgrid variability:\n");
                    firsttime = 0;
                }
                enkf_printf("    %s %d (%.2f%%)\n", ot->name, ot->nsubgrid, (double) ot->nsubgrid / (double) ot->ngood * 100.0);
            }
        }
    }

  finalise:

    if (describe_superob_id < 0) {
        enkf_printf("  printing observation summary:\n");
        print_obsstats(obs, sobs);
    }

    obs_destroy(obs);
    obs_destroy(sobs);
    model_destroy(m);
    enkfprm_destroy(prm);

    enkf_finish();

    return 0;
}
