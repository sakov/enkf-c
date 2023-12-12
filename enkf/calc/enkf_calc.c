/******************************************************************************
 *
 * File:        enkf_calc.c        
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

#define EPS_IJ 1.0e-6

observation* singleob = NULL;
char* singleobtype = NULL;
int ignorenoobs = 0;
int use_rmsd = 0;
int plogs_only = 0;
int skip_transforms = 0;
int print_mem = 0;
int write_HE = 0;
int strict_time_matching = 0;
int skip_bad_fc_obs = 0;

/**
 */
static void usage(void)
{
    enkf_printf("  Usage: enkf_calc <prm file> [<options>]\n");
    enkf_printf("  Options:\n");
    enkf_printf("  --allow-logspace-with-static-ens\n");
    enkf_printf("      confirm that static ensemble is conditioned for using log space\n");
    enkf_printf("  --describe-prm-format [main|model|grid|obstypes]\n");
    enkf_printf("      describe format of a parameter file and exit\n");
    enkf_printf("  --forecast-stats-only\n");
    enkf_printf("      calculate and print forecast observation stats only\n");
    enkf_printf("  --ignore-no-obs\n");
    enkf_printf("      proceed even if there are no observations\n");
    enkf_printf("  --point-logs-only\n");
    enkf_printf("      skip calculating transforms for the whole grid and observation stats\n");
    enkf_printf("  --print-memory-usage\n");
    enkf_printf("      print memory usage by each process\n");
    enkf_printf("  --single-observation <lon> <lat> <depth> <type> <inn> <std>\n");
    enkf_printf("      assimilate single observation with these parameters\n");
    enkf_printf("  --skip-bad-forecast-obs\n");
    enkf_printf("      skip observations with invalid forecasts\n");
    enkf_printf("  --strict-time-matching\n");
    enkf_printf("      when assimilating asynchronously -- check that the time of model dumps\n");
    enkf_printf("      matches centres of the corresponding time bins\n");
    enkf_printf("  --use-existing-transforms\n");
    enkf_printf("      skip calculating ensemble transforms; use existing transforms*.nc files\n");
    enkf_printf("  --use-rmsd-for-obsstats\n");
    enkf_printf("      use RMSD instead of MAD when printing observation stats\n");
    enkf_printf("  --use-these-obs <obs file>\n");
    enkf_printf("      assimilate observations from this file; the file format must be compatible\n");
    enkf_printf("      with that of observations.nc produced by `enkf_prep'\n");
    enkf_printf("  --version\n");
    enkf_printf("      print version and exit\n");
    enkf_printf("  --write-HE\n");
    enkf_printf("      write ensemble observations to file \"%s\"\n", FNAME_HE);

    exit(0);
}

/**
 */
static void parse_commandline(int argc, char* argv[], char** fname_prm, char** fname_obs)
{
    int i;

    if (argc < 2)
        usage();

    i = 1;
    while (i < argc) {
        if (argv[i][0] != '-') {
            if (*fname_prm == NULL) {
                *fname_prm = argv[i];
                i++;
                continue;
            } else
                usage();
        } else if (strcmp(argv[i], "--allow-logspace-with-static-ens") == 0) {
            enkf_allowenoilog = 1;
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
                else
                    usage();
            } else
                enkfprm_describeprm();
            exit(0);
        } else if (strcmp(argv[i], "--ignore-no-obs") == 0) {
            ignorenoobs = 1;
            i++;
            continue;
        } else if (strcmp(argv[i], "--print-memory-usage") == 0) {
            print_mem = 1;
            i++;
            continue;
        } else if (strcmp(argv[i], "--point-logs-only") == 0) {
            plogs_only = 1;
            i++;
            continue;
        } else if (strcmp(argv[i], "--single-observation") == 0) {
            singleob = calloc(1, sizeof(observation));
            i++;
            if (i >= argc)
                usage();
            if (!str2float(argv[i], &singleob->lon))
                enkf_quit("command line: could not convert \"%s\" to float\n", argv[i]);
            i++;
            if (i >= argc)
                usage();
            if (!str2float(argv[i], &singleob->lat))
                enkf_quit("command line: could not convert \"%s\" to float\n", argv[i]);
            i++;
            if (i >= argc)
                usage();
            if (!str2float(argv[i], &singleob->depth))
                enkf_quit("command line: could not convert \"%s\" to float", argv[i]);
            i++;
            if (i >= argc)
                usage();
            singleobtype = argv[i];
            i++;
            if (i >= argc)
                usage();
            if (!str2float(argv[i], &singleob->value))
                enkf_quit("command line: could not convert \"%s\" to float", argv[i]);
            i++;
            if (i >= argc)
                usage();
            if (!str2float(argv[i], &singleob->estd))
                enkf_quit("command line: could not convert \"%s\" to float", argv[i]);
            i++;
            continue;
        } else if (strcmp(argv[i], "--skip-bad-forecast-obs") == 0) {
            skip_bad_fc_obs = 1;
            i++;
            continue;
        } else if (strcmp(argv[i], "--strict-time-matching") == 0) {
            strict_time_matching = 1;
            i++;
            continue;
        } else if (strcmp(argv[i], "--use-existing-transforms") == 0) {
            skip_transforms = 1;
            i++;
            continue;
        } else if (strcmp(argv[i], "--use-these-obs") == 0) {
            i++;
            if (i >= argc)
                usage();
            if (*fname_obs == NULL) {
                *fname_obs = strdup(argv[i]);
                enkf_noobsdatecheck = 1;
                i++;
                continue;
            } else
                usage();
        } else if (strcmp(argv[i], "--use-rmsd-for-obsstats") == 0) {
            use_rmsd = 1;
            i++;
            continue;
        } else if (strcmp(argv[i], "--forecast-stats-only") == 0) {
            enkf_fstatsonly = 1;
            i++;
            continue;
        } else if (strcmp(argv[i], "--version") == 0) {
            enkf_printversion();
            exit(0);
        } else if (strcmp(argv[i], "--write-HE") == 0) {
            write_HE = 1;
            i++;
            continue;
        } else
            enkf_quit("command line: option \"%s\" not recognised", argv[i]);
    }

    if (singleob != NULL && enkf_fstatsonly != 0)
        enkf_quit("command line: \"--forecast-stats-only\" is not compatible with \"--single-observation\"");
    if (*fname_obs != NULL && singleob != NULL)
        enkf_quit("command line: \"--use-these-obs\" is not compatible with \"--single-observation\"");
    if (plogs_only != 0 && enkf_fstatsonly != 0)
        enkf_quit("command line: \"--point-logs-only\" is not compatible with \"--forecast-stats-only\"");
    if (plogs_only != 0 && singleob != NULL)
        enkf_quit("command line: \"--point-logs-only\" is not compatible with \"--single-observation\"");

    if (*fname_prm == NULL)
        enkf_quit("command line: parameter file not specified");
}

/**
 */
static observations* obs_create_fromsingleob(enkfprm* prm, dasystem* das)
{
    model* m = das->m;
    observations* obs = obs_create();
    observation* o = singleob;
    int vid = -1;
    void* g = NULL;

    enkf_printf("  reading observation type specs from \"%s\":\n", prm->obstypeprm);
    obstypes_read(prm, prm->obstypeprm, &obs->nobstypes, &obs->obstypes);

    obs->da_time = date2day(prm->fname, prm->date);
    obs->datestr = strdup(prm->date);

    o->type = obstype_getid(obs->nobstypes, obs->obstypes, singleobtype, 1);

    vid = model_getvarid(m, obs->obstypes[o->type].varnames[0], 1);
    g = model_getvargrid(m, vid);

    obs->obstypes[o->type].gridid = model_getvargridid(das->m, vid);

    obs->products = st_create("products");
    st_add_ifabsent(obs->products, "Synthetic", -1);
    obs->instruments = st_create("instruments");
    st_add_ifabsent(obs->instruments, "Virtual", -1);

    obs->nobs = 1;
    obs->data = o;

    if (obs->obstypes[o->type].issurface)
        o->depth = 0.0;

    o->status = model_xy2fij(m, vid, o->lon, o->lat, o->fij);
    if (o->status == STATUS_OK)
        o->status = model_z2fk_f(m, vid, o->fij, o->depth, &o->fk);
    else
        o->fk = NAN;
    if (o->status != STATUS_OK)
        enkf_quit("command line: could not map the observation");

    enkf_printf("  assimilating single observation:\n");
    enkf_printf("    type = %s\n", singleobtype);
    enkf_printf("    inn  = %.3f\n", singleob->value);
    enkf_printf("    estd = %.3f\n", singleob->estd);
    enkf_printf("    lon  = %.3f\n", o->lon);
    enkf_printf("    lon  = %.3f\n", o->lat);
    if (grid_isstructured(g)) {
        enkf_printf("    i    = %.3f\n", o->fij[0]);
        enkf_printf("    j    = %.3f\n", o->fij[1]);
    } else {
        enkf_printf("    i0   = %.3f\n", o->fij[0]);
        enkf_printf("    i1   = %.3f\n", o->fij[1]);
        enkf_printf("    i2   = %.3f\n", o->fij[2]);
    }
    if (!obs->obstypes[o->type].issurface)
        enkf_printf("    k    = %.3f\n", o->fk);

    obs_calcstats(obs);
    return obs;
}

/**
 */
int main(int argc, char* argv[])
{
    char* fname_prm = NULL;
    char* fname_obs = NULL;
    enkfprm* prm = NULL;
    dasystem* das = NULL;

    enkf_init(&argc, &argv);

    parse_commandline(argc, argv, &fname_prm, &fname_obs);

    enkf_printf("  running CALC for EnKF-C version %s:\n", ENKF_VERSION);
    print_commandinfo(argc, argv);
    enkf_printtime("  ");

    enkf_printf("  reading system specs from \"%s\":\n", fname_prm);
    prm = enkfprm_read(fname_prm);
    enkfprm_print(prm, "    ");

    if (enkf_fstatsonly)
        enkf_doplogs = 0;

    enkf_printf("  initialising the system:\n");
    das = das_create(prm);
    das->strict_time_matching = strict_time_matching;

    if (print_mem)
        print_memory_usage();

    if (singleob == NULL) {
        if (fname_obs == NULL)
            fname_obs = strdup(FNAME_SOBS);

        enkf_printf("  reading observations from \"%s\":\n", fname_obs);
        obs_read(das->obs, fname_obs);
    } else {
        das->obs = obs_create_fromsingleob(prm, das);
        enkf_obstype = OBSTYPE_INNOVATION;
    }
    enkfprm_destroy(prm);

    if (print_mem)
        print_memory_usage();

    if (model_destroygxytrees(das->m))
        if (print_mem) {
            enkf_printf("  (destroyed grid kd-trees)\n");
            print_memory_usage();
        }

    enkf_printf("    creating kd-trees for observations:\n");
    obs_createkdtrees(das->obs);

    if (print_mem)
        print_memory_usage();

    if (das->obs->nobs == 0 && !ignorenoobs)
        enkf_quit("nothing to do! (nobs = 0). Use \"--ignore-no-obs\" to proceed cleanly");

    enkf_printf("  calculating ensemble observations:\n");
    enkf_printtime("  ");
    das_getHE(das);
    if (write_HE)
        das_writeHE(das);
    enkf_printf("  calculating innovations and spread:\n");
    enkf_printtime("  ");
    das_calcinnandspread(das);

    if (print_mem)
        print_memory_usage();

    if (singleob == NULL) {
        enkf_printf("  writing forecast innovations and spread to \"%s\":\n", fname_obs);
        enkf_printtime("  ");
        das_writeforecastobs(das, fname_obs);
    }

    /*
     * calculate stats for observation batches and write to FNAME_BATCHES
     */
#if defined(USE_SHMEM)
    if (singleob == NULL && das->obs->nobs > 0 && sm_comm_rank == 0) {
#else
    if (singleob == NULL && das->obs->nobs > 0) {
#endif
        hashtable* batches = das_getbatches(das);
        hashtable* badbatches = das_processbatches(das, batches);

        obs_markbadbatches(das->obs, badbatches);

        ht_destroy(batches);
        ht_destroy(badbatches);

        /*
         * update observation status in observation.nc
         */
        if (rank == 0)
            obs_writeobsstatus(das->obs, fname_obs);
    }

    if (!enkf_fstatsonly) {
#if defined(MPI)
        MPI_Barrier(MPI_COMM_WORLD);
#endif
        if (singleob == NULL) {
            int done;

            enkf_printf("  moderating observations:\n");
            done = obs_modifiederrors_alreadywritten(das->obs, fname_obs);
#if defined(MPI)
            MPI_Barrier(MPI_COMM_WORLD);
#endif
            if (done)
                enkf_printf("    already done\n");
            else {
                das_moderateobs(das);
                enkf_printf("  writing modified obs errors to \"%s\":\n", fname_obs);
                if (rank == 0)
                    das_writemoderatedobs(das, fname_obs);
            }
        }

        if (!plogs_only && !skip_transforms) {
            enkf_printf("  calculating transforms:\n");
            enkf_printtime("  ");
            das_calctransforms(das);
        } else
            das_standardise(das);

        if (das->nplog > 0 && das->obs->loctrees == NULL)
            obs_createkdtrees(das->obs);
        if (rank == 0) {
            enkf_printf("  writing point logs:\n");
            das_calcpointlogtransforms(das);
        }

        if (print_mem)
            print_memory_usage();

        obs_destroykdtrees(das->obs);

        if (!plogs_only) {
            /*
             * the following is an optional bit - updating ensemble
             * observations and generating report 
             */
            enkf_printf("  calculating analysed observations:\n");
            enkf_printtime("  ");
            das_updateHE(das);

            if (singleob == NULL) {
                enkf_printf("  adding analysis innovations and spread to \"%s\":\n", fname_obs);
                das_writeanalysisobs(das, fname_obs);
            }

            enkf_printf("  printing observation statistics:\n");
            das_printobsstats(das, use_rmsd);
        }
    } else {
        enkf_printf("  printing observation statistics:\n");
        das_printfobsstats(das, use_rmsd);
    }

    das_destroy(das);
    free(fname_obs);
    enkf_finish();

    return 0;
}
