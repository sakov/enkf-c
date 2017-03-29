/******************************************************************************
 *
 * File:        enkfprm.c        
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
#include <float.h>
#include <math.h>
#include "ncw.h"
#include "definitions.h"
#include "utils.h"
#include "enkfprm.h"

#define NINC 10

/**
 */
static void enkfprm_check(enkfprm* prm)
{
    if (prm->mode == MODE_NONE)
        enkf_quit("%s: MODE not specified", prm->fname);
#if defined(ENKF_PREP)
    if (prm->obsprm == NULL)
        enkf_quit("%s: OBS not specified", prm->fname);
    if (prm->date == NULL)
        enkf_quit("%s: DATE not specified", prm->fname);
#endif
    if (prm->modelprm == NULL)
        enkf_quit("%s: MODEL not specified", prm->fname);
    if (prm->gridprm == NULL)
        enkf_quit("%s: GRID not specified", prm->fname);
#if defined(ENKF_PREP) || defined(ENKF_CALC)
    if (prm->obstypeprm == NULL)
        enkf_quit("%s: OBSTYPES not specified", prm->fname);
#endif
#if defined(ENKF_CALC) || defined(ENKF_UPDATE)
    if (prm->ensdir == NULL && (prm->mode == MODE_ENKF || !enkf_fstatsonly))
        enkf_quit("%s: ENSDIR not specified", prm->fname);
#endif
#if defined(ENKF_CALC)
    if (prm->mode == MODE_ENOI && prm->bgdir == NULL)
        enkf_quit("%s: BGDIR must be specified for MODE = ENOI", prm->fname);
    /*
     * (we skip the test for ENKF_UPDATE because (1) there are cases when BGDIR
     * is not requiredt, and (2) normally the same parameter file is used for
     * ENKF_CALC and ENKF_UPDATE)
     */
#endif
#if defined(ENKF_CALC)
    if (!isfinite(prm->locrad) && !enkf_fstatsonly)
        enkf_quit("%s: LOCRAD not specified", prm->fname);
#endif
    if (prm->ncformat != NC_CLASSIC_MODEL && prm->ncformat != NC_64BIT_OFFSET && prm->ncformat != NC_NETCDF4)
        enkf_quit("programming error");
    if (prm->nccompression > 0 && prm->ncformat != NC_NETCDF4)
        enkf_quit("%s: NCFORMAT must be set to NETCDF4 when NCCOMPRESSION > 0\n", prm->fname);
    if (prm->nccompression < 0 || prm->nccompression > 9)
        enkf_quit("%s: NCCOMPRESSION must be in the interval between 0 and 9\n", prm->fname);
}

/**
 */
enkfprm* enkfprm_read(char fname[])
{
    enkfprm* prm = NULL;
    FILE* f = NULL;
    char buf[MAXSTRLEN];
    int line;
    int i;

    prm = calloc(1, sizeof(enkfprm));   /* everything is initialised to 0 */
    prm->fname = fname;
    prm->mode = MODE_NONE;
    prm->scheme = SCHEME_DEFAULT;
    prm->alpha = ALPHA_DEFAULT;
    prm->date = NULL;
    prm->windowmin = NAN;
    prm->windowmax = NAN;
    prm->modelprm = NULL;
    prm->gridprm = NULL;

    prm->obstypeprm = NULL;
    prm->obsprm = NULL;
    prm->enssize = -1;
    prm->kfactor = NAN;
    prm->rfactor_base = 1.0;
    prm->inflation_base = 1.0;
    prm->inf_ratio = INFRATIO_DEFAULT;
    prm->locrad = NAN;
    prm->stride = 1;
    prm->fieldbufsize = 1;
    prm->sob_stride = 1;
    prm->ncformat = NC_64BIT_OFFSET;
    prm->nccompression = 0;

    f = enkf_fopen(fname, "r");
    line = 0;
    while (fgets(buf, MAXSTRLEN, f) != NULL) {
        char seps[] = " =\t\n";
        char* token;

        line++;
        if (buf[0] == '#')
            continue;
        if ((token = strtok(buf, seps)) == NULL)
            continue;
        if (strcasecmp(token, "MODE") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: MODE not specified", fname, line);
            else {
                if (strcasecmp(token, "ENKF") == 0)
                    prm->mode = MODE_ENKF;
                else if (strcasecmp(token, "ENOI") == 0)
                    prm->mode = MODE_ENOI;
                else
                    enkf_quit("%s, l.%d: mode \"%s\" is not known", fname, line, token);
            }
        } else if (strcasecmp(token, "SCHEME") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: SCHEME not specified", fname, line);
            else {
                if (strcasecmp(token, "DENKF") == 0)
                    prm->scheme = SCHEME_DENKF;
                else if (strcasecmp(token, "ETKF") == 0)
                    prm->scheme = SCHEME_ETKF;
                else
                    enkf_quit("%s, l.%d: scheme \"%s\" is not implemented", fname, line, token);
            }
        } else if (strcasecmp(token, "ALPHA") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: ALPHA not specified", fname, line);
            else if (!str2double(token, &prm->alpha))
                enkf_quit("%s, l.%d: could not convert \"%s\" to double", fname, line, token);
        } else if (strcasecmp(token, "DATE") == 0) {
            char seps_date[] = "=\n";

            if ((token = strtok(NULL, seps_date)) == NULL)
                enkf_quit("%s, l.%d: DATE not specified", fname, line);
            else if (prm->date != NULL)
                enkf_quit("%s, l.%d: DATE specified twice", fname, line);
            else {
                while (*token == ' ')
                    token++;
                prm->date = strdup(token);
            }
        } else if (strcasecmp(token, "WINDOWMIN") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: WINDOWMIN not specified", fname, line);
            else if (!str2double(token, &prm->windowmin))
                enkf_quit("%s, l.%d: could convert WINDOWMIN entry", fname, line);
        } else if (strcasecmp(token, "WINDOWMAX") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: WINDOWMAX not specified", fname, line);
            else if (!str2double(token, &prm->windowmax))
                enkf_quit("%s, l.%d: could convert WINDOWMAX entry", fname, line);
        } else if (strcasecmp(token, "MODEL") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: MODEL not specified", fname, line);
            else if (prm->modelprm != NULL)
                enkf_quit("%s, l.%d: MODEL specified twice", fname, line);
            else
                prm->modelprm = strdup(token);
        } else if (strcasecmp(token, "GRID") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: GRID not specified", fname, line);
            else if (prm->gridprm !=NULL)
                enkf_quit("%s, l.%d: GRID specified twice", fname, line);
            else
                prm->gridprm = strdup(token);
        } else if (strcasecmp(token, "OBS") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: OBS not specified", fname, line);
            else if (prm->obsprm != NULL)
                enkf_quit("%s, l.%d: OBS specified twice", fname, line);
            else
                prm->obsprm = strdup(token);
        } else if (strcasecmp(token, "OBSTYPES") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: OBSTYPES not specified", fname, line);
            else if (prm->obstypeprm != NULL)
                enkf_quit("%s, l.%d: OBSTYPES specified twice", fname, line);
            else
                prm->obstypeprm = strdup(token);
        } else if (strcasecmp(token, "ENSDIR") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: ENSDIR not specified", fname, line);
            else if (prm->ensdir != NULL)
                enkf_quit("%s, l.%d: ENSDIR specified twice", fname, line);
            else
                prm->ensdir = strdup(token);
        } else if (strcasecmp(token, "ENSSIZE") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: ENSSIZE not specified", fname, line);
            else if (prm->enssize >= 0)
                enkf_quit("%s, l.%d: ENSSIZE specified twice", fname, line);
            else if (!str2int(token, &prm->enssize))
                enkf_quit("%s, l.%d: could convert ENSSIZE entry", fname, line);
        } else if (strcasecmp(token, "BGDIR") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: BGDIR not specified", fname, line);
            else if (prm->bgdir != NULL)
                enkf_quit("%s, l.%d: BGDIR specified twice", fname, line);
            else
                prm->bgdir = strdup(token);
        } else if (strcasecmp(token, "KFACTOR") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: KFACTOR not specified", fname, line);
            else if (isfinite(prm->kfactor))
                enkf_quit("%s, l.%d: KFACTOR specified twice", fname, line);
            else if (!str2double(token, &prm->kfactor))
                enkf_quit("%s, l.%d: could convert KFACTOR entry", fname, line);
        } else if (strcasecmp(token, "RFACTOR") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: observation type for RFACTOR not specified", fname, line);
            if (!str2double(token, &prm->rfactor_base))
                enkf_quit("%s, l.%d: could not convert RFACTOR \"%s\" to double", fname, line, token);
        } else if (strcasecmp(token, "LOCRAD") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: LOCRAD not specified", fname, line);
            if (isfinite(prm->locrad))
                enkf_quit("%s, l.%d: LOCRAD specified twice", fname, line);
            if (!str2double(token, &prm->locrad))
                enkf_quit("%s, l.%d: could not convert LOCRAD value", fname, line);
        } else if (strcasecmp(token, "STRIDE") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: STRIDE not specified", fname, line);
            if (prm->stride != 1)
                enkf_quit("%s, l.%d: STRIDE specified twice", fname, line);
            if (!str2int(token, &prm->stride))
                enkf_quit("%s, l.%d: could not convert STRIDE value", fname, line);
        } else if (strcasecmp(token, "SOBSTRIDE") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: SOBSTRIDE not specified", fname, line);
            if (prm->sob_stride != 1)
                enkf_quit("%s, l.%d: SOBSTRIDE specified twice", fname, line);
            if (!str2int(token, &prm->sob_stride))
                enkf_quit("%s, l.%d: could not convert SOBSTRIDE value", fname, line);
        } else if (strcasecmp(token, "FIELDBUFFERSIZE") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: FIELDBUFFERSIZE not specified", fname, line);
            if (prm->fieldbufsize != 1)
                enkf_quit("%s, l.%d: FIELDBUFFERSIZE specified twice", fname, line);
            if (!str2int(token, &prm->fieldbufsize))
                enkf_quit("%s, l.%d: could not convert FIELDBUFFERSIZE value", fname, line);
        } else if (strcasecmp(token, "INFLATION") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: INFLATION not specified", fname, line);
            if (!str2double(token, &prm->inflation_base))
                enkf_quit("%s, l.%d: could not convert \"%s\" to double", fname, line, token);
            if ((token = strtok(NULL, seps)) != NULL) {
                if (strcasecmp(token, "PLAIN") == 0)
                    prm->inf_ratio = NAN;
                else if (!str2double(token, &prm->inf_ratio))
                    enkf_quit("%s, l.%d: could not convert \"%s\" to double", fname, line, token);
            }
        } else if (strcasecmp(token, "REGION") == 0) {
            char zseps[] = " =\t\n[](){}";
            char* space;
            char* newtoken;
            region* r = NULL;

            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: REGION name not specified", fname, line);
            if (prm->nregions % NINC == 0)
                prm->regions = realloc(prm->regions, (prm->nregions + NINC) * sizeof(region));

            r = &prm->regions[prm->nregions];
            r->nzints = 0;
            r->zints = NULL;

            newtoken = token;
            while ((space = strchr(newtoken, '_')) != NULL) {
                *space = ' ';
                newtoken = space + 1;
            }
            r->name = strdup(token);
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: minimal longitude not specified", fname, line);
            if (!str2double(token, &r->x1))
                enkf_quit("%s, l.%d: could not convert \"%s\" to double", fname, line, token);
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: maximal longitude not specified", fname, line);
            if (!str2double(token, &r->x2))
                enkf_quit("%s, l.%d: could not convert \"%s\" to double", fname, line, token);
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: minimal latitude not specified", fname, line);
            if (!str2double(token, &r->y1))
                enkf_quit("%s, l.%d: could not convert \"%s\" to double", fname, line, token);
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: maximal latitude not specified", fname, line);
            if (!str2double(token, &r->y2))
                enkf_quit("%s, l.%d: could not convert \"%s\" to double", fname, line, token);
            while ((token = strtok(NULL, zseps)) != NULL) {
                if (r->nzints % NINC == 0)
                    r->zints = realloc(r->zints, (r->nzints + NINC) * sizeof(zint));
                if (!str2double(token, &r->zints[r->nzints].z1))
                    enkf_quit("%s, l.%d: could not convert \"%s\" to double", fname, line, token);
                if ((token = strtok(NULL, zseps)) == NULL)
                    enkf_quit("%s, l.%d: maximal depth/height for an interval not specified", fname, line);
                if (!str2double(token, &r->zints[r->nzints].z2))
                    enkf_quit("%s, l.%d: could not convert \"%s\" to double", fname, line, token);
                r->nzints++;
            }
            prm->nregions++;
        } else if (strcasecmp(token, "ZSTATINTS") == 0) {
            char zseps[] = " =\t\n[](){}";

            while ((token = strtok(NULL, zseps)) != NULL) {
                if (prm->nzints % NINC == 0)
                    prm->zints = realloc(prm->zints, (prm->nzints + NINC) * sizeof(zint));
                if (!str2double(token, &prm->zints[prm->nzints].z1))
                    enkf_quit("%s, l.%d: could not convert \"%s\" to double", fname, line, token);
                if ((token = strtok(NULL, zseps)) == NULL)
                    enkf_quit("%s, l.%d: maximal depth/height for an interval not specified", fname, line);
                if (!str2double(token, &prm->zints[prm->nzints].z2))
                    enkf_quit("%s, l.%d: could not convert \"%s\" to double", fname, line, token);
                prm->nzints++;
            }
        } else if (strcasecmp(token, "POINTLOG") == 0) {
            pointlog* plog = NULL;

            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: \"I\" coordinate for POINTLOG not specified", fname, line);
            if (prm->nplogs % NINC == 0)
                prm->plogs = realloc(prm->plogs, (prm->nplogs + NINC) * sizeof(pointlog));
            plog = &prm->plogs[prm->nplogs];
            plog->id = prm->nplogs;
            if (!str2int(token, &plog->i))
                enkf_quit("%s, l.%d: could convert \"I\" entry", fname, line);
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: \"J\" coordinate for POINTLOG not specified", fname, line);
            if (!str2int(token, &plog->j))
                enkf_quit("%s, l.%d: could convert \"J\" entry", fname, line);
            if ((token = strtok(NULL, seps)) != NULL && token[0] != '#')
                plog->gridname = strdup(token);
            else
                plog->gridname = NULL;
            plog->lon = NAN;
            plog->lat = NAN;
            prm->nplogs++;
        } else if (strcasecmp(token, "EXITACTION") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: EXITACTION not specified", fname, line);
            if (strcasecmp(token, "BACKTRACE") == 0)
                enkf_exitaction = EXITACTION_BACKTRACE;
            else if (strcasecmp(token, "SEGFAULT") == 0)
                enkf_exitaction = EXITACTION_SEGFAULT;
            else
                enkf_quit("%s, l.%d: EXITACTION \"%s\" not known", fname, line, token);
        } else if (strcasecmp(token, "BADBATCHES") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: BADBATCHES not specified", fname, line);
            if (prm->nbadbatchspecs % NINC == 0)
                prm->badbatchspecs = realloc(prm->badbatchspecs, (prm->nbadbatchspecs + NINC) * sizeof(badbatchspec));
            prm->badbatchspecs[prm->nbadbatchspecs].obstype = strdup(token);

            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: maximal allowed bias magnitude not defined", fname, line);
            if (!str2double(token, &prm->badbatchspecs[prm->nbadbatchspecs].maxbias))
                enkf_quit("%s, l.%d: could not convert \"%s\" to double", fname, line, token);
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: maximal allowed MAD not defined", fname, line);
            if (!str2double(token, &prm->badbatchspecs[prm->nbadbatchspecs].maxmad))
                enkf_quit("%s, l.%d: could not convert \"%s\" to double", fname, line, token);
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: minimal # obs. not defined", fname, line);
            if (!str2int(token, &prm->badbatchspecs[prm->nbadbatchspecs].minnobs))
                enkf_quit("%s, l.%d: could not convert \"%s\" to int", fname, line, token);
            if ((token = strtok(NULL, seps)) != NULL) {
                if (!str2int(token, &prm->badbatchspecs[prm->nbadbatchspecs].minnobs))
                    enkf_quit("%s, l.%d: could convert \"%s\" to integer", fname, line, token);
            }
            prm->nbadbatchspecs++;
        } else if (strcasecmp(token, "NCFORMAT") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: NCFORMAT not specified", fname, line);
            if (strcasecmp(token, "CLASSIC") == 0)
                prm->ncformat = NC_CLASSIC_MODEL;
            else if (strcasecmp(token, "64BIT") == 0)
                prm->ncformat = NC_64BIT_OFFSET;
            else if (strcasecmp(token, "NETCDF4") == 0)
                prm->ncformat = NC_NETCDF4;
            else
                enkf_quit("%s, l.%d: could not recognise specified NetCDF format\"%s\"", fname, line, token);
        } else if (strcasecmp(token, "NCCOMPRESSION") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: NCCOMPRESSION not specified", fname, line);
            if (!str2int(token, &prm->nccompression))
                enkf_quit("%s, l.%d: could not convert NCCOMPRESSION value", fname, line);
        } else
            enkf_quit("%s, l.%d: unknown token \"%s\"", fname, line, token);
    }                           /* while */

    if (prm->alpha < 0.0 || prm->alpha > 1.0)
        enkf_quit("ALPHA = %.3g outside [0,1] interval", prm->alpha);

    if (prm->nregions == 0) {
        region* r = malloc(sizeof(region));

        prm->nregions = 1;
        prm->regions = r;
        r->name = strdup("Global");
        r->x1 = -999;
        r->x2 = 999;
        r->y1 = -999;
        r->y2 = 999;
        r->nzints = 0;
        r->zints = NULL;
    }

    if (prm->nzints == 0) {
        prm->nzints = 3;
        prm->zints = malloc(3 * sizeof(zint));
        prm->zints[0].z1 = 0.0;
        prm->zints[0].z2 = DEPTH_SHALLOW;
        prm->zints[1].z1 = DEPTH_SHALLOW;
        prm->zints[1].z2 = DEPTH_DEEP;
        prm->zints[2].z1 = DEPTH_DEEP;
        prm->zints[2].z2 = DEPTH_MAX;
    }

    for (i = 0; i < prm->nregions; ++i) {
        region* r = &prm->regions[i];

        if (r->nzints == 0) {
            int j;

            r->nzints = prm->nzints;
            r->zints = malloc(r->nzints * sizeof(zint));
            for (j = 0; j < r->nzints; ++j) {
                r->zints[j].z1 = prm->zints[j].z1;
                r->zints[j].z2 = prm->zints[j].z2;
            }
        }
    }

    for (i = 0; i < prm->nbadbatchspecs; ++i) {
        badbatchspec* bb = &prm->badbatchspecs[i];

        if (bb->maxmad < fabs(bb->maxbias))
            bb->maxmad = fabs(bb->maxbias);
    }

    fclose(f);
    enkfprm_check(prm);

    return prm;
}

/**
 */
void enkfprm_destroy(enkfprm* prm)
{
    int i;

    free(prm->obsprm);
    free(prm->modelprm);
    free(prm->gridprm);

    free(prm->obstypeprm);

    free(prm->date);
    free(prm->ensdir);
    if (prm->bgdir != NULL)
        free(prm->bgdir);
    if (prm->nregions > 0) {
        for (i = 0; i < prm->nregions; ++i) {
            free(prm->regions[i].name);
            free(prm->regions[i].zints);
        }
        free(prm->regions);
    }
    if (prm->nzints != 0)
        free(prm->zints);
    if (prm->nplogs > 0) {
        for (i = 0; i < prm->nplogs; ++i)
            if (prm->plogs[i].gridname != NULL)
                free(prm->plogs[i].gridname);
        free(prm->plogs);
    }
    if (prm->nbadbatchspecs > 0) {
        for (i = 0; i < prm->nbadbatchspecs; ++i)
            free(prm->badbatchspecs[i].obstype);

        free(prm->badbatchspecs);
    }

    free(prm);
}

/**
 */
void enkfprm_print(enkfprm* prm, char offset[])
{
    int i;

    if (prm->mode == MODE_NONE)
        enkf_printf("%sMODE = NONE\n", offset);
    else if (prm->mode == MODE_ENKF)
        enkf_printf("%sMODE = EnKF\n", offset);
    else if (prm->mode == MODE_ENOI)
        enkf_printf("%sMODE = EnOI\n", offset);
    if (prm->mode == MODE_ENKF) {
        if (prm->scheme == SCHEME_NONE)
            prm->scheme = SCHEME_DENKF;
        if (prm->scheme == SCHEME_DENKF)
            enkf_printf("%sSCHEME = DEnKF\n", offset);
        else if (prm->scheme == SCHEME_ETKF)
            enkf_printf("%sSCHEME = ETKF\n", offset);
        if (1.0 - prm->alpha > 1.0e-6)
            enkf_printf("%s  ALPHA = %.3g\n", offset, prm->alpha);
    }
    enkf_printf("%sMODEL PRM = \"%s\"\n", offset, prm->modelprm);
    enkf_printf("%sGRID PRM = \"%s\"\n", offset, prm->gridprm);

    enkf_printf("%sOBS TYPES PRM = \"%s\"\n", offset, prm->obstypeprm);

    enkf_printf("%sDATE = \"%s\"\n", offset, prm->date);
    enkf_printf("%sOBS PRM = \"%s\"\n", offset, prm->obsprm);

    if (prm->mode == MODE_ENOI)
        enkf_printf("%sBGDIR = \"%s\"\n", offset, prm->bgdir);
    if (prm->mode == MODE_ENKF || !enkf_fstatsonly) {
        enkf_printf("%sENSEMBLE DIR = \"%s\"\n", offset, prm->ensdir);
        if (prm->enssize > 0)
            enkf_printf("%sENSEMBLE SIZE = %d\n", offset, prm->enssize);
        else
            enkf_printf("%sENSEMBLE SIZE = <FULL>\n", offset);
    }
    if (!enkf_fstatsonly) {
        enkf_printf("%sRFACTOR BASE = %.1f\n", offset, prm->rfactor_base);
        enkf_printf("%sINFLATION BASE = %.4f\n", offset, prm->inflation_base);
        if (!isnan(prm->inf_ratio))
            enkf_printf("%sINFLATION MODE = SPREAD LIMITED, MAX RATIO = %.2f\n", offset, prm->inf_ratio);
        else
            enkf_printf("%sINFLATION MODE = PLAIN\n", offset);
        if (isfinite(prm->kfactor))
            enkf_printf("%sKFACTOR = %.1f\n", offset, prm->kfactor);
        else
            enkf_printf("%sKFACTOR = n/a\n", offset);
        enkf_printf("%sLOCRAD = %.0f\n", offset, prm->locrad);
        enkf_printf("%sSTRIDE = %d\n", offset, prm->stride);
        enkf_printf("%sFIELDBUFFERSIZE = %d\n", offset, prm->fieldbufsize);
    }
    for (i = 0; i < prm->nregions; ++i) {
        region* r = &prm->regions[i];
        int j;

        enkf_printf("%sREGION %s: x = [%.1f, %.1f], y = [%.1f, %.1f], z intervals = ", offset, r->name, r->x1, r->x2, r->y1, r->y2);
        for (j = 0; j < r->nzints; ++j)
            enkf_printf("[%.0f %.0f] ", r->zints[j].z1, r->zints[j].z2);
        enkf_printf("\n");
    }
    if (!enkf_fstatsonly) {
        for (i = 0; i < prm->nplogs; ++i) {
            pointlog* plog = &prm->plogs[i];

            if (plog->gridname == NULL)
                enkf_printf("%sPOINTLOG %d %d\n", offset, plog->i, plog->j);
            else
                enkf_printf("%sPOINTLOG %d %d %s\n", offset, plog->i, plog->j, plog->gridname);
        }
    }
    enkf_printf("%sSOBSTRIDE = %d\n", offset, prm->sob_stride);
    for (i = 0; i < prm->nbadbatchspecs; ++i) {
        badbatchspec* bb = &prm->badbatchspecs[i];

        enkf_printf("%sBADBATCHES = %s %.3f %.3f %d\n", offset, bb->obstype, bb->maxbias, bb->maxmad, bb->minnobs);
    }
    if (prm->ncformat == NC_CLASSIC_MODEL)
        enkf_printf("%sNCFORMAT = CLASSIC\n", offset);
    else if (prm->ncformat == NC_64BIT_OFFSET)
        enkf_printf("%sNCFORMAT = 64BIT\n", offset);
    else if (prm->ncformat == NC_NETCDF4)
        enkf_printf("%sNCFORMAT = NETCDF4\n", offset);
    else
        enkf_quit("programming error");
    enkf_printf("%sNCCOMPRESSION = %d\n, prm->nccompression\n", offset, prm->nccompression);
    enkf_printflags(offset);
}

/**
 */
void enkfprm_describeprm(void)
{
    enkf_printf("\n");
    enkf_printf("  Main parameter file format:\n");
    enkf_printf("\n");
    enkf_printf("    MODE            = { ENKF | ENOI }\n");
    enkf_printf("  [ SCHEME          = { DENKF* | ETKF } ]\n");
    enkf_printf("  [ ALPHA           = <alpha> ]                              (1*)\n");
    enkf_printf("    MODEL           = <model prm file>\n");
    enkf_printf("    GRID            = <grid prm file>\n");
    enkf_printf("    OBSTYPES        = <obs. types prm file>\n");
    enkf_printf("    OBS             = <obs. data prm file>\n");
    enkf_printf("    DATE            = <julian day of analysis>\n");
    enkf_printf("  [ WINDOWMIN       = <start of obs window in days from analysis> ] (-inf*)\n");
    enkf_printf("  [ WINDOWMAX       = <end of obs window in days from analysis> ]   (+inf*)\n");
    enkf_printf("    ENSDIR          = <ensemble directory>\n");
    enkf_printf("    BGDIR           = <background directory>                 (MODE = ENOI)\n");
    enkf_printf("  [ KFACTOR         = <kfactor> ]                            (NaN*)\n");
    enkf_printf("  [ RFACTOR         = <rfactor> ]                            (1*)\n");
    enkf_printf("    ...\n");
    enkf_printf("    LOCRAD          = <loc. radius in km>\n");
    enkf_printf("  [ STRIDE          = <stride> ]                             (1*)\n");
    enkf_printf("  [ SOBSTRIDE       = <stride> ]                             (1*)\n");
    enkf_printf("  [ FIELDBUFFERSIZE = <fieldbuffersize> ]                    (1*)\n");
    enkf_printf("  [ INFLATION       = <inflation> [ <VALUE>* | PLAIN ]       (1*)\n");
    enkf_printf("    ...\n");
    enkf_printf("  [ ZSTATINTS       = [<z1> <z2>] ... ]\n");
    enkf_printf("  [ REGION          = <name> <lon1> <lon2> <lat1> <lat2> [[<z1> <z2>] ... ]\n");
    enkf_printf("    ...\n");
    enkf_printf("  [ POINTLOG        <i> <j> [grid name]]\n");
    enkf_printf("    ...\n");
    enkf_printf("  [ EXITACTION      = { BACKTRACE* | SEGFAULT } ]\n");
    enkf_printf("  [ BADBATCHES      = <obstype> <max. bias> <max. mad> <min # obs.> ]\n");
    enkf_printf("  [ NCFORMAT        = { CLASSIC | 64BIT | NETCDF4 } ]        (64BIT*)\n");
    enkf_printf("  [ NCCOMPRESSION   = <compression level>                    (0*) ]\n");
    enkf_printf("    ...\n");
    enkf_printf("\n");
    enkf_printf("  Notes:\n");
    enkf_printf("    1. { ... | ... | ... } denotes the list of possible choices\n");
    enkf_printf("    2. [ ... ] denotes an optional input\n");
    enkf_printf("    3. ( ... ) is a note\n");
    enkf_printf("    4. * denotes the default value\n");
    enkf_printf("    5. < ... > denotes a description of an entry\n");
    enkf_printf("    6. ... denotes repeating the previous item an arbitrary number of times\n");
    enkf_printf("\n");
}
