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
#include <math.h>
#include <values.h>
#include "nan.h"
#include "definitions.h"
#include "utils.h"
#include "enkfprm.h"

#define NASYNC_INC 10
#define NVAR_INC 10
#define NREGIONS_INC 10
#define NPLOGS_INC 10
#define NBBSPECS_INC 10

/**
 */
static void enkfprm_check(enkfprm* prm)
{
    if (prm->mode == MODE_NONE)
        enkf_quit("%s: MODE not specified", prm->fname);
#if defined(ENKF_PREP)
    if (prm->obsprm == NULL)
        enkf_quit("%s: OBS not specified", prm->fname);
#endif
    if (prm->date == NULL)
        enkf_quit("%s: DATE not specified", prm->fname);
    if (prm->modelprm == NULL)
        enkf_quit("%s: MODEL not specified", prm->fname);
    if (prm->gridprm == NULL)
        enkf_quit("%s: GRID not specified", prm->fname);
#if defined(ENKF_PREP) || defined(ENKF_CALC)
    if (prm->obstypeprm == NULL)
        enkf_quit("%s: OBSTYPES not specified", prm->fname);
#endif
#if defined(ENKF_CALC) || defined(ENKF_POST)
    if (prm->ensdir == NULL && (prm->mode == MODE_ENKF || !enkf_fstatsonly))
        enkf_quit("%s: ENSDIR not specified", prm->fname);
#endif
    if (prm->mode == MODE_ENOI && prm->bgdir == NULL)
        enkf_quit("%s: BGDIR must be specified for MODE = ENOI", prm->fname);
#if defined(ENKF_CALC)
    if (!isfinite(prm->locrad) && !enkf_fstatsonly)
        enkf_quit("%s: LOCRAD not specified", prm->fname);
#endif
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
    prm->target = TARGET_DEFAULT;
    prm->date = NULL;
    prm->modelprm = NULL;
    prm->gridprm = NULL;

    prm->obstypeprm = NULL;
    prm->obsprm = NULL;
    prm->enssize = -1;
    prm->kfactor = NaN;
    prm->inflations = NULL;
    prm->inflation_base = 1.0;
    prm->rfactor_base = 1.0;
    prm->inflation_base = 1.0;
    prm->locrad = NaN;
    prm->stride = 1;
    prm->fieldbufsize = 1;
    prm->sob_stride = 1;

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
        } else if (strcasecmp(token, "TARGET") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: TARGET not specified", fname, line);
            else {
                if (strcasecmp(token, "ANALYSIS") == 0)
                    prm->target = TARGET_ANALYSIS;
                else if (strcasecmp(token, "INCREMENT") == 0)
                    prm->target = TARGET_INCREMENT;
                else
                    enkf_quit("%s, l.%d: target \"%s\" is not implemented", fname, line, token);
            }
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
        } else if (strcasecmp(token, "VARNAMES") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: VARNAMES not specified", fname, line);
            do {
                if (prm->nvar % NVAR_INC == 0) {
                    prm->varnames = realloc(prm->varnames, (prm->nvar + NVAR_INC) * sizeof(char*));
                    prm->inflations = realloc(prm->inflations, (prm->nvar + NVAR_INC) * sizeof(double));
                    for (i = prm->nvar; i < prm->nvar + NVAR_INC; ++i)
                        prm->inflations[i] = NaN;
                }
                prm->varnames[prm->nvar] = strdup(token);
                prm->nvar++;
            } while ((token = strtok(NULL, seps)) != NULL);
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
                enkf_quit("%s, l.%d: could not convert RFACTOR value", fname, line);
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
            if (strcasecmp(token, "BASE") == 0) {
                if ((token = strtok(NULL, seps)) == NULL)
                    enkf_quit("%s, l.%d: INFLATION BASE not specified", fname, line);
                if (!str2double(token, &prm->inflation_base))
                    enkf_quit("%s, l.%d: could not convert \"%s\" to double", fname, line, token);
            } else {
                for (i = 0; i < prm->nvar; ++i)
                    if (strcmp(token, prm->varnames[i]) == 0)
                        break;
                if (i == prm->nvar)
                    enkf_quit("%s, l.%d: could not identify the variable \"%s\". (Make sure that the entries for INFLATION appear after VARNAMES entry.)", fname, line, token);
                if ((token = strtok(NULL, seps)) == NULL)
                    enkf_quit("%s, l.%d: INFLATION for \"%s\" not specified", fname, line, prm->varnames[i]);
                if (!str2double(token, &prm->inflations[i]))
                    enkf_quit("%s, l.%d: could not convert \"%s\" to double", fname, line, token);
            }
        } else if (strcasecmp(token, "REGION") == 0) {
            char* space;
            char* newtoken;

            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: REGION name not specified", fname, line);
            if (prm->nregions % NREGIONS_INC == 0)
                prm->regions = realloc(prm->regions, (prm->nregions + NREGIONS_INC) * sizeof(region));
            newtoken = token;
            while ((space = strchr(newtoken, '_')) != NULL) {
                *space = ' ';
                newtoken = space + 1;
            }
            prm->regions[prm->nregions].name = strdup(token);
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: minimal longitude not specified", fname, line);
            if (!str2double(token, &prm->regions[prm->nregions].x1))
                enkf_quit("%s, l.%d: could not convert \"%s\" to double", fname, line, token);
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: maximal longitude not specified", fname, line);
            if (!str2double(token, &prm->regions[prm->nregions].x2))
                enkf_quit("%s, l.%d: could not convert \"%s\" to double", fname, line, token);
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: minimal latitude not specified", fname, line);
            if (!str2double(token, &prm->regions[prm->nregions].y1))
                enkf_quit("%s, l.%d: could not convert \"%s\" to double", fname, line, token);
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: maximal latitude not specified", fname, line);
            if (!str2double(token, &prm->regions[prm->nregions].y2))
                enkf_quit("%s, l.%d: could not convert \"%s\" to double", fname, line, token);
            prm->nregions++;
        } else if (strcasecmp(token, "POINTLOG") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: \"I\" coordinate for POINTLOG not specified", fname, line);
            if (prm->nplogs % NPLOGS_INC == 0)
                prm->plogs = realloc(prm->plogs, (prm->nplogs + NPLOGS_INC) * sizeof(pointlog));
            if (!str2int(token, &prm->plogs[prm->nplogs].i))
                enkf_quit("%s, l.%d: could convert \"I\" entry", fname, line);
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: \"J\" coordinate for POINTLOG not specified", fname, line);
            if (!str2int(token, &prm->plogs[prm->nplogs].j))
                enkf_quit("%s, l.%d: could convert \"J\" entry", fname, line);
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
            if (prm->nbadbatchspecs % NBBSPECS_INC == 0)
                prm->badbatchspecs = realloc(prm->badbatchspecs, (prm->nbadbatchspecs + NBBSPECS_INC) * sizeof(badbatchspec));
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
        } else
            enkf_quit("%s, l.%d: unknown token \"%s\"", fname, line, token);
    }                           /* while */

    for (i = 0; i < prm->nvar; ++i)
        if (isnan(prm->inflations[i]))
            prm->inflations[i] = prm->inflation_base;
        else
            prm->inflations[i] *= prm->inflation_base;
    prm->inflation_base = NaN;  /* not to be used */

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
    if (prm->nvar > 0) {
        for (i = 0; i < prm->nvar; ++i)
            free(prm->varnames[i]);
        free(prm->varnames);
        free(prm->inflations);
    }
    if (prm->nregions > 0) {
        for (i = 0; i < prm->nregions; ++i)
            free(prm->regions[i].name);
        free(prm->regions);
    }
    if (prm->nplogs > 0)
        free(prm->plogs);
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
    }
    enkf_printf("%sMODEL PRM = \"%s\"\n", offset, prm->modelprm);
    enkf_printf("%sGRID PRM = \"%s\"\n", offset, prm->gridprm);

    enkf_printf("%sOBS TYPES PRM = \"%s\"\n", offset, prm->obsprm);

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
    if (prm->nvar > 0) {
        enkf_printf("%sVARNAMES =", offset);
        for (i = 0; i < prm->nvar; ++i)
            enkf_printf(" %s", prm->varnames[i]);
        enkf_printf("\n");
    }
    if (!enkf_fstatsonly) {
        enkf_printf("%sRFACTOR BASE = %.1f\n", offset, prm->rfactor_base);
        if (isfinite(prm->kfactor))
            enkf_printf("%sKFACTOR = %.1f\n", offset, prm->kfactor);
        else
            enkf_printf("%sKFACTOR = n/a\n", offset);
        enkf_printf("%sLOCRAD = %.0f\n", offset, prm->locrad);
        enkf_printf("%sSTRIDE = %d\n", offset, prm->stride);
        enkf_printf("%sFIELDBUFFERSIZE = %d\n", offset, prm->fieldbufsize);
        for (i = 0; i < prm->nvar; ++i)
            enkf_printf("%sINFLATION %s = %.3f\n", offset, prm->varnames[i], prm->inflations[i]);
    }
    for (i = 0; i < prm->nregions; ++i) {
        region* r = &prm->regions[i];

        enkf_printf("%sREGION %s: x = [%.1f, %.1f], y = [%.1f, %.1f]\n", offset, r->name, r->x1, r->x2, r->y1, r->y2);
    }
    if (!enkf_fstatsonly) {
        for (i = 0; i < prm->nplogs; ++i) {
            pointlog* plog = &prm->plogs[i];

            enkf_printf("%sPOINTLOG %d %d\n", offset, plog->i, plog->j);
        }
    }
    enkf_printf("%sSOBSTRIDE = %d\n", offset, prm->sob_stride);
    for (i = 0; i < prm->nbadbatchspecs; ++i) {
        badbatchspec* bb = &prm->badbatchspecs[i];

        enkf_printf("%sBADBATCHES = %s %.3f %.3f %d\n", offset, bb->obstype, bb->maxbias, bb->maxmad, bb->minnobs);
    }
    enkf_printflags(offset);
}

/**
 */
void enkfprm_describeprm(void)
{
    enkf_printf("\n");
    enkf_printf("  Main parameter file format:\n");
    enkf_printf("\n");
    enkf_printf("    MODE                = { ENKF | ENOI }\n");
    enkf_printf("  [ SCHEME              = { DENKF* | ETKF | EnKF-N } ]\n");
    enkf_printf("  [ TARGET              = { ANALYSIS* | INCREMENT } ]\n");
    enkf_printf("    MODEL               = <model prm file>\n");
    enkf_printf("    GRID                = <grid prm file>\n");
    enkf_printf("    OBSTYPES            = <obs. types prm file>\n");
    enkf_printf("    OBS                 = <obs. data prm file>\n");
    enkf_printf("    DATE                = <julian day of analysis>\n");
    enkf_printf("    ENSDIR              = <ensemble directory>\n");
    enkf_printf("    BGDIR               = <background directory>                 (MODE = ENOI)\n");
    enkf_printf("  [ KFACTOR             = <kfactor> ]                            (1*)\n");
    enkf_printf("  [ RFACTOR             = <rfactor> ]                            (1*)\n");
    enkf_printf("    ...\n");
    enkf_printf("    LOCRAD              = <locrad>\n");
    enkf_printf("  [ STRIDE              = <stride> ]                             (1*)\n");
    enkf_printf("  [ SOBSTRIDE           = <stride> ]                             (1*)\n");
    enkf_printf("  [ FIELDBUFFERSIZE     = <fieldbuffersize> ]                    (1*)\n");
    enkf_printf("  [ INFLATION BASE      = <inflation> ]                          (1*)\n");
    enkf_printf("  [ INFLATION <VARNAME> = <inflation> ]                          (1*)\n");
    enkf_printf("    ...\n");
    enkf_printf("  [ REGION              <name> { <lon1> <lon2> <lat1> <lat2> } ]\n");
    enkf_printf("    ...\n");
    enkf_printf("  [ POINTLOG            { <i> <j> } ]\n");
    enkf_printf("    ...\n");
    enkf_printf("  [ EXITACTION          = { BACKTRACE* | SEGFAULT } ]\n");
    enkf_printf("  [ BADBATCHES          = <obstype> <max. bias> <max. mad> <min # obs.> ]\n");
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
