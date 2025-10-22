/******************************************************************************
 *
 * File:        obsprm.c        
 *
 * Created:     12/2012
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description:
 *
 * Revisions:   18062018 PS Renamed from obsmeta.c

 *
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "definitions.h"
#include "utils.h"
#include "obsprm.h"

#define NINC 10

/**
 */
static void obssection_init(obssection* section)
{
    memset(section, 0, sizeof(obssection));
}

/**
 */
static void obssection_addfname(obssection* section, char fname[])
{
    if (section->nfiles % NINC == 0)
        section->fnames = realloc(section->fnames, (section->nfiles + NINC) * sizeof(char*));
    section->fnames[section->nfiles] = strdup(fname);
    section->nfiles++;
}

/**
 */
void obsprm_read(char fname[], int* nsection, obssection** sections, int* nexclude, obsregion** exclude)
{
    FILE* f = NULL;
    char* buf = NULL;
    size_t bufsize = 0;
    obssection* section = NULL;
    int line;
    int i, j;

    *nsection = 0;
    *sections = NULL;

    *nexclude = 0;
    *exclude = NULL;

    f = enkf_fopen(fname, "r");

    line = 0;
    while (getline(&buf, &bufsize, f) >= 0) {
        char seps[] = " =\t\n";
        char* token;

        line++;
        if (buf[0] == '#')
            continue;
        if ((token = strtok(buf, seps)) == NULL)
            continue;
        if (strcasecmp(token, "PRODUCT") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: PRODUCT not specified", fname, line);

            *sections = realloc(*sections, (*nsection + 1) * sizeof(obssection));
            section = &(*sections)[*nsection];
            obssection_init(section);
            section->id = *nsection;
            section->prmfname = strdup(fname);
            section->product = strdup(token);
            (*nsection)++;
            continue;
        }

        if (strcasecmp(token, "EXCLUDE") == 0) {
            obsregion* r = NULL;

            if (*nexclude % NINC == 0)
                *exclude = realloc(*exclude, (*nexclude + NINC) * sizeof(obsregion));
            r = &(*exclude)[*nexclude];
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: EXCLUDE observation type not specified", fname, line);
            r->otname = strdup(token);
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
            (*nexclude)++;
            continue;
        }

        if (section == NULL)
            enkf_quit("%s, l.%d: expected entry PRODUCT", fname, line);

        if (strcasecmp(token, "READER") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: READER not specified", fname, line);
            section->reader = strdup(token);
        } else if (strcasecmp(token, "TYPE") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: TYPE not specified", fname, line);
            section->type = strdup(token);
        } else if (strcasecmp(token, "FILE") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: FILE not specified", fname, line);
            obssection_addfname(section, token);
        } else if (strcasecmp(token, "ERROR_STD") == 0) {
            std_entry* now = NULL;
            double std;

            section->estds = realloc(section->estds, (section->nestds + 1) * sizeof(std_entry));
            now = &section->estds[section->nestds];

            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: STD not specified", fname, line);

            if (str2double(token, &std)) {
                if (!isfinite(std) || std < 0.0)
                    enkf_quit("%s, l.%d: invalid STD value", fname, line);
                now->type = STDTYPE_VALUE;
                now->data = malloc(sizeof(double));
                ((double*) now->data)[0] = std;
                now->varname = NULL;
            } else {
                now->type = STDTYPE_FILE;
                now->data = strdup(token);
                if ((token = strtok(NULL, seps)) == NULL)
                    enkf_quit("%s, l.%d: variable name not specified", fname, line);
                now->varname = strdup(token);
            }
            if ((token = strtok(NULL, seps)) == NULL)
                now->op = ARITHMETIC_EQ;
            else if (strncasecmp(token, "EQ", 2) == 0)
                now->op = ARITHMETIC_EQ;
            else if (strncasecmp(token, "PL", 2) == 0)
                now->op = ARITHMETIC_PLUS;
            else if (strncasecmp(token, "MU", 2) == 0)
                now->op = ARITHMETIC_MULT;
            else if (strncasecmp(token, "MI", 2) == 0)
                now->op = ARITHMETIC_MIN;
            else if (strncasecmp(token, "MA", 2) == 0)
                now->op = ARITHMETIC_MAX;
            else
                enkf_quit("%s, l.%d:, unknown operation", fname, line);
            section->nestds++;
        } else if (strcasecmp(token, "PARAMETER") == 0) {
            par_entry* now = NULL;
            int p;

            section->pars = realloc(section->pars, (section->npars + 1) * sizeof(par_entry));
            now = &section->pars[section->npars];

            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: parameter name not specified (expected: PARAMETER <name> = <value>)", fname, line);
            /*
             * check that this parameter has not already been set in this
             * section, to avoid ambiguity
             */
            /*
             * skip this check for reader "z" to allow alternatives for some
             * parameters, such as e.g. VARNAME and ZNAME -- the first one found
             * in the data file will be used
             */
            if (strcasecmp(section->reader, "z") != 0) {
                if (strcasecmp(token, "QCFLAGNAME") != 0 && strcasecmp(token, "QCFLAGVARNAME") != 0 && strcasecmp(token, "QCFLAGVALS") != 0 && strcasecmp(token, "EXCLUDEINST") != 0) {
                    for (p = 0; p < section->npars; ++p) {
                        if (strncasecmp(token, section->pars[p].name, MAXSTRLEN - 1) == 0)
                            enkf_quit("%s: l.%d: parameter \"%s\" has already been set in this section", fname, line, token);
                    }
                }
            }
            now->name = strdup(token);
            token = strtok(NULL, seps);
            if (token != NULL)
                now->value = strdup(token);
            else
                enkf_quit("%s, l.%d: parameter value not specified (expected: PARAMETER <name> = <value>)", fname, line);
            while ((token = strtok(NULL, seps)) != NULL) {
                now->value = realloc(now->value, strlen(now->value) + strlen(token) + 2);
                strcat(now->value, " ");
                strcat(now->value, token);
            }
            section->npars++;
        } else
            enkf_quit("%s, l.%d: unknown token \"%s\" (have you missed \"PARAMETER\"?)", fname, line, token);
    }

    fclose(f);
    if (buf != NULL)
        free(buf);

    /*
     * print summary 
     */
    for (i = 0; i < *nsection; ++i) {
        section = &(*sections)[i];

        enkf_printf("    PRODUCT = %s\n", section->product);
        if (section->reader == NULL) {
            section->reader = strdup("standard");
            enkf_printf("      (assumed) READER = %s\n", section->reader);
        } else
            enkf_printf("      READER = %s\n", section->reader);
        if (section->type == NULL)
            enkf_quit("%s: observation type not specified for product \"%s\"", fname, section->product);
        enkf_printf("      TYPE = %s\n", section->type);
        for (j = 0; j < section->nfiles; ++j)
            enkf_printf("        File: %s\n", section->fnames[j]);
        for (j = 0; j < section->nestds; ++j) {
            char operstr[MAXSTRLEN] = "";
            std_entry* std = &section->estds[j];

            if (std->op == ARITHMETIC_EQ)
                strcpy(operstr, "EQUAL");
            else if (std->op == ARITHMETIC_PLUS)
                strcpy(operstr, "PLUS");
            else if (std->op == ARITHMETIC_MULT)
                strcpy(operstr, "MULT");
            else if (std->op == ARITHMETIC_MIN)
                strcpy(operstr, "MIN");
            else if (std->op == ARITHMETIC_MAX)
                strcpy(operstr, "MAX");

            if (std->type == STDTYPE_VALUE)
                enkf_printf("      ERROR_STD = %.3g, operation = %s\n", ((double*) std->data)[0], operstr);
            else if (std->type == STDTYPE_FILE)
                enkf_printf("      ERROR_STD = %s %s, operation = %s\n", (char*) std->data, std->varname, operstr);
        }
        for (j = 0; j < section->npars; ++j)
            enkf_printf("      PARAMETER %s = %s\n", section->pars[j].name, section->pars[j].value);
    }

    for (i = 0; i < *nexclude; ++i) {
        obsregion* r = &(*exclude)[i];

        enkf_printf("    EXCLUDE: TYPE = %s, %.3f <= lon <= %.3f, %.3f <= lat <= %.3f\n", r->otname, r->x1, r->x2, r->y1, r->y2);
    }
}

/**
 */
void obsprm_destroy(int nsection, obssection sections[], int nexclude, obsregion exclude[])
{
    int i, j;

    for (i = 0; i < nsection; ++i) {
        obssection* section = &sections[i];

        free(section->prmfname);
        free(section->product);
        free(section->type);
        free(section->reader);
        if (section->nfiles > 0) {
            for (j = 0; j < section->nfiles; ++j)
                free(section->fnames[j]);
            free(section->fnames);
        }
        if (section->nestds > 0) {
            for (j = 0; j < section->nestds; ++j) {
                std_entry* std = &section->estds[j];

                if (std->data != NULL)
                    free(std->data);
                if (std->varname != NULL)
                    free(std->varname);
            }
            free(section->estds);
        }
        if (section->npars > 0) {
            for (j = 0; j < section->npars; ++j) {
                free(section->pars[j].name);
                free(section->pars[j].value);
            }
            free(section->pars);
        }
    }
    if (sections != NULL)
        free(sections);

    for (i = 0; i < nexclude; ++i)
        free(exclude[i].otname);
    if (exclude != NULL)
        free(exclude);
}

/**
 */
void obsprm_describeprm(void)
{
    enkf_printf("\n");
    enkf_printf("  Observation data parameter file format:\n");
    enkf_printf("\n");
    enkf_printf("    PRODUCT   = <product>\n");
    enkf_printf("    TYPE      = <observation type>\n");
    enkf_printf("    READER    = <reader>\n");
    enkf_printf("    FILE      = <data file wildcard> \n");
    enkf_printf("    ...\n");
    enkf_printf("  [ PARAMETER <name> = <value> ]\n");
    enkf_printf("    ...\n");
    enkf_printf("  [ ERROR_STD = { <value> | <data file> <varname> } [ EQ* | PL | MU | MI | MA ] ]\n");
    enkf_printf("    ...\n");
    enkf_printf("\n");
    enkf_printf("  [ <more of the above blocks> ]\n");
    enkf_printf("\n");
    enkf_printf("  [ EXCLUDE   = { <observation type> | ALL } <lon1> <lon2> <lat1> <lat2> ]\n");
    enkf_printf("    ...\n");
    enkf_printf("\n");
    enkf_printf("  Notes:\n");
    enkf_printf("    1. { ... | ... | ... } denotes the list of possible choices\n");
    enkf_printf("    2. [ ... ] denotes an optional input\n");
    enkf_printf("    3. * denotes the default value\n");
    enkf_printf("    4. < ... > denotes a description of an entry\n");
    enkf_printf("    5. ... denotes repeating the previous item an arbitrary number of times\n");
    enkf_printf("    6. Use \"--list-readers\" to list available readers\n");
    enkf_printf("    7. Use \"--describe-reader <reader>\" to list available parameters\n");
    enkf_printf("\n");
}
