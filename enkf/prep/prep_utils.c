/******************************************************************************
 *
 * File:        prep_utils.c        
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
#include <assert.h>
#include "stringtable.h"
#include "ncw.h"
#include "kdtree.h"
#include "definitions.h"
#include "utils.h"
#include "ncutils.h"
#include "obsprm.h"
#include "grid.h"
#include "model.h"
#include "observations.h"
#include "allreaders.h"
#include "prep_utils.h"

#define NINC 5
#define QCFLAGVALMAX 31

/**
 */
static int obs_badob(observations* obs, int i)
{
    observation* o = &obs->data[i];

    if (o->status != STATUS_OK)
        return 0;
    if (o->type < 0 || o->product < 0 || o->instrument < 0 || o->fid < 0 || o->batch < 0 || !isfinite(o->value) || fabs(o->value) > MAXOBSVAL || isnan(o->estd) || o->estd <= 0.0 || !isfinite(o->fi) || !isfinite(o->fj) || !isfinite(o->fk) || !isfinite(o->lon) || !isfinite(o->lat) || !isfinite(o->depth))
        return 1;
    return 0;
}

/**
 */
static void readobs(obsmeta* meta, model* m, obsread_fn reader, observations* obs)
{
    int otid = obstype_getid(obs->nobstypes, obs->obstypes, meta->type, 1);
    obstype* ot = &obs->obstypes[otid];
    grid* g = model_getgridbyid(m, ot->gridid);
    int nfiles;
    char** fnames = NULL;
    int i;

    nfiles = 0;
    get_obsfiles(meta, &nfiles, &fnames);
    for (i = 0; i < nfiles; ++i) {
        int nobs0 = obs->nobs;
        int fid;

        enkf_printf("      reading %s:\n", fnames[i]);
        fid = st_add_ifabsent(obs->datafiles, fnames[i], -1);
        reader(fnames[i], fid, meta, g, obs);

        enkf_printf("        # good obs = %d\n", obs->nobs - nobs0);
        enkf_flush();
        free(fnames[i]);
    }
    free(fnames);
}

/** Add observations from a certain provider.
 *  This procedure contains generic/common operations done after reading the
 *  data.
 */
void obs_add(observations* obs, model* m, obsmeta* meta, int nexclude, obsregion* exclude)
{
    int nobs0 = obs->nobs;
    int otid = obstype_getid(obs->nobstypes, obs->obstypes, meta->type, 1);
    obstype* ot = &obs->obstypes[otid];
    int applylog = model_getvarislog(m, model_getvarid(m, ot->varnames[0], 1));

    grid* g = model_getgridbyid(m, ot->gridid);
    float** depth = grid_getdepth(g);
    int** numlevels = grid_getnumlevels(g);
    int isperiodic_i = grid_isperiodic_i(g);
    double lonbase = grid_getlonbase(g);

    double mindepth = NAN;
    double maxdepth = NAN;
    double footprint = 0.0;
    double varshift = 0.0;
    int thin = 0;
    obsread_fn reader;
    int i, ngood, npars;

    enkf_printf("    PRODUCT = %s, TYPE = %s, reader = %s\n", meta->product, meta->type, meta->reader);
    st_add_ifabsent(obs->products, meta->product, -1);

    /*
     * search for and take notice of parameters common for all readers, and
     * remove them from the parameter list
     */
    for (i = 0; i < meta->npars; ++i) {
        if (strcasecmp(meta->pars[i].name, "MINDEPTH") == 0) {
            if (!str2double(meta->pars[i].value, &mindepth))
                enkf_quit("observation prm file: can not convert MINDEPTH = \"%s\" to double\n", meta->pars[i].value);
            else {
                free(meta->pars[i].name);
                free(meta->pars[i].value);
                meta->pars[i].name = NULL;
                meta->pars[i].value = NULL;
            }
        } else if (strcasecmp(meta->pars[i].name, "MAXDEPTH") == 0) {
            if (!str2double(meta->pars[i].value, &maxdepth))
                enkf_quit("observation prm file: can not convert MAXDEPTH = \"%s\" to double\n", meta->pars[i].value);
            else {
                free(meta->pars[i].name);
                free(meta->pars[i].value);
                meta->pars[i].name = NULL;
                meta->pars[i].value = NULL;
            }
        } else if (strcasecmp(meta->pars[i].name, "FOOTPRINT") == 0) {
            if (!str2double(meta->pars[i].value, &footprint))
                enkf_quit("observation prm file: can not convert FOOTPRINT = \"%s\" to double\n", meta->pars[i].value);
            else {
                free(meta->pars[i].name);
                free(meta->pars[i].value);
                meta->pars[i].name = NULL;
                meta->pars[i].value = NULL;
            }
        } else if (strcasecmp(meta->pars[i].name, "VARSHIFT") == 0) {
            if (!str2double(meta->pars[i].value, &varshift))
                enkf_quit("observation prm file: can not convert VARSHIFT = \"%s\" to double\n", meta->pars[i].value);
            else {
                free(meta->pars[i].name);
                free(meta->pars[i].value);
                meta->pars[i].name = NULL;
                meta->pars[i].value = NULL;
            }
        } else if (strcasecmp(meta->pars[i].name, "THIN") == 0) {
            if (!str2int(meta->pars[i].value, &thin))
                enkf_quit("observation prm file: can not convert THIN = \"%s\" to double\n", meta->pars[i].value);
            else {
                free(meta->pars[i].name);
                free(meta->pars[i].value);
                meta->pars[i].name = NULL;
                meta->pars[i].value = NULL;
            }
        }
    }
    /*
     * compact parameters
     */
    npars = 0;
    for (i = 0; i < meta->npars; ++i) {
        if (meta->pars[i].name != NULL) {
            meta->pars[npars].name = meta->pars[i].name;
            meta->pars[npars].value = meta->pars[i].value;
            npars++;
        }
    }
    meta->npars = npars;

    reader = get_obsreadfn(meta);
    readobs(meta, m, reader, obs);      /* add the data to obs */

    if (!isnan(lonbase)) {
        for (i = nobs0; i < obs->nobs; ++i) {
            observation* o = &obs->data[i];

            o->id_orig = i;
            if (o->lon < lonbase)
                o->lon += 360.0;
            else if (o->lon >= lonbase + 360.0)
                o->lon -= 360.0;
        }
    }

    /*
     * common checks
     */
    if (obs->nobs - nobs0 > 0) {
        int ninf = 0;
        int nmin = 0;
        int nmax = 0;
        int noutow = 0;
        int noutod = 0;
        int nland = 0;
        int nshallow = 0;
        int nthin = 0;
        int nexcluded = 0;
        int ni, nj, ksurf, n;

        enkf_printf("      id = %d - %d\n", nobs0, obs->nobs - 1);
        grid_getsize(g, &ni, &nj, NULL);
        ksurf = grid_getsurflayerid(g);
        for (i = nobs0; i < obs->nobs; ++i) {
            observation* o = &obs->data[i];

            if (o->status != STATUS_OK)
                continue;
            if (o->time - obs->da_time < ot->obswindow_min) {
                o->status = STATUS_OUTSIDEOBSWINDOW;
                noutow++;
                continue;
            }
            if (o->time - obs->da_time >= ot->obswindow_max) {
                o->status = STATUS_OUTSIDEOBSWINDOW;
                noutow++;
                continue;
            }
            o->value += varshift;
            if (!isfinite(o->value)) {
                o->status = STATUS_RANGE;
                ninf++;
                continue;
            }
            if (o->value < ot->allowed_min) {
                o->status = STATUS_RANGE;
                nmin++;
                continue;
            }
            if (o->value > ot->allowed_max) {
                o->status = STATUS_RANGE;
                nmax++;
                continue;
            }
            if (o->lon <= ot->xmin || o->lon >= ot->xmax || o->lat <= ot->ymin || o->lat >= ot->ymax || o->depth < ot->zmin || o->depth > ot->zmax) {
                o->status = STATUS_OUTSIDEOBSDOMAIN;
                noutod++;
                continue;
            }
            if (depth != NULL) {
                o->model_depth = (double) interpolate2d(o->fi, o->fj, ni, nj, depth, numlevels, isperiodic_i);
                if (!isfinite(o->model_depth) || o->model_depth == 0) {
                    o->status = STATUS_LAND;
                    nland++;
                    continue;
                }
            } else if (island(o->fi, o->fj, o->fk, ni, nj, ksurf, numlevels, isperiodic_i)) {
                o->status = STATUS_LAND;
                nland++;
                continue;
            }
            if (isfinite(mindepth)) {
                if (depth == NULL)
                    enkf_quit("MINDEPTH specified for the obs reader, but no depth specified for grid \"%s\"", grid_getname(g));
                if (o->model_depth < mindepth || !isfinite(o->model_depth)) {
                    o->status = STATUS_SHALLOW;
                    nshallow++;
                }
            }
            if (isfinite(maxdepth)) {
                if (depth == NULL)
                    enkf_quit("MAXDEPTH specified for the obs reader, but no depth specified for grid \"%s\"", grid_getname(g));
                if (o->model_depth > maxdepth || !isfinite(o->model_depth)) {
                    o->status = STATUS_SHALLOW;
                    nshallow++;
                }
            }
            if (thin > 1 && i % thin != 0) {
                o->status = STATUS_THINNED;
                nthin++;
                continue;
            }

            if (footprint > 0.0)
                o->footprint = footprint;
            else
                o->footprint = 0.0;
        }

        for (n = 0; n < nexclude; ++n) {
            obsregion* r = &exclude[n];

            if (r->otid == -1 || r->otid == otid) {
                for (i = nobs0; i < obs->nobs; ++i) {
                    observation* o = &obs->data[i];

                    if (o->status != STATUS_OK)
                        continue;
                    if (o->lon >= r->x1 && o->lon <= r->x2 && o->lat >= r->y1 && o->lat <= r->y2) {
                        o->status = STATUS_EXCLUDED;
                        nexcluded++;
                    }
                }
            }
        }

        if (applylog) {
            for (i = nobs0; i < obs->nobs; ++i) {
                observation* o = &obs->data[i];

                if (o->status != STATUS_OK)
                    continue;
                o->value = log10(o->value);
                if (!isnormal(o->value)) {
                    o->status = STATUS_RANGE;
                    ninf++;
                }
            }
        }

        if (noutow > 0)
            enkf_printf("        %d observations outside obs. window\n", noutow);
        if (noutod > 0)
            enkf_printf("        %d obbservations outside obs. domain\n", noutod);
        if (ninf > 0)
            enkf_printf("        %d observations not finite\n", ninf);
        if (nmin > 0)
            enkf_printf("        %d observations below allowed minimum of %.4g\n", nmin, ot->allowed_min);
        if (nmax > 0)
            enkf_printf("        %d observations above allowed maximum of %.4g\n", nmax, ot->allowed_max);
        if (nland > 0)
            enkf_printf("        %d observations on land\n", nland);
        if (nshallow > 0)
            enkf_printf("        %d observations in shallow areas\n", nshallow);
        if (nthin > 0)
            enkf_printf("        %d observations thinned\n", nthin);
        if (nexcluded > 0)
            enkf_printf("        %d observations in excluded regions\n", nexcluded);
    }

    obs->compacted = 0;
    enkf_printf("      total %d observations\n", obs->nobs - nobs0);
    for (ngood = 0, i = nobs0; i < obs->nobs; ++i)
        if (obs->data[i].status == STATUS_OK)
            ngood++;
    enkf_printf("      %d valid observations\n", ngood);

    if (obs->nobs - nobs0 > 0 && applylog && meta->nestds == 0)
        enkf_quit("%s: observation error must be specified explicitly for observations of type associated with log-transformed model variables", ot->name);

    /*
     * add specified errors 
     */
    if (obs->nobs - nobs0 > 0 && meta->nestds > 0) {
        int i, o;

        for (i = 0; i < meta->nestds; ++i) {
            metastd* estd = &meta->estds[i];

            if (estd->type == STDTYPE_VALUE) {
                double v = ((double*) estd->data)[0];

                if (estd->op == ARITHMETIC_EQ)
                    enkf_printf("      setting error_std to %.3g\n", v);
                else if (estd->op == ARITHMETIC_PLUS)
                    enkf_printf("      adding error_std of %.3g\n", v);
                else if (estd->op == ARITHMETIC_MULT)
                    enkf_printf("      multiplying error_std by %.3g\n", v);
                else if (estd->op == ARITHMETIC_MIN)
                    enkf_printf("      setting minimum error_std to %.3g\n", v);
                else if (estd->op == ARITHMETIC_MAX)
                    enkf_printf("      setting maximum error_std to %.3g\n", v);
                else
                    enkf_quit("programming error");

                for (o = nobs0; o < obs->nobs; ++o) {
                    observation* oo = &obs->data[o];

                    if (oo->status != STATUS_OK)
                        continue;
                    if (estd->op == ARITHMETIC_EQ)
                        oo->estd = v;
                    else if (estd->op == ARITHMETIC_PLUS)
                        oo->estd = sqrt(oo->estd * oo->estd + v * v);
                    else if (estd->op == ARITHMETIC_MULT)
                        oo->estd *= v;
                    else if (estd->op == ARITHMETIC_MIN)
                        oo->estd = (oo->estd < v) ? v : oo->estd;
                    else if (estd->op == ARITHMETIC_MAX)
                        oo->estd = (oo->estd > v) ? v : oo->estd;
                    else
                        enkf_quit("programming error");
                }
            } else if (estd->type == STDTYPE_FILE) {
                char* fname = (char*) estd->data;
                int ni, nj, nk;

                if (estd->op == ARITHMETIC_EQ)
                    enkf_printf("      setting error_std to %s from %s\n", estd->varname, fname);
                else if (estd->op == ARITHMETIC_PLUS)
                    enkf_printf("      adding error_std of %s from %s\n", estd->varname, fname);
                else if (estd->op == ARITHMETIC_MULT)
                    enkf_printf("      multiplying error_std by %s from %s\n", estd->varname, fname);
                else if (estd->op == ARITHMETIC_MIN)
                    enkf_printf("      setting minimum error_std to %s from %s\n", estd->varname, fname);
                else if (estd->op == ARITHMETIC_MAX)
                    enkf_printf("      setting maximum error_std to %s from %s\n", estd->varname, fname);
                else
                    enkf_quit("programming error");

                grid_getsize(g, &ni, &nj, &nk);

                if (ot->issurface) {
                    float** v = alloc2d(nj, ni, sizeof(float));

                    ncu_readfield(fname, estd->varname, 0, ni, nj, nk, v[0]);
                    for (o = nobs0; o < obs->nobs; ++o) {
                        observation* oo = &obs->data[o];
                        float vv;

                        if (oo->status != STATUS_OK)
                            continue;

                        vv = (float) interpolate2d(oo->fi, oo->fj, ni, nj, v, numlevels, isperiodic_i);
                        if (estd->op == ARITHMETIC_EQ)
                            oo->estd = vv;
                        else if (estd->op == ARITHMETIC_PLUS)
                            oo->estd = sqrt(oo->estd * oo->estd + vv * vv);
                        else if (estd->op == ARITHMETIC_MULT)
                            oo->estd *= vv;
                        else if (estd->op == ARITHMETIC_MIN)
                            oo->estd = (oo->estd < vv) ? vv : oo->estd;
                        else if (estd->op == ARITHMETIC_MAX)
                            oo->estd = (oo->estd > vv) ? vv : oo->estd;
                        else
                            enkf_quit("programming error");
                    }
                    free(v);
                } else {
                    float*** v = alloc3d(nk, nj, ni, sizeof(float));
                    int ksurf = grid_getsurflayerid(g);

                    ncu_read3dfield(fname, estd->varname, ni, nj, nk, v[0][0]);
                    for (o = nobs0; o < obs->nobs; ++o) {
                        observation* oo = &obs->data[o];
                        float vv;

                        if (oo->status != STATUS_OK)
                            continue;

                        vv = (float) interpolate3d(oo->fi, oo->fj, oo->fk, ni, nj, nk, ksurf, v, numlevels, isperiodic_i);
                        if (estd->op == ARITHMETIC_EQ)
                            oo->estd = vv;
                        else if (estd->op == ARITHMETIC_PLUS)
                            oo->estd = sqrt(oo->estd * oo->estd + vv * vv);
                        else if (estd->op == ARITHMETIC_MULT)
                            oo->estd *= vv;
                        else if (estd->op == ARITHMETIC_MIN)
                            oo->estd = (oo->estd < vv) ? vv : oo->estd;
                        else if (estd->op == ARITHMETIC_MAX)
                            oo->estd = (oo->estd > vv) ? vv : oo->estd;
                        else
                            enkf_quit("programming error");
                    }
                    free(v);
                }
            }
        }
    }

    /*
     * report time range 
     */
    if (obs->nobs - nobs0 > 0) {
        double day_min = DBL_MAX;
        double day_max = -DBL_MAX;
        int i;

        for (i = nobs0; i < obs->nobs; ++i) {
            observation* o = &obs->data[i];

            if (!isnan(o->time))
                o->time -= obs->da_time;
            else
                o->time = 0.0;
            if (o->status != STATUS_OK)
                continue;
            if (o->time < day_min)
                day_min = o->time;
            if (o->time > day_max)
                day_max = o->time;
        }
        if (day_min <= day_max) {
            enkf_printf("      min day = %.3f\n", day_min);
            enkf_printf("      max day = %.3f\n", day_max);
        }
    }
    fflush(stdout);

    for (i = nobs0; i < obs->nobs; ++i) {
        if (obs_badob(obs, i)) {
            enkf_printf("        bad observation detected: ");
            obs_printob(obs, i);
            enkf_quit("bad observation");
        }
    }
}

/**
 */
int obs_checkforland(observations* obs, model* m)
{
    int hasland = 0;
    int ni = -1, nj = -1, ksurf = -1, isperiodic = -1;
    int** numlevels = NULL;
    int type_prev, i;

    for (i = 0, type_prev = -1; i < obs->nobs; ++i) {
        observation* o = &obs->data[i];

        if (o->type != type_prev) {
            obstype* ot = &obs->obstypes[o->type];
            int vid = model_getvarid(m, ot->varnames[0], 1);
            grid* g = model_getvargrid(m, vid);

            grid_getsize(g, &ni, &nj, NULL);
            ksurf = grid_getsurflayerid(g);
            numlevels = grid_getnumlevels(g);
            isperiodic = grid_isperiodic_i(g);
        }

        if (island(o->fi, o->fj, o->fk, ni, nj, ksurf, numlevels, isperiodic)) {
            o->status = STATUS_LAND;
            hasland = 1;
        }
        type_prev = o->type;
    }

    return hasland;
}

/**
 */
void get_obsfiles(obsmeta* meta, int* nfiles, char*** fnames)
{
    int i;

    for (i = 0; i < meta->nfiles; ++i)
        find_files(meta->fnames[i], nfiles, fnames);
}

/**
 */
void print_obsstats(observations* obs, observations* sobs)
{
    int i;

    enkf_printf("    type    #used    #dropped #out_grd #out_obs #out_wnd #land    #shallow #badbatch#badvalue#thinned #excluded#superobs\n");
    enkf_printf("    -----------------------------------------------------------------------------------------------------------\n");
    for (i = 0; i < obs->nobstypes; ++i) {
        obstype* ot = &obs->obstypes[i];

        enkf_printf("    %-7s %-8d %-8d %-8d %-8d %-8d %-8d %-8d %-8d %-8d %-8d %-8d %-8d\n", ot->name, ot->ngood, ot->nobs - ot->ngood, ot->noutside_grid, ot->noutside_obsdomain, ot->noutside_obswindow, ot->nland, ot->nshallow, ot->nbadbatch, ot->nrange, ot->nthinned, ot->nexcluded, sobs->obstypes[i].nobs);
    }
    enkf_printf("    total   %-8d %-8d %-8d %-8d %-8d %-8d %-8d %-8d %-8d %-8d %-8d %-8d\n", obs->ngood, obs->nobs - obs->ngood, obs->noutside_grid, obs->noutside_obsdomain, obs->noutside_obswindow, obs->nland, obs->nshallow, obs->nbadbatch, obs->nrange, obs->nthinned, obs->nexcluded, sobs->nobs);
}

#define NLONNAMES 3
char* LONNAMES[] = { "lon",
    "longitude",
    "LONGITUDE"
};

/**
 */
char* get_lonname(int ncid, char* lonname)
{
    int i;

    if (lonname != NULL) {
        ncw_inq_varid(ncid, lonname, &i);
        return lonname;
    }
    for (i = 0; i < NLONNAMES; ++i)
        if (ncw_var_exists(ncid, LONNAMES[i]))
            return LONNAMES[i];

    return NULL;
}

#define NLATNAMES 3
static char* LATNAMES[] = { "lat",
    "latitude",
    "LATITUDE"
};

/**
 */
char* get_latname(int ncid, char* latname)
{
    int i;

    if (latname != NULL) {
        ncw_inq_varid(ncid, latname, &i);
        return latname;
    }
    for (i = 0; i < NLATNAMES; ++i)
        if (ncw_var_exists(ncid, LATNAMES[i]))
            return LATNAMES[i];

    return NULL;
}

#define NZNAMES 3
static char* ZNAMES[] = { "z",
    "depth",
    "DEPTH"
};

/**
 */
char* get_zname(int ncid, char* zname)
{
    int i;

    if (zname != NULL) {
        ncw_inq_varid(ncid, zname, &i);
        return zname;
    }
    for (i = 0; i < NZNAMES; ++i)
        if (ncw_var_exists(ncid, ZNAMES[i]))
            return ZNAMES[i];

    return NULL;
}

/**
 */
static void get_timenames(int ncid, int* varids, char** timenames)
{
    int i, nfound, nvar;

    if (timenames[0] != NULL) {
        if (timenames[1] != NULL)
            ncw_inq_varid(ncid, timenames[1], &varids[1]);
        ncw_inq_varid(ncid, timenames[0], &varids[0]);
        nfound = (timenames[1] == NULL) ? 1 : 2;
    } else {
        ncw_inq_nvars(ncid, &nvar);
        for (i = 0, nfound = 0; i < nvar; ++i) {
            char varname[NC_MAX_NAME];

            ncw_inq_varname(ncid, i, varname);
            if (strncasecmp(varname, "time", MAXSTRLEN - 1) == 0) {
                if (nfound == 1)
                    enkf_quit("%s: more than 1 possible time variables found: %s, %s. Use entry TIMENAMES to specify offset and difference time variable names", ncw_get_path(ncid), timenames[0], varname);
                varids[nfound] = i;
                timenames[nfound] = strdup(varname);
                nfound++;
            }
        }
    }

    if (nfound == 0)
        enkf_quit("%s: no time variable found", ncw_get_path(ncid));
    if (nfound == 2) {
        size_t size[2];

        for (i = 0; i < 2; ++i)
            ncw_inq_varsize(ncid, varids[i], &size[i]);
        if (size[1] == 1)
            return;
        if (size[1] > 1 && size[0] == 1) {
            char* tmpc = timenames[0];
            int tmpi = varids[0];

            timenames[0] = timenames[1];
            timenames[1] = tmpc;
            varids[0] = varids[1];
            varids[1] = tmpi;

            return;
        }
        enkf_quit("%s: both possible time variables \"%s\" and \"%s\" have size greater than 1", ncw_get_path(ncid), timenames[0], timenames[1]);
    }
}

/** Returns scaling multiple to convert time units to days.
 */
static double tunits2days(char* tunits)
{
    double mult = NAN;

    if (strchr(tunits, ' ') != NULL)
        enkf_quit("tunits2days(): the input = \"%s\" is not supposed to contain spaces", tunits);
    if (strncasecmp(tunits, "sec", 3) == 0)
        mult = 1.0 / 86400.0;
    else if (strncasecmp(tunits, "hou", 3) == 0)
        mult = 1.0 / 24.0;
    else if (strncasecmp(tunits, "day", 3) == 0)
        mult = 1.0;
    else
        enkf_quit("can not convert \"%s\" to days", tunits);

    return mult;
}

/** Reads observation time. Can handle compound time represented as sum of
 ** a base time (of size 1) and offset time.
 * @param meta - meta data for the block of observations
 * @param ncid - input data file
 * @param size - output: length of the time vector
 * @param time - output: time vector
 */
void get_time(obsmeta* meta, int ncid, size_t* size, double** time)
{
    char* timenames[2] = { NULL, NULL };
    int varids[2] = { -1, -1 };
    char tunits[MAXSTRLEN];
    double tunits_multiple = NAN, tunits_offset = 0.0;
    int i;

    for (i = 0; i < meta->npars; ++i) {
        if (strcasecmp(meta->pars[i].name, "TIMENAME") == 0)
            timenames[0] = strdup(meta->pars[i].value);
        else if (strcasecmp(meta->pars[i].name, "TIMENAMES") == 0) {
            char seps[] = " ,";
            char* token;

            token = strtok(meta->pars[i].value, seps);
            if (token == NULL)
                enkf_quit("%s: no entries found after TIMENAMES", meta->prmfname);
            timenames[0] = strdup(token);
            token = strtok(NULL, seps);
            if (token != NULL)
                timenames[1] = strdup(token);
        }
    }
    get_timenames(ncid, varids, timenames);
    if (timenames[1] != NULL) {
        double offset = 0.0;

        ncw_get_att_text(ncid, varids[1], "units", tunits);
        tunits_convert(tunits, &tunits_multiple, &tunits_offset);
        ncu_readvardouble(ncid, varids[1], 1, &offset);
        tunits_offset = offset * tunits_multiple + tunits_offset;
    }
    ncw_get_att_text(ncid, varids[0], "units", tunits);
    if (timenames[1] != NULL)
        tunits_multiple = tunits2days(tunits);
    else
        tunits_convert(tunits, &tunits_multiple, &tunits_offset);
    ncw_inq_varsize(ncid, varids[0], size);
    assert(*size > 0);
    *time = malloc(*size * sizeof(size_t));
    ncu_readvardouble(ncid, varids[0], *size, *time);
    for (i = 0; i < *size; ++i)
        (*time)[i] = (*time)[i] * tunits_multiple + tunits_offset;

    if (timenames[0] != NULL)
        free(timenames[0]);
    if (timenames[1] != NULL)
        free(timenames[1]);
}

#define NINSTNAMES 5
static char* INSTNAMES[] = { "instrument",
    "Instrument",
    "INSTRUMENT",
    "platform",
    "Platform",
    "PLATFORM"
};

/**
 */
int get_insttag(int ncid, char* varname, char* insttag)
{
    int i;

    if (varname != NULL) {
        int varid;

        ncw_inq_varid(ncid, varname, &varid);
        for (i = 0; i < NINSTNAMES; ++i) {
            if (ncw_att_exists(ncid, varid, INSTNAMES[i])) {
                nc_type type;
                size_t len;

                ncw_inq_att(ncid, varid, INSTNAMES[i], &type, &len);
                if (type == NC_CHAR && len < MAXSTRLEN) {
                    ncw_get_att_text(ncid, varid, INSTNAMES[i], insttag);
                    return 1;
                }
            }
        }
    }
    for (i = 0; i < NINSTNAMES; ++i) {
        if (ncw_att_exists(ncid, NC_GLOBAL, INSTNAMES[i])) {
            nc_type type;
            size_t len;

            ncw_inq_att(ncid, NC_GLOBAL, INSTNAMES[i], &type, &len);
            if (type == NC_CHAR && len < MAXSTRLEN) {
                ncw_get_att_text(ncid, NC_GLOBAL, INSTNAMES[i], insttag);
                return 1;
            }
        }
    }

    return 0;
}

/**
 */
void get_qcflags(obsmeta* meta, int* nqcflagvars, char*** qcflagvarnames, uint32_t** qcflagmasks)
{
    int i;

    assert(*nqcflagvars == 0);
    assert(*qcflagvarnames == NULL);
    assert(*qcflagmasks == NULL);

    for (i = 0; i < meta->npars; ++i) {
        if (strcasecmp(meta->pars[i].name, "QCFLAGVARNAME") == 0 || strcasecmp(meta->pars[i].name, "QCFLAGNAME") == 0) {
            if (*nqcflagvars % NINC == 0)
                *qcflagvarnames = realloc(*qcflagvarnames, (*nqcflagvars + NINC) * sizeof(char*));
            (*qcflagvarnames)[*nqcflagvars] = meta->pars[i].value;

            /*
             * make sure that the entry "QCFLAGVARNAME" or "QCFLAGNAME" is
             * followed by "QCFLAGVALS"
             */
            if (meta->npars == i + 1 || strcasecmp(meta->pars[i + 1].name, "QCFLAGVALS") != 0)
                enkf_quit("%s: parameter \"%s\" must be followed by parameter \"QCFLAGVALS\"\n", meta->prmfname, meta->pars[i].name);
        } else if (strcasecmp(meta->pars[i].name, "QCFLAGVALS") == 0) {
            char seps[] = " ,";
            char* line = meta->pars[i].value;
            char* token;
            int val;
            int ii;

            if (i == 0 || (strcasecmp(meta->pars[i - 1].name, "QCFLAGVARNAME") != 0 && strcasecmp(meta->pars[i - 1].name, "QCFLAGNAME") != 0))
                enkf_quit("%s: parameter \"QCFLAGVALS\" must be preceeded by parameter \"QCFLAGVARNAME\"\n", meta->prmfname);

            if (*nqcflagvars % NINC == 0)
                *qcflagmasks = realloc(*qcflagmasks, (*nqcflagvars + NINC) * sizeof(uint32_t));

            (*qcflagmasks)[*nqcflagvars] = 0;
            while ((token = strtok(line, seps)) != NULL) {
                if (!str2int(token, &val))
                    enkf_quit("%s: could not convert QCFLAGVALS entry \"%s\" to integer", meta->prmfname, token);
                if (val < 0 || val > 31)
                    enkf_quit("%s: QCFLAGVALS entry = %d (supposed to be in [0,31] interval", meta->prmfname, val);
                (*qcflagmasks)[*nqcflagvars] |= 1 << val;
                line = NULL;
            }
            if ((*qcflagmasks)[*nqcflagvars] == 0)
                enkf_quit("%s: no valid flag entries found after QCFLAGVALS\n", meta->prmfname);
            enkf_printf("        QCFLAGVARNAME = %s\n", (*qcflagvarnames)[*nqcflagvars]);
            enkf_printf("        QCFLAG values used =");
            for (ii = 0; ii <= QCFLAGVALMAX; ++ii)
                if ((*qcflagmasks)[*nqcflagvars] & (1 << ii))
                    enkf_printf(" %d", ii);
            enkf_printf("\n");
            (*nqcflagvars)++;
        }
    }
}
