/******************************************************************************
 *
 * File:        pointlog.c        
 *
 * Created:     7/10/2013
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
#include <assert.h>
#include <math.h>
#include "definitions.h"
#include "distribute.h"
#include "utils.h"
#include "ncutils.h"
#include "grid.h"
#include "dasystem.h"
#include "pointlog.h"
#include "version.h"

#define NPLOG_INC 10

/**
 */
void plogs_add(int* nplog, pointlog** plogs, double lon, double lat, char* gridname)
{
    pointlog* plog;

    if (*nplog % NPLOG_INC == 0)
        *plogs = realloc(*plogs, (*nplog + NPLOG_INC) * sizeof(pointlog));

    plog = &(*plogs)[*nplog];
    plog->id = *nplog;
    plog->lon = lon;
    plog->lat = lat;
    plog->gridid = -1;
    plog->gridname = (gridname != NULL) ? strdup(gridname) : NULL;
    plog->fij = NULL;
    (*nplog)++;
}

/**
 */
void plogs_destroy(int nplog, pointlog plogs[])
{
    int i;

    for (i = 0; i < nplog; ++i) {
        pointlog* plog = &plogs[i];

        if (plog->gridname != NULL)
            free(plog->gridname);
        if (plog->fij != NULL)
            free(plog->fij);
    }
    if (nplog > 0)
        free(plogs);
}

/**
 */
static void get_gridstr(dasystem* das, int gid, char str[])
{
    if (model_getngrid(das->m) == 1)
        strcpy(str, "");
    else
        sprintf(str, "-%d", gid);
}

#if defined(ENKF_CALC)

/** Create a pointlog file and write local obs to it.
 */
void das_createplog(dasystem* das, int plogid, int ploc, int* lobs, double* lcoeffs)
{
    pointlog* plog = &das->plogs[plogid];
    observations* obs = das->obs;
    int ngrid = model_getngrid(das->m);
    char fname[MAXSTRLEN];
    int ncid;
    int dimid;
    float* olon;
    float* olat;
    float* odepth;
    float* oval;
    float* oestd;
    float* ofij0;
    float* ofij1;
    float* ofij2;
    float* ofk;
    int* otype;
    int* oinst;
    float* otime;
    int vid_ids, vid_lcoeffs, vid_lon, vid_lat, vid_depth, vid_val, vid_estd, vid_fij0, vid_fij1, vid_fij2, vid_fk, vid_type, vid_inst, vid_time;
    char tunits[MAXSTRLEN];
    int otid, gid, iid;

    assert(das->s_mode == S_MODE_S_f);

    das_getfname_plog(das, plog, fname);
    ncw_create(fname, NC_CLOBBER | das->ncformat, &ncid);
    if (das->mode == MODE_ENOI)
        ncw_def_dim(ncid, "m", das->nmem, NULL);
    else if (das->mode == MODE_ENKF || das->mode == MODE_HYBRID) {
        ncw_def_dim(ncid, "m1", das->nmem_dynamic, NULL);
        ncw_def_dim(ncid, "m2", das->nmem, NULL);
    }
    ncw_def_dim(ncid, "p", ploc, &dimid);
    ncw_def_var(ncid, "obs_ids", NC_INT, 1, &dimid, &vid_ids);
    ncw_def_var(ncid, "lcoeffs", NC_FLOAT, 1, &dimid, &vid_lcoeffs);
    ncw_def_var(ncid, "lon", NC_FLOAT, 1, &dimid, &vid_lon);
    ncw_def_var(ncid, "lat", NC_FLOAT, 1, &dimid, &vid_lat);
    ncw_def_var(ncid, "depth", NC_FLOAT, 1, &dimid, &vid_depth);
    ncw_def_var(ncid, "obs_val", NC_FLOAT, 1, &dimid, &vid_val);
    ncw_def_var(ncid, "obs_estd", NC_FLOAT, 1, &dimid, &vid_estd);
    ncw_def_var(ncid, "obs_fij0", NC_FLOAT, 1, &dimid, &vid_fij0);
    ncw_def_var(ncid, "obs_fij1", NC_FLOAT, 1, &dimid, &vid_fij1);
    ncw_def_var(ncid, "obs_fij2", NC_FLOAT, 1, &dimid, &vid_fij2);
    ncw_def_var(ncid, "obs_fk", NC_FLOAT, 1, &dimid, &vid_fk);
    ncw_def_var(ncid, "obs_type", NC_INT, 1, &dimid, &vid_type);
    ncw_def_var(ncid, "obs_inst", NC_INT, 1, &dimid, &vid_inst);
    ncw_def_var(ncid, "obs_time", NC_FLOAT, 1, &dimid, &vid_time);
    snprintf(tunits, MAXSTRLEN, "days from %s", obs->datestr);
    ncw_put_att_text(ncid, vid_time, "units", tunits);

    /*
     * observation types
     */
    /*
     * first write type IDs
     */
    for (otid = 0; otid < obs->nobstypes; ++otid) {
        if (obs->obstypes[otid].statsonly)
            continue;

        ncw_put_att_int(ncid, vid_type, obs->obstypes[otid].name, 1, &otid);
    }
    /*
     * now write type related settings
     */
    for (otid = 0; otid < obs->nobstypes; ++otid) {
        char name[NC_MAX_NAME];
        char domains[MAXSTRLEN];
        int i;

        if (obs->obstypes[otid].statsonly)
            continue;

        snprintf(name, NC_MAX_NAME, "RFACTOR_%s", obs->obstypes[otid].name);
        ncw_put_att_double(ncid, vid_type, name, 1, &obs->obstypes[otid].rfactor);
        snprintf(name, NC_MAX_NAME, "LOCRAD_%s", obs->obstypes[otid].name);
        ncw_put_att_double(ncid, vid_type, name, obs->obstypes[otid].nlocrad, obs->obstypes[otid].locrad);
        snprintf(name, NC_MAX_NAME, "WEIGHT_%s", obs->obstypes[otid].name);
        ncw_put_att_double(ncid, vid_type, name, obs->obstypes[otid].nlocrad, obs->obstypes[otid].locweight);
        snprintf(name, NC_MAX_NAME, "GRIDID_%s", obs->obstypes[otid].name);
        ncw_put_att_int(ncid, vid_type, name, 1, &obs->obstypes[otid].gridid);
        domains[0] = 0;
        for (i = 0; i < obs->obstypes[otid].ndomains; ++i) {
            if (i > 0)
                strncat(domains, " ", NC_MAX_NAME);
            strncat(domains, obs->obstypes[otid].domainnames[i], NC_MAX_NAME);
        }
        if (i > 0) {
            snprintf(name, NC_MAX_NAME, "DOMAINS_%s", obs->obstypes[otid].name);
            ncw_put_att_text(ncid, vid_type, name, domains);
        }
    }

    /*
     * grids
     */
    for (gid = 0; gid < ngrid; ++gid) {
        grid* g = model_getgridbyid(das->m, gid);
        int aliasid = grid_getaliasid(g);
        char varname[NC_MAX_NAME];
        int vid_grid;
        int ni, nj, nk;
        float depth;
        char gridstr[SHORTSTRLEN - 5];

        if (plog->gridid >= 0 && plog->gridid != gid)
            continue;

        get_gridstr(das, gid, gridstr);
        snprintf(varname, SHORTSTRLEN - 1, "grid%s", gridstr);
        ncw_def_var(ncid, varname, NC_INT, 0, NULL, &vid_grid);
        ncw_put_att_int(ncid, vid_grid, "id", 1, &gid);
        if (aliasid >= 0)
            ncw_put_att_int(ncid, vid_grid, "hid", 1, &aliasid);
        ncw_put_att_text(ncid, vid_grid, "name", grid_getname(g));
        ncw_put_att_text(ncid, vid_grid, "domain", grid_getdomainname(g));
        if (grid_isstructured(g)) {
            ncw_put_att_double(ncid, vid_grid, "fi", 1, &plog->fij[gid][0]);
            ncw_put_att_double(ncid, vid_grid, "fj", 1, &plog->fij[gid][1]);
        } else {
            ncw_put_att_double(ncid, vid_grid, "fi0", 1, &plog->fij[gid][0]);
            ncw_put_att_double(ncid, vid_grid, "fi1", 1, &plog->fij[gid][1]);
            ncw_put_att_double(ncid, vid_grid, "fi2", 1, &plog->fij[gid][2]);
        }
        grid_getsize(g, &ni, &nj, &nk);
        ncw_put_att_int(ncid, vid_grid, "nk", 1, &nk);
        depth = grid_interpolate2d(g, plog->fij[gid], grid_getdepth(g));
        ncw_put_att_float(ncid, vid_grid, "model_depth", 1, &depth);
    }

    /*
     * instruments
     */
    for (iid = 0; iid < st_getsize(obs->instruments); ++iid)
        ncw_put_att_int(ncid, vid_inst, st_findstringbyindex(obs->instruments, iid), 1, &iid);

    /*
     * global atts
     */
    ncw_put_att_text(ncid, NC_GLOBAL, "version", ENKF_VERSION);
    ncw_put_att_double(ncid, NC_GLOBAL, "lon", 1, &plog->lon);
    ncw_put_att_double(ncid, NC_GLOBAL, "lat", 1, &plog->lat);
    if (das->mode == MODE_ENKF || das->mode == MODE_HYBRID) {
        ncw_put_att_text(ncid, NC_GLOBAL, "MODE", (das->mode == MODE_ENKF) ? "EnKF" : "Hybrid");
        ncw_put_att_text(ncid, NC_GLOBAL, "SCHEME", (das->scheme == SCHEME_DENKF) ? "DEnKF" : "ETKF");
        ncw_put_att_double(ncid, NC_GLOBAL, "ALPHA", 1, &das->alpha);
    } else
        ncw_put_att_text(ncid, NC_GLOBAL, "MODE", "EnOI");
    if (das->mode == MODE_HYBRID)
        ncw_put_att_double(ncid, NC_GLOBAL, "GAMMA", 1, &das->gamma);
    ncw_put_att_int(ncid, NC_GLOBAL, "ngrids", 1, &ngrid);

    if (das->nccompression > 0)
        ncw_def_deflate(ncid, 0, 1, das->nccompression);
    ncw_enddef(ncid);

    /*
     * obs
     */
    if (ploc > 0) {
        int oid;

        ncw_put_var_int(ncid, vid_ids, lobs);
        ncw_put_var_double(ncid, vid_lcoeffs, lcoeffs);

        olon = malloc(ploc * sizeof(float));
        olat = malloc(ploc * sizeof(float));
        odepth = malloc(ploc * sizeof(float));
        oval = malloc(ploc * sizeof(float));
        oestd = malloc(ploc * sizeof(float));
        ofij0 = malloc(ploc * sizeof(float));
        ofij1 = malloc(ploc * sizeof(float));
        ofij2 = malloc(ploc * sizeof(float));
        ofk = malloc(ploc * sizeof(float));
        otype = malloc(ploc * sizeof(int));
        oinst = malloc(ploc * sizeof(int));
        otime = malloc(ploc * sizeof(float));

        for (oid = 0; oid < ploc; ++oid) {
            observation* o = &obs->data[lobs[oid]];

            olon[oid] = o->lon;
            olat[oid] = o->lat;
            odepth[oid] = o->depth;
            oval[oid] = o->value;
            oestd[oid] = o->estd;
            ofij0[oid] = o->fij[0];
            ofij1[oid] = o->fij[1];
            ofij2[oid] = o->fij[2];
            ofk[oid] = o->fk;
            otype[oid] = o->type;
            oinst[oid] = o->instrument;
            otime[oid] = o->time;
        }

        ncw_put_var_float(ncid, vid_lon, olon);
        ncw_put_var_float(ncid, vid_lat, olat);
        ncw_put_var_float(ncid, vid_depth, odepth);
        ncw_put_var_float(ncid, vid_val, oval);
        ncw_put_var_float(ncid, vid_estd, oestd);
        ncw_put_var_float(ncid, vid_fij0, ofij0);
        ncw_put_var_float(ncid, vid_fij1, ofij1);
        ncw_put_var_float(ncid, vid_fij2, ofij2);
        ncw_put_var_float(ncid, vid_fk, ofk);
        ncw_put_var_int(ncid, vid_type, otype);
        ncw_put_var_int(ncid, vid_inst, oinst);
        ncw_put_var_float(ncid, vid_time, otime);

        free(olon);
        free(olat);
        free(odepth);
        free(oval);
        free(oestd);
        free(ofij0);
        free(ofij1);
        free(ofij2);
        free(ofk);
        free(otype);
        free(oinst);
        free(otime);
    }

    ncw_close(ncid);
}

/** Creates a point log file and writes ensemble observations, transforms, and
 * accompanying information.
 */
void das_writeplogtransform(dasystem* das, int plogid, int gid, int ploc, double* s, double* S, double* w, double* T)
{
    pointlog* plog = &das->plogs[plogid];
    grid* g = model_getgridbyid(das->m, gid);
    char* gridname = grid_getname(g);

    char fname[MAXSTRLEN];
    int ncid;
    int dimids[3];
    char name[NC_MAX_NAME];
    int vid_S, vid_s, vid_w, vid_T;
    char gridstr[SHORTSTRLEN];

    assert(das->s_mode == S_MODE_S_f);

    get_gridstr(das, gid, gridstr);

    das_getfname_plog(das, plog, fname);
    ncw_open(fname, NC_WRITE, &ncid);
    ncw_redef(ncid);
    if (das->mode == MODE_ENOI)
        ncw_inq_dimid(ncid, "m", &dimids[1]);
    else if (das->mode == MODE_ENKF || das->mode == MODE_HYBRID) {
        ncw_inq_dimid(ncid, "m1", &dimids[0]);
        ncw_inq_dimid(ncid, "m2", &dimids[1]);
    }

    if (ploc > 0) {
        snprintf(name, NC_MAX_NAME, "p%s", gridstr);
        if (!ncw_dim_exists(ncid, name))
            ncw_def_dim(ncid, name, ploc, &dimids[2]);
        else
            ncw_inq_dimid(ncid, name, &dimids[2]);
        snprintf(name, NC_MAX_NAME, "s%s", gridstr);
        ncw_def_var(ncid, name, NC_FLOAT, 1, &dimids[2], &vid_s);
        snprintf(name, NC_MAX_NAME, "S%s", gridstr);
        ncw_def_var(ncid, name, NC_FLOAT, 2, &dimids[1], &vid_S);
    }
    {
        char attstr[MAXSTRLEN];
        char varname[NC_MAX_NAME];

        snprintf(varname, NC_MAX_NAME, "w%s", gridstr);
        ncw_def_var(ncid, varname, NC_DOUBLE, 1, &dimids[1], &vid_w);
        if (grid_isstructured(g))
            snprintf(attstr, MAXSTRLEN, "ensemble coefficients for location (fi,fj)=(%.3f,%.3f) on grid %d (\"%s\")", plog->fij[gid][0], plog->fij[gid][1], gid, gridname);
        else
            snprintf(attstr, MAXSTRLEN, "ensemble coefficients for location (fi0,fi1,fi2)=(%.3f,%.3f,%.3f) on grid %d (\"%s\")", plog->fij[gid][0], plog->fij[gid][1], plog->fij[gid][2], gid, gridname);
        ncw_put_att_text(ncid, vid_w, "long_name", attstr);
    }
    if (T != NULL) {
        char attstr[MAXSTRLEN];

        snprintf(name, NC_MAX_NAME, "T%s", gridstr);
        ncw_def_var(ncid, name, NC_DOUBLE, 2, dimids, &vid_T);
        if (grid_isstructured(g))
            snprintf(attstr, MAXSTRLEN, "ensemble anomalies transform for location (fi,fj)=(%.3f,%.3f) on grid %d (\"%s\")", plog->fij[gid][0], plog->fij[gid][1], gid, gridname);
        else
            snprintf(attstr, MAXSTRLEN, "ensemble anomalies transform for location (fi0,fi1,fi2)=(%.3f,%.3f,%.3f) on grid %d (\"%s\")", plog->fij[gid][0], plog->fij[gid][1], plog->fij[gid][2], gid, gridname);
        ncw_put_att_text(ncid, vid_T, "long_name", attstr);
    }
    ncw_enddef(ncid);

    if (ploc > 0) {
        ncw_put_var_double(ncid, vid_s, s);
        ncw_put_var_double(ncid, vid_S, S);
        ncw_put_var_double(ncid, vid_w, w);
        if (T != NULL)
            ncw_put_var_double(ncid, vid_T, T);
    }

    ncw_close(ncid);
}
#endif

#if defined(ENKF_UPDATE)
/**
 */
void plogs_definestatevars(dasystem* das)
{
    int nvar = model_getnvar(das->m);
    int plogid;

    if (rank != 0)
        return;
    if (das->nplog == 0)
        return;

    for (plogid = 0; plogid < das->nplog; ++plogid) {
        pointlog* plog = &das->plogs[plogid];
        char fname[MAXSTRLEN];
        int ncid, vid;

        das_getfname_plog(das, plog, fname);
        ncw_open(fname, NC_WRITE, &ncid);
        ncw_redef(ncid);

        for (vid = 0; vid < nvar; ++vid) {
            int gid = model_getvargridid(das->m, vid);
            grid* g = model_getgridbyid(das->m, gid);
            char* varname = model_getvarname(das->m, vid);
            char gridstr[NC_MAX_NAME - 2];
            char varname_an[NC_MAX_NAME];
            int varid, varid_an;
            int nk;

            if (plog->gridid >= 0 && plog->gridid != gid)
                continue;

            get_gridstr(das, gid, gridstr);

            if (!(das->updatespec & UPDATE_OUTPUTINC))
                snprintf(varname_an, NC_MAX_NAME, "%s_an", varname);
            else
                snprintf(varname_an, NC_MAX_NAME, "%s_inc", varname);

            {
                char memfname[MAXSTRLEN];

                das_getmemberfname(das, varname, 1, memfname);
                nk = ncu_getnlevels(memfname, varname, grid_isstructured(g));
            }
            {
                char nkname[NC_MAX_NAME];
                int dimids[2];

                snprintf(nkname, NC_MAX_NAME, "nk%s", gridstr);
                if (nk > 1) {
                    if (!ncw_dim_exists(ncid, nkname))
                        ncw_def_dim(ncid, nkname, nk, &dimids[0]);
                    else
                        ncw_inq_dimid(ncid, nkname, &dimids[0]);
                }
                if (das->mode == MODE_ENOI)
                    ncw_inq_dimid(ncid, "m", &dimids[1]);
                else if (das->mode == MODE_ENKF || das->mode == MODE_HYBRID)
                    ncw_inq_dimid(ncid, "m2", &dimids[1]);
                if (!ncw_var_exists(ncid, varname)) {
                    if (nk > 1)
                        ncw_def_var(ncid, varname, NC_FLOAT, 2, dimids, &varid);
                    else
                        ncw_def_var(ncid, varname, NC_FLOAT, 1, &dimids[1], &varid);
                } else
                    ncw_inq_varid(ncid, varname, &varid);

                if (das->mode == MODE_ENKF || das->mode == MODE_HYBRID)
                    ncw_inq_dimid(ncid, "m1", &dimids[1]);
                if ((das->updatespec & UPDATE_DOPLOGSAN) || das->haveanalysis) {
                    if (!ncw_var_exists(ncid, varname_an)) {
                        if (nk > 1) {
                            if (das->mode == MODE_ENKF || das->mode == MODE_HYBRID)
                                ncw_def_var(ncid, varname_an, NC_FLOAT, 2, dimids, &varid_an);
                            else if (das->mode == MODE_ENOI)
                                ncw_def_var(ncid, varname_an, NC_FLOAT, 1, dimids, &varid_an);
                        } else {
                            if (das->mode == MODE_ENKF || das->mode == MODE_HYBRID)
                                ncw_def_var(ncid, varname_an, NC_FLOAT, 1, &dimids[1], &varid_an);
                            else if (das->mode == MODE_ENOI)
                                ncw_def_var(ncid, varname_an, NC_FLOAT, 0, NULL, &varid_an);
                        }
                    } else
                        ncw_inq_varid(ncid, varname_an, &varid_an);
                }

                if (das->mode == MODE_ENOI) {
                    char varname_bg[NC_MAX_NAME];
                    int varid_bg;

                    snprintf(varname_bg, NC_MAX_NAME, "%s_bg", varname);
                    if (!ncw_var_exists(ncid, varname_bg))
                        ncw_def_var(ncid, varname_bg, NC_FLOAT, (nk > 1) ? 1 : 0, (nk > 1) ? dimids : NULL, &varid_bg);
                    else
                        ncw_inq_varid(ncid, varname_bg, &varid_bg);
                    ncw_put_att_int(ncid, varid_bg, "gridid", 1, &gid);
                }
            }
            ncw_put_att_int(ncid, varid, "gridid", 1, &gid);
            ncw_put_att_int(ncid, varid_an, "gridid", 1, &gid);

            if (das->mode == MODE_ENKF || das->mode == MODE_HYBRID) {
                float inflation[2];
                double tmp;

                model_getvarinflation(das->m, vid, &inflation[0], &tmp);
                inflation[1] = (float) tmp;
                ncw_put_att_float(ncid, varid_an, "INFLATION", 2, inflation);
            }
        }
        /*
         * (putting this attribute should have been done in das_createplog(),
         * but it is called in CALC, which knows nothing about das->updatespec)
         */
        ncw_put_att_text(ncid, NC_GLOBAL, "output", (das->updatespec & UPDATE_OUTPUTINC) ? "increment" : "analysis");
        ncw_close(ncid);
        enkf_writeinfo(fname);
    }
}

static void plogs_writeens(dasystem* das, int isanalysis)
{
    model* m = das->m;
    int nvar = model_getnvar(m);
    int nmem = das->nmem;
    int vid, plogid;

    if (das->nplog == 0)
        return;
    if (isanalysis && das->mode == MODE_ENOI)
        return;

    if (isanalysis && (das->mode == MODE_ENKF || das->mode == MODE_HYBRID))
        nmem = das->nmem_dynamic;

    distribute_iterations(0, nmem - 1, nprocesses, "    ");
    for (vid = 0; vid < nvar; ++vid) {
        char* varname = model_getvarname(m, vid);
        grid* g = model_getvargrid(m, vid);
        int gid = grid_getid(g);
        float*** dst = NULL;
        float** src = NULL;
        int ni, nj, nk, e;

        enkf_printf("    %s:", varname);

        model_getvargridsize(m, vid, &ni, &nj, NULL);
        {
            char memfname[MAXSTRLEN];

            das_getmemberfname(das, varname, 1, memfname);
            nk = ncu_getnlevels(memfname, varname, nj > 0);
        }
        if (my_number_of_iterations > 0) {
            dst = alloc3d(das->nplog, (rank == 0) ? nmem : my_number_of_iterations, nk, sizeof(float));
            src = alloc2d((nj > 0) ? nj : 1, ni, sizeof(float));
        }

        for (e = my_first_iteration; e <= my_last_iteration; ++e) {
            char fname_src[MAXSTRLEN];
            int k;

            das_getmemberfname(das, varname, e + 1, fname_src);
            if (isanalysis) {
                if ((das->updatespec & UPDATE_OUTPUTINC))
                    strncat(fname_src, ".increment", MAXSTRLEN - 1);
                else
                    strncat(fname_src, ".analysis", MAXSTRLEN - 1);
            }

            for (k = 0; k < nk; ++k) {
                model_readfield(m, fname_src, varname, 0, (nk > 1) ? k : -1, src[0], 1);
                for (plogid = 0; plogid < das->nplog; ++plogid) {
                    pointlog* plog = &das->plogs[plogid];

                    if (plog->gridid >= 0 && plog->gridid != gid)
                        continue;
                    if (isnan(plog->fij[gid][0] + plog->fij[gid][1]))
                        continue;

                    dst[plogid][e - my_first_iteration][k] = grid_interpolate2d(g, plog->fij[gid], src);
                }
                enkf_printf(".");
                enkf_flush();
            }
        }

#if defined(MPI)
        if (nprocesses > 1) {
            int* recvcounts = NULL;
            int* displs = NULL;
            int p;

            if (rank == 0) {
                recvcounts = malloc(nprocesses * sizeof(int));
                displs = malloc(nprocesses * sizeof(int));
                for (p = 0; p < nprocesses; ++p)
                    recvcounts[p] = number_of_iterations[p] * nk;
                displs[0] = 0;
                for (p = 1; p < nprocesses; ++p)
                    displs[p] = displs[p - 1] + recvcounts[p];
            }

            for (plogid = 0; plogid < das->nplog; ++plogid) {
                MPI_Gatherv((my_number_of_iterations > 0) ? dst[plogid][0] : NULL, my_number_of_iterations * nk, MPI_FLOAT, (my_number_of_iterations > 0) ? dst[plogid][0] : NULL, recvcounts, displs, MPI_FLOAT, 0, MPI_COMM_WORLD);
            }

            if (rank == 0) {
                free(displs);
                free(recvcounts);
            }
        }
        if (rank == 0) {
            float** dst2 = alloc2d(nk, nmem, sizeof(float));

            for (plogid = 0; plogid < das->nplog; ++plogid) {
                pointlog* plog = &das->plogs[plogid];
                float** dst1 = dst[plogid];
                char fname_dst[MAXSTRLEN];
                int ncid_dst, vid_dst;
                int k;

                for (k = 0; k < nk; ++k)
                    for (e = 0; e < nmem; ++e)
                        dst2[k][e] = dst1[e][k];

                das_getfname_plog(das, plog, fname_dst);
                ncw_open(fname_dst, NC_WRITE, &ncid_dst);
                if (!isanalysis)
                    ncw_inq_varid(ncid_dst, varname, &vid_dst);
                else {
                    char varname_dst[MAXSTRLEN];

                    strncpy(varname_dst, varname, MAXSTRLEN - 1);
                    strncat(varname_dst, !(das->updatespec & UPDATE_OUTPUTINC) ? "_an" : "_inc", NC_MAX_NAME - 5);
                    ncw_inq_varid(ncid_dst, varname_dst, &vid_dst);
                }
                ncw_put_var_float(ncid_dst, vid_dst, dst2[0]);
                ncw_close(ncid_dst);
            }

            free(dst2);
        }
#endif
        if (src != NULL)
            free(src);
        if (dst != NULL)
            free(dst);
        enkf_printf("\n");
        enkf_flush();
    }
}

static void plogs_writebg(dasystem* das, int isanalysis)
{
    model* m = das->m;
    int nvar = model_getnvar(m);
    int vid, plogid;

    if (das->nplog == 0)
        return;
    if (rank != 0 || das->mode != MODE_ENOI)
        return;

    for (vid = 0; vid < nvar; ++vid) {
        char* varname = model_getvarname(m, vid);
        grid* g = model_getvargrid(m, vid);
        int gid = grid_getid(g);
        float** dst = NULL;
        float** src = NULL;
        char fname_src[MAXSTRLEN];
        char varname_dst[MAXSTRLEN];
        char fname_dst[MAXSTRLEN];
        int ni, nj, nk, k;

        enkf_printf("    %s:", varname);

        model_getvargridsize(m, vid, &ni, &nj, NULL);
        das_getbgfname(das, varname, fname_src);
        nk = ncu_getnlevels(fname_src, varname, nj > 0);

        dst = alloc2d(das->nplog, nk, sizeof(float));
        src = alloc2d((nj > 0) ? nj : 1, ni, sizeof(float));

        if (!isanalysis)
            snprintf(varname_dst, MAXSTRLEN, "%s_bg", varname);
        else {
            snprintf(varname_dst, MAXSTRLEN, "%s_%s", varname, (das->updatespec & UPDATE_OUTPUTINC) ? "inc" : "an");
            if ((das->updatespec & UPDATE_OUTPUTINC))
                strncat(fname_src, ".increment", MAXSTRLEN - 1);
            else
                strncat(fname_src, ".analysis", MAXSTRLEN - 1);
        }

        for (k = 0; k < nk; ++k) {
            model_readfield(m, fname_src, varname, 0, (nk > 1) ? k : -1, src[0], 1);
            for (plogid = 0; plogid < das->nplog; ++plogid) {
                pointlog* plog = &das->plogs[plogid];

                if (plog->gridid >= 0 && plog->gridid != gid)
                    continue;
                if (isnan(plog->fij[gid][0] + plog->fij[gid][1]))
                    continue;

                dst[plogid][k] = grid_interpolate2d(g, plog->fij[gid], src);
            }
            enkf_printf(".");
            enkf_flush();
        }

        for (plogid = 0; plogid < das->nplog; ++plogid) {
            pointlog* plog = &das->plogs[plogid];
            int ncid_dst, vid_dst;

            das_getfname_plog(das, plog, fname_dst);
            ncw_open(fname_dst, NC_WRITE, &ncid_dst);
            ncw_inq_varid(ncid_dst, varname_dst, &vid_dst);
            ncw_put_var_float(ncid_dst, vid_dst, dst[plogid]);
            ncw_close(ncid_dst);
        }

        free(src);
        free(dst);
        enkf_printf("\n");
        enkf_flush();
    }
}

void plogs_writestatevars(dasystem* das, int isanalysis)
{
    plogs_writeens(das, isanalysis);
    if (das->mode == MODE_ENOI && rank == 0)
        plogs_writebg(das, isanalysis);
}
#endif                          /* ENKF_UPDATE */
