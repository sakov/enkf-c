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
    plog->fi = NULL;
    plog->fj = NULL;
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
        if (plog->fi != NULL) {
            free(plog->fi);
            free(plog->fj);
        }
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
void plog_create(dasystem* das, int plogid, int ploc, int* lobs, double* lcoeffs)
{
    pointlog* plog = &das->plogs[plogid];
    observations* obs = das->obs;
    int ngrid = model_getngrid(das->m);
    char fname[MAXSTRLEN];
    int ncid;
    int dimids[2];
    float* olon;
    float* olat;
    float* odepth;
    float* oval;
    float* oestd;
    float* ofi;
    float* ofj;
    float* ofk;
    int* otype;
    int* oinst;
    float* otime;
    int vid_ids, vid_lcoeffs, vid_lon, vid_lat, vid_depth, vid_val, vid_estd, vid_fi, vid_fj, vid_fk, vid_type, vid_inst, vid_time;
    char tunits[MAXSTRLEN];
    int otid, gid, iid;

    assert(das->s_mode == S_MODE_S_f);

    das_getfname_plog(das, plog, fname);
    ncw_create(fname, NC_CLOBBER | das->ncformat, &ncid);
    ncw_def_dim(ncid, "m", das->nmem, &dimids[0]);
    ncw_def_dim(ncid, "p", ploc, &dimids[1]);
    ncw_def_var(ncid, "obs_ids", NC_INT, 1, &dimids[1], &vid_ids);
    ncw_def_var(ncid, "lcoeffs", NC_FLOAT, 1, &dimids[1], &vid_lcoeffs);
    ncw_def_var(ncid, "lon", NC_FLOAT, 1, &dimids[1], &vid_lon);
    ncw_def_var(ncid, "lat", NC_FLOAT, 1, &dimids[1], &vid_lat);
    ncw_def_var(ncid, "depth", NC_FLOAT, 1, &dimids[1], &vid_depth);
    ncw_def_var(ncid, "obs_val", NC_FLOAT, 1, &dimids[1], &vid_val);
    ncw_def_var(ncid, "obs_estd", NC_FLOAT, 1, &dimids[1], &vid_estd);
    ncw_def_var(ncid, "obs_fi", NC_FLOAT, 1, &dimids[1], &vid_fi);
    ncw_def_var(ncid, "obs_fj", NC_FLOAT, 1, &dimids[1], &vid_fj);
    ncw_def_var(ncid, "obs_fk", NC_FLOAT, 1, &dimids[1], &vid_fk);
    ncw_def_var(ncid, "obs_type", NC_INT, 1, &dimids[1], &vid_type);
    ncw_def_var(ncid, "obs_inst", NC_INT, 1, &dimids[1], &vid_inst);
    ncw_def_var(ncid, "obs_time", NC_FLOAT, 1, &dimids[1], &vid_time);
    snprintf(tunits, MAXSTRLEN, "days from %s", obs->datestr);
    ncw_put_att_text(ncid, vid_time, "units", tunits);

    /*
     * observation types
     */
    for (otid = 0; otid < obs->nobstypes; ++otid) {
        char name[NC_MAX_NAME];
        char domains[MAXSTRLEN];
        int i;

        if (obs->obstypes[otid].statsonly)
            continue;

        ncw_put_att_int(ncid, vid_type, obs->obstypes[otid].name, 1, &otid);
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
        char varname[NC_MAX_NAME];
        int vid_grid;
        int ni, nj, nk;
        float** depths;
        float depth = NAN;
        char gridstr[MAXSTRLEN];

        if (plog->gridid >= 0 && plog->gridid != gid)
            continue;

        get_gridstr(das, gid, gridstr);
        snprintf(varname, NC_MAX_NAME, "grid%s", gridstr);
        ncw_def_var(ncid, varname, NC_INT, 0, NULL, &vid_grid);
        ncw_put_att_int(ncid, vid_grid, "id", 1, &gid);
        ncw_put_att_text(ncid, vid_grid, "NAME", grid_getname(g));
        ncw_put_att_text(ncid, vid_grid, "DOMAIN", grid_getdomainname(g));
        ncw_put_att_double(ncid, vid_grid, "fi", 1, &plog->fi[gid]);
        ncw_put_att_double(ncid, vid_grid, "fj", 1, &plog->fj[gid]);
        grid_getsize(g, &ni, &nj, &nk);
        ncw_put_att_int(ncid, vid_grid, "nk", 1, &nk);
        depths = grid_getdepth(g);
        if (depths != NULL) {
            depth = interpolate2d(plog->fi[gid], plog->fj[gid], ni, nj, depths, grid_getnumlevels(g), grid_isperiodic_i(g));
            ncw_put_att_float(ncid, vid_grid, "model_depth", 1, &depth);
        }
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
    ncw_put_att_text(ncid, NC_GLOBAL, "date", obs->datestr);
    ncw_put_att_double(ncid, NC_GLOBAL, "lon", 1, &plog->lon);
    ncw_put_att_double(ncid, NC_GLOBAL, "lat", 1, &plog->lat);
    if (das->mode == MODE_ENKF) {
        ncw_put_att_text(ncid, NC_GLOBAL, "MODE", "EnKF");
        ncw_put_att_text(ncid, NC_GLOBAL, "SCHEME", (das->scheme == SCHEME_DENKF) ? "DEnKF" : "ETKF");
        ncw_put_att_double(ncid, NC_GLOBAL, "ALPHA", 1, &das->alpha);
    } else
        ncw_put_att_text(ncid, NC_GLOBAL, "MODE", "EnOI");
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
        ofi = malloc(ploc * sizeof(float));
        ofj = malloc(ploc * sizeof(float));
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
            ofi[oid] = o->fi;
            ofj[oid] = o->fj;
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
        ncw_put_var_float(ncid, vid_fi, ofi);
        ncw_put_var_float(ncid, vid_fj, ofj);
        ncw_put_var_float(ncid, vid_fk, ofk);
        ncw_put_var_int(ncid, vid_type, otype);
        ncw_put_var_int(ncid, vid_inst, oinst);
        ncw_put_var_float(ncid, vid_time, otime);

        free(olon);
        free(olat);
        free(odepth);
        free(oval);
        free(oestd);
        free(ofi);
        free(ofj);
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
void plog_writetransform(dasystem* das, int plogid, int gid, int ploc, double* s, double* S, double* transform)
{
    pointlog* plog = &das->plogs[plogid];
    char* gridname = grid_getname(model_getgridbyid(das->m, gid));

    char fname[MAXSTRLEN];
    int ncid;
    int dimids[2];
    char name[NC_MAX_NAME];
    int vid_S, vid_s, vid_transform;
    char gridstr[MAXSTRLEN];

    assert(das->s_mode == S_MODE_S_f);

    get_gridstr(das, gid, gridstr);

    das_getfname_plog(das, plog, fname);
    ncw_open(fname, NC_WRITE, &ncid);
    ncw_redef(ncid);
    ncw_inq_dimid(ncid, "m", &dimids[0]);
    if (ploc > 0) {
        snprintf(name, NC_MAX_NAME, "p%s", gridstr);
        if (!ncw_dim_exists(ncid, name))
            ncw_def_dim(ncid, name, ploc, &dimids[1]);
        else
            ncw_inq_dimid(ncid, name, &dimids[1]);
        snprintf(name, NC_MAX_NAME, "s%s", gridstr);
        ncw_def_var(ncid, name, NC_FLOAT, 1, &dimids[1], &vid_s);
        snprintf(name, NC_MAX_NAME, "S%s", gridstr);
        ncw_def_var(ncid, name, NC_FLOAT, 2, dimids, &vid_S);
    }
    if (das->mode == MODE_ENKF) {
        char attstr[MAXSTRLEN];

        dimids[1] = dimids[0];
        snprintf(name, NC_MAX_NAME, "X5%s", gridstr);
        ncw_def_var(ncid, name, NC_DOUBLE, 2, dimids, &vid_transform);
        snprintf(attstr, MAXSTRLEN, "ensemble transform calculated for location (fi,fj)=(%.3f,%.3f) on grid %d (\"%s\")", plog->fi[gid], plog->fj[gid], gid, gridname);
        ncw_put_att_text(ncid, vid_transform, "long_name", attstr);

    } else if (das->mode == MODE_ENOI) {
        char attstr[MAXSTRLEN];
        char varname[NC_MAX_NAME];

        snprintf(varname, NC_MAX_NAME, "w%s", gridstr);
        ncw_def_var(ncid, varname, NC_DOUBLE, 1, &dimids[0], &vid_transform);
        snprintf(attstr, MAXSTRLEN, "ensemble coefficients calculated for location (fi,fj)=(%.3f,%.3f) on grid %d (\"%s\")", plog->fi[gid], plog->fj[gid], gid, gridname);
        ncw_put_att_text(ncid, vid_transform, "long_name", attstr);
    }
    ncw_enddef(ncid);

    if (ploc > 0) {
        ncw_put_var_double(ncid, vid_s, s);
        ncw_put_var_double(ncid, vid_S, S);
        ncw_put_var_double(ncid, vid_transform, transform);
    }

    ncw_close(ncid);
}
#endif

#if defined(ENKF_UPDATE)
/**
 */
void plog_definestatevars(dasystem* das)
{
    int nvar = model_getnvar(das->m);
    int vid, plogid;

    if (das->nplog == 0)
        return;

    for (plogid = 0; plogid < das->nplog; ++plogid) {
        pointlog* plog = &das->plogs[plogid];
        char fname[MAXSTRLEN];
        int ncid;

        das_getfname_plog(das, plog, fname);
        ncw_open(fname, NC_WRITE, &ncid);
        ncw_redef(ncid);

        for (vid = 0; vid < nvar; ++vid) {
            int gid = model_getvargridid(das->m, vid);
            char* varname = model_getvarname(das->m, vid);
            char gridstr[MAXSTRLEN];
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
                char fname[MAXSTRLEN];

                model_getmemberfname(das->m, das->ensdir, varname, 1, fname);
                nk = getnlevels(fname, varname);
            }
            if (nk > 1) {
                char gridstr[NC_MAX_NAME];
                char nkname[NC_MAX_NAME];
                int dimids[2];

                get_gridstr(das, gid, gridstr);
                snprintf(nkname, NC_MAX_NAME, "nk%s", gridstr);
                if (!ncw_dim_exists(ncid, nkname))
                    ncw_def_dim(ncid, nkname, nk, &dimids[0]);
                else
                    ncw_inq_dimid(ncid, nkname, &dimids[0]);
                ncw_inq_dimid(ncid, "m", &dimids[1]);

                if (!ncw_var_exists(ncid, varname)) {
                    ncw_def_var(ncid, varname, NC_FLOAT, 2, dimids, &varid);
                    ncw_def_var(ncid, varname_an, NC_FLOAT, 2, dimids, &varid_an);
                } else {
                    ncw_inq_varid(ncid, varname, &varid);
                    ncw_inq_varid(ncid, varname_an, &varid_an);
                }
            } else {
                int dimid;

                ncw_inq_dimid(ncid, "m", &dimid);
                if (!ncw_var_exists(ncid, varname)) {
                    ncw_def_var(ncid, varname, NC_FLOAT, 1, &dimid, &varid);
                    ncw_def_var(ncid, varname_an, NC_FLOAT, 1, &dimid, &varid_an);
                } else {
                    ncw_inq_varid(ncid, varname, &varid);
                    ncw_inq_varid(ncid, varname_an, &varid_an);
                }
            }
            ncw_put_att_int(ncid, varid, "gridid", 1, &gid);

            if (das->mode == MODE_ENKF) {
                float inflation[2];
                double tmp;

                model_getvarinflation(das->m, vid, &inflation[0], &tmp);
                inflation[1] = (float) tmp;
                ncw_put_att_float(ncid, varid, "INFLATION", 2, inflation);
            }
        }
        /*
         * (putting this attribute should have been done in plog_create(),
         * but it is called in CALC, which knows nothing about das->updatespec)
         */
        ncw_put_att_text(ncid, NC_GLOBAL, "output", (das->updatespec & UPDATE_OUTPUTINC) ? "increment" : "analysis");
        ncw_close(ncid);
    }
}

/** Writes state variables directly to the pointlogs (without tiles).
 */
static void plog_writestatevars_direct(dasystem* das, int nfields, void** fieldbuffer, field* fields, int isanalysis)
{
    int p, fid, e;
    float*** v_src = NULL;
    float* v = NULL;
    size_t start[2] = { 0, 0 };
    size_t count[2] = { 1, das->nmem };

    v = malloc(das->nmem * sizeof(float));

    for (p = 0; p < das->nplog; ++p) {
        pointlog* plog = &das->plogs[p];
        char fname[MAXSTRLEN];
        int ncid;

        das_getfname_plog(das, plog, fname);
        ncw_open(fname, NC_WRITE, &ncid);

        for (fid = 0; fid < nfields; ++fid) {
            field* f = &fields[fid];
            grid* g = model_getvargrid(das->m, f->varid);
            int gid = grid_getid(g);
            int** mask = model_getnumlevels(das->m, f->varid);
            int periodic_i = grid_isperiodic_i(g);
            char varname[NC_MAX_NAME];
            int vid, ndims;
            int ni, nj;

            if (plog->gridid >= 0 && plog->gridid != gid)
                continue;

            v_src = (float***) fieldbuffer[fid];
            grid_getsize(g, &ni, &nj, NULL);

            if (das->mode == MODE_ENKF) {
                for (e = 0; e < das->nmem; ++e)
                    v[e] = interpolate2d(plog->fi[gid], plog->fj[gid], ni, nj, v_src[e], mask, periodic_i);
            } else if (das->mode == MODE_ENOI) {
                float bg = interpolate2d(plog->fi[gid], plog->fj[gid], ni, nj, v_src[das->nmem], mask, periodic_i);

                if (isanalysis && (das->updatespec & UPDATE_OUTPUTINC)) {
                    /*
                     * same increments for all ensemble members
                     */
                    for (e = 0; e < das->nmem; ++e)
                        v[e] = bg;
                } else {
                    for (e = 0; e < das->nmem; ++e)
                        v[e] = bg + interpolate2d(plog->fi[gid], plog->fj[gid], ni, nj, v_src[e], mask, periodic_i);
                }
            }

            snprintf(varname, NC_MAX_NAME, "%s", f->varname);
            if (isanalysis)
                strncat(varname, !(das->updatespec & UPDATE_OUTPUTINC) ? "_an" : "_inc", NC_MAX_NAME);
            ncw_inq_varid(ncid, varname, &vid);
            ncw_inq_varndims(ncid, vid, &ndims);
            if (ndims == 1)
                ncw_put_var_float(ncid, vid, v);
            else if (ndims == 2) {
                start[0] = f->level;
                ncw_put_vara_float(ncid, vid, start, count, v);
            }
        }

        ncw_close(ncid);
    }

    free(v);
}

/** Writes state variables to tiles that are to be assembled into the pointlogs
 *  later.
 */
static void plog_writestatevars_toassemble(dasystem* das, int nfields, void** fieldbuffer, field* fields, int isanalysis)
{
    float* v = malloc(das->nmem * sizeof(float));
    int vid;
    float*** v_src = NULL;
    int fid;

    for (fid = 0; fid < nfields; ++fid) {
        field* f = &fields[fid];
        grid* g = model_getvargrid(das->m, f->varid);
        int gid = grid_getid(g);
        int** mask = model_getnumlevels(das->m, f->varid);
        int periodic_i = grid_isperiodic_i(g);
        char fname[MAXSTRLEN];
        int ncid, dimid;
        int plogid;
        int ni, nj;
        int e;

        v_src = (float***) fieldbuffer[fid];
        grid_getsize(g, &ni, &nj, NULL);

        for (plogid = 0; plogid < das->nplog; ++plogid) {
            pointlog* plog = &das->plogs[plogid];
            char varname[NC_MAX_NAME];

            if (plog->gridid >= 0 && plog->gridid != gid)
                continue;

            snprintf(fname, MAXSTRLEN, "%s/pointlog-%d_%s-%03d.nc", DIRNAME_TMP, plogid, f->varname, f->level);
            if (!isanalysis) {
                ncw_create(fname, NC_CLOBBER | das->ncformat, &ncid);
                if (das->nccompression > 0)
                    ncw_def_deflate(ncid, 0, 1, das->nccompression);
                ncw_def_dim(ncid, "m", das->nmem, &dimid);
            } else {
                ncw_open(fname, NC_WRITE, &ncid);
                ncw_inq_dimid(ncid, "m", &dimid);
                ncw_redef(ncid);
            }

            snprintf(varname, NC_MAX_NAME, "%s", f->varname);
            if (isanalysis)
                strncat(varname, !(das->updatespec & UPDATE_OUTPUTINC) ? "_an" : "_inc", NC_MAX_NAME);
            ncw_def_var(ncid, varname, NC_FLOAT, 1, &dimid, &vid);
            ncw_enddef(ncid);

            if (das->mode == MODE_ENKF) {
                for (e = 0; e < das->nmem; ++e)
                    v[e] = interpolate2d(plog->fi[gid], plog->fj[gid], ni, nj, v_src[e], mask, periodic_i);
            } else if (das->mode == MODE_ENOI) {
                float bg = interpolate2d(plog->fi[gid], plog->fj[gid], ni, nj, v_src[das->nmem], mask, periodic_i);

                if (isanalysis && (das->updatespec & UPDATE_OUTPUTINC)) {
                    /*
                     * same increments for all ensemble members
                     */
                    for (e = 0; e < das->nmem; ++e)
                        v[e] = bg;
                } else {
                    for (e = 0; e < das->nmem; ++e)
                        v[e] = bg + interpolate2d(plog->fi[gid], plog->fj[gid], ni, nj, v_src[e], mask, periodic_i);
                }
            }

            ncw_put_var_float(ncid, vid, v);
            ncw_close(ncid);
        }
    }

    free(v);
}

/**
 */
void plog_writestatevars(dasystem* das, int nfields, void** fieldbuffer, field* fields, int isanalysis)
{
    if (das->nplog == 0)
        return;

    if (das->updatespec & UPDATE_DIRECTWRITE)
        plog_writestatevars_direct(das, nfields, fieldbuffer, fields, isanalysis);
    else
        plog_writestatevars_toassemble(das, nfields, fieldbuffer, fields, isanalysis);
}

/**
 */
void plog_assemblestatevars(dasystem* das)
{
    float* v = NULL;
    int nfields = 0;
    field* fields = NULL;
    int plogid, fid;

    if (das->nplog == 0)
        return;

    v = malloc(das->nmem * sizeof(float));

    das_getfields(das, -1, &nfields, &fields);

#if defined(MPI)
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    distribute_iterations(0, das->nplog - 1, nprocesses, rank, "    ");
    for (plogid = my_first_iteration; plogid <= my_last_iteration; ++plogid) {
        pointlog* plog = &das->plogs[plogid];
        char fname_dst[MAXSTRLEN];
        int ncid_dst;
        size_t start[2] = { 0, 0 };
        size_t count[2] = { 1, das->nmem };

        das_getfname_plog(das, plog, fname_dst);
        ncw_open(fname_dst, NC_WRITE, &ncid_dst);

        for (fid = 0; fid < nfields; ++fid) {
            field* f = &fields[fid];
            int gid = model_getvargridid(das->m, f->varid);
            char fname_src[MAXSTRLEN];
            int ii;

            if (plog->gridid >= 0 && plog->gridid != gid)
                continue;

            snprintf(fname_src, MAXSTRLEN, "%s/pointlog-%d_%s-%03d.nc", DIRNAME_TMP, plogid, f->varname, f->level);

            for (ii = 0; ii < 2; ++ii) {
                char varname[NC_MAX_NAME];
                int ncid_src, vid_src, vid_dst, ndims_dst;

                snprintf(varname, NC_MAX_NAME, "%s", f->varname);
                if (ii == 1)
                    strncat(varname, !(das->updatespec & UPDATE_OUTPUTINC) ? "_an" : "_inc", NC_MAX_NAME);

                ncw_open(fname_src, NC_NOWRITE, &ncid_src);
                ncw_inq_varid(ncid_src, varname, &vid_src);
                ncw_get_var_float(ncid_src, vid_src, v);
                ncw_close(ncid_src);

                ncw_inq_varid(ncid_dst, varname, &vid_dst);
                ncw_inq_varndims(ncid_dst, vid_dst, &ndims_dst);
                if (ndims_dst == 1)
                    ncw_put_var_float(ncid_dst, vid_dst, v);
                else if (ndims_dst == 2) {
                    start[0] = f->level;
                    ncw_put_vara_float(ncid_dst, vid_dst, start, count, v);
                }
            }
        }
        ncw_close(ncid_dst);
    }

#if defined(MPI)
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    enkf_printf("    deleting tiles:\n");
    distribute_iterations(0, nfields - 1, nprocesses, rank, "      ");
    for (fid = my_first_iteration; fid <= my_last_iteration; ++fid) {
        field* f = &fields[fid];
        char fname[MAXSTRLEN];
        int plogid;

        for (plogid = 0; plogid < das->nplog; ++plogid) {
            snprintf(fname, MAXSTRLEN, "%s/pointlog-%d_%s-%03d.nc", DIRNAME_TMP, plogid, f->varname, f->level);
            file_delete(fname);
        }
    }

    free(v);
    free(fields);
}
#endif
