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
#include "definitions.h"
#include "distribute.h"
#include "utils.h"
#include "dasystem.h"
#include "pointlog.h"

/** Creates a point log file and writes ensemble observations, transforms, and
 * accompanying information.
 */
void plog_write(dasystem* das, int id, double depth, int p, int* lobs, double* lcoeffs, double* s, double* S, double* transform)
{
    pointlog* plog = &das->plogs[id];
    observations* obs = das->obs;

    char fname[MAXSTRLEN];
    int ncid;
    int dimids[2];
    int vid_ids, vid_S, vid_s, vid_lcoeffs, vid_transform, vid_transform_actual, vid_lon, vid_lat, vid_depth, vid_val, vid_std, vid_fi, vid_fj, vid_fk, vid_type, vid_date;
    char tunits[MAXSTRLEN];

    float* olon;
    float* olat;
    float* odepth;
    float* val;
    float* std;
    float* fi;
    float* fj;
    float* fk;
    int* otype;
    float* odate;

    int ot, oid;

    if (das->nplogs == 0)
        return;

    assert(das->s_mode == S_MODE_S_f);

    das_getfname_plog(das, plog, fname);
    ncw_create(fname, NC_CLOBBER | NETCDF_FORMAT, &ncid);
    ncw_def_dim(fname, ncid, "m", das->nmem, &dimids[0]);
    ncw_def_dim(fname, ncid, "p", p, &dimids[1]);
    if (p > 0) {
        ncw_def_var(fname, ncid, "obs_ids", NC_INT, 1, &dimids[1], &vid_ids);
        ncw_def_var(fname, ncid, "lcoeffs", NC_FLOAT, 1, &dimids[1], &vid_lcoeffs);
        ncw_def_var(fname, ncid, "lon", NC_FLOAT, 1, &dimids[1], &vid_lon);
        ncw_def_var(fname, ncid, "lat", NC_FLOAT, 1, &dimids[1], &vid_lat);
        ncw_def_var(fname, ncid, "depth", NC_FLOAT, 1, &dimids[1], &vid_depth);
        ncw_def_var(fname, ncid, "obs_val", NC_FLOAT, 1, &dimids[1], &vid_val);
        ncw_def_var(fname, ncid, "obs_std", NC_FLOAT, 1, &dimids[1], &vid_std);
        ncw_def_var(fname, ncid, "obs_fi", NC_FLOAT, 1, &dimids[1], &vid_fi);
        ncw_def_var(fname, ncid, "obs_fj", NC_FLOAT, 1, &dimids[1], &vid_fj);
        ncw_def_var(fname, ncid, "obs_fk", NC_FLOAT, 1, &dimids[1], &vid_fk);
        ncw_def_var(fname, ncid, "obs_type", NC_INT, 1, &dimids[1], &vid_type);
        ncw_def_var(fname, ncid, "obs_date", NC_FLOAT, 1, &dimids[1], &vid_date);
        ncw_def_var(fname, ncid, "s", NC_FLOAT, 1, &dimids[1], &vid_s);
        ncw_def_var(fname, ncid, "S", NC_FLOAT, 2, dimids, &vid_S);
        snprintf(tunits, MAXSTRLEN, "days from %s", obs->datestr);
        ncw_put_att_text(fname, ncid, vid_date, "units", tunits);
        for (ot = 0; ot < obs->nobstypes; ++ot) {
            char attname[NC_MAX_NAME];

            ncw_put_att_int(fname, ncid, vid_type, obs->obstypes[ot].name, 1, &ot);
            snprintf(attname, NC_MAX_NAME, "RFACTOR_%s", obs->obstypes[ot].name);
            ncw_put_att_double(fname, ncid, vid_type, attname, 1, &obs->obstypes[ot].rfactor);
            snprintf(attname, NC_MAX_NAME, "LOCRAD_%s", obs->obstypes[ot].name);
            ncw_put_att_double(fname, ncid, vid_type, attname, obs->obstypes[ot].nscale, obs->obstypes[ot].locrad);
            snprintf(attname, NC_MAX_NAME, "WEIGHT_%s", obs->obstypes[ot].name);
            ncw_put_att_double(fname, ncid, vid_type, attname, obs->obstypes[ot].nscale, obs->obstypes[ot].weight);
        }
    }
    if (das->mode == MODE_ENKF) {
        char attstr[MAXSTRLEN];

        dimids[1] = dimids[0];
        ncw_def_var(fname, ncid, "X5", NC_DOUBLE, 2, dimids, &vid_transform);
        snprintf(attstr, MAXSTRLEN, "ensemble transform calculated for this (i,j) location in grid-0 (%s)", grid_getname(model_getgridbyid(das->m, 0)));
        ncw_put_att_text(fname, ncid, vid_transform, "long_name", attstr);
        ncw_def_var(fname, ncid, "X5_actual", NC_FLOAT, 2, dimids, &vid_transform_actual);
        snprintf(attstr, MAXSTRLEN, "the actual (interpolated) ensemble transform used at this (i,j) location in grid-0 (%s)", grid_getname(model_getvargrid(das->m, 0)));
        ncw_put_att_text(fname, ncid, vid_transform_actual, "long_name", attstr);
    } else if (das->mode == MODE_ENOI) {
        ncw_def_var(fname, ncid, "w", NC_DOUBLE, 1, &dimids[0], &vid_transform);
        ncw_def_var(fname, ncid, "w_actual", NC_FLOAT, 1, &dimids[0], &vid_transform_actual);
    }
    ncw_put_att_text(fname, ncid, NC_GLOBAL, "date", obs->datestr);
    ncw_put_att_text(fname, ncid, NC_GLOBAL, "grid_name", plog->gridname);
    ncw_put_att_int(fname, ncid, NC_GLOBAL, "i", 1, &plog->i);
    ncw_put_att_int(fname, ncid, NC_GLOBAL, "j", 1, &plog->j);
    ncw_put_att_double(fname, ncid, NC_GLOBAL, "lon", 1, &plog->lon);
    ncw_put_att_double(fname, ncid, NC_GLOBAL, "lat", 1, &plog->lat);
    ncw_put_att_double(fname, ncid, NC_GLOBAL, "depth", 1, &depth);

    ncw_enddef(fname, ncid);

    if (p > 0) {
        olon = malloc(p * sizeof(float));
        olat = malloc(p * sizeof(float));
        odepth = malloc(p * sizeof(float));
        val = malloc(p * sizeof(float));
        std = malloc(p * sizeof(float));
        fi = malloc(p * sizeof(float));
        fj = malloc(p * sizeof(float));
        fk = malloc(p * sizeof(float));
        otype = malloc(p * sizeof(int));
        odate = malloc(p * sizeof(float));

        for (oid = 0; oid < p; ++oid) {
            observation* o = &obs->data[lobs[oid]];

            olon[oid] = o->lon;
            olat[oid] = o->lat;
            odepth[oid] = o->depth;
            val[oid] = o->value;
            std[oid] = o->std;
            fi[oid] = o->fi;
            fj[oid] = o->fj;
            fk[oid] = o->fk;
            otype[oid] = o->type;
            odate[oid] = o->date;
        }

        ncw_put_var_float(fname, ncid, vid_lon, olon);
        ncw_put_var_float(fname, ncid, vid_lat, olat);
        ncw_put_var_float(fname, ncid, vid_depth, odepth);
        ncw_put_var_float(fname, ncid, vid_val, val);
        ncw_put_var_float(fname, ncid, vid_std, std);
        ncw_put_var_float(fname, ncid, vid_fi, fi);
        ncw_put_var_float(fname, ncid, vid_fj, fj);
        ncw_put_var_float(fname, ncid, vid_fk, fk);
        ncw_put_var_int(fname, ncid, vid_type, otype);
        ncw_put_var_float(fname, ncid, vid_date, odate);

        ncw_put_var_int(fname, ncid, vid_ids, lobs);
        ncw_put_var_double(fname, ncid, vid_lcoeffs, lcoeffs);
        ncw_put_var_double(fname, ncid, vid_s, s);
        ncw_put_var_double(fname, ncid, vid_S, S);

        free(olon);
        free(olat);
        free(odepth);
        free(val);
        free(std);
        free(fi);
        free(fj);
        free(fk);
        free(otype);
        free(odate);
    }

    ncw_put_var_double(fname, ncid, vid_transform, transform);

    ncw_close(fname, ncid);
}

/** Writes the "actual" (interpolated) ensemble transform to a pointlog file.
 */
void plog_writeactualtransform(dasystem* das, int id, float* transform)
{
    pointlog* plog = &das->plogs[id];
    char fname[MAXSTRLEN];
    int ncid;
    int vid = -1;

    das_getfname_plog(das, plog, fname);
    ncw_open(fname, NC_WRITE, &ncid);
    if (das->mode == MODE_ENKF)
        ncw_inq_varid(fname, ncid, "X5_actual", &vid);
    else if (das->mode == MODE_ENOI)
        ncw_inq_varid(fname, ncid, "w_actual", &vid);
    else
        enkf_quit("programming error");
    ncw_put_var_float(fname, ncid, vid, transform);
    ncw_close(fname, ncid);
}

/**
 */
static void get_nkname(dasystem* das, void* grid, char name[])
{
    if (model_getngrid(das->m) == 1)
        strcpy(name, "nk");
    else
        sprintf(name, "nk-%d", grid_getid(grid));
}

/**
 */
void plog_definestatevars(dasystem* das)
{
    int nvar = model_getnvar(das->m);
    int vid, p;

    if (das->nplogs == 0)
        return;

    for (vid = 0; vid < nvar; ++vid) {
        char* varname = model_getvarname(das->m, vid);
        char fname[MAXSTRLEN];
        char nkname[NC_MAX_NAME];
        char varname_an[NC_MAX_NAME];
        int nk;

        model_getmemberfname(das->m, das->ensdir, varname, 1, fname);
        nk = getnlevels(fname, varname);
        get_nkname(das, model_getvargrid(das->m, vid), nkname);

        for (p = 0; p < das->nplogs; ++p) {
            pointlog* plog = &das->plogs[p];
            int ncid;
            int dimid_nk, dimid_m;
            int dimids[2];
            int varid;

            das_getfname_plog(das, plog, fname);
            ncw_open(fname, NC_WRITE, &ncid);
            if (!ncw_var_exists(ncid, varname)) {
                ncw_redef(fname, ncid);
                ncw_inq_dimid(fname, ncid, "m", &dimid_m);
                if (!(das->updatespec & UPDATE_OUTPUTINC))
                    snprintf(varname_an, NC_MAX_NAME, "%s_an", varname);
                else
                    snprintf(varname_an, NC_MAX_NAME, "%s_inc", varname);
                if (nk > 1) {
                    if (!ncw_dim_exists(ncid, nkname))
                        ncw_def_dim(fname, ncid, nkname, nk, &dimid_nk);
                    else
                        ncw_inq_dimid(fname, ncid, nkname, &dimid_nk);
                    dimids[0] = dimid_nk;
                    dimids[1] = dimid_m;
                    ncw_def_var(fname, ncid, varname, NC_FLOAT, 2, dimids, &varid);
                    ncw_def_var(fname, ncid, varname_an, NC_FLOAT, 2, dimids, &varid);
                } else {
                    dimids[0] = dimid_m;
                    ncw_def_var(fname, ncid, varname, NC_FLOAT, 1, dimids, &varid);
                    ncw_def_var(fname, ncid, varname_an, NC_FLOAT, 1, dimids, &varid);
                }

                {
                    float inflation[2];
                    double tmp;

                    model_getvarinflation(das->m, vid, &inflation[0], &tmp);
                    inflation[1] = (float) tmp;
                    ncw_put_att_float(fname, ncid, varid, "INFLATION", 2, inflation);
                }
            }
            ncw_close(fname, ncid);
        }
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

    for (p = 0; p < das->nplogs; ++p) {
        pointlog* plog = &das->plogs[p];
        char fname[MAXSTRLEN];
        int ncid;

        das_getfname_plog(das, plog, fname);
        ncw_open(fname, NC_WRITE, &ncid);

        for (fid = 0; fid < nfields; ++fid) {
            field* f = &fields[fid];
            char varname[NC_MAX_NAME];
            int vid, ndims;

            v_src = (float***) fieldbuffer[fid];
            for (e = 0; e < das->nmem; ++e)
                v[e] = v_src[e][plog->j][plog->i];
            if (das->mode == MODE_ENOI)
                for (e = 0; e < das->nmem; ++e)
                    v[e] += v_src[das->nmem][plog->j][plog->i];

            if (!isanalysis)
                strncpy(varname, f->varname, NC_MAX_NAME);
            else {
                if (!(das->updatespec & UPDATE_OUTPUTINC))
                    snprintf(varname, NC_MAX_NAME, "%s%s", f->varname, "_an");
                else
                    snprintf(varname, NC_MAX_NAME, "%s%s", f->varname, "_inc");
            }
            ncw_inq_varid(fname, ncid, varname, &vid);
            ncw_inq_varndims(fname, ncid, vid, &ndims);
            if (ndims == 1)
                ncw_put_var_float(fname, ncid, vid, v);
            else if (ndims == 2) {
                start[0] = f->level;
                ncw_put_vara_float(fname, ncid, vid, start, count, v);
            }
        }

        ncw_close(fname, ncid);
    }

    free(v);
}

/** Writes state variables to tiles that are to be assembled into the pointlogs
 *  later.
 */
static void plog_writestatevars_toassemble(dasystem* das, int nfields, void** fieldbuffer, field* fields, int isanalysis)
{
    float* v = malloc(das->nmem * sizeof(float));
    int* vid = malloc(das->nplogs * sizeof(int));
    float*** v_src = NULL;
    int fid;

    for (fid = 0; fid < nfields; ++fid) {
        field* f = &fields[fid];
        char fname[MAXSTRLEN];
        int ncid, dimid;
        int p;

        snprintf(fname, MAXSTRLEN, "pointlog_%s-%03d.nc", f->varname, f->level);
        if (!isanalysis) {
            ncw_create(fname, NC_CLOBBER | NETCDF_FORMAT, &ncid);
            ncw_def_dim(fname, ncid, "m", das->nmem, &dimid);
        } else {
            ncw_open(fname, NC_WRITE, &ncid);
            ncw_inq_dimid(fname, ncid, "m", &dimid);
            ncw_redef(fname, ncid);
        }

        for (p = 0; p < das->nplogs; ++p) {
            pointlog* plog = &das->plogs[p];
            char varname[NC_MAX_NAME];

            if (!isanalysis)
                snprintf(varname, NC_MAX_NAME, "%s_%d_%d-%d", f->varname, plog->i, plog->j, plog->gridid);
            else {
                if (!(das->updatespec & UPDATE_OUTPUTINC))
                    snprintf(varname, NC_MAX_NAME, "%s%s_%d_%d-%d", f->varname, "_an", plog->i, plog->j, plog->gridid);
                else
                    snprintf(varname, NC_MAX_NAME, "%s%s_%d_%d-%d", f->varname, "_inc", plog->i, plog->j, plog->gridid);
            }
                    
            ncw_def_var(fname, ncid, varname, NC_FLOAT, 1, &dimid, &vid[p]);
        }
        ncw_enddef(fname, ncid);

        v_src = (float***) fieldbuffer[fid];

        for (p = 0; p < das->nplogs; ++p) {
            pointlog* plog = &das->plogs[p];
            int e;

            for (e = 0; e < das->nmem; ++e)
                v[e] = v_src[e][plog->j][plog->i];
            if (das->mode == MODE_ENOI)
                for (e = 0; e < das->nmem; ++e)
                    v[e] += v_src[das->nmem][plog->j][plog->i];

            ncw_put_var_float(fname, ncid, vid[p], v);
        }
        ncw_close(fname, ncid);
    }

    free(vid);
    free(v);
}

/**
 */
void plog_writestatevars(dasystem* das, int nfields, void** fieldbuffer, field* fields, int isanalysis)
{
    if (das->nplogs == 0)
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
    int p, fid;

    if (das->nplogs == 0)
        return;

    v = malloc(das->nmem * sizeof(float));

    das_getfields(das, -1, &nfields, &fields);

#if defined(MPI)
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    distribute_iterations(0, das->nplogs - 1, nprocesses, rank, "    ");
    for (p = my_first_iteration; p <= my_last_iteration; ++p) {
        pointlog* plog = &das->plogs[p];
        char fname_dst[MAXSTRLEN];
        int ncid_dst;
        size_t start[2] = { 0, 0 };
        size_t count[2] = { 1, das->nmem };

        das_getfname_plog(das, plog, fname_dst);
        ncw_open(fname_dst, NC_WRITE, &ncid_dst);

        for (fid = 0; fid < nfields; ++fid) {
            field* f = &fields[fid];
            char fname_src[MAXSTRLEN];
            int ii;

            snprintf(fname_src, MAXSTRLEN, "pointlog_%s-%03d.nc", f->varname, f->level);

            for (ii = 0; ii < 2; ++ii) {
                char varname[NC_MAX_NAME];
                char suffix[5] = "";
                int ncid_src, vid_src, vid_dst, ndims_dst;

                if (ii > 0) {
                    if (!(das->updatespec & UPDATE_OUTPUTINC))
                        strcpy(suffix, "_an");
                    else
                        strcpy(suffix, "_inc");
                }

                snprintf(varname, NC_MAX_NAME, "%s%s_%d_%d-%d", f->varname, suffix, plog->i, plog->j, plog->gridid);
                ncw_open(fname_src, NC_NOWRITE, &ncid_src);
                ncw_inq_varid(fname_src, ncid_src, varname, &vid_src);
                ncw_get_var_float(fname_src, ncid_src, vid_src, v);
                ncw_close(fname_src, ncid_src);

                snprintf(varname, NC_MAX_NAME, "%s%s", f->varname, suffix);
                ncw_inq_varid(fname_dst, ncid_dst, varname, &vid_dst);
                ncw_inq_varndims(fname_dst, ncid_dst, vid_dst, &ndims_dst);
                if (ndims_dst == 1)
                    ncw_put_var_float(fname_dst, ncid_dst, vid_dst, v);
                else if (ndims_dst == 2) {
                    start[0] = f->level;
                    ncw_put_vara_float(fname_dst, ncid_dst, vid_dst, start, count, v);
                }
            }
        }
        ncw_close(fname_dst, ncid_dst);
    }

#if defined(MPI)
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    enkf_printf("    deleting tiles:\n");
    distribute_iterations(0, nfields - 1, nprocesses, rank, "      ");
    for (fid = my_first_iteration; fid <= my_last_iteration; ++fid) {
        field* f = &fields[fid];
        char fname[MAXSTRLEN];

        snprintf(fname, MAXSTRLEN, "pointlog_%s-%03d.nc", f->varname, f->level);
        file_delete(fname);
    }

    free(v);
    free(fields);
}
