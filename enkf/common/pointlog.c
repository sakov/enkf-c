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
#include <assert.h>
#include "definitions.h"
#include "utils.h"
#include "dasystem.h"
#include "pointlog.h"

/**
 */
void plog_write(dasystem* das, int i, int j, double lon, double lat, double depth, int p, int* lobs, double* lcoeffs, double* s, double* S, double* transform)
{
    observations* obs = das->obs;
    char fname[MAXSTRLEN];
    int ncid;
    int dimids[2];
    int vid_ids, vid_S, vid_s, vid_lcoeffs, vid_transform, vid_lon, vid_lat, vid_depth, vid_val, vid_std, vid_fi, vid_fj, vid_fk, vid_type, vid_date;
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

    assert(das->s_mode == S_MODE_S_f);

    sprintf(fname, "pointlog_%d,%d.nc", i, j);
    ncw_create(fname, NC_CLOBBER | NC_64BIT_OFFSET, &ncid);
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
        sprintf(tunits, "days from %s", obs->datestr);
        ncw_put_att_text(fname, ncid, vid_date, "units", tunits);
        for (ot = 0; ot < obs->types->n; ++ot) {
            char attname[MAXSTRLEN];

            ncw_put_att_int(fname, ncid, vid_type, obs->obstypes[ot].name, 1, &ot);
            sprintf(attname, "RFACTOR_%s", obs->obstypes[ot].name);
            ncw_put_att_double(fname, ncid, vid_type, attname, 1, &obs->obstypes[ot].rfactor);
        }
    }
    if (das->mode == MODE_ENKF) {
        dimids[1] = dimids[0];
        ncw_def_var(fname, ncid, "X5", NC_DOUBLE, 2, dimids, &vid_transform);
    } else if (das->mode == MODE_ENOI) {
        ncw_def_var(fname, ncid, "w", NC_DOUBLE, 1, &dimids[0], &vid_transform);
    }
    ncw_put_att_text(fname, ncid, NC_GLOBAL, "date", obs->datestr);
    ncw_put_att_int(fname, ncid, NC_GLOBAL, "i", 1, &i);
    ncw_put_att_int(fname, ncid, NC_GLOBAL, "j", 1, &j);
    ncw_put_att_double(fname, ncid, NC_GLOBAL, "lon", 1, &lon);
    ncw_put_att_double(fname, ncid, NC_GLOBAL, "lat", 1, &lat);
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
    }

    ncw_put_var_double(fname, ncid, vid_transform, transform);

    ncw_close(fname, ncid);
}

/**
 */
void plog_definestatevars(dasystem* das)
{
    int nvar = model_getnvar(das->m);
    int vid, p;

    for (vid = 0; vid < nvar; ++vid) {
        char* varname = model_getvarname(das->m, vid);
        char fname[MAXSTRLEN];
        int nk;

        model_getmemberfname(das->m, das->ensdir, varname, 1, fname);
        nk = getnlevels(fname, varname);

        for (p = 0; p < das->nplogs; ++p) {
            int i = das->plogs[p].i;
            int j = das->plogs[p].j;
            int ncid;
            int dimid_nk, dimid_m;
            int dimids[2];
            int vid;

            sprintf(fname, "pointlog_%d,%d.nc", i, j);
            ncw_open(fname, NC_WRITE, &ncid);
            if (!ncw_var_exists(ncid, varname)) {
                ncw_redef(fname, ncid);
                ncw_inq_dimid(fname, ncid, "m", &dimid_m);
                if (nk > 1) {
                    if (!ncw_dim_exists(ncid, "nk"))
                        ncw_def_dim(fname, ncid, "nk", nk, &dimid_nk);
                    else
                        ncw_inq_dimid(fname, ncid, "nk", &dimid_nk);
                    dimids[0] = dimid_nk;
                    dimids[1] = dimid_m;
                    ncw_def_var(fname, ncid, varname, NC_FLOAT, 2, dimids, &vid);
                } else {
                    dimids[0] = dimid_m;
                    ncw_def_var(fname, ncid, varname, NC_FLOAT, 1, dimids, &vid);
                }
            }
            ncw_close(fname, ncid);
        }
    }
}

/**
 */
void plog_writestatevars(dasystem* das, int nfields, void** fieldbuffer, field* fields)
{
    int ni, nj, nk;
    int p, fid, e;
    float*** v_src = NULL;
    float* v = NULL;
    size_t start[2] = { 0, 0 };
    size_t count[2] = { 1, das->nmem };

    model_getdims(das->m, &ni, &nj, &nk);
    v = malloc(das->nmem * sizeof(float));

    for (p = 0; p < das->nplogs; ++p) {
        int i = das->plogs[p].i;
        int j = das->plogs[p].j;
        char fname[MAXSTRLEN];
        int ncid;

        sprintf(fname, "pointlog_%d,%d.nc", i, j);
        ncw_open(fname, NC_WRITE, &ncid);

        for (fid = 0; fid < nfields; ++fid) {
            field* f = &fields[fid];
            int vid, ndims;

            v_src = (float***) fieldbuffer[fid];
            for (e = 0; e < das->nmem; ++e)
                v[e] = v_src[e][j][i];
            if (das->mode == MODE_ENOI)
                for (e = 0; e < das->nmem; ++e)
                    v[e] += v_src[das->nmem][j][i];

            ncw_inq_varid(fname, ncid, f->varname, &vid);
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
