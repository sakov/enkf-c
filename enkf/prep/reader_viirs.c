/******************************************************************************
 *
 * File:        reader_viirs.c        
 *
 * Created:     08/03/2017
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description: Reader for pre-preprocessed (L3C) SST from VIIRS.
 *
 * Revisions:  
 *
 *****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "ncw.h"
#include "nan.h"
#include "definitions.h"
#include "utils.h"
#include "obsmeta.h"
#include "model.h"
#include "grid.h"
#include "observations.h"
#include "prep_utils.h"

#define ERRORSTD_DEF 0.5

void reader_viirs_standard(char* fname, int fid, obsmeta* meta, model* m, observations* obs)
{
    int ncid;
    int dimid_ni, dimid_nj;
    size_t ni, nj;
    int varid_lon, varid_lat, varid_sst, varid_npoints, varid_std, varid_estd, varid_time;
    double* lon;
    double* lat;
    short* sst;
    double sst_add_offset, sst_scale_factor;
    short *std;
    double std_add_offset, std_scale_factor;
    short* estd;
    double estd_add_offset, estd_scale_factor;
    short* npoints;
    unsigned char* time;
    double time_add_offset, time_scale_factor;
    char tunits[MAXSTRLEN];
    double tunits_multiple, tunits_offset;
    int mvid;
    float** depth;
    int ktop;
    int i, nobs;

    for (i = 0; i < meta->npars; ++i)
        enkf_quit("unknown PARAMETER \"%s\"\n", meta->pars[i].name);

    ncw_open(fname, NC_NOWRITE, &ncid);
    ncw_inq_dimid(ncid, "lon", &dimid_ni);
    ncw_inq_dimlen(ncid, dimid_ni, &ni);
    ncw_inq_dimid(ncid, "lat", &dimid_nj);
    ncw_inq_dimlen(ncid, dimid_nj, &nj);
    enkf_printf("        (ni, nj) = (%u, %u)\n", ni, nj);
    lon = malloc(ni * sizeof(double));
    lat = malloc(nj * sizeof(double));
    ncw_inq_varid(ncid, "lon", &varid_lon);
    ncw_inq_varid(ncid, "lat", &varid_lat);
    ncw_get_var_double(ncid, varid_lon, lon);
    ncw_get_var_double(ncid, varid_lat, lat);

    ncw_inq_varid(ncid, "sst", &varid_sst);
    sst = malloc(ni * nj * sizeof(short));
    ncw_get_var_short(ncid, varid_sst, sst);
    ncw_get_att_double(ncid, varid_sst, "add_offset", &sst_add_offset);
    ncw_get_att_double(ncid, varid_sst, "scale_factor", &sst_scale_factor);

    ncw_inq_varid(ncid, "std", &varid_std);
    std = malloc(ni * nj * sizeof(short));
    ncw_get_var_short(ncid, varid_std, std);
    ncw_get_att_double(ncid, varid_std, "add_offset", &std_add_offset);
    ncw_get_att_double(ncid, varid_std, "scale_factor", &std_scale_factor);

    ncw_inq_varid(ncid, "error_std", &varid_estd);
    estd = malloc(ni * nj * sizeof(short));
    ncw_get_var_short(ncid, varid_estd, estd);
    ncw_get_att_double(ncid, varid_estd, "add_offset", &estd_add_offset);
    ncw_get_att_double(ncid, varid_estd, "scale_factor", &estd_scale_factor);

    ncw_inq_varid(ncid, "npoints", &varid_npoints);
    npoints = malloc(ni * nj * sizeof(short));
    ncw_get_var_short(ncid, varid_npoints, npoints);

    ncw_inq_varid(ncid, "time", &varid_time);
    time = malloc(ni * nj);
    ncw_get_var_uchar(ncid, varid_time, time);
    ncw_get_att_double(ncid, varid_time, "add_offset", &time_add_offset);
    ncw_get_att_double(ncid, varid_time, "scale_factor", &time_scale_factor);
    ncw_get_att_text(ncid, varid_time, "units", tunits);

    ncw_close(ncid);

    tunits_convert(tunits, &tunits_multiple, &tunits_offset);

    mvid = model_getvarid(m, obs->obstypes[obstype_getid(obs->nobstypes, obs->obstypes, meta->type, 1)].varnames[0], 1);
    ktop = grid_gettoplayerid(model_getvargrid(m, mvid));
    depth = model_getdepth(m, mvid, 0);

    nobs = 0;
    for (i = 0; i < (int) (ni * nj); ++i) {
        observation* o;
        obstype* ot;

        if (npoints[i] == 0)
            continue;

        nobs++;
        obs_checkalloc(obs);
        o = &obs->data[obs->nobs];

        o->product = st_findindexbystring(obs->products, meta->product);
        assert(o->product >= 0);
        o->type = obstype_getid(obs->nobstypes, obs->obstypes, meta->type, 1);
        ot = &obs->obstypes[o->type];
        o->instrument = st_add_ifabscent(obs->instruments, "VIIRS", -1);
        o->id = obs->nobs;
        o->fid = fid;
        o->batch = 0;
        o->value = (double) sst[i] * sst_scale_factor + sst_add_offset - 273.15;
        {
            double std1 = (double) std[i] * std_scale_factor + std_add_offset;
            double std2 = (double) estd[i] * estd_scale_factor + estd_add_offset;
            
            o->std = (std1 > std2) ? std1 : std2;
        }
        o->lon = lon[i % ni];
        o->lat = lat[i / ni];
        o->depth = 0.0;
        o->fk = (double) ktop;
        o->status = model_xy2fij(m, mvid, o->lon, o->lat, &o->fi, &o->fj);
        if (!obs->allobs && o->status == STATUS_OUTSIDEGRID)
            continue;
        if ((o->status == STATUS_OK) && (o->lon <= ot->xmin || o->lon >= ot->xmax || o->lat <= ot->ymin || o->lat >= ot->ymax))
            o->status = STATUS_OUTSIDEOBSDOMAIN;
        o->model_depth = (depth == NULL || isnan(o->fi + o->fj)) ? NaN : depth[(int) (o->fj + 0.5)][(int) (o->fi + 0.5)];
        o->date = ((double) time[i] * time_scale_factor + time_add_offset) * tunits_multiple + tunits_offset;
        o->aux = -1;

        obs->nobs++;
    }
    enkf_printf("        nobs = %d\n", nobs);

    free(lon);
    free(lat);
    free(sst);
    free(std);
    free(estd);
    free(npoints);
    free(time);
}
