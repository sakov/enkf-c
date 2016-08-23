/******************************************************************************
 *
 * File:        reader_aquarius.c        
 *
 * Created:     6/2016
 *
 * Author:      Paul Sandery
 *              Bureau of Meteorology
 *
 * Description: Contains 1 reader reader_aquarius_standard() for preprocessed
 *              SSS data from AQUARIUS.
 *              Parameters:
 *               
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

#define ERRORSTD_DEF 0.25
#define MINDEPTH_DEF 200.0

void reader_aquarius_standard(char* fname, int fid, obsmeta* meta, model* m, observations* obs)
{
    double mindepth = MINDEPTH_DEF;
    char* addname = NULL;
    int ncid;
    int dimid_ni, dimid_nj;
    size_t ni, nj;
    int varid_lon, varid_lat, varid_sss;
    double* lon;
    double* lat;
    double** sss;
    int year, month, day;
    char tunits[MAXSTRLEN];
    double tunits_multiple, tunits_offset;
    char* basename;
    int model_vid;
    int i, j;
    int ktop;
    float** depth;
    double missval;

    for (i = 0; i < meta->npars; ++i) {
        if (strcasecmp(meta->pars[i].name, "MINDEPTH") == 0) {
            if (!str2double(meta->pars[i].value, &mindepth))
                enkf_quit("observation prm file: can not convert MINDEPTH = \"%s\" to double\n", meta->pars[i].value);
        } else if (strcasecmp(meta->pars[i].name, "ADD") == 0) {
            addname = meta->pars[i].value;
            enkf_printf("        ADDING \"%s\"\n", addname);
        } else
            enkf_quit("unknown PARAMETER \"%s\"\n", meta->pars[i].name);
    }
    enkf_printf("        MINDEPTH = %.0f\n", mindepth);

    basename = strrchr(fname, '/');
    if (basename == NULL)
        basename = fname;
    else
        basename += 1;

    ncw_open(fname, NC_NOWRITE, &ncid);
    ncw_inq_dimid(ncid, "longitude", &dimid_ni);
    ncw_inq_dimlen(ncid, dimid_ni, &ni);
    ncw_inq_dimid(ncid, "latitude", &dimid_nj);
    ncw_inq_dimlen(ncid, dimid_nj, &nj);
    enkf_printf("        (ni, nj) = (%u, %u)\n", ni, nj);

    ncw_inq_varid(ncid, "longitude", &varid_lon);
    lon = malloc(ni * sizeof(double));
    ncw_get_var_double(ncid, varid_lon, lon);

    ncw_inq_varid(ncid, "latitude", &varid_lat);
    lat = malloc(nj * sizeof(double));
    ncw_get_var_double(ncid, varid_lat, lat);

    ncw_inq_varid(ncid, "SSS", &varid_sss);
    sss = alloc2d(nj, ni, sizeof(double));
    ncw_get_var_double(ncid, varid_sss, sss[0]);
    ncw_get_att_double(ncid, varid_sss, "_FillValue", &missval);
    ncw_close(ncid);

    basename[9] = 0;
    if (!str2int(&basename[7], &day))
        enkf_quit("AQUARIUS reader: could not convert file name \"%s\" to date", fname);
    basename[7] = 0;
    if (!str2int(&basename[5], &month))
        enkf_quit("AQUARIUS reader: could not convert file name \"%s\" to date", fname);
    basename[5] = 0;
    if (!str2int(&basename[1], &year))
        enkf_quit("AQUARIUS reader: could not convert file name \"%s\" to date", fname);
    snprintf(tunits, MAXSTRLEN, "days since %4d-%02d-%02d", year, month, day);

    tunits_convert(tunits, &tunits_multiple, &tunits_offset);

    model_vid = model_getvarid(m, obs->obstypes[obstype_getid(obs->nobstypes, obs->obstypes, meta->type)].varnames[0], 1);
    depth = model_getdepth(m, model_vid, 1);
    ktop = grid_gettoplayerid(model_getvargrid(m, model_vid));

    for (j = 0; j < (int) nj; ++j) {
        for (i = 0; i < (int) ni; ++i) {
            observation* o;
            obstype* ot;

            obs_checkalloc(obs);
            o = &obs->data[obs->nobs];

            o->product = st_findindexbystring(obs->products, meta->product);
            assert(o->product >= 0);
            o->type = obstype_getid(obs->nobstypes, obs->obstypes, meta->type);
            assert(o->type >= 0);
            ot = &obs->obstypes[o->type];
            o->instrument = st_add_ifabscent(obs->instruments, "Aquarius", -1);
            o->id = obs->nobs;
            o->fid = fid;
            o->batch = 0;
            o->value = sss[j][i];
            if (o->value == missval)
                o->status = STATUS_MISSING;
            o->lon = lon[i];
            o->lat = lat[j];
            o->std = ERRORSTD_DEF;
            o->depth = 0.0;
            if (o->status == STATUS_OK)
                o->status = model_xy2fij(m, model_vid, o->lon, o->lat, &o->fi, &o->fj);
            if (!obs->allobs && o->status == STATUS_OUTSIDEGRID)
                continue;
            o->model_depth = (isnan(o->fi + o->fj)) ? NaN : depth[(int) (o->fj + 0.5)][(int) (o->fi + 0.5)];
            o->fk = 0.0;
            if (o->status == STATUS_OK && o->model_depth < mindepth)
                o->status = STATUS_SHALLOW;
            if ((o->status == STATUS_OK) && (o->lon <= ot->xmin || o->lon >= ot->xmax || o->lat <= ot->ymin || o->lat >= ot->ymax))
                o->status = STATUS_OUTSIDEOBSDOMAIN;

            o->fk = (double) ktop;
            o->date = tunits_offset + 0.5;
            o->aux = -1;

            obs->nobs++;
        }
    }

    free(lon);
    free(lat);
    free(sss);
}
