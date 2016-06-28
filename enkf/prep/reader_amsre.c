/******************************************************************************
 *
 * File:        reader_amsre.c        
 *
 * Created:     09/07/2015
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description: AMSRE reader. Note that (at the moment) this is the only
 *              reader for grided data in EnKF-C.
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

#define ORBITS_ALL        0
#define ORBITS_DESCENDING 1
#define ORBITS_ASCENDING  2
#define ORBITS_DEF        ORBITS_ALL

#define ERRORSTD_DEF 0.25
#define MINWIND_DEF 5.0

void reader_amsre_standard(char* fname, int fid, obsmeta* meta, model* m, observations* obs)
{
    double minwind = MINWIND_DEF;
    int orbits = ORBITS_DEF;
    int ncid;
    int dimid_ni, dimid_nj;
    size_t ni, nj;
    int varid_lon, varid_lat, varid_sst_a, varid_sst_d, varid_wind_a, varid_wind_d, varid_time_a, varid_time_d;
    double* lon;
    double* lat;
    double** sst_a;
    double** sst_d;
    double** wind_a;
    double** wind_d;
    double** time_a;
    double** time_d;
    int year, month, day;
    char tunits[MAXSTRLEN];
    size_t tunits_len;
    double tunits_multiple, tunits_offset;
    char* basename;
    int model_vid;
    int i, j, k;
    int channel;

    for (i = 0; i < meta->npars; ++i) {
        if (strcasecmp(meta->pars[i].name, "MINWIND") == 0) {
            if (!str2double(meta->pars[i].value, &minwind))
                enkf_quit("observation prm file: can not convert MINWIND = \"%s\" to double\n", meta->pars[i].value);
        } else if (strcasecmp(meta->pars[i].name, "ORBITS") == 0) {
            if (strcasecmp(meta->pars[i].value, "ALL") == 0)
                orbits = ORBITS_ALL;
            else if (strcasecmp(meta->pars[i].value, "DESCENDING") == 0)
                orbits = ORBITS_DESCENDING;
            else if (strcasecmp(meta->pars[i].value, "ASCENDING") == 0)
                orbits = ORBITS_ASCENDING;
            else
                enkf_quit("observation prm file: parameter \"ORBITS\": value \"%s\" not understood: expected either \"ALL\", \"DESCENDING\", or \"ASCENDING\"\n", meta->pars[i].value);
        } else
            enkf_quit("unknown PARAMETER \"%s\"\n", meta->pars[i].name);
    }
    enkf_printf("        MINWIND = %.0f\n", minwind);
    enkf_printf("        ORBITS = ");
    if (orbits == ORBITS_ALL)
        enkf_printf("ALL\n");
    else if (orbits == ORBITS_DESCENDING)
        enkf_printf("DESCENDING\n");
    else if (orbits == ORBITS_ASCENDING)
        enkf_printf("ASCENDING\n");
    else
        enkf_quit("programming error");

    basename = strrchr(fname, '/');
    if (basename == NULL)
        basename = fname;
    else
        basename += 1;

    ncw_open(fname, NC_NOWRITE, &ncid);
    ncw_inq_dimid(fname, ncid, "lon", &dimid_ni);
    ncw_inq_dimlen(fname, ncid, dimid_ni, &ni);
    ncw_inq_dimid(fname, ncid, "lat", &dimid_nj);
    ncw_inq_dimlen(fname, ncid, dimid_nj, &nj);
    enkf_printf("        (ni, nj) = (%u, %u)\n", ni, nj);

    ncw_inq_varid(fname, ncid, "lon", &varid_lon);
    lon = malloc(ni * sizeof(double));
    ncw_get_var_double(fname, ncid, varid_lon, lon);

    ncw_inq_varid(fname, ncid, "lat", &varid_lat);
    lat = malloc(nj * sizeof(double));
    ncw_get_var_double(fname, ncid, varid_lat, lat);

    ncw_inq_varid(fname, ncid, "a_sst", &varid_sst_a);
    sst_a = alloc2d(nj, ni, sizeof(double));
    ncw_get_var_double(fname, ncid, varid_sst_a, sst_a[0]);

    ncw_inq_varid(fname, ncid, "d_sst", &varid_sst_d);
    sst_d = alloc2d(nj, ni, sizeof(double));
    ncw_get_var_double(fname, ncid, varid_sst_d, sst_d[0]);

    ncw_inq_varid(fname, ncid, "a_wind", &varid_wind_a);
    wind_a = alloc2d(nj, ni, sizeof(double));
    ncw_get_var_double(fname, ncid, varid_wind_a, wind_a[0]);

    ncw_inq_varid(fname, ncid, "d_wind", &varid_wind_d);
    wind_d = alloc2d(nj, ni, sizeof(double));
    ncw_get_var_double(fname, ncid, varid_wind_d, wind_d[0]);

    ncw_inq_varid(fname, ncid, "a_time", &varid_time_a);
    time_a = alloc2d(nj, ni, sizeof(double));
    ncw_get_var_double(fname, ncid, varid_time_a, time_a[0]);

    ncw_inq_varid(fname, ncid, "d_time", &varid_time_d);
    time_d = alloc2d(nj, ni, sizeof(double));
    ncw_get_var_double(fname, ncid, varid_time_d, time_d[0]);

    ncw_inq_attlen(fname, ncid, varid_time_a, "units", &tunits_len);
    ncw_get_att_text(fname, ncid, varid_time_a, "units", tunits);
    basename[14] = 0;
    if (!str2int(&basename[12], &day))
        enkf_quit("AMSRE reader: could not convert file name \"%s\" to date", fname);
    basename[12] = 0;
    if (!str2int(&basename[10], &month))
        enkf_quit("AMSRE reader: could not convert file name \"%s\" to date", fname);
    basename[10] = 0;
    if (!str2int(&basename[6], &year))
        enkf_quit("AMSRE reader: could not convert file name \"%s\" to date", fname);
    snprintf(&tunits[tunits_len], MAXSTRLEN - tunits_len, " since %4d-%02d-%02d", year, month, day);

    ncw_close(fname, ncid);

    tunits_convert(tunits, &tunits_multiple, &tunits_offset);

    model_vid = model_getvarid(m, obs->obstypes[obstype_getid(obs->nobstypes, obs->obstypes, meta->type)].varnames[0], 1);
    k = grid_gettoplayerid(model_getvargrid(m, model_vid));

    for (channel = 0; channel < 2; ++channel) {
        double** data = (channel == 0) ? sst_a : sst_d;
        double** wind = (channel == 0) ? wind_a : wind_d;
        double** time = (channel == 0) ? time_a : time_d;

        if ((orbits == ORBITS_DESCENDING && channel == 0) || (orbits == ORBITS_ASCENDING && channel == 1))
            continue;

        for (j = 0; j < (int) nj; ++j) {
            for (i = 0; i < (int) ni; ++i) {
                observation* o;
                obstype* ot;

                if (fabs(data[j][i]) > MAXOBSVAL)
                    continue;
                if (wind[j][i] > MAXOBSVAL || wind[j][i] < minwind)
                    continue;   /* (do not need fabs, as wind non-negative) */

                obs_checkalloc(obs);
                o = &obs->data[obs->nobs];

                o->product = st_findindexbystring(obs->products, meta->product);
                assert(o->product >= 0);
                o->type = obstype_getid(obs->nobstypes, obs->obstypes, meta->type);
                assert(o->type >= 0);
                ot = &obs->obstypes[o->type];
                o->instrument = st_add_ifabscent(obs->instruments, "AMSRE", -1);
                o->id = obs->nobs;
                o->fid = fid;
                o->batch = channel;
                o->value = data[j][i];
                o->std = ERRORSTD_DEF;
                o->lon = lon[i];
                o->lat = lat[j];
                o->depth = 0.0;
                o->status = model_xy2fij(m, model_vid, o->lon, o->lat, &o->fi, &o->fj);
                if (!obs->allobs && o->status == STATUS_OUTSIDEGRID)
                    continue;
                if ((o->status == STATUS_OK) && (o->lon <= ot->xmin || o->lon >= ot->xmax || o->lat <= ot->ymin || o->lat >= ot->ymax || o->depth <= ot->zmin || o->depth >= ot->zmax))
                    o->status = STATUS_OUTSIDEOBSDOMAIN;
                o->fk = (double) k;
                o->date = time[j][i] * tunits_multiple + tunits_offset;
                o->aux = -1;

                obs->nobs++;
            }
        }
    }

    free(lon);
    free(lat);
    free(sst_a);
    free(sst_d);
    free(wind_a);
    free(wind_d);
    free(time_a);
    free(time_d);
}
