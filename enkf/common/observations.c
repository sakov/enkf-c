/******************************************************************************
 *
 * File:        observations.c        
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <values.h>
#include <assert.h>
#include "nan.h"
#include "stringtable.h"
#include "ncw.h"
#include "kdtree.h"
#include "definitions.h"
#include "utils.h"
#include "enkfprm.h"
#include "observations.h"

#define NOBSTYPES_INC 10
#define KD_INC 10000
#define PLAINDISTANCE 1

/* allowed range for obs of each type */
obstypedesc otdescs[] = {
    {"SST", 1, -2.5, 50.0},
    {"SLA", 1, -2.0, 2.0},
    {"TEM", 0, -2.5, 50.0},
    {"SAL", 0, 2.0, 41.0}
};

int notdescs = sizeof(otdescs) / sizeof(obstypedesc);

/**
 */
static double distance(double lon1, double lat1, double lon2, double lat2)
{
    double londiff, costheta;

    londiff = fabs(lon1 - lon2);
    if (londiff > abs(londiff - TWOPI))
        londiff -= TWOPI;

    costheta = sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(londiff);
    /*
     * around for round-off error 
     */
    if (fabs(costheta) > 1.0)
        costheta = (costheta > 0.0) ? 1.0 : -1.0;

    if (PLAINDISTANCE)
        return REARTH * sqrt(2.0 - 2.0 * costheta);
    else
        return REARTH * acos(costheta);
}

/**
 */
static double locfun(double x)
{
    double x2, x3;

    assert(x >= 0 && x <= 1.0 + 1.0e-8);

    if (x >= 1.0)               /* handle possible round-up error */
        return 0.0;

    x *= 2.0;
    x2 = x * x;
    x3 = x2 * x;
    if (x < 1.0)
        return 1.0 + x2 * (-x3 / 4.0 + x2 / 2.0) + x3 * (5.0 / 8.0) - x2 * (5.0 / 3.0);
    return x2 * (x3 / 12.0 - x2 / 2.0) + x3 * (5.0 / 8.0) + x2 * (5.0 / 3.0) - x * 5.0 + 4.0 - (2.0 / 3.0) / x;
}

/**
 */
void obs_addtype(observations* obs, char name[], int issurface, char varname[], char hfunction[], double rfactor, int isasync, double async_tstep)
{
    obstype* ot;

    if (obs->nobstypes % NOBSTYPES_INC == 0)
        obs->obstypes = realloc(obs->obstypes, (obs->nobstypes + NOBSTYPES_INC) * sizeof(obstype));
    ot = &obs->obstypes[obs->nobstypes];
    ot->id = obs->nobstypes;
    ot->name = strdup(name);
    ot->issurface = issurface;
    ot->varname = strdup(varname);
    ot->hfunction = strdup(hfunction);
    ot->rfactor = rfactor;
    ot->isasync = isasync;
    ot->async_tstep = async_tstep;
    ot->nobs = -1;
    ot->ngood = -1;
    ot->noutside = -1;
    ot->nland = -1;
    ot->nshallow = -1;
    ot->nrange = -1;
    ot->nmodified = 0;
    ot->date_min = NaN;
    ot->date_max = NaN;

    obs->nobstypes++;

    st_add_ifabscent(obs->types, name, ot->id);
}

/**
 */
observations* obs_create(void)
{
    observations* obs = malloc(sizeof(observations));

    obs->nobstypes = 0;
    obs->obstypes = NULL;
    obs->types = NULL;
    obs->products = NULL;
    obs->instruments = NULL;
    obs->nobstypes = 0;
    obs->obstypes = NULL;
    obs->da_date = NaN;
    obs->datestr = NULL;
    obs->allobs = 0;
    obs->nobs = 0;
    obs->data = NULL;
    obs->compacted = 0;
    obs->stride = 1;
    obs->hasstats = 0;
    obs->ngood = 0;
    obs->noutside = 0;
    obs->nland = 0;
    obs->nshallow = 0;
    obs->nrange = 0;
    obs->nmodified = 0;
    obs->tree = NULL;

    return obs;
}

/**
 */
observations* obs_create_fromprm(enkfprm* prm)
{
    observations* obs = obs_create();
    int i, j;

    obs->types = st_create("types");
    obs->products = st_create("products");
    obs->instruments = st_create("instruments");

    for (i = 0; i < prm->ntypes; ++i) {
        int isasync = 0;
        int issurface = -1;
        double tstep = NaN;

        for (j = 0; j < prm->nasync; ++j) {
            if (strcmp(prm->types[i], prm->async_types[j]) == 0) {
                isasync = 1;
                tstep = prm->async_timesteps[j];

                break;
            }
        }

        for (j = 0; j < notdescs; ++j) {
            if (strcmp(prm->types[i], otdescs[j].typename) == 0) {
                issurface = otdescs[j].issurface;
                break;
            }
        }
        assert(issurface >= 0);

        obs_addtype(obs, prm->types[i], issurface, prm->typevars[i], prm->hfunctions[i], prm->rfactors[i], isasync, tstep);
    }

    obs->da_date = date_str2dbl(prm->date);
    obs->datestr = strdup(prm->date);

    obs->stride = prm->sob_stride;

    return obs;
}

/**
 */
observations* obs_create_fromdata(observations* parentobs, int nobs, measurement data[])
{
    observations* obs = obs_create();
    int i;

    obs->types = st_copy(parentobs->types);
    obs->products = st_copy(parentobs->products);
    obs->instruments = st_copy(parentobs->instruments);

    for (i = 0; i < parentobs->nobstypes; ++i) {
        obstype* ot = &parentobs->obstypes[i];

        obs_addtype(obs, ot->name, ot->issurface, ot->varname, ot->hfunction, ot->rfactor, ot->isasync, ot->async_tstep);
    }

    obs->da_date = parentobs->da_date;
    obs->datestr = strdup(parentobs->datestr);

    obs->nobs = nobs;
    obs->data = data;

    return obs;
}

/**
 */
void obs_destroy(observations* obs)
{
    int i;

    st_destroy(obs->types);
    st_destroy(obs->products);
    st_destroy(obs->instruments);
    if (obs->nobstypes > 0) {
        for (i = 0; i < obs->nobstypes; ++i) {
            obstype* ot = &obs->obstypes[i];

            free(ot->name);
            free(ot->varname);
            free(ot->hfunction);
        }
        free(obs->obstypes);
    }
    if (obs->nobs > 0)
        free(obs->data);
    if (obs->datestr != NULL)
        free(obs->datestr);
    if (obs->tree != NULL)
        kd_destroy(obs->tree);
    free(obs);
}

/**
 */
void obs_checklon(observations* obs)
{
    int i;

    for (i = 0; i < obs->nobs; ++i) {
        measurement* o = &obs->data[i];

        if (o->lon < 0)
            o->lon += 360.0;
    }
}

/**
 */
static int comp_obsstatus(const void* p1, const void* p2)
{
    measurement* m1 = (measurement*) p1;
    measurement* m2 = (measurement*) p2;

    if (m1->status > m2->status)
        return 1;
    if (m1->status < m2->status)
        return -1;
    return 0;
}

/**
 */
void obs_compact(observations* obs)
{
    if (obs->compacted)
        return;

    enkf_printf("    compacting obs:");
    assert(STATUS_OK == 0);
    qsort(obs->data, obs->nobs, sizeof(measurement), comp_obsstatus);
    enkf_printf("\n");
    obs->compacted = 1;
}

/**
 */
void obs_calcstats(observations* obs)
{
    int i;

    if (obs->hasstats)
        return;

    enkf_printf("    calculating obs stats:");

    obs->ngood = 0;
    obs->noutside = 0;
    obs->nland = 0;
    obs->nshallow = 0;
    obs->nrange = 0;
    for (i = 0; i < obs->nobstypes; ++i) {
        obstype* ot = &obs->obstypes[i];

        ot->nobs = 0;
        ot->ngood = 0;
        ot->noutside = 0;
        ot->nland = 0;
        ot->nshallow = 0;
        ot->nrange = 0;
        ot->date_min = DBL_MAX;
        ot->date_max = -DBL_MAX;
    }

    for (i = 0; i < obs->nobs; ++i) {
        measurement* m = &obs->data[i];
        obstype* ot = &obs->obstypes[m->type];

        ot->nobs++;
        if (m->status == STATUS_OK) {
            obs->ngood++;
            ot->ngood++;
        } else if (m->status == STATUS_OUTSIDE) {
            obs->noutside++;
            ot->noutside++;
        } else if (m->status == STATUS_LAND) {
            obs->nland++;
            ot->nland++;
        } else if (m->status == STATUS_SHALLOW) {
            obs->nshallow++;
            ot->nshallow++;
        } else if (m->status == STATUS_RANGE) {
            obs->nrange++;
            ot->nrange++;
        }

        if (m->date < ot->date_min)
            ot->date_min = m->date;
        if (m->date > ot->date_max)
            ot->date_max = m->date;
    }
    obs->hasstats = 1;
    enkf_printf("\n");
}

/**
 */
void obs_write(observations* obs, char fname[])
{
    int nobs = obs->nobs;
    char tunits[MAXSTRLEN];

    int ncid;
    int dimid_nobs[1];
    int varid_type, varid_product, varid_instrument, varid_id, varid_idorig, varid_value, varid_std, varid_lon, varid_lat, varid_depth, varid_fi, varid_fj, varid_fk, varid_date, varid_status, varid_aux;

    int* type;
    int* product;
    int* instrument;
    int* id;
    int* id_orig;
    double* value;
    double* std;
    double* lon;
    double* lat;
    double* depth;
    double* fi;
    double* fj;
    double* fk;
    double* date;
    int* status;
    int* aux;

    int i, ii;

    if (rank != 0)
        return;

    if (file_exists(fname))
        enkf_quit("file \"%s\" already exists", fname);
    ncw_create(fname, NC_NOCLOBBER, &ncid);

    ncw_def_dim(fname, ncid, "nobs", nobs, dimid_nobs);
    ncw_def_var(fname, ncid, "type", NC_SHORT, 1, dimid_nobs, &varid_type);
    ncw_def_var(fname, ncid, "product", NC_SHORT, 1, dimid_nobs, &varid_product);
    ncw_def_var(fname, ncid, "instrument", NC_SHORT, 1, dimid_nobs, &varid_instrument);
    ncw_def_var(fname, ncid, "id", NC_INT, 1, dimid_nobs, &varid_id);
    ncw_def_var(fname, ncid, "id_orig", NC_INT, 1, dimid_nobs, &varid_idorig);
    ncw_def_var(fname, ncid, "value", NC_FLOAT, 1, dimid_nobs, &varid_value);
    ncw_def_var(fname, ncid, "std", NC_FLOAT, 1, dimid_nobs, &varid_std);
    ncw_def_var(fname, ncid, "lon", NC_FLOAT, 1, dimid_nobs, &varid_lon);
    ncw_def_var(fname, ncid, "lat", NC_FLOAT, 1, dimid_nobs, &varid_lat);
    ncw_def_var(fname, ncid, "depth", NC_FLOAT, 1, dimid_nobs, &varid_depth);
    ncw_def_var(fname, ncid, "fi", NC_FLOAT, 1, dimid_nobs, &varid_fi);
    ncw_def_var(fname, ncid, "fj", NC_FLOAT, 1, dimid_nobs, &varid_fj);
    ncw_def_var(fname, ncid, "fk", NC_FLOAT, 1, dimid_nobs, &varid_fk);
    ncw_def_var(fname, ncid, "date", NC_FLOAT, 1, dimid_nobs, &varid_date);
    ncw_def_var(fname, ncid, "status", NC_BYTE, 1, dimid_nobs, &varid_status);
    i = STATUS_OK;
    ncw_put_att_int(fname, ncid, varid_status, "STATUS_OK", 1, &i);
    i = STATUS_OUTSIDE;
    ncw_put_att_int(fname, ncid, varid_status, "STATUS_OUSIDE", 1, &i);
    i = STATUS_LAND;
    ncw_put_att_int(fname, ncid, varid_status, "STATUS_LAND", 1, &i);
    i = STATUS_SHALLOW;
    ncw_put_att_int(fname, ncid, varid_status, "STATUS_SHALLOW", 1, &i);
    i = STATUS_RANGE;
    ncw_put_att_int(fname, ncid, varid_status, "STATUS_RANGE", 1, &i);
    ncw_def_var(fname, ncid, "aux", NC_INT, 1, dimid_nobs, &varid_aux);
    sprintf(tunits, "days from %s", obs->datestr);
    ncw_put_att_text(fname, ncid, varid_date, "units", tunits);

    for (i = 0; i < obs->types->n; ++i)
        ncw_put_att_int(fname, ncid, varid_type, st_findstringbyindex(obs->types, i), 1, &i);

    for (i = 0; i < obs->products->n; ++i)
        ncw_put_att_int(fname, ncid, varid_product, st_findstringbyindex(obs->products, i), 1, &i);

    for (i = 0; i < obs->instruments->n; ++i)
        ncw_put_att_int(fname, ncid, varid_instrument, st_findstringbyindex(obs->instruments, i), 1, &i);

    ncw_enddef(fname, ncid);

    type = malloc(nobs * sizeof(int));
    product = malloc(nobs * sizeof(int));
    instrument = malloc(nobs * sizeof(int));
    id = malloc(nobs * sizeof(int));
    id_orig = malloc(nobs * sizeof(int));
    value = malloc(nobs * sizeof(double));
    std = malloc(nobs * sizeof(double));
    lon = malloc(nobs * sizeof(double));
    lat = malloc(nobs * sizeof(double));
    depth = malloc(nobs * sizeof(double));
    fi = malloc(nobs * sizeof(double));
    fj = malloc(nobs * sizeof(double));
    fk = malloc(nobs * sizeof(double));
    date = malloc(nobs * sizeof(double));
    status = malloc(nobs * sizeof(int));
    aux = malloc(nobs * sizeof(int));

    for (i = 0, ii = 0; i < obs->nobs; ++i) {
        measurement* m = &obs->data[i];

        if (!isfinite(m->value) || fabs(m->value) > FLT_MAX || !isfinite((float) m->value))
            enkf_quit("bad value");

        type[ii] = m->type;
        product[ii] = m->product;
        instrument[ii] = m->instrument;
        id[ii] = i;
        /*
         * id of the first ob contributed to this sob 
         */
        id_orig[ii] = m->id_orig;
        value[ii] = m->value;
        std[ii] = m->std;
        lon[ii] = m->lon;
        lat[ii] = m->lat;
        depth[ii] = m->depth;
        fi[ii] = m->fi;
        fj[ii] = m->fj;
        fk[ii] = m->fk;
        date[ii] = m->date;
        status[ii] = m->status;
        aux[ii] = m->aux;
        ii++;
    }
    assert(ii == nobs);

    ncw_put_var_int(fname, ncid, varid_type, type);
    ncw_put_var_int(fname, ncid, varid_product, product);
    ncw_put_var_int(fname, ncid, varid_instrument, instrument);
    ncw_put_var_int(fname, ncid, varid_id, id);
    ncw_put_var_int(fname, ncid, varid_idorig, id_orig);
    ncw_put_var_double(fname, ncid, varid_value, value);
    ncw_put_var_double(fname, ncid, varid_std, std);
    ncw_put_var_double(fname, ncid, varid_lon, lon);
    ncw_put_var_double(fname, ncid, varid_lat, lat);
    ncw_put_var_double(fname, ncid, varid_depth, depth);
    ncw_put_var_double(fname, ncid, varid_fi, fi);
    ncw_put_var_double(fname, ncid, varid_fj, fj);
    ncw_put_var_double(fname, ncid, varid_fk, fk);
    ncw_put_var_double(fname, ncid, varid_date, date);
    ncw_put_var_int(fname, ncid, varid_status, status);
    ncw_put_var_int(fname, ncid, varid_aux, aux);

    ncw_close(fname, ncid);
    free(type);
    free(product);
    free(instrument);
    free(id);
    free(id_orig);
    free(value);
    free(std);
    free(lon);
    free(lat);
    free(depth);
    free(fi);
    free(fj);
    free(fk);
    free(date);
    free(status);
    free(aux);
}

/**
 */
void obs_writeaux(observations* obs, char fname[])
{
    int ncid;
    int dimid_nobs;
    int varid_aux;
    size_t n;
    int i, ii;
    int* aux;

    if (rank != 0)
        return;

    ncw_open(fname, NC_WRITE, &ncid);
    ncw_inq_varid(fname, ncid, "aux", &varid_aux);
    ncw_inq_dimid(fname, ncid, "nobs", &dimid_nobs);
    ncw_inq_dimlen(fname, ncid, dimid_nobs, &n);
    assert((int) n == obs->nobs);

    aux = malloc(n * sizeof(int));
    for (i = 0, ii = 0; i < obs->nobs; ++i) {
        aux[ii] = obs->data[i].aux;
        ii++;
    }

    ncw_put_var_int(fname, ncid, varid_aux, aux);
    ncw_close(fname, ncid);
    free(aux);
}

/**
 */
void obs_superob(observations* obs, __compar_d_fn_t cmp_obs, observations** sobs, int sobid)
{
    int i1 = 0, i2 = 0;
    int nsobs = 0;
    measurement* data = obs->data;
    measurement* sdata = NULL;

    obs_calcstats(obs);
    obs_compact(obs);

    qsort_r(obs->data, obs->ngood, sizeof(measurement), cmp_obs, obs);

    while (i2 < obs->ngood) {
        measurement* so;
        measurement* o;
        double n;
        double lon_min, lon_max;
        int ii;

        /*
         * identify obs that will be combined into this superob 
         */
        while (i2 + 1 < obs->nobs && cmp_obs(&data[i1], &data[i2 + 1], obs) == 0)
            i2++;
        if (nsobs % NOBS_INC == 0)
            sdata = realloc(sdata, (nsobs + NOBS_INC) * sizeof(measurement));

        if (sobid == nsobs) {
            int i;

            enkf_printf("    sob # %d is formed by the following obs:\n", sobid);
            for (i = i1; i <= i2; ++i) {
                enkf_printf("      ");
                obs_printob(obs, i);
            }
        }

        so = &sdata[nsobs];
        o = &data[i1];

        assert(o->status == STATUS_OK);

        so->type = o->type;
        so->product = o->product;
        so->instrument = o->instrument;
        so->id = nsobs;
        so->id_orig = o->id;
        so->value = o->value;
        so->std = 1.0 / (o->std * o->std);
        so->lon = o->lon;
        so->lat = o->lat;
        so->depth = o->depth;
        so->fi = o->fi;
        so->fj = o->fj;
        so->fk = o->fk;
        so->date = o->date;
        so->status = o->status;
        so->aux = 1;
        n = 1.0;

        lon_min = o->lon;
        lon_max = o->lon;

        for (ii = i1; ii <= i2; ++ii)
            data[ii].aux = nsobs;

        for (ii = i1 + 1; ii <= i2; ++ii) {
            o = &data[ii];
            if (so->product != o->product)
                so->product = -1;
            if (so->instrument != o->instrument)
                so->instrument = -1;

            so->value = so->value * so->std + o->value / o->std / o->std;
            so->lon = so->lon * so->std + o->lon / o->std / o->std;
            so->lat = so->lat * so->std + o->lat / o->std / o->std;
            so->depth = so->depth * so->std + o->depth / o->std / o->std;
            so->fi = so->fi * so->std + o->fi / o->std / o->std;
            so->fj = so->fj * so->std + o->fj / o->std / o->std;
            so->fk = so->fk * so->std + o->fk / o->std / o->std;
            so->date = so->date * so->std + o->date / o->std / o->std;

            so->std += 1.0 / (o->std * o->std);

            so->value /= so->std;
            so->lon /= so->std;
            so->lat /= so->std;
            so->depth /= so->std;
            so->fi /= so->std;
            so->fj /= so->std;
            so->fk /= so->std;
            so->date /= so->std;
            n++;

            if (o->lon < lon_min)
                lon_min = o->lon;
            if (o->lon > lon_max)
                lon_max = o->lon;

            assert(o->status == STATUS_OK);
        }

        so->std = sqrt(1.0 / so->std);
        if (lon_max - lon_min > 180.0) {
            /*
             * (there is a possibility of merging observations separated by
             * more than 360 degrees in longitude) 
             */
            so->lon = 0.0;
            for (ii = i1; ii < i2; ++ii) {
                o = &data[ii];
                if (o->lon - lon_min > 180.0)
                    o->lon -= 360.0;
                so->lon += o->lon / o->std / o->std;
            }
            so->lon *= so->std * so->std;
            if (so->lon < 0.0)
                so->lon += 360.0;
        }
        so->aux = n;

        nsobs += 1;

        i1 = i2 + 1;
        i2 = i1;
    }
    enkf_printf("    %d superoobservations\n", nsobs);

    *sobs = obs_create_fromdata(obs, nsobs, sdata);
    obs_calcstats(*sobs);
}

/**
 */
void obs_find_bytype(observations* obs, int type, int* nobs, int** obsids, int fstatsonly)
{
    int i;

    if (!fstatsonly)
        assert(obs->obstypes[type].nobs == obs->obstypes[type].ngood);

    *nobs = 0;
    if (obs->obstypes[type].nobs == 0) {
        *obsids = NULL;
        return;
    }
    *obsids = malloc(obs->obstypes[type].nobs * sizeof(int));
    for (i = 0; i < obs->nobs; ++i) {
        measurement* o = &obs->data[i];

        if (o->type == type && o->status == STATUS_OK) {
            (*obsids)[*nobs] = i;
            (*nobs)++;
        }
    }
    if (*nobs == 0) {
        free(*obsids);
        *obsids = NULL;
    }
}

/**
 */
void obs_find_bytypeandtime(observations* obs, int type, int time, int* nobs, int** obsids, int fstatsonly)
{
    double tstep = obs->obstypes[type].async_tstep;
    int i;

    if (!fstatsonly)
        assert(obs->obstypes[type].nobs == obs->obstypes[type].ngood);

    *nobs = 0;
    if (obs->obstypes[type].ngood == 0) {
        *obsids = NULL;
        return;
    }
    *obsids = malloc(obs->obstypes[type].ngood * sizeof(int));
    for (i = 0; i < obs->nobs; ++i) {
        measurement* o = &obs->data[i];

        if (o->type == type && o->status == STATUS_OK && get_tshift(o->date, tstep) == time) {
            (*obsids)[*nobs] = i;
            (*nobs)++;
        }
    }
    if (*nobs == 0) {
        free(*obsids);
        *obsids = NULL;
    }
}

/**
 */
void obs_printob(observations* obs, int i)
{
    measurement* o = &obs->data[i];

    enkf_printf("type = %s, product = %s, instrument = %s, id = %d, original id = %d, value = %.3g, std = %.3g, lon = %.3f, lat = %.3f, depth = %.1f, fi = %.3f, fj = %.3f, fk = %.3f, date = %.3g, status = %d\n", obs->obstypes[o->type].name, st_findstringbyindex(obs->products, o->product), st_findstringbyindex(obs->instruments, o->instrument), o->id, o->id_orig, o->value, o->std, o->lon, o->lat, o->depth, o->fi, o->fj, o->fk, o->date, o->status);
}

/**
 */
void obs_createkdtree(observations* obs)
{
    int i;

    if (obs->tree != NULL)
        kd_destroy(obs->tree);

    obs->tree = kd_create(3);
    for (i = 0; i < obs->nobs; ++i) {
        measurement* o = &obs->data[i];
        double lat = o->lat * DEG2RAD;
        double lon = o->lon * DEG2RAD;
        double coslat = cos(lat);
        double point[3];

        point[0] = REARTH * sin(lon) * coslat;
        point[1] = REARTH * cos(lon) * coslat;
        point[2] = REARTH * sin(lat);
        kd_insert(obs->tree, point);
    }
}

/**
 */
void obs_findlocal(observations* obs, double lon, double lat, double r, int* n, int** ids, double** lcoeffs)
{
    double point[3];
    kdset* set = NULL;
    int i, id;

    lon *= DEG2RAD;
    lat *= DEG2RAD;

    point[0] = REARTH * sin(lon) * cos(lat);
    point[1] = REARTH * cos(lon) * cos(lat);
    point[2] = REARTH * sin(lat);

    set = kd_nearest_range(obs->tree, point, r, 1);
    for (i = 0; (id = kd_res_item_getid(set)) >= 0; kd_res_next(set), ++i) {
        if (i % KD_INC == 0) {
            *ids = realloc(*ids, (i + KD_INC) * sizeof(int));
            *lcoeffs = realloc(*lcoeffs, (i + KD_INC) * sizeof(double));
        }
        (*ids)[i] = id;
    }
    kd_res_free(set);
    *n = i;

    for (i = 0; i < *n; ++i) {
        measurement* o = &obs->data[(*ids)[i]];

        (*lcoeffs)[i] = locfun(distance(lon, lat, o->lon * DEG2RAD, o->lat * DEG2RAD) / r);
    }
}
