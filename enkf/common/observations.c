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
#include <float.h>
#include <limits.h>
#include <assert.h>
#include <stdint.h>
#include "ncw.h"
#include "definitions.h"
#include "utils.h"
#include "observations.h"

#define NOBSTYPES_INC 10
#define PLOC_INC 10000
#define HT_SIZE 500
#define EPSF 1.0e-5

typedef struct {
    int obstypeid;
    char fname[MAXSTRLEN];
    int fid;
    int batch;
} badbatch;

typedef union {
    int key_int[2];
    short int key_short[4];
} keydata;

/**
 */
void obs_addtype(observations* obs, obstype* src)
{
    obstype* ot;
    int i;

    if (obs->nobstypes % NOBSTYPES_INC == 0)
        obs->obstypes = realloc(obs->obstypes, (obs->nobstypes + NOBSTYPES_INC) * sizeof(obstype));
    ot = &obs->obstypes[obs->nobstypes];
    ot->id = obs->nobstypes;
    ot->name = strdup(src->name);
    ot->issurface = src->issurface;
    ot->statsonly = src->statsonly;
    ot->nvar = src->nvar;
    ot->varnames = malloc(src->nvar * sizeof(char*));
    for (i = 0; i < src->nvar; ++i)
        ot->varnames[i] = strdup(src->varnames[i]);
    ot->alias = strdup(src->alias);
    ot->logapplied = src->logapplied;
    ot->offset_fname = (src->offset_fname != NULL) ? strdup(src->offset_fname) : NULL;
    ot->offset_varname = (src->offset_varname != NULL) ? strdup(src->offset_varname) : NULL;
    ot->mld_varname = (src->mld_varname != NULL) ? strdup(src->mld_varname) : NULL;
    ot->mld_threshold = src->mld_threshold;
    ot->hfunction = strdup(src->hfunction);
    ot->allowed_min = src->allowed_min;
    ot->allowed_max = src->allowed_max;
    ot->isasync = src->isasync;
    ot->async_tstep = src->async_tstep;
    ot->async_centred = src->async_centred;
    ot->async_tname = (src->async_tname == NULL) ? NULL : strdup(src->async_tname);
    ot->nlocrad = src->nlocrad;
    ot->locrad = malloc(sizeof(double) * ot->nlocrad);
    ot->locweight = malloc(sizeof(double) * ot->nlocrad);
    for (i = 0; i < ot->nlocrad; ++i) {
        ot->locrad[i] = src->locrad[i];
        ot->locweight[i] = src->locweight[i];
    }
    ot->rfactor = src->rfactor;
    ot->nlobsmax = src->nlobsmax;
    ot->estdmin = src->estdmin;
    ot->vid = src->vid;
    ot->gridid = src->gridid;
    ot->sob_stride = src->sob_stride;

    ot->obswindow_min = src->obswindow_min;
    ot->obswindow_max = src->obswindow_max;
    ot->time_min = src->time_min;
    ot->time_max = src->time_max;

    ot->xmin = src->xmin;
    ot->xmax = src->xmax;
    ot->ymin = src->ymin;
    ot->ymax = src->ymax;
    ot->zmin = src->zmin;
    ot->zmax = src->zmax;

    ot->ndomains = src->ndomains;
    if (ot->ndomains > 0)
        ot->domainnames = malloc(ot->ndomains * sizeof(char*));
    for (i = 0; i < ot->ndomains; ++i)
        ot->domainnames[i] = strdup(src->domainnames[i]);

    ot->nsubgrid = src->nsubgrid;
    ot->nmodified = 0;
    /*
     * (these fields are set by obs_calcstats())
     */
    ot->nobs = -1;
    ot->ngood = -1;
    ot->noutside_grid = -1;
    ot->noutside_obsdomain = -1;
    ot->noutside_obswindow = -1;
    ot->nland = -1;
    ot->nshallow = -1;
    ot->nbadbatch = -1;
    ot->nbadfc = -1;
    ot->nrange = -1;
    ot->nthinned = -1;
    ot->nexcluded = -1;

    obs->nobstypes++;
}

/**
 */
observations* obs_create(void)
{
    observations* obs = malloc(sizeof(observations));

    obs->products = NULL;
    obs->instruments = NULL;
    obs->datafiles = NULL;
    obs->nobstypes = 0;
    obs->obstypes = NULL;
#if defined(ENKF_CALC)
    obs->loctrees = NULL;
    obs->obsids = NULL;
#if defined(USE_SHMEM)
    obs->sm_comm_wins_kd = NULL;
    obs->sm_comm_win_data = MPI_WIN_NULL;
#endif
#endif
    obs->da_time = NAN;
    obs->datestr = NULL;
    obs->allobs = 0;
    obs->nallocated = 0;
    obs->nobs = 0;
    obs->data = NULL;
    obs->compacted = 0;
    obs->has_nonpointobs = 0;
    obs->ngood = 0;
    obs->noutside_grid = 0;
    obs->noutside_obsdomain = 0;
    obs->noutside_obswindow = 0;
    obs->nland = 0;
    obs->nshallow = 0;
    obs->nthinned = 0;
    obs->nbadbatch = 0;
    obs->nrange = 0;
    obs->nmodified = 0;
    obs->badbatches = NULL;
    obs->ncformat = NETCDF_FORMAT;
    obs->nccompression = 0;
#if defined(ENKF_PREP)
    obs->model = NULL;
#endif

    return obs;
}

#if defined(ENKF_PREP) || defined(ENKF_CALC)
/**
 */
observations* obs_create_fromprm(enkfprm* prm)
{
    observations* obs = obs_create();

    obs->products = st_create("products");
    obs->instruments = st_create("instruments");
    obs->datafiles = st_create("datafiles");
    obs->da_time = date2day(prm->fname, prm->date);
    obs->datestr = strdup(prm->date);

    enkf_printf("  reading observation type specs from \"%s\":\n", prm->obstypeprm);
    obstypes_read(prm, prm->obstypeprm, &obs->nobstypes, &obs->obstypes);

    {
        int otid;

        for (otid = 0; otid < obs->nobstypes; ++otid) {
            obstype* ot = &obs->obstypes[otid];

            if (isnan(ot->obswindow_min))
                ot->obswindow_min = prm->obswindow_min;
            if (isnan(ot->obswindow_max))
                ot->obswindow_max = prm->obswindow_max;
        }
    }

    obs->ncformat = prm->ncformat;
    obs->nccompression = prm->nccompression;

    return obs;
}
#endif

/**
 */
observations* obs_create_fromdata(observations* parentobs, int nobs, observation data[])
{
    observations* obs = obs_create();
    int i;

    obs->products = st_copy(parentobs->products);
    obs->instruments = st_copy(parentobs->instruments);
    obs->datafiles = st_copy(parentobs->datafiles);

    for (i = 0; i < parentobs->nobstypes; ++i)
        obs_addtype(obs, &parentobs->obstypes[i]);

    obs->da_time = parentobs->da_time;
    obs->datestr = strdup(parentobs->datestr);

    obs->nobs = nobs;
    obs->data = data;

    obs->ncformat = parentobs->ncformat;
    obs->nccompression = parentobs->nccompression;

    obs_calcstats(obs);

    return obs;
}

/**
 */
void obs_destroy(observations* obs)
{
    st_destroy(obs->products);
    st_destroy(obs->instruments);
    st_destroy(obs->datafiles);
    obstypes_destroy(obs->nobstypes, obs->obstypes);
#if defined(ENKF_CALC)
    obs_destroykdtrees(obs);
    if (obs->obsids != NULL) {
        int i;

        for (i = 0; i < obs->nobstypes; ++i)
            if (obs->obsids[i] != NULL)
                free(obs->obsids[i]);
        free(obs->obsids);
    }
#endif
#if defined(USE_SHMEM)
    if (obs->sm_comm_win_data != MPI_WIN_NULL) {
        MPI_Win_free(&obs->sm_comm_win_data);
        assert(obs->sm_comm_win_data == MPI_WIN_NULL);
    }
#else
    if (obs->nobs > 0)
        free(obs->data);
#endif
    if (obs->datestr != NULL)
        free(obs->datestr);
    if (obs->badbatches != NULL) {
        ht_process(obs->badbatches, free);
        ht_destroy(obs->badbatches);
    }
    free(obs);
}

#if defined(ENKF_PREP)
/**
 */
void obs_checkalloc(observations* obs)
{
    if (obs->nobs == obs->nallocated) {
        obs->data = realloc(obs->data, (obs->nobs + NOBS_INC) * sizeof(observation));
        if (obs->data == NULL)
            enkf_quit("obs_checkalloc(): not enough memory");
        obs->nallocated += NOBS_INC;
        memset(&obs->data[obs->nobs], 0, NOBS_INC * sizeof(observation));
    }
}
#endif

#if defined(ENKF_PREP)
/** Move good observations to the head of the observation array.
 */
void obs_compact(observations* obs)
{
    int i, ii;

    if (obs->compacted)
        return;

    enkf_flush();
    assert(STATUS_OK == 0);

    for (i = 0, ii = obs->nobs - 1; i < ii; ++i) {
        observation tmp;

        if (obs->data[i].status == 0)
            continue;
        while (obs->data[ii].status != 0 && ii > i)
            ii--;
        tmp = obs->data[i];
        obs->data[i] = obs->data[ii];
        obs->data[ii] = tmp;
        obs->data[i].id = i;
        obs->data[ii].id = ii;
        ii--;
    }
    obs->compacted = 1;
}
#endif

/**
 */
static int comp_obsid(const void* p1, const void* p2)
{
    observation* m1 = (observation*) p1;
    observation* m2 = (observation*) p2;

    if (m1->id > m2->id)
        return 1;
    if (m1->id < m2->id)
        return -1;
    return 0;
}

/** Sort observations by id.
 */
void obs_inorder(observations* obs)
{
    qsort(obs->data, obs->nobs, sizeof(observation), comp_obsid);
}

/**
 */
void obs_calcstats(observations* obs)
{
    int i;

    obs->ngood = 0;
    obs->noutside_grid = 0;
    obs->noutside_obsdomain = 0;
    obs->noutside_obswindow = 0;
    obs->nland = 0;
    obs->nbadbatch = 0;
    obs->nbadfc = 0;
    obs->nshallow = 0;
    obs->nrange = 0;
    obs->nthinned = 0;
    obs->nexcluded = 0;
    for (i = 0; i < obs->nobstypes; ++i) {
        obstype* ot = &obs->obstypes[i];

        ot->nobs = 0;
        ot->ngood = 0;
        ot->noutside_grid = 0;
        ot->noutside_obsdomain = 0;
        ot->noutside_obswindow = 0;
        ot->nland = 0;
        ot->nshallow = 0;
        ot->nbadbatch = 0;
        ot->nbadfc = 0;
        ot->nrange = 0;
        ot->nthinned = 0;
        ot->nexcluded = 0;
        ot->time_min = DBL_MAX;
        ot->time_max = -DBL_MAX;
    }

    for (i = 0; i < obs->nobs; ++i) {
        observation* o = &obs->data[i];
        obstype* ot = &obs->obstypes[o->type];

        ot->nobs++;
        if (o->status == STATUS_OK) {
            obs->ngood++;
            ot->ngood++;
        } else if (o->status == STATUS_OUTSIDEGRID) {
            obs->noutside_grid++;
            ot->noutside_grid++;
        } else if (o->status == STATUS_OUTSIDEOBSDOMAIN) {
            obs->noutside_obsdomain++;
            ot->noutside_obsdomain++;
        } else if (o->status == STATUS_OUTSIDEOBSWINDOW) {
            obs->noutside_obswindow++;
            ot->noutside_obswindow++;
        } else if (o->status == STATUS_LAND) {
            obs->nland++;
            ot->nland++;
        } else if (o->status == STATUS_SHALLOW) {
            obs->nshallow++;
            ot->nshallow++;
        } else if (o->status == STATUS_BADBATCH) {
            obs->nbadbatch++;
            ot->nbadbatch++;
        } else if (o->status == STATUS_BADFC) {
            obs->nbadfc++;
            ot->nbadfc++;
        } else if (o->status == STATUS_RANGE) {
            obs->nrange++;
            ot->nrange++;
        } else if (o->status == STATUS_THINNED) {
            obs->nthinned++;
            ot->nthinned++;
        } else if (o->status == STATUS_EXCLUDED) {
            obs->nexcluded++;
            ot->nexcluded++;
        }

        if (o->time < ot->time_min)
            ot->time_min = o->time;
        if (o->time > ot->time_max)
            ot->time_max = o->time;
    }
}

/** Reads observations from "observations.nc".
 */
void obs_read(observations* obs, char fname[])
{
    int ncid;
    double da_time = NAN;
    size_t nobs;

    ncw_open(fname, NC_NOWRITE, &ncid);

    if (ncw_att_exists(ncid, NC_GLOBAL, "DA_DAY"))
        ncw_get_att_double(ncid, NC_GLOBAL, "DA_DAY", &da_time);
    if (!enkf_noobsdatecheck && !isnan(da_time)) {
        obs->da_time = da_time;
        if (fabs(obs->da_time - da_time) > 1e-6)
            enkf_quit("observation data file \"%s\" from a different cycle", fname);
    }

    {
        int dimid_nobs[1];

        ncw_inq_dimid(ncid, "nobs", dimid_nobs);
        ncw_inq_dimlen(ncid, dimid_nobs[0], &nobs);
    }

    obs->nobs = nobs;
    enkf_printf("    %zu observations\n", nobs);
    if (nobs == 0) {
        obs->data = NULL;
        ncw_close(ncid);
        goto finish;
    }

    enkf_printf("    allocating %zu bytes for array of observations\n", nobs * sizeof(observation));
#if defined(USE_SHMEM)
    {
        MPI_Aint size;
        int ierror;

        assert(sizeof(MPI_Aint) == sizeof(size_t));
        size = nobs * sizeof(observation);
        ierror = MPI_Win_allocate_shared((sm_comm_rank == 0) ? size : 0, sizeof(observation), MPI_INFO_NULL, sm_comm, &obs->data, &obs->sm_comm_win_data);
        assert(ierror == MPI_SUCCESS);
        if (sm_comm_rank != 0) {
            int disp_unit;
            MPI_Aint my_size;

            ierror = MPI_Win_shared_query(obs->sm_comm_win_data, 0, &my_size, &disp_unit, &obs->data);
            assert(ierror == MPI_SUCCESS);
            assert(my_size == size);
            assert(disp_unit == sizeof(observation));
        } else
            memset(obs->data, 0, size);
        MPI_Win_fence(0, obs->sm_comm_win_data);
        MPI_Barrier(sm_comm);
    }
#else
    obs->data = calloc(nobs, sizeof(observation));
    assert(obs->data != NULL);
#endif

    /*
     * "type"
     */
    {
        int varid, natts, i, nobstypes_read;

        ncw_inq_varid(ncid, "type", &varid);
        ncw_inq_varnatts(ncid, varid, &natts);
        /*
         * check consistency of "type" attributes
         */
        for (i = 0, nobstypes_read = 0; i < natts; ++i) {
            char attname[NC_MAX_NAME];
            int typeid;

            ncw_inq_attname(ncid, varid, i, attname);
            typeid = obstype_getid(obs->nobstypes, obs->obstypes, attname, 0);
            if (typeid >= 0) {
                obstype* ot = &obs->obstypes[typeid];
                int typeid_read;

                ncw_check_attlen(ncid, varid, attname, 1);
                ncw_get_att_int(ncid, varid, attname, &typeid_read);
                assert(typeid == typeid_read);
                nobstypes_read++;

                if (ot->logapplied) {
                    char logattname[NC_MAX_NAME];
                    char logattval[MAXSTRLEN];

                    snprintf(logattname, NC_MAX_NAME, "%s:LOGAPPLIED", ot->name);
                    if (!ncw_att_exists(ncid, varid, logattname))
                        enkf_quit("%s: variable = \"type\": expected attribute \"%s\" to be present for a log-transformed obs. type", fname, logattname);
                    ncw_get_att_text(ncid, varid, logattname, logattval);
                    if (strncmp(logattval, "true", 4) != 0)
                        enkf_quit("%s: variable = \"type\": expected attribute \"%s\" to have value \"true\" a log-transformed obs. type", fname, logattname);
                }
            }
        }
        if (nobstypes_read != obs->nobstypes)
            enkf_quit("number of observation types in observations.nc = %d does not match that in the observation types parameter file = %d", nobstypes_read, obs->nobstypes);

        /*
         * "product"
         */
        ncw_inq_varid(ncid, "product", &varid);
        ncw_inq_varnatts(ncid, varid, &natts);
        for (i = 0; i < natts; ++i) {
            char name[NC_MAX_NAME];
            nc_type nctype;
            size_t len;

            ncw_inq_attname(ncid, varid, i, name);
            ncw_inq_att(ncid, varid, name, &nctype, &len);
            if (nctype == NC_INT && len == 1) {
                int productid;

                ncw_get_att_int(ncid, varid, name, &productid);
                st_add(obs->products, name, productid);
            }
        }

        /*
         * "instrument"
         */
        ncw_inq_varid(ncid, "instrument", &varid);
        ncw_inq_varnatts(ncid, varid, &natts);
        for (i = 0; i < natts; ++i) {
            char name[NC_MAX_NAME];
            nc_type nctype;
            size_t len;

            ncw_inq_attname(ncid, varid, i, name);
            ncw_inq_att(ncid, varid, name, &nctype, &len);
            if (nctype == NC_INT && len == 1) {
                int instid;

                ncw_get_att_int(ncid, varid, name, &instid);
                st_add(obs->instruments, name, instid);
            }
        }

        /*
         * fill datafile stringtable
         */
        ncw_inq_varid(ncid, "fid", &varid);
        ncw_inq_varnatts(ncid, varid, &natts);
        for (i = 0; i < natts; ++i) {
            char name[NC_MAX_NAME];
            char attstr[MAXSTRLEN];
            size_t len;
            int fileid;

            ncw_inq_attname(ncid, varid, i, name);
            if (!str2int(name, &fileid))
                continue;
            ncw_inq_attlen(ncid, varid, name, &len);
            assert(len < MAXSTRLEN);
            ncw_get_att_text(ncid, varid, name, attstr);
            attstr[len] = 0;
            st_add_ifabsent(obs->datafiles, attstr, fileid);
        }

        /*
         * get time units
         */
        ncw_inq_varid(ncid, "time", &varid);
        if (obs->datestr == NULL) {
            size_t len = 0;

            ncw_inq_attlen(ncid, varid, "units", &len);
            obs->datestr = malloc(len + 1);
            ncw_get_att_text(ncid, varid, "units", obs->datestr);
            obs->datestr[len] = 0;
        }
    }

#if defined(USE_SHMEM)
    if (sm_comm_rank == 0) {
#endif
        void* v = NULL;
        int structuredonly = ncw_var_exists(ncid, "fi");
        int varid;
        size_t i;

        if (structuredonly)
            v = malloc(nobs * ((sizeof(float) >= sizeof(int)) ? sizeof(float) : sizeof(int)));
        else
            v = malloc(nobs * sizeof(double));

        ncw_inq_varid(ncid, "id", &varid);
        ncw_get_var_int(ncid, varid, v);
        for (i = 0; i < nobs; ++i)
            obs->data[i].id = ((int*) v)[i];

        ncw_inq_varid(ncid, "id_orig", &varid);
        ncw_get_var_int(ncid, varid, v);
        for (i = 0; i < nobs; ++i)
            obs->data[i].id_orig = ((int*) v)[i];

        ncw_inq_varid(ncid, "status", &varid);
        ncw_get_var_short(ncid, varid, v);
        for (i = 0; i < nobs; ++i)
            obs->data[i].status = ((short int*) v)[i];

        ncw_inq_varid(ncid, "type", &varid);
        ncw_get_var_short(ncid, varid, v);
        for (i = 0; i < nobs; ++i)
            obs->data[i].type = ((short int*) v)[i];

        ncw_inq_varid(ncid, "product", &varid);
        ncw_get_var_short(ncid, varid, v);
        for (i = 0; i < nobs; ++i)
            obs->data[i].product = ((short int*) v)[i];

        ncw_inq_varid(ncid, "instrument", &varid);
        ncw_get_var_short(ncid, varid, v);
        for (i = 0; i < nobs; ++i)
            obs->data[i].instrument = ((short int*) v)[i];

        ncw_inq_varid(ncid, "fid", &varid);
        ncw_get_var_short(ncid, varid, v);
        for (i = 0; i < nobs; ++i)
            obs->data[i].fid = ((short int*) v)[i];

        ncw_inq_varid(ncid, "batch", &varid);
        ncw_get_var_int(ncid, varid, v);
        for (i = 0; i < nobs; ++i)
            obs->data[i].batch = ((int*) v)[i];

        ncw_inq_varid(ncid, "value", &varid);
        ncw_get_var_float(ncid, varid, v);
        for (i = 0; i < nobs; ++i)
            obs->data[i].value = ((float*) v)[i];

        if (!ncw_var_exists(ncid, "estd_orig"))
            ncw_inq_varid(ncid, "estd", &varid);
        else
            ncw_inq_varid(ncid, "estd_orig", &varid);
        ncw_get_var_float(ncid, varid, v);
        for (i = 0; i < nobs; ++i)
            obs->data[i].estd = ((float*) v)[i];

        if (ncw_var_exists(ncid, "footprint")) {
            obs->has_nonpointobs = 1;
            ncw_inq_varid(ncid, "footprint", &varid);
            ncw_get_var_float(ncid, varid, v);
            for (i = 0; i < nobs; ++i)
                obs->data[i].footprint = ((float*) v)[i];
        }

        ncw_inq_varid(ncid, "lon", &varid);
        ncw_get_var_float(ncid, varid, v);
        for (i = 0; i < nobs; ++i)
            obs->data[i].lon = ((float*) v)[i];

        ncw_inq_varid(ncid, "lat", &varid);
        ncw_get_var_float(ncid, varid, v);
        for (i = 0; i < nobs; ++i)
            obs->data[i].lat = ((float*) v)[i];

        ncw_inq_varid(ncid, "depth", &varid);
        ncw_get_var_float(ncid, varid, v);
        for (i = 0; i < nobs; ++i)
            obs->data[i].depth = ((float*) v)[i];

        ncw_inq_varid(ncid, "model_depth", &varid);
        ncw_get_var_float(ncid, varid, v);
        for (i = 0; i < nobs; ++i)
            obs->data[i].model_depth = ((float*) v)[i];

        if (structuredonly) {
            ncw_inq_varid(ncid, "fi", &varid);
            ncw_get_var_float(ncid, varid, v);
            for (i = 0; i < nobs; ++i)
                obs->data[i].fij[0] = ((float*) v)[i];

            ncw_inq_varid(ncid, "fj", &varid);
            ncw_get_var_float(ncid, varid, v);
            for (i = 0; i < nobs; ++i)
                obs->data[i].fij[1] = ((float*) v)[i];
        } else {
            ncw_inq_varid(ncid, "fi0", &varid);
            ncw_get_var_double(ncid, varid, v);
            for (i = 0; i < nobs; ++i)
                obs->data[i].fij[0] = ((double*) v)[i];

            ncw_inq_varid(ncid, "fi1", &varid);
            ncw_get_var_double(ncid, varid, v);
            for (i = 0; i < nobs; ++i)
                obs->data[i].fij[1] = ((double*) v)[i];

            ncw_inq_varid(ncid, "fi2", &varid);
            ncw_get_var_double(ncid, varid, v);
            for (i = 0; i < nobs; ++i)
                obs->data[i].fij[2] = ((double*) v)[i];
        }

        ncw_inq_varid(ncid, "fk", &varid);
        ncw_get_var_float(ncid, varid, v);
        for (i = 0; i < nobs; ++i)
            obs->data[i].fk = ((float*) v)[i];

        ncw_inq_varid(ncid, "time", &varid);
        ncw_get_var_float(ncid, varid, v);
        for (i = 0; i < nobs; ++i)
            obs->data[i].time = ((float*) v)[i];

        ncw_inq_varid(ncid, "aux", &varid);
        ncw_get_var_int(ncid, varid, v);
        for (i = 0; i < nobs; ++i)
            obs->data[i].aux = ((int*) v)[i];

        free(v);

        /*
         * if because of the roundup error time gets outside the allowed
         * range, then correct it
         */
        for (i = 0; i < nobs; ++i) {
            observation* o = &obs->data[i];
            obstype* ot = &obs->obstypes[o->type];

            if (o->time <= ot->obswindow_min)
                o->time = ot->obswindow_min + EPSF;
            else if (o->time >= ot->obswindow_max)
                o->time = ot->obswindow_max - EPSF;
        }
#if defined(USE_SHMEM)
    }
#endif

    ncw_close(ncid);

#if defined(MPI)
    MPI_Barrier(MPI_COMM_WORLD);
#endif

  finish:
    /*
     * clearing STATUS_BADBATCH is necessary for replicating the original run
     * of CALC for repeat runs
     */
    {
        size_t i;

        for (i = 0; i < nobs; ++i)
            if (obs->data[i].status == STATUS_BADBATCH)
                obs->data[i].status = STATUS_OK;
    }
    obs_calcstats(obs);
}

#if defined(ENKF_PREP)
/**
 */
void obs_write(observations* obs, char fname[])
{
    int nobs = obs->nobs;
    char tunits[MAXSTRLEN];
    int structuredonly = 1;

    int ncid;
    int dimid_nobs[1];
    int varid;
    void* v;
    int i;

    if (rank != 0)
        return;

    for (i = 0; i < obs->nobs; ++i)
        if (!isnan(obs->data[i].fij[2])) {
            structuredonly = 0;
            break;
        }

    if (file_exists(fname))
        enkf_quit("file \"%s\" already exists", fname);
    ncw_create(fname, NC_NOCLOBBER | obs->ncformat, &ncid);

    ncw_put_att_double(ncid, NC_GLOBAL, "DA_DAY", 1, &obs->da_time);

    ncw_def_dim(ncid, "nobs", nobs, dimid_nobs);
    /*
     * id
     */
    ncw_def_var(ncid, "id", NC_INT, 1, dimid_nobs, &varid);
    ncw_put_att_text(ncid, varid, "long_name", "observation ID");
    /*
     * id_orig
     */
    ncw_def_var(ncid, "id_orig", NC_INT, 1, dimid_nobs, &varid);
    ncw_put_att_text(ncid, varid, "long_name", "original observation ID");
    ncw_put_att_text(ncid, varid, "description", "for primary observations - the serial number of the primary observation during the reading of data files; for superobs - the original ID of the very first observation collated into this observation");
    /*
     * status
     */
    ncw_def_var(ncid, "status", NC_BYTE, 1, dimid_nobs, &varid);
    ncw_put_att_text(ncid, varid, "long_name", "observation status");
    i = STATUS_OK;
    ncw_put_att_int(ncid, varid, "STATUS_OK", 1, &i);
    i = STATUS_OUTSIDEGRID;
    ncw_put_att_int(ncid, varid, "STATUS_OUTSIDEGRID", 1, &i);
    i = STATUS_LAND;
    ncw_put_att_int(ncid, varid, "STATUS_LAND", 1, &i);
    i = STATUS_SHALLOW;
    ncw_put_att_int(ncid, varid, "STATUS_SHALLOW", 1, &i);
    i = STATUS_RANGE;
    ncw_put_att_int(ncid, varid, "STATUS_RANGE", 1, &i);
    i = STATUS_BADBATCH;
    ncw_put_att_int(ncid, varid, "STATUS_BADBATCH", 1, &i);
    i = STATUS_OUTSIDEOBSDOMAIN;
    ncw_put_att_int(ncid, varid, "STATUS_OUTSIDEOBSDOMAIN", 1, &i);
    i = STATUS_OUTSIDEOBSWINDOW;
    ncw_put_att_int(ncid, varid, "STATUS_OUTSIDEOBSWINDOW", 1, &i);
    i = STATUS_THINNED;
    ncw_put_att_int(ncid, varid, "STATUS_THINNED", 1, &i);
    i = STATUS_EXCLUDED;
    ncw_put_att_int(ncid, varid, "STATUS_EXCLUDED", 1, &i);
    i = STATUS_BADFC;
    ncw_put_att_int(ncid, varid, "STATUS_BADFC", 1, &i);
    /*
     * type
     */
    ncw_def_var(ncid, "type", NC_SHORT, 1, dimid_nobs, &varid);
    ncw_put_att_text(ncid, varid, "long_name", "observation type ID");
    for (i = 0; i < obs->nobstypes; ++i) {
        ncw_put_att_int(ncid, varid, obs->obstypes[i].name, 1, &i);
        if (obs->obstypes[i].logapplied) {
            char attname[NC_MAX_NAME];

            snprintf(attname, NC_MAX_NAME, "%s:LOGAPPLIED", obs->obstypes[i].name);
            ncw_put_att_text(ncid, varid, attname, "true");
        }
    }
    /*
     * product
     */
    ncw_def_var(ncid, "product", NC_SHORT, 1, dimid_nobs, &varid);
    ncw_put_att_text(ncid, varid, "long_name", "observation product ID");
    for (i = 0; i < st_getsize(obs->products); ++i)
        ncw_put_att_int(ncid, varid, st_findstringbyindex(obs->products, i), 1, &i);
    /*
     * instrument
     */
    ncw_def_var(ncid, "instrument", NC_SHORT, 1, dimid_nobs, &varid);
    ncw_put_att_text(ncid, varid, "long_name", "observation instrument ID");
    for (i = 0; i < st_getsize(obs->instruments); ++i)
        ncw_put_att_int(ncid, varid, st_findstringbyindex(obs->instruments, i), 1, &i);
    /*
     * fid
     */
    ncw_def_var(ncid, "fid", NC_SHORT, 1, dimid_nobs, &varid);
    ncw_put_att_text(ncid, varid, "long_name", "observation data file ID");
    for (i = 0; i < st_getsize(obs->datafiles); ++i) {
        char attname[NC_MAX_NAME];
        char* datafile = st_findstringbyindex(obs->datafiles, i);

        snprintf(attname, NC_MAX_NAME, "%d", i);
        ncw_put_att_text(ncid, varid, attname, datafile);
    }
    /*
     * batch
     */
    ncw_def_var(ncid, "batch", NC_INT, 1, dimid_nobs, &varid);
    ncw_put_att_text(ncid, varid, "long_name", "observation batch ID");
    /*
     * vallue
     */
    ncw_def_var(ncid, "value", NC_FLOAT, 1, dimid_nobs, &varid);
    ncw_put_att_text(ncid, varid, "long_name", "observation value");
    /*
     * estd
     */
    ncw_def_var(ncid, "estd", NC_FLOAT, 1, dimid_nobs, &varid);
    ncw_put_att_text(ncid, varid, "long_name", "standard deviation of observation error used in DA");
    /*
     * footprint
     */
    if (obs->has_nonpointobs) {
        ncw_def_var(ncid, "footprint", NC_FLOAT, 1, dimid_nobs, &varid);
        ncw_put_att_text(ncid, varid, "long_name", "observation footprint in km");
    }
    /*
     * lon
     */
    ncw_def_var(ncid, "lon", NC_FLOAT, 1, dimid_nobs, &varid);
    ncw_put_att_text(ncid, varid, "long_name", "observation longitude");
    /*
     * lat
     */
    ncw_def_var(ncid, "lat", NC_FLOAT, 1, dimid_nobs, &varid);
    ncw_put_att_text(ncid, varid, "long_name", "observation latitude");
    /*
     * depth
     */
    ncw_def_var(ncid, "depth", NC_FLOAT, 1, dimid_nobs, &varid);
    ncw_put_att_text(ncid, varid, "long_name", "observation depth/height");
    /*
     * model_depth
     */
    ncw_def_var(ncid, "model_depth", NC_FLOAT, 1, dimid_nobs, &varid);
    ncw_put_att_text(ncid, varid, "long_name", "model bottom depth at the observation location");
    if (structuredonly) {
        /*
         * fi
         */
        ncw_def_var(ncid, "fi", NC_FLOAT, 1, dimid_nobs, &varid);
        ncw_put_att_text(ncid, varid, "long_name", "fractional grid index i of the observation");
        /*
         * fj
         */
        ncw_def_var(ncid, "fj", NC_FLOAT, 1, dimid_nobs, &varid);
        ncw_put_att_text(ncid, varid, "long_name", "fractional grid index j of the observation");
    } else {
        /*
         * fi0
         */
        ncw_def_var(ncid, "fi0", NC_DOUBLE, 1, dimid_nobs, &varid);
        ncw_put_att_text(ncid, varid, "long_name", "fractional grid index i0 (i) of the observation");
        /*
         * fi1
         */
        ncw_def_var(ncid, "fi1", NC_DOUBLE, 1, dimid_nobs, &varid);
        ncw_put_att_text(ncid, varid, "long_name", "fractional grid index i1 (j) of the observation");
        /*
         * fi2
         */
        ncw_def_var(ncid, "fi2", NC_DOUBLE, 1, dimid_nobs, &varid);
        ncw_put_att_text(ncid, varid, "long_name", "fractional grid index i2 of the observation");
    }
    /*
     * fk
     */
    ncw_def_var(ncid, "fk", NC_FLOAT, 1, dimid_nobs, &varid);
    ncw_put_att_text(ncid, varid, "long_name", "fractional grid index k of the observation");
    /*
     * time
     */
    ncw_def_var(ncid, "time", NC_FLOAT, 1, dimid_nobs, &varid);
    ncw_put_att_text(ncid, varid, "long_name", "observation time");
    snprintf(tunits, MAXSTRLEN, "days from %s", obs->datestr);
    ncw_put_att_text(ncid, varid, "units", tunits);
    /*
     * aux
     */
    ncw_def_var(ncid, "aux", NC_INT, 1, dimid_nobs, &varid);
    ncw_put_att_text(ncid, varid, "long_name", "auxiliary information");
    ncw_put_att_text(ncid, varid, "description", "for primary observations - the ID of the superobservation it is collated into; for superobservations - the number of primary observations collated");

    if (obs->nccompression > 0)
        ncw_def_deflate(ncid, 0, 1, obs->nccompression);
    ncw_enddef(ncid);

    if (nobs == 0) {
        ncw_close(ncid);
        return;
    }

    for (i = 0; i < obs->nobs; ++i) {
        observation* o = &obs->data[i];

        if (o->status == STATUS_OK && !isfinite(o->value)) {
            obs_printob(obs, i);
            enkf_quit("obs_write(): bad value");
        }
    }

    if (structuredonly)
        v = malloc(nobs * ((sizeof(float) >= sizeof(int)) ? sizeof(float) : sizeof(int)));
    else
        v = malloc(nobs * sizeof(double));

    ncw_inq_varid(ncid, "id", &varid);
    for (i = 0; i < obs->nobs; ++i)
        ((int*) v)[i] = obs->data[i].id;
    ncw_put_var_int(ncid, varid, v);

    ncw_inq_varid(ncid, "id_orig", &varid);
    for (i = 0; i < obs->nobs; ++i)
        ((int*) v)[i] = obs->data[i].id_orig;
    ncw_put_var_int(ncid, varid, v);

    ncw_inq_varid(ncid, "status", &varid);
    for (i = 0; i < obs->nobs; ++i)
        ((unsigned char*) v)[i] = obs->data[i].status;
    ncw_put_var_uchar(ncid, varid, v);

    ncw_inq_varid(ncid, "type", &varid);
    for (i = 0; i < obs->nobs; ++i)
        ((short int*) v)[i] = obs->data[i].type;
    ncw_put_var_short(ncid, varid, v);

    ncw_inq_varid(ncid, "product", &varid);
    for (i = 0; i < obs->nobs; ++i)
        ((short int*) v)[i] = obs->data[i].product;
    ncw_put_var_short(ncid, varid, v);

    ncw_inq_varid(ncid, "instrument", &varid);
    for (i = 0; i < obs->nobs; ++i)
        ((short int*) v)[i] = obs->data[i].instrument;
    ncw_put_var_short(ncid, varid, v);

    ncw_inq_varid(ncid, "fid", &varid);
    for (i = 0; i < obs->nobs; ++i)
        ((short int*) v)[i] = obs->data[i].fid;
    ncw_put_var_short(ncid, varid, v);

    ncw_inq_varid(ncid, "batch", &varid);
    for (i = 0; i < obs->nobs; ++i)
        ((int*) v)[i] = obs->data[i].batch;
    ncw_put_var_int(ncid, varid, v);

    ncw_inq_varid(ncid, "value", &varid);
    for (i = 0; i < obs->nobs; ++i)
        ((float*) v)[i] = obs->data[i].value;
    ncw_put_var_float(ncid, varid, v);

    ncw_inq_varid(ncid, "estd", &varid);
    for (i = 0; i < obs->nobs; ++i)
        ((float*) v)[i] = obs->data[i].estd;
    ncw_put_var_float(ncid, varid, v);

    if (obs->has_nonpointobs) {
        ncw_inq_varid(ncid, "footprint", &varid);
        for (i = 0; i < obs->nobs; ++i)
            ((float*) v)[i] = obs->data[i].footprint;
        ncw_put_var_float(ncid, varid, v);
    }

    ncw_inq_varid(ncid, "lon", &varid);
    for (i = 0; i < obs->nobs; ++i)
        ((float*) v)[i] = obs->data[i].lon;
    ncw_put_var_float(ncid, varid, v);

    ncw_inq_varid(ncid, "lat", &varid);
    for (i = 0; i < obs->nobs; ++i)
        ((float*) v)[i] = obs->data[i].lat;
    ncw_put_var_float(ncid, varid, v);

    ncw_inq_varid(ncid, "depth", &varid);
    for (i = 0; i < obs->nobs; ++i)
        ((float*) v)[i] = obs->data[i].depth;
    ncw_put_var_float(ncid, varid, v);

    ncw_inq_varid(ncid, "model_depth", &varid);
    for (i = 0; i < obs->nobs; ++i)
        ((float*) v)[i] = obs->data[i].model_depth;
    ncw_put_var_float(ncid, varid, v);

    if (structuredonly) {
        ncw_inq_varid(ncid, "fi", &varid);
        for (i = 0; i < obs->nobs; ++i)
            ((float*) v)[i] = obs->data[i].fij[0];
        ncw_put_var_float(ncid, varid, v);

        ncw_inq_varid(ncid, "fj", &varid);
        for (i = 0; i < obs->nobs; ++i)
            ((float*) v)[i] = obs->data[i].fij[1];
        ncw_put_var_float(ncid, varid, v);
    } else {
        ncw_inq_varid(ncid, "fi0", &varid);
        for (i = 0; i < obs->nobs; ++i)
            ((double*) v)[i] = obs->data[i].fij[0];
        ncw_put_var_double(ncid, varid, v);

        ncw_inq_varid(ncid, "fi1", &varid);
        for (i = 0; i < obs->nobs; ++i)
            ((double*) v)[i] = obs->data[i].fij[1];
        ncw_put_var_double(ncid, varid, v);

        ncw_inq_varid(ncid, "fi2", &varid);
        for (i = 0; i < obs->nobs; ++i)
            ((double*) v)[i] = obs->data[i].fij[2];
        ncw_put_var_double(ncid, varid, v);
    }

    ncw_inq_varid(ncid, "fk", &varid);
    for (i = 0; i < obs->nobs; ++i)
        ((float*) v)[i] = obs->data[i].fk;
    ncw_put_var_float(ncid, varid, v);

    ncw_inq_varid(ncid, "time", &varid);
    for (i = 0; i < obs->nobs; ++i)
        ((float*) v)[i] = obs->data[i].time;
    ncw_put_var_float(ncid, varid, v);

    ncw_inq_varid(ncid, "aux", &varid);
    for (i = 0; i < obs->nobs; ++i)
        ((int*) v)[i] = obs->data[i].aux;
    ncw_put_var_int(ncid, varid, v);

    ncw_close(ncid);
    free(v);
}
#endif

#if defined(ENKF_PREP)
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
    ncw_inq_varid(ncid, "aux", &varid_aux);
    ncw_inq_dimid(ncid, "nobs", &dimid_nobs);
    ncw_inq_dimlen(ncid, dimid_nobs, &n);
    assert((int) n == obs->nobs);

    aux = malloc(n * sizeof(int));
    for (i = 0, ii = 0; i < obs->nobs; ++i) {
        aux[ii] = obs->data[i].aux;
        ii++;
    }

    ncw_put_var_int(ncid, varid_aux, aux);
    ncw_close(ncid);
    free(aux);
}
#endif

#if defined(ENKF_CALC)
/**
 */
int obs_modifiederrors_alreadywritten(observations* obs, char fname[])
{
    int ncid;
    int iswritten;

    ncw_open(fname, NC_NOWRITE, &ncid);
    iswritten = ncw_var_exists(ncid, "estd_orig");
    ncw_close(ncid);

    return iswritten;
}
#endif

#if defined(ENKF_PREP)
/**
 */
static int cmp_xyz(const void* p1, const void* p2)
{
    observation* o1 = (observation*) p1;
    observation* o2 = (observation*) p2;

    if (o1->lon > o2->lon)
        return 1;
    else if (o1->lon < o2->lon)
        return -1;
    else if (o1->lat > o2->lat)
        return 1;
    else if (o1->lat < o2->lat)
        return -1;
    else if (o1->depth > o2->depth)
        return 1;
    else if (o1->depth < o2->depth)
        return -1;
    else if (o1->instrument > o2->instrument)
        return 1;
    else if (o1->instrument < o2->instrument)
        return -1;
    return 0;
}

/**
 */
static int cmp_xy(const void* p1, const void* p2)
{
    observation* o1 = (observation*) p1;
    observation* o2 = (observation*) p2;

    if (o1->lon > o2->lon)
        return 1;
    else if (o1->lon < o2->lon)
        return -1;
    else if (o1->lat > o2->lat)
        return 1;
    else if (o1->lat < o2->lat)
        return -1;
    else if (o1->instrument > o2->instrument)
        return 1;
    else if (o1->instrument < o2->instrument)
        return -1;
    return 0;
}

/**
 */
void obs_superob(observations* obs, __compar_d_fn_t cmp_obs, observations** sobs, int sobid, int do_thin)
{
    int i1 = 0, i2 = 0;
    int nsobs = 0;
    int nthinned = 0;
    observation* data = obs->data;
    observation* sdata = NULL;
    int has_nonpointobs = 0;

    while (i2 < obs->ngood) {
        observation* so;
        observation* o;
        double lon_min, lon_max;
        int ii;
        double sum, sumsq, subvar;
        double evar;

        /*
         * identify obs that will be combined into this superob 
         */
        while (i2 + 1 < obs->ngood && cmp_obs(&data[i1], &data[i2 + 1], obs) == 0)
            i2++;

        /*
         * Thin observations with identical positions (these are supposedly
         * coming from high-frequency instruments)
         */
        if (i2 > i1 && do_thin && obs->obstypes[data[i1].type].can_thin) {
            int i11 = i1;
            int i22 = i1;
            int (*cmp) (const void*, const void*) = (do_thin == THIN_XYZ) ? cmp_xyz : cmp_xy;

            qsort(&data[i1], i2 - i1 + 1, sizeof(observation), cmp);
            while (i22 <= i2) {
                while (i22 + 1 <= i2 && cmp(&data[i11], &data[i22 + 1]) == 0)
                    i22++;
                /*
                 * replace the value of the kept observation by the average
                 * of the batch
                 */
                if (i22 > i11) {
                    if (do_thin == THIN_XYZ) {
                        float sum_value = data[i11].value;

                        for (ii = i11 + 1; ii <= i22; ++ii) {
                            data[ii].status = STATUS_THINNED;
                            sum_value += data[ii].value;
                            nthinned++;
                        }
                        data[i11].value = sum_value / (float) (i22 - i11 + 1);
                    } else if (do_thin == THIN_XY) {
                        float sum_value = data[i11].value;
                        float sum_depth = data[i11].depth;
                        float sum_fk = data[i11].fk;

                        for (ii = i11 + 1; ii <= i22; ++ii) {
                            data[ii].status = STATUS_THINNED;
                            sum_value += data[ii].value;
                            sum_depth += data[ii].depth;
                            sum_fk += data[ii].fk;
                            nthinned++;
                        }
                        data[i11].value = sum_value / (float) (i22 - i11 + 1);
                        data[i11].depth = sum_depth / (float) (i22 - i11 + 1);
                        data[i11].fk = sum_fk / (float) (i22 - i11 + 1);
                    } else
                        enkf_quit("programming error");
                }
                i22++;
                i11 = i22;
            }
        }

        if (nsobs % NOBS_INC == 0)
            sdata = realloc(sdata, (nsobs + NOBS_INC) * sizeof(observation));

        if (sobid == nsobs) {
            int i;

            enkf_printf("    sob # %d is formed by the following obs:\n", sobid);
            for (i = i1; i <= i2; ++i) {
                enkf_printf("      ");
                obs_printob(obs, i);
            }
        }

        if (!enkf_considersubgridvar)
            subvar = 0.0;
        else {
            /*
             * Calculate subgrid variance. For now assume equal weights.
             */
            sum = 0.0;
            sumsq = 0.0;
            for (ii = i1; ii <= i2; ++ii) {
                sum += data[ii].value;
                sumsq += data[ii].value * data[ii].value;
            }
            subvar = (i2 - i1 > 0) ? (sumsq - sum * sum / (double) (i2 - i1 + 1)) / (double) (i2 - i1) : 0.0;
        }

        /*
         * set the "aux" field in the original obs to the id of the superob
         * they are merged into
         */
        for (ii = i1; ii <= i2; ++ii)
            data[ii].aux = nsobs;

        o = &data[i1];
        assert(o->status == STATUS_OK);
        so = &sdata[nsobs];

        so->type = o->type;
        so->product = o->product;
        so->instrument = o->instrument;
        so->id = nsobs;
        so->id_orig = o->id_orig;
        so->fid = o->fid;
        so->batch = o->batch;

        evar = o->estd * o->estd;
        if (subvar > evar) {
            evar = subvar;
            obs->obstypes[o->type].nsubgrid++;
        }
        so->estd = 1.0 / evar;
        so->value = o->value;
        so->footprint = o->footprint;
        so->lon = o->lon;
        so->lat = o->lat;
        so->depth = o->depth;
        so->model_depth = o->model_depth;
        so->fij[0] = o->fij[0];
        so->fij[1] = o->fij[1];
        so->fij[2] = o->fij[2];
        so->fk = o->fk;
        so->time = o->time;
        so->status = o->status;
        so->aux = 1;

        lon_min = o->lon;
        lon_max = o->lon;
        for (ii = i1 + 1; ii <= i2; ++ii) {
            o = &data[ii];
            if (o->status == STATUS_THINNED)
                continue;
            if (so->product != o->product)
                so->product = -1;
            if (so->instrument != o->instrument)
                so->instrument = -1;
            if (so->fid != o->fid)
                so->fid = -1;
            if (so->batch != o->batch)
                so->batch = -1;

            /*
             * find average of collated obs weighted by inverse error variance
             */
            evar = o->estd * o->estd;
            evar = (subvar > evar) ? subvar : evar;
            so->value = so->value * so->estd + o->value / evar;
            so->lon = so->lon * so->estd + o->lon / evar;
            so->lat = so->lat * so->estd + o->lat / evar;
            so->depth = so->depth * so->estd + o->depth / evar;
            so->model_depth = so->model_depth * so->estd + o->model_depth / evar;
            so->fij[0] = so->fij[0] * so->estd + o->fij[0] / evar;
            so->fij[1] = so->fij[1] * so->estd + o->fij[1] / evar;
            so->fij[2] = so->fij[2] * so->estd + o->fij[2] / evar;
            so->fk = so->fk * so->estd + o->fk / evar;
            so->time = so->time * so->estd + o->time / evar;

            so->estd += 1.0 / evar;

            so->value /= so->estd;
            so->lon /= so->estd;
            so->lat /= so->estd;
            so->depth /= so->estd;
            so->model_depth /= so->estd;
            so->fij[0] /= so->estd;
            so->fij[1] /= so->estd;
            so->fij[2] /= so->estd;
            so->fk /= so->estd;
            so->time /= so->estd;
            so->aux++;

            if (o->lon < lon_min)
                lon_min = o->lon;
            if (o->lon > lon_max)
                lon_max = o->lon;

            assert(o->status == STATUS_OK);
        }

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
                evar = o->estd * o->estd;
                evar = (subvar > evar) ? subvar : evar;
                so->lon += o->lon / evar;
            }
            so->lon *= so->estd;
            if (so->lon < 0.0)
                so->lon += 360.0;
        }
        so->estd = sqrt(1.0 / so->estd);
        if (so->estd < obs->obstypes[so->type].estdmin)
            so->estd = obs->obstypes[so->type].estdmin;
        if (so->footprint > 0.0)
            has_nonpointobs = 1;

        nsobs++;

        i1 = i2 + 1;
        i2 = i1;
    }                           /* main cycle */
    if (nthinned > 0) {
        enkf_printf("    %d observations thinned\n", nthinned);
        obs_calcstats(obs);
    }

    enkf_printf("    %d superobservations\n", nsobs);

    *sobs = obs_create_fromdata(obs, nsobs, sdata);
    if (sobid >= 0) {
        enkf_printf("    sob # %d info:\n", sobid);
        enkf_printf("      ");
        obs_printob(*sobs, sobid);
    }
    (*sobs)->has_nonpointobs = has_nonpointobs;
}
#endif

#if defined (ENKF_CALC)
/**
 */
void obs_find_bytype(observations* obs, int type, int* nobs, int** obsids)
{
    int i;

    /*
     * it is likely that this check can be safely removed
     */
    if (!enkf_fstatsonly)
        assert(obs->obstypes[type].nobs == obs->obstypes[type].ngood);

    *nobs = 0;
    if (obs->obstypes[type].nobs == 0) {
        *obsids = NULL;
        return;
    }
    *obsids = malloc(obs->obstypes[type].nobs * sizeof(int));
    for (i = 0; i < obs->nobs; ++i) {
        observation* o = &obs->data[i];

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
#endif

#if defined (ENKF_CALC)
/**
 */
void obs_find_bytypeandtime(observations* obs, int type, int time, int* nobs, int** obsids)
{
    obstype* ot = &obs->obstypes[type];
    int i;

    /*
     * it is likely that this check can be safely removed
     */
    if (!enkf_fstatsonly)
        assert(ot->nobs == ot->ngood);

    *nobs = 0;
    if (ot->ngood == 0) {
        *obsids = NULL;
        return;
    }
    *obsids = malloc(obs->obstypes[type].ngood * sizeof(int));
    for (i = 0; i < obs->nobs; ++i) {
        observation* o = &obs->data[i];

        if (o->type == type && o->status == STATUS_OK && get_tshift(o->time, ot->async_tstep, ot->async_centred) == time) {
            (*obsids)[*nobs] = i;
            (*nobs)++;
        }
    }
    if (*nobs == 0) {
        free(*obsids);
        *obsids = NULL;
    }
}
#endif

/**
 */
void obs_printob(observations* obs, int i)
{
    observation* o = &obs->data[i];

    enkf_printf("type = %s, product = %s, instrument = %s, datafile = %s, id = %d, original id = %d, batch = %d, value = %.3g, estd = %.3g, footprint = %.3g, lon = %.3f, lat = %.3f, depth = %.1f, model_depth = %.1f, ", obs->obstypes[o->type].name, st_findstringbyindex(obs->products, o->product), st_findstringbyindex(obs->instruments, o->instrument), st_findstringbyindex(obs->datafiles, o->fid), o->id, o->id_orig, (int) o->batch, o->value, o->estd, o->footprint, o->lon, o->lat, o->depth, o->model_depth);
    if (!isfinite(o->fij[2]))
        enkf_printf("fi = %.3f, fj = %.3f, ", o->fij[0], o->fij[1]);
    else
        enkf_printf("fi0 = %.3f, fi1 = %.3f, fi2 = %.3f, ", o->fij[0], o->fij[1], o->fij[2]);
    enkf_printf("fk = %.3f, time = %.3g, status = %d\n", o->fk, o->time, o->status);
}

#if defined(ENKF_CALC)
/** Create kd-trees for (3D locations of) observations of each type. If
 ** (USE_SHMEM) then put them on core #0 of each compute node and access from
 ** other cores using shared memory machinery of MPI-3.
 */
void obs_createkdtrees(observations* obs)
{
    int otid;

    assert(obs->loctrees == NULL && obs->obsids == NULL);
    obs->loctrees = calloc(obs->nobstypes, sizeof(kdtree*));
    obs->obsids = calloc(obs->nobstypes, sizeof(int*));
#if defined(USE_SHMEM)
    assert(obs->sm_comm_wins_kd == NULL);
    obs->sm_comm_wins_kd = calloc(obs->nobstypes, sizeof(MPI_Win));
    for (otid = 0; otid < obs->nobstypes; ++otid)
        obs->sm_comm_wins_kd[otid] = MPI_WIN_NULL;
#endif

    for (otid = 0; otid < obs->nobstypes; ++otid) {
        obstype* ot = &obs->obstypes[otid];
        kdtree** tree = &obs->loctrees[otid];
        int nobs = 0;
        int* obsids = NULL;

#if defined(OBS_SHUFFLE)
        size_t* ids = NULL;
#endif
        int i;

        if (ot->statsonly)
            continue;

        obs_find_bytype(obs, otid, &nobs, &obsids);
        if (nobs == 0)
            continue;
        ot->nobs = nobs;
        if (obs->obsids[otid] != NULL)
            free(obs->obsids[otid]);
        obs->obsids[otid] = obsids;

        assert(*tree == NULL);

#if defined(OBS_SHUFFLE)
        ids = malloc(nobs * sizeof(size_t));
        for (i = 0; i < nobs; ++i)
            ids[i] = i;
        shuffle(nobs, ids);
#endif

        *tree = kd_create(ot->name, 3);
#if defined(USE_SHMEM)
        {
            MPI_Win* sm_comm_win = &obs->sm_comm_wins_kd[otid];
            MPI_Aint size;
            int ierror;

            size = kd_getstoragesize(*tree, nobs);
            if (sm_comm_rank == 0) {
                void* storage = NULL;

                assert(sizeof(MPI_Aint) == sizeof(size_t));
                ierror = MPI_Win_allocate_shared(size, sizeof(double), MPI_INFO_NULL, sm_comm, &storage, sm_comm_win);
                assert(ierror == MPI_SUCCESS);
                kd_setstorage(*tree, nobs, storage, 1);
            } else {
                MPI_Aint my_size;
                void* storage = NULL;
                int disp_unit;

                ierror = MPI_Win_allocate_shared(0, sizeof(double), MPI_INFO_NULL, sm_comm, &storage, sm_comm_win);
                assert(ierror == MPI_SUCCESS);
                ierror = MPI_Win_shared_query(*sm_comm_win, 0, &my_size, &disp_unit, &storage);
                assert(ierror == MPI_SUCCESS);
                assert(my_size = size);
                kd_setstorage(*tree, nobs, storage, 0);
            }
        }
#else
        kd_setstorage(*tree, nobs, NULL, -1);
#endif
#if defined(USE_SHMEM)
        if (sm_comm_rank == 0) {
#endif
            for (i = 0; i < nobs; ++i) {
#if defined(OBS_SHUFFLE)
                int id = ids[i];
                observation* o = &obs->data[obsids[id]];
#else
                observation* o = &obs->data[obsids[i]];
#endif
                double ll[2] = { o->lon, o->lat };
                double xyz[3];

                ll2xyz(ll, xyz);
#if defined(OBS_SHUFFLE)
                kd_insertnode(*tree, xyz, id);
#else
                kd_insertnode(*tree, xyz, i);
#endif
            }
#if defined(USE_SHMEM)
        } else
            kd_syncsize(*tree);
        MPI_Win_fence(0, obs->sm_comm_wins_kd[otid]);
        MPI_Barrier(sm_comm);
#endif
        kd_printinfo(*tree, "      ");

#if defined(OBS_SHUFFLE)
        free(ids);
#endif
    }
}
#endif

#if defined(ENKF_CALC)
/**
 */
void obs_destroykdtrees(observations* obs)
{
    int otid;

    if (obs->loctrees == NULL)
        return;

    for (otid = 0; otid < obs->nobstypes; ++otid) {
        kdtree** tree = &obs->loctrees[otid];

        if (*tree != NULL)
            kd_destroy(*tree);
        *tree = NULL;
    }
#if defined(USE_SHMEM)
    if (obs->sm_comm_wins_kd != NULL) {
        for (otid = 0; otid < obs->nobstypes; ++otid) {
            if (obs->sm_comm_wins_kd[otid] == MPI_WIN_NULL)
                continue;
            MPI_Win_free(&obs->sm_comm_wins_kd[otid]);
            assert(obs->sm_comm_wins_kd[otid] == MPI_WIN_NULL);
        }
        free(obs->sm_comm_wins_kd);
        obs->sm_comm_wins_kd = NULL;
    }
#endif
    free(obs->loctrees);
    obs->loctrees = NULL;
}
#endif

#if defined(ENKF_CALC)
/**
 */
static void obs_getglobal(observations* obs, char* domainname, int* n, int** ids, double** lcoeffs, int* ploc_allocated)
{
    int otid;

    for (otid = 0; otid < obs->nobstypes; ++otid) {
        obstype* ot = &obs->obstypes[otid];
        int* obsids = obs->obsids[otid];
        int i;

        if (ot->nobs == 0 || ot->statsonly)
            continue;
        /*
         * check whether observations of this type are visible from the
         * specified domain
         */
        if (domainname != NULL && ot->ndomains > 0) {
            /*
             * (if ot->ndomains = 0 then observations of this type are visible
             * from all grids)
             */
            int d;

            for (d = 0; d < ot->ndomains; ++d)
                if (strcasecmp(domainname, ot->domainnames[d]) == 0)
                    break;
            if (d == ot->ndomains)
                continue;
        }
        if (ploc_allocated != NULL) {
            if (*n + ot->nobs >= *ploc_allocated) {
                *ploc_allocated = *n + ot->nobs;
                *ids = realloc(*ids, *ploc_allocated * sizeof(int));
                *lcoeffs = realloc(*lcoeffs, *ploc_allocated * sizeof(double));
            }
        } else {
            *ids = realloc(*ids, (*n + ot->nobs) * sizeof(int));
            *lcoeffs = realloc(*lcoeffs, (*n + ot->nobs) * sizeof(double));
        }

        for (i = 0; i < ot->nobs; ++i, ++(*n)) {
            (*ids)[*n] = obsids[i];
            (*lcoeffs)[*n] = 1.0;
        }
    }
}

/** For each observation type find observations within localisation radius from
 ** the specified location (lon,lat) and calculate the corresponding taper
 ** coefficients. If there is a limit on the number of local observations then
 ** sort observations according to distance and keep the specified number of the
 ** closest observations. If *ploc_allocated is NULL then allocate arrays of
 ** observation ids and taper coefficients, if not -- then assume that these
 ** arrays are pre-allocated.
 * @param obs - observations
 * @param lon - longitude
 * @param lat - latitude
 * @param domainname - the domain
 * @param n - number of local obs. found
 * @param ids - array of the local obs. ids
 * @param lcoeffs - array of taper coefficients
 * @ploc_allocated - (pointer to) the number of local observations the output
 *                   arrays have been pre-allocated for; not preallocated if
 *                   NULL
 */
void obs_findlocal(observations* obs, double lon, double lat, char* domainname, int* n, int** ids, double** lcoeffs, int* ploc_allocated)
{
    double ll[2] = { lon, lat };
    double xyz[3];
    int otid;
    int i;

    if (obs->nobstypes == 0)
        return;

    /*
     * This is supposed to happen if horizontal grid has type GRIDHTYPE_NONE.
     * The treatment could be very simple (just adding all obs), but it is
     * still necessary to check each obs. type if it belongs to the grid's
     * domain.
     */
    if (isnan(lon + lat)) {
        obs_getglobal(obs, domainname, n, ids, lcoeffs, ploc_allocated);
        return;
    }

    ll2xyz(ll, xyz);

    if (obs->loctrees == NULL)
        obs_createkdtrees(obs);

    for (otid = 0, i = 0; otid < obs->nobstypes; ++otid) {
        obstype* ot = &obs->obstypes[otid];
        kdtree* tree = obs->loctrees[otid];
        int* obsids = obs->obsids[otid];
        size_t nloc;
        kdresult* results;
        int iloc;

        if (ot->nobs == 0 || ot->statsonly)
            continue;

        /*
         * check whether observations of this type are visible from the
         * specified domain
         */
        if (domainname != NULL && ot->ndomains > 0) {
            /*
             * (if ot->ndomains = 0 then observations of this type are visible
             * from all grids)
             */
            int d;

            for (d = 0; d < ot->ndomains; ++d)
                if (strcasecmp(domainname, ot->domainnames[d]) == 0)
                    break;
            if (d == ot->ndomains)
                continue;
        }

        kd_findnodeswithinrange(tree, xyz, obstype_getmaxlocrad(ot), (ot->nlobsmax == INT_MAX) ? 0 : 1, &nloc, &results);
        if (nloc > ot->nlobsmax)
            nloc = ot->nlobsmax;
        for (iloc = 0; iloc < nloc; ++iloc) {
            size_t id_orig = kd_getnodedata(tree, results[iloc].id);
            observation* o = &obs->data[obsids[id_orig]];

            if (o->status != STATUS_OK)
                continue;

            if (ploc_allocated != NULL) {
                if (i >= *ploc_allocated) {
                    *ploc_allocated += PLOC_INC;
                    *ids = realloc(*ids, *ploc_allocated * sizeof(int));
                    *lcoeffs = realloc(*lcoeffs, *ploc_allocated * sizeof(double));
                }
            } else {
                if (i % PLOC_INC == 0) {
                    *ids = realloc(*ids, (i + PLOC_INC) * sizeof(int));
                    *lcoeffs = realloc(*lcoeffs, (i + PLOC_INC) * sizeof(double));
                }
            }
            (*ids)[i] = obsids[id_orig];
            (*lcoeffs)[i] = obstype_calclcoeff(ot, sqrt(results[iloc].distsq));
            i++;
        }
    }
    *n = i;
}
#endif

#if defined(ENKF_CALC)
#define BYTE_PER_SHORT 2
#define BYTE_PER_INT 4
/**
 */
void obs_markbadbatches(observations* obs, hashtable* badbatches)
{
    int i;

    assert(BYTE_PER_SHORT == sizeof(short));
    assert(BYTE_PER_INT == sizeof(int));

    if (badbatches == NULL)
        return;

    for (i = 0; i < obs->nobs; ++i) {
        observation* o = &obs->data[i];
        int key_int[2] = { -1, -1 };
        short* key_short = (short*) key_int;

        key_int[0] = o->batch;
        key_short[2] = o->type;
        key_short[3] = o->fid;

        if (ht_find(badbatches, key_int) != NULL)
            o->status = STATUS_BADBATCH;
    }
    obs_calcstats(obs);
}
#endif

#if defined(ENKF_CALC)
/**
 */
void obs_writeobsstatus(observations* obs, char fname[])
{
    int ncid;
    int dimid_nobs[1];
    size_t nobs;
    int varid;
    unsigned char* v;
    int o;

    if (obs->nobs == 0)
        return;
    if (rank != 0)
        return;

    ncw_open(fname, NC_WRITE, &ncid);
    ncw_inq_dimid(ncid, "nobs", dimid_nobs);
    ncw_inq_dimlen(ncid, dimid_nobs[0], &nobs);
    assert(nobs == obs->nobs);

    v = malloc(nobs);
    for (o = 0; o < (int) nobs; ++o)
        v[o] = obs->data[o].status;
    ncw_inq_varid(ncid, "status", &varid);
    ncw_put_var_uchar(ncid, varid, v);
    free(v);

    ncw_close(ncid);
}
#endif
