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
#include <assert.h>
#include <stdint.h>
#include "ncw.h"
#include "definitions.h"
#include "utils.h"
#include "observations.h"

#define NOBSTYPES_INC 10
#define PLOC_INC 10000
#define HT_SIZE 500
#define EPSD 1.0e-10

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

/**
 */
observations* obs_create_fromprm(enkfprm* prm)
{
    observations* obs = obs_create();

    obs->products = st_create("products");
    obs->instruments = st_create("instruments");
    obs->datafiles = st_create("datafiles");

    enkf_printf("  reading observation type specs from \"%s\":\n", prm->obstypeprm);
    obstypes_read(prm, prm->obstypeprm, &obs->nobstypes, &obs->obstypes);

#if defined(ENKF_PREP)
    obs->da_time = date2day(prm->date);
    obs->datestr = strdup(prm->date);

    if (file_exists(FNAME_BADBATCHES)) {
        FILE* f = NULL;
        char buf[MAXSTRLEN];
        int line;

        enkf_printf("  reading bad batches:\n");

        obs->badbatches = ht_create_i1s2(HT_SIZE);
        f = enkf_fopen(FNAME_BADBATCHES, "r");
        line = 0;
        while (fgets(buf, MAXSTRLEN, f) != NULL) {
            char otname[MAXSTRLEN];
            badbatch* bb;
            keydata key;

            line++;
            if (buf[0] == '#')
                continue;
            bb = malloc(sizeof(badbatch));
            if (sscanf(buf, "%s %s %d %d", otname, bb->fname, &bb->fid, &bb->batch) != 4)
                enkf_quit("%s, l.%d: wrong bad batch specification (expected \"%s %s %d %d\"\n", FNAME_BADBATCHES, line);

            bb->obstypeid = obstype_getid(obs->nobstypes, obs->obstypes, otname, 1);

            key.key_int[0] = bb->batch;
            key.key_short[2] = bb->obstypeid;
            key.key_short[3] = bb->fid;
            ht_insert(obs->badbatches, &key, bb);
            enkf_printf("    %s %s %d %d\n", otname, bb->fname, bb->fid, bb->batch);
        }
        fclose(f);
    }
#endif

    {
        int otid;

        for (otid = 0; otid <= obs->nobstypes; ++otid) {
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

#if defined(ENKF_PREP)
/**
 */
#define BYTE_PER_SHORT 2
#define BYTE_PER_INT 4
void obs_markbadbatches(observations* obs)
{
    int i;

    assert(BYTE_PER_SHORT == sizeof(short));
    assert(BYTE_PER_INT == sizeof(int));

    if (obs->badbatches == NULL || ht_getnentries(obs->badbatches) == 0)
        return;

    for (i = 0; i < obs->nobs; ++i) {
        observation* o = &obs->data[i];
        keydata key;
        badbatch* bb;

        key.key_int[0] = o->batch;
        key.key_short[2] = o->type;
        key.key_short[3] = o->fid;
        bb = ht_find(obs->badbatches, &key);

        if (bb != NULL && bb->obstypeid == o->type) {
            if (strcmp(bb->fname, st_findstringbyindex(obs->datafiles, o->fid)) != 0)
                enkf_quit("bad batch processing: file name for fid = %d in \"%s\" does not match the data file name. Check that \"%s\" is the right file.", o->fid, FNAME_BADBATCHES, FNAME_BADBATCHES);
            o->status = STATUS_BADBATCH;
        }
    }
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

#if defined(ENKF_PREP)
/**
 */
static int comp_obsstatus(const void* p1, const void* p2)
{
    observation* m1 = (observation*) p1;
    observation* m2 = (observation*) p2;

    if (m1->status > m2->status)
        return 1;
    if (m1->status < m2->status)
        return -1;
    return 0;
}

/** Move good observations to the head of the observation array.
 */
void obs_compact(observations* obs)
{
    int i;

    if (obs->compacted)
        return;

    enkf_flush();
    assert(STATUS_OK == 0);
    qsort(obs->data, obs->nobs, sizeof(observation), comp_obsstatus);
    for (i = 0; i < obs->nobs; ++i) {
        obs->data[i].id_orig = obs->data[i].id;
        obs->data[i].id = i;
    }
    enkf_printf("\n");
    enkf_flush();
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
    int dimid_nobs[1];
    size_t nobs;
    int varid_type, varid_product, varid_instrument, varid_id, varid_idorig, varid_fid, varid_batch, varid_value, varid_estd, varid_footprint, varid_lon, varid_lat, varid_depth, varid_mdepth, varid_fi, varid_fj, varid_fk, varid_time, varid_status, varid_aux;
    int* id = NULL;
    int* id_orig = NULL;
    short int* type = NULL;
    short int* product = NULL;
    short int* instrument = NULL;
    short int* fid = NULL;
    int* batch = NULL;
    double* value = NULL;
    double* estd = NULL;
    double* footprint = NULL;
    double* lon = NULL;
    double* lat = NULL;
    double* depth = NULL;
    double* model_depth = NULL;
    double* fi = NULL;
    double* fj = NULL;
    double* fk = NULL;
    double* time = NULL;
    int* status = NULL;
    int* aux = NULL;
    int natts;
    int i;

    ncw_open(fname, NC_NOWRITE, &ncid);

    if (ncw_att_exists(ncid, NC_GLOBAL, "DA_DAY"))
        ncw_get_att_double(ncid, NC_GLOBAL, "DA_DAY", &da_time);
    if (!enkf_noobsdatecheck && !isnan(da_time) && fabs(obs->da_time - da_time) > 1e-6)
        enkf_quit("observation data file \"%s\" from a different cycle");

    ncw_inq_dimid(ncid, "nobs", dimid_nobs);
    ncw_inq_dimlen(ncid, dimid_nobs[0], &nobs);

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

    ncw_inq_varid(ncid, "type", &varid_type);
    ncw_inq_varid(ncid, "product", &varid_product);
    ncw_inq_varid(ncid, "instrument", &varid_instrument);
    ncw_inq_varid(ncid, "id", &varid_id);
    ncw_inq_varid(ncid, "id_orig", &varid_idorig);
    ncw_inq_varid(ncid, "fid", &varid_fid);
    ncw_inq_varid(ncid, "batch", &varid_batch);
    ncw_inq_varid(ncid, "value", &varid_value);
    ncw_inq_varid(ncid, "estd", &varid_estd);
    if (ncw_var_exists(ncid, "footprint")) {
        ncw_inq_varid(ncid, "footprint", &varid_footprint);
        obs->has_nonpointobs = 1;
    }
    ncw_inq_varid(ncid, "lon", &varid_lon);
    ncw_inq_varid(ncid, "lat", &varid_lat);
    ncw_inq_varid(ncid, "depth", &varid_depth);
    ncw_inq_varid(ncid, "model_depth", &varid_mdepth);
    ncw_inq_varid(ncid, "fi", &varid_fi);
    ncw_inq_varid(ncid, "fj", &varid_fj);
    ncw_inq_varid(ncid, "fk", &varid_fk);
    ncw_inq_varid(ncid, "time", &varid_time);
    ncw_inq_varid(ncid, "status", &varid_status);
    ncw_inq_varid(ncid, "aux", &varid_aux);

    /*
     * date
     */
    if (obs->datestr == NULL) {
        size_t len = 0;

        ncw_inq_attlen(ncid, varid_time, "units", &len);
        obs->datestr = malloc(len + 1);
        ncw_get_att_text(ncid, varid_time, "units", obs->datestr);
        obs->datestr[len] = 0;
    }

    /*
     * type (basically, just a check)
     */
    ncw_inq_varnatts(ncid, varid_type, &natts);
    for (i = 0; i < natts; ++i) {
        char attname[NC_MAX_NAME];
        int typeid;

        ncw_inq_attname(ncid, varid_type, i, attname);
        typeid = obstype_getid(obs->nobstypes, obs->obstypes, attname, 0);
        if (typeid >= 0) {
            obstype* ot = &obs->obstypes[typeid];
            int typeid_read;

            ncw_check_attlen(ncid, varid_type, attname, 1);
            ncw_get_att_int(ncid, varid_type, attname, &typeid_read);
            assert(typeid == typeid_read);

            if (ot->logapplied) {
                char logattname[NC_MAX_NAME];
                char logattval[MAXSTRLEN];

                snprintf(logattname, NC_MAX_NAME, "%s:LOGAPPLIED", ot->name);
                if (!ncw_att_exists(ncid, varid_type, logattname))
                    enkf_quit("%s: variable = \"type\": expected attribute \"%s\" to be present for a log-transformed obs. type", fname, logattname);
                ncw_get_att_text(ncid, varid_type, logattname, logattval);
                if (strncmp(logattval, "true", 4) != 0)
                    enkf_quit("%s: variable = \"type\": expected attribute \"%s\" to have value \"true\" a log-transformed obs. type", fname, logattname);
            }
        }
    }

    /*
     * product 
     */
    ncw_inq_varnatts(ncid, varid_product, &natts);
    for (i = 0; i < natts; ++i) {
        char name[NC_MAX_NAME];
        nc_type nctype;
        size_t len;

        ncw_inq_attname(ncid, varid_product, i, name);
        ncw_inq_att(ncid, varid_product, name, &nctype, &len);
        if (nctype == NC_INT && len == 1) {
            int productid;

            ncw_get_att_int(ncid, varid_product, name, &productid);
            st_add(obs->products, name, productid);
        }
    }

    /*
     * instrument 
     */
    ncw_inq_varnatts(ncid, varid_instrument, &natts);
    for (i = 0; i < natts; ++i) {
        char name[NC_MAX_NAME];
        nc_type nctype;
        size_t len;

        ncw_inq_attname(ncid, varid_instrument, i, name);
        ncw_inq_att(ncid, varid_instrument, name, &nctype, &len);
        if (nctype == NC_INT && len == 1) {
            int instid;

            ncw_get_att_int(ncid, varid_instrument, name, &instid);
            st_add(obs->instruments, name, instid);
        }
    }

    /*
     * datafiles
     */
    ncw_inq_varnatts(ncid, varid_fid, &natts);
    for (i = 0; i < natts; ++i) {
        char name[NC_MAX_NAME];
        char attstr[MAXSTRLEN];
        size_t len;
        int fileid;

        ncw_inq_attname(ncid, varid_fid, i, name);
        if (!str2int(name, &fileid))
            continue;
        ncw_inq_attlen(ncid, varid_fid, name, &len);
        assert(len < MAXSTRLEN);
        ncw_get_att_text(ncid, varid_fid, name, attstr);
        attstr[len] = 0;
        st_add_ifabsent(obs->datafiles, attstr, fileid);
    }

#if defined(USE_SHMEM)
    if (sm_comm_rank == 0) {
#endif
        id = malloc(nobs * sizeof(int));
        id_orig = malloc(nobs * sizeof(int));
        type = malloc(nobs * sizeof(short int));
        product = malloc(nobs * sizeof(short int));
        instrument = malloc(nobs * sizeof(short int));
        fid = malloc(nobs * sizeof(short int));
        batch = malloc(nobs * sizeof(int));
        value = malloc(nobs * sizeof(double));
        estd = malloc(nobs * sizeof(double));
        if (obs->has_nonpointobs)
            footprint = malloc(nobs * sizeof(double));
        lon = malloc(nobs * sizeof(double));
        lat = malloc(nobs * sizeof(double));
        depth = malloc(nobs * sizeof(double));
        model_depth = malloc(nobs * sizeof(double));
        fi = malloc(nobs * sizeof(double));
        fj = malloc(nobs * sizeof(double));
        fk = malloc(nobs * sizeof(double));
        time = malloc(nobs * sizeof(double));
        status = malloc(nobs * sizeof(int));
        aux = malloc(nobs * sizeof(int));

        ncw_get_var_int(ncid, varid_id, id);
        ncw_get_var_int(ncid, varid_idorig, id_orig);
        ncw_get_var_short(ncid, varid_type, type);
        ncw_get_var_short(ncid, varid_product, product);
        ncw_get_var_short(ncid, varid_instrument, instrument);
        ncw_get_var_short(ncid, varid_fid, fid);
        ncw_get_var_int(ncid, varid_batch, batch);
        ncw_get_var_double(ncid, varid_value, value);
        ncw_get_var_double(ncid, varid_estd, estd);
        if (obs->has_nonpointobs)
            ncw_get_var_double(ncid, varid_footprint, footprint);
        ncw_get_var_double(ncid, varid_lon, lon);
        ncw_get_var_double(ncid, varid_lat, lat);
        ncw_get_var_double(ncid, varid_depth, depth);
        ncw_get_var_double(ncid, varid_mdepth, model_depth);
        ncw_get_var_double(ncid, varid_fi, fi);
        ncw_get_var_double(ncid, varid_fj, fj);
        ncw_get_var_double(ncid, varid_fk, fk);
        ncw_get_var_double(ncid, varid_time, time);
        ncw_get_var_int(ncid, varid_status, status);
        ncw_get_var_int(ncid, varid_aux, aux);
#if defined(USE_SHMEM)
    }
#endif

    ncw_close(ncid);

#if defined(USE_SHMEM)
    if (sm_comm_rank == 0) {
#endif
        for (i = 0; i < (int) nobs; ++i) {
            observation* o = &obs->data[i];

            o->type = type[i];
            o->product = product[i];
            o->instrument = instrument[i];
            o->id = id[i];
            o->id_orig = id_orig[i];
            o->fid = fid[i];
            o->batch = batch[i];
            o->value = value[i];
            o->estd = estd[i];
            if (obs->has_nonpointobs)
                o->footprint = footprint[i];
            o->lon = lon[i];
            o->lat = lat[i];
            o->depth = depth[i];
            o->model_depth = model_depth[i];
            o->fi = fi[i];
            o->fj = fj[i];
            o->fk = fk[i];
            /*
             * if because of the roundup error time gets outside the allowed
             * range, then correct it
             */
            {
                obstype* ot = &obs->obstypes[o->type];

                if (time[i] <= ot->obswindow_min)
                    o->time = ot->obswindow_min + EPSD;
                else if (time[i] >= ot->obswindow_max)
                    o->time = ot->obswindow_max - EPSD;
                else
                    o->time = time[i];
            }
            o->time = time[i];
            o->status = status[i];
            o->aux = aux[i];
        }

        free(type);
        free(product);
        free(instrument);
        free(id);
        free(id_orig);
        free(fid);
        free(batch);
        free(value);
        free(estd);
        if (obs->has_nonpointobs)
            free(footprint);
        free(lon);
        free(lat);
        free(depth);
        free(model_depth);
        free(fi);
        free(fj);
        free(fk);
        free(time);
        free(status);
        free(aux);
#if defined(USE_SHMEM)
    }
#endif
#if defined(MPI)
    MPI_Barrier(MPI_COMM_WORLD);
#endif

  finish:

    obs_calcstats(obs);
}

/**
 */
void obs_write(observations* obs, char fname[])
{
    int nobs = obs->nobs;
    char tunits[MAXSTRLEN];

    int ncid;
    int dimid_nobs[1];
    int varid_type, varid_product, varid_instrument, varid_id, varid_idorig, varid_fid, varid_batch, varid_value, varid_estd, varid_footprint, varid_lon, varid_lat, varid_depth, varid_mdepth, varid_fi, varid_fj, varid_fk, varid_time, varid_status, varid_aux;

    int* id;
    int* id_orig;
    short int* type;
    short int* product;
    short int* instrument;
    short int* fid;
    int* batch;
    double* value;
    double* estd;
    double* footprint = NULL;
    double* lon;
    double* lat;
    double* depth;
    double* model_depth;
    double* fi;
    double* fj;
    double* fk;
    double* time;
    int* status;
    int* aux;

    int i, ii;

    if (rank != 0)
        return;

    if (file_exists(fname))
        enkf_quit("file \"%s\" already exists", fname);
    ncw_create(fname, NC_NOCLOBBER | obs->ncformat, &ncid);

    ncw_put_att_double(ncid, NC_GLOBAL, "DA_DAY", 1, &obs->da_time);

    ncw_def_dim(ncid, "nobs", nobs, dimid_nobs);
    ncw_def_var(ncid, "id", NC_INT, 1, dimid_nobs, &varid_id);
    ncw_put_att_text(ncid, varid_id, "long_name", "observation ID");
    ncw_def_var(ncid, "id_orig", NC_INT, 1, dimid_nobs, &varid_idorig);
    ncw_put_att_text(ncid, varid_idorig, "long_name", "original observation ID");
    ncw_put_att_text(ncid, varid_idorig, "description", "for primary observations - the serial number of the primary observation during the reading of data files; for superobs - the original ID of the very first observation collated into this observation");
    ncw_def_var(ncid, "type", NC_SHORT, 1, dimid_nobs, &varid_type);
    ncw_put_att_text(ncid, varid_type, "long_name", "observation type ID");
    ncw_def_var(ncid, "product", NC_SHORT, 1, dimid_nobs, &varid_product);
    ncw_put_att_text(ncid, varid_product, "long_name", "observation product ID");
    ncw_def_var(ncid, "instrument", NC_SHORT, 1, dimid_nobs, &varid_instrument);
    ncw_put_att_text(ncid, varid_instrument, "long_name", "observation instrument ID");
    ncw_def_var(ncid, "fid", NC_SHORT, 1, dimid_nobs, &varid_fid);
    ncw_put_att_text(ncid, varid_fid, "long_name", "observation data file ID");
    ncw_def_var(ncid, "batch", NC_INT, 1, dimid_nobs, &varid_batch);
    ncw_put_att_text(ncid, varid_batch, "long_name", "observation batch ID");
    ncw_def_var(ncid, "value", NC_FLOAT, 1, dimid_nobs, &varid_value);
    ncw_put_att_text(ncid, varid_value, "long_name", "observation value");
    ncw_def_var(ncid, "estd", NC_FLOAT, 1, dimid_nobs, &varid_estd);
    ncw_put_att_text(ncid, varid_estd, "long_name", "standard deviation of observation error used in DA");
    if (obs->has_nonpointobs)
        ncw_def_var(ncid, "footprint", NC_FLOAT, 1, dimid_nobs, &varid_footprint);
    ncw_def_var(ncid, "lon", NC_FLOAT, 1, dimid_nobs, &varid_lon);
    ncw_put_att_text(ncid, varid_lon, "long_name", "observation longitude");
    ncw_def_var(ncid, "lat", NC_FLOAT, 1, dimid_nobs, &varid_lat);
    ncw_put_att_text(ncid, varid_lat, "long_name", "observation latitude");
    ncw_def_var(ncid, "depth", NC_FLOAT, 1, dimid_nobs, &varid_depth);
    ncw_put_att_text(ncid, varid_depth, "long_name", "observation depth/height");
    ncw_def_var(ncid, "model_depth", NC_FLOAT, 1, dimid_nobs, &varid_mdepth);
    ncw_put_att_text(ncid, varid_mdepth, "long_name", "model bottom depth at the observation location");
    ncw_def_var(ncid, "fi", NC_FLOAT, 1, dimid_nobs, &varid_fi);
    ncw_put_att_text(ncid, varid_fi, "long_name", "fractional grid index i of the observation");
    ncw_def_var(ncid, "fj", NC_FLOAT, 1, dimid_nobs, &varid_fj);
    ncw_put_att_text(ncid, varid_fj, "long_name", "fractional grid index j of the observation");
    ncw_def_var(ncid, "fk", NC_FLOAT, 1, dimid_nobs, &varid_fk);
    ncw_put_att_text(ncid, varid_fk, "long_name", "fractional grid index k of the observation");
    ncw_def_var(ncid, "time", NC_FLOAT, 1, dimid_nobs, &varid_time);
    ncw_put_att_text(ncid, varid_time, "long_name", "observation time");
    ncw_def_var(ncid, "status", NC_BYTE, 1, dimid_nobs, &varid_status);
    ncw_put_att_text(ncid, varid_status, "long_name", "observation status");
    i = STATUS_OK;
    ncw_put_att_int(ncid, varid_status, "STATUS_OK", 1, &i);
    i = STATUS_OUTSIDEGRID;
    ncw_put_att_int(ncid, varid_status, "STATUS_OUTSIDEGRID", 1, &i);
    i = STATUS_LAND;
    ncw_put_att_int(ncid, varid_status, "STATUS_LAND", 1, &i);
    i = STATUS_SHALLOW;
    ncw_put_att_int(ncid, varid_status, "STATUS_SHALLOW", 1, &i);
    i = STATUS_RANGE;
    ncw_put_att_int(ncid, varid_status, "STATUS_RANGE", 1, &i);
    i = STATUS_BADBATCH;
    ncw_put_att_int(ncid, varid_status, "STATUS_BADBATCH", 1, &i);
    i = STATUS_OUTSIDEOBSDOMAIN;
    ncw_put_att_int(ncid, varid_status, "STATUS_OUTSIDEOBSDOMAIN", 1, &i);
    i = STATUS_OUTSIDEOBSWINDOW;
    ncw_put_att_int(ncid, varid_status, "STATUS_OUTSIDEOBSWINDOW", 1, &i);
    i = STATUS_THINNED;
    ncw_put_att_int(ncid, varid_status, "STATUS_THINNED", 1, &i);
    i = STATUS_EXCLUDED;
    ncw_put_att_int(ncid, varid_status, "STATUS_EXCLUDED", 1, &i);
    ncw_def_var(ncid, "aux", NC_INT, 1, dimid_nobs, &varid_aux);
    ncw_put_att_text(ncid, varid_aux, "long_name", "auxiliary information");
    ncw_put_att_text(ncid, varid_aux, "description", "for primary observations - the ID of the superobservation it is collated into; for superobservations - the number of primary observations collated");
    snprintf(tunits, MAXSTRLEN, "days from %s", obs->datestr);
    ncw_put_att_text(ncid, varid_time, "units", tunits);

    for (i = 0; i < obs->nobstypes; ++i) {
        ncw_put_att_int(ncid, varid_type, obs->obstypes[i].name, 1, &i);
        if (obs->obstypes[i].logapplied) {
            char attname[NC_MAX_NAME];

            snprintf(attname, NC_MAX_NAME, "%s:LOGAPPLIED", obs->obstypes[i].name);
            ncw_put_att_text(ncid, varid_type, attname, "true");
        }
    }

    for (i = 0; i < st_getsize(obs->products); ++i)
        ncw_put_att_int(ncid, varid_product, st_findstringbyindex(obs->products, i), 1, &i);

    for (i = 0; i < st_getsize(obs->instruments); ++i)
        ncw_put_att_int(ncid, varid_instrument, st_findstringbyindex(obs->instruments, i), 1, &i);

    for (i = 0; i < st_getsize(obs->datafiles); ++i) {
        char attname[NC_MAX_NAME];
        char* datafile = st_findstringbyindex(obs->datafiles, i);

        snprintf(attname, NC_MAX_NAME, "%d", i);
        ncw_put_att_text(ncid, varid_fid, attname, datafile);
    }

    if (obs->nccompression > 0)
        ncw_def_deflate(ncid, 0, 1, obs->nccompression);
    ncw_enddef(ncid);

    if (nobs == 0) {
        ncw_close(ncid);
        return;
    }

    id = malloc(nobs * sizeof(int));
    id_orig = malloc(nobs * sizeof(int));
    type = malloc(nobs * sizeof(short int));
    product = malloc(nobs * sizeof(short int));
    instrument = malloc(nobs * sizeof(short int));
    fid = malloc(nobs * sizeof(short int));
    batch = malloc(nobs * sizeof(int));
    value = malloc(nobs * sizeof(double));
    estd = malloc(nobs * sizeof(double));
    if (obs->has_nonpointobs)
        footprint = malloc(nobs * sizeof(footprint));
    lon = malloc(nobs * sizeof(double));
    lat = malloc(nobs * sizeof(double));
    depth = malloc(nobs * sizeof(double));
    model_depth = malloc(nobs * sizeof(double));
    fi = malloc(nobs * sizeof(double));
    fj = malloc(nobs * sizeof(double));
    fk = malloc(nobs * sizeof(double));
    time = malloc(nobs * sizeof(double));
    status = malloc(nobs * sizeof(int));
    aux = malloc(nobs * sizeof(int));

    for (i = 0, ii = 0; i < obs->nobs; ++i) {
        observation* o = &obs->data[i];

        if (!isfinite(o->value) || fabs(o->value) > FLT_MAX || !isfinite((float) o->value))
            enkf_quit("bad value");

        type[ii] = o->type;
        product[ii] = o->product;
        instrument[ii] = o->instrument;
        id[ii] = o->id;
        fid[ii] = o->fid;
        batch[ii] = o->batch;
        /*
         * id of the first ob contributed to this sob 
         */
        id_orig[ii] = o->id_orig;
        value[ii] = o->value;
        estd[ii] = o->estd;
        if (obs->has_nonpointobs)
            footprint[ii] = o->footprint;
        lon[ii] = o->lon;
        lat[ii] = o->lat;
        depth[ii] = o->depth;
        model_depth[ii] = o->model_depth;
        fi[ii] = o->fi;
        fj[ii] = o->fj;
        fk[ii] = o->fk;
        time[ii] = o->time;
        status[ii] = o->status;
        aux[ii] = o->aux;
        ii++;
    }
    assert(ii == nobs);

    ncw_put_var_int(ncid, varid_id, id);
    ncw_put_var_int(ncid, varid_idorig, id_orig);
    ncw_put_var_short(ncid, varid_type, type);
    ncw_put_var_short(ncid, varid_product, product);
    ncw_put_var_short(ncid, varid_instrument, instrument);
    ncw_put_var_short(ncid, varid_fid, fid);
    ncw_put_var_int(ncid, varid_batch, batch);
    ncw_put_var_double(ncid, varid_value, value);
    ncw_put_var_double(ncid, varid_estd, estd);
    if (obs->has_nonpointobs)
        ncw_put_var_double(ncid, varid_footprint, footprint);
    ncw_put_var_double(ncid, varid_lon, lon);
    ncw_put_var_double(ncid, varid_lat, lat);
    ncw_put_var_double(ncid, varid_depth, depth);
    ncw_put_var_double(ncid, varid_mdepth, model_depth);
    ncw_put_var_double(ncid, varid_fi, fi);
    ncw_put_var_double(ncid, varid_fj, fj);
    ncw_put_var_double(ncid, varid_fk, fk);
    ncw_put_var_double(ncid, varid_time, time);
    ncw_put_var_int(ncid, varid_status, status);
    ncw_put_var_int(ncid, varid_aux, aux);

    ncw_close(ncid);
    free(type);
    free(product);
    free(instrument);
    free(id);
    free(id_orig);
    free(fid);
    free(batch);
    free(value);
    free(estd);
    if (obs->has_nonpointobs)
        free(footprint);
    free(lon);
    free(lat);
    free(depth);
    free(model_depth);
    free(fi);
    free(fj);
    free(fk);
    free(time);
    free(status);
    free(aux);
#if defined(MPI)
    MPI_Barrier(MPI_COMM_WORLD);
#endif
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
void obs_superob(observations* obs, __compar_d_fn_t cmp_obs, observations** sobs, int sobid, int do_thin)
{
    int i1 = 0, i2 = 0;
    int nsobs = 0;
    int nthinned = 0;
    observation* data = obs->data;
    observation* sdata = NULL;
    int has_nonpointobs = 0;

    qsort_r(obs->data, obs->ngood, sizeof(observation), cmp_obs, obs);

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
         * thin observations with identical positions (these are supposedly
         * coming from high-frequency instruments)
         */
        if (do_thin && obs->obstypes[data[i1].type].can_thin) {
            int i11 = i1;
            int i22 = i1;

            qsort(&data[i1], i2 - i1 + 1, sizeof(observation), cmp_xyz);
            while (i22 <= i2) {
                while (i22 + 1 <= i2 && cmp_xyz(&data[i11], &data[i22 + 1]) == 0)
                    i22++;
                for (ii = i11 + 1; ii <= i22; ++ii) {
                    data[ii].status = STATUS_THINNED;
                    nthinned++;
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

        so = &sdata[nsobs];
        o = &data[i1];

        assert(o->status == STATUS_OK);

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
        so->fi = o->fi;
        so->fj = o->fj;
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
            so->fi = so->fi * so->estd + o->fi / evar;
            so->fj = so->fj * so->estd + o->fj / evar;
            so->fk = so->fk * so->estd + o->fk / evar;
            so->time = so->time * so->estd + o->time / evar;

            so->estd += 1.0 / evar;

            so->value /= so->estd;
            so->lon /= so->estd;
            so->lat /= so->estd;
            so->depth /= so->estd;
            so->model_depth /= so->estd;
            so->fi /= so->estd;
            so->fj /= so->estd;
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

/**
 */
void obs_find_bytype(observations* obs, int type, int* nobs, int** obsids)
{
    int i;

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

/**
 */
void obs_find_bytypeandtime(observations* obs, int type, int time, int* nobs, int** obsids)
{
    obstype* ot = &obs->obstypes[type];
    int i;

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

/**
 */
void obs_printob(observations* obs, int i)
{
    observation* o = &obs->data[i];

    enkf_printf("type = %s, product = %s, instrument = %s, datafile = %s, id = %d, original id = %d, batch = %d, value = %.3g, estd = %.3g, footprint = %.3g, ", obs->obstypes[o->type].name, st_findstringbyindex(obs->products, o->product), st_findstringbyindex(obs->instruments, o->instrument), st_findstringbyindex(obs->datafiles, o->fid), o->id, o->id_orig, (int) o->batch, o->value, o->estd, o->footprint);
    enkf_printf("lon = %.3f, lat = %.3f, depth = %.1f, model_depth = %.1f, fi = %.3f, fj = %.3f, fk = %.3f, day = %.3g, status = %d\n", o->lon, o->lat, o->depth, o->model_depth, o->fi, o->fj, o->fk, o->time, o->status);
}

#if defined(ENKF_CALC)
/**
 */
static double distance(double xyz1[3], double xyz2[3])
{
    return sqrt((xyz1[0] - xyz2[0]) * (xyz1[0] - xyz2[0]) + (xyz1[1] - xyz2[1]) * (xyz1[1] - xyz2[1]) + (xyz1[2] - xyz2[2]) * (xyz1[2] - xyz2[2]));
}
#endif

#if defined(ENKF_CALC)
/**
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
        int* ids = NULL;
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
        ids = malloc(nobs * sizeof(int));
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
void obs_findlocal(observations* obs, double lon, double lat, char* domainname, int* n, int** ids, double** lcoeffs, int* ploc_allocated)
{
    double ll[2] = { lon, lat };
    double xyz[3];
    int otid;
    int i, ntot, ngood;

    ll2xyz(ll, xyz);

    if (obs->nobstypes == 0)
        return;
    if (obs->loctrees == NULL)
        obs_createkdtrees(obs);

    for (otid = 0, i = 0; otid < obs->nobstypes; ++otid) {
        obstype* ot = &obs->obstypes[otid];
        kdtree* tree = obs->loctrees[otid];
        int* obsids = obs->obsids[otid];
        kdset* set = NULL;
        size_t id;
        int iloc;

        if (ot->nobs == 0 || ot->statsonly)
            continue;

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

        set = kd_findnodeswithinrange(tree, xyz, obstype_getmaxlocrad(ot), 1);
        for (iloc = 0; iloc < ot->nlobsmax && (id = kdset_readnext(set, NULL)) != SIZE_MAX; ++i, ++iloc) {
            size_t id_orig = kd_getnodedata(tree, id);

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
        }
        kdset_free(set);
    }

    /*
     * compact the result, calculate taper coefficients
     */
    ntot = i;
    for (i = 0, ngood = 0; i < ntot; ++i) {
        int id = (*ids)[i];
        observation* o = &obs->data[id];
        double ll2[2] = { o->lon, o->lat };
        double xyz2[3];

        if (o->status != STATUS_OK)
            continue;

        ll2xyz(ll2, xyz2);
        (*ids)[ngood] = id;
        (*lcoeffs)[ngood] = obstype_calclcoeff(&obs->obstypes[o->type], distance(xyz, xyz2));
        ngood++;
    }
    *n = ngood;
#if defined (MINIMISE_ALLOC)
    if (ploc_allocated == NULL && ngood == 0 && *ids != NULL) {
#else
    if (ngood == 0 && *ids != NULL) {
#endif
        free(*ids);
        free(*lcoeffs);
    }
}
#endif
