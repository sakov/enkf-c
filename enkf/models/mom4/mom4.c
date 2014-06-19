/******************************************************************************
 *
 * File:        mom4.c        
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "definitions.h"
#include "utils.h"
#include "grid.h"
#include "enkfprm.h"
#include "model.h"
#include "mom4.h"
#include "ncw.h"

#define EPS 1.0e-6
#define EPSLON 1.0e-3

typedef struct {
    char* gridspec;
    char* xdimname;
    char* ydimname;
    char* zdimname;
    char* xvarname;
    char* yvarname;
    char* zvarname;
    char* depthvarname;
    char* numlevelsvarname;
} gridprm;

/**
 */
static gridprm* gridprm_create(char* fname)
{
    gridprm* prm = NULL;

    FILE* f = NULL;
    char buf[MAXSTRLEN];
    int line;

    prm = calloc(1, sizeof(gridprm));

    f = enkf_fopen(fname, "r");

    line = 0;
    while (fgets(buf, MAXSTRLEN, f) != NULL) {
        char seps[] = " =\t\n";
        char* token;

        line++;
        if (buf[0] == '#')
            continue;
        if ((token = strtok(buf, seps)) == NULL)
            continue;
        if (strcasecmp(token, "GRIDSPEC") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: GRIDSPEC not specified", fname, line);
            else if (prm->gridspec != NULL)
                enkf_quit("%s, l.%d: GRIDSPEC specified twice", fname, line);
            else
                prm->gridspec = strdup(token);
        } else if (strcasecmp(token, "XDIMNAME") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: XDIMNAME not specified", fname, line);
            else if (prm->xdimname != NULL)
                enkf_quit("%s, l.%d: XDIMNAME specified twice", fname, line);
            else
                prm->xdimname = strdup(token);
        } else if (strcasecmp(token, "YDIMNAME") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: YDIMNAME not specified", fname, line);
            else if (prm->ydimname != NULL)
                enkf_quit("%s, l.%d: YDIMNAME specified twice", fname, line);
            else
                prm->ydimname = strdup(token);
        } else if (strcasecmp(token, "ZDIMNAME") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: ZDIMNAME not specified", fname, line);
            else if (prm->zdimname != NULL)
                enkf_quit("%s, l.%d: ZDIMNAME specified twice", fname, line);
            else
                prm->zdimname = strdup(token);
        } else if (strcasecmp(token, "XVARNAME") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: XVARNAME not specified", fname, line);
            else if (prm->xvarname != NULL)
                enkf_quit("%s, l.%d: XVARNAME specified twice", fname, line);
            else
                prm->xvarname = strdup(token);
        } else if (strcasecmp(token, "YVARNAME") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: YVARNAME not specified", fname, line);
            else if (prm->yvarname != NULL)
                enkf_quit("%s, l.%d: YVARNAME specified twice", fname, line);
            else
                prm->yvarname = strdup(token);
        } else if (strcasecmp(token, "ZVARNAME") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: ZVARNAME not specified", fname, line);
            else if (prm->zvarname != NULL)
                enkf_quit("%s, l.%d: ZVARNAME specified twice", fname, line);
            else
                prm->zvarname = strdup(token);
        } else if (strcasecmp(token, "DEPTHVARNAME") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: DEPTHVARNAME not specified", fname, line);
            else if (prm->depthvarname != NULL)
                enkf_quit("%s, l.%d: DEPTHVARNAME specified twice", fname, line);
            else
                prm->depthvarname = strdup(token);
        } else if (strcasecmp(token, "NUMLEVELSVARNAME") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: NUMLEVELSVARNAME not specified", fname, line);
            else if (prm->numlevelsvarname != NULL)
                enkf_quit("%s, l.%d: NUMLEVELSVARNAME specified twice", fname, line);
            else
                prm->numlevelsvarname = strdup(token);
        }
    }

    fclose(f);

    if (prm->gridspec == NULL)
        enkf_quit("%s: GRIDSPEC not specified", fname);
    if (prm->xdimname == NULL)
        enkf_quit("%s: XDIMNAME not specified", fname);
    if (prm->ydimname == NULL)
        enkf_quit("%s: YDIMNAME not specified", fname);
    if (prm->zdimname == NULL)
        enkf_quit("%s: ZDIMNAME not specified", fname);
    if (prm->xvarname == NULL)
        enkf_quit("%s: XVARNAME not specified", fname);
    if (prm->yvarname == NULL)
        enkf_quit("%s: YVARNAME not specified", fname);
    if (prm->zvarname == NULL)
        enkf_quit("%s: ZVARNAME not specified", fname);
    if (prm->depthvarname == NULL)
        enkf_quit("%s: DEPTHVARNAME not specified", fname);
    if (prm->numlevelsvarname == NULL)
        enkf_quit("%s: NUMLEVELSVARNAME not specified", fname);

    return prm;
}

/**
 */
static void gridprm_destroy(gridprm* prm)
{
    free(prm->gridspec);
    free(prm->xdimname);
    free(prm->ydimname);
    free(prm->zdimname);
    free(prm->xvarname);
    free(prm->yvarname);
    free(prm->zvarname);
    free(prm->depthvarname);
    free(prm->numlevelsvarname);
    free(prm);
}

/**
 */
static void mom4_getmemberfname(model* m, char ensdir[], char varname[], int mem, char fname[])
{
    snprintf(fname, MAXSTRLEN, "%s/mem%03d_%s.nc", ensdir, mem, varname);
}

/**
 */
static int mom4_getmemberfname_async(model* m, char ensdir[], char varname[], char otname[], int mem, int t, char fname[])
{
    snprintf(fname, MAXSTRLEN, "%s/mem%03d_%s_%d.nc", ensdir, mem, varname, t);
    if (!file_exists(fname)) {
        snprintf(fname, MAXSTRLEN, "%s/mem%03d_%s.nc", ensdir, mem, varname);
        return 0;
    }
    return 1;
}

/**
 */
static void mom4_getbgfname(model* m, char ensdir[], char varname[], char fname[])
{
    snprintf(fname, MAXSTRLEN, "%s/bg_%s.nc", ensdir, varname);
}

/**
 */
static int mom4_getbgfname_async(model* m, char bgdir[], char varname[], char otname[], int t, char fname[])
{
    snprintf(fname, MAXSTRLEN, "%s/bg_%s_%d.nc", bgdir, varname, t);
    if (!file_exists(fname)) {
        snprintf(fname, MAXSTRLEN, "%s/bg_%s.nc", bgdir, varname);
        return 0;
    }
    return 1;
}

/**
 */
static void mom4_readfield(model* m, char fname[], int mem, int time, char varname[], int k, float* v)
{
    readfield(fname, k, varname, v);
}

/**
 */
static void mom4_read3dfield(model* m, char fname[], int mem, int time, char varname[], float* v)
{
    read3dfield(fname, varname, v);
}

/**
 */
static void mom4_writefield(model* m, char fname[], int time, char varname[], int k, float* v)
{
    writefield(fname, k, varname, v);
}

/**
 */
void mom4_setup(model* m, char gridspec[])
{
    char* fname;
    gridprm* prm;
    int ncid;
    int dimid_x, dimid_y, dimid_z;
    int varid_x, varid_y, varid_z;
    int ndims_x, ndims_y, ndims_z;
    size_t nx, ny, nz;
    int varid_depth, varid_numlevels;

    prm = gridprm_create(gridspec);
    fname = prm->gridspec;
    enkf_printf("    grid file = \"%s\"\n", fname);

    ncw_open(fname, NC_NOWRITE, &ncid);
    ncw_inq_dimid(fname, ncid, prm->xdimname, &dimid_x);
    ncw_inq_dimid(fname, ncid, prm->ydimname, &dimid_y);
    ncw_inq_dimid(fname, ncid, prm->zdimname, &dimid_z);
    ncw_inq_dimlen(fname, ncid, dimid_x, &nx);
    ncw_inq_dimlen(fname, ncid, dimid_y, &ny);
    ncw_inq_dimlen(fname, ncid, dimid_z, &nz);

    enkf_printf("    grid dimensions = %u x %u x %u\n", (unsigned int) nx, (unsigned int) ny, (unsigned int) nz);

    ncw_inq_varid(fname, ncid, prm->xvarname, &varid_x);
    ncw_inq_varid(fname, ncid, prm->yvarname, &varid_y);
    ncw_inq_varid(fname, ncid, prm->zvarname, &varid_z);

    ncw_inq_varndims(fname, ncid, varid_x, &ndims_x);
    ncw_inq_varndims(fname, ncid, varid_y, &ndims_y);
    ncw_inq_varndims(fname, ncid, varid_z, &ndims_z);

    if (ndims_x == 1 && ndims_y == 1) {
        double* x;
        double* y;
        double* z;
        int i;
        double dx, dy;
        int periodic_x;

        x = malloc(nx * sizeof(double));
        y = malloc(ny * sizeof(double));
        z = malloc(nz * sizeof(double));

        ncw_get_var_double(fname, ncid, varid_x, x);
        ncw_get_var_double(fname, ncid, varid_y, y);
        ncw_get_var_double(fname, ncid, varid_z, z);

        periodic_x = fabs(fmod(2.0 * x[nx - 1] - x[nx - 2], 360.0) - x[0]) < EPSLON;

        dx = (x[nx - 1] - x[0]) / (double) (nx - 1);
        for (i = 1; i < nx; ++i)
            if (fabs(x[i] - x[i - 1] - dx) / fabs(dx) > EPS)
                break;
        if (i != nx)
            grid_set(model_getgrid(m), GRIDTYPE_LATLON_IRREGULAR, periodic_x, 0, nx, ny, nz, x, y, z);
        else {
            dy = (y[ny - 1] - y[0]) / (double) (ny - 1);
            for (i = 1; i < ny; ++i)
                if (fabs(y[i] - y[i - 1] - dy) / fabs(dy) > EPS)
                    break;
            if (i != ny)
                grid_set(model_getgrid(m), GRIDTYPE_LATLON_IRREGULAR, periodic_x, 0, nx, ny, nz, x, y, z);
            else
                grid_set(model_getgrid(m), GRIDTYPE_LATLON_REGULAR, periodic_x, 0, nx, ny, nz, x, y, z);
        }
    } else if (ndims_x == 2 && ndims_y == 2) {
        double** x;
        double** y;
        double* z;

        x = alloc2d(ny, nx, sizeof(double));
        y = alloc2d(ny, nx, sizeof(double));
        z = malloc(nz * sizeof(double));

        ncw_get_var_double(fname, ncid, varid_x, x[0]);
        ncw_get_var_double(fname, ncid, varid_y, y[0]);
        ncw_get_var_double(fname, ncid, varid_z, z);

        grid_set(model_getgrid(m), GRIDTYPE_CURVILINEAR, 0, 0, nx, ny, nz, x, y, z);
    } else
        enkf_quit("%s: could not determine the grid type", fname);

    ncw_inq_varid(fname, ncid, prm->depthvarname, &varid_depth);
    {
        float** depth = alloc2d(ny, nx, sizeof(float));

        ncw_get_var_float(fname, ncid, varid_depth, depth[0]);
        model_setdepth(m, depth);
    }

    ncw_inq_varid(fname, ncid, prm->numlevelsvarname, &varid_numlevels);
    {
        int** numlevels = alloc2d(ny, nx, sizeof(int));

        ncw_get_var_int(fname, ncid, varid_numlevels, numlevels[0]);
        model_setnumlevels(m, numlevels);
    }

    ncw_close(fname, ncid);
    gridprm_destroy(prm);

    model_setgetmemberfname_fn(m, mom4_getmemberfname);
    model_setgetmemberfnameasync_fn(m, mom4_getmemberfname_async);
    model_setbgfname_fn(m, mom4_getbgfname);
    model_setbgfnameasync_fn(m, mom4_getbgfname_async);
    model_setreadfield_fn(m, mom4_readfield);
    model_setread3dfield_fn(m, mom4_read3dfield);
    model_setwritefield_fn(m, mom4_writefield);
}
