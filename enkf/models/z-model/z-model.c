/******************************************************************************
 *
 * File:        z-model.c        
 *
 * Created:     12/2012
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Description:
 *
 * Revisions:   PS 1.7.2014 -- renamed from mom4.c
 *
 *****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "definitions.h"
#include "utils.h"
#include "grid.h"
#include "enkfprm.h"
#include "model.h"
#include "z-model.h"
#include "ncw.h"

#define EPSLON 1.0e-3

typedef struct {
    char* name;
    char* fname;
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
static void gridprm_create(char* fname, int* ngrid, gridprm** prm)
{
    gridprm* now = NULL;
    FILE* f = NULL;
    char buf[MAXSTRLEN];
    int line;
    int i;

    assert(*ngrid == 0 && *prm == NULL);

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
        if (strcasecmp(token, "NAME") == 0) {
            (*ngrid)++;
            *prm = realloc(*prm, sizeof(gridprm) * (*ngrid));
            now = &(*prm)[*ngrid - 1];
            memset(now, 0, sizeof(gridprm));
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: NAME not specified", fname, line);
            else
                now->name = strdup(token);
        } else if (strcasecmp(token, "DATA") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: DATA not specified", fname, line);
            else if (now->fname != NULL)
                enkf_quit("%s, l.%d: DATA specified twice", fname, line);
            else
                now->fname = strdup(token);
        } else if (strcasecmp(token, "XDIMNAME") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: XDIMNAME not specified", fname, line);
            else if (now->xdimname != NULL)
                enkf_quit("%s, l.%d: XDIMNAME specified twice", fname, line);
            else
                now->xdimname = strdup(token);
        } else if (strcasecmp(token, "YDIMNAME") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: YDIMNAME not specified", fname, line);
            else if (now->ydimname != NULL)
                enkf_quit("%s, l.%d: YDIMNAME specified twice", fname, line);
            else
                now->ydimname = strdup(token);
        } else if (strcasecmp(token, "ZDIMNAME") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: ZDIMNAME not specified", fname, line);
            else if (now->zdimname != NULL)
                enkf_quit("%s, l.%d: ZDIMNAME specified twice", fname, line);
            else
                now->zdimname = strdup(token);
        } else if (strcasecmp(token, "XVARNAME") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: XVARNAME not specified", fname, line);
            else if (now->xvarname != NULL)
                enkf_quit("%s, l.%d: XVARNAME specified twice", fname, line);
            else
                now->xvarname = strdup(token);
        } else if (strcasecmp(token, "YVARNAME") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: YVARNAME not specified", fname, line);
            else if (now->yvarname != NULL)
                enkf_quit("%s, l.%d: YVARNAME specified twice", fname, line);
            else
                now->yvarname = strdup(token);
        } else if (strcasecmp(token, "ZVARNAME") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: ZVARNAME not specified", fname, line);
            else if (now->zvarname != NULL)
                enkf_quit("%s, l.%d: ZVARNAME specified twice", fname, line);
            else
                now->zvarname = strdup(token);
        } else if (strcasecmp(token, "DEPTHVARNAME") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: DEPTHVARNAME not specified", fname, line);
            else if (now->depthvarname != NULL)
                enkf_quit("%s, l.%d: DEPTHVARNAME specified twice", fname, line);
            else
                now->depthvarname = strdup(token);
        } else if (strcasecmp(token, "NUMLEVELSVARNAME") == 0) {
            if ((token = strtok(NULL, seps)) == NULL)
                enkf_quit("%s, l.%d: NUMLEVELSVARNAME not specified", fname, line);
            else if (now->numlevelsvarname != NULL)
                enkf_quit("%s, l.%d: NUMLEVELSVARNAME specified twice", fname, line);
            else
                now->numlevelsvarname = strdup(token);
        }
    }

    fclose(f);

    for (i = 0; i < *ngrid; ++i) {
        gridprm* now  = &(*prm)[i];

        if (now->name == NULL)
            enkf_quit("%s: NAME not specified", fname);
        if (now->fname == NULL)
            enkf_quit("%s: DATA not specified", fname);
        if (now->xdimname == NULL)
            enkf_quit("%s: XDIMNAME not specified", fname);
        if (now->ydimname == NULL)
            enkf_quit("%s: YDIMNAME not specified", fname);
        if (now->zdimname == NULL)
            enkf_quit("%s: ZDIMNAME not specified", fname);
        if (now->xvarname == NULL)
            enkf_quit("%s: XVARNAME not specified", fname);
        if (now->yvarname == NULL)
            enkf_quit("%s: YVARNAME not specified", fname);
        if (now->zvarname == NULL)
            enkf_quit("%s: ZVARNAME not specified", fname);
        if (now->depthvarname == NULL)
            enkf_quit("%s: DEPTHVARNAME not specified", fname);
        if (now->numlevelsvarname == NULL)
            enkf_quit("%s: NUMLEVELSVARNAME not specified", fname);
    }
}

/**
 */
static void gridprm_destroy(int ngrid, gridprm prm[])
{
    int i;

    for (i = 0; i < ngrid; ++i) {
        gridprm* now = &prm[i];

        free(now->name);
        free(now->fname);
        free(now->xdimname);
        free(now->ydimname);
        free(now->zdimname);
        free(now->xvarname);
        free(now->yvarname);
        free(now->zvarname);
        free(now->depthvarname);
        free(now->numlevelsvarname);
    }
    free(prm);
}

/**
 */
static void gridprm_describe(gridprm* prm, char offset[])
{
    enkf_printf("%sgrid prm info:\n", offset);
    enkf_printf("%s  name = \"%s\"\n", offset, prm->name);
    enkf_printf("%s  file = \"%s\"\n", offset, prm->fname);
    enkf_printf("%s  xdimname = \"%s\"\n", offset, prm->xdimname);
    enkf_printf("%s  ydimname = \"%s\"\n", offset, prm->ydimname);
    enkf_printf("%s  zdimname = \"%s\"\n", offset, prm->zdimname);
    enkf_printf("%s  xvarname = \"%s\"\n", offset, prm->xvarname);
    enkf_printf("%s  yvarname = \"%s\"\n", offset, prm->yvarname);
    enkf_printf("%s  zvarname = \"%s\"\n", offset, prm->zvarname);
    enkf_printf("%s  depthvarname = \"%s\"\n", offset, prm->depthvarname);
    enkf_printf("%s  numlevelsvarname = \"%s\"\n", offset, prm->numlevelsvarname);
}


/**
 */
static void ll2xyz(double in[2], double out[3])
{
    double lon = in[0] * DEG2RAD;
    double lat = in[1] * DEG2RAD;
    double coslat = cos(lat);

    out[0] = REARTH * sin(lon) * coslat;
    out[1] = REARTH * cos(lon) * coslat;
    out[2] = REARTH * sin(lat);
}

/**
 */
void zmodel_setgrids(model* m, char gfname[])
{
    int ngrid = 0;
    gridprm* prm = NULL;
    int ncid;
    int dimid_x, dimid_y, dimid_z;
    int varid_x, varid_y, varid_z;
    int ndims_x, ndims_y, ndims_z;
    size_t nx, ny, nz;
    int varid_depth, varid_numlevels;
    int i;

    gridprm_create(gfname, &ngrid, &prm);
    assert(ngrid > 0);

    for (i = 0; i < ngrid; ++i) {
        gridprm* now = &prm[i];
        grid* g = NULL;
        char* fname = now->fname;

        g = grid_create(now->name, i);

        ncw_open(fname, NC_NOWRITE, &ncid);
        ncw_inq_dimid(fname, ncid, now->xdimname, &dimid_x);
        ncw_inq_dimid(fname, ncid, now->ydimname, &dimid_y);
        ncw_inq_dimid(fname, ncid, now->zdimname, &dimid_z);
        ncw_inq_dimlen(fname, ncid, dimid_x, &nx);
        ncw_inq_dimlen(fname, ncid, dimid_y, &ny);
        ncw_inq_dimlen(fname, ncid, dimid_z, &nz);

        ncw_inq_varid(fname, ncid, now->xvarname, &varid_x);
        ncw_inq_varid(fname, ncid, now->yvarname, &varid_y);
        ncw_inq_varid(fname, ncid, now->zvarname, &varid_z);

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
            for (i = 1; i < (int) nx; ++i)
                if (fabs(x[i] - x[i - 1] - dx) / fabs(dx) > EPSLON)
                    break;
            if (i != (int) nx)
                grid_setcoords(g, GRIDHTYPE_LATLON_IRREGULAR, periodic_x, 0, nx, ny, nz, x, y, z);
            else {
                dy = (y[ny - 1] - y[0]) / (double) (ny - 1);
                for (i = 1; i < (int) ny; ++i)
                    if (fabs(y[i] - y[i - 1] - dy) / fabs(dy) > EPSLON)
                        break;
                if (i != (int) ny)
                    grid_setcoords(g, GRIDHTYPE_LATLON_IRREGULAR, periodic_x, 0, nx, ny, nz, x, y, z);
                else
                grid_setcoords(g, GRIDHTYPE_LATLON_REGULAR, periodic_x, 0, nx, ny, nz, x, y, z);
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

            grid_setcoords(g, GRIDHTYPE_CURVILINEAR, 0, 0, nx, ny, nz, x, y, z);
        } else
            enkf_quit("%s: could not determine the grid type", fname);

        model_setgrid(m, g);
        grid_settocartesian_fn(g, ll2xyz);

        ncw_inq_varid(fname, ncid, now->depthvarname, &varid_depth);
        {
            float** depth = alloc2d(ny, nx, sizeof(float));

            ncw_get_var_float(fname, ncid, varid_depth, depth[0]);
            grid_setdepth(g, depth);
        }

        ncw_inq_varid(fname, ncid, now->numlevelsvarname, &varid_numlevels);
        {
            int** numlevels = alloc2d(ny, nx, sizeof(int));

            ncw_get_var_int(fname, ncid, varid_numlevels, numlevels[0]);
            grid_setnumlevels(g, numlevels);
        }
        ncw_close(fname, ncid);

        gridprm_describe(now, "    ");
        grid_print(g, "    ");
    }

    gridprm_destroy(ngrid, prm);
}

/**
 */
static void zmodel_getmemberfname(model* m, char ensdir[], char varname[], int mem, char fname[])
{
    snprintf(fname, MAXSTRLEN, "%s/mem%03d_%s.nc", ensdir, mem, varname);
}

/**
 */
static int zmodel_getmemberfname_async(model* m, char ensdir[], char varname[], char otname[], int mem, int t, char fname[])
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
static void zmodel_getbgfname(model* m, char ensdir[], char varname[], char fname[])
{
    snprintf(fname, MAXSTRLEN, "%s/bg_%s.nc", ensdir, varname);
}

/**
 */
static int zmodel_getbgfname_async(model* m, char bgdir[], char varname[], char otname[], int t, char fname[])
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
static void zmodel_adddata(model* m, char* token, char* fname, int line)
{
    char seps[] = " =\t\n";

    if (strcasecmp(token, "MSL") == 0) {
        char mslfname[MAXSTRLEN];
	float** msl = NULL;
        int mvid;
	int nx, ny;

        enkf_printf("    adding custom data \"%s\":\n", token);
        if ((token = strtok(NULL, seps)) == NULL)
            enkf_quit("%s, l.%d: MSL file not specified", fname, line);
        enkf_printf("      file = %s\n", token);
        strcpy(mslfname, token);
        if ((token = strtok(NULL, seps)) == NULL)
            enkf_quit("%s, l.%d: MSL variable not specified", fname, line);
        enkf_printf("      variable = %s\n", token);

        mvid = model_getvarid(m, token);
        assert(mvid >= 0);
	model_getvardims(m, mvid, &nx, &ny, NULL);
	msl = alloc2d(ny, nx, sizeof(float));
	readfield(mslfname, 0, token, msl[0]);

	model_addmodeldata(m, "MSL", ALLOCTYPE_2D, msl);
    } else
        enkf_quit("%s, l.%d: data tag \"%s\" not handled by z-model", fname, line, token);
}

/**
 */
static void zmodel_readfield(model* m, char fname[], int mem, int time, char varname[], int k, float* v)
{
    readfield(fname, k, varname, v);
}

/**
 */
static void zmodel_read3dfield(model* m, char fname[], int mem, int time, char varname[], float* v)
{
    read3dfield(fname, varname, v);
}

/**
 */
static void zmodel_writefield(model* m, char fname[], int time, char varname[], int k, float* v)
{
    writefield(fname, k, varname, v);
}

/**
 */
void zmodel_setup(model* m, char fname[])
{
    model_setgetmemberfname_fn(m, zmodel_getmemberfname);
    model_setgetmemberfnameasync_fn(m, zmodel_getmemberfname_async);
    model_setbgfname_fn(m, zmodel_getbgfname);
    model_setbgfnameasync_fn(m, zmodel_getbgfname_async);
    model_setreadfield_fn(m, zmodel_readfield);
    model_setread3dfield_fn(m, zmodel_read3dfield);
    model_setwritefield_fn(m, zmodel_writefield);
    model_setadddata_fn(m, zmodel_adddata);
}
