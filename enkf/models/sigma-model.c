/******************************************************************************
 *
 * File:        sigma-model.c        
 *
 * Created:     03/12/2014
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
#if !defined(NO_GRIDUTILS)
#include <gridnodes.h>
#endif
#include <assert.h>
#include "ncw.h"
#include "definitions.h"
#include "utils.h"
#include "grid.h"
#include "enkfprm.h"
#include "model.h"
#include "standard-model.h"
#include "z-model.h"

/**
 */
void smodel_setgrids(model* m, char gfname[])
{
    int ngrid = 0;
    gridprm* prm = NULL;
    int gid;

    gridprm_create(gfname, &ngrid, &prm, "MASKVARNAME");
    assert(ngrid > 0);

    for (gid = 0; gid < ngrid; ++gid) {
        grid* g = NULL;
        int** numlevels = NULL;
        int nx, ny, nz;
        int i, j;

        g = grid_create(&prm[gid], gid, GRIDVTYPE_SIGMA);
        grid_settocartesian_fn(g, ll2xyz);
        model_setgrid(m, g);

        numlevels = grid_getnumlevels(g);
        grid_getdims(g, &nx, &ny, &nz);
        for (j = 0; j < ny; ++j)
            for (i = 0; i < nx; ++i)
                numlevels[j][i] *= nz;
    }

    gridprm_destroy(ngrid, prm);
}

/**
 */
static void smodel_adddata(model* m, char* token, char* fname, int line)
{
    if (strcasecmp(token, "MSL") != 0)
        enkf_quit("%s, l.%d: data tag \"%s\" not handled by sigma-model", fname, line, token);

    standardmodel_adddata_2D(m, "MSL", fname, line);
}

/**
 */
void smodel_setup(model* m, char fname[])
{
    model_setgetmemberfname_fn(m, standardmodel_getmemberfname);
    model_setgetmemberfnameasync_fn(m, standardmodel_getmemberfname_async);
    model_setbgfname_fn(m, standardmodel_getbgfname);
    model_setbgfnameasync_fn(m, standardmodel_getbgfname_async);
    model_setreadfield_fn(m, standardmodel_readfield);
    model_setread3dfield_fn(m, standardmodel_read3dfield);
    model_setwritefield_fn(m, standardmodel_writefield);
    model_setadddata_fn(m, smodel_adddata);
}
