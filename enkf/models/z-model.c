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
#if !defined(NO_GRIDUTILS)
#include <gridnodes.h>
#endif
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
void zmodel_setgrids(model* m, char gfname[])
{
    int ngrid = 0;
    gridprm* prm = NULL;
    int i;

    gridprm_create(gfname, &ngrid, &prm, "NUMLEVELSVARNAME");
    assert(ngrid > 0);

    for (i = 0; i < ngrid; ++i) {
        grid* g = NULL;

        g = grid_create(&prm[i], i, GRIDVTYPE_Z);
        grid_settocartesian_fn(g, ll2xyz);
        model_setgrid(m, g);
    }

    gridprm_destroy(ngrid, prm);
}

/**
 */
static void zmodel_adddata(model* m, char* token, char* fname, int line)
{
    if (strcasecmp(token, "MSL") == 0)
        standardmodel_adddata_2D(m, token, fname, line);
    else
        enkf_quit("%s, l.%d: data tag \"%s\" not handled by z-model", fname, line, token);
}

/**
 */
void zmodel_setup(model* m, char fname[])
{
    model_setgetmemberfname_fn(m, standardmodel_getmemberfname);
    model_setgetmemberfnameasync_fn(m, standardmodel_getmemberfname_async);
    model_setbgfname_fn(m, standardmodel_getbgfname);
    model_setbgfnameasync_fn(m, standardmodel_getbgfname_async);
    model_setreadfield_fn(m, standardmodel_readfield);
    model_setread3dfield_fn(m, standardmodel_read3dfield);
    model_setwritefield_fn(m, standardmodel_writefield);
    model_setadddata_fn(m, zmodel_adddata);
}
