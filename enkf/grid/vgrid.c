/******************************************************************************
 *
 * File:        vgrid.c
 *
 * Created:     25/01/2022
 *
 * Author:      Pavel Sakov
 *              Bureau of Meteorology
 *
 * Purpose:     Code for `vgrid' object (generic vertical grid) and specific
 *              types of vertical grids.
 *
 * Description: Specific grid types are handled via field `gz'.
 *
 * Revisions:
 *
 *****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "ncw.h"
#include "global.h"
#include "ncutils.h"
#include "utils.h"
#include "grid.h"
#include "gridprm.h"
#include "hgrid.h"
#include "vgrid.h"

#define EPS_IJ 1.0e-3

typedef struct {
    int nk;
    int direction;
    int sign_changed;
    double* zt;
    double* zc;
} gz_z;

typedef struct {
    int nk;
    int direction;
    int sign_changed;
    double* zt;
} gz_zt;

/*
 * This structure manages the "new" ROMS vertical coordinate as described in
 * https://www.myroms.org/wiki/Vertical_S-coordinate, Eq. 2 and by Shchepetkin
 * in https://www.myroms.org/forum/viewtopic.php?f=20&t=2189. With hc = 0 it
 * turns into "pure" sigma. To be consistent with the rest of the code we use
 * ct and cc for the stretching function C at the layer centres and corners
 * (these arrays correspond to variables Cs_r and Cs_w in ROMS restarts).
 *
 * Some recent stretching options in ROMS also require stretching of the sigma
 * coordinate, which previously used to be uniform -1 <= s < = 0. This why there
 * are st and sc arrays added.
 */
typedef struct {
    int nk;
    int direction;
    double hc;
    double* ct;
    double* cc;
    double* st;
    double* sc;

    double fij_prev[3];
    double* zt;
    double* zc;

    /*
     * separate cache for grid_fk2z()
     */
    int ij_prev[2];
    double* zt_fk2z;
    double* zc_fk2z;
} gz_sigma;

/*
 * See e.g. Eckerman (2009), Eq. 6 
 * https://journals.ametsoc.org/doi/pdf/10.1175/2008MWR2537.1
 */
typedef struct {
    int nk;
    int direction;
    double* at;
    double* ac;
    double* bt;
    double* bc;
    void* p1;
    void* p2;

    double fij_prev[3];
    double* pt;
    double* pc;

    /*
     * separate cache for grid_fk2z()
     */
    int ij_prev[2];
    double* pt_fk2z;
    double* pc_fk2z;
} gz_hybrid;

/*
 * Structure for vertical grid of general type.
 * Straightforward approach, expensive memory-wise, we will see how it goes.
 * 3D/2D arrays of vertical coordinates for grid centres and corners are read
 * into arrays of floats zt and zc.
 */
typedef struct {
    int ni, nj, nk;
    int direction;
    void* zt;
    void* zc;

    double fij_prev[3];
    double* zt1d;
    double* zc1d;
} gz_numeric;

#define EPSZ 0.001

/** Builds coordinates of layer edges zc[nk + 1] from those of layer centres
 ** zt[nk]. The surface edge is set to 0. All following edges are placed based
 ** on specified layer centres. E.g. if zt = {0.5, 5} then zc = {0, 1, 9}.
 */
static void zt2zc(int nk, int direction, double* zt, double* zc)
{
    int k;

    if (direction == GRIDVDIR_FROMSURF) {
        zc[0] = 0.0;
        for (k = 1; k <= nk; ++k)
            zc[k] = 2.0 * zt[k - 1] - zc[k - 1];
    } else {
        zc[nk] = 0.0;
        for (k = nk - 1; k >= 0; --k)
            zc[k] = 2.0 * zt[k] - zc[k + 1];
    }
}

/** Sets the surface edge when edges were specified on one side of each layer
 ** only (the side away from the surface): it is set to 0 if 0 is on the
 ** surface side of the surface layer centre, or reflected otherwise.
 */
static double zc_surface(double zt_surf, double zc_adjacent)
{
    if (zt_surf * (zc_adjacent - zt_surf) > 0.0)
        return 0.0;
    return 2.0 * zt_surf - zc_adjacent;
}

/** Checks that layer centres and edges are monotonic and consistent.
 */
static void check_zc(char* gname, char* what, int nk, int direction, double* zt, double* zc)
{
    int i;

    if (direction == GRIDVDIR_FROMSURF) {
        if (!(zt[0] <= zt[nk - 1] && zc[0] < zc[nk]))
            enkf_quit("%s: %s: coordinates must increase from surface for VDIR = FROMSURF", gname, what);
    } else {
        if (!(zt[0] >= zt[nk - 1] && zc[0] > zc[nk]))
            enkf_quit("%s: %s: coordinates must decrease towards surface for VDIR = TOSURF", gname, what);
    }
    for (i = 0; i < nk; ++i)
        if (!((zt[i] - zc[i]) * (zc[i + 1] - zt[i]) > 0.0))
            enkf_quit("%s: %s: layer centre %d is not strictly inside its layer edges (zt = %f, zc = %f, %f)", gname, what, i, zt[i], zc[i], zc[i + 1]);
}

/**
 */
static gz_z* gz_z_create(char* gname, int nk, double* z, int nkc, double* zc, char* direction)
{
    gz_z* gz = malloc(sizeof(gz_z));
    int i;

    if (strncasecmp(direction, "FROMSURF", 8) == 0)
        gz->direction = GRIDVDIR_FROMSURF;
    else if (strncasecmp(direction, "TOSURF", 6) == 0)
        gz->direction = GRIDVDIR_TOSURF;
    else
        enkf_quit("programming error");

    gz->zt = z;
    gz->nk = nk;
    gz->sign_changed = 0;
    /*
     * check monotonicity
     */
    for (i = 0; i < nk - 1; ++i)
        if ((z[i + 1] - z[i]) * (z[nk - 1] - z[0]) <= 0.0)
            enkf_quit("%s: non-monotonic Z grid", gname);
    /*
     * change sign if negative
     */
    {
        double sum = 0.0;

        for (i = 0; i < nk; ++i)
            sum += z[i];
        if (sum < 0.0) {
            gz->sign_changed = 1;
            for (i = 0; i < nk; ++i)
                z[i] = -z[i];
        }
    }

    if ((gz->direction == GRIDVDIR_FROMSURF && z[0] == 0) || (gz->direction == GRIDVDIR_TOSURF && z[nk - 1] == 0))
        enkf_quit("%s: the depth of surface layer must be non-zero", gname);

    gz->zc = malloc((nk + 1) * sizeof(double));
    if (zc == NULL)
        zt2zc(nk, gz->direction, z, gz->zc);
    else {
        /*
         * change sign if negative (before using zc)
         */
        {
            double sum = 0.0;

            for (i = 0; i < nkc; ++i)
                sum += zc[i];
            if (sum < 0.0) {
                if (!gz->sign_changed)
                    enkf_quit("%s: signs inconsistent in zt and zc vertical grids", gname);
                for (i = 0; i < nkc; ++i)
                    zc[i] = -zc[i];
            } else if (gz->sign_changed)
                enkf_quit("%s: signs inconsistent in zt and zc vertical grids", gname);
        }

        if (nkc == nk + 1)
            memcpy(gz->zc, zc, nkc * sizeof(double));
        else if (nkc == nk) {
            /*
             * (e.g. in MOM) edges are given on the side of each layer away
             * from the surface (as for NUMERIC grids). For backward
             * compatibility, edges on the surface side of each layer are
             * still accepted for FROMSURF grids.
             */
            if (gz->direction == GRIDVDIR_FROMSURF) {
                if (nk > 1 && (z[0] - zc[0]) * (z[nk - 1] - z[0]) > 0.0) {
                    memcpy(gz->zc, zc, nk * sizeof(double));
                    gz->zc[nk] = 2.0 * z[nk - 1] - zc[nk - 1];
                } else {
                    memcpy(&gz->zc[1], zc, nk * sizeof(double));
                    gz->zc[0] = zc_surface(z[0], zc[0]);
                }
            } else {
                memcpy(gz->zc, zc, nk * sizeof(double));
                gz->zc[nk] = zc_surface(z[nk - 1], zc[nk - 1]);
            }
        } else
            enkf_quit("%s: ZC vertical dimension (%d) must be equal to or larger by one than that of Z (%d)", gname, nkc, nk);
        free(zc);
    }

    check_zc(gname, "Z grid", nk, gz->direction, gz->zt, gz->zc);

    return gz;
}

/**
 */
static gz_zt* gz_zt_create(char* gname, int nk, double* z, char* direction)
{
    gz_zt* gz = malloc(sizeof(gz_zt));
    int i;

    if (strncasecmp(direction, "FROMSURF", 8) == 0)
        gz->direction = GRIDVDIR_FROMSURF;
    else if (strncasecmp(direction, "TOSURF", 6) == 0)
        gz->direction = GRIDVDIR_TOSURF;
    else
        enkf_quit("programming error");

    gz->zt = z;
    gz->nk = nk;
    gz->sign_changed = 0;
    /*
     * check monotonicity
     */
    for (i = 0; i < nk - 1; ++i)
        if ((z[i + 1] - z[i]) * (z[nk - 1] - z[0]) <= 0.0)
            enkf_quit("%s: non-monotonic Z grid", gname);
    /*
     * change sign if negative
     */
    {
        double sum = 0.0;

        for (i = 0; i < nk; ++i)
            sum += z[i];
        if (sum < 0.0) {
            gz->sign_changed = 1;
            for (i = 0; i < nk; ++i)
                z[i] = -z[i];
        }
    }

    if (nk < 2)
        enkf_quit("%s: ZT grid must have at least 2 layers", gname);
    if (gz->direction == GRIDVDIR_FROMSURF && gz->zt[0] > gz->zt[nk - 1])
        enkf_quit("%s: ZT grid: coordinates must increase from surface for VDIR = FROMSURF", gname);
    else if (gz->direction == GRIDVDIR_TOSURF && gz->zt[0] < gz->zt[nk - 1])
        enkf_quit("%s: ZT grid: coordinates must decrease towards surface for VDIR = TOSURF", gname);

    return gz;
}

/**
 */
static gz_sigma* gz_sigma_create(char* gname, int nk, double* st, double* sc, double* ct, double* cc, double hc, char* direction)
{
    gz_sigma* gz = malloc(sizeof(gz_sigma));
    int i;

    gz->nk = nk;
    gz->hc = hc;
    if (strncasecmp(direction, "FROMSURF", 8) == 0)
        gz->direction = GRIDVDIR_FROMSURF;
    else if (strncasecmp(direction, "TOSURF", 6) == 0)
        gz->direction = GRIDVDIR_TOSURF;
    else
        enkf_quit("programming error");

    assert(ct != NULL);
    {
        for (i = 0; i < nk; ++i)
            if (isnan(ct[i]))
                break;
        if (i < nk)
            for (i = 0; i < nk; ++i)
                ct[i] = (0.5 + (double) i) / (double) nk;
    }
    /*
     * check monotonicity
     */
    for (i = 0; i < nk - 2; ++i)
        if ((ct[i + 1] - ct[i]) * (ct[i + 2] - ct[i + 1]) < 0.0)
            enkf_quit("%s: non-monotonic Cs_t", gname);

    gz->ct = ct;

    if (cc != NULL)
        gz->cc = cc;
    else
        gz->cc = malloc((nk + 1) * sizeof(double));
    /*
     * change sign if negative
     */
    {
        double sum = 0.0;

        for (i = 0; i < nk; ++i)
            sum += ct[i];
        if (sum < 0.0) {
            for (i = 0; i < nk; ++i)
                ct[i] = -ct[i];
            if (cc != NULL)
                for (i = 0; i <= nk; ++i)
                    cc[i] = -cc[i];
        }
    }
    /*
     * build Cs_w if necessary
     */
    if (cc == NULL) {
        if (ct[0] > ct[nk - 1]) {
            gz->cc[0] = 1.0;
            gz->cc[nk] = 0.0;
        } else {
            gz->cc[0] = 0.0;
            gz->cc[nk] = 1.0;
        }
        for (i = 1; i < nk; ++i)
            gz->cc[i] = (gz->ct[i - 1] + gz->ct[i]) / 2.0;
    }

    /*
     * read sigma arrays if necessary (assume uniform stretching by default)
     */
    if (st == NULL) {
        if (sc != NULL)
            enkf_quit("%s: SCVARNAME can only be used together with SVARNAME", gname);
        gz->st = malloc(nk * sizeof(double));
        if (gz->direction == GRIDVDIR_FROMSURF)
            for (i = 0; i < nk; ++i)
                gz->st[i] = ((double) i + 0.5) / (double) gz->nk;
        else
            for (i = 0; i < nk; ++i)
                gz->st[i] = 1.0 - ((double) i + 0.5) / (double) gz->nk;
        gz->sc = malloc((nk + 1) * sizeof(double));
        if (gz->direction == GRIDVDIR_FROMSURF)
            for (i = 0; i <= nk; ++i)
                gz->sc[i] = (double) i / (double) gz->nk;
        else
            for (i = 0; i <= nk; ++i)
                gz->sc[i] = 1.0 - (double) i / (double) gz->nk;
    } else {
        double sum = 0.0;

        gz->st = st;
        gz->sc = (sc != NULL) ? sc : malloc((nk + 1) * sizeof(double));
        /*
         * check monotonicity
         */
        for (i = 0; i < nk - 2; ++i)
            if ((st[i + 1] - st[i]) * (st[i + 2] - st[i + 1]) < 0.0)
                enkf_quit("%s: non-monotonic s_rho coordinate (entry SVARNAME in the grid parameter file)", gname);
        /*
         * change sign if negative
         */
        for (i = 0; i < nk; ++i)
            sum += st[i];
        if (sum < 0.0)
            for (i = 0; i < nk; ++i)
                st[i] = -st[i];
    }
    /*
     * build s_w if necessary
     */
    if (sc == NULL) {
        if (gz->st[0] > gz->st[nk - 1]) {
            gz->sc[0] = 1.0;
            gz->sc[nk] = 0.0;
        } else {
            gz->sc[0] = 0.0;
            gz->sc[nk] = 1.0;
        }
        for (i = 1; i < nk; ++i)
            gz->sc[i] = (gz->st[i - 1] + gz->st[i]) / 2.0;
    } else {
        double sum = 0.0;

        /*
         * check monotonicity
         */
        for (i = 0; i < nk - 1; ++i)
            if ((sc[i + 1] - sc[i]) * (sc[i + 2] - sc[i + 1]) < 0.0)
                enkf_quit("%s: non-monotonic s_w coordinate (entry SCVARNAME in the grid parameter file)", gname);
        /*
         * change sign if negative
         */
        for (i = 0; i < nk; ++i)
            sum += sc[i];
        if (sum < 0.0)
            for (i = 0; i <= nk; ++i)
                sc[i] = -sc[i];
    }

    gz->zt = malloc(nk * sizeof(double));
    gz->zc = malloc((nk + 1) * sizeof(double));
    gz->fij_prev[0] = NAN;
    gz->zt_fk2z = malloc(nk * sizeof(double));
    gz->zc_fk2z = malloc((nk + 1) * sizeof(double));
    gz->ij_prev[0] = -1;
    gz->ij_prev[1] = -1;

    check_zc(gname, "SIGMA grid (Cs_r, Cs_w)", nk, gz->direction, gz->ct, gz->cc);
    check_zc(gname, "SIGMA grid (s_rho, s_w)", nk, gz->direction, gz->st, gz->sc);

    return gz;
}

/**
 */
gz_hybrid* gz_hybrid_create(int nk, double* a, double* b, double* ac, double* bc, void* p1, void* p2, char* direction)
{
    gz_hybrid* gz = malloc(sizeof(gz_hybrid));
    int i;

    if (strncasecmp(direction, "FROMSURF", 8) == 0)
        gz->direction = GRIDVDIR_FROMSURF;
    else if (strncasecmp(direction, "TOSURF", 6) == 0)
        gz->direction = GRIDVDIR_TOSURF;
    else
        enkf_quit("programming error");

    gz->nk = nk;
    gz->at = a;
    gz->bt = b;
    gz->p1 = p1;
    gz->p2 = p2;

    if (ac != NULL) {
        assert(bc != NULL);

        gz->ac = ac;
        gz->bc = bc;
    } else {
        assert(bc == NULL);

        gz->ac = malloc((nk + 1) * sizeof(double));
        gz->bc = malloc((nk + 1) * sizeof(double));

        gz->ac[0] = 1.5 * gz->at[0] - 0.5 * gz->at[1];
        if (gz->ac[0] < 0.0)
            gz->ac[0] = 0.0;
        gz->ac[nk] = 1.5 * gz->at[nk - 1] - 0.5 * gz->at[nk - 2];
        if (gz->ac[nk] < 0.0)
            gz->ac[nk] = 0.0;
        for (i = 1; i < nk; ++i)
            gz->ac[i] = (gz->at[i - 1] + gz->at[i]) / 2.0;
        gz->bc[0] = 1.5 * gz->bt[0] - 0.5 * gz->bt[1];
        if (gz->bc[0] < 0.0)
            gz->bc[0] = 0.0;
        gz->bc[nk] = 1.5 * gz->bt[nk - 1] - 0.5 * gz->bt[nk - 2];
        if (gz->bc[nk] < 0.0)
            gz->bc[nk] = 0.0;
        for (i = 1; i < nk; ++i)
            gz->bc[i] = (gz->bt[i - 1] + gz->bt[i]) / 2.0;
    }

    gz->pt = malloc(nk * sizeof(double));
    gz->pc = malloc((nk + 1) * sizeof(double));
    gz->fij_prev[0] = NAN;
    gz->pt_fk2z = malloc(nk * sizeof(double));
    gz->pc_fk2z = malloc((nk + 1) * sizeof(double));
    gz->ij_prev[0] = -1;
    gz->ij_prev[1] = -1;

    return gz;
}

/**
 */
static gz_numeric* gz_numeric_create(int ni, int nj, int nk, char* direction, void* zt, void* zc, int zc_short)
{
    gz_numeric* gz = malloc(sizeof(gz_numeric));
    int nij = (nj > 0) ? ni * nj : ni;
    float* zt0 = (nj > 0) ? ((float***) zt)[0][0] : ((float**) zt)[0];
    float* zc0;
    int ij, k;

    gz->ni = ni;
    gz->nj = nj;
    gz->nk = nk;
    if (strncasecmp(direction, "FROMSURF", 8) == 0)
        gz->direction = GRIDVDIR_FROMSURF;
    else if (strncasecmp(direction, "TOSURF", 6) == 0)
        gz->direction = GRIDVDIR_TOSURF;
    else
        enkf_quit("programming error");
    gz->zt = zt;
    if (zc != NULL) {
        gz->zc = zc;
        zc0 = (nj > 0) ? ((float***) zc)[0][0] : ((float**) zc)[0];
        if (zc_short) {
            /*
             * edges were specified on the side of each layer away from the
             * surface only; set the surface edge
             */
            int ks = (gz->direction == GRIDVDIR_FROMSURF) ? 0 : nk;
            int ka = (gz->direction == GRIDVDIR_FROMSURF) ? 1 : nk - 1;
            int kt = (gz->direction == GRIDVDIR_FROMSURF) ? 0 : nk - 1;

            for (ij = 0; ij < nij; ++ij)
                zc0[ks * nij + ij] = (float) zc_surface(zt0[kt * nij + ij], zc0[ka * nij + ij]);
        }
    } else {
        double* zt1d = malloc(nk * sizeof(double));
        double* zc1d = malloc((nk + 1) * sizeof(double));

        if (nj > 0)
            gz->zc = alloc3d(nk + 1, nj, ni, sizeof(float));
        else
            gz->zc = alloc2d(nk + 1, ni, sizeof(float));
        zc0 = (nj > 0) ? ((float***) gz->zc)[0][0] : ((float**) gz->zc)[0];
        for (ij = 0; ij < nij; ++ij) {
            for (k = 0; k < nk; ++k)
                zt1d[k] = zt0[k * nij + ij];
            zt2zc(nk, gz->direction, zt1d, zc1d);
            for (k = 0; k <= nk; ++k)
                zc0[k * nij + ij] = (float) zc1d[k];
        }
        free(zt1d);
        free(zc1d);
    }

    gz->zt1d = malloc(gz->nk * sizeof(double));
    gz->zc1d = malloc((gz->nk + 1) * sizeof(double));
    gz->fij_prev[0] = NAN;

    return gz;
}

/**
 * @param p - struct gridprm
 * @param g - struct grid
 */
vgrid* vgrid_create(void* p, void* g)
{
    gridprm* prm = p;
    vgrid* vg;
    int ncid;
    size_t nk = 0;

    vg = calloc(1, sizeof(vgrid));
    vg->parent = g;
    vg->type = gridprm_getvtype(prm);

    if (vg->type == GRIDVTYPE_NONE) {
        vg->nk = 1;
        vg->direction = GRIDVDIR_FROMSURF;

        return vg;
    }

    ncw_open(prm->gdatafname, NC_NOWRITE, &ncid);

    if (vg->type == GRIDVTYPE_Z) {
        int varid;
        double* z = NULL;
        double* zc = NULL;
        size_t nkc = 0;

        ncw_inq_varid(ncid, prm->zvarname, &varid);
        ncw_inq_vardims(ncid, varid, 1, NULL, &nk);
        z = malloc(nk * sizeof(double));
        ncu_readvardouble(ncid, varid, nk, z);
        if (prm->zcvarname != NULL) {
            ncw_inq_varid(ncid, prm->zcvarname, &varid);
            ncw_inq_vardims(ncid, varid, 1, NULL, &nkc);
            /*
             * (nkc = nk in MOM)
             */
            if (nkc != nk && nkc != nk + 1)
                enkf_quit("%s: %s: %s: vertical dimension must be equal to or larger by one than that of \"%s\"", prm->gdatafname, prm->name, prm->zcvarname, prm->zvarname);
            zc = malloc(nkc * sizeof(double));
            ncu_readvardouble(ncid, varid, nkc, zc);
        }
        vg->gz = gz_z_create(prm->name, nk, z, nkc, zc, prm->vdirection);
    } else if (vg->type == GRIDVTYPE_ZT) {
        int varid;
        double* z = NULL;

        ncw_inq_varid(ncid, prm->zvarname, &varid);
        ncw_inq_vardims(ncid, varid, 1, NULL, &nk);
        z = malloc(nk * sizeof(double));
        ncu_readvardouble(ncid, varid, nk, z);
        vg->gz = gz_zt_create(prm->name, nk, z, prm->vdirection);
    } else if (vg->type == GRIDVTYPE_SIGMA) {
        int varid;
        double* ct = NULL;
        double* cc = NULL;
        double* st = NULL;
        double* sc = NULL;
        double hc = 0.0;

        ncw_inq_varid(ncid, prm->cvarname, &varid);
        ncw_inq_vardims(ncid, varid, 1, NULL, &nk);
        ct = malloc(nk * sizeof(double));
        ncu_readvardouble(ncid, varid, nk, ct);
        if (prm->ccvarname != NULL) {
            size_t nkc = nk + 1;

            ncw_inq_varid(ncid, prm->ccvarname, &varid);
            ncw_check_vardims(ncid, varid, 1, &nkc);
            cc = malloc(nkc * sizeof(double));
            ncu_readvardouble(ncid, varid, nkc, cc);
        }
        if (prm->hcvarname != NULL) {
            ncw_inq_varid(ncid, prm->hcvarname, &varid);
            ncw_check_varndims(ncid, varid, 0);
            ncu_readvardouble(ncid, varid, 1, &hc);
        }
        if (prm->svarname != NULL) {
            st = malloc(nk * sizeof(double));
            ncw_inq_varid(ncid, prm->svarname, &varid);
            ncw_check_vardims(ncid, varid, 1, &nk);
            ncu_readvardouble(ncid, varid, nk, st);
        }
        if (prm->scvarname != NULL) {
            size_t nkc = nk + 1;

            sc = malloc((nk + 1) * sizeof(double));
            ncw_inq_varid(ncid, prm->scvarname, &varid);
            ncw_check_vardims(ncid, varid, 1, &nkc);
            ncu_readvardouble(ncid, varid, nkc, sc);
        }

        vg->gz = gz_sigma_create(prm->name, nk, st, sc, ct, cc, hc, prm->vdirection);
    } else if (vg->type == GRIDVTYPE_HYBRID) {
        double* a = NULL;
        double* b = NULL;
        double* ac = NULL;
        double* bc = NULL;
        void* p1 = NULL;
        void* p2 = NULL;
        int ni, nj, nij;
        int varid;
        size_t dimlen[2];

        grid_getsize(g, &ni, &nj, NULL);
        if (nj > 0) {
            nij = ni * nj;
            p1 = alloc2d(nj, ni, sizeof(float));
            p2 = alloc2d(nj, ni, sizeof(float));
        } else {
            nij = ni;
            p1 = calloc(ni, sizeof(float));
            p2 = calloc(ni, sizeof(float));
        }

        ncw_inq_varid(ncid, prm->avarname, &varid);
        ncw_inq_vardims(ncid, varid, 1, NULL, &nk);
        a = malloc(nk * sizeof(double));
        ncu_readvardouble(ncid, varid, nk, a);

        ncw_inq_varid(ncid, prm->bvarname, &varid);
        ncw_check_vardims(ncid, varid, 1, &nk);
        b = malloc(nk * sizeof(double));
        ncu_readvardouble(ncid, varid, nk, b);

        if (prm->acvarname != NULL) {
            size_t nkc = nk + 1;

            if (prm->bcvarname == NULL)
                enkf_quit("%s: %s: BCVARNAME must be specified together with ACVARNAME", prm->gdatafname, prm->name);

            ncw_inq_varid(ncid, prm->acvarname, &varid);
            ncw_check_vardims(ncid, varid, 1, &nkc);
            ac = malloc(nkc * sizeof(double));
            ncu_readvardouble(ncid, varid, nkc, ac);

            ncw_inq_varid(ncid, prm->bcvarname, &varid);
            ncw_check_vardims(ncid, varid, 1, &nkc);
            bc = malloc(nkc * sizeof(double));
            ncu_readvardouble(ncid, varid, nkc, bc);
        }

        if (prm->acvarname == NULL && prm->bcvarname != NULL)
            enkf_quit("%s: %s: ACVARNAME must be specified together with BCVARNAME", prm->gdatafname, prm->name);

        ncw_inq_varid(ncid, prm->p1varname, &varid);
        if (nj > 0) {
            dimlen[0] = nj;
            dimlen[1] = ni;
        } else
            dimlen[0] = ni;
        ncw_check_vardims(ncid, varid, (nj > 0) ? 2 : 1, dimlen);
        ncu_readvarfloat(ncid, varid, nij, (nj > 0) ? ((float**) p1)[0] : p1);

        ncw_inq_varid(ncid, prm->p2varname, &varid);
        ncw_check_vardims(ncid, varid, (nj > 0) ? 2 : 1, dimlen);
        ncu_readvarfloat(ncid, varid, nij, (nj > 0) ? ((float**) p2)[0] : p2);

        vg->gz = gz_hybrid_create(nk, a, b, ac, bc, p1, p2, prm->vdirection);
    } else if (vg->type == GRIDVTYPE_NUMERIC) {
        int varid;
        int ndimt, ndimc;
        size_t dimlent[3], dimlenc[3];
        int ni = 0, nj = 0;
        void* zt = NULL;
        void* zc = NULL;
        int kc_offset = 0;
        int zc_short = 0;

        /*
         * zt
         */
        ncw_inq_varid(ncid, prm->zvarname, &varid);
        if (ncw_var_hasunlimdim(ncid, varid))
            enkf_quit("%s: %s: %s: can not have unlimited dimension", prm->gdatafname, prm->name, prm->zvarname);
        ncw_inq_vardims(ncid, varid, 3, &ndimt, dimlent);
        nk = dimlent[0];
        if (ndimt == 3) {
            nj = dimlent[1];
            ni = dimlent[2];
            zt = alloc3d(dimlent[0], dimlent[1], dimlent[2], sizeof(float));
            ncw_get_var_float(ncid, varid, ((float***) zt)[0][0]);
        } else if (ndimt == 2) {
            nj = 0;
            ni = dimlent[1];
            zt = alloc2d(dimlent[0], dimlent[1], sizeof(float));
            ncw_get_var_float(ncid, varid, ((float**) zt)[0]);
        } else
            enkf_quit("%s: %s: %s: can not have single dimension", prm->gdatafname, prm->name, prm->zvarname);

        {
            int ni2, nj2;

            grid_getsize(g, &ni2, &nj2, NULL);
            if (ni != ni2 || nj != nj2)
                enkf_quit("%s: %s: horizontal dimension of variable \"%s\" must be equal to that of the grid", prm->gdatafname, prm->name, prm->zvarname);
        }

        /*
         * zc
         */
        if (prm->zcvarname != NULL) {
            ncw_inq_varid(ncid, prm->zcvarname, &varid);
            if (ncw_var_hasunlimdim(ncid, varid))
                enkf_quit("%s: %s: %s: can not have unlimited dimension", prm->gdatafname, prm->name, prm->zcvarname);
            ncw_inq_vardims(ncid, varid, 3, &ndimc, dimlenc);
            if (dimlenc[0] == nk) {
                /*
                 * edges on the side of each layer away from the surface
                 */
                dimlenc[0]++;
                zc_short = 1;
                if (strncasecmp(prm->vdirection, "FROMSURF", 8) == 0)
                    kc_offset = 1;
            }
            if (dimlenc[0] != nk + 1)
                enkf_quit("%s: %s: %s: vertical dimension must be equal or larger by one than that of \"%s\"", prm->gdatafname, prm->name, prm->zcvarname, prm->zvarname);
            if (ndimc != ndimt)
                enkf_quit("%s: %s: \"%s\" and \"%s\" must have equal number of dimensions", prm->gdatafname, prm->name, prm->zcvarname, prm->zvarname);
            if (dimlenc[1] != dimlent[1] || (ndimc == 3 && dimlenc[2] != dimlent[2]))
                enkf_quit("%s: %s: \"%s\" and \"%s\" must have equal horizontal dimensions", prm->gdatafname, prm->name, prm->zcvarname, prm->zvarname);
            if (ndimc == 3) {
                zc = alloc3d(dimlenc[0], dimlenc[1], dimlenc[2], sizeof(float));
                ncw_get_var_float(ncid, varid, ((float***) zc)[kc_offset][0]);
            } else if (ndimc == 2) {
                zc = alloc2d(dimlenc[0], dimlenc[1], sizeof(float));
                ncw_get_var_float(ncid, varid, ((float**) zc)[kc_offset]);
            } else
                enkf_quit("%s: %s: %s: can not have single dimension", prm->gdatafname, prm->name, prm->zcvarname);
        }

        vg->gz = gz_numeric_create(ni, nj, nk, prm->vdirection, zt, zc, zc_short);
    } else
        enkf_quit("programming error");

    ncw_close(ncid);

    vg->nk = nk;
    if (strncasecmp(prm->vdirection, "FROMSURF", 8) == 0)
        vg->direction = GRIDVDIR_FROMSURF;
    else if (strncasecmp(prm->vdirection, "TOSURF", 6) == 0)
        vg->direction = GRIDVDIR_TOSURF;
    else
        enkf_quit("programming error");

    return vg;
}

/**
 */
static void gz_z_destroy(gz_z* gz)
{
    free(gz->zt);
    free(gz->zc);
    free(gz);
}

/**
 */
static void gz_zt_destroy(gz_zt* gz)
{
    free(gz->zt);
    free(gz);
}

/**
 */
static void gz_sigma_destroy(gz_sigma* gz)
{
    free(gz->ct);
    free(gz->cc);
    free(gz->zt);
    free(gz->zc);
    free(gz->st);
    free(gz->sc);
    free(gz->zt_fk2z);
    free(gz->zc_fk2z);
    free(gz);
}

/**
 */
static void gz_hybrid_destroy(gz_hybrid* gz)
{
    free(gz->at);
    free(gz->bt);
    free(gz->ac);
    free(gz->bc);
    free(gz->p1);
    free(gz->p2);
    free(gz->pt);
    free(gz->pc);
    free(gz->pt_fk2z);
    free(gz->pc_fk2z);
    free(gz);
}

/**
 */
static void gz_numeric_destroy(gz_numeric* gz)
{
    free(gz->zt);
    free(gz->zc);
    free(gz->zt1d);
    free(gz->zc1d);
    free(gz);
}

/**
 */
void vgrid_destroy(vgrid* vg)
{
    if (vg == NULL)
        return;
    if (vg->type == GRIDVTYPE_NONE);
    else if (vg->type == GRIDVTYPE_Z)
        gz_z_destroy(vg->gz);
    else if (vg->type == GRIDVTYPE_ZT)
        gz_zt_destroy(vg->gz);
    else if (vg->type == GRIDVTYPE_SIGMA)
        gz_sigma_destroy(vg->gz);
    else if (vg->type == GRIDVTYPE_HYBRID)
        gz_hybrid_destroy(vg->gz);
    else if (vg->type == GRIDVTYPE_NUMERIC)
        gz_numeric_destroy(vg->gz);
    else
        enkf_quit("programming error");
    free(vg);
}

/**
 */
void vgrid_describe(vgrid* vg, char* offset)
{
    switch (vg->type) {
    case GRIDVTYPE_NONE:
        enkf_printf("%s  v type = NONE\n", offset);
        break;
    case GRIDVTYPE_Z:
        enkf_printf("%s  v type = Z\n", offset);
        break;
    case GRIDVTYPE_ZT:
        enkf_printf("%s  v type = ZT\n", offset);
        break;
    case GRIDVTYPE_SIGMA:
        enkf_printf("%s  v type = SIGMA\n", offset);
        enkf_printf("%s  hc = %f\n", offset, ((gz_sigma*) vg->gz)->hc);
        break;
    case GRIDVTYPE_HYBRID:
        enkf_printf("%s  v type = HYBRID\n", offset);
        break;
    case GRIDVTYPE_NUMERIC:
        enkf_printf("%s  v type = NUMERIC\n", offset);
        break;
    default:
        enkf_printf("%s  v type = UNDEFINED\n", offset);
    }
    if (vg->type != GRIDVTYPE_UNDEFINED)
        enkf_printf("%s  v dir = %s\n", offset, (vg->direction == GRIDVDIR_FROMSURF) ? "FROMSURF" : "TOSURF");
}

/** Maps vertical coordinate to fractional cell index.
 * Both zt and zc are assumed monotonous. zt[i] maps to (double) i, zc[i] maps
 * to (i - 0.5) for ascending arrays, and to (i + 0.5) for descending arrays.
 * More generally, the fractional indices are calculated by linear interpolation
 * of indices in the interval betweenn the points in arrays zt and zc closest to
 * z. z values above surface are considered to be at surface, z values below
 * bottom are mapped to NaN.
 * @param n Number of layers
 * @param zt Coordinates of layer centres [n]
 * @param zc Coordinates of layer edges [n + 1]
 * @param z Depth/height
 * @return Fractional layer index
 */
static double z2fk_basic(int n, double* zt, double* zc, double z)
{
    int ascending, i1, i2, imid;

    ascending = (zc[n] > zc[0]) ? 1 : 0;

    if (ascending) {
        if (z <= zc[0])
            return -0.5;
        if (z > zc[n])
            return NAN;
    } else {
        if (z <= zc[n])
            return (double) n - 0.5;
        if (z > zc[0])
            return NAN;
    }

    i1 = 0;
    i2 = n;
    if (ascending) {
        while (1) {
            imid = (i1 + i2) / 2;
            if (imid == i1)
                break;
            if (z > zc[imid])
                i1 = imid;
            else
                i2 = imid;
        }
        if (z > zc[i2])
            i1++;
        if (z < zt[i1])
            return (double) i1 + 0.5 * (z - zt[i1]) / (zt[i1] - zc[i1]);
        else
            return (double) i1 + 0.5 * (z - zt[i1]) / (zc[i1 + 1] - zt[i1]);
    } else {
        while (1) {
            imid = (i1 + i2) / 2;
            if (imid == i1)
                break;
            if (z > zc[imid])
                i2 = imid;
            else
                i1 = imid;
        }
        if (z < zc[i2])
            i1++;
        if (z > zt[i1])
            return (double) i1 + 0.5 * (z - zt[i1]) / (zt[i1] - zc[i1]);
        else
            return (double) i1 + 0.5 * (z - zt[i1]) / (zc[i1 + 1] - zt[i1]);
    }
}

/**
 */
static void gz_z_z2fk(vgrid* vg, double z, double* fk)
{
    gz_z* gz = vg->gz;

    *fk = z2fk_basic(gz->nk, gz->zt, gz->zc, z);
}

/**
 */
static void gz_zt_z2fk(vgrid* vg, double z, double* fk)
{
    gz_zt* gz = vg->gz;
    int n = gz->nk;
    double* zt = gz->zt;

    int ascending = (zt[n - 1] > zt[0]) ? 1 : 0;
    int i1 = 0;
    int i2 = n - 1;
    int imid;

    if (ascending) {
        if (z <= zt[0])
            *fk = (z - zt[0]) / (zt[1] - zt[0]);
        else if (z >= zt[n - 1])
            *fk = n - 1 + (z - zt[n - 1]) / (zt[n - 1] - zt[n - 2]);
        else {
            while (1) {
                imid = (i1 + i2) / 2;
                if (imid == i1)
                    break;
                if (z > zt[imid])
                    i1 = imid;
                else
                    i2 = imid;
            }
            *fk = i1 + (z - zt[i1]) / (zt[i1 + 1] - zt[i1]);
        }
    } else {
        if (z >= zt[0])
            *fk = (z - zt[0]) / (zt[1] - zt[0]);
        else if (z <= zt[n - 1])
            *fk = n - 1 + (z - zt[n - 1]) / (zt[n - 1] - zt[n - 2]);
        else {
            while (1) {
                imid = (i1 + i2) / 2;
                if (imid == i1)
                    break;
                if (z > zt[imid])
                    i2 = imid;
                else
                    i1 = imid;
            }
            *fk = i1 + (z - zt[i1]) / (zt[i1 + 1] - zt[i1]);
        }
    }
}

/** Calculates depths of layer centres and edges for a sigma grid.
 */
static void gz_sigma_calcz(gz_sigma* gz, double h, double* zt, double* zc)
{
    int i;

    if (gz->hc != 0.0) {
        for (i = 0; i < gz->nk; ++i)
            zt[i] = h * (gz->hc * gz->st[i] + h * gz->ct[i]) / (gz->hc + h);
        for (i = 0; i <= gz->nk; ++i)
            zc[i] = h * (gz->hc * gz->sc[i] + h * gz->cc[i]) / (gz->hc + h);
    } else {
        for (i = 0; i < gz->nk; ++i)
            zt[i] = h * gz->ct[i];
        for (i = 0; i <= gz->nk; ++i)
            zc[i] = h * gz->cc[i];
    }
}

/** Calculates pressure of layer centres and edges for a hybrid grid.
 */
static void gz_hybrid_calcp(gz_hybrid* gz, double p1, double p2, double* pt, double* pc)
{
    int i;

    for (i = 0; i < gz->nk; ++i)
        pt[i] = gz->at[i] + gz->bt[i] * (p1 - p2);
    for (i = 0; i <= gz->nk; ++i)
        pc[i] = gz->ac[i] + gz->bc[i] * (p1 - p2);
}

/**
 */
static void gz_sigma_z2fk(vgrid* vg, double* fij, double z, double* fk)
{
    gz_sigma* gz = vg->gz;
    void* g = vg->parent;
    void* depth = grid_getdepth(g);
    int htype = grid_gethtype(g);

    if (depth == NULL)
        enkf_quit("%s: DEPTHVARNAME must be entered for SIGMA grids to assimilate subsurface obs", grid_getname(g));

    if (htype == GRIDHTYPE_RECTANGULAR || htype == GRIDHTYPE_CURVILINEAR || htype == GRIDHTYPE_CURVILINEAR2) {
        if (isnan(gz->fij_prev[0]) || fabs(fij[0] - gz->fij_prev[0]) > EPS_IJ || fabs(fij[1] - gz->fij_prev[1]) > EPS_IJ) {
            double h = (double) grid_interpolate2d(g, fij, depth);

            gz_sigma_calcz(gz, h, gz->zt, gz->zc);
            gz->fij_prev[0] = fij[0];
            gz->fij_prev[1] = fij[1];
        }
    } else if (htype == GRIDHTYPE_UNSTRUCTURED) {
        if (isnan(gz->fij_prev[0]) || fabs(fij[0] - gz->fij_prev[0]) > EPS_IJ || fabs(fij[1] - gz->fij_prev[1]) > EPS_IJ || fabs(fij[2] - gz->fij_prev[2]) > EPS_IJ) {
            double h = (double) grid_interpolate2d(g, fij, depth);

            gz_sigma_calcz(gz, h, gz->zt, gz->zc);
            gz->fij_prev[0] = fij[0];
            gz->fij_prev[1] = fij[1];
            gz->fij_prev[2] = fij[2];
        }
    }

    *fk = z2fk_basic(gz->nk, gz->zt, gz->zc, z);
}

/**
 */
static void gz_hybrid_z2fk(vgrid* vg, double* fij, double z, double* fk)
{
    void* g = vg->parent;
    gz_hybrid* gz = vg->gz;
    int htype = grid_gethtype(g);

    if (htype == GRIDHTYPE_RECTANGULAR || htype == GRIDHTYPE_CURVILINEAR || htype == GRIDHTYPE_CURVILINEAR2) {
        if (isnan(gz->fij_prev[0]) || fabs(fij[0] - gz->fij_prev[0]) > EPS_IJ || fabs(fij[1] - gz->fij_prev[1]) > EPS_IJ) {
            double p1 = grid_interpolate2d(g, fij, gz->p1);
            double p2 = grid_interpolate2d(g, fij, gz->p2);

            gz_hybrid_calcp(gz, p1, p2, gz->pt, gz->pc);
            gz->fij_prev[0] = fij[0];
            gz->fij_prev[1] = fij[1];
        }
    } else if (htype == GRIDHTYPE_UNSTRUCTURED) {
        if (isnan(gz->fij_prev[0]) || fabs(fij[0] - gz->fij_prev[0]) > EPS_IJ || fabs(fij[1] - gz->fij_prev[1]) > EPS_IJ || fabs(fij[2] - gz->fij_prev[2]) > EPS_IJ) {
            double p1 = grid_interpolate2d(g, fij, gz->p1);
            double p2 = grid_interpolate2d(g, fij, gz->p2);

            gz_hybrid_calcp(gz, p1, p2, gz->pt, gz->pc);
            gz->fij_prev[0] = fij[0];
            gz->fij_prev[1] = fij[1];
            gz->fij_prev[2] = fij[2];
        }
    }

    *fk = z2fk_basic(gz->nk, gz->pt, gz->pc, z);
}

/**
 */
static void gz_numeric_z2fk(vgrid* vg, double* fij, double z, double* fk)
{
    gz_numeric* gz = vg->gz;
    void* g = vg->parent;
    int htype = grid_gethtype(g);

    if (htype == GRIDHTYPE_RECTANGULAR || htype == GRIDHTYPE_CURVILINEAR || htype == GRIDHTYPE_CURVILINEAR2) {
        if (isnan(gz->fij_prev[0]) || fabs(fij[0] - gz->fij_prev[0]) > EPS_IJ || fabs(fij[1] - gz->fij_prev[1]) > EPS_IJ) {
            int k;

            for (k = 0; k < gz->nk; ++k)
                gz->zt1d[k] = grid_interpolate2d(g, fij, ((float***) gz->zt)[k]);
            for (k = 0; k <= gz->nk; ++k)
                gz->zc1d[k] = grid_interpolate2d(g, fij, ((float***) gz->zc)[k]);
            gz->fij_prev[0] = fij[0];
            gz->fij_prev[1] = fij[1];
        }
    } else if (htype == GRIDHTYPE_UNSTRUCTURED) {
        if (isnan(gz->fij_prev[0]) || fabs(fij[0] - gz->fij_prev[0]) > EPS_IJ || fabs(fij[1] - gz->fij_prev[1]) > EPS_IJ || fabs(fij[2] - gz->fij_prev[2]) > EPS_IJ) {
            int k;

            for (k = 0; k < gz->nk; ++k)
                gz->zt1d[k] = grid_interpolate2d(g, fij, ((float**) gz->zt)[k]);
            for (k = 0; k <= gz->nk; ++k)
                gz->zc1d[k] = grid_interpolate2d(g, fij, ((float**) gz->zc)[k]);
            gz->fij_prev[0] = fij[0];
            gz->fij_prev[1] = fij[1];
            gz->fij_prev[2] = fij[2];
        }
    }

    *fk = z2fk_basic(gz->nk, gz->zt1d, gz->zc1d, z);
}

/** This procedure is called from grid_z2fk(), which determines and returns the
 ** status. 
 */
void vgrid_z2fk(vgrid* vg, double* fij, double z, double* fk)
{
    if (vg->type == GRIDVTYPE_Z)
        gz_z_z2fk(vg, z, fk);
    else if (vg->type == GRIDVTYPE_ZT)
        gz_zt_z2fk(vg, z, fk);
    else if (vg->type == GRIDVTYPE_SIGMA)
        gz_sigma_z2fk(vg, fij, z, fk);
    else if (vg->type == GRIDVTYPE_HYBRID)
        gz_hybrid_z2fk(vg, fij, z, fk);
    else if (vg->type == GRIDVTYPE_NUMERIC)
        gz_numeric_z2fk(vg, fij, z, fk);
    else if (vg->type == GRIDVTYPE_NONE)
        *fk = INFINITY;         /* not NAN to pass the check in grid_z2fk() */
    else
        enkf_quit("programming error");
}

/** Maps vertical coordinate to fractional layer index at horizontal node `ij'
 ** (ij = j * ni + i for structured grids). Used for calculating the number of
 ** levels from depth; not applicable to SIGMA and HYBRID grids.
 */
void vgrid_z2fk_node(vgrid* vg, int ij, double z, double* fk)
{
    if (vg->type == GRIDVTYPE_NUMERIC) {
        gz_numeric* gz = vg->gz;
        int k;

        for (k = 0; k < gz->nk; ++k)
            gz->zt1d[k] = (gz->nj > 0) ? ((float***) gz->zt)[k][0][ij] : ((float**) gz->zt)[k][ij];
        for (k = 0; k <= gz->nk; ++k)
            gz->zc1d[k] = (gz->nj > 0) ? ((float***) gz->zc)[k][0][ij] : ((float**) gz->zc)[k][ij];
        gz->fij_prev[0] = NAN;  /* invalidate cache */
        *fk = z2fk_basic(gz->nk, gz->zt1d, gz->zc1d, z);
    } else if (vg->type == GRIDVTYPE_Z || vg->type == GRIDVTYPE_ZT || vg->type == GRIDVTYPE_NONE)
        vgrid_z2fk(vg, NULL, z, fk);
    else
        enkf_quit("programming error");
}

/** Maps fractional index to vertical coordinate.
 * @param n Number of layers
 * @param zt Coordinates of layer centres [n]
 * @param zc Coordinates of layer edges [n + 1]
 * @param fk Fractional layer index
 * @return z
 */
#if defined(ENKF_PREP) || defined(ENKF_CALC)
static double fk2z_basic(int n, double* zt, double* zc, double fk)
{
    int k;
    double dk;

    fk += 0.5;

    if (fk <= 0.0)
        return zc[0];
    if (fk >= n)
        return zc[n];

    k = (int) floor(fk);
    dk = fk - (double) k;
    if (dk < 0.5)
        return zc[k] + dk * (zt[k] - zc[k]) / 0.5;
    else
        return zt[k] + (dk - 0.5) * (zc[k + 1] - zt[k]) / 0.5;
}

static double gzt_fk2z(int n, double* zt, double fk)
{
    int k;

    if (fk <= 0.0)
        return zt[0] + fk * (zt[1] - zt[0]);
    else if (fk >= (double) (n - 1))
        return zt[n - 1] + (fk - (double) (n - 1)) * (zt[n - 1] - zt[n - 2]);
    k = fk;
    return zt[k] + (fk - (double) k) * (zt[k + 1] - zt[k]);
}

/**
 */
int grid_fk2z(grid* g, int* ij, double fk, double* z)
{
    int htype = grid_gethtype(g);
    int vtype = grid_getvtype(g);

    if (htype == GRIDHTYPE_RECTANGULAR || htype == GRIDHTYPE_CURVILINEAR || htype == GRIDHTYPE_CURVILINEAR2) {
        int i = ij[0];
        int j = ij[1];
        int ni, nj;

        grid_getsize(g, &ni, &nj, NULL);
        if (i < 0 || j < 0 || i >= ni || j >= nj) {
            *z = NAN;
            return STATUS_OUTSIDEGRID;
        }

        if (vtype == GRIDVTYPE_Z) {
            gz_z* gz = ((vgrid*) grid_getvgrid(g))->gz;

            *z = fk2z_basic(gz->nk, gz->zt, gz->zc, fk);
        } else if (vtype == GRIDVTYPE_ZT) {
            gz_zt* gz = ((vgrid*) grid_getvgrid(g))->gz;

            *z = gzt_fk2z(gz->nk, gz->zt, fk);
        } else if (vtype == GRIDVTYPE_SIGMA) {
            gz_sigma* gz = ((vgrid*) grid_getvgrid(g))->gz;
            float** depth = grid_getdepth(g);

            if (depth == NULL)
                enkf_quit("%s: DEPTHVARNAME must be entered for SIGMA grids to assimilate subsurface obs", grid_getname(g));

            if (gz->ij_prev[0] != i || gz->ij_prev[1] != j) {
                gz_sigma_calcz(gz, depth[j][i], gz->zt_fk2z, gz->zc_fk2z);
                gz->ij_prev[0] = i;
                gz->ij_prev[1] = j;
            }

            *z = fk2z_basic(gz->nk, gz->zt_fk2z, gz->zc_fk2z, fk);
        } else if (vtype == GRIDVTYPE_HYBRID) {
            gz_hybrid* gz = ((vgrid*) grid_getvgrid(g))->gz;

            if (gz->ij_prev[0] != i || gz->ij_prev[1] != j) {
                gz_hybrid_calcp(gz, ((float**) gz->p1)[j][i], ((float**) gz->p2)[j][i], gz->pt_fk2z, gz->pc_fk2z);
                gz->ij_prev[0] = i;
                gz->ij_prev[1] = j;
            }

            *z = fk2z_basic(gz->nk, gz->pt_fk2z, gz->pc_fk2z, fk);
        } else
            enkf_quit("programming error");
    } else if (htype == GRIDHTYPE_UNSTRUCTURED) {
        int i = ij[0];
        int n;

        grid_getsize(g, &n, NULL, NULL);
        if (i < 0 || i >= n) {
            *z = NAN;
            return STATUS_OUTSIDEGRID;
        }
        if (vtype == GRIDVTYPE_Z) {
            gz_z* gz = ((vgrid*) grid_getvgrid(g))->gz;

            *z = fk2z_basic(gz->nk, gz->zt, gz->zc, fk);
        } else if (vtype == GRIDVTYPE_SIGMA) {
            gz_sigma* gz = ((vgrid*) grid_getvgrid(g))->gz;
            float* depth = grid_getdepth(g);

            if (depth == NULL)
                enkf_quit("%s: DEPTHVARNAME must be entered for SIGMA grids to assimilate subsurface obs", grid_getname(g));

            if (gz->ij_prev[0] != i || gz->ij_prev[1] != -1) {
                gz_sigma_calcz(gz, depth[i], gz->zt_fk2z, gz->zc_fk2z);
                gz->ij_prev[0] = i;
                gz->ij_prev[1] = -1;
            }

            *z = fk2z_basic(gz->nk, gz->zt_fk2z, gz->zc_fk2z, fk);
        } else if (vtype == GRIDVTYPE_HYBRID) {
            gz_hybrid* gz = ((vgrid*) grid_getvgrid(g))->gz;

            if (gz->ij_prev[0] != i || gz->ij_prev[1] != -1) {
                gz_hybrid_calcp(gz, ((float*) gz->p1)[i], ((float*) gz->p2)[i], gz->pt_fk2z, gz->pc_fk2z);
                gz->ij_prev[0] = i;
                gz->ij_prev[1] = -1;
            }

            *z = fk2z_basic(gz->nk, gz->pt_fk2z, gz->pc_fk2z, fk);
        } else
            enkf_quit("programming error");
    }

    return STATUS_OK;
}
#endif

/**
 */
int vgrid_signchanged(vgrid* vg)
{
    if (vg->type == GRIDVTYPE_Z)
        return ((gz_z*) vg->gz)->sign_changed;
    else if (vg->type == GRIDVTYPE_ZT)
        return ((gz_zt*) vg->gz)->sign_changed;
    return 0;
}
