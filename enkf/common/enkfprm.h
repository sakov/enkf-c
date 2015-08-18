/******************************************************************************
 *
 * File:        enkfprm.h        
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

#if !defined(_ENKFPRM_H)

typedef struct {
    double z1, z2;
} zint;

typedef struct {
    char* name;
    double x1, x2, y1, y2;
    int nzints;
    zint* zints;
} region;

typedef struct {
    double x1, x2, y1, y2, z1, z2;
} obsdomain;

typedef struct {
    int i, j;
} pointlog;

typedef struct {
    char* obstype;
    double maxbias;
    double maxmad;
    int minnobs;
} badbatchspec;

typedef struct {
    char* fname;
    int mode;
    int scheme;
    int target;
    char* date;

    char* modelprm;
    char* gridprm;
    char* obstypeprm;
    char* obsprm;

    char* ensdir;
    int enssize;
    char* bgdir;

    int nasync;
    char** async_types;
    double* async_timesteps;    /* in days */

    double kfactor;
    double rfactor_base;
    double inflation_base;
    double inf_ratio;
    double locrad;
    int stride;
    int fieldbufsize;
    int nzints;
    zint* zints;
    int nregions;
    region* regions;
    int nplogs;
    pointlog* plogs;
    int sob_stride;
    int nbadbatchspecs;
    badbatchspec* badbatchspecs;
} enkfprm;

enkfprm* enkfprm_read(char fname[]);
void enkfprm_destroy(enkfprm* prm);
void enkfprm_print(enkfprm* prm, char offset[]);
void enkfprm_describeprm(void);

#define _ENKFPRM_H
#endif
