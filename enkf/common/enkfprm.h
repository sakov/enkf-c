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
    char* name;
    double x1, x2, y1, y2;
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
    char* obsspec;
    char* date;
    int nasync;
    char** async_types;
    double* async_timesteps;    /* in days */
    char* modelprm;
    char* gridprm;
    char* ensdir;
    int enssize;
    char* bgdir;
    int nvar;
    char** varnames;
    double* inflations;
    int ntypes;
    char** types;
    char** typevars;
    char** hfunctions;
    double* rfactors;
    obsdomain* obsdomains;
    char* msl_fname;
    char* msl_varname;
    double kfactor;
    double rfactor_base;
    double inflation_base;
    double locrad;
    int stride;
    int fieldbufsize;
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
void enkfprm_describe(void);

#define _ENKFPRM_H
#endif
