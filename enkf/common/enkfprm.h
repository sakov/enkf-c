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
    int id;
    int i, j;
    double lon, lat;
    char* gridname;
    int gridid;
} pointlog;

typedef struct {
    char* obstype;
    double maxbias;
    double maxmad;
    int minnobs;
} badbatchspec;

#if !defined(_ENKFPRM_TYPEDEF)
struct enkfprm;
typedef struct enkfprm enkfprm;
#define _ENKFPRM_TYPEDEF
#endif

struct enkfprm {
    char* fname;
    int mode;
    int scheme;
    double alpha;
    char* date;
    double windowmin;
    double windowmax;

    char* modelprm;
    char* gridprm;
    char* obstypeprm;
    char* obsprm;

    char* ensdir;
    int enssize;
    char* bgdir;

    double kfactor;
    /*
     * Unlike other parameters defined in the main parameter file and obstypes
     * parameter file, rfactor_base does not provide the default common value,
     * but a COMMON MULTIPLE for rfactors defined for each observation type.
     */
    double rfactor_base;
    double inflation;
    double inf_ratio;
    int nlocrad;
    double* locrad;
    double* locweight;
    int nlobsmax;
    int stride;
    int fieldbufsize;
    int nregions;
    region* regions;
    int nplogs;
    pointlog* plogs;
    int sob_stride;
    int nbadbatchspecs;
    badbatchspec* badbatchspecs;

    int ncformat;
    int nccompression;
};

enkfprm* enkfprm_read(char fname[]);
void enkfprm_destroy(enkfprm* prm);
void enkfprm_print(enkfprm* prm, char offset[]);
void enkfprm_describeprm(void);

#define _ENKFPRM_H
#endif
