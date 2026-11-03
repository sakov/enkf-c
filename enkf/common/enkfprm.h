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

struct enkfprm {
    char* fname;
    int mode;
    int scheme;
    double alpha;
    char* time;                 /* analysis time */
    double obswindow_min;
    double obswindow_max;

    char* modelprm;
    char* gridprm;
    char* obstypeprm;
    char* obsprm;

    /*
     * directory of the dynamic ensemble
     */
    char* ensdir_dynamic;
    /*
     * size of the dynamic ensemble
     */
    int enssize_dynamic;
    /*
     * number of static ensembles
     */
    int nstatic;
    /*
     * directories of static ensembles
     */
    char** ensdir_static;
    /*
     * sizes of static ensembles
     */
    int* enssize_static;
    /*
     * scaling coefficients for static covariances
     */
    double* gamma;
    /*
     * adaptive QC
     */
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
    int nzints;
    zint* zints;
    int sob_stride;
    int fieldbufsize;
    int nfieldsplit;
    int nregions;
    region* regions;
    int nplog;
    pointlog* plogs;
    int nbadbatchspecs;
    badbatchspec* badbatchspecs;

    int ncformat;
    int nccompression;
};

enkfprm* enkfprm_read(char fname[]);
void enkfprm_destroy(enkfprm* prm);
void enkfprm_print(enkfprm* prm, char offset[]);
void enkfprm_describeprm(void);
void enkfprm_describeprm_ensdiag(void);

#define _ENKFPRM_H
#endif
