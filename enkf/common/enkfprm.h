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
    char* date;
    double obswindow_min;
    double obswindow_max;

    char* modelprm;
    char* gridprm;
    char* obstypeprm;
    char* obsprm;

    /*
     * Ensemble directory for mode = MODE_ENKF or mode = MODE_ENOI. Directory of
     * the dynamic ensemble for mode = MODE_HYBRID.
     */
    char* ensdir;
    /*
     * directory of the static ensemble for mode = MODE_HYBRID
     */
    char* ensdir2;
    /*
     * total ensemble size
     */
    int enssize;
    /*
     * size of the dynamic ensemble
     */
    int enssize_dynamic;
    /*
     * size of the static ensemble
     */
    int enssize_static;
    /*
     * mixing coefficient for hybrid covariance
     */
    double gamma;
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
    int nzints;
    zint* zints;
    int sob_stride;
    int fieldbufsize;
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
