/******************************************************************************
 *
 * File:        calcs.h        
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

#if !defined(_CALCS_H)

double traceprod(int transposeA, int transposeB, int m, int n, double** A, double** B, int iscompact);
void calc_G(int m, int p, double** M, double** S, double** G);
void calc_w(int m, int p, double** G, double* s, double* w);
void calc_wT_denkf(int m, int mout, int p, double* s, double** S, double** Sa, double** M, double** G, double* w, double** T);
void calc_wT_etkf(int m, int mout, int p, double* s, double** S, double** Sa, double** M, double** G, double* w, double** T);

/*
 * redundant
 */
void calc_GT_etkf(int m, int p, double** Min, double** S, double** G, double** T);
void calc_X5(int m, int mout, double alpha, double* w, double** X5);

#define _CALCS_H
#endif
