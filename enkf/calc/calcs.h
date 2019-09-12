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

double traceprod(int transposeA, int transposeB, int m, int n, double** A, double** B);
void calc_G(int m, int p, double** S, int i, int j, double** G);
void calc_T_denkf(int m, int p, double** G, double** S, double** T);
void calc_GT_etkf(int m, int p, double** S, int i, int j, double** G, double** T);
void calc_w(int m, int p, double** G, double* s, double* w);
void calc_X5(int m, double alpha, double* w, double** X5);

#define _CALCS_H
#endif
