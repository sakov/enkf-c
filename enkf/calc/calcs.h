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

void shuffle(int n, int ids[]);
double traceprod(int transposeA, int transposeB, int m, int n, double** A, double** B);
void calc_G_denkf(int m, int p, double** S, int i, int j, double** G);
void calc_X5_denkf(int m, int p, double** G, double** S, double* s, double alpha, int ii, int jj, double** X5);
void calc_G_etkf(int m, int p, double** S, double alpha, int i, int j, double** G, double** T);
void calc_X5_etkf(int m, int p, double** G, double* s, int ii, int jj, double** X5);
void calc_w(int m, int p, double** G, double* s, double* w);

#define _CALCS_H
#endif
