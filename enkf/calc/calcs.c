/******************************************************************************
 *
 * File:        calcs.c        
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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include "definitions.h"
#include "lapack.h"
#include "utils.h"
#include "calcs.h"

#define EPS 1.0e-8

static char doT = 'T';
static char noT = 'N';

/** Calculates trace of the product of two matrices.
 * @param transposeA - flag: whether to transpose A
 * @param transposeB - flag: whether to transpose B
 * @param m          - number of rows/columns in A/B (after transposition)
 * @param n          - number of columns/rows in A/B (after transposition)
 * @param A          - n x m matrix (A[m][n]) (after transposition)
 * @param B          - m x n matrix (B[n][m]) (after transposition)
 * @param iscompact  - whether A and B are allocated as continuous blocks
 */
double traceprod(int transposeA, int transposeB, int m, int n, double** A, double** B, int iscompact)
{
    double trace = 0.0;
    int i, j;

    if (!iscompact) {
        if (!transposeA && !transposeB) {
            for (i = 0; i < n; ++i) {
                double* Bi = B[i];

                for (j = 0; j < m; ++j)
                    trace += A[j][i] * Bi[j];
            }
        } else if (transposeA && transposeB) {
            for (i = 0; i < n; ++i) {
                double* Ai = A[i];

                for (j = 0; j < m; ++j)
                    trace += Ai[j] * B[j][i];
            }
        } else if (transposeA || transposeB) {
            for (j = 0; j < m; ++j)
                for (i = 0; i < n; ++i)
                    trace += A[j][i] * B[j][i];
        }
    } else {
        if (!transposeA && !transposeB) {
            for (i = 0; i < n; ++i) {
                double* Aji = &A[0][i];
                double* Bij = &B[i][0];

                for (j = 0; j < m; ++j, Aji += n, Bij++)
                    trace += *Aji * *Bij;
            }
        } else if (transposeA && transposeB) {
            for (i = 0; i < n; ++i) {
                double* Bji = &B[0][i];
                double* Aij = &A[i][0];

                for (j = 0; j < m; ++j, Bji += n, Aij++)
                    trace += *Bji * *Aij;
            }
        } else if (transposeA || transposeB) {
            double* Aij = &A[0][0];
            double* Bij = &B[0][0];

            for (i = 0; i < n * m; ++i, Aij++, Bij++)
                trace += *Aij * *Bij;
        }
    }

    return trace;
}

/** Calculates inverse of a symmetric matrix via Cholesky decomposition.
 * @param m - matrix size
 * @param S - input: matrix; output: S^-1
 * @return lapack_info from dgesvd_()
 */
static int invm(int m, double** S)
{
    char uplo = 'U';
    int lapack_info;
    int i, j;

    dpotrf_(&uplo, &m, S[0], &m, &lapack_info);
    if (lapack_info != 0)
        return lapack_info;
    dpotri_(&uplo, &m, S[0], &m, &lapack_info);
    if (lapack_info != 0)
        return lapack_info;
#if 0
    for (j = 1; j < m; ++j)
        for (i = 0; i < j; ++i)
            S[i][j] = S[j][i];
#else
    for (j = 1; j < m; ++j) {
        double* colj_rowi = S[j];
        double* rowj_coli = &S[0][j];

        for (i = 0; i < j; ++i, colj_rowi++, rowj_coli += m)
            *rowj_coli = *colj_rowi;
    }
#endif
    return 0;
}

/** Calculates M = inv(I + S' * S) [transpose = 0] or M = inv(I + S * S')
 ** [transpose = 1].
 * @param p - size(S, 1)
 * @param m - size(S, 2)
 * @param transpose - flag: whether to transpose S
 * @param S - input, p x m
 * @param M - output, m x m [transpose = 0], p x p [transpose = 1]
 * @return lapack_info (0 = success)
 */
static int calc_M(int p, int m, int transpose, double** S, double** M)
{
    int i;
    double a = 1.0;
    double b = 0.0;
    int lapack_info;

    assert(p > 0);

    if (!transpose) {
        /*
         * M = S' * S 
         */
        dgemm_(&doT, &noT, &m, &m, &p, &a, S[0], &p, S[0], &p, &b, M[0], &m);

        /*
         * M = I + S' * S
         */
        for (i = 0; i < m; ++i)
            M[i][i] += 1.0;

        /*
         * M = inv(I + S' * S) 
         */
        lapack_info = invm(m, M);
    } else {
        /*
         * M = S * S' 
         */
        dgemm_(&noT, &doT, &p, &p, &m, &a, S[0], &p, S[0], &p, &b, M[0], &p);
        /*
         * M = I + S * S'
         */
        for (i = 0; i < p; ++i)
            M[i][i] += 1.0;

        lapack_info = invm(p, M);
    }

    return lapack_info;
}

/** Calculates inverse _and_ inverse square root of a square matrix via SVD.
 * @param m - matrix size
 * @param S - input: matrix; output: S^-1
 * @param D - output: S^-1/2
 * @return lapack_info from dgesvd_() (0 = success)
 * Note: the actual space allocated for S is [2 * m + 11][m]; the surplus is
 * used for work arrays and matrices.
 */
static int invsqrtm2(int m, double** S, double** D)
{
    double** U = &S[m];
    double* sigmas = S[2 * m];
    double* work = S[2 * m + 1];
    int lwork = 10 * m;
    char specU = 'A';           /* all M columns of U are returned in array U 
                                 */
    char specV = 'N';           /* no rows of V**T are computed */
    int lapack_info;
    double a = 1.0;
    double b = 0.0;
    int i, j;

    dgesvd_(&specU, &specV, &m, &m, S[0], &m, sigmas, U[0], &m, NULL, &m, work, &lwork, &lapack_info);
    if (lapack_info != 0)
        return lapack_info;

    for (i = 0; i < m; ++i) {
        double* Ui = U[i];
        double s = sqrt(sqrt(sigmas[i]));

        for (j = 0; j < m; ++j)
            Ui[j] /= s;
    }
    dgemm_(&noT, &doT, &m, &m, &m, &a, U[0], &m, U[0], &m, &b, D[0], &m);

    for (i = 0; i < m; ++i) {
        double* Ui = U[i];
        double s = sqrt(sqrt(sigmas[i]));

        for (j = 0; j < m; ++j)
            Ui[j] /= s;
    }
    dgemm_(&noT, &doT, &m, &m, &m, &a, U[0], &m, U[0], &m, &b, S[0], &m);

    return 0;
}

/** Calculates G = inv(I + S' * S) * S' = S' * inv(I + S * S').
 */
void calc_G(int m, int p, double** Min, double** S, int i, int j, double** G)
{
    double** M;
    double a = 1.0;
    double b = 0.0;
    int lapack_info;

    if (p < m) {
        M = (Min != NULL) ? cast2d(Min, p, p, sizeof(double)) : alloc2d(p, p, sizeof(double));
        /*
         * M = inv(I + S * S')
         */
        lapack_info = calc_M(p, m, 1, S, M);
        if (lapack_info != 0)
            enkf_quit("dpotrf() or dpotri(): lapack_info = %d at (i, j) = (%d, %d)", lapack_info, i, j);

        /*
         * G = S' * inv(I + S * S') 
         */
        dgemm_(&doT, &noT, &m, &p, &p, &a, S[0], &p, M[0], &p, &b, G[0], &m);
    } else {
        /*
         * M = inv(I + S' * S)
         */
        M = (Min != NULL) ? cast2d(Min, m, m, sizeof(double)) : alloc2d(m, m, sizeof(double));
        lapack_info = calc_M(p, m, 0, S, M);
        if (lapack_info != 0)
            enkf_quit("dpotrf() or dpotri(): lapack_info = %d at (i, j) = (%d, %d)", lapack_info, i, j);

        /*
         * G = inv(I + S * S') * S'
         */
        dgemm_(&noT, &doT, &m, &p, &m, &a, M[0], &m, S[0], &p, &b, G[0], &m);
    }

    if (Min == NULL)
        free(M);
}

/** Calculates X5 = G * s * 1' + T.
 * G is [m x p]
 * S is [p x m]
 * s is [p]
 * X5 is [m x m]
 */
void calc_T_denkf(int m, int p, double** G, double** S, double** T)
{
    double a, b;
    int i;

    /*
     * T <- -1/2 * G * S 
     */
    a = -0.5;
    b = 0.0;
    dgemm_(&noT, &noT, &m, &m, &p, &a, G[0], &m, S[0], &p, &b, T[0], &m);

    /*
     * T = I - 1/2 * G * S
     */
    for (i = 0; i < m; ++i)
        T[i][i] += 1.0;
}

/** Calculates G = inv(I + S' * S) * S' and T = (I + S' * S)^-1/2.
 */
 void calc_GT_etkf(int m, int p, double** Min, double** S, int ii, int jj, double** G, double** T)
{
    double** M;
    int lapack_info;
    int i;

    /*
     * dgemm stuff 
     */
    double a = 1.0;
    double b = 0.0;

    M = (Min != NULL) ? cast2d(Min, 2 * m + 11, m, sizeof(double)) : alloc2d(2 * m + 11, m, sizeof(double));
    /*
     * M = S' * S 
     */
    dgemm_(&doT, &noT, &m, &m, &p, &a, S[0], &p, S[0], &p, &b, M[0], &m);

    /*
     * M = I + S' * S
     */
    for (i = 0; i < m; ++i)
        M[i][i] += 1.0;

    lapack_info = invsqrtm2(m, M, T);   /* M = M^-1, T = M^-1/2 */
    if (lapack_info != 0)
        enkf_quit("dgesvd(): lapack_info = %d at (i, j) = (%d, %d)", lapack_info, ii, jj);

    /*
     * G = inv(I + S * S') * S'
     */
    dgemm_(&noT, &doT, &m, &p, &m, &a, M[0], &m, S[0], &p, &b, G[0], &m);

    if (Min == NULL)
        free(M);
}

/** Calculates X5 = w * 1' + I + alpha * (T - I).
 * @param m - size
 * @param alpha - relaxation parameter (0 <= alpha <= 1; 
 *        alpha = 1 -- no relaxation, alpha = 0 -- no anomalies update)
 * @param w[m] - weights
 * @param X5[m][m] - input: T, output: X5
 */
void calc_X5(int m, double alpha, double* w, double** X5)
{
    int i, j;

    /*
     * X5 = w * 1^T + I + alpha * (T - I)
     */
    for (i = 0; i < m; ++i) {
        double* X5i = X5[i];

        X5i[i] -= 1.0;
        for (j = 0; j < m; ++j)
            X5i[j] *= alpha;
        X5i[i] += 1.0;

        for (j = 0; j < m; ++j)
            X5i[j] += w[j];
    }
}

/** Calculates w = G * s.
 */
void calc_w(int m, int p, double** G, double* s, double* w)
{
    double a = 1.0;
    int inc = 1;
    double b = 0.0;

    /*
     * (no need to initialize w because b = 0, and w <- a * G * s + b * w)
     */
    dgemv_(&noT, &m, &p, &a, G[0], &m, s, &inc, &b, w, &inc);
}
