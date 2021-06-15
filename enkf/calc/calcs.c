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
            for (j = 0; j < m; ++j) {
                double* Aj = A[j];
                double* Bj = B[j];

                for (i = 0; i < n; ++i)
                    trace += Aj[i] * Bj[i];
            }
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
 */
static void invm(int m, double** S)
{
    char uplo = 'U';
    int lapack_info;
    int i, j;

    dpotrf_(&uplo, &m, S[0], &m, &lapack_info);
    if (lapack_info != 0)
        /*
         * the failures of these procedured below are extremely rare, basically
         * unheard of, so we do not bother to elaborate
         */
        enkf_quit("dpotrf(): lapack_info = %d", lapack_info);
    dpotri_(&uplo, &m, S[0], &m, &lapack_info);
    if (lapack_info != 0)
        enkf_quit("dpotri(): lapack_info = %d", lapack_info);
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
}

/** Calculates Minv = inv(I + S' * S) if m <= p, or M = inv(I + S * S') if
 ** m > p.
 * @param p - size(S, 1)
 * @param m - size(S, 2)
 * @param S - input, p x m
 * @param M - output, m x m (m <= p), p x p (m > p)
 */
static void calc_Minv(int p, int m, double** S, double** M)
{
    int i;
    double a = 1.0;
    double b = 0.0;

    assert(p > 0);

    if (m <= p) {
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
        invm(m, M);
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

        invm(p, M);
    }
}

/** Calculates inverse _and_ inverse square root of a square matrix via SVD.
 * @param m - matrix size
 * @param S - input: matrix; output: S^-1
 * @param D - output: S^-1/2
 *
 * Note: the actual space allocated for S has to be [2 * m + 11][m]; the
 * surplus is used for work arrays and matrices.
 */
static void invsqrtm2(int m, double** S, double** D)
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
        /*
         * the failures are extremely rare so we do not bother to elaborate
         */
        enkf_quit("dgesvd(): lapack_info = %d", lapack_info);

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
}

/** Calculates G = inv(I + S' * S) * S' = S' * inv(I + S * S').
 */
void calc_G(int m, int p, double** Min, double** S, double** G)
{
    double** M;
    double a = 1.0;
    double b = 0.0;

    if (p < m) {
        M = (Min != NULL) ? cast2d(Min, p, p, sizeof(double)) : alloc2d(p, p, sizeof(double));
        /*
         * M = inv(I + S * S')
         */
        calc_Minv(p, m, S, M);

        /*
         * G = S' * inv(I + S * S') 
         */
        dgemm_(&doT, &noT, &m, &p, &p, &a, S[0], &p, M[0], &p, &b, G[0], &m);
    } else {
        /*
         * M = inv(I + S' * S)
         */
        M = (Min != NULL) ? cast2d(Min, m, m, sizeof(double)) : alloc2d(m, m, sizeof(double));
        calc_Minv(p, m, S, M);

        /*
         * G = inv(I + S * S') * S'
         */
        dgemm_(&noT, &doT, &m, &p, &m, &a, M[0], &m, S[0], &p, &b, G[0], &m);
    }

    if (Min == NULL)
        free(M);
}

/** Calculates w = G * s.
 */
void calc_w(int m, int p, double** G, double* s, double* w)
{
    double a = 1.0;
    int inc = 1;
    double b = 0.0;

    dgemv_(&noT, &m, &p, &a, G[0], &m, s, &inc, &b, w, &inc);
}

/** Calculate ensemble transform for the DEnKF.
 * @param m - full ensemble size
 * @param mout - actual (dynamic) ensemble size
 * @param p - (local) number of obs
 * @param s - s (see sec. 1.3.3 of the User Guide), calculated over the dynamic
 *            ensemble
 * @param S - S full S, p x m (see sec. 1.3.3 of the User Guide)
 * @param Sa - actual S, p x mout, can be same as S, or the leading part of S
 *             (for hybrid systems), or different from S (for modulated systems)
 * @param M - storage for intermediate matrix M (can be NULL)
 * @param G - storage for intermediate matrix G
 * @return w - mean update coefficients
 * @return T - ensemble anomalies transform matrix
 */
void calc_wT_denkf(int m, int mout, int p, double* s, double** S, double** Sa, double** M, double** G, double* w, double** T)
{
    calc_G(m, p, M, S, G);
    calc_T_denkf(m, mout, p, G, Sa, T);
    calc_w(m, p, G, s, w);
}

/** Calculate ensemble transform for the ETKF.
 * @param m - full ensemble size
 * @param mout - actual (dynamic) ensemble size
 * @param p - (local) number of obs
 * @param s - s (see sec. 1.3.3 of the User Guide), calculated over the dynamic
 *            ensemble
 * @param S - S full S, p x m (see sec. 1.3.3 of the User Guide)
 * @param Sa - actual S, p x mout, can be same as S, or the leading part of S
 *             (for hybrid systems), or different from S (for modulated systems)
 * @param M - storage for intermediate matrix M (can be NULL)
 * @param G - storage for intermediate matrix G
 * @return w - mean update coefficients
 * @return T - ensemble anomalies transform matrix
 */
void calc_wT_etkf(int m, int mout, int p, double* s, double** S, double** Sa, double** Min, double** G, double* w, double** T)
{
    if (p < m) {
        double** M = (Min != NULL) ? cast2d(Min, 2 * p + 11, p, sizeof(double)) : alloc2d(2 * p + 11, p, sizeof(double));
        double** U = &M[p];

        double* sigmas = M[2 * p];
        double* work = M[2 * p + 1];

        int lwork = 10 * p;
        char specU = 'A';       /* all M columns of U are returned in array U 
                                 */
        char specV = 'N';       /* no rows of V**T are computed */
        double a = 1.0;
        double b = 0.0;
        int lapack_info;
        int i;

        M = cast2d(M, 2 * p + 11, p, sizeof(double));
        U = &M[p];
        work = M[2 * p + 1];

        /*
         * M = S * S' 
         */
        dgemm_(&noT, &doT, &p, &p, &m, &a, S[0], &p, S[0], &p, &b, M[0], &p);
        /*
         * M = I + S * S'
         */
        for (i = 0; i < p; ++i)
            M[i][i] += 1.0;

        /*
         * svd(M): M = U * diag(sigmas) * U'
         */
        dgesvd_(&specU, &specV, &p, &p, M[0], &p, sigmas, U[0], &p, NULL, &p, work, &lwork, &lapack_info);
        if (lapack_info != 0)
            /*
             * the failures are extremely rare so we do not bother to elaborate
             */
            enkf_quit("dgesvd(): lapack_info = %d", lapack_info);

        /*
         * calculate Mp = (M + M^1/2)^-1 and store in M
         */
        for (i = 0; i < p; ++i) {
            double* Ui = U[i];
            double s = sqrt(sigmas[i] + sqrt(sigmas[i]));
            int j;

            for (j = 0; j < p; ++j)
                Ui[j] /= s;
        }
        dgemm_(&noT, &doT, &p, &p, &p, &a, U[0], &p, U[0], &p, &b, M[0], &p);
        /*
         * calculate S' * (Mp + Mp^1/2)^-1 and store in G
         */
        dgemm_(&doT, &noT, &m, &p, &p, &a, S[0], &p, M[0], &p, &b, G[0], &m);
        /*
         * calculate S' * (Mp + Mp^1/2)^-1 * Sa and store in T
         */
        dgemm_(&noT, &noT, &m, &mout, &p, &a, G[0], &m, Sa[0], &p, &b, T[0], &m);
        /*
         * calculate T = I - S' * (Mp + Mp^1/2)^-1 * Sa
         */
        for (i = 0; i < mout; ++i) {
            double* Ti = T[i];
            int j;

            for (j = 0; j < m; ++j)
                Ti[j] = -Ti[j];
            Ti[i] += 1.0;
        }

        /*
         * calculate M^-1 and store in M
         */
        for (i = 0; i < p; ++i) {
            double* Ui = U[i];
            double s = sqrt(1.0 + 1.0 / sqrt(sigmas[i]));
            int j;

            for (j = 0; j < p; ++j)
                Ui[j] *= s;
        }
        dgemm_(&noT, &doT, &p, &p, &p, &a, U[0], &p, U[0], &p, &b, M[0], &p);
        /*
         * calculate G = S' * (I + S * S')^-1
         */
        dgemm_(&doT, &noT, &m, &p, &p, &a, S[0], &p, M[0], &p, &b, G[0], &m);

        if (Min == NULL)
            free(M);
    } else {
        double** M = (Min != NULL) ? cast2d(Min, 2 * m + 11, m, sizeof(double)) : alloc2d(2 * m + 11, m, sizeof(double));
        double** U = &M[m];
        double* sigmas = M[2 * m];
        double* work = M[2 * m + 1];

        int lwork = 10 * m;
        char specU = 'A';       /* all M columns of U are returned in array U 
                                 * no rows of V**T are computed */
        char specV = 'N';       /* no rows of V**T are computed */
        double a = 1.0;
        double b = 0.0;
        int lapack_info;
        int i;

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
         * svd(M): M = U * diag(sigmas) * U'
         */
        dgesvd_(&specU, &specV, &m, &m, M[0], &m, sigmas, U[0], &m, NULL, &m, work, &lwork, &lapack_info);
        if (lapack_info != 0)
            /*
             * the failures are extremely rare so we do not bother to elaborate
             */
            enkf_quit("dgesvd(): lapack_info = %d", lapack_info);

        /*
         * calculate Mm = (M + M^1/2)^-1 and store in M
         */
        for (i = 0; i < m; ++i) {
            double* Ui = U[i];
            double s = sqrt(sigmas[i] + sqrt(sigmas[i]));
            int j;

            for (j = 0; j < m; ++j)
                Ui[j] /= s;
        }
        dgemm_(&noT, &doT, &m, &m, &m, &a, U[0], &m, U[0], &m, &b, M[0], &m);
        /*
         * calculate (Mm + Mm^1/2)^-1 * S' and store in G
         */
        dgemm_(&noT, &doT, &m, &p, &m, &a, M[0], &m, S[0], &p, &b, G[0], &m);
        /*
         * calculate (Mm + Mm^1/2)^-1 * S' * Sa and store in T
         */
        dgemm_(&noT, &noT, &m, &mout, &p, &a, G[0], &m, Sa[0], &p, &b, T[0], &m);
        /*
         * calculate T = I - (Mm + Mm^1/2)^-1 * S' * Sa
         */
        for (i = 0; i < mout; ++i) {
            double* Ti = T[i];
            int j;

            for (j = 0; j < m; ++j)
                Ti[j] = -Ti[j];
            Ti[i] += 1.0;
        }

        /*
         * calculate M^-1 and store in M
         */
        for (i = 0; i < m; ++i) {
            double* Ui = U[i];
            double s = sqrt(1.0 + 1.0 / sqrt(sigmas[i]));
            int j;

            for (j = 0; j < m; ++j)
                Ui[j] *= s;
        }
        dgemm_(&noT, &doT, &m, &m, &m, &a, U[0], &m, U[0], &m, &b, M[0], &m);
        /*
         * calculate G = (I + S' * S)^-1 * S'
         */
        dgemm_(&noT, &doT, &m, &p, &m, &a, M[0], &m, S[0], &p, &b, G[0], &m);

        if (Min == NULL)
            free(M);
    }

    calc_w(m, p, G, s, w);
}

/*
 * Below are now redundant procedures (on the way out).
 */

/** Calculates G = inv(I + S' * S) * S' and T = (I + S' * S)^-1/2.
 */
void calc_GT_etkf(int m, int p, double** Min, double** S, double** G, double** T)
{
    double** M;
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

    invsqrtm2(m, M, T);         /* M = M^-1, T = M^-1/2 */

    /*
     * G = inv(I + S * S') * S'
     */
    dgemm_(&noT, &doT, &m, &p, &m, &a, M[0], &m, S[0], &p, &b, G[0], &m);

    if (Min == NULL)
        free(M);
}

/** Calculates T = I - 1/2 * G * S.
 */
void calc_T_denkf(int m, int mout, int p, double** G, double** S, double** T)
{
    double a, b;
    int i;

    /*
     * T <- -1/2 * G * S 
     */
    a = -0.5;
    b = 0.0;
    dgemm_(&noT, &noT, &m, &mout, &p, &a, G[0], &m, S[0], &p, &b, T[0], &m);

    /*
     * T = I - 1/2 * G * S
     */
    for (i = 0; i < mout; ++i)
        T[i][i] += 1.0;
}

/** Calculates X5 = w * 1' + I + alpha * (T - I).
 * @param m - full ensemble size
 * @param mout - actual ensemble size (m_dynamic for hybrid mode)
 * @param alpha - relaxation parameter (0 <= alpha <= 1; 
 *        alpha = 1 -- no relaxation, alpha = 0 -- no anomalies update)
 * @param w[m] - weights
 * @param X5[m][m] - input: T, output: X5
 */
void calc_X5(int m, int mout, double alpha, double* w, double** X5)
{
    int i, j;

    /*
     * X5 = w * 1^T + I + alpha * (T - I)
     */
    for (i = 0; i < mout; ++i) {
        double* X5i = X5[i];

        X5i[i] -= 1.0;
        for (j = 0; j < m; ++j)
            X5i[j] *= alpha;
        X5i[i] += 1.0;

        for (j = 0; j < m; ++j)
            X5i[j] += w[j];
    }
}
