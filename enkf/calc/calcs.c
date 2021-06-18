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
 * @return           - trace(op(A) * op(B))
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
 * @param M - input: m x m matrix; output: M^-1
 */
static void invm(int m, double** M)
{
    char uplo = 'U';
    int lapack_info;
    int i, j;

    dpotrf_(&uplo, &m, M[0], &m, &lapack_info);
    if (lapack_info != 0)
        /*
         * the failures of these procedured below are extremely rare, basically
         * unheard of, so we do not bother to elaborate
         */
        enkf_quit("dpotrf(): lapack_info = %d", lapack_info);
    dpotri_(&uplo, &m, M[0], &m, &lapack_info);
    if (lapack_info != 0)
        enkf_quit("dpotri(): lapack_info = %d", lapack_info);
#if 0
    for (j = 1; j < m; ++j)
        for (i = 0; i < j; ++i)
            M[i][j] = M[j][i];
#else
    for (j = 1; j < m; ++j) {
        double* colj_rowi = M[j];
        double* rowj_coli = &M[0][j];

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

/** Calculates G = inv(I + S' * S) * S' = S' * inv(I + S * S'), eqs. (1.37)
 ** and (1.38) of the User Guide.
 * @param m - ensemble size
 * @param p - number of obs.
 * @param Min - storage space for M of size >= n x n x sizeof(double),
 *              n = min(m, p); or NULL
 * @param S - p x m matrix of normalised ensemble observation anomalies
 * @param G - output, m x p
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

/** Calculates ETM for the DEnKF, T = I - 1/2 * G * S, eq. (1.36) of the User
 ** Guide.
 * @param m - full (expanded) ensemble size
 * @param mout - actual (dynamic) ensemble size
 * @param p - number of observations
 * @param G - m x p matrix, G = inv(I + S' * S) * S' = S' * inv(I + S * S')
 * @param S - p x m matrix of normalised ensemble observation anomalies
 * @param T - output ETM, m x m
 */
static void calc_T_denkf(int m, int mout, int p, double** G, double** S, double** T)
{
    double a = -0.5;
    double b = 0.0;
    int i;

    /*
     * T <- -1/2 * G * S 
     */
    dgemm_(&noT, &noT, &m, &mout, &p, &a, G[0], &m, S[0], &p, &b, T[0], &m);

    /*
     * T = I - 1/2 * G * S
     */
    for (i = 0; i < mout; ++i)
        T[i][i] += 1.0;
}

/** Calculates increment coefficients w = G * s, eq. (1.34) of the User Guide.
 * @param m - ensemble size = size(G, 1) 
 * @param p - number of obs. = size(G, 2)
 * @param G - m x p matrix, G = inv(I + S' * S) * S' = S' * inv(I + S * S').
 * @param s - normalised innovations
 */
void calc_w(int m, int p, double** G, double* s, double* w)
{
    double a = 1.0;
    double b = 0.0;
    int inc = 1;

    dgemv_(&noT, &m, &p, &a, G[0], &m, s, &inc, &b, w, &inc);
}

/** Calculate ensemble transform for the DEnKF, eq. (1.36) of the User Guide.
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
 * @param w - (output) mean update coefficients
 * @param T - (output) ensemble anomalies transform matrix
 */
void calc_wT_denkf(int m, int mout, int p, double* s, double** S, double** Sa, double** M, double** G, double* w, double** T)
{
    calc_G(m, p, M, S, G);
    calc_T_denkf(m, mout, p, G, Sa, T);
    calc_w(m, p, G, s, w);
}

/** Calculate ensemble transform for the ETKF, eqs. (1.29) of the User Guide.
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
 * @param w - (output) mean update coefficients
 * @param T - (output) ensemble anomalies transform matrix
 */
void calc_wT_etkf(int m, int mout, int p, double* s, double** S, double** Sa, double** Min, double** G, double* w, double** T)
{
    int n = (p < m) ? p : m;

    double** M = (Min != NULL) ? cast2d(Min, 2 * n + 11, n, sizeof(double)) : alloc2d(2 * n + 11, n, sizeof(double));
    double** U = &M[n];
    double* sigmas = M[2 * n];
    double a = 1.0;
    double b = 0.0;
    int i;

    if (p < m)
        /*
         * calculate S * S' and store in M
         */
        dgemm_(&noT, &doT, &p, &p, &m, &a, S[0], &p, S[0], &p, &b, M[0], &p);
    else
        /*
         * calculate S' * S and store in M
         */
        dgemm_(&doT, &noT, &m, &m, &p, &a, S[0], &p, S[0], &p, &b, M[0], &m);
    /*
     * Mp = I + S * S' or Mm = I + S' * S
     */
    for (i = 0; i < n; ++i)
        M[i][i] += 1.0;

    /*
     * svd(M): Mp or Mm = U * diag(sigmas) * U'
     */
    {
        double* work = M[2 * n + 1];
        int lwork = 10 * n;
        char specU = 'A';       /* all M columns of U are returned in array U 
                                 */
        char specV = 'N';       /* no rows of V**T are computed */
        int lapack_info;

        dgesvd_(&specU, &specV, &n, &n, M[0], &n, sigmas, U[0], &n, NULL, &n, work, &lwork, &lapack_info);
        if (lapack_info != 0)
            /*
             * the failures are extremely rare so we do not bother to elaborate
             */
            enkf_quit("dgesvd(): lapack_info = %d", lapack_info);
    }

    /*
     * calculate (Mp + Mp^1/2)^-1 or (Mm + Mm^1/2)^-1 and store in M
     */
    for (i = 0; i < n; ++i) {
        double* Ui = U[i];
        double s = sqrt(sigmas[i] + sqrt(sigmas[i]));
        int j;

        for (j = 0; j < n; ++j)
            Ui[j] /= s;
    }
    dgemm_(&noT, &doT, &n, &n, &n, &a, U[0], &n, U[0], &n, &b, M[0], &n);

    if (p < m)
        /*
         * calculate S' * (Mp + Mp^1/2)^-1 and store in G
         */
        dgemm_(&doT, &noT, &m, &p, &p, &a, S[0], &p, M[0], &p, &b, G[0], &m);
    else
        /*
         * calculate (Mm + Mm^1/2)^-1 * S' and store in G
         */
        dgemm_(&noT, &doT, &m, &p, &m, &a, M[0], &m, S[0], &p, &b, G[0], &m);

    /*
     * calculate -S' * (Mp + Mp^1/2)^-1 * Sa or -(Mm + Mm^1/2)^-1 * S' * Sa
     * and store in T
     */
    a = -1.0;
    dgemm_(&noT, &noT, &m, &mout, &p, &a, G[0], &m, Sa[0], &p, &b, T[0], &m);
    a = 1.0;
    /*
     * calculate T = I - (Mm + Mm^1/2)^-1 * S' * Sa
     *             = I - S' * (Mp + Mp^1/2)^-1 * Sa
     */
    for (i = 0; i < mout; ++i)
        T[i][i] += 1.0;

    /*
     * calculate M^-1 and store in M
     */
    for (i = 0; i < n; ++i) {
        double* Ui = U[i];
        double s = sqrt(1.0 + 1.0 / sqrt(sigmas[i]));
        int j;

        for (j = 0; j < n; ++j)
            Ui[j] *= s;
    }
    dgemm_(&noT, &doT, &n, &n, &n, &a, U[0], &n, U[0], &n, &b, M[0], &n);

    if (p < m)
        /*
         * calculate G = S' * (I + S * S')^-1
         */
        dgemm_(&doT, &noT, &m, &p, &p, &a, S[0], &p, M[0], &p, &b, G[0], &m);
    else
        /*
         * calculate G = (I + S' * S)^-1 * S'
         */
        dgemm_(&noT, &doT, &m, &p, &m, &a, M[0], &m, S[0], &p, &b, G[0], &m);

    calc_w(m, p, G, s, w);

    if (Min == NULL)
        free(M);
}
