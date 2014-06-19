/******************************************************************************
 *
 * File:        lapack.h        
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

#if !defined(_LAPACK_H)

/* C = alpha * A * B + beta * C */
int dgemm_(char* transa, char* transb, int* m, int* n, int* k, double* alpha, double* A, int* lda, double* B, int* ldb, double* beta, double* C, int* ldc);

/* y = alpha * A * x + beta * y */
int dgemv_(char* trans, int* m, int* n, double* alpha, double* A, int* lda, double* x, int* incx, double* beta, double* y, int* incy);
int sgemv_(char* trans, int* m, int* n, float* alpha, float* A, int* lda, float* x, int* incx, float* beta, float* y, int* incy);

int dpotrf_(char* uplo, int* n, double* A, int* lda, int* info);
int dpotri_(char* uplo, int* n, double* A, int* lda, int* info);

int dgesvd_(char* jobu, char* jobvt, int* m, int* n, double* A, int* lda, double* S, double* U, int* ldu, double* VT, int* ldvt, double* work, int* lwork, int* Info);

#define _LAPACK_H
#endif
