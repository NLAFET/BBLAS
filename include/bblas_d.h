/**
 * @file bblas_d.h
 *
 * @brief BBLAS header file for double routines.
 *
 * BBLAS is a software package provided by Univ. of Manchester,
 * Univ. of Tennessee.
 *
 * @version 1.0.0
 * @author  Samuel  D. Relton
 * @author  Pedro   V. Lara
 * @author  Mawussi Zounon
 * @date    2016-02-20
 *
 **/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/**
 * Code generation
 * @generated from ../include/bblas_z.h normal z -> d, Mon Jun  6 09:44:13 2016
 **/
#endif

#ifndef BBLAS_D_H
#define BBLAS_D_H


#include "bblas_types.h"
#include "bblas_macros.h"
#include "auxiliary.h"


#define REAL

/*
 *  Declarations of level 3 BATCH BLAS  - alphabetical order
 */
void bblas_dgemm_batch(
    const enum BBLAS_TRANS *transA, const enum BBLAS_TRANS *transB,
    const int *M,  const int *N, const int *K,
    const double *alpha,
    const double **arrayA, const int *lda,
    const double **arrayB,
    const int *ldb, const double *beta,
    double **arrayC, const int *ldc, const int batch_count,
    const enum BBLAS_OPTS batch_opts, int *info);

#ifdef COMPLEX
void bblas_dsymm_batch(
    const enum BBLAS_SIDE *side, const enum BBLAS_UPLO *uplo,
    const int *M, const int *N, const double *alpha,
    const double **arrayA, const int *lda,
    const double **arrayB, const int *ldb,
    const double *beta, double **arrayC,
    const int *ldc, const int batch_count, const enum BBLAS_OPTS batch_opts, int *info);
#endif

void bblas_dsymm_batch(
    const enum BBLAS_SIDE *side, const enum BBLAS_UPLO *uplo,
    const int *M, const int *N, const double *alpha,
    const double **arrayA, const int *lda,
    const double **arrayB, const int *ldb,
    const double *beta, double **arrayC,
    const int *ldc,  const int batch_count, const enum BBLAS_OPTS batch_opts, int *info);

void bblas_dsyr2k_batch(
    const enum BBLAS_UPLO *uplo, const enum BBLAS_TRANS *trans,
    const int *N, const int *K, const double *alpha,
    const double **arrayA, const int *lda,
    const double **arrayB, const int *ldb,
    const double *beta, double **arrayC,
    const int *ldc, const int batch_count, const enum BBLAS_OPTS batch_opts, int *info);

#ifdef COMPLEX
void bblas_dsyr2k_batch(
    const enum BBLAS_UPLO *uplo, const enum  BBLAS_TRANS *trans,
    const int *N, const int *K, const double *alpha,
    const double **arrayA, const int *lda,
    const double **arrayB, const int *ldb,
    const double  *beta, double **arrayC,
    const int *ldc, const int batch_count, const enum BBLAS_OPTS batch_opts, int *info);
#endif

void bblas_dsyrk_batch(
    const enum BBLAS_UPLO *uplo, const enum BBLAS_TRANS *trans,
    const int *N, const int *K, const double *alpha,
    const double **arrayA, const int *lda,
    const double *beta, double **arrayC,
    const int *ldc, const int batch_count, const enum BBLAS_OPTS batch_opts, int *info);

#ifdef COMPLEX
void bblas_dsyrk_batch(
    const enum BBLAS_UPLO *uplo, const enum  BBLAS_TRANS *trans,
    const int *N, const int *K, const double *alpha,
    const double **arrayA, const int *lda,
    const double  *beta, double **arrayC,
    const int *ldc, const int batch_count, const  enum BBLAS_OPTS batch_opts, int *info);
#endif

void bblas_dtrmm_batch(
    const enum BBLAS_SIDE *side, const enum  BBLAS_UPLO *uplo,
    const enum BBLAS_TRANS *transA, const enum BBLAS_DIAG *diag,
    const int *M, const int *N, const double *alpha,
    const double **arrayA, const int *lda,
    double **arrayB, const int *ldb,
    const int batch_count, const  enum BBLAS_OPTS batch_opts, int *info);

void bblas_dtrsm_batch(
    const enum BBLAS_SIDE *side, const enum  BBLAS_UPLO *uplo,
    const enum BBLAS_TRANS *transA, const enum BBLAS_DIAG *diag,
    const int *M, const int *N, const double *alpha,
    const double **arrayA, const int *lda,
    double **arrayB, const int *ldb,
    const int batch_count, const  enum BBLAS_OPTS batch_opts, int *info);

int xerbla_batch(char *func_name, int error, int subproblem_ind);

#undef REAL
#endif /* BBLAS_D_H */
