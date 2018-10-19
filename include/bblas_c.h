/**
 * @file bblas_c.h
 *
 * @brief BBLAS header file for float _Complex routines.
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
 * @generated from ../include/bblas_z.h normal z -> c, Mon Jun  6 09:44:13 2016
 **/
#endif

#ifndef BBLAS_C_H
#define BBLAS_C_H


#include "bblas_types.h"
#include "bblas_macros.h"
#include "auxiliary.h"


#define COMPLEX

/*
 *  Declarations of level 3 BATCH BLAS  - alphabetical order
 */
void bblas_cgemm_batch(
    const enum BBLAS_TRANS *transA, const enum BBLAS_TRANS *transB,
    const int *M,  const int *N, const int *K,
    const BBLAS_Complex32_t *alpha,
    const BBLAS_Complex32_t **arrayA, const int *lda,
    const BBLAS_Complex32_t **arrayB,
    const int *ldb, const BBLAS_Complex32_t *beta,
    BBLAS_Complex32_t **arrayC, const int *ldc, const int batch_count,
    const enum BBLAS_OPTS batch_opts, int *info);

#ifdef COMPLEX
void bblas_chemm_batch(
    const enum BBLAS_SIDE *side, const enum BBLAS_UPLO *uplo,
    const int *M, const int *N, const BBLAS_Complex32_t *alpha,
    const BBLAS_Complex32_t **arrayA, const int *lda,
    const BBLAS_Complex32_t **arrayB, const int *ldb,
    const BBLAS_Complex32_t *beta, BBLAS_Complex32_t **arrayC,
    const int *ldc, const int batch_count, const enum BBLAS_OPTS batch_opts, int *info);
#endif

void bblas_csymm_batch(
    const enum BBLAS_SIDE *side, const enum BBLAS_UPLO *uplo,
    const int *M, const int *N, const BBLAS_Complex32_t *alpha,
    const BBLAS_Complex32_t **arrayA, const int *lda,
    const BBLAS_Complex32_t **arrayB, const int *ldb,
    const BBLAS_Complex32_t *beta, BBLAS_Complex32_t **arrayC,
    const int *ldc,  const int batch_count, const enum BBLAS_OPTS batch_opts, int *info);

void bblas_csyr2k_batch(
    const enum BBLAS_UPLO *uplo, const enum BBLAS_TRANS *trans,
    const int *N, const int *K, const BBLAS_Complex32_t *alpha,
    const BBLAS_Complex32_t **arrayA, const int *lda,
    const BBLAS_Complex32_t **arrayB, const int *ldb,
    const BBLAS_Complex32_t *beta, BBLAS_Complex32_t **arrayC,
    const int *ldc, const int batch_count, const enum BBLAS_OPTS batch_opts, int *info);

#ifdef COMPLEX
void bblas_cher2k_batch(
    const enum BBLAS_UPLO *uplo, const enum  BBLAS_TRANS *trans,
    const int *N, const int *K, const BBLAS_Complex32_t *alpha,
    const BBLAS_Complex32_t **arrayA, const int *lda,
    const BBLAS_Complex32_t **arrayB, const int *ldb,
    const float  *beta, BBLAS_Complex32_t **arrayC,
    const int *ldc, const int batch_count, const enum BBLAS_OPTS batch_opts, int *info);
#endif

void bblas_csyrk_batch(
    const enum BBLAS_UPLO *uplo, const enum BBLAS_TRANS *trans,
    const int *N, const int *K, const BBLAS_Complex32_t *alpha,
    const BBLAS_Complex32_t **arrayA, const int *lda,
    const BBLAS_Complex32_t *beta, BBLAS_Complex32_t **arrayC,
    const int *ldc, const int batch_count, const enum BBLAS_OPTS batch_opts, int *info);

#ifdef COMPLEX
void bblas_cherk_batch(
    const enum BBLAS_UPLO *uplo, const enum  BBLAS_TRANS *trans,
    const int *N, const int *K, const float *alpha,
    const BBLAS_Complex32_t **arrayA, const int *lda,
    const float  *beta, BBLAS_Complex32_t **arrayC,
    const int *ldc, const int batch_count, const  enum BBLAS_OPTS batch_opts, int *info);
#endif

void bblas_ctrmm_batch(
    const enum BBLAS_SIDE *side, const enum  BBLAS_UPLO *uplo,
    const enum BBLAS_TRANS *transA, const enum BBLAS_DIAG *diag,
    const int *M, const int *N, const BBLAS_Complex32_t *alpha,
    const BBLAS_Complex32_t **arrayA, const int *lda,
    BBLAS_Complex32_t **arrayB, const int *ldb,
    const int batch_count, const  enum BBLAS_OPTS batch_opts, int *info);

void bblas_ctrsm_batch(
    const enum BBLAS_SIDE *side, const enum  BBLAS_UPLO *uplo,
    const enum BBLAS_TRANS *transA, const enum BBLAS_DIAG *diag,
    const int *M, const int *N, const BBLAS_Complex32_t *alpha,
    const BBLAS_Complex32_t **arrayA, const int *lda,
    BBLAS_Complex32_t **arrayB, const int *ldb,
    const int batch_count, const  enum BBLAS_OPTS batch_opts, int *info);

int xerbla_batch(char *func_name, int error, int subproblem_ind);

#undef COMPLEX
#endif /* BBLAS_C_H */
