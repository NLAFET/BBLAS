/**
 *
 * @file
 *
 *  PLASMA header.
 *  PLASMA is a software package provided by Univ. of Tennessee,
 *  Univ. of Manchester, Univ. of California Berkeley and
 *  Univ. of Colorado Denver.
 *
 * @precisions normal z -> s d c
 *
 **/

#ifndef BBLAS_Z_H
#define BBLAS_Z_H


#include "bblas_types.h"
#include "bblas_macros.h"
#include "auxiliary.h"


#define COMPLEX

void blas_zgemm_batch(int group_count, const int *group_sizes,
		bblas_enum_t layout, const bblas_enum_t *transa, const bblas_enum_t *transb,
		const int *m, const int *n, const int *k,
		const bblas_complex64_t *alpha, bblas_complex64_t const *const *A, const int *lda,
		bblas_complex64_t const* const *B, const int *ldb, 
		const bblas_complex64_t *beta,  bblas_complex64_t** C, const int *ldc, 
		int *info);


#undef COMPLEX
#endif /* BBLAS_Z_H */
