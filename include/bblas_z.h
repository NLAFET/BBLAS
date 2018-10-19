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

#include "core_z.h"

#define COMPLEX

void blas_zgemm_batch(int group_count, const int *group_sizes,
		bblas_enum_t layout, const bblas_enum_t *transa, const bblas_enum_t *transb,
		const int *m, const int *n, const int *k,
		const bblas_complex64_t *alpha, bblas_complex64_t const *const *A, const int *lda,
		bblas_complex64_t const* const *B, const int *ldb, 
		const bblas_complex64_t *beta,  bblas_complex64_t** C, const int *ldc, 
		int *info);


void blas_zhemm_batch( int group_count, const int *group_sizes,
			bblas_enum_t layout, const bblas_enum_t *side, const bblas_enum_t *uplo,
    			const int *m, const int *n, 
			const bblas_complex64_t *alpha, bblas_complex64_t const *const *A, const int *lda, 
    			 		 	 bblas_complex64_t const* const *B, const int *ldb, 
			const bblas_complex64_t *beta,  bblas_complex64_t** C, const int *ldc, 
    			int *info);


void blas_zher2k_batch( int group_count, const int *group_sizes,
			bblas_enum_t layout, const bblas_enum_t *uplo, const bblas_enum_t *trans,
    			const int *n, const int *k, 
			const bblas_complex64_t *alpha, bblas_complex64_t const *const *A, const int *lda, 
    			 		 	 bblas_complex64_t const* const *B, const int *ldb, 
			const double  *beta,  bblas_complex64_t** C, const int *ldc, 
    			int *info);



void blas_zherk_batch( int group_count, const int *group_sizes,
		      bblas_enum_t layout, const bblas_enum_t *uplo, const bblas_enum_t *trans,
		      const int *n, const int *k, 
		      const double *alpha, bblas_complex64_t const *const *A, const int *lda, 
		      const double  *beta,  bblas_complex64_t** C, const int *ldc, 
    		      int *info);



void blas_zsymm_batch( int group_count, const int *group_sizes,
		        bblas_enum_t layout, const bblas_enum_t *side, const bblas_enum_t *uplo,
    			const int *m, const int *n, 
			const bblas_complex64_t *alpha, bblas_complex64_t const *const *A, const int *lda, 
    			 		 	        bblas_complex64_t const* const *B, const int *ldb, 
			const bblas_complex64_t *beta,  bblas_complex64_t** C, const int *ldc, 
    			 int *info);



void blas_zsyr2k_batch( int group_count, const int *group_sizes, 
			bblas_enum_t layout, const bblas_enum_t *uplo, const bblas_enum_t *trans,
    			const int *n, const int *k, 
			const bblas_complex64_t *alpha, bblas_complex64_t const *const *A, const int *lda, 
    			 		 	        bblas_complex64_t const* const *B, const int *ldb, 
			const double  *beta,  bblas_complex64_t** C, const int *ldc, 
    			 int *info);



void blas_zsyrk_batch( int group_count, const int *group_sizes,
		      bblas_enum_t layout, const bblas_enum_t *uplo, const bblas_enum_t *trans,
		      const int *n, const int *k, 
		      const double *alpha, bblas_complex64_t const *const *A, const int *lda, 
		      const double  *beta,  bblas_complex64_t** C, const int *ldc, 
    		      int *info);



void blas_ztrmm_batch( int group_count, const int *group_sizes,
			bblas_enum_t layout, const bblas_enum_t *side, const bblas_enum_t *uplo,
    			const bblas_enum_t *transa, const bblas_enum_t *diag,
    			const int *m, const int *n, 
			const bblas_complex64_t *alpha, bblas_complex64_t const *const *A, const int *lda,
    			bblas_complex64_t **B, int const *ldb,
			 int *info);



void blas_ztrsm_batch(int group_count, const int *group_sizes,
			bblas_enum_t layout, const bblas_enum_t *side, const bblas_enum_t *uplo,
    			const bblas_enum_t *transa, const bblas_enum_t *diag,
    			const int *m, const int *n, 
			const bblas_complex64_t *alpha, bblas_complex64_t const *const *A, const int *lda,
    			bblas_complex64_t **B, const int *ldb,
			int *info);



#undef COMPLEX
#endif /* BBLAS_Z_H */
