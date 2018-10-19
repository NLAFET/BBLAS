/**
 *
 * @file
 *
 *  BBLAS is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @precisions normal z -> s d c
 *
 **/
#ifndef ICL_CORE_Z_H
#define ICL_CORE_Z_H

#include "bblas_types.h"
#include "bblas_macros.h"
#include "auxiliary.h"

#ifdef __cplusplus
extern "C" {
#endif

#define COMPLEX


void blas_zgemm_batchf (int group_size, 
			bblas_enum_t layout, bblas_enum_t transa, bblas_enum_t transb,
			int m,  int n, int k,
			bblas_complex64_t alpha, bblas_complex64_t const *const *A, int lda,
    			 		 	 bblas_complex64_t const* const *B, int ldb, 
			bblas_complex64_t beta,  bblas_complex64_t** C, int ldc, 
			int *info);

#ifdef COMPLEX
void blas_zhemm_batchf(int group_size,
			bblas_enum_t layout, bblas_enum_t side, bblas_enum_t uplo,
    			int m, int n, 
			bblas_complex64_t alpha, bblas_complex64_t const *const *A, int lda, 
    			 		 	 bblas_complex64_t const* const *B, int ldb, 
			bblas_complex64_t beta,  bblas_complex64_t** C, int ldc, 
    			int *info);
#endif

void blas_zsymm_batchf( int group_size,
		        bblas_enum_t layout, bblas_enum_t side, bblas_enum_t uplo,
    			int m,  int n, 
			bblas_complex64_t alpha, bblas_complex64_t const *const *A, int lda, 
    			 		 	 bblas_complex64_t const* const *B, int ldb, 
			bblas_complex64_t beta,  bblas_complex64_t** C, int ldc, 
    			 int *info);

void blas_zsyr2k_batchf( int group_size, 
			bblas_enum_t layout, bblas_enum_t uplo, bblas_enum_t trans,
    			int n, int k, 
			bblas_complex64_t alpha, bblas_complex64_t const *const *A, int lda, 
    			 		 	 bblas_complex64_t const* const *B, int ldb, 
			const double  beta,  bblas_complex64_t** C, int ldc, 
    			 int *info);

void blas_zsyrk_batchf( int group_size,
		      bblas_enum_t layout, bblas_enum_t uplo, bblas_enum_t trans,
		      int n, int k, 
		      const double alpha, bblas_complex64_t const *const *A, int lda, 
		      const double  beta,  bblas_complex64_t** C, int ldc, 
    		      int *info);


#ifdef COMPLEX
void blas_zher2k_batchf(int group_size,
			bblas_enum_t layout, bblas_enum_t uplo, bblas_enum_t trans,
    			int n, int k, 
			bblas_complex64_t alpha, bblas_complex64_t const *const *A, int lda, 
    			 		 	 bblas_complex64_t const* const *B, int ldb, 
			const double  beta,  bblas_complex64_t** C, int ldc, 
    			int *info);
#endif

#ifdef COMPLEX
void blas_zherk_batchf( int group_size,
		      bblas_enum_t layout, bblas_enum_t uplo, bblas_enum_t trans,
		      int n, int k, 
		      const double alpha, bblas_complex64_t const *const *A, int lda, 
		      const double  beta,  bblas_complex64_t** C, int ldc, 
    		      int *info);
#endif

void blas_ztrmm_batchf( int group_size,
			bblas_enum_t layout, bblas_enum_t side, bblas_enum_t uplo,
    			bblas_enum_t transa, bblas_enum_t diag,
    			int m, int n, 
			bblas_complex64_t alpha, bblas_complex64_t const *const *A, int lda,
    			bblas_complex64_t **B, int ldb,
			 int *info);

void blas_ztrsm_batchf( int group_size,
			bblas_enum_t layout, bblas_enum_t side, bblas_enum_t uplo,
    			bblas_enum_t transa, bblas_enum_t diag,
    			int m, int n, 
			bblas_complex64_t alpha, bblas_complex64_t const *const *A, int lda,
    			bblas_complex64_t **B, int ldb,
			int *info);

int xerbla_batch(char *func_name, int error, int subproblem_ind);

#undef COMPLEX

#ifdef __cplusplus
}  // extern "C"
#endif    
#endif // CORE_Z_H 
