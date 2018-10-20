/**
 *
 * @file
 *
 *  BBLAS is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @precisions normal z -> c
 *
 **/

#include "cblas.h"
#include "bblas.h"

/***************************************************************************//**
 *
 * @ingroup her2k_batch
 *
 *  Performs one of the batch Hermitian rank 2k operations 
 *
 *    \f[ C[i] = \alpha[i] A[i] \times B[i]^H + conjg( \alpha[i] ) B \times A[i]^H + \beta[i] C[i], \f]
 *    or
 *    \f[ C[i] = \alpha[i] A[i]^H \times B[i] + conjg( \alpha[i] ) B[i]^H \times A[i] + \beta[i] C[i], \f]
 *
 *  for a group of matrices, where alpha[i] is a complex scalar, beta[i] is a real scalar,
 *  C[i]-s are n-by-n Hermitian matrices, and A[i]-s and B[i]-s are n-by-k matrices
 *  in the first case and k-by-n matrices in the second case.
 *
 *******************************************************************************
 * @param[in] group_count
 * 	    The number groups of matrices.	  
 *
 * @param[in] group_sizes
 * 	    An array of integers of length group_count-1, where group_sizes[i] denotes
 * 	    the number of matrices in i-th group.	
 * 
 * @param[in] layout
 * 	    Specifies if the matrix is stored in row major or column major
 * 	    format:
 * 	    - BblasRowMajor: Row major format
 * 	    - BblasColMajor: Column major format
 *
 * @param[in] uplo
 * 	    An array of length group_count-1, where uplo[i]
 *          specifies whether the upper or lower triangular part of
 *          the Hermitian matrices C[j]-s of i-th group are to be stored
 *  
 *          - BblasLower:     Only the lower triangular part of the
 *                             Hermitian matrices C[j] are to be stored.
 *          - BblasUpper:     Only the upper triangular part of the
 *                             Hermitian matrices C[j] are to be stored.
 *
 * @param[in] trans
 * 	    An array of size group_count-1, where 
 * 	    for j-th matrix in i-th group
 *          - BblasNoTrans:
 *            \f[ C[j] = \alpha[i] A[j] \times B[j]^H
 *                  + conjg( \alpha[i] ) B[j] \times A[j]^H + \beta[i] C[j]; \f]
 *          - BblasConjTrans:
 *            \f[ C[j] = \alpha[i] A[j]^H \times B[j]
 *                  + conjg( \alpha[i] ) B[j]^H \times A[j] + \beta[i] C[j]. \f]
 *
 * @param[in] n
 *          An array of integers of size group count-1, where n[i] is 
 *          the order of the matrices C[j] in i-th group. n[i] >= 0.
 *
 * @param[in] k
 * 	    An array of integers of length group_count-1. For matrices
 * 	    in i-th group
 *          If trans = BblasNoTrans, number of columns of A[j]-s and B[j]-s matrices;
 *          if trans = BblasConjTrans, number of rows of A[j]-s and B[j]-s matrices.
 *
 * @param[in] alpha
 *          An array of scalars of length group_count-1.
 *
 * @param[in] A
 * 	    A is an array of pointers to matrices A[0], A[1] .. A[batch_count-1]. 
 * 	    In i-th group each element A[j] is a pointer to a matrix of 
 * 	    A[j] of size lda[i]-by-ka.
 *    	    If trans[i] = BblasNoTrans,   ka = k[i];
 *          if trans[i] = BblasConjTrans, ka = n[i].
 *
 * @param[in] lda
 *          An array of integers of length group_count-1, 
 *          are the leading dimension of the arrays A[j]
 *          in i-th group.
 *          If trans[i] = BblasNoTrans,   lda[i] >= max(1, n[i]);
 *          if trans[i] = BblasConjTrans, lda[i] >= max(1, k[i]).
 *
 * @param[in] B
 * 		B is an array of pointers to matrices B[0], B[1] .. B[batch_count-1]. 
 * 		In i-th group each element B[j] is a pointer to a matrix B[j] of size 
 *          	ldb[i]-by-kb.
 *          	If trans[i] = BblasNoTrans,   kb = k[i];
 *          	if trans[i] = BblasConjTrans, kb = n[i].
 *
 * @param[in] ldb
 *          An array of integers of size group_count-1, are 
 *          the leading dimension of the arrays B[j] in i-th group.
 *          If trans[i] = BblasNoTrans,   ldb[i] >= max(1, n[i]);
 *          if trans[i] = BblasConjTrans, ldb[i] >= max(1, k[i]).
 *
 * @param[in] beta
 *          An array of scalars of size group_count-1.
 *
 * @param[in,out] C
 * 		C is an array of pointers to matrices C[0], C[1] .. C[batc_count-1].
 * 		In i-th group each element C[j] is a pointer to a matrix C[j] of size
 *          	ldc[i]-by-n[i].
 *          	On exit, the uplo[i] part of the matrix is overwritten
 *          	by the uplo[i] part of the updated matrix.
 *
 * @param[in] ldc
 *          An array of integers of length group_count-1. Where ldc[i]
 *          is the leading dimension of the arrays C[j] in i-th group. 
 *          ldc[i] >= max(1, n[i]).
 *
 * @param[in,out] info
 * 		Array of int for error handling. On entry info[0] should have one of the 
 * 		following values
 *			- BblasErrorsReportAll    :  All errors will be specified on output.
 *						     Length of the array should be atleast
 *						     \sum_{i=1}^{group_count-1}group_sizes[i]..
 *			- BblasErrorsReportGroup  :  Single error from each group will be 
 *						     reported. Length of the array should 
 *						     be atleast (group_count).
 *			- BblasErrorsReportAny    :  Occurence of an error will be indicated
 *						     by a single integer value, and length 
 *						     of the array should be atleast 1.
 *			- BblasErrorsReportNone   :  No error will be reported on output, and
 *						     length of the array should be atleast 1.
 *******************************************************************************
 *
 * @retval BblasSuccess successful exit
 *
 *******************************************************************************
 *
 * @sa zher2k_batch
 * @sa cher2k_batch
 *
 ******************************************************************************/
void blas_zher2k_batch( int group_count, const int *group_sizes,
			bblas_enum_t layout, const bblas_enum_t *uplo, const bblas_enum_t *trans,
    			const int *n, const int *k, 
			const bblas_complex64_t *alpha, bblas_complex64_t const *const *A, const int *lda, 
    			 		 	 bblas_complex64_t const* const *B, const int *ldb, 
			const double  *beta,  bblas_complex64_t** C, const int *ldc, 
    			int *info)

{
	// Local variables 
	int group_iter;
	int offset = 0;
	int info_offset = offset;

	// Check input arguments 
	if (group_count < 0) {
		bblas_error("Illegal value of group_count");
		info[0] = 1;
		return;
	}

	// Check group_size and call fixed batch computation 
	for (group_iter = 0; group_iter < group_count; group_iter++) {
		if (group_sizes[group_iter] < 0) {
			bblas_error("Illegal values of group_sizes");
			if (info[0] != BblasErrorsReportNone) {
				bblas_set_info(info[0], &info[0], group_sizes[group_iter], 2);
			}
			return;
		}

		if (group_iter != 0) {
			if (info[0] == BblasErrorsReportAll) 
				info_offset = offset;
			else
				info_offset = group_iter;
		}
		info[info_offset] = info[0];	

		// Call to bblas_zher2k_batchf 
		blas_zher2k_batchf (group_sizes[group_iter], 
				    layout,
				    uplo[group_iter],
				    trans[group_iter],
				    n[group_iter],
				    k[group_iter],
				    alpha[group_iter],
				    A+offset,
				    lda[group_iter],
				    B+offset,
				    ldb[group_iter],
				    beta[group_iter],
				    C+offset,
				    ldc[group_iter],
				    &info[info_offset]);    

		offset += group_sizes[group_iter];    
	}
}
#undef COMPLEX
