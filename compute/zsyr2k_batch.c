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

#include "bblas.h"

/***************************************************************************//**
 *
 * @ingroup syr2k_batch 
 *
 *  Performs one of the symmetric rank 2k operations on a group of matrices, .
 *
 *    \f[ C[i] = \alpha A[i] \times B[i]^T + \alpha B[i] \times A[i]^T + \beta C[i], \f]
 *    or
 *    \f[ C[i] = \alpha A[i]^T \times B[i] + \alpha B[i]^T \times A[i] + \beta C[i], \f]
 *
 *  where alpha[i]-s and beta[i]-s are scalars,
 *  C[i]-s are n-by-n symmetric matrices and A[i]-s and B[i]-s are n-by-k matrices
 *  in the first case and k-by-n matrices in the second case.
 *
 *******************************************************************************
 *
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
 *                            symmetric matrices C[j] are to be stored.
 *          - BblasUpper:     Only the upper triangular part of the
 *                            symmetric matrices C[j] are to be stored.
 *
 * @param[in] trans
 * 	    An array of length group_count-1, where 
 * 	    for j-th matrix in i-th group
 *          - BblasNoTrans:
 *            \f[ C[i] = \alpha[i] A[j] \times B[j]^T + \alpha[i] B[j] \times A[j]^T + \beta[i] C[j]; \f]
 *          - BblasTrans:
 *            \f[ C[j] = \alpha[i] A[j]^T \times B[j] + \alpha[i] B[j]^T \times A[j] + \beta[i] C[j]. \f]
 *
 * @param[in] n
 *          An array of integers of length group count-1, where n[i] is 
 *          the order of the matrices C[j] in i-th group. n[i] >= 0.
 *
 * @param[in] k
 * 	    An array of integers of length group_count-1. For matrices
 * 	    in i-th group
 *          If trans = BblasNoTrans, number of columns of A[j]-s and B[j]-s matrices;
 *          if trans = BblasTrans, number of rows of A[j]-s and B[j]-s matrices.
 *
 * @param[in] alpha
 *          An array of scalars of length group_count-1.
 *
 * @param[in] A
 * 	    A is an array of pointers to matrices A[0], A[1] .. A[batch_count-1]. 
 * 	    In i-th group each element A[j] is a pointer to a matrix of 
 * 	    A[j] of size lda[i]-by-ka.
 *    	    If trans[i] = BblasNoTrans,   ka = k[i];
 *          if trans[i] = BblasTrans, ka = n[i].
 * 	    batch_count=\sum_{i=1}^{group_count}group_sizes[i].
 *
 *
 * @param[in] lda
 *          An array of integers of length group_count-1, 
 *          are the leading dimension of the arrays A[j]
 *          in i-th group.
 *          If trans[i] = BblasNoTrans,   lda[i] >= max(1, n[i]);
 *          if trans[i] = BblasTrans, lda[i] >= max(1, k[i]).
 *
 * @param[in] B
 * 		B is an array of pointers to matrices B[0], B[1] .. B[batch_count-1]. 
 * 		In i-th group each element B[j] is a pointer to a matrix B[j] of size 
 *          	ldb[i]-by-kb.
 *          	If trans[i] = BblasNoTrans,   kb = k[i];
 *          	if trans[i] = BblasTrans, kb = n[i].
 * 		batch_count=\sum_{i=1}^{group_count}group_sizes[i].
 *
 *
 * @param[in] ldb
 *          An array of integers of size group_count-1, are 
 *          the leading dimension of the arrays B[j] in i-th group.
 *          If trans[i] = BblasNoTrans,   ldb[i] >= max(1, n[i]);
 *          if trans[i] = BblasTrans, ldb[i] >= max(1, k[i]).
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
 * 		batch_count=\sum_{i=1}^{group_count}group_sizes[i].
 *
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
 *						     (group_count*group_sizes).
 *			- BblasErrorsReportGroup  :  Single error from each group will be 
 *						     reported. Length of the array should 
 *						     be atleast (group_count).
 *			- BblasErrorsReportAny    :  Occurence of an error will be indicated
 *						     by a single integer value, and length 
 *						     of the array should be atleast 1.
 *			- BblasErrorsReportNone   :  No error will be reported on output, and
 *						     length of the array should be atleast 1.
 ******************************************************************************
 *
 * @retval  BblasSuccess successful exit
 *
 *******************************************************************************
 *
 * @sa zsyr2k_batch
 * @sa csyr2k_batch
 * @sa dsyr2k_batch
 * @sa ssyr2k_batch
 *
 ******************************************************************************/
void blas_zsyr2k_batch(int group_count, const int *group_sizes, 
		       bblas_enum_t layout, const bblas_enum_t *uplo, const bblas_enum_t *trans,
		       const int *n, const int *k, 
		       const bblas_complex64_t *alpha, bblas_complex64_t const *const *A, const int *lda, 
		       				       bblas_complex64_t const* const *B, const int *ldb, 
		       const bblas_complex64_t  *beta, bblas_complex64_t	    ** C, const int *ldc, 
		       int *info)

{

	// Check input arguments 
	if (group_count < 0) {
		bblas_error("Illegal value of group_count");
		info[0] = -1;
		return;
	}

	int offset = 0;
	int info_offset = 0;
	int flag = 0;
	int info_option = info[0];
	// Check group_size and call fixed batch computation 
	for (int group_iter = 0; group_iter < group_count; group_iter++) {

		if (info_option == BblasErrorsReportAll) 
			info_offset = offset+1;
		else if (info_option == BblasErrorsReportGroup)
			info_offset = group_iter+1;	
		else 
			info_offset = 0;
		info[info_offset] = info_option;	

		if (group_sizes[group_iter] < 0) {
			bblas_error("Illegal values of group_sizes");
			info[0] = -2;
			return;
		}

		// Skip the group where nothing needs to be done
		if (n[group_iter] == 0 || k[group_iter] == 0 ||
				(alpha[group_iter] == (bblas_complex64_t)0.0 || 
				 beta[group_iter] == (bblas_complex64_t)1.0) || 
				group_sizes[group_iter] == 0) {
			bblas_success(info_option, &info[info_offset], group_sizes[group_iter]);
			continue;
		}

		// Call to blas_zher2k_batchf 
		blas_zsyr2k_batchf(group_sizes[group_iter], 
				   layout, uplo[group_iter], trans[group_iter],
				   n[group_iter], k[group_iter],
				   alpha[group_iter], A+offset, lda[group_iter],
				   		      B+offset, ldb[group_iter],
				   beta[group_iter],  C+offset, ldc[group_iter],
				   &info[info_offset]);    

		// check for errors in batchf function
		if (info[info_offset] != 0 && flag == 0) {
			info[0] = info[info_offset];	
			flag = 1;
		}

		offset += group_sizes[group_iter];    
	}
}
