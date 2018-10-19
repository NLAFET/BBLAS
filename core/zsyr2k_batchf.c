/**
 * @file zsyr2k_batchf.c
 *
 *  @brief BBLAS zsyr2k_batchf  for double _Complex routine.
 *
 *  BBLAS is a software package provided by 
 *  Univ. of Manchester,
 *  Univ. of Tennessee.
 *
 * @author  Srikara Pranesh
 * @author  Mawussi Zounon
 * @date    2018-09-22
 *
 **/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
/**
 * Code generation
 * @precisions normal z -> c d s
 **/
#endif

#include "cblas.h"
#include "bblas_z.h"

#define COMPLEX

/***************************************************************************//**
 *
 * @ingroup syr2k_batchf
 *
 *  Performs one of the batch symmetric rank 2k operations
 *
 *    \f[ C[i] = \alpha A[i] \times B[i]^T + \alpha B[i] \times A[i]^T + \beta C[i], \f]
 *    or
 *    \f[ C[i] = \alpha A[i]^T \times B[i] + \alpha B[i]^T \times A[i] + \beta C[i], \f]
 *
 *  where alpha and beta are scalars,
 *  C[i]-s are n-by-n symmetric matrix, and A[i]-s and B[i]-s are n-by-k matrices
 *  in the first case and k-by-n matrices in the second case.
 *
 *******************************************************************************
 *
 * @param[in] group_size
 * 	    The number of matrices to operate on
 *
 * @param[in] layout
 * 	    Specifies if the matrix is stored in row major or column major
 * 	    format:
 * 	    - BblasRowMajor: Row major format
 * 	    - BblasColMajor: Column major format
 *
 * @param[in] uplo
 *          - BblasUpper: Upper triangle of C[i]-s are stored;
 *          - BblasLower: Lower triangle of C[i]-s are stored.
 *
 * @param[in] trans
 *          - BblasNoTrans:
 *            \f[ C[i] = \alpha A[i] \times B[i]^T + \alpha B[i] \times A[i]^T + \beta C[i]; \f]
 *          - BblasTrans:
 *            \f[ C[i] = \alpha A[i]^T \times B[i] + \alpha B[i]^T \times A[i] + \beta C[i]. \f]
 *
 * @param[in] n
 *          The order of the matrices C[i]. n >= zero.
 *
 * @param[in] k
 *          If trans = BblasNoTrans, number of columns of the A[i] and B[i] matrices;
 *          if trans = BblasTrans, number of rows of the A[i] and B[i] matrices.
 *
 * @param[in] alpha
 *          The scalar alpha.
 *
 * @param[in] A
 * 	    A is an array of pointers to matrices A[0], A[1] .. A[group_size-1], 
 * 	    where each element A[i] is a pointer to a matrix A[i] of size 
 *          lda-by-ka.
 *    	    If trans = BblasNoTrans,   ka = k;
 *          if trans = BblasTrans, ka = n.
 *
 * @param[in] lda
 *          The leading dimension of the arrays A[i].
 *          If trans = BblasNoTrans, lda >= max(1, n);
 *          if trans = BblasTrans,   lda >= max(1, k).
 *
 * @param[in] B
 * 		B is an array of pointers to matrices B[0], B[1] .. B[group_size-1], 
 * 		where each element B[i] is a pointer to a matrix B[i] of size 
 *          	ldb-by-kb.
 *          	If trans = BblasNoTrans,   kb = k;
 *          	if trans = BblasTrans, kb = n.
 *
 * @param[in] ldb
 *          The leading dimension of the arrays B[i].
 *          If trans = BblasNoTrans, ldb >= max(1, n);
 *          if trans = BblasTrans,   ldb >= max(1, k).
 *
 * @param[in] beta
 *          The scalar beta.
 *
 * @param[in,out] C
 * 		C is an array of pointers to matrices C[0], C[1] .. C[group_size-1], 
 * 		where each element C[i] is a pointer to a matrix C[i] of size 
 *          	ldc-by-n.
 *          	On exit, the uplo part of the matrix is overwritten
 *          	by the uplo part of the updated matrix.
 *
 * @param[in] ldc
 *          The leading dimension of the arrays C[i]. ldc >= max(1, n).
 *
 * @param[in,out] info
 * 		Array of int for error handling. On entry info[0] should have one of the 
 * 		following values
 *			- BblasErrorsReportAll    :  All errors will be specified on output.
 *						     Length of the array should be atleast
 *						     (group_count*group_size).
 *			- BblasErrorsReportGroup  :  Single error from each group will be 
 *						     reported. Length of the array should 
 *						     be atleast (group_count).
 *			- BblasErrorsReportAny    :  Occurence of an error will be indicated
 *						     by a single integer value, and length 
 *						     of the array should be atleast 1.
 *			- BblasErrorsReportNone   :  No error will be reported on output, and
 *						     length of the array should be atleast 1.
 ******************************************************************************/

void blas_zsyr2k_batchf( int group_size, 
			bblas_enum_t layout, bblas_enum_t uplo, bblas_enum_t trans,
    			int n, int k, 
			bblas_complex64_t alpha, bblas_complex64_t const *const *A, int lda, 
    			 		 	 bblas_complex64_t const* const *B, int ldb, 
			const double  beta,  bblas_complex64_t** C, int ldc, 
    			 int *info)
{
	// Local variables 
	int first_index = 0;
	int batch_iter;
	int LDA,  LDB;
	char func_name[15] = "zsyr2k_batchf";

	// Check input arguments 
	if (group_size < 0) {
		xerbla_batch(func_name, BblasErrorBatchCount, -1);
	}
	else {
		if ((uplo != BblasUpper) && (uplo != BblasLower)) {
			xerbla_batch(func_name, BblasErrorUplo, first_index);
			for (batch_iter = 0; batch_iter < group_size; batch_iter++) {
				info[batch_iter] = BblasErrorUplo;
			}
			return;
		}
		if ((trans != BblasNoTrans) &&
				(trans != BblasTrans) && (trans != BblasConjTrans)) {
			xerbla_batch(func_name, BblasErrorTrans, first_index);
			for (batch_iter = 0; batch_iter < group_size; batch_iter++) {
				info[batch_iter]  = BblasErrorTrans;
			}
			return;
		}
		if (n < 0) {
			xerbla_batch(func_name, BblasErrorN, first_index);
			for (batch_iter = 0; batch_iter < group_size; batch_iter++) {
				info[batch_iter] = BblasErrorN;
			}
			return;
		}
		if (k < 0) {
			xerbla_batch(func_name, BblasErrorK, first_index);
			for (batch_iter = 0; batch_iter < group_size; batch_iter++) {
				info[batch_iter] = BblasErrorK;
			}
			return;
		}
		if (trans == BblasNoTrans) {
			LDA = n;
			LDB = n;
		} 
		else {
			LDA = k;
			LDB = k;
		}
		if (lda < max(1,LDA)) {
			xerbla_batch(func_name, BblasErrorlda, first_index);
			for (batch_iter = 0; batch_iter < group_size; batch_iter++) {
				info[batch_iter] =  BblasErrorlda;
			}
			return;
		}
		if (ldb < max(1, LDB)) {
			xerbla_batch(func_name, BblasErrorldb, first_index);
			for (batch_iter = 0; batch_iter < group_size; batch_iter++) {
				info[batch_iter] = BblasErrorldb;
			}
			return;
		}
		if (ldc < max(1, n)) {
			xerbla_batch(func_name, BblasErrorldc, first_index);
			for (batch_iter = 0; batch_iter < group_size; batch_iter++) {
				info[batch_iter] = BblasErrorldc;
			}
			return;
		}
		// Skip subproblems where nothing needs to be done
		if (n == 0 || k == 0 ||
				(alpha == (bblas_complex64_t)0.0 || 
				 beta == (bblas_complex64_t)1.0)) {
			for (batch_iter = 0; batch_iter < group_size; batch_iter++) {
				info[batch_iter] =  BblasSuccess;
			}
			return;
		}
		for (batch_iter = 0; batch_iter < group_size; batch_iter++) {
			// Call to cblas_zsyr2k
			cblas_zsyr2k(layout,
					uplo, trans,
					n, k,
					CBLAS_SADDR(alpha),
					A[batch_iter], lda,
					B[batch_iter], ldb,
					CBLAS_SADDR(beta),
					C[batch_iter], ldc);
			// Successful 
			info[batch_iter] = BblasSuccess;
		} // END FIXED SIZE FOR LOOP 
	}
}
#undef COMPLEX
