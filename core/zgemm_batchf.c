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

#include "cblas.h"
#include "bblas.h"

/***************************************************************************//**
 *
 * @ingroup gemm_batchf
 *
 * zgemm_batchf is a batch version of zgemm. It performs 
 * matrix-matrix multiplication of matrices, where all the
 * matrices of the batch have a fixed size.
 *
 *  \f[ C[i] = \alpha [op( A[i] )\times op( B[i] )] + \beta C[i], \f]
 *
 *  where op( X ) is one of:
 *    \f[ op( X ) = X,   \f]
 *    \f[ op( X ) = X^T, \f]
 *    \f[ op( X ) = X^H, \f]
 *
 *  alpha and beta are scalars, and A[i], B[i] and C[i] are matrices, with 
 *  op( A[i] ) an m-by-k matrix, op( B[i] ) a k-by-n matrix and C[i] an m-by-n matrix.
 *
 *******************************************************************************i
 * @param[in] group_size
 * 	    The number of matrices to operate on
 *
 * @param[in] layout
 * 	    Specifies if the matrix is stored in row major or column major
 * 	    format:
 * 	    - BblasRowMajor: Row major format
 * 	    - BblasColMajor: Column major format
 *
 * @param[in] transa
 *          - BblasNoTrans:   A[i] is not transposed,
 *          - BblasTrans:     A[i] is transposed,
 *          - BblasConjTrans: A[i] is conjugate transposed.
 *
 * @param[in] transb
 *          - BblasNoTrans:   B[i] is not transposed,
 *          - BblasTrans:     B[i] is transposed,
 *          - BblasConjTrans: B[i] is conjugate transposed.
 *
 * @param[in] m
 *          The number of rows of the matrix op( A[i] ) and of the matrix C[i].
 *          m >= 0.
 *
 * @param[in] n
 *          The number of columns of the matrix op( B[i] ) and of the matrix C[i].
 *          n >= 0.
 *
 * @param[in] k
 *          The number of columns of the matrix op( A[i] ) and the number of rows
 *          of the matrix op( B[i] ). k >= 0.
 *
 * @param[in] alpha
 *          The scalar alpha.
 *
 * @param[in] A
 * 		A is an array of pointers to matrices A[0], A[1] .. A[group_size-1], 
 * 		where each element A[i] is a pointer to a matrix of 
 * 		dimension lda-by-ka, where ka is k when transa = 
 * 		BblasNoTrans, and is m otherwise. When using transa = 
 * 		BblasNoTrans the leading m-by-k part of A[i] 
 * 		must contain the matrix elements, otherwise the leading  
 * 		k-by-m part of A[i] must contain the matrix elements.
 *
 * @param[in] lda
 *          The leading dimension of the array A[i].
 *          When transa = BblasNoTrans, lda >= max(1,m),
 *          otherwise, lda >= max(1,k).
 *
 * @param[in] B
 * 		B is an array of pointers to matrices B[0], B[1],..,B[group_size-1],
 * 		where each element B[i] is a pointer to a matrix of 
 * 		dimension lda-by-kb, where kb is n when transb = 
 * 		BblasNoTrans, and is k otherwise. When using transb = 
 * 		BblasNoTrans the leading k-by-n part of B[i] 
 * 		must contain the matrix elements, otherwise the leading  
 * 		n-by-k part of B[i] must contain the matrix elements.
 *
 * @param[in] ldb
 *          The leading dimension of the array B[i].
 *          When transb = BblasNoTrans, ldb >= max(1,k),
 *          otherwise, ldb >= max(1,n).
 *
 * @param[in] beta
 *          The scalar beta.
 *
 * @param[in,out] C
 * 		C is an array of pointers to matrices C[0], C[1],...,C[group_size-1],
 *		where each element of C[i] is a pointer to a matrix of dimension 
 *		ldc-by-n. On exit, each array C[i] is overwritten by the m-by-n
 *		matrix ( alpha*op(A[i] )*op( B[i] ) + beta*C[i] ).
 *
 * @param[in] ldc
 *          The leading dimension of the array C[i]. ldc >= max(1,m).
 *
 * @param[in,out] info
 * 		Array of int for error handling. On entry info[0] should have one of the 
 * 		following values
 *			- BblasErrorsReportAll    :  All errors will be specified on output.
 *						     Length of the array should be atleast
 *						     \sum_{i=0}^{group_count-1}group_size[i].
 *			- BblasErrorsReportGroup  :  Single error from each group will be 
 *						     reported. Length of the array should 
 *						     be atleast group_count.
 *			- BblasErrorsReportAny    :  Occurence of an error will be indicated
 *						     by a single integer value, and length 
 *						     of the array should be atleast 1.
 *			- BblasErrorsReportNone   :  No error will be reported on output, and
 *						     length of the array should be atleast 1.
 ********************************************************************************
 * @retval BblasSuccess successful exit
 *
 *******************************************************************************
 *
 * @sa zgemm_batchf
 * @sa cgemm_batchf
 * @sa dgemm_batchf
 * @sa sgemm_batchf
 *
 ******************************************************************************/
void blas_zgemm_batchf(int group_size, 
		       bblas_enum_t layout, bblas_enum_t transa, bblas_enum_t transb,
		       int m,  int n, int k,
		       bblas_complex64_t alpha, bblas_complex64_t const *const *A, int lda,
		       				bblas_complex64_t const* const *B, int ldb, 
		       bblas_complex64_t beta,  bblas_complex64_t** C, int ldc, 
		       int *info)
{
	// Local variables
	int iter;

	// Check input arguments 
	if ((layout != BblasRowMajor) &&
			(layout != BblasColMajor)) {
		bblas_error("Illegal value of layout");
		if (info[0] != BblasErrorsReportNone) {
			bblas_set_info(info[0], &info[0], group_size, 3);
		}
		return;
	}
	if ((transa != BblasNoTrans) &&
			(transa != BblasTrans) &&
			(transa != BblasConjTrans)) {
		bblas_error("Illegal value of transa");
		if (info[0] != BblasErrorsReportNone) {
			bblas_set_info(info[0], &info[0], group_size, 4);
		}
		return;
	}
	if ((transb != BblasNoTrans) &&
			(transb != BblasTrans) &&
			(transb != BblasConjTrans)) {
		bblas_error("Illegal value of transb");
		if (info[0] != BblasErrorsReportNone) {
			bblas_set_info(info[0], &info[0], group_size, 5);
		}
		return;
	}
	if ( transa == BblasNoTrans ) { 
		lda = m;
	} 
	else { 
		lda = k;
	}
	if ( transb == BblasNoTrans ) { 
		ldb = k;
	} 
	else { 
		ldb = n;
	}
	if (m < 0) {
		bblas_error("Illegal value of m");
		if (info[0] != BblasErrorsReportNone) {
			bblas_set_info(info[0], &info[0], group_size, 6);
		}
		return;
	}
	if (n < 0) {
		bblas_error("Illegal value of n");
		if (info[0] != BblasErrorsReportNone) {
			bblas_set_info(info[0], &info[0], group_size, 7);
		}
		return;
	}
	if (k < 0) {
		bblas_error("Illegal value of k");
		if (info[0] != BblasErrorsReportNone) {
			bblas_set_info(info[0], &info[0], group_size, 8);
		}
		return;
	}
	if (lda < imax(1, lda)) {
		bblas_error("Illegal value of lda");
		if (info[0] != BblasErrorsReportNone) {
			bblas_set_info(info[0], &info[0], group_size, 9);
		}
		return;
	}
	if (ldb < imax(1, ldb)) {
		bblas_error("Illegal value of ldb");
		if (info[0] != BblasErrorsReportNone) {
			bblas_set_info(info[0], &info[0], group_size, 10);
		}
		return;
	}
	if (ldc < imax(1, m)) {
		bblas_error("Illegal value of ldc");
		if (info[0] != BblasErrorsReportNone) {
			bblas_set_info(info[0], &info[0], group_size, 11);
		}
		return;
	}
	// Skip subproblems where nothing needs to be done
	if (m == 0 || n == 0 ||
			((alpha == (bblas_complex64_t)0.0 || k == 0) &&
			 beta == (bblas_complex64_t)1.0 )) {
		for (iter = 0; iter < group_size; iter++) {
			info[iter] =  0;
		}
		return;
	}
	for (iter = 0; iter < group_size; iter++) {
		// Call to cblas_zgemm 
		cblas_zgemm(layout,
			    transa, transb,
			    m, n, k,
			    CBLAS_SADDR(alpha),
			    A[iter], lda,
			    B[iter], ldb,
			    CBLAS_SADDR(beta),
			    C[iter], ldc);
		// Successful 
		info[iter] = 0;
	} // END FIXED SIZE FOR LOOP 

}
#undef COMPLEX
