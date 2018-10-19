/**
 * @file zherk_batchf.c
 *
 * @brief BBLAS zherk_batchf double _Complex routine.
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
 * @precisions normal z -> c
 **/
#endif


#include "cblas.h"
#include "bblas.h"

#define COMPLEX

/***************************************************************************//**
 *
 * @ingroup herk_batchf
 *
 *  Performs one of the batch Hermitian rank k operations
 *
 *    \f[ C[i] = \alpha A[i] \times A[i]^H + \beta C[i], \f]
 *    or
 *    \f[ C[i] = \alpha A[i]^H \times A[i] + \beta C[i], \f]
 *
 *  where alpha and beta are real scalars, C[i]-s are n-by-n Hermitian
 *  matrices, and A[i]-s are n-by-k matrices in the first case and k-by-n
 *  matrices in the second case.
 *
 *******************************************************************************
 * @param[in] group_size
 * 	    The number of matrices to operate on
 *
 * @param[in] layout
 * 	    Specifies if the matrix is stored in row major or column major
 * 	    format:
 * 	    - BblasRowMajor: Row major format
 * 	    - BblasColMajor: Column major format
 * @param[in] uplo
 *          - BblasUpper: Upper triangle of C[i]-s are stored;
 *          - BblasLower: Lower triangle of C[i]-s are stored.
 *
 * @param[in] trans
 *          - BblasNoTrans:   \f[ C[i] = \alpha A[i] \times A[i]^H + \beta C[i]. \f]
 *          - BblasConjTrans: \f[ C[i] = \alpha A[i]^H \times A[i] + \beta C[i]. \f]
 *
 * @param[in] n
 *          The order of the matrices C[i]. n >= 0.
 *
 * @param[in] k
 *          If trans = BblasNoTrans, number of columns of the matrices A[i];
 *          if trans = BblasConjTrans, number of rows of the matrices A[i].
 *
 * @param[in] alpha
 *          The scalar alpha.
 *
 * @param[in] A
 * 	    A is an array of pointers to matrices A[0], A[1] .. A[group_size-1], 
 * 	    where each element A[i] is a pointer to a matrix A[i] of size 
 *          lda-by-ka.
 *    	    If trans = BblasNoTrans,   ka = k;
 *          if trans = BblasConjTrans, ka = n.
 *
 * @param[in] lda
 *          The leading dimension of the arrays A[i].
 *          If trans = BblasNoTrans,   lda >= max(1, n);
 *          if trans = BblasConjTrans, lda >= max(1, k).
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
 *						     be atleast to (group_count).
 *			- BblasErrorsReportAny    :  Occurence of an error will be indicated
 *						     by a single integer value, and length 
 *						     of the array should be atleast 1.
 *			- BblasErrorsReportNone   :  No error will be reported on output, and
 *						     length of the array should be atleast 1.
 ******************************************************************************/

void blas_zherk_batchf( int group_size,
		      bblas_enum_t layout, bblas_enum_t uplo, bblas_enum_t trans,
		      int n, int k, 
		      const double alpha, bblas_complex64_t const *const *A, int lda, 
		      const double  beta,  bblas_complex64_t** C, int ldc, 
    		      int *info)
{
	// Local variables 
	int iter;
	int LDA;

	// Check input arguments 
	if ((layout != CblasRowMajor) &&
			(layout != CblasColMajor)) {
		bblas_error("Illegal value of layout");
		if (info[0] != BblasErrorsReportNone) {
			bblas_set_info(info[0], &info[0], group_size, 3);
		}
		return;
	}
	if ((uplo != BblasUpper) && (uplo != BblasLower)) {
		bblas_error("Illegal value of uplo");
		if (info[0] != BblasErrorsReportNone) {
			bblas_set_info(info[0], &info[0], group_size, 4);
		}
		return;
	}
	if ((trans != BblasNoTrans) &&
			(trans != BblasTrans) && (trans != BblasConjTrans)) {
		bblas_error("Illegal value of trans");
		if (info[0] != BblasErrorsReportNone) {
			bblas_set_info(info[0], &info[0], group_size, 5);
		}
		return;
	}
	if (n < 0) {
		bblas_error("Illegal value of n");
		if (info[0] != BblasErrorsReportNone) {
			bblas_set_info(info[0], &info[0], group_size, 6);
		}
		return;
	}
	if (k < 0) {
		bblas_error("Illegal value of k");
		if (info[0] != BblasErrorsReportNone) {
			bblas_set_info(info[0], &info[0], group_size, 7);
		}
		return;
	}
	if (trans == BblasNoTrans) {
		LDA = n;
	} 
	else {
		LDA = k;
	}
	if (lda < imax(1, LDA)) {
		bblas_error("Illegal value of lda");
		if (info[0] != BblasErrorsReportNone) {
			bblas_set_info(info[0], &info[0], group_size, 8);
		}
		return;
	}
	if (ldc < imax(1, n))
	{
		bblas_error("Illegal value of ldc");
		if (info[0] != BblasErrorsReportNone) {
			bblas_set_info(info[0], &info[0], group_size, 9);
		}
		return;
	}
	// Skip subproblems where nothing needs to be done 
	if (n == 0   || ((k == 0 || alpha == (double)0.0) &&
				(beta == (double)1.0))) {
		for (iter = 0; iter < group_size; iter++) {
			info[iter] =  0;
		}
		return;
	}
	for (iter = 0; iter < group_size; iter++) {
		// Call to cblas_zherk 
		cblas_zherk(layout,
				uplo, trans,
				n, k,
				alpha,
				A[iter], lda,
				beta,
				C[iter], ldc);
		// Successful 
		info[iter] = 0;
	} // END FIXED SIZE FOR LOOP 
}
#undef COMPLEX
