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
 * @ingroup symm_batchf
 *
 *  Performs one of the matrix-matrix operations
 *
 *     \f[ C[i] = \alpha \times A[i] \times B[i] + \beta \times C[i] \f]
 *  or
 *     \f[ C[i] = \alpha \times B[i] \times A[i] + \beta \times C[i] \f]
 *
 *  where alpha and beta are scalars, A[i]-s are symmetric matrices and B[i]-s are
 *  C[i] are m-by-n matrices.
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
 * @param[in] side
 *          Specifies whether the symmetric matrices A[i] appears on the
 *          left or right in the operation as follows:
 *          - BblasLeft:  \f[ C[i] = \alpha \times A[i] \times B[i] + \beta \times C[i] \f]
 *          - BblasRight: \f[ C[i] = \alpha \times B[i] \times A[i] + \beta \times C[i] \f]
 *
 * @param[in] uplo
 *          Specifies whether the upper or lower triangular part of
 *          the symmetric matrices A[i] are to be referenced as follows:
 *          - BblasLower:     Only the lower triangular part of the
 *                             symmetric matrices A[i] is to be referenced.
 *          - BblasUpper:     Only the upper triangular part of the
 *                             symmetric matrices A[i] is to be referenced.
 *
 * @param[in] m
 *          The number of rows of the matrices C[i]. m >= 0.
 *
 * @param[in] n
 *          The number of columns of the matrices C[i]. n >= 0.
 *
 * @param[in] alpha
 *          The scalar alpha.
 *
 * @param[in] A
 * 		A is an array of pointers to matrices A[0], A[1] .. A[group_size-1], 
 * 		where each element A[i] is a pointer to a matrix A[i] of size
 *		lda-by-ka, where ka is m when side = BblasLeft,
 *          	and is n otherwise. Only the uplo triangular part is referenced.
 *
 * @param[in] lda
 *          The leading dimension of the arrays A[i]. lda >= max(1,ka).
 *
 * @param[in] B
 * 	    B is an array of pointers to matrices B[0], B[1] .. B[group_size-1], 
 * 	    where each element B[i] is a pointer to a matrix B[i] of size
 *          ldb-by-n matrix, where the leading m-by-n part of
 *          the array B[i] must contain the matrix B[i].
 *
 * @param[in] ldb
 *          The leading dimension of the arrays B[i]. ldb >= max(1,m).
 *
 * @param[in] beta
 *          The scalar beta.
 *
 * @param[in,out] C
 * 	    C is an array of pointers to matrices C[0], C[1] .. C[group_size-1], 
 * 	    where each element C[i] is a pointer to a matrix C[i].
 *          On exit, the array is overwritten by the m-by-n updated matrix.
 *
 * @param[in] ldc
 *          The leading dimension of the arrays C[i]. ldc >= max(1,m).
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
 *
 * @retval  BblasSuccess successful exit
 *
 *******************************************************************************
 *
 * @sa zsymm_batchf
 * @sa csymm_batchf
 * @sa dsymm_batchf
 * @sa ssymm_batchf
 *
 ******************************************************************************/
void blas_zsymm_batchf( int group_size,
		        bblas_enum_t layout, bblas_enum_t side, bblas_enum_t uplo,
    			int m,  int n, 
			bblas_complex64_t alpha, bblas_complex64_t const *const *A, int lda, 
    			 		 	 bblas_complex64_t const* const *B, int ldb, 
			bblas_complex64_t beta,  bblas_complex64_t** C, int ldc, 
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
	if ((side != BblasLeft) && (side != BblasRight)) {
		bblas_error("Illegal value of side");
		if (info[0] != BblasErrorsReportNone) {
			bblas_set_info(info[0], &info[0], group_size, 4);
		}
		return;
	}
	if ((uplo != BblasUpper) && (uplo != BblasLower)) {
		bblas_error("Illegal value of uplo");
		if (info[0] != BblasErrorsReportNone) {
			bblas_set_info(info[0], &info[0], group_size, 5);
		}
		return;
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
	if (side == BblasLeft) {
		LDA = m;
	} 
	else {
		LDA = n;
	}
	if (lda < LDA) {
		bblas_error("Illegal value of lda");
		if (info[0] != BblasErrorsReportNone) {
			bblas_set_info(info[0], &info[0], group_size, 8);
		}
		return;
	}
	if (ldb < imax(1, m)) {
		bblas_error("Illegal value of ldb");
		if (info[0] != BblasErrorsReportNone) {
			bblas_set_info(info[0], &info[0], group_size, 9);
		}
		return;
	}
	if (ldc < imax(1, m)) {
		bblas_error("Illegal value of ldc");
		if (info[0] != BblasErrorsReportNone) {
			bblas_set_info(info[0], &info[0], group_size, 10);
		}
		return;
	} 
	// Skip subproblems where nothing needs to be done 
	if (m == 0 || n == 0 ||
			(alpha == (bblas_complex64_t)0.0 &&
			 beta == (bblas_complex64_t)1.0)) {
		for (iter = 0; iter < group_size; iter++) {
			info[iter] = 0;
		}
		return;
	}
	for (iter = 0; iter < group_size; iter++) {
		// Call to cblas_zsymm 
		cblas_zsymm(layout,
			    side, uplo,
			    m, n,
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
