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
#include "bblas.h"


/***************************************************************************//**
 *
 * @ingroup her2k_batchf
 *
 *  Performs one of the batch Hermitian rank 2k operations
 *
 *    \f[ C[i] = \alpha A[i] \times B[i]^H + conjg( \alpha ) B \times A[i]^H + \beta C[i], \f]
 *    or
 *    \f[ C[i] = \alpha A[i]^H \times B[i] + conjg( \alpha ) B[i]^H \times A[i] + \beta C[i], \f]
 *
 *  where alpha is a complex scalar, beta is a real scalar,
 *  C[i]-s are n-by-n Hermitian matrices, and A[i]-s and B[i]-s are n-by-k matrices
 *  in the first case and k-by-n matrices in the second case.
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
 *
 * @param[in] uplo
 *          - BblasUpper: Upper triangle of C[i]-s is stored;
 *          - BblasLower: Lower triangle of C[i]-s is stored.
 *
 * @param[in] trans
 *          - BblasNoTrans:
 *            \f[ C[i] = \alpha A[i] \times B[i]^H
 *                  + conjg( \alpha ) B[i] \times A[i]^H + \beta C[i]; \f]
 *          - BblasConjTrans:
 *            \f[ C[i] = \alpha A[i]^H \times B[i]
 *                  + conjg( \alpha ) B[i]^H \times A[i] + \beta C[i]. \f]
 *
 * @param[in] n
 *          The order of the matrix C[i]. n >= zero.
 *
 * @param[in] k
 *          If trans = BblasNoTrans, number of columns of the A[i] and B[i] matrices;
 *          if trans = BblasConjTrans, number of rows of the A[i] and B[i] matrices.
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
 * @param[in] B
 * 		B is an array of pointers to matrices B[0], B[1] .. B[group_size-1], 
 * 		where each element B[i] is a pointer to a matrix B[i] of size 
 *          	ldb-by-kb.
 *          	If trans = BblasNoTrans,   kb = k;
 *          	if trans = BblasConjTrans, kb = n.
 *
 * @param[in] ldb
 *          The leading dimension of the arrays B[i].
 *          If trans = BblasNoTrans,   ldb >= max(1, n);
 *          if trans = BblasConjTrans, ldb >= max(1, k).
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
 *						     Length of the array should be at least
 *						     (group_count*group_size).
 *			- BblasErrorsReportGroup  :  Single error from each group will be 
 *						     reported. Length of the array should 
 *						     be at least (group_count).
 *			- BblasErrorsReportAny    :  Occurence of an error will be indicated
 *						     by a single integer value, and length 
 *						     of the array should be at least 1.
 *			- BblasErrorsReportNone   :  No error will be reported on output, and
 *						     length of the array should be at least 1.
 *******************************************************************************
 *
 * @retval BblasSuccess successful exit
 *
 *******************************************************************************
 *
 * @sa zher2k_batchf
 * @sa cher2k_batchf
 *
 ******************************************************************************/
void blas_zher2k_batchf(int group_size, bblas_enum_t layout, bblas_enum_t uplo,
                        bblas_enum_t trans, int n, int k,
                        bblas_complex64_t alpha, bblas_complex64_t const *const *A, int lda,
                                                 bblas_complex64_t const* const *B, int ldb,
                        const double  beta,      bblas_complex64_t             **C, int ldc,
                        int *info)
{
	// Check input arguments 
	if ((layout != BblasRowMajor) &&
        (layout != BblasColMajor)) {
		bblas_error("Illegal value of layout");
		if (info[0] != BblasErrorsReportNone) {
			bblas_set_info(info[0], &info[0], group_size, 1);
		}
		return;
	}
	if ((uplo != BblasUpper) && (uplo != BblasLower)) {
		bblas_error("Illegal value of uplo");
		if (info[0] != BblasErrorsReportNone) {
			bblas_set_info(info[0], &info[0], group_size, 2);
		}
		return;
	}
	if ((trans != BblasNoTrans) &&
        (trans != BblasTrans) && (trans != BblasConjTrans)) {
		bblas_error("Illegal value of trans");
		if (info[0] != BblasErrorsReportNone) {
			bblas_set_info(info[0], &info[0], group_size, 3);
		}
		return;
	}
	if (n < 0) {
		bblas_error("Illegal value of n");
		if (info[0] != BblasErrorsReportNone) {
			bblas_set_info(info[0], &info[0], group_size, 4);
		}
		return;
	}
	if (k < 0) {
		bblas_error("Illegal value of k");
		if (info[0] != BblasErrorsReportNone) {
			bblas_set_info(info[0], &info[0], group_size, 5);
		}
		return;
	}
    int am, bm;
	if (trans == BblasNoTrans) {
		am = n;
		bm = n;
	} 
	else {
		am = k;
		bm = k;
	}
	if (lda < imax(1, am)) {
		bblas_error("Illegal value of lda");
		if (info[0] != BblasErrorsReportNone) {
			bblas_set_info(info[0], &info[0], group_size, 6);
		}
		return;
	}
	if (ldb < imax(1, bm)) {
		bblas_error("Illegal value of ldb");
		if (info[0] != BblasErrorsReportNone) {
			bblas_set_info(info[0], &info[0], group_size, 7);
		}
		return;
	}
	if (ldc < imax(1, n)) {
		bblas_error("Illegal value of ldb");
		if (info[0] != BblasErrorsReportNone) {
			bblas_set_info(info[0], &info[0], group_size, 8);
		}
		return;
	}
	for (int iter = 0; iter < group_size; iter++) {
		cblas_zher2k(layout, uplo, trans,
                     	     n, k,
			     CBLAS_SADDR(alpha), A[iter], lda,
			     			 B[iter], ldb,
			     beta, C[iter], ldc);
		// BblasSuccess
		if (info[0] == BblasErrorsReportAll)
			info[iter] = 0;
	}
	// BblasSuccess
	if (info[0] != BblasErrorsReportAll)
		info[0] = 0;
}
