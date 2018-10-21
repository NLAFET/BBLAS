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
 * @ingroup trsm_batchf
 *
 *  Solves one of the batch matrix equations
 *
 *    \f[ op( A[i] )\times X[i]  = \alpha B[i], \f] or
 *    \f[ X[i] \times op( A[i] ) = \alpha B[i], \f]
 *
 *  where op( A[i] ) is one of:
 *    \f[ op( A[i] ) = A[i],   \f]
 *    \f[ op( A[i] ) = A[i]^T, \f]
 *    \f[ op( A[i] ) = A[i]^H, \f]
 *
 *  alpha is a scalar, X[i]-s and B[i]-s are m-by-n matrices, and
 *  A[i]-s are unit or non-unit, upper or lower triangular matrices.
 *  The matrix X[i] overwrites B[i].
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
 *          Specifies whether op( A[i] ) appears on the left or on the right of X[i]:
 *          - BblasLeft:  op( A[i] )*X[i] = B[i],
 *          - BblasRight: X[i]*op( A[i] ) = B[i].
 *
 * @param[in] uplo
 *          Specifies whether the matrices A[i]-s are upper triangular or lower
 *          triangular:
 *          - BblasUpper: Upper triangle of A[i] is stored;
 *          - BblasLower: Lower triangle of A[i] is stored.
 *
 * @param[in] transa
 *          Specifies whether the matrices A[i] are transposed, not transposed or
 *          conjugate transposed:
 *          - BblasNoTrans:   A[i]-s are transposed;
 *          - BblasTrans:     A[i]-s are not transposed;
 *          - BblasConjTrans: A[i]-s are conjugate transposed.
 *
 * @param[in] diag
 *          Specifies whether or not A[i]-s are unit triangular:
 *          - BblasNonUnit: A[i]-s are non-unit triangular;
 *          - BblasUnit:    A[i]-s are unit triangular.
 *
 * @param[in] m
 *          The number of rows of matrices B[i].
 *          m >= 0.
 *
 * @param[in] n
 *          The number of columns of matrices B[i].
 *          n >= 0.
 *
 * @param[in] alpha
 *          The scalar alpha.
 *
 * @param[in] A
 * 	    A is an array of pointers to matrices A[0], A[1] .. A[group_size-1], 
 * 	    where each element A[i] is a pointer to a triangular matrix of dimension
 *          lda-by-ka triangular,
 *          where ka = m if side = BblasLeft,
 *            and ka = n if side = BblasRight.
 *          If uplo = BblasUpper, the leading k-by-k upper triangular part
 *          of the array A[i] contains the upper triangular matrix, and the
 *          strictly lower triangular part of A[i] is not referenced.
 *          If uplo = BblasLower, the leading k-by-k lower triangular part
 *          of the array A[i] contains the lower triangular matrix, and the
 *          strictly upper triangular part of A[i] is not referenced.
 *          If diag = BblasUnit, the diagonal elements of A[i] are also not
 *          referenced and are assumed to be 1.
 *
 * @param[in] lda
 *          The leading dimension of the array A. lda >= max(1,k).
 *
 * @param[in,out] B
 * 	    B is an array of pointers to matrices B[0], B[1] .. B[group_size-1], 
 *          On entry, each B[i]-s are ldb-by-n right hand side matrix.
 *          On exit, if return value = 0, the ldb-by-n solution matrix X.
 *
 * @param[in] ldb
 *          The leading dimension of the array B. ldb >= max(1,m).
 *
 * @param[in,out] info
 * 		Array of int for error handling. On entry info[0] should have one of the 
 * 		following values
 *			- BblasErrorsReportAll    :  All errors will be specified on output.
 *						     Length of the array should be atleast
 *						     \sum_{i=1}^{group_count-1}group_sizes[i].
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
 * @retval BblasSuccess successful exit
 *
 *******************************************************************************
 *
 * @sa ztrsm_batchf
 * @sa ctrsm_batchf
 * @sa dtrsm_batchf
 * @sa strsm_batchf
 *
 ******************************************************************************/
void blas_ztrsm_batchf(int group_size, bblas_enum_t layout, bblas_enum_t side,
                       bblas_enum_t uplo, bblas_enum_t transa, bblas_enum_t diag,
                       int m, int n, 
                       bblas_complex64_t alpha, bblas_complex64_t const *const *A, int lda,
                                                bblas_complex64_t             **B, int ldb,
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
	if ((side != BblasLeft) && (side != BblasRight)) {
		bblas_error("Illegal value of side");
		if (info[0] != BblasErrorsReportNone) {
			bblas_set_info(info[0], &info[0], group_size, 2);
		}
		return;
	}
	if ((uplo != BblasUpper) && (uplo != BblasLower)) {
		bblas_error("Illegal value of uplo");
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
	if ((diag != BblasNonUnit) && (diag != BblasUnit)) {
		bblas_error("Illegal value of diag");
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
    int an;
	if (side == BblasLeft) {
		an = m;
	} 
	else {
		an = n;
	}
    if (lda < imax(1, an)) {
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
    for (int iter = 0; iter < group_size; iter++) {
	    cblas_ztrsm(layout, side, uplo,
			transa, diag,
			m, n,
			CBLAS_SADDR(alpha), A[iter], lda,
					    B[iter], ldb);
		// BblasSuccess
		if (info[0] == BblasErrorsReportAll)
			info[iter] = 0;
	}
	// BblasSuccess
	if (info[0] != BblasErrorsReportAll)
		info[0] = 0;
}
