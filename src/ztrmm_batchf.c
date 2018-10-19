/**
 * @file ztrmm_batchf.c
 *
 * @brief BBLAS ztrmm_batchf double _Complex routine.
 *
 * BBLAS is a software package provided by 
 * Univ. of Manchester,
 * Univ. of Tennessee.
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
 * @ingroup trmm_batchf
 *
 *  Performs a triangular batch matrix-matrix multiply of the form
 *
 *          \f[B[i] = \alpha [op( A[i] ) \times B[i]] \f], if side = BblasLeft  or
 *          \f[B[i] = \alpha [B[i] \times op( A[i]) ] \f], if side = BblasRight
 *
 *  where op( X ) is one of:
 *
 *          - op(A[i]) = A[i]   or
 *          - op(A[i]) = A[i]^T or
 *          - op(A[i]) = A[i]^H
 *
 *  alpha is a scalar, B[i]-s are m-by-n matrices and A[i]-s are a unit or non-unit, upper
 *  or lower triangular matrix.
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
 * @param[in] side
 *          Specifies whether op( A[i] ) appears on the left or on the right of B[i]:
 *          - BblasLeft:  alpha*op( A[i] )*B[i]
 *          - BblasRight: alpha*B[i]*op( A[i] )
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
 * 		A is an array of pointers to matrices A[0], A[1] .. A[group_size-1], 
 * 		where each element A[i] is a pointer to a triangular matrix 
 *          	of dimension lda-by-k, where k is m when
 *          	side='L' or 'l' and k is n when when side='R' or 'r'. If uplo =
 *          	BblasUpper, the leading k-by-k upper triangular part of the array
 *          	A[i] contains the upper triangular matrix, and the strictly lower
 *          	triangular part of A[i] is not referenced. If uplo = BblasLower, the
 *          	leading k-by-k lower triangular part of the array A[i] contains the
 *          	lower triangular matrix, and the strictly upper triangular part of
 *          	A[i] is not referenced. If diag = BblasUnit, the diagonal elements of
 *          	A[i] are also not referenced and are assumed to be 1.
 *
 * @param[in] lda
 *          The leading dimension of the arrays A[i]. When side='L' or 'l',
 *          lda >= max(1,m), when side='R' or 'r' then lda >= max(1,n).
 *
 * @param[in] B
 * 		B is an array of pointers to matrices B[0], B[1],..,B[group_size-1],
 * 		where each element B[i] is a pointer to a matrix. 
 *          	On entry, the matrices B[i] are of dimension ldb-by-n.
 *          	On exit, the result of a triangular matrix-matrix multiply
 *          	( alpha*op(A[i])*B[i] ) or ( alpha*B[i]*op(A[i]) ).
 *
 * @param[in] ldb
 *          The leading dimension of the arrays B[i]. ldb >= max(1,m).
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

void blas_ztrmm_batchf( int group_size,
			bblas_enum_t layout, bblas_enum_t side, bblas_enum_t uplo,
    			bblas_enum_t transa, bblas_enum_t diag,
    			int m, int n, 
			bblas_complex64_t alpha, bblas_complex64_t const *const *A, int lda,
    			bblas_complex64_t **B, int ldb,
			 int *info)
{

	// Local variables 
	int first_index = 0;
	int batch_iter;
	int LDA;
	char func_name[15] = "ztrmm_batchf";

	// Check input arguments 
	if (group_size < 0)
	{
		xerbla_batch(func_name, BblasErrorBatchCount, -1);
	}
	else {
		if ((side != BblasLeft) && (side != BblasRight)) {
			xerbla_batch(func_name, BblasErrorSide, first_index);
			for (batch_iter = 0; batch_iter < group_size; batch_iter++) {
				info[batch_iter]  = BblasErrorSide;
			}
			return;
		}
		if ((uplo != BblasUpper) && (uplo != BblasLower)) {
			xerbla_batch(func_name, BblasErrorUplo, first_index);
			for (batch_iter = 0; batch_iter < group_size; batch_iter++) {
				info[batch_iter] = BblasErrorUplo;
			}
			return;
		}
		if ((transa != BblasNoTrans) &&
				(transa != BblasTrans) &&
				(transa != BblasConjTrans)) {
			xerbla_batch(func_name, BblasErrorTransa, first_index);
			for (batch_iter = 0; batch_iter < group_size; batch_iter++) {
				info[batch_iter]  = BblasErrorTransa;
			}
			return;
		}
		if ((diag != BblasNonUnit) && (diag != BblasUnit)) {
			xerbla_batch(func_name, BblasErrorDiag, first_index);
			for (batch_iter = 0; batch_iter < group_size; batch_iter++) {
				info[batch_iter]  = BblasErrorDiag;
			}
			return;
		}
		if (m < 0) {
			xerbla_batch(func_name, BblasErrorM, first_index);
			for (batch_iter = 0; batch_iter < group_size; batch_iter++) {
				info[batch_iter] = BblasErrorM;
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
		if (side == BblasLeft) {
			LDA = m;
		} 
		else {
			LDA = n;
		}
		if (lda < max(1, LDA)) {
			xerbla_batch(func_name, BblasErrorlda, first_index);
			for (batch_iter = 0; batch_iter < group_size; batch_iter++) {
				info[batch_iter] =  BblasErrorlda;
			}
			return;
		}
		if (ldb < max(1, m)) {
			xerbla_batch(func_name, BblasErrorldb, first_index);
			for (batch_iter = 0; batch_iter < group_size; batch_iter++) {
				info[batch_iter] = BblasErrorldb;
			}
			return;
		}
		// Skip subproblems where nothing needs to be done 
		if (min(m, n) == 0) {
			for (batch_iter = 0; batch_iter < group_size; batch_iter++) {
				info[batch_iter] =  BblasSuccess;
			}
			return;
		}
		for (batch_iter = 0; batch_iter < group_size; batch_iter++) {
			// Call to cblas_ztrmm 
			cblas_ztrmm(layout,
					side, uplo,
					transa, diag,
					m, n,
					CBLAS_SADDR(alpha), 
					A[batch_iter], lda,
					B[batch_iter], ldb);
			// Successful 
			info[batch_iter] = BblasSuccess;
		} // END FIXED SIZE FOR LOOP 
	}
}
#undef COMPLEX
