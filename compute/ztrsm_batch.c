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
 * @ingroup trsm_batch
 *
 *  Solves one of the batch matrix equations
 *
 *    \f[ op( A[i] )\times X[i]  = \alpha[i] B[i], \f] or
 *    \f[ X[i] \times op( A[i] ) = \alpha[i] B[i], \f]
 *
 *  for a group of matrices, and op( A[i] ) is one of:
 *
 *    \f[ op( A[i] ) = A[i],   \f]
 *    \f[ op( A[i] ) = A[i]^T, \f]
 *    \f[ op( A[i] ) = A[i]^H, \f]
 *
 *  alpha[i] is a scalar, X[i]-s and B[i]-s are m-by-n matrices, and
 *  A[i]-s are unit or non-unit, upper or lower triangular matrices.
 *  The matrix X[i] overwrites B[i].
 *
 *******************************************************************************
 *
 * @param[in] group_count
 * 	    The number groups of matrices.	  
 *
 * @param[in] group_sizes
 * 	    An array of integers of length group_count-1
 * 	    denoting the number of matrices in each group.	
 * 
 * @param[in] layout
 * 	    Specifies if the matrix is stored in row major or column major
 * 	    format:
 * 	    - BblasRowMajor: Row major format
 * 	    - BblasColMajor: Column major format
 *
 * @param[in] side
 * 	    An array of length group_count-1, where for matrices of i-th group
 *          specifies whether op( A[j] ) appears on the left or on the right of X[j]:
 *          - BblasLeft:  op( A[j] )*X[j] = B[j],
 *          - BblasRight: X[j]*op( A[j] ) = B[j].
 *
 * @param[in] uplo
 * 	    An array of length group_count-1, where for matrices of i-th group
 *          specifies whether the matrices A[j]-s are upper triangular or lower
 *          triangular:
 *          - BblasUpper: Upper triangle of A[j] is stored;
 *          - BblasLower: Lower triangle of A[j] is stored.
 *
 * @param[in] transa
 * 	    An array of length group_count-1, where
 *          - BblasNoTrans:   A[j]-s in i-th group are not transposed,
 *          - BblasTrans:     A[j]-s in i-th group are transposed,
 *          - BblasConjTrans: A[j]-s in i-th group are conjugate transposed.
 *
 * @param[in] diag
 *          An array of length group_count-1, which specifies 
 *          whether or not A[j]-s of i-th group are unit triangular:
 *          - BblasNonUnit: A[j]-s are non-unit triangular;
 *          - BblasUnit:    A[j]-s are unit triangular.
 *
 * @param[in] m
 *          An array integers of length group_count-1,
 *          which specified the number of rows of matrices B[j]
 *          in i-th group. m[i] >= 0.
 *
 * @param[in] n
 *          An array integers of length group_count-1,
 *          which specified the number of columns of matrices B[j]
 *          in i-th group. n[i] >= 0.
 *
 * @param[in] alpha
 *          An array of length group_count-1, where alpha[i] is
 *          a scalar.
 *
 * @param[in] A
 * 	    A is an array of pointers to matrices A[0], A[1] .. A[batch_count-1], 
 * 	    where for i-th group each element A[j] is a pointer to a triangular matrix of dimension
 *          lda[i]-by-ka triangular,
 *          where ka = m[i] if side = BblasLeft,
 *            and ka = n[i] if side = BblasRight.
 *          If uplo[i] = BblasUpper, the leading k-by-k upper triangular part
 *          of the array A[j] contains the upper triangular matrix, and the
 *          strictly lower triangular part of A[j] is not referenced.
 *          If uplo[i] = BblasLower, the leading k-by-k lower triangular part
 *          of the array A[j] contains the lower triangular matrix, and the
 *          strictly upper triangular part of A[j] is not referenced.
 *          If diag[i] = BblasUnit, the diagonal elements of A[j] are also not
 *          referenced and are assumed to be 1.
 * 	    batch_count=\sum_{i=1}^{group_count}group_sizes[i].
 *
 * @param[in] lda
 * 	    An arrray of integers of length group_count-1, where
 * 	    lda[i] is the leading dimension of the arrays A[j]
 * 	    in i-th group. lda[i] >= max(1,k[i]).
 *
 * @param[in,out] B
 * 	    B is an array of pointers to matrices B[0], B[1] .. B[batch_count-1], 
 *          On entry, for i-th group each B[j]-s are ldb[i]-by-n[i] right hand side matrix.
 *          On exit, if return value = 0, the ldb[i]-by-n[i] solution matrix X.
 *	    batch_count=\sum_{i=1}^{group_count}group_sizes[i].
 *
 * @param[in] ldb
 * 	    An array of integers of length group_count-1, where
 *	    ldb[i] is the leading dimension of the arrays B[j]
 *	    in i-th group. ldb[i] >= max(1,m[i]).
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
 * @sa ztrsm_batch
 * @sa ctrsm_batch
 * @sa dtrsm_batch
 * @sa strsm_batch
 *
 ******************************************************************************/
void blas_ztrsm_batch(int group_count, const int *group_sizes,
		      bblas_enum_t layout, const bblas_enum_t *side, const bblas_enum_t *uplo,
		      const bblas_enum_t *transa, const bblas_enum_t *diag,
		      const int *m, const int *n, 
		      const bblas_complex64_t *alpha, bblas_complex64_t const *const *A, const int *lda,
		      				      bblas_complex64_t 	    **B, const int *ldb,
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
		if (imin(m[group_iter], n[group_iter]) == 0 || 
				group_sizes[group_iter] == 0) {
			bblas_success(info_option, &info[info_offset], group_sizes[group_iter]);
			continue;
		}

		// Call to blas_ztrsm_batchf 
		blas_ztrsm_batchf(group_sizes[group_iter], 
				  layout, side[group_iter], uplo[group_iter],
				  transa[group_iter], diag[group_iter],
				  m[group_iter], n[group_iter],
				  alpha[group_iter], A+offset, lda[group_iter],
				  		     B+offset, ldb[group_iter],
				  &info[info_offset]);    

		// check for errors in batchf function
		if (info[info_offset] != 0 && flag == 0) {
			info[0] = info[info_offset];	
			flag = 1;
		}

		offset += group_sizes[group_iter];    
	}
}
