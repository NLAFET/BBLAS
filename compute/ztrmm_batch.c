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
 * @ingroup trmm_batch
 *
 *  Performs a triangular batch matrix-matrix multiply of the form 
 *
 *          \f[B[i] = \alpha [op( A[i] ) \times B[i]] \f], if side = BblasLeft  or
 *          \f[B[i] = \alpha [B[i] \times op( A[i]) ] \f], if side = BblasRight
 *
 *  for a group of matrices, and op( X ) is one of:
 *
 *          - op(A[i]) = A[i]   or
 *          - op(A[i]) = A[i]^T or
 *          - op(A[i]) = A[i]^H
 *
 *  alpha[i]-s are scalars, B[i]-s are m-by-n matrices and A[i]-s are a unit or non-unit, upper
 *  or lower triangular matrix.
 *
 *******************************************************************************
 * @param[in] group_count
 * 	    The number groups of matrices.	  
 *
 * @param[in] group_sizes
 * 	    The number of matrices in each group.	
 * 
 * @param[in] layout
 * 	    Specifies if the matrix is stored in row major or column major
 * 	    format:
 * 	    - BblasRowMajor: Row major format
 * 	    - BblasColMajor: Column major format
 *
 * @param[in] side
 * 	    An array of length group_count-1, for matrices of i-th group it
 *          specifies whether op( A[j] ) appears on the left or on the right of B[j]:
 *          - BblasLeft:  alpha[i]*op( A[j] )*B[j]
 *          - BblasRight: alpha[i]*B[j]*op( A[j] )
 *
 * @param[in] uplo
 * 	    An array of length group_count-1, where uplo[i]
 *          specifies whether the upper or lower triangular part of
 *          the symmetric matrices A[j]-s of i-th group are to be referenced 
 *          
 *          - BblasLower:     Only the lower triangular part of the
 *                            symmetric matrices A[j] is to be referenced.
 *          - BblasUpper:     Only the upper triangular part of the
 *                            symmetric matrices A[j] is to be referenced.
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
 * 		A is an array of pointers to matrices A[0], A[1] .. A[batch_count-1], 
 * 		where for i-th group each element A[j] is a pointer to a triangular matrix 
 *          	of dimension lda[i]-by-k, where k is m[i] when
 *          	side='L' or 'l' and k is n[i] when when side='R' or 'r'. If uplo =
 *          	BblasUpper, the leading k-by-k upper triangular part of the array
 *          	A[j] contains the upper triangular matrix, and the strictly lower
 *          	triangular part of A[j] is not referenced. If uplo = BblasLower, the
 *          	leading k-by-k lower triangular part of the array A[j] contains the
 *          	lower triangular matrix, and the strictly upper triangular part of
 *          	A[j] is not referenced. If diag = BblasUnit, the diagonal elements of
 *          	A[j] are also not referenced and are assumed to be 1.
 * 		batch_count=\sum_{i=1}^{group_count}group_sizes[i].
 *
 *
 * @param[in] lda
 * 	    An array of integers of length group_count-1, where lda[i]
 *          denotes the leading dimension of the arrays A[j] of i-th group. 
 *          When side='L' or 'l', lda[i] >= max(1,m[i]), 
 *          when side='R' or 'r' then lda[i] >= max(1,n[i]).
 *
 * @param[in] B
 * 		B is an array of pointers to matrices B[0], B[1],..,B[batch_count-1],
 * 		where for i-th group each element B[j] is a pointer to a matrix. 
 *          	On entry, the matrices B[j] are of dimension ldb[i]-by-n[i].
 *          	On exit, the result of a triangular matrix-matrix multiply
 *          	( alpha[i]*op(A[j])*B[j] ) or ( alpha[i]*B[j]*op(A[j]) ).
 * 		batch_count=\sum_{i=1}^{group_count}group_sizes[i].
 *
 *
 * @param[in] ldb
 * 	    An array of integers of length group_count-1, where ldb[i]
 *          is the leading dimension of the arrays B[j] of i-th group. 
 *          ldb[i] >= max(1,m[i]).
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
 * @sa ztrmm_batch
 * @sa ctrmm_batch
 * @sa dtrmm_batch
 * @sa strmm_batch
 *
 ******************************************************************************/
void blas_ztrmm_batch(int group_count, const int *group_sizes,
		      bblas_enum_t layout, const bblas_enum_t *side, const bblas_enum_t *uplo,
		      const bblas_enum_t *transa, const bblas_enum_t *diag,
		      const int *m, const int *n, 
		      const bblas_complex64_t *alpha, bblas_complex64_t const *const *A, const int *lda,
		      				      bblas_complex64_t 	    **B, int const *ldb,
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

		// Call to blas_zher2k_batchf 
		blas_ztrmm_batchf(group_sizes[group_iter], 
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

