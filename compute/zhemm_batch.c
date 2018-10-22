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
 * @ingroup hemm_batch
 *
 *  Performs one of the batch matrix-matrix operations on each group of matrices
 *
 *     \f[ C[i] = \alpha[i] \times A[i] \times B[i] + \beta[i] \times C[i] \f]
 *  or
 *     \f[ C[i] = \alpha[i] \times B[i] \times A[i] + \beta[i] \times C[i] \f]
 *
 *  where alpha[i] and beta[i] are scalars, A[i]-s are Hermitian matrices B[i]-s and
 *  C[i]-s are m-by-n matrices.
 *
 *******************************************************************************
 * @param[in] group_count
 * 	    The number groups of matrices.	  
 *
 * @param[in] group_sizes
 * 	    An array of integers of length group_count, where 
 * 	    group_sizes[i] denotes the number of matrices in i-th group.	
 *
 * @param[in] layout
 * 	    Specifies if the matrix is stored in row major or column major
 * 	    format:
 * 	    - BblasRowMajor: Row major format
 * 	    - BblasColMajor: Column major format
 *
 * @param[in] side
 * 	    An array of length group_count. side[i]
 *          Specifies whether the Hermitian matrices A[j]-s of i-th 
 *          group appear on the left or right in the operation as follows:
 *          - BblasLeft: \f[ C[j] = \alpha[i] \times A[j] \times B[j] + \beta[i] \times C[j] \f] 
 *          - BblasRight: \f[ C[j] = \alpha[i] \times B[j] \times A[j] + \beta[i] \times C[j] \f]
 *
 * @param[in] uplo
 * 	    An array of length group_count, where uplo[i]
 *          specifies whether the upper or lower triangular part of
 *          the Hermitian matrices A[j]-s of i-th group are to be referenced 
 *          
 *          - BblasLower:     Only the lower triangular part of the
 *                             Hermitian matrices A[j] is to be referenced.
 *          - BblasUpper:     Only the upper triangular part of the
 *                             Hermitian matrices A[j] is to be referenced.
 *
 * @param[in] m
 * 	    An of array of integers of length group_count,
 *          where m[i] denotes the number of rows in the matrices C[j]
 *          of i-th group. m[i] >= 0.
 *
 * @param[in] n
 * 	    An of array of integers of length group_count,
 *          where n[i] denotes the number of columns in the matrices C[j]
 *          of i-th group. n[i] >= 0.
 *
 * @param[in] alpha
 *          An array of scalars of length group-count.
 *
 * @param[in] A
 * 		A is an array of pointers to matrices A[0], A[1] .. A[batch_count-1], 
 * 		where each element A[j] of i-th group is a pointer to a 
 * 		matrix A[j] of size lda[i]-by-ka, where ka is m[i] when 
 * 		side[i] = BblasLeft, and is n[i] otherwise. Only the uplo 
 * 		triangular part is referenced.
 *              batch_count=\sum_{i=1}^{group_count}group_sizes[i].
 *
 *
 * @param[in] lda
 * 	    An array of length group_count, where lda[i] is
 *          the leading dimension of the arrays A[j] of i-th group. 
 *          lda[i] >= max(1,ka).
 *
 * @param[in] B
 * 	    B is an array of pointers to matrices B[0], B[1] .. B[batch_count-1], 
 * 	    where each element B[j] of i-th group is a pointer to a 
 * 	    matrix B[j] of size ldb[i]-by-n[i] matrix, where the leading 
 * 	    m[i]-by-n[i] part of the array B[j] must contain the matrix B[j].
 *          batch_count=\sum_{i=1}^{group_count}group_sizes[i].
 *
 * @param[in] ldb
 * 	    An array of length group_count, where ldb[i] is
 *          the leading dimension of the arrays B[j] of i-th
 *          group. ldb[i] >= max(1,m[i]).
 *
 * @param[in] beta
 *          An array of scalars of length group_count.
 *
 * @param[in,out] C
 * 		C is an array of pointers to matrices C[0], C[1],...,C[batch_count-1].
 * 		In i-th group each element C[j] is a pointer to a matrix C[j].
 *		On exit, each array C[j] of i-th group is overwritten 
 *              by the m[i]-by-n[i] updated matrix.
 *              batch_count=\sum_{i=1}^{group_count}group_sizes[i].
 *
 * @param[in] ldc
 * 	    An array of integers of size group_count, which
 *          denotes the leading dimension of the arrays C[j]
 *          in i-th group. ldc[i] >= max(1,m[i]).
 *
 * @param[in,out] info
 * 		Array of int for error handling. On entry info[0] should have one of the 
 * 		following values
 *			- BblasErrorsReportAll    :  All errors will be specified on output.
 *						     Length of the array should be atleast
 *						     (group_count*batch_count).
 *			- BblasErrorsReportGroup  :  Single error from each group will be 
 *						     reported. Length of the array should 
 *						     be atleast to (group_count).
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
 * @sa zhemm_batch
 * @sa chemm_batch
 *
 ******************************************************************************/
void blas_zhemm_batch(int group_count, const int *group_sizes,
		      bblas_enum_t layout, const bblas_enum_t *side, const bblas_enum_t *uplo,
		      const int *m, const int *n, 
		      const bblas_complex64_t *alpha, bblas_complex64_t const *const *A, const int *lda, 
		      				      bblas_complex64_t const* const *B, const int *ldb, 
		      const bblas_complex64_t *beta,  bblas_complex64_t		   ** C, const int *ldc, 
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
		if (m[group_iter] == 0 || n[group_iter]== 0 ||
				((alpha[group_iter] == (bblas_complex64_t)0.0) &&
				 (beta[group_iter] == (bblas_complex64_t)1.0)) || 
				group_sizes[group_iter] == 0) {
			bblas_success(info_option, &info[info_offset], group_sizes[group_iter]);
			continue;
		}

		// Call to blas_zhemm_batchf 
		blas_zhemm_batchf(group_sizes[group_iter], 
				  layout, side[group_iter], uplo[group_iter],
				  m[group_iter], n[group_iter],
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

