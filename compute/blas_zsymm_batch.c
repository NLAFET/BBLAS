/**
 * @file blas_zsymm_batch.c
 *
 * @brief BBLAS zsymm_batch double _Complex routine.
 *
 *  BBLAS is a software package provided by 
 *  Univ. of Manchester,
 *  Univ. of Tennessee.
 *
 * @author  Srikara Pranesh
 * @author  Mawussi Zounon
 * @date    2018-10-08
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
 * @ingroup symm_batch
 *
 *  Performs one of the batch matrix-matrix operations 
 *
 *     \f[ C[i] = \alpha \times A[i] \times B[i] + \beta \times C[i] \f]
 *  or
 *     \f[ C[i] = \alpha \times B[i] \times A[i] + \beta \times C[i] \f]
 *
 *
 *  for a group of matrices, where alpha[i] and beta[i] are scalars, A[i]-s are 
 *  symmetric matrices and B[i]-s are C[i] are m[i]-by-n[i] matrices.
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
 * 	    An array of size group_count. side[i]
 *          Specifies whether the Hermitian matrices A[j]-s of i-th 
 *          group appear on the left or right in the operation as follows:
 *          - BblasLeft: \f[ C[j] = \alpha[i] \times A[j] \times B[j] + \beta[i] \times C[j] \f] 
 *          - BblasRight: \f[ C[j] = \alpha[i] \times B[j] \times A[j] + \beta[i] \times C[j] \f]
 *
 * @param[in] uplo
 * 	    An array of length group_count, where uplo[i]
 *          specifies whether the upper or lower triangular part of
 *          the symmetric matrices A[j]-s of i-th group are to be referenced 
 *          
 *          - BblasLower:     Only the lower triangular part of the
 *                            symmetric matrices A[j] is to be referenced.
 *          - BblasUpper:     Only the upper triangular part of the
 *                            symmetric matrices A[j] is to be referenced.
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

void blas_zsymm_batch( int group_count, const int *group_sizes,
		        bblas_enum_t layout, const bblas_enum_t *side, const bblas_enum_t *uplo,
    			const int *m, const int *n, 
			const bblas_complex64_t *alpha, bblas_complex64_t const *const *A, const int *lda, 
    			 		 	 bblas_complex64_t const* const *B, const int *ldb, 
			const bblas_complex64_t *beta,  bblas_complex64_t** C, const int *ldc, 
    			 int *info)


{
	// Local variables 
	int group_iter;
	int offset = 0;
	char func_name[15] = "batch_zsymm";

	// Check input arguments 
	if (group_count < 0) {
		xerbla_batch(func_name, BblasErrorGroupCount, -1);
		return;
	}

	// Check group_size and call fixed batch computation 
	for (group_iter = 0; group_iter < group_count; group_iter++) {
		if (group_sizes[group_iter] < 0) {
			xerbla_batch(func_name, BblasErrorGroupSize, group_iter);
			info[group_iter] = BblasErrorGroupSize;
			continue;
		}

		// Call to bblas_zsymm_batchf 
		blas_zsymm_batchf (group_sizes[group_iter], 
				layout,
				side[group_iter],
				uplo[group_iter],
				m[group_iter],
				n[group_iter],
				alpha[group_iter],
				A+offset,
				lda[group_iter],
				B+offset,
				ldb[group_iter],
				beta[group_iter],
				C+offset,
				ldc[group_iter],
				&info[group_iter]);    

		offset += group_sizes[group_iter];    
	}
}
#undef COMPLEX