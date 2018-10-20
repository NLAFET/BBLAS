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
 * @ingroup gemm_batch
 *
 * blas_zgemm_batch is a batch version of zgemm. It performs 
 * matrix-matrix multiplication
 *
 *  \f[ C[i] = \alpha[i] [op( A[i] )\times op( B[i] )] + \beta[i] C[i], \f]
 *
 *  where op( X ) is one of:
 *    \f[ op( X ) = X,   \f]
 *    \f[ op( X ) = X^T, \f]
 *    \f[ op( X ) = X^H, \f]
 *
 *  where alpha[i] and beta[i] are scalars, and A[i], B[i] and C[i] are matrices, with 
 *  op( A[i] ) an m[i]-by-k[i] matrix, op( B[i] ) a k[i]-by-n[i] matrix and C[i] an 
 *  m[i]-by-n[i] matrix.
 *
 ******************************************************************************
 * @param[in] group_count
 * 	    The number groups of matrices.	  
 *
 * @param[in] group_sizes
 * 	    An array of integers of length group_count, where group_sizes[i] denotes
 * 	    the number of matrices in i-th group.	
 *
 * @param[in] layout
 * 	    Specifies if the matrix is stored in row major or column major
 * 	    format:
 * 	    - BblasRowMajor: Row major format
 * 	    - BblasColMajor: Column major format
 *
 * @param[in] transa
 * 	    An array of size group_count-1, where
 *          - BblasNoTrans:   A[j]-s in i-th group are not transposed,
 *          - BblasTrans:     A[j]-s in i-th group are transposed,
 *          - BblasConjTrans: A[j]-s in i-th group are conjugate transposed.
 *
 * @param[in] transb
 * 	    An array of size group_count, where
 *          - BblasNoTrans:   B[j]-s in the i-th group are not transposed,
 *          - BblasTrans:     B[j]-s in the i-th group are  transposed,
 *          - BblasConjTrans: B[j]-s in the i-th group are conjugate transposed.
 *
 * @param[in] m
 *          An array of integers of size group count, where m[i] is 
 *          the number of rows of matrices op( A[j] ) and of matrices 
 *          C[j] in i-th group. m[i] >= 0.
 *
 * @param[in] n
 *          An array of integers of size group count, where n[i] is 
 *          the number of columns of matrices op( B[j] ) and C[j] in 
 *          i-th group. n[i] >= 0.
 *
 * @param[in] k
 *          An array of integers, where k[i] is the number of columns of  
 *          matrices op( A[j] ) and number of rows of matrices op( B[j] )
 *          in i-th group. k[i] >= 0.
 *
 * @param[in] alpha
 *          An array of scalars of length group-count.
 *
 * @param[in] A
 * 		A is an array of pointers to matrices A[0], A[1] .. A[batch_count-1]. 
 * 		In i-th group each element A[j] is a pointer to a matrix of 
 * 		dimension lda[i]-by-ka[i], where ka[i] is k[i] when transa[i] = 
 * 		BblasNoTrans, and is m[i] otherwise. When using transa[i] = 
 * 		BblasNoTrans the leading m[i]-by-k[i] part of A[j] 
 * 		must contain the matrix elements, otherwise the leading  
 * 		k[i]-by-m[i] part of A[j] must contain the matrix elements.
 * 		batch_count=\sum_{i=1}^{group_count}group_sizes[i].
 *
 * @param[in] lda
 * 	    An array of integers of size group_count, which
 * 	    denotes the leading dimension of the arrays A[j]-s
 * 	    in i-th group. When transa[i] = BblasNoTrans, 
 * 	    lda[i] >= max(1,m[i]), otherwise, lda[i] >= max(1,k[i]).
 *
 * @param[in] B
 * 		B is an array of pointers to matrices B[0], B[1],..,B[batch_count-1].
 * 		In i-th group each element B[j] is a pointer to a matrix
 * 		of dimension lda[i]-by-kb, where kb is n[i] when transb[i] = 
 * 		BblasNoTrans, and is k[i] otherwise. When using transb[i] = 
 * 		BblasNoTrans the leading k[i]-by-n[i] part of B[j] 
 * 		must contain the matrix elements, otherwise the leading  
 * 		n[i]-by-k[i] part of B[j] must contain the matrix elements.
 ** 		batch_count=\sum_{i=1}^{group_count}group_sizes[i].
 *
 * @param[in] ldb
 * 	    An array of integers of size group_count, which 
 * 	    denotes the leading dimension of the array B[j]-s
 * 	    in i-th group. When transb[i] = BblasNoTrans, 
 * 	    ldb[i] >= max(1,k[i]), otherwise, ldb[i] >= max(1,n[i]).
 *
 * @param[in] beta
 *          An array of scalars of length group_count.
 *
 * @param[in,out] C
 * 		C is an array of pointers to matrices C[0], C[1],...,C[batch_count-1].
 * 		In i-th group each element C[j] is a pointer to a matrix of dimension 
 *		ldc[j]-by-n[j]. On exit, each array C[j] of i-th group is overwritten 
 *		by the m[i]-by-n[i] matrix ( alpha[i]*op(A[j] )*op( B[j] ) + beta[i]*C[j] ),
 *		where j=0,1,...,group_sizes[i-1].
 *
 * @param[in] ldc
 * 	    An array of integers of size group_count, which
 *          denotes the leading dimension of the arrays C[j]
 *          in i-th group. ldc[i] >= max(1,m[i]).
 *
 *
 * @param[in,out] info
 * 		Array of int for error handling. On entry info[0] should have one of the 
 * 		following values
 *			- BblasErrorsReportAll    :  All errors will be specified on output.
 *						     Length of the array should be atleast
 *						     \sum_{i=0}^{group_count-1}group_sizes[i].
 *			- BblasErrorsReportGroup  :  Single error from each group will be 
 *						     reported. Length of the array should 
 *						     be atleast group_count.
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
 * @sa zgemm_batch
 * @sa cgemm_batch
 * @sa dgemm_batch
 * @sa sgemm_batch
 *
 ******************************************************************************/
void blas_zgemm_batch( int group_count, const int *group_sizes,
		bblas_enum_t layout, const bblas_enum_t *transa, const bblas_enum_t *transb,
		const int *m, const int *n, const int *k,
		const bblas_complex64_t *alpha, bblas_complex64_t const *const *A, const int *lda,
		                                bblas_complex64_t const* const *B, const int *ldb, 
		const bblas_complex64_t *beta,  bblas_complex64_t** C, const int *ldc, 
		int *info)
{
	// Local variables 
	int group_iter;
	int offset = 0;
	int info_offset = offset;

	// Check input arguments 
	if (group_count < 0) {
		bblas_error("Illegal value of group_count");
		info[0] = 1;
		return;
	}

	// Check group_size and call fixed batch computation 
	for (group_iter = 0; group_iter < group_count; group_iter++) {
		if (group_sizes[group_iter] < 0) {
			bblas_error("Illegal values of group_sizes");
			if (info[0] != BblasErrorsReportNone) {
				bblas_set_info(info[0], &info[0], group_sizes[group_iter], 2);
			}
			return;
		}

		if (group_iter != 0) {
			if (info[0] == BblasErrorsReportAll) 
				info_offset = offset;
			else
				info_offset = group_iter;
		}
		info[info_offset] = info[0];	

		// Call to bblas_zgemm_batchf 
		blas_zgemm_batchf (group_sizes[group_iter], 
				   layout,
				   transa[group_iter],
				   transb[group_iter],
				   m[group_iter],
				   n[group_iter],
				   k[group_iter],
				   alpha[group_iter],
				   A+offset,
				   lda[group_iter],
				   B+offset,
				   ldb[group_iter],
				   beta[group_iter],
				   C+offset,
				   ldc[group_iter],
				   &info[info_offset]);    

		offset += group_sizes[group_iter];    
	}
}
#undef COMPLEX
