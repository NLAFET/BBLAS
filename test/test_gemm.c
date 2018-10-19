/**
 *
 * @file test_gemm.c
 *
 * @author Samuel  D. Relton
 * @author Pedro   V. Lara
 * @author Mawussi Zounon
 * @date 2016-06-15
 * @precisions normal z -> c
 *
 **/
#include "bblas_z.h"
#include "example_z.h"
#include <stdlib.h>
#include <stdio.h>

int main(int argc, char *argv[])
{
    int i;

    printf("Testing gemm using some very simple tests.\n");

	/* Create simple test matrices */
	printf("Creating inputs for testing grouped zgemm...\n");

	/* Create two fixed size batches to use for grouped API */
	struct zgemm_batchf_example *zgemmf_example =
		(struct zgemm_batchf_example*) malloc(sizeof(struct zgemm_batchf_example));
	set_params_fixed_zgemm(zgemmf_example);

	/*struct zgemm_batchf_example *zgemmf_example2 =
		(struct zgemm_batchf_example*) malloc(sizeof(struct zgemm_batchf_example));
	set_params_fixed_zgemm(zgemmf_example2);*/


	printf("Computing with the group API\n");
	
    /* Group two fixed batched together */
	int group_count = 2;
	int infog[2];
	int group_size[2] = {zgemmf_example->batch_count, zgemmf_example->batch_count};
	int group_transA[2] = {zgemmf_example->transA, zgemmf_example->transA};
	int group_transB[2] = {zgemmf_example->transB, zgemmf_example->transB};
	int group_m[2] = {zgemmf_example->m, zgemmf_example->m};
	int group_n[2] = {zgemmf_example->n, zgemmf_example->n};
	int group_k[2] = {zgemmf_example->k, zgemmf_example->k};
	bblas_complex64_t group_alpha[2] = {zgemmf_example->alpha, zgemmf_example->alpha};
	bblas_complex64_t group_beta[2] = {zgemmf_example->beta, zgemmf_example->beta};
	int group_lda[2] = {zgemmf_example->lda, zgemmf_example->lda};
	int group_ldb[2] = {zgemmf_example->ldb, zgemmf_example->ldb};
	int group_ldc[2] = {zgemmf_example->ldc, zgemmf_example->ldc};
	
    /* ArrayA */
    bblas_complex64_t* group_arrayA[zgemmf_example->batch_count + zgemmf_example->batch_count];
	for (i = 0; i < zgemmf_example->batch_count; i++)
	{
		group_arrayA[i] = zgemmf_example->arrayA[i];
	}
	for (i = 0; i < zgemmf_example->batch_count; i++)
	{
		group_arrayA[i+zgemmf_example->batch_count] = zgemmf_example->arrayA[i];
	}
    
    /* ArrayB */
	bblas_complex64_t* group_arrayB[zgemmf_example->batch_count + zgemmf_example->batch_count];
	for (i = 0; i < zgemmf_example->batch_count; i++)
	{
		group_arrayB[i] = zgemmf_example->arrayB[i];
	}
	for (i = 0; i < zgemmf_example->batch_count; i++)
	{
		group_arrayB[i+zgemmf_example->batch_count] = zgemmf_example->arrayB[i];
	}
    
    /* ArrayC */
	bblas_complex64_t* group_arrayC[zgemmf_example->batch_count + zgemmf_example->batch_count];
	for (i = 0; i < zgemmf_example->batch_count; i++)
	{
		group_arrayC[i] = zgemmf_example->arrayC[i];
	}
	for (i = 0; i < zgemmf_example->batch_count; i++)
	{
		group_arrayC[i+zgemmf_example->batch_count] = zgemmf_example->arrayC[i];
	}


void blas_zgemm_batch( (int) group_count, (const int*) group_sizes,
		(bblas_enum_t) layout, (const bblas_enum_t*) group_transA, (const bblas_enum_t*) group_transB,
		(const int*) group_m, (const int*) group_n, (const int*) group_k,
		(const bblas_complex64_t*) group_alpha, (const bblas_complex64_t**) const group_arrayA, (const int*) lda,
		(const bblas_complex64_t**) group_arrayB, (const int*) group_ldb, 
		(const bblas_complex64_t*) group_beta, group_arrayC, (const int*) group_ldc, 
		info);


/*	batchg_zgemm(
		(const enum BBLAS_TRANS*) group_transA,
		(const enum BBLAS_TRANS*) group_transB,
		(const int*) group_m, 
        (const int*) group_n,
        (const int*) group_k,
		(const BBLAS_Complex64_t*) group_alpha,
		(const BBLAS_Complex64_t**) group_arrayA, (const int*) group_lda,
		(const BBLAS_Complex64_t**) group_arrayB, (const int*) group_ldb,
		(const BBLAS_Complex64_t*) group_beta,
		group_arrayC, (const int*) group_ldc,
		(const int) group_count, 
        (const int*) group_size,
		infog);*/


	printf("Completed computations with the group API\n");

    return 0;
}
