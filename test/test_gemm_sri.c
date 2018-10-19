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

	/* Create two fixed size batches to use for fixed/grouped API */
	struct zgemm_batchf_example *zgemmf_example =
		(struct zgemm_batchf_example*) malloc(sizeof(struct zgemm_batchf_example));
	set_params_fixed_zgemm(zgemmf_example);


	struct zgemm_batchf_example *zgemmf_example2 =
		(struct zgemm_batchf_example*) malloc(sizeof(struct zgemm_batchf_example));
	set_params_fixed_zgemm(zgemmf_example2);

    /* Create one fixed batch for the one pointer (stride) approach */
    struct zgemm_batchf_example_stride *zgemmf_example_stride =
		(struct zgemm_batchf_example_stride*) malloc(sizeof(struct zgemm_batchf_example_stride));
	set_params_fixed_zgemm_stride(zgemmf_example_stride);

	/* Creat one variable size batch for variable API */
	struct zgemm_batchv_example *zgemmv_example =
		(struct zgemm_batchv_example*) malloc(sizeof(struct zgemm_batchv_example));
    set_params_variable_zgemm(zgemmv_example);

	printf("Finished creating inputs, running each API...\n");

	int *infof = malloc(sizeof(int)*1);
	int *infof_stride = malloc(sizeof(int)*1);
	int *infov = malloc(sizeof(int)*zgemmv_example->batch_count);

	printf("Computing with first API...\n");
	/*
	   First API: Use separate calls for fixed/variable

        batchf_zgemm(
        const enum BBLAS_TRANS transA, const enum BBLAS_TRANS transB,
        const int M,  const int N, const int K,
        const BBLAS_Complex64_t alpha,
        const BBLAS_Complex64_t **arrayA, const int lda,
        const BBLAS_Complex64_t **arrayB, const int ldb, 
        const BBLAS_Complex64_t beta,
        BBLAS_Complex64_t **arrayC, const int ldc, 
        const int batch_count, int info);
	*/
	
    batchf_zgemm(
		(const enum BBLAS_TRANS) zgemmf_example->transA,
		(const enum BBLAS_TRANS) zgemmf_example->transB,
		(const int) zgemmf_example->m,
		(const int) zgemmf_example->n,
		(const int) zgemmf_example->k,
		(const BBLAS_Complex64_t) zgemmf_example->alpha,
		(const BBLAS_Complex64_t**) zgemmf_example->arrayA,
		(const int) zgemmf_example->lda,
		(const BBLAS_Complex64_t **) zgemmf_example->arrayB,
		(const int) zgemmf_example->ldb,
		(const BBLAS_Complex64_t) zgemmf_example->beta,
		(BBLAS_Complex64_t **) zgemmf_example->arrayC,
		(const int) zgemmf_example->ldc,
		(const int) zgemmf_example->batch_count, 
        infof);
    /*
	    Use one pointer approach for fixed size

        batchf_zgemm_stride(
        const enum BBLAS_TRANS transA, const enum BBLAS_TRANS transB,
        const int M,  const int N, const int K,
        const BBLAS_Complex64_t alpha,
        const BBLAS_Complex64_t *arrayA, const int lda, const int strideA,
        const BBLAS_Complex64_t *arrayB, const int ldb, const int strideB, 
        const BBLAS_Complex64_t beta,
        BBLAS_Complex64_t *arrayC, const int ldc, const int strideC, 
        const int batch_count, int info);
	*/
	
    batchf_zgemm_stride(
		(const enum BBLAS_TRANS) zgemmf_example_stride->transA,
		(const enum BBLAS_TRANS) zgemmf_example_stride->transB,
		(const int) zgemmf_example_stride->m,
		(const int) zgemmf_example_stride->n,
		(const int) zgemmf_example_stride->k,
		(const BBLAS_Complex64_t) zgemmf_example_stride->alpha,
		(const BBLAS_Complex64_t *) zgemmf_example_stride->arrayA,
		(const int) zgemmf_example_stride->lda,
		(const int) zgemmf_example_stride->strideA,
		(const BBLAS_Complex64_t *) zgemmf_example_stride->arrayB,
		(const int) zgemmf_example_stride->ldb,
		(const int) zgemmf_example_stride->strideB,
		(const BBLAS_Complex64_t) zgemmf_example_stride->beta,
		(BBLAS_Complex64_t *) zgemmf_example_stride->arrayC,
		(const int) zgemmf_example_stride->ldc,
		(const int) zgemmf_example_stride->strideC,
		(const int) zgemmf_example_stride->batch_count, 
        infof_stride);

	/*
        batchv_zgemm(
        const enum BBLAS_TRANS *transA, const enum BBLAS_TRANS *transB,
        const int *M,  const int *N, const int *K,
        const BBLAS_Complex64_t *alpha,
        const BBLAS_Complex64_t **arrayA, const int *lda,
        const BBLAS_Complex64_t **arrayB, const int *ldb, 
        const BBLAS_Complex64_t *beta,
        BBLAS_Complex64_t **arrayC, const int *ldc, 
        const int batch_count, int *info);
	*/
	batchv_zgemm(
		(const enum BBLAS_TRANS*) zgemmv_example->transA,
		(const enum BBLAS_TRANS*) zgemmv_example->transB,
		(const int*) zgemmv_example->m,
		(const int*) zgemmv_example->n,
		(const int*) zgemmv_example->k,
		(const BBLAS_Complex64_t*) zgemmv_example->alpha,
		(const BBLAS_Complex64_t**) zgemmv_example->arrayA,
		(const int*) zgemmv_example->lda,
		(const BBLAS_Complex64_t **) zgemmv_example->arrayB,
		(const int*) zgemmv_example->ldb,
		(const BBLAS_Complex64_t*) zgemmv_example->beta,
		(BBLAS_Complex64_t **) zgemmv_example->arrayC,
		(const int*) zgemmv_example->ldc,
		(const int) zgemmv_example->batch_count, 
        infov);

	printf("Computing with second API...\n");
	
    /* Group two fixed batched together */
	int group_count = 2;
	int infog[2];
	int group_size[2] = {zgemmf_example->batch_count, zgemmf_example2->batch_count};
	int group_transA[2] = {zgemmf_example->transA, zgemmf_example2->transA};
	int group_transB[2] = {zgemmf_example->transB, zgemmf_example2->transB};
	int group_m[2] = {zgemmf_example->m, zgemmf_example2->m};
	int group_n[2] = {zgemmf_example->n, zgemmf_example2->n};
	int group_k[2] = {zgemmf_example->k, zgemmf_example2->k};
	BBLAS_Complex64_t group_alpha[2] = {zgemmf_example->alpha, zgemmf_example2->alpha};
	BBLAS_Complex64_t group_beta[2] = {zgemmf_example->beta, zgemmf_example2->beta};
	int group_lda[2] = {zgemmf_example->lda, zgemmf_example2->lda};
	int group_ldb[2] = {zgemmf_example->ldb, zgemmf_example2->ldb};
	int group_ldc[2] = {zgemmf_example->ldc, zgemmf_example2->ldc};
	
    /* ArrayA */
    BBLAS_Complex64_t* group_arrayA[zgemmf_example->batch_count + zgemmf_example2->batch_count];
	for (i = 0; i < zgemmf_example->batch_count; i++)
	{
		group_arrayA[i] = zgemmf_example->arrayA[i];
	}
	for (i = 0; i < zgemmf_example2->batch_count; i++)
	{
		group_arrayA[i+zgemmf_example->batch_count] = zgemmf_example2->arrayA[i];
	}
    
    /* ArrayB */
	BBLAS_Complex64_t* group_arrayB[zgemmf_example->batch_count + zgemmf_example2->batch_count];
	for (i = 0; i < zgemmf_example->batch_count; i++)
	{
		group_arrayB[i] = zgemmf_example->arrayB[i];
	}
	for (i = 0; i < zgemmf_example2->batch_count; i++)
	{
		group_arrayB[i+zgemmf_example->batch_count] = zgemmf_example2->arrayB[i];
	}
    
    /* ArrayC */
	BBLAS_Complex64_t* group_arrayC[zgemmf_example->batch_count + zgemmf_example2->batch_count];
	for (i = 0; i < zgemmf_example->batch_count; i++)
	{
		group_arrayC[i] = zgemmf_example->arrayC[i];
	}
	for (i = 0; i < zgemmf_example2->batch_count; i++)
	{
		group_arrayC[i+zgemmf_example->batch_count] = zgemmf_example2->arrayC[i];
	}

	/*
        batchg_zgemm(
        const enum BBLAS_TRANS *transA, const enum BBLAS_TRANS *transB,
        const int *M,  const int *N, const int *K,
        const BBLAS_Complex64_t *alpha,
        const BBLAS_Complex64_t **arrayA, const int *lda,
        const BBLAS_Complex64_t **arrayB, const int *ldb, 
        const BBLAS_Complex64_t *beta,
        BBLAS_Complex64_t **arrayC, const int *ldc, 
	    const int group_count, const int *group_size,
        int *info);
	 */
	batchg_zgemm(
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
		infog);

	printf("Computing with third API...\n");
	// Fixed batch
	enum BBLAS_OPTS batch_opts = BBLAS_FIXED;
	batch_zgemm(
		(const enum BBLAS_TRANS*) &zgemmf_example->transA,
		(const enum BBLAS_TRANS*) &zgemmf_example->transB,
		(const int*) &zgemmf_example->m, 
        (const int*) &zgemmf_example->n,
        (const int*) &zgemmf_example->k,
		(const BBLAS_Complex64_t*) &zgemmf_example->alpha,
		(const BBLAS_Complex64_t**) zgemmf_example->arrayA,
		(const int*) &zgemmf_example->lda,
		(const BBLAS_Complex64_t**) zgemmf_example->arrayB,
		(const int*) &zgemmf_example->ldc,
		(const BBLAS_Complex64_t*) &zgemmf_example->beta,
		(BBLAS_Complex64_t**) zgemmf_example->arrayC,
		(const int*) &zgemmf_example->ldc,
		(const int) zgemmf_example->batch_count, 
        (const enum BBLAS_OPTS) batch_opts,
		infof);

	batch_opts = BBLAS_VARIABLE;

	batch_zgemm(
		(const enum BBLAS_TRANS*) zgemmv_example->transA,
		(const enum BBLAS_TRANS*) zgemmv_example->transA,
		(const int*) zgemmv_example->m, 
        (const int*) zgemmv_example->n,
        (const int*) zgemmv_example->k,
		(const BBLAS_Complex64_t*) zgemmv_example->alpha,
		(const BBLAS_Complex64_t**) zgemmv_example->arrayA,
		(const int*) zgemmv_example->lda,
		(const BBLAS_Complex64_t **) zgemmv_example->arrayB,
		(const int*) zgemmv_example->ldb,
		(const BBLAS_Complex64_t*) zgemmv_example->beta,
		(BBLAS_Complex64_t **) zgemmv_example->arrayC,
		(const int*) zgemmv_example->ldc,
		(const int) zgemmv_example->batch_count, 
        (const enum BBLAS_OPTS) batch_opts,
		infov);

	printf("Completed computations with all APIs...\n");

    return 0;
}
