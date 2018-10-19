#include "example_z.h"
#include <stdlib.h>
#include <mkl_lapacke.h>

void random_zvec(int len, BBLAS_Complex64_t *vec)
{
//	BBLAS_Complex64_t *work =
//		(BBLAS_Complex64_t*) malloc(sizeof(BBLAS_Complex64_t)*len);
//	int info = 0;
	double scalar = 1.0;
	int seed = 2;
	int zero = 0;
	int one = 1;
	int lenm1 = len-1;
	int colmaj = BblasColMajor;
	LAPACKE_zlagge(colmaj, len, one, lenm1, zero, &scalar, vec, len, &seed);
}

void random_mat(int m, int n, BBLAS_Complex64_t *mat)
{
//	BBLAS_Complex64_t *work = malloc(sizeof(BBLAS_Complex64_t)*m*n);
//	int info = 0;
	double scalar = 1.0;
	int seed = 2;
	int mm1 = m-1;
	int nn1 = n-1;
	int colmaj = BblasColMajor;
	LAPACKE_zlagge(colmaj, m, n, mm1, nn1, &scalar, mat, m, &seed);
}

BBLAS_Complex64_t randz()
{
	BBLAS_Complex64_t val;
	random_mat(1, 1, &val);
	return val;
}

void set_params_fixed_zgemv(struct zgemv_batchf_example *zgemv_example)
{
	int batch_count = rand() % 50;
	zgemv_example->batch_count = batch_count;
	zgemv_example->trans = rand() % 3+111;
	int m = rand() % 200;
	zgemv_example->m = m;
	int n = rand() % 200;
	zgemv_example->n = n;
	zgemv_example->alpha = randz();
	/* arrayA */
	BBLAS_Complex64_t **arrayA =
		(BBLAS_Complex64_t**) malloc(sizeof(BBLAS_Complex64_t*)*batch_count);
	for (int i = 0; i < batch_count; i++)
	{
		arrayA[i] = (BBLAS_Complex64_t*) malloc(sizeof(BBLAS_Complex64_t)*m*n);
		random_mat(m, n, arrayA[i]);
	}
	zgemv_example->arrayA = arrayA;
	zgemv_example->lda = zgemv_example->m;
	/* arrayx */
	BBLAS_Complex64_t **arrayx =
		(BBLAS_Complex64_t**) malloc(sizeof(BBLAS_Complex64_t*)*batch_count);
	int curdim;
	if (zgemv_example->trans == BblasNoTrans)
	{
		curdim = n;
	}
	else
	{
		curdim = m;
	}
	for (int i = 0; i < batch_count; i++)
	{
		arrayx[i] =
			(BBLAS_Complex64_t*) malloc(sizeof(BBLAS_Complex64_t)*curdim);
		random_zvec(curdim, arrayx[i]);
	}
	zgemv_example->arrayx = arrayx;
	zgemv_example->incx = 1;
	zgemv_example->beta = randz();
	/* arrayy */
	BBLAS_Complex64_t **arrayy =
		(BBLAS_Complex64_t**) malloc(sizeof(BBLAS_Complex64_t*)*batch_count);
	if (zgemv_example->trans == BblasNoTrans)
	{
		curdim = m;
	}
	else
	{
		curdim = n;
	}
	for (int i = 0; i < batch_count; i++)
	{
		arrayy[i] =
			(BBLAS_Complex64_t*) malloc(sizeof(BBLAS_Complex64_t)*curdim);
		random_zvec(curdim, arrayy[i]);
	}
	zgemv_example->arrayy = arrayy;
	zgemv_example->incy = 1;
}

void set_params_variable_zgemv(struct zgemv_batchv_example *zgemv_example)
{
	int batch_count = rand() % 50;
	zgemv_example->batch_count = batch_count;
	int *trans = (int*) malloc(sizeof(int)*batch_count);
	for(int i = 0; i < batch_count; i++)
	{
		trans[i] = rand() % 3 + 111;
	}
	zgemv_example->trans = (enum BBLAS_TRANS*) trans;
	int *m = (int*) malloc(sizeof(int)*batch_count);
	for (int i = 0; i < batch_count; i++)
	{
		m[i] = rand() % 200;
	}
	zgemv_example->m = m;
	int *n = (int*) malloc(sizeof(int)*batch_count);
	for (int i = 0; i < batch_count; i++)
	{
		n[i] = rand() % 200;
	}
	zgemv_example->n = n;
	BBLAS_Complex64_t *alpha = malloc(sizeof(BBLAS_Complex64_t)*batch_count);
	for (int i = 0; i < batch_count; i++)
	{
		alpha[i] = randz();
	}
	zgemv_example->alpha = alpha;
	/* arrayA */
	BBLAS_Complex64_t **arrayA =
		(BBLAS_Complex64_t**) malloc(sizeof(BBLAS_Complex64_t*)*batch_count);
	for (int i = 0; i < batch_count; i++)
	{
		arrayA[i] = (BBLAS_Complex64_t*) malloc(sizeof(BBLAS_Complex64_t)*m[i]*n[i]);
		random_mat(m[i], n[i], arrayA[i]);
	}
	zgemv_example->arrayA = arrayA;
	int *lda = (int*) malloc(sizeof(int)*batch_count);
	for (int i = 0; i < batch_count; i++)
	{
		lda[i] = m[i];
	}
	zgemv_example->lda = lda;
	/* arrayx */
	BBLAS_Complex64_t **arrayx =
		(BBLAS_Complex64_t**) malloc(sizeof(BBLAS_Complex64_t*)*batch_count);
	int curdim;
	for (int i = 0; i < batch_count; i++)
	{
		if (trans[i] == BblasNoTrans)
		{
			curdim = n[i];
		}
		else
		{
			curdim = m[i];
		}
		arrayx[i] =
			(BBLAS_Complex64_t*) malloc(sizeof(BBLAS_Complex64_t)*curdim);
		random_zvec(curdim, arrayx[i]);
	}
	zgemv_example->arrayx = arrayx;
	int *incx = (int*) malloc(sizeof(int)*batch_count);
	for (int i = 0; i < batch_count; i++)
	{
		incx[i] = 1;
	}
	zgemv_example->incx = incx;
	BBLAS_Complex64_t *beta = malloc(sizeof(BBLAS_Complex64_t)*batch_count);
	for (int i = 0; i < batch_count; i++)
	{
		beta[i] = randz();
	}
	zgemv_example->beta = beta;
	/* arrayy */
	BBLAS_Complex64_t **arrayy =
		(BBLAS_Complex64_t**) malloc(sizeof(BBLAS_Complex64_t*)*batch_count);
	for (int i = 0; i < batch_count; i++)
	{
		if (trans[i] == BblasNoTrans)
		{
			curdim = m[i];
		}
		else
		{
			curdim = n[i];
		}
		arrayy[i] =
			(BBLAS_Complex64_t*) malloc(sizeof(BBLAS_Complex64_t)*curdim);
		random_zvec(curdim, arrayy[i]);
	}
	zgemv_example->arrayy = arrayy;
	int *incy = (int*) malloc(sizeof(int)*batch_count);
	for (int i = 0; i < batch_count; i++)
	{
		incy[i] = 1;
	}
	zgemv_example->incy = incy;
}

void set_params_fixed_zgemm(struct zgemm_batchf_example *zgemm_example)
{
    int i;
    int batch_count;
    int m;
    int n;
    int k;
    int lda;
    int ldb;
 
	batch_count = rand() % 5;
	zgemm_example->batch_count = batch_count;
	zgemm_example->transA = rand() % 3+111;
	zgemm_example->transB = rand() % 3+111;
	m = (rand()+25) % 20;
	zgemm_example->m = m;
	n = (rand()+25) % 20;
	zgemm_example->n = n;
	k = (rand()+25) % 20;
	zgemm_example->k = k;
	zgemm_example->alpha = randz();
	zgemm_example->beta = randz();


    /* arrayA */
	BBLAS_Complex64_t **arrayA =
		(BBLAS_Complex64_t**) malloc(sizeof(BBLAS_Complex64_t*)*batch_count);
	for (i = 0; i < batch_count; i++)
	{
		arrayA[i] = (BBLAS_Complex64_t*) malloc(sizeof(BBLAS_Complex64_t)*m*k);
		random_mat(m, k, arrayA[i]);
	}
	zgemm_example->arrayA = arrayA;
	if (zgemm_example->transA == BblasNoTrans)
	{
		lda = m;
	}
	else
	{
		lda = k;
	}
	zgemm_example->lda = lda;

    /* arrayB */
	BBLAS_Complex64_t **arrayB =
		(BBLAS_Complex64_t**) malloc(sizeof(BBLAS_Complex64_t*)*batch_count);
	for (i = 0; i < batch_count; i++)
	{
		arrayB[i] = (BBLAS_Complex64_t*) malloc(sizeof(BBLAS_Complex64_t)*k*n);
		random_mat(k, n, arrayB[i]);
	}
	zgemm_example->arrayB = arrayB;
	if (zgemm_example->transB == BblasNoTrans)
	{
		ldb = k;
	}
	else
	{
		ldb = n;
	}
	zgemm_example->ldb = ldb;

    /* arrayC */
	BBLAS_Complex64_t **arrayC =
		(BBLAS_Complex64_t**) malloc(sizeof(BBLAS_Complex64_t*)*batch_count);
	for (i = 0; i < batch_count; i++)
	{
		arrayC[i] = (BBLAS_Complex64_t*) malloc(sizeof(BBLAS_Complex64_t)*m*n);
		random_mat(m, n, arrayC[i]);
	}
	zgemm_example->arrayC = arrayC;
	zgemm_example->ldc = m;
}

// One pointer (op) approach

void set_params_fixed_zgemm_stride(struct zgemm_batchf_example_stride *zgemm_example)
{
    int i;
    int batch_count;
    int m;
    int n;
    int k;
    int lda;
    int ldb;
    enum BBLAS_TRANS transA;
    enum BBLAS_TRANS transB;
 
	batch_count = rand() % 5;
	zgemm_example->batch_count = batch_count;
    transA = BblasNoTrans;
    zgemm_example->transA = transA;
	transB = BblasNoTrans;
    zgemm_example->transB = transB;
	m = rand() % 20;
	zgemm_example->m = m;
	n = rand() % 20;
	zgemm_example->n = n;
	k = rand() % 20;
	zgemm_example->k = k;
	zgemm_example->alpha = randz();
	zgemm_example->beta = randz();

    /* arrayA */
	BBLAS_Complex64_t *arrayA =
		(BBLAS_Complex64_t*) malloc(sizeof(BBLAS_Complex64_t)*batch_count*m*k);
	for (i = 0; i < batch_count; i++)
	{
		random_mat(m, k, &arrayA[i*m*k]);
	}
	zgemm_example->arrayA = arrayA;
	if (zgemm_example->transA == BblasNoTrans)
	{
		lda = m;
	}
	else
	{
		lda = k;
	}
	zgemm_example->lda = lda;
	zgemm_example->strideA = m*k;

    /* arrayB */
	BBLAS_Complex64_t *arrayB =
		(BBLAS_Complex64_t*) malloc(sizeof(BBLAS_Complex64_t*)*batch_count*k*n);
	for (i = 0; i < batch_count; i++)
	{
		random_mat(k, n, &arrayB[i*k*n]);
	}
	zgemm_example->arrayB = arrayB;
	if (zgemm_example->transB == BblasNoTrans)
	{
		ldb = k;
	}
	else
	{
		ldb = n;
	}
	zgemm_example->ldb = ldb;
	zgemm_example->strideB = k*n;

    /* arrayC */
	BBLAS_Complex64_t *arrayC =
		(BBLAS_Complex64_t*) malloc(sizeof(BBLAS_Complex64_t*)*batch_count*m*n);
	for (i = 0; i < batch_count; i++)
	{
		random_mat(m, n, &arrayC[i*m*n]);
	}
	zgemm_example->arrayC = arrayC;
	zgemm_example->ldc = m;
	zgemm_example->strideC = m*n;
}



void set_params_variable_zgemm(struct zgemm_batchv_example *zgemm_example)
{    
    int i;
    int batch_count;
    int lda;
    int ldb;
 
    batch_count = rand() % 5;
    zgemm_example->batch_count = batch_count;
	for (i = 0; i < batch_count; i++)
	{	
	    zgemm_example->transA[i] = rand() % 3+111;
	    zgemm_example->transB[i] = rand() % 3+111;
	    zgemm_example->m[i] = rand() % 20;
	    zgemm_example->n[i] = rand() % 20;
	    zgemm_example->k[i] = rand() % 20;
	    zgemm_example->alpha[i] = randz();
	    zgemm_example->beta[i] = randz();
    }

    /* arrayA */
	BBLAS_Complex64_t **arrayA =
		(BBLAS_Complex64_t**) malloc(sizeof(BBLAS_Complex64_t*)*batch_count);
	for (i = 0; i < batch_count; i++)
	{
		arrayA[i] = (BBLAS_Complex64_t*) malloc(sizeof(BBLAS_Complex64_t)*
                zgemm_example->m[i]*zgemm_example->k[i]);
		random_mat(zgemm_example->m[i], zgemm_example->k[i], arrayA[i]);
	    if (zgemm_example->transA[i] == BblasNoTrans)
	    {
		    lda = zgemm_example->m[i];
	    }
	    else
	    {
		    lda = zgemm_example->k[i];
	    }
	    zgemm_example->lda[i] = lda;
	}
	zgemm_example->arrayA = arrayA;

    /* arrayB */
	BBLAS_Complex64_t **arrayB =
		(BBLAS_Complex64_t**) malloc(sizeof(BBLAS_Complex64_t*)*batch_count);
	for (i = 0; i < batch_count; i++)
	{
		arrayB[i] = (BBLAS_Complex64_t*) malloc(sizeof(BBLAS_Complex64_t)*
                zgemm_example->k[i]*zgemm_example->n[i]);
		random_mat(zgemm_example->k[i], zgemm_example->n[i], arrayB[i]);
        if (zgemm_example->transB[i] == BblasNoTrans)
	    {
		    ldb = zgemm_example->k[i];
	    }
	    else
	    {
		    ldb = zgemm_example->n[i];
	    }
	    zgemm_example->ldb[i] = ldb;
	}
	zgemm_example->arrayB = arrayB;
	
    /* arrayC */
	BBLAS_Complex64_t **arrayC =
		(BBLAS_Complex64_t**) malloc(sizeof(BBLAS_Complex64_t*)*batch_count);
	for (i = 0; i < batch_count; i++)
	{
		arrayC[i] = (BBLAS_Complex64_t*) malloc(sizeof(BBLAS_Complex64_t)*
                zgemm_example->m[i]*zgemm_example->n[i]);
		random_mat(zgemm_example->m[i], zgemm_example->n[i], arrayC[i]);
	    zgemm_example->ldc[i] = zgemm_example->m[i];
	}
	zgemm_example->arrayC = arrayC;

}
