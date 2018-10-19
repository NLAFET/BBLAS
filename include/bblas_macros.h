/**
 * @file bblas_macros.h
 *
 * @brief BBLAS macro definitions.
 *
 *  BBLAS is a software package provided by 
 *  Univ. of Manschester,
 *  Univ. of Tennessee.
 *
 * @version 1.0.0
 * @author Samuel  D. Relton
 * @author Pedro   V. Lara
 * @author Mawussi Zounon
 * @date 2018-09-21
 *
 * Contains macros for success and error codes, max, min, and passing complex values to CBLAS.
 *
 **/
#ifndef BBLAS_MACROS_H
#define BBLAS_MACROS_H


/*
 * BBLAS Return Codes
 */
#define BblasSuccess                0 //!< Computation completed successfully
/* Batch error codes */
#define BblasErrorBatchCount       -1 //!< Error in batch_count
/* Subproblem error codes */
#define BblasErrorM                 -3 //!< Error in M
#define BblasErrorN                 -4 //!< Error in N
#define BblasErrorK                 -5 //!< Error in K
#define BblasErrorlda               -6 //!< Error in lda
#define BblasErrorldb               -7 //!< Error in ldb
#define BblasErrorldc               -8 //!< Error in ldc
#define BblasErrorUplo              -9 //!< Error in uplo
#define BblasErrorTransa            -10 //!< Error in transA
#define BblasErrorTransb            -11 //!< Error in transB
#define BblasErrorTrans             -12 //!< Error in trans
#define BblasErrorSide              -13 //!< Error in side
#define BblasErrorDiag              -14 //!< Error in diag
#define BblasErrorGroupCount        -15 //!< Error in group count
#define BblasErrorGroupSize         -16 //!< Error in group size

#define NAME_LENGTH 30 //!< Maximum length of routine name e.g. zgemm_batch

#ifndef imax
#define imax(a, b) ((a) > (b) ? (a) : (b)) //!< Take max of two numbers
#endif

#ifndef imin
#define imin(a, b) ((a) < (b) ? (a) : (b)) //!< Take min of two numbers
#endif

/**
 * CBLAS requires for complex scalar arguments to be passed by address rather than by value
 **/

#endif /* BBLAS_MACROS_H */
