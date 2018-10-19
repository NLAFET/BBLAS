/**
 *
 * @file
 *
 *  BBLAS is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 **/

#ifndef BBLAS_TYPES_H
#define BBLAS_TYPES_H

#include <complex.h>

#ifdef __cplusplus
extern "C" {
#endif

/******************************************************************************/
#if defined(HAVE_MKL) || defined(BBLAS_WITH_MKL)
#define lapack_complex_float plasma_complex32_t
#define lapack_complex_double plasma_complex64_t
#endif

/***************************************************************************//**
 *
 *  Some CBLAS routines take scalars by value in real arithmetic
 *  and by pointer in complex arithmetic.
 *  In precision generation, CBLAS_SADDR is removed from real arithmetic files.
 *
 **/
#ifndef CBLAS_SADDR
#define CBLAS_SADDR(var) &(var)
#endif

/*
 *  BBLAS enumerates for parameter values. During 
 *  generating the code for other precisions, Bblas_ConjTrans 
 *  is converted to BblasTrans, where BblasConjTrans is
 *  preserved.
 */


enum {
    BblasInvalid       = -1,

    BblasNoTrans       = 111,
    BblasTrans         = 112,
    BblasConjTrans     = 113,
    Bblas_ConjTrans    = BblasConjTrans,

    BblasUpper         = 121,
    BblasLower         = 122,
    BblasGeneral       = 123,
    BblasGeneralBand   = 124,

    BblasNonUnit       = 131,
    BblasUnit          = 132,

    BblasLeft          = 141,
    BblasRight         = 142,

    BblasOneNorm       = 171,
    BblasRealOneNorm   = 172,
    BblasTwoNorm       = 173,
    BblasFrobeniusNorm = 174,
    BblasInfNorm       = 175,
    BblasRealInfNorm   = 176,
    BblasMaxNorm       = 177,
    BblasRealMaxNorm   = 178,

    BblasForward       = 391,
    BblasBackward      = 392,

    BblasColumnwise    = 401,
    BblasRowwise       = 402,

    BblasW             = 501,
    BblasA2            = 502
};

enum {
    BblasSuccess = 0,
    BblasFail
};

enum {
    BblasErrorsReportAll = 0,
    BblasErrorsReportGroup,
    BblasErrorsReportAny,
    BblasErrorsReportNone 
};

/******************************************************************************/
typedef int bblas_enum_t;

typedef float  _Complex bblas_complex32_t;
typedef double _Complex bblas_complex64_t;

/******************************************************************************/
bblas_enum_t bblas_diag_const(char lapack_char);
bblas_enum_t bblas_direct_const(char lapack_char);
bblas_enum_t bblas_norm_const(char lapack_char);
bblas_enum_t bblas_side_const(char lapack_char);
bblas_enum_t bblas_storev_const(char lapack_char);
bblas_enum_t bblas_trans_const(char lapack_char);
bblas_enum_t bblas_uplo_const(char lapack_char);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif  // BBLAS_TYPES_H 
