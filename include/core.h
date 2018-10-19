/**
 *
 * @file
 *
 *  BBLAS is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 **/
#ifndef CORE_BBLAS_H
#define CORE_BBLAS_H

#include <stdio.h>
#include "bblas_error.h"
#include "bblas_types"

#ifdef __cplusplus
extern "C" {
#endif

/******************************************************************************/
static inline int imin(int a, int b)
{
    if (a < b)
        return a;
    else
        return b;
}

/******************************************************************************/
static inline int imax(int a, int b)
{
    if (a > b)
        return a;
    else
        return b;
}

    
static const char *lapack_constants[] = {
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",

    "", "", "", "", "", "", "", "", "", "",
    "",
    "NoTrans",                              ///< 111: BblasNoTrans
    "Trans",                                ///< 112: BblasTrans
    "ConjTrans",                            ///< 113: BblasConjTrans

    "", "", "", "", "", "", "",
    "Upper",                                ///< 121: BblasUpper
    "Lower",                                ///< 122: BblasLower
    "General",                              ///< 123: BblasGeneral

    "", "", "", "", "", "", "",
    "NonUnit",                              ///< 131: BblasNonUnit
    "Unit",                                 ///< 132: BblasUnit

    "", "", "", "", "", "", "", "",
    "Left",                                 ///< 141: BblasLeft
    "Right",                                ///< 142: BblasRight

    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "",
    "One",                                  ///< 171: BblasOneNorm
    "",                                     ///< 172: BblasRealOneNorm
    "Two",                                  ///< 173: BblasTwoNorm
    "Frobenius",                            ///< 174: BblasFrobeniusNorm
    "Infinity",                             ///< 175: BblasInfNorm
    "",                                     ///< 176: BblasRealInfNorm
    "Maximum",                              ///< 177: BblasMaxNorm
    "",                                     ///< 178: BblasRealMaxNorm

    "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "", "", "",
    "", "", "", "", "", "", "", "",
    "Forward",                             ///< 391: BblasForward
    "Backward",                            ///< 392: BblasBackward
    "", "", "", "", "", "", "", "",
    "Columnwise",                          ///< 401: BblasColumnwise
    "Rowwise"                              ///< 402: BblasRowwise
};

/***************************************************************************//**
 * @retval LAPACK character constant corresponding to PLASMA constant
 * @ingroup plasma_const
 ******************************************************************************/
static inline char lapack_const(int plasma_const) {
    return lapack_constants[plasma_const][0];
}

#ifdef __cplusplus
}  // extern "C"
#endif

#include "core_s.h"
#include "core_d.h"
#include "core_c.h"
#include "core_z.h"


#endif // ICL_CORE_BLAS_H
