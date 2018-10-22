/**
 *
 * @file
 *
 *  BBLAS is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 **/

#include "bblas_types.h"

#include <assert.h>
#include <stdbool.h>

/***************************************************************************//**
 * @addtogroup bblas_const
 * Convert LAPACK character constants to BBLAS constants.
 * This is a one-to-many mapping, requiring multiple translators
 * (e.g., "N" can be NoTrans or NonUnit or NoVec).
 * Matching is case-insensitive.
 * @{
 ******************************************************************************/

// These functions and cases are in the same order as the constants are
// declared in bblas_types.h

/***************************************************************************//**
 * @retval BblasNonUnit if lapack_char = 'N'
 * @retval BblasUnit    if lapack_char = 'U'
 ******************************************************************************/
bblas_enum_t bblas_diag_const(char lapack_char)
{
    switch (lapack_char) {
    case 'N': case 'n': return BblasNonUnit;
    case 'U': case 'u': return BblasUnit;
    default:            return BblasInvalid;
    }
}

/***************************************************************************//**
 * @retval BblasErrorsReportAll       if lapack_char = 'a'
 * @retval BblasErrorsReportGroup     if lapack_char = 'g'
 * @retval BblasErrorsReportAny       if lapack_char = 'o'
 * @retval BblasErrorsReportNone      if lapack_char = 'n'
 ******************************************************************************/
bblas_enum_t bblas_info_const(char lapack_char)
{
    switch (lapack_char) {
    case 'A': case 'a': return BblasErrorsReportAll;
    case 'G': case 'g': return BblasErrorsReportGroup;
    case 'O': case 'o': return BblasErrorsReportAny;
    case 'N': case 'n': return BblasErrorsReportNone;
    default:            return BblasInvalid;
    }
}


/***************************************************************************//**
 * @retval BblasForward  if lapack_char = 'F'
 * @retval BblasBackward if lapack_char = 'B'
 ******************************************************************************/
bblas_enum_t bblas_direct_const(char lapack_char)
{
    switch (lapack_char) {
    case 'F': case 'f': return BblasForward;
    case 'B': case 'b': return BblasBackward;
    default:            return BblasInvalid;
    }
}

/***************************************************************************//**
 * @retval BblasOneNorm       if lapack_char = 'O|o|1'
 * @retval BblasTwoNorm       if lapack_char = '2'
 * @retval BblasFrobeniusNorm if lapack_char = 'F|f|E|e'
 * @retval BblasInfNorm       if lapack_char = 'I|i'
 * @retval BblasMaxNorm       if lapack_char = 'M|m'
 ******************************************************************************/
bblas_enum_t bblas_norm_const(char lapack_char)
{
    switch (lapack_char) {
    case 'O': case 'o': case '1':           return BblasOneNorm;
    case '2':                               return BblasTwoNorm;
    case 'F': case 'f': case 'E': case 'e': return BblasFrobeniusNorm;
    case 'I': case 'i':                     return BblasInfNorm;
    case 'M': case 'm':                     return BblasMaxNorm;
    // MagmaRealOneNorm
    // MagmaRealInfNorm
    // MagmaRealMaxNorm
    default:                                return BblasInvalid;
    }
}

/***************************************************************************//**
 * @retval BblasLeft      if lapack_char = 'L'
 * @retval BblasRight     if lapack_char = 'R'
 ******************************************************************************/
// @retval BblasBothSides if lapack_char = 'B'  // for trevc
bblas_enum_t bblas_side_const(char lapack_char)
{
    switch (lapack_char) {
    case 'L': case 'l': return BblasLeft;
    case 'R': case 'r': return BblasRight;
    //case 'B': case 'b': return BblasBothSides;  // for trevc
    default:            return BblasInvalid;
    }
}

/***************************************************************************//**
 * @retval BblasColumnwise if lapack_char = 'C'
 * @retval BblasRowwise    if lapack_char = 'R'
 ******************************************************************************/
bblas_enum_t bblas_storev_const(char lapack_char)
{
    switch (lapack_char) {
    case 'C': case 'c': return BblasColumnwise;
    case 'R': case 'r': return BblasRowwise;
    default:            return BblasInvalid;
    }
}

/***************************************************************************//**
 * @retval BblasNoTrans   if lapack_char = 'N'
 * @retval BblasTrans     if lapack_char = 'T'
 * @retval BblasConjTrans if lapack_char = 'C'
 ******************************************************************************/
bblas_enum_t bblas_trans_const(char lapack_char)
{
    switch (lapack_char) {
    case 'N': case 'n': return BblasNoTrans;
    case 'T': case 't': return BblasTrans;
    case 'C': case 'c': return BblasConjTrans;
    default:            return BblasInvalid;
    }
}

/***************************************************************************//**
 * @retval BblasUpper   if lapack_char = 'U'
 * @retval BblasLower   if lapack_char = 'L'
 * @retval BblasGeneral otherwise
 ******************************************************************************/
bblas_enum_t bblas_uplo_const(char lapack_char)
{
    switch (lapack_char) {
    case 'U': case 'u': return BblasUpper;
    case 'L': case 'l': return BblasLower;
    default:            return BblasGeneral;
    }
}

/***************************************************************************//**
 * @}
 * end group bblas_const
 ******************************************************************************/
