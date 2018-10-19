/**
 *
 * @file bblas_types.h
 *
 * @brief BBLAS typedefs and enumerates.
 *
 *  BBLAS is a software package provided by 
 *  Univ. of Manschester,
 *  Univ. of Tennessee.
 *
 * @version 1.0.0
 * @author Srikara Pranesh
 * @author Samuel  D. Relton
 * @author Pedro   V. Lara
 * @author Mawussi Zounon
 * @date 2016-02-20
 *
 * Contains enumerates for parameter values and typedefs to support various platforms.
 *
 **/
#ifndef BBLAS_TYPES_H
#define BBLAS_TYPES_H

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


/*
 *  BBLAS typedefs
 */
typedef int  bblas_bool; //!< Define bblas_bool
typedef long bblas_index; //!< Define bblas_index
typedef long bblas_size; //!< Define bblas_size
typedef double bblas_double_t; //!< Define bblas_double_t
typedef int bblas_enum_t;


/*
 * BBLAS Complex numbers
 */


#define BBLAS_HAS_COMPLEX_H 1 //!< Is complex arithmetic is available?

#if defined(_WIN32)
# include <float.h>
# if defined(__INTEL_COMPILER)
    /* Fix name conflict within the cabs prototype (_Complex) that
     * conflicts with a C99 keyword.  */
    #define _Complex __ConflictingComplex
    #include <math.h>
    #undef _Complex
    #undef complex
typedef float  _Complex bblas_complex32_t; //!< Define bblas_complex32_t
    typedef double _Complex bblas_complex64_t; //!< Define bblas_complex64_t
# else
    /* Use MS VC complex class */
    #include <complex>
    typedef std::complex<float> bblas_complex32_t; //!< Define bblas_complex32_t
    typedef std::complex<double> bblas_complex64_t; //!< Define bblas_complex64_t
    #undef BBLAS_HAS_COMPLEX_H
# endif
# define isnan _isnan
# define isinf !_finite

#else /* defined(_WIN32) */

    typedef float  _Complex bblas_complex32_t; //!< Define bblas_complex32_t
    typedef double _Complex bblas_complex64_t; //!< Define bblas_complex64_t

#endif /* defined(_WIN32) */

/* Sun doesn't ship the complex.h header. Sun Studio doesn't have it and older GCC compilers don't have it either. */
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC) || defined(sun) || defined(__sun)
#undef BBLAS_HAS_COMPLEX_H
#endif

#if (__GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 1))
#define BBLAS_DEPRECATED  __attribute__((__deprecated__)) //!< Define BBLAS_DEPRECIATED
#else
#define BBLAS_DEPRECATED //!< Define BBLAS_DEPRECIATED
#endif /* __GNUC__ */

#ifdef BBLAS_HAS_COMPLEX_H
#include <complex.h>

#else

#ifdef __cplusplus
extern "C" {
#endif

/* These declarations will not clash with what C++ provides because the names in C++ are name-mangled. */
#if !defined(_WIN32)
extern double cabs(bblas_complex64_t z);
extern bblas_complex64_t conj(bblas_complex64_t z);
#endif
extern float cabsf(bblas_complex32_t z);
extern double cimag(bblas_complex64_t z);
extern double creal(bblas_complex64_t z);

#ifdef __cplusplus
}
#endif
#endif /* defined(BBLAS_HAS_COMPLEX_H) */

#endif  /* BBLAS_TYPES_H */
