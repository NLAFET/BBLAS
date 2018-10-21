/**
 *
 * @file
 *
 *  BBLAS is a software package provided by:
 *  University of Manchester, UK,
 *  University of Tennessee, US.
 *
 **/
#ifndef BBLAS_FLOPS_H
#define BBLAS_FLOPS_H

#include "bblas_types.h"

#ifdef __cplusplus
extern "C" {
#endif

//==============================================================================
// Generic formulas come from LAWN 41
// BLAS formulas generally assume alpha == 1 or -1, and beta == 1, -1, or 0;
// otherwise add some smaller order term.
// Some formulas are wrong when m, n, or k == 0; flops should be 0
// (e.g., syr2k, unmqr).
// Formulas may give negative results for invalid combinations of m, n, k
// (e.g., ungqr, unmqr).

//==============================================================================
// Level 2 BLAS
//==============================================================================

//------------------------------------------------------------ gemv
static double fmuls_gemv(double m, double n)
    { return m*n; }

static double fadds_gemv(double m, double n)
    { return m*n; }

static double  flops_zgemv(double m, double n)
    { return 6.*fmuls_gemv(m, n) + 2.*fadds_gemv(m, n); }

static double  flops_cgemv(double m, double n)
    { return 6.*fmuls_gemv(m, n) + 2.*fadds_gemv(m, n); }

static double  flops_dgemv(double m, double n)
    { return    fmuls_gemv(m, n) +    fadds_gemv(m, n); }

static double  flops_sgemv(double m, double n)
    { return    fmuls_gemv(m, n) +    fadds_gemv(m, n); }

//------------------------------------------------------------ symv/hemv
static double fmuls_symv(double n)
    { return fmuls_gemv(n, n); }

static double fadds_symv(double n)
    { return fadds_gemv(n, n); }

static double fmuls_hemv(double n)
    { return fmuls_symv(n); }

static double fadds_hemv(double n)
    { return fadds_symv(n); }

static double  flops_zhemv(double n)
    { return 6.*fmuls_hemv(n) + 2.*fadds_hemv(n); }

static double  flops_chemv(double n)
    { return 6.*fmuls_hemv(n) + 2.*fadds_hemv(n); }

static double  flops_zsymv(double n)
    { return 6.*fmuls_symv(n) + 2.*fadds_symv(n); }

static double  flops_csymv(double n)
    { return 6.*fmuls_symv(n) + 2.*fadds_symv(n); }

static double  flops_dsymv(double n)
    { return    fmuls_symv(n) +    fadds_symv(n); }

static double  flops_ssymv(double n)
    { return    fmuls_symv(n) +    fadds_symv(n); }

//==============================================================================
// Level 3 BLAS
//==============================================================================

//------------------------------------------------------------ gemm
static double fmuls_gemm(double m, double n, double k)
    { return m*n*k; }

static double fadds_gemm(double m, double n, double k)
    { return m*n*k; }

static double  flops_zgemm(double m, double n, double k)
    { return 6.*fmuls_gemm(m, n, k) + 2.*fadds_gemm(m, n, k); }

static double  flops_cgemm(double m, double n, double k)
    { return 6.*fmuls_gemm(m, n, k) + 2.*fadds_gemm(m, n, k); }

static double  flops_dgemm(double m, double n, double k)
    { return    fmuls_gemm(m, n, k) +    fadds_gemm(m, n, k); }

static double  flops_sgemm(double m, double n, double k)
    { return    fmuls_gemm(m, n, k) +    fadds_gemm(m, n, k); }

//------------------------------------------------------------ symm/hemm
static double fmuls_symm(bblas_enum_t side, double m, double n)
{
    return (side == BblasLeft)
        ? fmuls_gemm(m, m, n)
        : fmuls_gemm(m, n, n);
}

static double fadds_symm(bblas_enum_t side, double m, double n)
{
    return (side == BblasLeft)
        ? fadds_gemm(m, m, n)
        : fadds_gemm(m, n, n);
}

static double fmuls_hemm(bblas_enum_t side, double m, double n)
    { return fmuls_symm(side, m, n); }

static double fadds_hemm(bblas_enum_t side, double m, double n)
    { return fadds_symm(side, m, n); }

static double  flops_zhemm(bblas_enum_t side, double m, double n)
    { return 6.*fmuls_hemm(side, m, n) + 2.*fadds_hemm(side, m, n); }

static double  flops_chemm(bblas_enum_t side, double m, double n)
    { return 6.*fmuls_hemm(side, m, n) + 2.*fadds_hemm(side, m, n); }

static double  flops_zsymm(bblas_enum_t side, double m, double n)
    { return 6.*fmuls_symm(side, m, n) + 2.*fadds_symm(side, m, n); }

static double  flops_csymm(bblas_enum_t side, double m, double n)
    { return 6.*fmuls_symm(side, m, n) + 2.*fadds_symm(side, m, n); }

static double  flops_dsymm(bblas_enum_t side, double m, double n)
    { return    fmuls_symm(side, m, n) +    fadds_symm(side, m, n); }

static double  flops_ssymm(bblas_enum_t side, double m, double n)
    { return    fmuls_symm(side, m, n) +    fadds_symm(side, m, n); }

//------------------------------------------------------------ syrk/herk
static double fmuls_syrk(double n, double k)
    { return 0.5*k*n*(n + 1); }

static double fadds_syrk(double n, double k)
    { return 0.5*k*n*(n + 1); }

static double fmuls_herk(double n, double k)
    { return fmuls_syrk(n, k); }

static double fadds_herk(double n, double k)
    { return fadds_syrk(n, k); }

static double  flops_zherk(double n, double k)
    { return 6.*fmuls_herk(n, k) + 2.*fadds_herk(n, k); }

static double  flops_cherk(double n, double k)
    { return 6.*fmuls_herk(n, k) + 2.*fadds_herk(n, k); }

static double  flops_zsyrk(double n, double k)
    { return 6.*fmuls_syrk(n, k) + 2.*fadds_syrk(n, k); }

static double  flops_csyrk(double n, double k)
    { return 6.*fmuls_syrk(n, k) + 2.*fadds_syrk(n, k); }

static double  flops_dsyrk(double n, double k)
    { return    fmuls_syrk(n, k) +    fadds_syrk(n, k); }

static double  flops_ssyrk(double n, double k)
    { return    fmuls_syrk(n, k) +    fadds_syrk(n, k); }

//------------------------------------------------------------ syr2k/her2k
static double fmuls_syr2k(double n, double k)
    { return k*n*n; }

static double fadds_syr2k(double n, double k)
    { return k*n*n + n; }

static double fmuls_her2k(double n, double k)
    { return fmuls_syr2k(n, k); }

static double fadds_her2k(double n, double k)
    { return fadds_syr2k(n, k); }

static double  flops_zher2k(double n, double k)
    { return 6.*fmuls_her2k(n, k) + 2.*fadds_her2k(n, k); }

static double  flops_cher2k(double n, double k)
    { return 6.*fmuls_her2k(n, k) + 2.*fadds_her2k(n, k); }

static double  flops_zsyr2k(double n, double k)
    { return 6.*fmuls_syr2k(n, k) + 2.*fadds_syr2k(n, k); }

static double  flops_csyr2k(double n, double k)
    { return 6.*fmuls_syr2k(n, k) + 2.*fadds_syr2k(n, k); }

static double  flops_dsyr2k(double n, double k)
    { return    fmuls_syr2k(n, k) +    fadds_syr2k(n, k); }

static double  flops_ssyr2k(double n, double k)
    { return    fmuls_syr2k(n, k) +    fadds_syr2k(n, k); }

//------------------------------------------------------------ trmm
static double fmuls_trmm_2(double m, double n)
    { return 0.5*n*m*(m + 1); }

static double fadds_trmm_2(double m, double n)
    { return 0.5*n*m*(m - 1); }

static double fmuls_trmm(bblas_enum_t side, double m, double n)
{
    return (side == BblasLeft)
        ? fmuls_trmm_2(m, n)
        : fmuls_trmm_2(n, m);
}

static double fadds_trmm(bblas_enum_t side, double m, double n)
{
    return (side == BblasLeft)
        ? fadds_trmm_2(m, n)
        : fadds_trmm_2(n, m);
}

static double  flops_ztrmm(bblas_enum_t side, double m, double n)
    { return 6.*fmuls_trmm(side, m, n) + 2.*fadds_trmm(side, m, n); }

static double  flops_ctrmm(bblas_enum_t side, double m, double n)
    { return 6.*fmuls_trmm(side, m, n) + 2.*fadds_trmm(side, m, n); }

static double  flops_dtrmm(bblas_enum_t side, double m, double n)
    { return    fmuls_trmm(side, m, n) +    fadds_trmm(side, m, n); }

static double  flops_strmm(bblas_enum_t side, double m, double n)
    { return    fmuls_trmm(side, m, n) +    fadds_trmm(side, m, n); }

//------------------------------------------------------------ trsm
static double fmuls_trsm(bblas_enum_t side, double m, double n)
    { return fmuls_trmm(side, m, n); }

static double fadds_trsm(bblas_enum_t side, double m, double n)
    { return fadds_trmm(side, m, n); }

static double  flops_ztrsm(bblas_enum_t side, double m, double n)
    { return 6.*fmuls_trsm(side, m, n) + 2.*fadds_trsm(side, m, n); }

static double  flops_ctrsm(bblas_enum_t side, double m, double n)
    { return 6.*fmuls_trsm(side, m, n) + 2.*fadds_trsm(side, m, n); }

static double  flops_dtrsm(bblas_enum_t side, double m, double n)
    { return    fmuls_trsm(side, m, n) +    fadds_trsm(side, m, n); }

static double  flops_strsm(bblas_enum_t side, double m, double n)
    { return    fmuls_trsm(side, m, n) +    fadds_trsm(side, m, n); }

#ifdef __cplusplus
}  // extern "C"
#endif

#endif // BBLAS_FLOPS_H
