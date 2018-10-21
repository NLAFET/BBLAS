/**
 *
 * @file
 *
 *  PLASMA is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 * @precisions normal z -> s d c
 *
 **/
#ifndef TEST_Z_H
#define TEST_Z_H

#include "test.h"

//==============================================================================
// test routines
//==============================================================================

void test_zgemm_batch(param_value_t param[], bool run);
void test_zhemm_batch(param_value_t param[], bool run);
void test_zher2k_batch(param_value_t param[], bool run);
void test_zherk_batch(param_value_t param[], bool run);
void test_zsymm_batch(param_value_t param[], bool run);
void test_zsyr2k_batch(param_value_t param[], bool run);
void test_zsyrk_batch(param_value_t param[], bool run);
void test_ztrmm_batch(param_value_t param[], bool run);
void test_ztrsm_batch(param_value_t param[], bool run);


#endif // TEST_Z_H
