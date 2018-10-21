/**
 *
 * @file
 *
 *  BBLAS is a software package provided by:
 *  University of Manchester, UK,
 *  University of Tennessee, US.
 *
 * @precisions normal z -> s d c
 *
 **/
#include "test.h"
#include "flops.h"
#include "bblas.h"

#include <assert.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#define COMPLEX

/***************************************************************************//**
 *
 * @brief Tests BATCHED ZGEMM.
 *
 * @param[in,out] param - array of parameters
 * @param[in]     run - whether to run test
 *
 * Sets used flags in param indicating parameters that are used.
 * If run is true, also runs test and stores output parameters.
 ******************************************************************************/
void test_zgemm_batch(param_value_t param[], bool run)
{
    //================================================================
    // Mark which parameters are used.
    //================================================================
    param[PARAM_NG     ].used = true;
    param[PARAM_GS     ].used = true;
    param[PARAM_INCM   ].used = true;
    param[PARAM_INCG   ].used = true;
    param[PARAM_TRANSA ].used = true;
    param[PARAM_TRANSB ].used = true;
    param[PARAM_DIM    ].used = PARAM_USE_M | PARAM_USE_N | PARAM_USE_K;
    param[PARAM_ALPHA  ].used = true;
    param[PARAM_BETA   ].used = true;
    param[PARAM_PADA   ].used = true;
    param[PARAM_PADB   ].used = true;
    param[PARAM_PADC   ].used = true;
    if (! run)
        return;

    /* 
      [1] to fill see plasma/test/test_zgemm.c for inspiration,
      [2] use ./test zgemm_batch -h for help,
      [3] we don't depend on openMP any more so I have added a new function called gettime() in place.
    */
}
