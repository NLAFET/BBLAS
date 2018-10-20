/**
 *
 * @file
 *
 *  BBLAS is a software package provided by:
 *  University of Tennessee, US,
 *  University of Manchester, UK.
 *
 **/
#ifndef ICL_BBLAS_ERROR_H
#define ICL_BBLAS_ERROR_H

#include <stdio.h>
#include <stdlib.h>
#include "bblas_types.h"
#include "bblas_error.h"
#ifdef __cplusplus
extern "C" {
#endif

/******************************************************************************/
#define bblas_warning(msg) \
    bblas_warning_func_line_file(__func__, __LINE__, __FILE__, msg)

#define bblas_error(msg) \
    bblas_error_func_line_file(__func__, __LINE__, __FILE__, msg)

#define bblas_error_with_code(msg, code) \
    bblas_error_func_line_file_code(__func__, __LINE__, __FILE__, msg, \
                                         code)

#define bblas_fatal_error(msg) \
    bblas_fatal_error_func_line_file(__func__, __LINE__, __FILE__, msg)

/******************************************************************************/
static inline void bblas_warning_func_line_file(
    char const *func, int line, const char *file, const char *msg)
{
    fprintf(stderr,
            "BBLAS WARNING at %d of %s() in %s: %s\n",
            line, func, file, msg);
}

/******************************************************************************/
static inline void bblas_error_func_line_file(
    char const *func, int line, const char *file, const char *msg)
{
    fprintf(stderr,
            "BBLAS ERROR at %d of %s() in %s: %s\n",
            line, func, file, msg);
}

/******************************************************************************/
static inline void bblas_error_func_line_file_code(
    char const *func, int line, const char *file, const char *msg, int code)
{
    fprintf(stderr,
            "BBLAS ERROR at %d of %s() in %s: %s %d\n",
            line, func, file, msg, code);
}

/******************************************************************************/
static inline void bblas_fatal_error_func_line_file(
    char const *func, int line, const char *file, const char *msg)
{
    fprintf(stderr,
            "BBLAS FATAL ERROR at %d of %s() in %s: %s\n",
            line, func, file, msg);
    exit(EXIT_FAILURE);
}

static inline void bblas_set_info(int error_flag, int *info,
                                  int batch_count, int code) {
    switch (error_flag) {
    case BblasErrorsReportAll :
        for( int i=0; i < batch_count; i++) {
            info[i] = code;
        }
        break;
    case BblasErrorsReportGroup :
        info[0] = code;
        break;
    case BblasErrorsReportAny :
        info[0] = code;
        break;
    default :
        bblas_error("illegal value of info");
        info[0] = -1;
    }
}

static inline void bblas_success(int error_flag, int *info,
                                  int batch_count) {
    switch (error_flag) {
    case BblasErrorsReportAll :
        for( int i=0; i < batch_count; i++) {
            info[i] = 0;
        }
        break;
    case BblasErrorsReportGroup :
        info[0] = 0;
        break;
    case BblasErrorsReportAny :
        info[0] = 0;
        break;
    default :
        bblas_error("illegal value of info");
        info[0] = -1;
    }
}
    
#ifdef __cplusplus
}  // extern "C"
#endif

#endif // ICL_BBLAS_ERROR_H
