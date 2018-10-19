/**
 * @file auxiliary.h
 *
 * @brief BBLAS auxiliary routines, should be unused.
 *
 * BBLAS is a software package provided by 
 * Univ. of Manchester
 * Univ. of Tennessee.
 *
 * @version 1.0.0
 * @date 2016-02-20
 *
 **/
#ifndef BBLAS_AUXILIARY_H
#define BBLAS_AUXILIARY_H


/*-------------------------
 *  Auxiliary routines
 *------------------------*/
/* Review*/
void bblas_warning(const char *func_name, char* msg_text);
void bblas_error(const char *func_name, char* msg_text);
void bblas_fatal_error(const char *func_name, char* msg_text);
void bblas_malloc_check(void *ptr, char* msg_text);
#endif
