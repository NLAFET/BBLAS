include ./make.inc

# BBLAS Source Code
BBLAS_BASE_DIR  = $(shell pwd)
BBLAS_UTILS_DIR = $(BBLAS_BASE_DIR)/utils
BBLAS_SRC_DIR   = $(BBLAS_BASE_DIR)/src
BBLAS_TEST_DIR  = $(BBLAS_BASE_DIR)/test
DEPS            = -I$(BBLAS_BASE_DIR)/include -I$(BBLAS_TEST_DIR)

#BBLAS_UTILS     = $(BBLAS_UTILS_DIR)/xerbla_batch.c
BBLAS_SRC_LIST   = blas_zgemm_batch.c 

BBLAS_SRC=$(addprefix $(BBLAS_SRC_DIR)/, $(BBLAS_SRC_LIST))

TEST_SRC      = testing_utils.c
TEST_SRC_LIST = $(addprefix $(BBLAS_TEST_DIR)/, $(TEST_SRC))

SOURCES       = $(BBLAS_UTILS) $(BBLAS_SRC) $(TEST_SRC_LIST)
SOURCES_Z     = $(SOURCES)
OBJECTS_Z     = $(SOURCES_Z:.c=.o)

all:
	$(MAKE) test_gemm
.DEFAULT_GOAL := all

test_gemm: $(OBJECTS_Z)
	$(CC) $(CFLAGS) $(DEPS) $(BBLAS_TEST_DIR)/test_gemm.c -o $(BBLAS_TEST_DIR)/test_gemm.o
	$(CC) $(OBJECTS_Z) $(BBLAS_TEST_DIR)/test_gemm.o $(LDFLAGS)   -o $(BBLAS_TEST_DIR)/$@

.c.o:
	$(CC) $(CFLAGS) $(DEPS) $<   -o $@

clean:
	rm */*.o
