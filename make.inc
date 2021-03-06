# BBLAS example make.inc, using Intel MKL and gcc
#
# BBLAS is a software package provided by:
# University of Tennessee, US,
# University of Manchester, UK.

# --------------------
# programs

CC        = gcc
AR        = ar
RANLIB    = ranlib


# --------------------
# flags

# Use -fPIC to make shared (.so) and static (.a) libraries;
# can be commented out if making only static libraries.
FPIC      = -fPIC

CFLAGS    =  $(FPIC) -O3 -std=c99   -Wall -pedantic -Wshadow -Wno-unused-function

LDFLAGS   =  $(FPIC)

# one of: aix bsd c89 freebsd generic linux macosx mingw posix solaris
# usually generic is fine
lua_platform = macosx


# --------------------
# libraries

# ----- with MKL
# This assumes $MKLROOT is set in your environment.
# Add these to your .cshrc or .bashrc, adjusting for where MKL is installed:
# in .cshrc:   source /opt/intel/bin/compilervars.csh intel64
# in .bashrc:  source /opt/intel/bin/compilervars.sh  intel64
# Currently, MKL on MacOS has only Intel threads, which may not work with gcc,
# so use sequential.
CFLAGS   += -DHAVE_MKL

LIBS      = -L$(MKLROOT)/lib -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lm
INC       = -I$(MKLROOT)/include


# ----- with MacOS Accelerate (Veclib)
# do not use Netlib LAPACK with Accelerate; zdotc crashes.
#LIBS      = -L$(LAPACKDIR) -llapacke -framework Accelerate
#INC       = -I$(LAPACKDIR)/LAPACKE/include -I$(CBLASDIR)
