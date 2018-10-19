# ------------------------------------------------------------------------------
# Usage:
#   make [all]      --  make lib test
#   make lib        --  make lib/libbblas.{a,so} lib/libcore.{a,so}
#   make test       --  make test/test
#   make docs       --  make docs/html
#   make generate   --  generate precisions
#   make clean      --  remove objects, libraries, and executables
#   make cleangen   --  remove generated precision files
#   make distclean  --  remove above, Makefile.*.gen, and anything else that can be generated
#
# Options:
#   quiet=1         --  abbreviate long compilation commands
#   static=0        --  build only shared libraries (no static)
#   static=1        --  build only static libraries (no shared)
#                       By default, builds both static and shared libraries.


# ------------------------------------------------------------------------------
# Define default rule first, before including Makefiles

all: lib test


# ------------------------------------------------------------------------------
# Tools and flags
# Variables defined in make.inc, or use make's defaults:
#   CC, CFLAGS, INC -- C compiler and flags
#   LDFLAGS, LIBS   -- Linker options, library paths, and libraries
#   AR, RANLIB      -- Archiver, ranlib updates library TOC
#   prefix          -- where to install BBLAS

include make.inc

# dependencies here interfere with manually edited make.inc
make.inc: #make.inc.in configure.py $(wildcard config/*.py)
	python configure.py

# GNU make doesn't have defaults for these
RANLIB    ?= ranlib

prefix    ?= /usr/local/bblas

# one of: aix bsd c89 freebsd generic linux macosx mingw posix solaris
# usually generic is fine

ifeq ($(FCFLAGS),)
    ifneq ($(FFLAGS),)
        $(warning Warning: FFLAGS renamed FCFLAGS, per autoconf. Update make.inc.)
    endif
endif


# ------------------------------------------------------------------------------
# Internal tools and flags

codegen     := ./tools/codegen.py

BBLAS_INC  := -Iinclude 
BBLAS_LIBS := -Llib -lbblas -lcore 

.DELETE_ON_ERROR:

.SUFFIXES:

ifeq ($(quiet),1)
   quiet_CC = @echo "$(CC) ... $@";
   quiet_FC = @echo "$(FC) ... $@";
   quiet_AR = @echo "$(AR) ... $@";
endif


# ------------------------------------------------------------------------------
# Define sources, objects, libraries, executables.
# These makefiles define lists of source and header files in
# $(bblas_all), $(core_all), and $(test_all).

makefiles_gen := \
	Makefile.bblas.gen   \
	Makefile.core.gen \
	Makefile.test.gen     

-include $(makefiles_gen)

bblas_hdr   := $(filter %.h, $(bblas_all))
core_hdr := $(filter %.h, $(core_all))
test_hdr     := $(filter %.h, $(test_all))
headers      := $(bblas_hdr) $(core_hdr) $(test_hdr)

bblas_obj   := $(addsuffix .o, $(basename $(filter-out %.h, $(bblas_all))))
core_obj := $(addsuffix .o, $(basename $(filter-out %.h, $(core_all))))
test_obj     := $(addsuffix .o, $(basename $(filter-out %.h, $(test_all))))

test_exe     := test/test


# ------------------------------------------------------------------------------
# Build static libraries

ifneq ($(static),0)
    libfiles += \
	lib/libbblas.a    \
	lib/libcore.a  \

    # In case changing Makefile.gen changes $(obj), also depend on it,
    # which recreates the library if a file is removed.
    lib/libbblas.a: $(bblas_obj) Makefile.bblas.gen
	-rm -f $@
	$(quiet_AR) $(AR) cr $@ $(bblas_obj)
	$(RANLIB) $@

    lib/libcore.a: $(core_obj) Makefile.core.gen
	-rm -f $@
	$(quiet_AR) $(AR) cr $@ $(core_obj)
	$(RANLIB) $@
endif


# ------------------------------------------------------------------------------
# Build shared libraries

# if all FLAGS have -fPIC, allow compiling shared libraries
have_fpic = $(and $(findstring -fPIC, $(CFLAGS)),   \
                  $(findstring -fPIC, $(LDFLAGS)))

ifneq ($(static),1)
ifneq ($(have_fpic),)
   libfiles += \
	lib/libbblas.so    \
	lib/libcore.so

   top = $(shell pwd)

   rpath = -Wl,-rpath,$(top)/lib

   # MacOS (darwin) needs shared library's path set
   ifneq ($(findstring darwin, $(OSTYPE)),)
       install_name = -install_name @rpath/$(notdir $@)
   else
       install_name =
   endif

   lib/libbblas.so: $(bblas_obj) Makefile.bblas.gen lib/libcore.so
	$(quiet_CC) $(CC) -shared $(LDFLAGS) -o $@ $(bblas_obj) \
	-Llib -lcore \
	$(install_name)

   lib/libcore.so: $(core_obj) Makefile.core.gen
	$(quiet_CC) $(CC) -shared $(LDFLAGS) -o $@ $(core_obj) \
	-Llib $(LIBS) \
	$(install_name)
endif
endif

# creare the lib
libdir:
	mkdir -p lib

$(libfiles): libdir

.PHONY: lib

lib: $(libfiles)


# ------------------------------------------------------------------------------
# Build tester

.PHONY: test

test: $(test_exe)

$(test_exe): $(test_obj) $(libfiles) Makefile.test.gen
	$(quiet_CC) $(CC) $(LDFLAGS) -o $@ $(test_obj) \
	$(BBLAS_LIBS) \
	$(LIBS) \
	$(rpath)


# ------------------------------------------------------------------------------
# Build objects
# Headers must exist before compiling, but use order-only prerequisite
# (after "|") so as not to force recompiling everything if a header changes.
# (Should use compiler's -MMD flag to create header dependencies.)

%.o: %.c | $(headers)
	$(quiet_CC) $(CC) $(CFLAGS) $(BBLAS_INC) $(INC) -c -o $@ $<

%.i: %.c | $(headers)
	# $(quiet_CC) $(CC) $(CFLAGS) $(BBLAS_INC) $(INC) -E -o $@ $<
	$(quiet_CC) $(CC) $(CFLAGS) $(BBLAS_INC) $(INC) -o $@ $<


# ------------------------------------------------------------------------------
# Build documentation

.PHONY: docs

docs: generate
	doxygen docs/doxygen/doxyfile.conf


# ------------------------------------------------------------------------------

.PHONY: install_dirs install

install_dirs:
	mkdir -p $(prefix)
	mkdir -p $(prefix)/include
	mkdir -p $(prefix)/lib
	mkdir -p $(prefix)/lib/pkgconfig

install: lib install_dirs
	cp include/*.h $(prefix)/include
	cp $(libfiles) $(prefix)/lib
	# pkgconfig
	cat lib/pkgconfig/bblas.pc.in         | \
	sed -e s:@INSTALL_PREFIX@:"$(prefix)": | \
	sed -e s:@CFLAGS@:"$(INC)":            | \
	sed -e s:@LIBS@:"$(LIBS)":             | \
	sed -e s:@REQUIRES@::                    \
	    > $(prefix)/lib/pkgconfig/bblas.pc

uninstall:
	rm -f $(prefix)/include/core*.h
	rm -f $(prefix)/include/bblas*.h
	rm -f $(prefix)/lib/libcore*
	rm -f $(prefix)/lib/libbblas*
	rm -f $(prefix)/lib/pkgconfig/bblas.pc


# ------------------------------------------------------------------------------
# Maintenance rules
# makefiles_gen define generate and cleangen.

.PHONY: clean distclean

clean:
	-rm -f $(bblas_obj) $(core_obj) $(test_obj) $(test_exe) $(libfiles)
	-rmdir lib

# cleangen removes generated files if the template still exists;
# grep for any stale generated files without a template.
distclean: clean cleangen
	grep -s -l @generated $(bblas_src) $(core_src) $(test_src) | xargs rm -f
	-rm -f compute/*.o control/*.o  core/*.o test/*.o
	-rm -f $(makefiles_gen)
	-rm -rf docs/html


# ------------------------------------------------------------------------------
# Create dependencies to do precision generation.

bblas_src   := $(wildcard compute/*.c include/bblas*.h)

core_src := $(wildcard core/*.c  control/*.c include/core*.h)

test_src     := $(wildcard test/*.c test/*.h)

Makefile.bblas.gen: $(codegen)
	$(codegen) --make --prefix bblas   $(bblas_src)   > $@

Makefile.core.gen: $(codegen)
	$(codegen) --make --prefix core $(core_src) > $@

Makefile.test.gen: $(codegen)
	$(codegen) --make --prefix test     $(test_src)     > $@

# --------------------
# If the list of src files changes, then force remaking Makefile.gen
# To reduce unnecesary remaking, don't remake if either:
# 1) src == old:
#    src has same files now as when Makefile.gen was generated, or
# 2) src - generated == templates:
#    src has all the templates from Makefile.gen, and no new non-generated files.
ifneq ($(bblas_src),$(bblas_old))
ifneq ($(filter-out $(bblas_generated),$(bblas_src)),$(bblas_templates))
    Makefile.bblas.gen: force_gen
endif
endif

ifneq ($(core_src),$(core_old))
ifneq ($(filter-out $(core_generated),$(core_src)),$(core_templates))
    Makefile.core.gen: force_gen
endif
endif

ifneq ($(test_src),$(test_old))
ifneq ($(filter-out $(test_generated),$(test_src)),$(test_templates))
    Makefile.test.gen: force_gen
endif
endif

# --------------------

force_gen: ;


# ------------------------------------------------------------------------------
# Debugging

echo:
	@echo "CC      $(CC)"
	@echo "CFLAGS  $(CFLAGS)"
	@echo "LDFLAGS $(LDFLAGS)"
	@echo
	@echo "test_exe           <$(test_exe)>"
	@echo "libfiles           <$(libfiles)>"
	@echo
	@echo "bblas_src         <$(bblas_src)>"
	@echo "bblas_old         <$(bblas_old)>"
	@echo "bblas_templates   <$(bblas_templates)>"
	@echo "bblas_filtered    <$(filter-out $(bblas_generated),$(bblas))>"
	@echo "bblas_hdr         <$(bblas_hdr)>"
	@echo "bblas_obj         <$(bblas_obj)>"
	@echo
	@echo "core_src       <$(core_src)>"
	@echo "core_old       <$(core_old)>"
	@echo "core_templates <$(core_templates)>"
	@echo "core_filtered  <$(filter-out $(core_generated),$(core))>"
	@echo "core_hdr       <$(core_hdr)>"
	@echo
	@echo "test_src           <$(test_src)>"
	@echo "test_old           <$(test_old)>"
	@echo "test_templates     <$(test_templates)>"
	@echo "test_filtered      <$(filter-out $(test_generated),$(test))>"
	@echo "test_hdr           <$(test_hdr)>"
	@echo
	@echo "headers            <$(headers)>"
