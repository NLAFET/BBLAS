**Batched Basic Linear Algebra Subroutines**

**University of Manchester (UK)**

**University of Tennessee (US)**


* * *

[Download BBLAS Software](https://github.com/NLAFET/BBLAS/archive/master.zip)

* * *

About
=====

BBLAS is the reference implementation of the Batched BLAS standard
specification.  A current trend in high-performance computing is to
decompose a large linear algebra problem into batches containing
thousands of smaller problems, which can be solved independently,
before collating the results. To standardize the interface to these
routines, the community  developed  an extension to the BLAS
standard (the batched BLAS), enabling users to perform thousands of
small BLAS operations in parallel whilst making efficient use of their
hardware. Please visit [BBLAS workshops](http://icl.utk.edu/bblas)
for more information on the standardization efforts. 

The main folders & files 
========================

* **compute**: contains the standard BBLAS group API functions 
* **core**:    contains the auxiliary batched BLAS functions
               which perform on groups of same size problems

* **control**: contains auxiliary functions for type conversions

* **test**: contains testing routines associated with the BBLAS functions
            and provides an insight on how the BBLAS functions should be called/used.

* **include**: contains header files

* _make.inc_: a configuration file to specify a C/C++ compiler,
                compilation flags and a BLAS library. The default
                configuration should work when MKL is installed.

* _Makefile_: the Makefile, normally, it should  not be modified.

Requirements
===========
#### BLAS & LAPACK
* **MKL** is now free for academics (students and researchers) available at [https://software.intel.com/en-us/articles/free-mkl](https://software.intel.com/en-us/articles/free-mkl)

OR

* **Netlib BLAS** no optimized BLAS routines, available at [BLAS-3.8.0.tgz](http://www.netlib.org/blas/blas-3.8.0.tgz)
* **Netlib LAPACK** no optimized LAPACK routines, available at [LAPACK-3.8.0.tgz](http://www.netlib.org/lapack/lapack-3.8.0.tar.gz)

#### Doxygen for documentation 
* **Doxygen** can be install on Unix systems by _sudo apt-get install doxygen_ or
downloaded on [the doxygen page](http://www.doxygen.org/download.html).

Compilation 
===========
After the configuration of **make.inc**, the compilation is very simple:

*  **make [all]**     --  make lib test
*  **make lib**      --  make lib/libbblas.{a,so} lib/libcore.{a,so}
*  **make test**      --  make test/test
*  **make docs**      --  make docs/html
*  **make generate**  --  generate precisions
*  **make clean**     --  remove objects, libraries, and executables
*  **make cleangen**  --  remove generated precision files
*  **make distclean** --  remove above, Makefile.*.gen, and anything else that can be generated 


Citing
======

Feel free to use the following publications to reference BBLAS:

* Jack  Dongarra, Sven Hammarling, Nicholas J. Higham,
  Samuel D. Relton, Mawussi Zounon:
  **Optimized Batched Linear Algebra for Modern Architectures.**
  *Euro-Par 2017: 511-522*.

* Jack  Dongarra, Sven Hammarling, Nicholas J. Higham,
  Samuel D. Relton, Pedro Valero-Lara, Mawussi Zounon:
  **The Design and Performance of Batched BLAS on Modern High-Performance Computing Systems**,
  *ICCS 2017: 495-504*


* Jack Dongarra, Iain Duff, Mark Gates, Azzam Haidar,
  Sven Hammarling, Nicholas J. Higham, Jonathan Hogg,
  Pedro Valero Lara, Mawussi Zounon, Samuel D. Relton,
  and Stanimire Tomov,
  **A Proposed API for Batched Basic Linear Algebra Subprograms**,
  [Draft Report, May 2016.](https://www.dropbox.com/s/olocmipyxfvcaui/batched_api_03_30_2016.pdf?dl=0)


Funding
=======

Primary funding for BBLAS was provided  the European Union grant:

* [NLAFET: Parallel Numerical Linear Algebra for Future Extreme Scale Systems](http://www.nlafet.eu), Grant Agreement no. 671633


People
======

The following people listed in alphabetical order contributed to the BBLAS reference implementation:

* Jack Dongarra
* Mark Gates
* Srikara Pranesh
* Samuel Relton
* Pedro Valero Lara
* Mawussi Zounon

