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

*  **make [all]**     --  compiles lib test
*  **make lib**       --  generates library files in lib/libbblas.{a,so} lib/libcore.{a,so}
*  **make test**      --  generates tester files in  test/test
*  **make docs**      --  generates documentation docs/html
*  **make generate**  --  generates routines of other precisions
*  **make clean**     --  removes objects, libraries, and executables
*  **make cleangen**  --  removes generated precision files
*  **make distclean** --  removes above, Makefile.*.gen, and anything else that can be generated 

For MacOS the libraries have to be linked manually. Use the following in the .profile file.

export DYLD_FALLBACK_LIBRARY_PATH=/path/to/BBLAS/lib:$DYLD_FALLBACK_LIBRARY_PATH

Testing
=======
At the end of the compilation, testing routines (in four precisions [c,d,s,z]) will be available in the folder BBLAS/test.
The main testing driver is the binary **./test**. It should be used with the kernel to test as the firt parameter followed 
by the kernel arguments. For example **./test dgemm_batch** will run the double precision of version of **dgemm_batch** with 
default arguments. For help **./test -h** will display the list of kernel availble for testing, while **./test dgemm_batch -h** will display all the possible more help and details on how to test **dgemm_batch**, this holds for all the kernels. 
We use random matrices for the test. To simplify the way to provide matrices in a batch and the size of each group, we provide 
the following arguments:
* **--ng**   : the number groups  [default: --ng=10].
* **--gs**   : the number matrices in the first group [default: --gs=100].
* **--incg** : the increment of group sizes. The size of the i^th group is gs + i * incg [default: --incg=10].
* **--dim**  : M x N x K dimensions of th matrices in the first group [default: --dim=1000 x 1000 x 1000].
* **--incm** : The increment of matrix size across the group. If the matrix size in the first group is M x N x K,
the matrix size in the i^th group will be (M + incm * i) x (N + incm * i) x (K + incm * i)  [default: --incm = 1].
* **--info[a|g|n|o]** : The parameter to set an error handling option [default: --info=a].
  - **a**: which indicates that all errors will be specified on output.
  - **g**: which indicates that only a single error will be reported for each group, independently.
  - **n**: which indicates that no errors will be reported on output.
  - **o**: which indicates that the occurrence of errors will be specified on output as a single integer value.

In addition to these options, arguments like **trans, transa, transb, side, diag, uplo, etc.**, can be set, and
for the sake of simplicity, they have the same value in all the groups.

* **--transa=[n|t|c]** : transposition      [default: --transa=n].
* **--transa=[n|t|c]** : transposition of A [default: --transa=n].
* **--transb=[n|t|c]** : transposition of B [default: --transb=n].
* **--side=[l|r]**     : left or right side application [default: --side=l].
* **--uplo=[g|u|l]**   : general rectangular or upper or lower triangular matrix [default: --uplo=l].
* **--diag=[n|u]**     : non-unit diagonal or unit diagonal [default: --diag=n].
                                                            
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

