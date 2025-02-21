
Download and Installation
==================================

=================
Download
=================

The FLAME software is available on the flame website http://flame-code.org, and the
latest development version can be found on GitHub https://github.com/flame-code
and on GitLab https://gitlab.com/flame-code.

Download the current version of FLAME with 

``git clone https://github.com/flame-code/FLAME.git``

or

``git clone https://gitlab.com/flame-code/FLAME.git``

=================
Prerequisites
=================

FLAME requires **autoconf** and **automake**.
IMPORTANT NOTE:currently only **automake** up to version 1.15.1
is supported due to changes introduced in later versions
that break the Makefile structure.

Any Fortran and C compiler should in principle work for compiling FLAME.
However, we recommend using the Intel Fortran and C compiler.


FLAME has to be linked to Blas, LaPack, and FFT libraries.  
They can be obtained as part of the Intel Math Kernel Library (MKL), 
which is the recommended route. In principle, other
implementations of the libraries should also work.


Linking to Atsushi Togo's SPGLIB is recommended. The currently supported
and well-tested version is 1.6.x and can be found here:

https://sourceforge.net/projects/spglib/files/spglib/


Linking to LAMMPS requires the installation of LAMMPS with 
the desired potentials. The best upported and well tested version is
r12824:

http://lammps.sandia.gov

Futile is required, a library of tools developed as part of the BigDFT project.

http://bigdft.org/

Installation of python is required. Currently,
only python 2.7 is supported. Future releases of FLAME will
support python 3


=========================
Installing FLAME 
=========================


On Linux
----------------

It is recommended to install FLAME in a different
directory than the source code.

#. First, install futile which is
   a set of utilities from the BigDFT project.
   Preferably, use the version provided with
   FLAME to avoid conflicts.
   Untar the included futile-suite.tar.gz (``tar -zxvf futile-suite.tar.gz``), then 
   create a new build directory (e.g., ``mkdir futile-build ; cd futile-build``), and from there run

   for GNU compilers:

      ``path_to_futile_source/Installer.py build futile -c 
      CC=gcc CXX=g++ FC=gfortran F77=gfortran \
      --with-ext-linalg="-llapack -lblas"``

   for Intel compilers:

      ``path_to_futile_source/Installer.py build futile -c FCFLAGS=-O2 \
      --with-ext-linalg="-L$MKLROOT/lib/intel64 \
      -lmkl_rt -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 -liomp5 -lm" \
      CC=icc CXX=icpc FC=ifort F77=ifort``

   for parallel compilation use the corresponding MPI wrappers:

      ``path_to_futile_source/Installer.py build futile -c FCFLAGS=-O2 \
      --with-ext-linalg="-L$MKLROOT/lib/intel64 \
      -lmkl_rt -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 -liomp5 -lm" \
      CC=mpicc CXX=mpicxx FC=mpif90 F77=mpif90``

   follwed by ``make build`` if necessary.
   Make sure to adapt the library locations and
   the linking flags appropriately.
   After the installation of futile, we need to link it
   with FLAME.
   To display details on the general linking procedure, run:

   ``path_to_futile_source/Installer.py link futile``

#. To compile FLAME, change into the main FLAME directory and run:

   ``autoreconf -fi``

#. Create a build directory for FLAME (e.g., ``mkdir build-FLAME ; cd build-FLAME``). 
   Explicitly replace $FUTILE with the full path of the futile build-directory during ``configure``, 
   or define it as an environmental variable:

   ``export FUTILE=path_to_futile_build``

   Note that providing the FUTILE variable is required to successfully compile FLAME, and is not optional.
   Then, run ``configure``. For Intel compilers and MPI parallelization:

   ``path_to_flame_source/configure FC=mpif90 F77=mpif90 CXX=mpicc CC=mpicc \
   FCFLAGS="-I$FUTILE/install/include -shared-intel -mcmodel=large  -mkl=sequential" \
   CFLAGS=-mcmodel=large "LIBS=-L$FUTILE/install/lib \
   -L/$MKLROOT/lib/intel64 -lfutile-1 -lmkl_rt \
   -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 -liomp5 -lm -lyaml -ldl -lrt -cxxlib"``


   For GNU compilers without MKL:

   ``path_to_flame_source/configure FC=mpif90 F77=mpif90 CXX=mpicc CC=mpicc \
   FCFLAGS="-I$FUTILE/install/include -mcmodel=large" \
   CFLAGS=-mcmodel=large "LIBS=-L$FUTILE/install/lib \
   -lfutile-1 -lm -lyaml -llapack -lfftw3 -ldl -cxxlib"``
   
   To link with SPGLIB, append
   ``--with-spglib SPGLIB_ROOT=path_to_spglib``

   To link with the BigDFT PSolver, append
   ``--with-bps BDIR=path_to_bigdft_root``

   To link with LAMMPS, append
   ``--with-lammps LAMMPS_ROOT=path_to_lammps_root``

#. Run ``make`` to compile the code. 
   Upon successful compilation, the executable can be found in ``src/flame``
