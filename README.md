<!--
<div align="center">
  <img src="https://www.iasbs.ac.ir/~aghasemi/images/logo.png"><br><br>
</div>
-->

[comment]: # (-----------------)

# FLAME: a library for atomistic modeling environments

[FLAME](flame-code.org) is a highly modular open source software package to perform atomistic simulations using a variety of techniques.


## Prerequisites

FLAME requires `autoconf` and `automake`.
IMPORTANT NOTE:currently only `automake` up to version 1.15.1
is supported due to changes introduced in later versions
that break the Makefile structure.

Any Fortran and C compiler should in principle work for compiling FLAME.
However, we recommend using the Intel Fortran and C compiler.


FLAME has to be linked to Blas, LaPack, and FFT libraries.  
They can be obtained as part of the Intel Math Kernel Library (MKL), 
which is the recommended route. In principle, other
implementations of the libraries should also work.


Linking to Atsushi Togo's [SPGLIB](https://sourceforge.net/projects/spglib/files/spglib/) is recommended. 
The currently supported
and well-tested version is 1.6.x and can be found here:

https://sourceforge.net/projects/spglib/files/spglib/


Linking to LAMMPS requires the installation of [LAMMPS](http://lammps.sandia.gov) with 
the desired potentials. The best upported and well tested version is
r12824:

http://lammps.sandia.gov

Futile is required, a library of tools developed as part of the [BigDFT](http://bigdft.org/) project.

http://bigdft.org/

Installation of python is required. Currently,
only python 2.7 is supported. Future releases of FLAME will
support python 3


### Installing FLAME on Linux

* It is recommended to install FLAME in a different
directory than the source code.

#### Here are steps:

1. First, install futile which is
   a set of utilities from the BigDFT project.
   Preferably, use the version provided with
   FLAME to avoid conflicts.
   Untar the included futile-suite.tar.gz (`tar -zxvf futile-suite.tar.gz`), then 
   create a new build directory (e.g., `mkdir futile-build ; cd futile-build`), and from there run

   - for GNU compilers:

      `path_to_futile_source/Installer.py build futile -c 
      CC=gcc CXX=g++ FC=gfortran F77=gfortran \
      --with-ext-linalg="-llapack -lblas"`

   - for Intel compilers:

      `path_to_futile_source/Installer.py build futile -c FCFLAGS=-O2 \
      --with-ext-linalg="-L$MKLROOT/lib/intel64 \
      -lmkl_rt -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 -liomp5 -lm" \
      CC=icc CXX=icpc FC=ifort F77=ifort`

   - for parallel compilation use the corresponding MPI wrappers:

      `path_to_futile_source/Installer.py build futile -c FCFLAGS=-O2 \
      --with-ext-linalg="-L$MKLROOT/lib/intel64 \
      -lmkl_rt -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 -liomp5 -lm" \
      CC=mpicc CXX=mpicxx FC=mpif90 F77=mpif90`

   follwed by `make build` if necessary.
   Make sure to adapt the library locations and
   the linking flags appropriately.
   After the installation of futile, we need to link it
   with FLAME.
   To display details on the general linking procedure, run:

   `path_to_futile_source/Installer.py link futile`



2. in FLAME source directory run:

   - `autoreconf -fi`

3. make a directory to build flame and then replace `$FUTILE` by the
full path of futile built directory or define it as an environment variable:

   - `export FUTILE=path_to_built_futile_directory`

   - IMPORTANT: exporting FUTILE variable is required, not optional.

4. `path_to_flame_source/configure FC=mpif90 F77=mpif90 CXX=mpicc CC=mpicc 
FCFLAGS="-I$FUTILE/install/include -shared-intel -mcmodel=large  -mkl=sequential" 
CFLAGS=-mcmodel=large "LIBS=-L$FUTILE/install/lib 
-L/opt/intel/composer_xe_2013.2.146/mkl/lib/intel64 -lfutile-1 -lmkl_rt 
-lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 -liomp5 -lm -lyaml -ldl -lrt -cxxlib"`

   Linking to exterior packages:

   - if you want to link with [SPGLIB](https://atztogo.github.io/spglib/) then add

     `--with-spglib SPGLIB_ROOT=path_to_spglib`

   - if you want to link with [BigDFT](http://bigdft.org) PSolver then add

     `--with-bps BDIR=path_to_bigdft_root`

   - if you want to link with [LAMMPS](https://lammps.sandia.gov) then add

     `--with-lammps LAMMPS_ROOT=path_to_lammps_root`

4. `make`

   You should be ready to use FLAME, the executable is at `src/flame`
