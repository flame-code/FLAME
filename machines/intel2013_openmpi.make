#compilers
MY_FC= mpif90
MY_CC= mpicc
MY_CXX= mpicxx
#MY_FC= ifort
#MY_CC= icc
#MY_CXX= icpc

#compiler flags
MY_FFLAGS=-O2 
MY_FFLAGS_DEBUG= -g -O0

BDIR = /home/ghasemi/build/bigdft/rev599

#includes
MY_MPI_INCLUDE=/home/ghasemi/Programs/openmpi-1.6.5/include

#libraries
MY_MKLROOT=/opt/intel/composer_xe_2011_sp1.7.256/mkl

MY_LAPACK_scalapack= $(MKLROOT)/lib/intel64/libmkl_scalapack_lp64.a \
    -Wl,--start-group  $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a \
    $(MKLROOT)/lib/intel64/libmkl_sequential.a $(MKLROOT)/lib/intel64/libmkl_core.a \
    $(MKLROOT)/lib/intel64/libmkl_blacs_openmpi_lp64.a -Wl,--end-group -lpthread -lm

MY_LAPACK_normal= -Wl,--start-group  $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a \
    $(MKLROOT)/lib/intel64/libmkl_sequential.a \
    $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm

FUTILE = -Ibuild/install/include -Lbuild/install/lib \
    -L/opt/intel/composer_xe_2013.2.146/mkl/lib/intel64 -lfutile-1 \
    -lmkl_rt -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 -liomp5 -lm -lyaml -lrt

MY_FFTW=
#MY_FFTW=/usr/local/lib/libfftw3.a
MY_LIBS= $(MY_FFTW) $(MY_LAPACK_normal) 
MY_LIBS_LACPACK= $(MY_FFTW) $(MY_LAPACK_scalapack) 
