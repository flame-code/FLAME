# To be installed:
noinst_LIBRARIES = libwrappers.a
module_DATA = $(modules)

SUBDIRS = mpi

# Source code.
if USE_DGEMMSY
DGEMMSY_SOURCES = interface_dgemmsy.f90 \
	dgemmsy/dgemmsy.c       dgemmsy/dgemmsy.h \
	dgemmsy/gemm_block_c.c  dgemmsy/gemm_block_c.h \
	dgemmsy/patterns.c      dgemmsy/patterns.h \
	dgemmsy/dgemmsy_utils.c dgemmsy/utils.h \
	dgemmsy/visitors.c      dgemmsy/visitors.h
else
DGEMMSY_SOURCES = interface_dgemmsy_fake.f90
endif
if USE_CUDA_GPU
CUBLAS_SOURCES = fortran.c
else
CUBLAS_SOURCES = cublas_fake.f90
endif
if USE_OCL
OCL_SOURCES =
else
OCL_SOURCES = interface_ocl_fake.f90
endif

# Hack for MIC compilation.
if USE_MIC
cublas_fake.o: cublas_fake.f90
	$(FCCOMPILE) -offload-attribute-target=mic -c -o $@ $<
endif

libwrappers_a_SOURCES = $(DGEMMSY_SOURCES) $(CUBLAS_SOURCES) $(OCL_SOURCES) \
	fft/fft2d.f90 \
	fft/fft3d.f90 \
	linalg.f90 \
	f_blas.f90 \
	op2p_module.f90

# Dependencies
$(modules): $(libwrappers_a_OBJECTS)
interface_cuda_fake.o: linalg.o
interface_dgemmsy_fake.o: linalg.o
f_blas.o: linalg.o
op2p_module.o: linalg.o
#fft2d.o: fft3d.o

AM_FCFLAGS = -I. -I$(top_builddir)/dicts/ -I$(top_builddir)/flib/ -I$(builddir)/mpi #@MPI_INCLUDE@
AM_CPPFLAGS = @DGEMMSY_CPPFLAGS@
if USE_CUDA_GPU
AM_CPPFLAGS += -I@CUDA_PATH@/include
endif

# Module files.
if CAPITALIZE
modules = MODULE_FFT_SG.@MODULE_EXT@ \
	WRAPPER_LINALG.@MODULE_EXT@ \
	F_BLAS.@MODULE_EXT@ \
	OVERLAP_POINT_TO_POINT.@MODULE_EXT@
else
modules = module_fft_sg.@MODULE_EXT@ \
	wrapper_linalg.@MODULE_EXT@ \
	f_blas.@MODULE_EXT@ \
	overlap_point_to_point.@MODULE_EXT@
endif
CLEANFILES = *.@MODULE_EXT@
