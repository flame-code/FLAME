# To be installed:
noinst_LIBRARIES = libfmpi.a
module_DATA = $(modules)

ALLTOALL_INCLUDES = alltoallv-inc.f90 mpi_get_alltoallv-inc.f90

ALLREDUCE_INCLUDES = allreduce-inc.f90 allreduce-arr-inc.f90 allreduce-core-inc.f90 allreduce-decl-inc.f90 allreduce-multi-inc.f90

ALLGATHER_INCLUDES = allgather-inc.f90 control-routines-inc.f90

MAXDIFF_INCLUDES =  maxdiff-arr-inc.f90 maxdiff-inc.f90 \
	maxdiff-decl-inc.f90 maxdiff-end-inc.f90

BCAST_INCLUDES = bcast-inc.f90 bcast-decl-inc.f90 bcast-decl-arr-inc.f90 \
	maxdiff-routines-inc.f90 $(MAXDIFF_INCLUDES)

WIN_INCLUDES = win-create-inc.f90

GET_INCLUDES = get-inc.f90 get-decl-inc.f90 get-end-inc.f90

GATHER_INCLUDES = gather-inc.f90 gather-inner-inc.f90

MPI_INCLUDES = scatter-inc.f90 scatterv-inc.f90 scatterv-decl-inc.f90 \
	$(ALLGATHER_INCLUDES) $(ALLREDUCE_INCLUDES) $(ALLTOALL_INCLUDES) \
	$(WIN_INCLUDES) $(GET_INCLUDES) \
	$(BCAST_INCLUDES) $(GATHER_INCLUDES)

if USE_MPI
if HAVE_MPI3
mpi_source =
else
mpi_source = MPI3fake.c
endif
mpi_include = $(MPI_INCLUDES)
else
mpi_source = MPIfake.f90
mpi_include = $(MPI_INCLUDES) mpif.h
endif

EXTRA_DIST = $(MPI_INCLUDES)

#Fake mpif.h for serial compilation
mpif.h:
	touch mpif.h &&\
	echo "integer, parameter :: MPI_COMM_NULL=2, MPI_COMM_SELF=0, MPI_COMM_WORLD=1" >> mpif.h &&\
	echo "integer, parameter :: MPI_SUCCESS=0, MPI_SUM=1" >> mpif.h &&\
	echo "integer, parameter :: MPI_DOUBLE_PRECISION=1, MPI_REAL=1, MPI_INTEGER=1, MPI_INTEGER8=1" >> mpif.h &&\
	echo "integer, parameter :: MPI_STATUSES_IGNORE=1, MPI_LOGICAL=1" >> mpif.h &&\
	echo "integer, parameter :: MPI_MIN=1, MPI_MAX=1, MPI_LOR=1, MPI_CHARACTER=1, MPI_REAL8=1" >> mpif.h &&\
	echo "integer, parameter :: MPI_MAX_PROCESSOR_NAME=10, MPI_STATUS_SIZE=1, MPI_LAND=1" >> mpif.h &&\
	echo "integer, parameter :: MPI_REQUEST_NULL=1, MPI_STATUS_IGNORE=1" >> mpif.h &&\
	echo "integer, parameter :: MPI_INTEGER_KIND=4, MPI_OFFSET_KIND=8, MPI_ADDRESS_KIND=8" >> mpif.h &&\
	echo "integer, parameter :: MPI_TAG_UB=1, MPI_INFO_NULL=0" >> mpif.h &&\
	echo "integer, parameter :: MPI_MODE_NOPRECEDE=0, MPI_MODE_NOSUCCEED=0, MPI_MODE_APPEND=0" >> mpif.h &&\
	echo "integer, parameter :: MPI_MODE_DELETE_ON_CLOSE=0, MPI_MODE_EXCL=0, MPI_MODE_RDWR=0, MPI_MODE_SEQUENTIAL=0" >> mpif.h &&\
	echo "integer, parameter :: MPI_LOCK_EXCLUSIVE=0, MPI_LOCK_SHARED=0, MPI_LONG=8" >> mpif.h &&\
	echo "integer, parameter :: MPI_ANY_SOURCE=0, MPI_MODE_CREATE=0, MPI_MODE_WRONLY=0, MPI_ANY_TAG=0" >> mpif.h &&\
	echo "integer, parameter :: MPI_MODE_NOPUT=0, MPI_MODE_RDONLY=0, MPI_MODE_NOSTORE=0, MPI_MODE_NOCHECK=0" >> mpif.h &&\
	echo "integer, parameter :: MPI_GROUP_NULL=0, MPI_PROC_NULL=0, MPI_WIN_NULL=0" >> mpif.h &&\
	echo "integer, parameter :: MPI_IN_PLACE=0" >> mpif.h &&\
	echo "real(kind=8), external :: mpi_wtime" >> mpif.h

libfmpi_a_SOURCES = scalars.f90 allreduce.f90 alltoall.f90 onesided.f90 allgather.f90 bcast.f90 gather.f90 sendrecv.f90 \
	mpi.F90 $(mpi_source) test_mpi_wrappers.f90

$(modules): $(libfmpi_a_OBJECTS)
scalars.o: $(mpi_include)
mpi.o: $(mpi_include) alltoall.o allreduce.o allgather.o
test_mpi_wrappers.o: mpi.o
onesided.o:  $(GET_INCLUDES) scalars.o
allreduce.o: $(ALLREDUCE_INCLUDES) scalars.o
alltoall.o: $(ALLTOALL_INCLUDES) allreduce.o onesided.o
allgather.o: $(ALLGATHER_INCLUDES) scalars.o onesided.o bcast.o
bcast.o: $(BCAST_INCLUDES) scalars.o allreduce.o gather.o
gather.o: $(GATHER_INCLUDES) scalars.o
sendrecv.o: scalars.o
#control.o: scalars.o bcast.o allgather.o

AM_FCFLAGS = -I. -I$(top_builddir)/dicts/ -I$(top_builddir)/flib/ @MPI_INCLUDE@

# Module files.
if CAPITALIZE
modules = WRAPPER_MPI.@MODULE_EXT@ \
	MPIF_MODULE.@MODULE_EXT@ \
	TEST_MPI_WRAPPERS.@MODULE_EXT@\
	FMPI_TYPES.@MODULE_EXT@ \
	F_ALLREDUCE.@MODULE_EXT@ \
	F_ONESIDED.@MODULE_EXT@ \
	F_SENDRECV.@MODULE_EXT@ \
	F_ALLGATHER.@MODULE_EXT@ \
	F_BCAST.@MODULE_EXT@ \
	F_GATHER.@MODULE_EXT@ \
	F_ALLTOALL.@MODULE_EXT@
else
modules = wrapper_mpi.@MODULE_EXT@ \
	mpif_module.@MODULE_EXT@ \
	test_mpi_wrappers.@MODULE_EXT@ \
	fmpi_types.@MODULE_EXT@ \
	f_allreduce.@MODULE_EXT@ \
	f_onesided.@MODULE_EXT@ \
	f_sendrecv.@MODULE_EXT@ \
	f_allgather.@MODULE_EXT@ \
	f_bcast.@MODULE_EXT@ \
	f_gather.@MODULE_EXT@ \
	f_alltoall.@MODULE_EXT@
endif
CLEANFILES = *.@MODULE_EXT@ mpif.h



