include def_com.mk
include def_pot.mk
include def_opt.mk

LIB_MOD = modules/libmodules.a
LIBS = liball.a

INCLUDES = -I../modules

LIB_SRC = src/libsrc.a

ifdef BPS
	LIBS += $(BDIR)/PSolver/src/libPSolver-1.a \
		$(BDIR)/wrappers/libwrappers.a \
		$(BDIR)/flib/src/libflib-1.a \
		$(BDIR)/yaml-0.1.4/src/.libs/libyaml.a
	PRE_PROC += -DHAVE_BPS
	INCLUDES += -I$(BDIR)/includes -I$(BDIR)/PSolver/src
endif

ifdef SPGLIB
	ARGS+= SPGLIB=1
	LIB_SPGLIB =  $(SPGLIB_ROOT)/src/.libs/libsymspg.a
	LIBS+= $(LIB_SPGLIB)
endif
ifdef LAMMPS
	ARGS+= LAMMPS=1
	LAMMPS_SRC = $(LAMMPS_ROOT)/src
	LIB_MPI_STUBS = $(LAMMPS_SRC)/STUBS/libmpi_stubs.a
	LIB_LAMMPS = $(LAMMPS_SRC)/liblammps_serial_intel.a  $(LAMMPS_ROOT)/lib/reax/*.o $(LAMMPS_ROOT)/lib/meam/*.o $(LAMMPS_ROOT)/lib/poems/*.o
	LIBS+= $(LIB_MPI_STUBS) $(LIB_LAMMPS)
endif
ifdef TINKER
	LIB_TINKER = $(TINKER_ROOT)/source/libtinker.a
	LIB_FFTW3 = /usr/lib/libfftw_threads.a /usr/lib/libfftw.a
	ARGS+= TINKER=1
	LIBS+= $(LIB_TINKER) $(LIB_FFTW3)
endif

MINHOCAO = minhocao/libminhocao.a minhocao/parsestring.o minhocao/lenosky_tb/*.o

DIRS += modules
DIRS += src
DIRS += minhocao

all: build/install/lib/libfutile-1.a $(DIRS) liball.a flame
	@echo "POTENTIALS: $(POTENTIALS)"
	@echo "Pre-processing: $(PRE_PROC)"

MKL = '--with-ext-linalg=-L$(MKLPATH) -lmkl_rt -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 -liomp5 -lm'
COMPILERS = 'CXX=$(MY_CXX)' 'FC=$(F90)' 'F77=$(F90)' 'FCLIBS= '
build/install/lib/libfutile-1.a:
	cd build ; ../Installer.py build futile -v -c 'FCFLAGS=-O2' $(MKL) $(COMPILERS)

.NOTPARALLEL: modules

$(DIRS):
	$(MAKE) PRE_PROC="$(PRE_PROC)" FFLAGS="$(FFLAGS)" INCLUDES="$(INCLUDES)" CC=$(CC) F90=$(F90) $(ARGS) -C $@



LIBS_A = $(LIB_SRC) $(LIB_POT) $(LIB_MOD)
#LIBS += $(LIB_LAPACK)
#LIBS += /home/ghasemi/ghasemi/ghasemi/oldsilicon/Coulomb/EE2DP1DF/MMM2D/MMM2D-1.0/Linux/libMMM2D.a

liball.a: $(LIBS_A)
	ar -scru liball.a `ls -1 src/ofiles/*.o modules/ofiles/*.o |grep -v alborz.o`

#FFLAGS := $(filter-out -traceback,$(FFLAGS))
flame: liball.a src/ofiles/alborz.o
	$(F90) $(filter-out -traceback,$(FFLAGS)) -openmp src/ofiles/alborz.o $(MINHOCAO) $(LIBS) $(INCLUDES) $(FUTILE) -o flame

clean:
	rm -f flame liball.a

cleanall: clean
	for pot in $(DIRS); do printf "\n" ; $(MAKE) -C $$pot clean ; done
	find . -not -path "./build/*" -type f -iname '*.o'
	find . -not -path "./build/*" -type f -iname '*.a'
	find . -not -path "./build/*" -type f -iname '*.mod'

cleanall2: clean
	for pot in $(DIRS); do printf "\n" ; $(MAKE) -C $$pot clean ; done
	rm -rf build/buildrc build/install build/libyaml build/futile
	find . -type f -iname '*.o'
	find . -type f -iname '*.a'
	find . -type f -iname '*.mod'

.PHONY : all $(DIRS) clean cleanall
#*****************************************************************************************
