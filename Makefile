include def_com.mk
include def_pot.mk
include def_opt.mk

LIB_MOD = modules/libmodules.a
LIBS = liball.a

INCLUDES = -I../modules -I../build/install/include

LIB_SRC = src/libsrc.a

ifdef BPS
	LIBS += -L$(BDIR)/install/lib/ -lPSolver-1
	PRE_PROC += -DHAVE_BPS
	INCLUDES += -I$(BDIR)/install/include
endif

ifdef SPGLIB
	ARGS+= SPGLIB=1
	LIB_SPGLIB =  $(SPGLIB_ROOT)/src/.libs/libsymspg.a
	LIBS+= $(LIB_SPGLIB)
endif
ifdef LAMMPS
	ARGS+= LAMMPS=1
	LAMMPS_SRC = $(LAMMPS_ROOT)/src
	#LIB_MPI_STUBS = $(LAMMPS_SRC)/STUBS/libmpi_stubs.a
	#LIB_LAMMPS = $(LAMMPS_SRC)/liblammps_serial_intel.a  $(LAMMPS_ROOT)/lib/reax/*.o $(LAMMPS_ROOT)/lib/meam/*.o $(LAMMPS_ROOT)/lib/poems/*.o
	LIB_LAMMPS = $(LAMMPS_SRC)/liblammps_mpi.a
	#LIBS+= $(LIB_MPI_STUBS) $(LIB_LAMMPS) do not use, it duplicates MPI libraries
	LIBS+= $(LIB_LAMMPS)
	PRE_PROC += -DHAVE_LAMMPS
endif
ifdef TINKER
	LIB_TINKER = $(TINKER_ROOT)/source/libtinker.a
	LIB_FFTW3 = /usr/lib/libfftw_threads.a /usr/lib/libfftw.a
	ARGS+= TINKER=1
	LIBS+= $(LIB_TINKER) $(LIB_FFTW3)
endif

#DIR  +=MINHOCAO = minhocao/libminhocao.a minhocao/lenosky_tb/*.o

DIRS += modules
DIRS += src
#DIRS += minhocao

all: makedirs build/install/lib/libfutile-1.a $(DIRS) liball.a flame\
	vasp_recompute_kpt.x\
	vasp_recompute_kpt_odd.x vasp_recompute_cell.x expand_poslows.x \
	convex_hull.x binaries.x ascii2POSCAR.x POSCAR2ascii.x recompute_kpt.x\
	espresso_restruct.x ternaries.x
	@echo "POTENTIALS: $(POTENTIALS)"
	@echo "Pre-processing: $(PRE_PROC)"

makedirs:
	mkdir -p build
	mkdir -p docs/_build
	mkdir -p docs/_static
	mkdir -p modules/ofiles
	mkdir -p src/ofiles/lenosky_tb

#MKL = '--with-ext-linalg=-L$(MKLPATH) -lmkl_rt -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 -liomp5 -lm'
MKL = '--with-ext-linalg=-L${MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_blacs_openmpi_lp64 -lpthread -lm -ldl'
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
	ar -scru liball.a `ls -1 src/ofiles/*.o  src/ofiles/lenosky_tb/*.o modules/ofiles/*.o  |grep -v \<alborz.o\> | grep -v vasp_recompute_kpt.o | grep -v expand_poslows.o | grep -v convex_hull.o | grep -v vasp_recompute_kpt_odd.o | grep -v vasp_recompute_cell.o | grep -v binaries.o | grep -v ascii2POSCAR.o | grep -v POSCAR2ascii.o | grep -v recompute_kpt.o | grep -v PWSCF_restruct.o | grep -v ternaries.o`

#FFLAGS := $(filter-out -traceback,$(FFLAGS))
flame: liball.a src/ofiles/alborz.o
	$(F90) $(filter-out -traceback,$(FFLAGS)) -openmp src/ofiles/alborz.o $(MINHOCAO) $(LIBS) $(INCLUDES) $(FUTILE) -o flame

OBJDIR = src/ofiles
PARSER = $(OBJDIR)/parser_core_minhocao.o  modules/ofiles/minhocao_mod.o
EXEC2 = $(OBJDIR)/vasp_recompute_kpt.o 
EXEC3 = $(OBJDIR)/expand_poslows.o
EXEC4 = $(OBJDIR)/convex_hull.o
EXEC5 = $(OBJDIR)/vasp_recompute_kpt_odd.o 
EXEC6 = $(OBJDIR)/vasp_recompute_cell.o 
EXEC7 = $(OBJDIR)/binaries.o 
EXEC8 = $(OBJDIR)/ascii2POSCAR.o
EXEC9 = $(OBJDIR)/POSCAR2ascii.o 
EXEC10 = $(OBJDIR)/recompute_kpt.o 
EXEC11 = $(OBJDIR)/PWSCF_restruct.o
EXEC12 = $(OBJDIR)/ternaries.o
vasp_recompute_kpt.x: $(OBJ_MINHOCAO)  $(EXEC2)
	$(F90) $(FFLAGS)  -o vasp_recompute_kpt.x $(EXEC2) $(INCLUDES) $(PARSER) 

expand_poslows.x: $(OBJ_MINHOCAO) $(EXEC3)
	$(F90) $(FFLAGS)  -o expand_poslows.x $(EXEC3)  $(OBJ_MINHOCAO)

convex_hull.x: $(OBJ_MINHOCAO) $(EXEC4)
	$(F90) $(FFLAGS)  -o convex_hull.x $(EXEC4) $(INCLUDES) $(PARSER) src/ofiles/atoms_minhocao.o src/ofiles/envelope.o

vasp_recompute_kpt_odd.x: $(OBJ_MINHOCAO) $(EXEC5)
	$(F90) $(FFLAGS)  -o vasp_recompute_kpt_odd.x $(EXEC5) $(INCLUDES) $(PARSER) 

vasp_recompute_cell.x: $(OBJ_MINHOCAO)  $(EXEC6)
	$(F90) $(FFLAGS)  -o vasp_recompute_cell.x $(EXEC6)

binaries.x: $(OBJ_MINHOCAO)  $(EXEC7)
	$(F90) $(FFLAGS)  -o binaries.x $(EXEC7) 

ascii2POSCAR.x: $(OBJ_MINHOCAO)  $(EXEC8)
	$(F90) $(FFLAGS)  -o ascii2POSCAR.x $(EXEC8) 

POSCAR2ascii.x: $(OBJ_MINHOCAO) $(EXEC9)
	$(F90) $(FFLAGS)  -o POSCAR2ascii.x $(EXEC9) 

recompute_kpt.x: $(OBJ_MINHOCAO)  $(EXEC10)
	$(F90) $(FFLAGS)  -o recompute_kpt.x $(EXEC10) $(INCLUDES) $(PARSER) 

espresso_restruct.x: $(OBJ_MINHOCAO)  $(EXEC11)
	$(F90) $(FFLAGS)  -o espresso_restruct.x $(EXEC11) $(INCLUDES) $(PARSER) 

ternaries.x: $(OBJ_MINHOCAO)  $(EXEC12)
	$(F90) $(FFLAGS)  -o ternaries.x $(EXEC12) 

interface: 
	./utils/python/build_mod_interface.py
clean:
	rm -f flame liball.a vasp_recompute_kpt.x expand_poslows.x \
	convex_hull.x vasp_recompute_kpt_odd.x vasp_recompute_cell.x \
	binaries.x ascii2POSCAR.x POSCAR2ascii.x recompute_kpt.x \
	espresso_restruct.x ternaries.x

cleanall: clean
	for pot in $(DIRS); do printf "\n" ; $(MAKE) -C $$pot clean ; done
	find . -not -path "./build/*" -type f -iname '*.o'
	find . -not -path "./build/*" -type f -iname '*.a'
	find . -not -path "./build/*" -type f -iname '*.mod'

veryclean: clean
	for pot in $(DIRS); do printf "\n" ; $(MAKE) -C $$pot clean ; done
	rm -rf build/buildrc build/install build/libyaml build/futile
	find . -type f -iname '*.o'
	find . -type f -iname '*.a'
	find . -type f -iname '*.mod'

.PHONY : all $(DIRS) clean cleanall
#*****************************************************************************************
