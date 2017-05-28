include ../machines/thismachine.make

##LD     = gfortran 
##FC     = gfortran
##CC     = gcc
##FFLAGS  =  -O2 -Wno-error -fbounds-check -fbacktrace -ffpe-trap=invalid,zero,overflow 
##LD     = ifort  -m64 -I$(MKLROOT)/include
##FC     = ifort  -m64 -I$(MKLROOT)/include
#LD     = mpif90 -O2 ##-m64 -I$(MKLROOT)/include
LD     = mpif90 #-C -g -traceback ##-m64 -I$(MKLROOT)/include
#FC     = mpif90 -O2 ##-m64 -I$(MKLROOT)/include
FC     = mpif90 #-C -g -traceback  ##-m64 -I$(MKLROOT)/include

CC     = icc
LDFLAGS =  -shared-intel -mcmodel=large -I$(MKLROOT)/include -mkl=sequential
FFLAGS  =  -shared-intel -mcmodel=large # -fPIC
##LDFLAGS =  -C -fp-stack-check -check all -g -traceback -heap-arrays -shared-intel -mcmodel=large -I$(MKLROOT)/include -mkl=sequential
##FFLAGS  =  -C -fp-stack-check -check all -g -traceback -heap-arrays -shared-intel -mcmodel=large # -fPIC
CFLAGS  =  -mcmodel=large
##FFLAGS = -O2 -fPIC 
CLINK =  -cxxlib -lm
CXXFLAGS = -O0 #-fPIC
CPPFLAGS =
CXX = icpc   # replace with your C++ compiler
#CFLAGS  =  -O2
#FFLAGS  = -O2
#LDFLAGS  = $(MKLROOT)/lib/em64t/libmkl_solver_lp64_sequential.a -Wl,--start-group  $(MKLROOT)/lib/em64t/libmkl_gf_lp64.a $(MKLROOT)/lib/em64t/libmkl_sequential.a $(MKLROOT)/lib/em64t/libmkl_core.a -Wl,--end-group -lpthread -lm
 #-mkl=sequential 
LDFLAGS = -L$(MKLROOT)/lib/intel64 -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lpthread -lm
#LIBS=
#PRE_PROC=


LTB = lenosky_tb/clsfunc.o  lenosky_tb/minimize.o  lenosky_tb/nrutil.o  lenosky_tb/tblink.o  lenosky_tb/tbnorth.o
PARSER = parsestring.o
MODULES = modules.o module_sqn.o defs_basis.o atmdata.o

OBJ = $(MODULES) $(PARSER)

ifdef SPGLIB
	# Path to spglib
	LIB_SPGLIB =  $(SPGLIB_ROOT)/src/.libs/libsymspg.a
	SPGLIB_EX = $(SPGLIB_ROOT)/example/
	#SPGLIB_OBJ= $(SPGLIB_EX)/spglib_f08.o
	SPGLIB_OBJ= spglib_f08.o
	OBJ+= $(SPGLIB_OBJ) spglib_int.o
	LIBS+= $(LIB_SPGLIB)
	override PRE_PROC+= -DSPGLIB
endif

#print-%  : ; @echo $* = $($*)

ifdef LAMMPS
	# Path to LAMMPS extraction directory
	##LAMMPS_ROOT = /home/maxamsler/Homefolder/lammps-8Jul13
	#LAMMPS_ROOT = /home/maxamsler/Homefolder/lammps-30Oct14
	# Uncomment the line below if using the MPI stubs library
	MPI_STUBS = -I$(LAMMPS_SRC)/STUBS/*
	LAMMPS_SRC = $(LAMMPS_ROOT)/src
	OBJ+= lammps_int.o LAMMPS.o LAMMPS-wrapper.o
	LIBS+= $(LIB_LAMMPS)
	override INCLUDE+= $(MPI_STUBS)
	override PRE_PROC+= -DLAMMPS
endif

ifdef TINKER
	# Path to tinker library: compile energyandforces.f in the tinker source and include it in libtinker.a
	#TINKER_ROOT = /home/maxamsler/Homefolder/tinker
	OBJ+= tinker.o
	FLAG_TINKER = -no-prec-div -fno-omit-frame-pointer -recursive 
	LIBS+= $(LIB_TINKER) $(LIB_FFTW3)
	FFLAGS+= $(FLAG_TINKER)
	override PRE_PROC+= -DTINKER
endif

### ifdef ALBORZ
### 	OBJ+= alborz.o
### 	LIBS+= $(LIB_ALBORZ) $(FUTIL)
### 	PRE_PROC+= -DALBORZ
### endif

#OBJ+= lammps_int.o LAMMPS.o LAMMPS-wrapper.o
OBJ+= alborz.o fsockets.o sockets.o  corerepulsion.o confinement.o ipi.o msock.o abinit.o  lj_voids.o blj.o mlj.o lenosky_tb.o lenosky_tb_lj.o lenosky_meam.o espresso.o vasp.o siesta.o cp2k.o tersoff.o edip.o dftb.o mopac.o parser.o io_poscar.o io_ascii.o interface_code.o spher_harm_mathematica.o\
	md_fixlat.o find_symmetry.o sqnm.o qe_bfgs.o sd.o bfgs.o lib_lbfgs.o euler.o envelope.o niggli.o fingerprint_gom.o fingerprint_molgom.o fingerprint_oganov.o fingerprint_oganov_cont.o fingerprint_bcm.o fingerprint_atorb.o fingerprint_gaussmol.o fingerprint_xyz2sm.o MHMoF.o


EXEC2 = vasp_recompute_kpt.o 
EXEC3 = expand_poslows.o
EXEC4 = hull.o envelope.o
EXEC5 = vasp_recompute_kpt_odd.o 
EXEC6 = vasp_recompute_cell.o 
EXEC7 = binaries.o 
EXEC8 = ascii2POSCAR.o
EXEC9 = POSCAR2ascii.o 
EXEC10 = recompute_kpt.o 
EXEC11 = espresso_restruct.o
EXEC12 = ternaries.o



all: vasp_recompute_kpt.x vasp_recompute_kpt_odd.x vasp_recompute_cell.x expand_poslows.x \
	hull.x binaries.x ascii2POSCAR.x POSCAR2ascii.x recompute_kpt.x espresso_restruct.x ternaries.x
	
vasp_recompute_kpt.x: $(MODULES) $(OBJ) libminhocao.a $(EXEC2)
	$(LD) $(FFLAGS)  -o vasp_recompute_kpt.x $(EXEC2) $(MODULES) $(PARSER) $(LDFLAGS)

expand_poslows.x: $(MODULES) $(OBJ) libminhocao.a $(EXEC3)
	$(LD) $(FFLAGS)  -o expand_poslows.x $(EXEC3) $(LDFLAGS)

hull.x: $(MODULES) $(OBJ) libminhocao.a $(EXEC4)
	$(LD) $(FFLAGS)  -o hull.x $(EXEC4) $(MODULES) $(PARSER) $(LDFLAGS) 

vasp_recompute_kpt_odd.x: $(MODULES) $(OBJ) libminhocao.a $(EXEC5)
	$(LD) $(FFLAGS)  -o vasp_recompute_kpt_odd.x $(EXEC5) $(MODULES) $(PARSER) $(LDFLAGS)

vasp_recompute_cell.x: $(MODULES) $(OBJ) libminhocao.a $(EXEC6)
	$(LD) $(FFLAGS)  -o vasp_recompute_cell.x $(EXEC6) $(LDFLAGS)

binaries.x: $(MODULES) $(OBJ) libminhocao.a $(EXEC7)
	$(LD) $(FFLAGS)  -o binaries.x $(EXEC7) $(LDFLAGS)

ascii2POSCAR.x: $(MODULES) $(OBJ) libminhocao.a $(EXEC8)
	$(LD) $(FFLAGS)  -o ascii2POSCAR.x $(EXEC8) $(LDFLAGS)

POSCAR2ascii.x: $(MODULES) $(OBJ) libminhocao.a $(EXEC9)
	$(LD) $(FFLAGS)  -o POSCAR2ascii.x $(EXEC9) $(LDFLAGS)

recompute_kpt.x: $(MODULES) $(OBJ) libminhocao.a $(EXEC10)
	$(LD) $(FFLAGS)  -o recompute_kpt.x $(EXEC10) $(MODULES) $(PARSER) $(LDFLAGS)

espresso_restruct.x: $(MODULES) $(OBJ) libminhocao.a $(EXEC11)
	$(LD) $(FFLAGS)  -o espresso_restruct.x $(EXEC11) $(MODULES) $(PARSER) $(LDFLAGS)

ternaries.x: $(MODULES) $(OBJ) libminhocao.a $(EXEC12)
	$(LD) $(FFLAGS)  -o ternaries.x $(EXEC12) $(LDFLAGS)

lenosky_tn%.o: %.c
	$(CC) $(CFLAGS) -c $<
sockets.o: sockets.c
	$(CC) $(FLAGS) -c -o sockets.o sockets.c
%.o: %.F90
	$(FC) $(FFLAGS) $(PRE_PROC) -c $<

#MHMoF.o: $(MODULES) MHMoF.F90

libminhocao.a: $(OBJ) $(LTB)
	ar -scru libminhocao.a $(OBJ)

clean:
	rm -rf $(OBJ) $(LTB) libminhocao.a MHMoF.o vasp_recompute_kpt.o *.mod *~
	rm -rf vasp_recompute_kpt.x vasp_recompute_kpt_odd.x vasp_recompute_cell.x expand_poslows.x \
		hull.x binaries.x ascii2POSCAR.x POSCAR2ascii.x recompute_kpt.x espresso_restruct.x ternaries.x
	rm -rf  *.o *.mod liblammps_fortran.a liblammps_fortran.so  liblammps_serial_intel.a  libmpi_stubs.a 
#	touch $(SPGLIB_EX)/spglib_f08.f90


$(SPGLIB_OBJ): $(SPGLIB_EX)/spglib_f08.f90
	$(FC) $(FFLAGS) -c $< -o $@
spglib_f08.mod: $(SPGLIB_EX)/spglib_f08.f90 $(SPGLIB_OBJ)
	@true

#MODULE: $(MODULES) $(SPGLIB_OBJ)
MODULE: $(MODULES) $(SPGLIB_OBJ) spglib_f08.mod 
	@true



# Flags for Fortran compiler, C++ compiler, and C preprocessor, respectively

lammps_int.o : lammps_int.F90 LAMMPS.o LAMMPS-wrapper.o
	$(FC) $(FFLAGS) -c $<


LAMMPS.o lammps.mod : LAMMPS.F90
	$(FC) $(CPPFLAGS) $(FFLAGS)  -c $<

LAMMPS-wrapper.o : LAMMPS-wrapper.cpp LAMMPS-wrapper.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS)  -c $< -I$(LAMMPS_SRC)    $(MPI_STUBS)




###mpif90 -m64 -I/opt/intel/composer_xe_2011_sp1.6.233/mkl/include -O2 -fPIC  -o MHMoF.x MHMoF.o modules.o defs_basis.o atmdata.o  confinement.o abinit.o lammps_int.o blj.o lenosky_tb.o vasp.o siesta.o cp2k.o dftb.o mopac.o interface_code.o find_symmetry.o bfgs.o lib_lbfgs.o euler.o niggli.o fingerprint_oganov.o fingerprint_gaussmol.o LAMMPS-wrapper.o LAMMPS.o lenosky_tb/clsfunc.o  lenosky_tb/minimize.o  lenosky_tb/nrutil.o  lenosky_tb/tblink.o  lenosky_tb/tbnorth.o liblammps_serial_intel.a libmpi_stubs.a  -cxxlib -lm -I$MKLROOT/include -mkl=sequential
###mpif90 -O2 -fPIC   -o MHMoF.x MHMoF.o lenosky_tb/clsfunc.o  lenosky_tb/minimize.o  lenosky_tb/nrutil.o  lenosky_tb/tblink.o  lenosky_tb/tbnorth.o modules.o defs_basis.o atmdata.o  confinement.o abinit.o lammps_int.o blj.o lenosky_tb.o vasp.o siesta.o cp2k.o dftb.o mopac.o interface_code.o find_symmetry.o bfgs.o lib_lbfgs.o euler.o niggli.o fingerprint_oganov.o fingerprint_gaussmol.o LAMMPS.o LAMMPS-wrapper.o  -L/opt/intel/composer_xe_2011_sp1.6.233/mkl/lib/intel64 -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -I/home/maxamsler/Homefolder/lammps-8Jul13/src -I/home/maxamsler/Homefolder/lammps-8Jul13/src/STUBS/* -cxxlib -lm

