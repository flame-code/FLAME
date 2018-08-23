POTENTIALS += LJ

POTENTIALS += VCBLJ

POTENTIALS += VASP

POTENTIALS += LenoskyTB

POTENTIALS += MPMD

POTENTIALS += QSC
FFLAGS += -cxxlib

POTENTIALS += ANN

POTENTIALS += BigDFT

ifdef SIESTA
	POTENTIALS += SIESTA
	LIB_POT = ../packages/siesta/siesta.a \
		`../packages/siesta/siesta-3.1/Obj-serial/FoX/FoX-config --libs --wcml`
	PPFLAG += -DHAVE_SIESTA
	INC_POT = potentials/SIESTA/modfiles
endif

OPTIONS =

ifdef BPS
    OPTIONS += BPS
endif

ifdef LAMMPS
	# Path to LAMMPS extraction directory
	##LAMMPS_ROOT = /home/maxamsler/Homefolder/lammps-8Jul13
	#LAMMPS_ROOT = /home/maxamsler/Homefolder/lammps-30Oct14
	# Uncomment the line below if using the MPI stubs library
	MPI_STUBS = -I$(LAMMPS_SRC)/STUBS/*
	LAMMPS_SRC = $(LAMMPS_ROOT)/src
	OBJ+= lammps_int.o LAMMPS.o LAMMPS-wrapper.o
	LIBS+= $(LIB_LAMMPS)
	INCLUDE_C= -I$(MY_MPI_INCLUDE) -I$(LAMMPS_ROOT)/src
	override INCLUDE+= $(MPI_STUBS)
	#override PRE_PROC+= -DLAMMPS
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
