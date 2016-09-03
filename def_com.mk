include machines/thismachine.make
F90 = $(MY_FC)
ifdef NOMPI
	PRE_PROC+=
else
	PRE_PROC+= -DMPI
endif

ifdef DEBUG
	FFLAGS = -C -g -traceback -O0 -ftrapuv
else
	FFLAGS = $(MY_FFLAGS)
endif
LIB_LAPACK = $(MY_LIBS)
CC = icc
