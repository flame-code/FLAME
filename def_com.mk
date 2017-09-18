include machines/thismachine.make
F90 = $(MY_FC)
ifdef NOMPI
	PRE_PROC+=
else
	PRE_PROC+= -DMPI
endif

ifdef DEBUG
	#FFLAGS = -C -g -traceback -O0 -ftrapuv  -shared-intel -mcmodel=large  -mkl=sequential
	FFLAGS = $(MY_FFLAGS_DEBUG)  -shared-intel -mcmodel=large  -mkl=sequential
else
	FFLAGS = $(MY_FFLAGS)  -shared-intel -mcmodel=large  -mkl=sequential
endif
LIB_LAPACK = $(MY_LIBS)
CC = icc
