POTENTIALS += LJ

POTENTIALS += VCBLJ

POTENTIALS += VASP

POTENTIALS += LenoskyTB

POTENTIALS += MPMD

POTENTIALS += QSC
FFLAGS += -cxxlib

POTENTIALS += ANN

POTENTIALS += BigDFT

ifdef PLATO
	POTENTIALS += PLATO 
	LIB_POT += ../packages/plato/plato.a -L/home/ghasemi/Programs/gsl-1.15/lib -lgsl -lm
	PPFLAG += -DHAVE_PLATO
endif

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

