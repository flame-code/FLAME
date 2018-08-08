noinst_LIBRARIES = libmodules.a
.NOTPARALLEL:

libmodules_a_SOURCES =  task_mod.F90 unitsconversion_mod.F90 processors_mod.F90 atoms_mod.F90 \
		linked_lists_mod.F90 \
		constants_mod.F90 \
		spline_mod.F90 \
		electrostatics_mod.F90 \
		ann_mod.F90 \
		shortrange_mod.F90 \
		potential_mod.F90 opt_mod.F90 genconf_mod.F90 parini_mod.F90 \
		tightbinding_mod.F90 \
		bader_mod.F90 \
		minhocao_mod.F90 \
		interface_mod.F90 minhopp_mod.F90 \
		saddle_mod.F90 dynamics_mod.F90 \
		fsockets.F90 \
		alborz_as_potential_mod.F90 \
		timing_mod.F90 \
		constants_minhocao_mod.F90

#temporary compiling line for gfortran
AM_FCFLAGS = -I. -O
AM_CFLAGS=-I$(top_builddir)

if HAVE_BPS
AM_FCFLAGS += -I$(BDIR)/install/include
AM_FCFLAGS += -DHAVE_BPS
endif

CLEANFILES = *.mod *.MOD


