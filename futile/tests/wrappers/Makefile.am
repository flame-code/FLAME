EXTRA_DIST = \
	op2p_check.ref.yaml \
	example.f90 \
	mpi_check.f90 \
	mpi_layer.f90 \
	op2p_check.f90

check_PROGRAMS = op2p_check

check: op2p_check.report.yaml $(check_PROGRAMS)

report:
	@if test $(MAKELEVEL) = 0 ; then python $(top_builddir)/tests/report.py ; fi

EXTRA_PROGRAMS = mpi_check

CLEANFILES = op2p_check.out.yaml \
	mpi_layer.@MODULE_EXT@ MPI_LAYER.@MODULE_EXT@ \
	overlap_point_to_point.@MODULE_EXT@ OVERLAP_POINT_TO_POINT.@MODULE_EXT@

AM_FCFLAGS = -I$(top_builddir)/includes @LIB_YAML_CFLAGS@ @MPI_INCLUDE@ -I. -I$(srcdir)
AM_LDFLAGS = -L$(top_builddir)/src
if BUILD_DYNAMIC_LIBS
AM_LDFLAGS += -Wl,-rpath=$(abs_top_builddir)/src -Wl,-rpath=$(DESTDIR)$(libdir)
endif

op2p_check_SOURCES = op2p_check.f90
op2p_check_LDADD = $(bigdft_library) @LIB_FUTILE_LIBS@

mpi_check_SOURCES = mpi_check.f90
mpi_check_LDADD = $(bigdft_library)

op2p_check.o: $(top_builddir)/wrappers/libwrappers.a

op2p_check.out.yaml: op2p_check
	$(run_parallel) ./op2p_check -o [10,10,10,10] -n 100 > $@
	$(run_parallel) ./op2p_check -o [100,100] -n 1000 >> $@
	$(run_parallel) ./op2p_check -o [201,201] -n 2000 >> $@

op2p_check.report.yaml: op2p_check.out.yaml
	python $(pythondir)/fldiff_yaml.py -t $(srcdir)/tols.yaml \
			   -l op2p_check -r $(srcdir)/op2p_check.ref.yaml -d $< -o $@
