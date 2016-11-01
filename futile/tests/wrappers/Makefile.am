EXTRA_DIST = \
	op2p_check.ref.yaml \
	example.f90 \
	mpi_check.f90 \
	mpi_layer.f90 \
	op2p_check.f90 \
	linalg_check.f90

check_PROGRAMS = op2p_check linalg_check

check: regtest report #op2p_check.report.yaml $(check_PROGRAMS)

report:
	@if test $(MAKELEVEL) = 0 ; then python $(top_builddir)/tests/report.py ; fi

EXTRA_PROGRAMS = mpi_check

CLEANFILES = *.out.yaml *.report.yaml *.@MODULE_EXT@ report_remarks

AM_FCFLAGS = -I$(top_builddir)/includes @LIB_YAML_CFLAGS@ @MPI_INCLUDE@ -I. -I$(srcdir)
AM_LDFLAGS = -L$(top_builddir)/src
if BUILD_DYNAMIC_LIBS
AM_LDFLAGS += -Wl,-rpath=$(abs_top_builddir)/src -Wl,-rpath=$(DESTDIR)$(libdir)
else
if USE_SIMGRID_SHARED_ALLOCS
AM_LDFLAGS += -Wl,-rpath=$(abs_top_builddir)/src -Wl,-rpath=$(DESTDIR)$(libdir)
endif
endif

op2p_check_SOURCES = op2p_check.f90
op2p_check_LDADD = @LIB_FUTILE_LIBS@

linalg_check_SOURCES = linalg_check.f90
linalg_check_LDADD = @LIB_FUTILE_LIBS@


mpi_check_SOURCES = mpi_check.f90
mpi_check_LDADD =  @LIB_FUTILE_LIBS@

op2p_check.o linalg_check.o: $(top_builddir)/wrappers/libwrappers.a

op2p_check.out.yaml: op2p_check
	$(run_parallel) ./op2p_check -o [10,10,10,10] -n 100 > $@
	$(run_parallel) ./op2p_check -o [100,100] -n 1000 >> $@
	$(run_parallel) ./op2p_check -o [201,201] -n 2000 >> $@

op2p_check.report.yaml: op2p_check.out.yaml
	python $(pythondir)/fldiff_yaml.py -t $(srcdir)/tols.yaml \
			   -l op2p_check -r $(srcdir)/op2p_check.ref.yaml -d $< -o $@

regtest: op2p_check linalg_check
	python $(pythondir)/f_regtest.py --fldiff=$(pythondir)/fldiff_yaml.py --tols=$(srcdir)/tols.yaml --srcdir=$(srcdir)

define F_REGTEST_INSTRUCTIONS
- op2p_check:
   runs:
    - $(run_parallel) ./op2p_check -o [10,10,10,10] -n 100 > op2p_check.out.yaml
    - $(run_parallel) ./op2p_check -o [100,100] -n 1000 >> op2p_check.out.yaml
    - $(run_parallel) ./op2p_check -o [201,201] -n 2000 >> op2p_check.out.yaml
- linalg_check:
   runs:
    - OMP_NUM_THREADS=1 ./linalg_check -o 75 -v 10000 > linalg_check.out.yaml
    - ./linalg_check -o 75 -v 10000 >> linalg_check.out.yaml
    - ./linalg_check -o 75 -v 10000 -s Yes >> linalg_check.out.yaml
    - ./linalg_check -o 75 -v 10000 -c 2 >> linalg_check.out.yaml
    - ./linalg_check -o 75 -v 10000 -c 2 -s Yes >> linalg_check.out.yaml
endef
export F_REGTEST_INSTRUCTIONS
