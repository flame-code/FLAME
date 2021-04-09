EXTRA_DIST = \
	op2p_check.ref.yaml \
	linalg_check.ref.yaml \
	example.f90 \
	mpi_layer.f90

CHECK_PROGS = linalg_check
if USE_MPI
CHECK_PROGS += op2p_check mpi_check
endif

check_PROGRAMS = $(CHECK_PROGS)

check: regtest report #op2p_check.report.yaml $(check_PROGRAMS)

report:
	@if test $(MAKELEVEL) = 0 ; then $(PYTHON) $(top_builddir)/tests/report.py ; fi

CLEANFILES = *.out.yaml *.report.yaml *.@MODULE_EXT@ report_remarks report

AM_FCFLAGS = -I$(top_builddir)/wrappers -I$(top_builddir)/wrappers/mpi -I$(top_builddir)/flib -I$(top_builddir)/dicts -I$(top_builddir)/src @LIB_YAML_CFLAGS@ @MPI_INCLUDE@ -I. -I$(srcdir)
AM_LDFLAGS = -L$(top_builddir)/src
if BUILD_DYNAMIC_LIBS
AM_LDFLAGS += @FC_RPATH@$(abs_top_builddir)/src @FC_RPATH@$(DESTDIR)$(libdir)
else
if USE_SIMGRID_SHARED_ALLOCS
AM_LDFLAGS += @FC_RPATH@$(abs_top_builddir)/src @FC_RPATH@$(DESTDIR)$(libdir)
endif
endif

#note: this way of iluding LDADD does not force the
#dependency between the programs and the library.
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
	$(PYTHON) $(top_builddir)/tests/fldiff_yaml.py -t $(srcdir)/tols.yaml \
			   -l op2p_check -r $(srcdir)/op2p_check.ref.yaml -d $< -o $@

regtest: $(CHECK_PROGS)
	$(PYTHON) $(top_srcdir)/tests/f_regtest.py --fldiff=$(top_builddir)/tests/fldiff_yaml.py --tols=$(srcdir)/tols.yaml --srcdir=$(srcdir)

define F_REGTEST_INSTRUCTIONS
if USE_MPI
- op2p_check:
   runs:
    - $(run_parallel) ./op2p_check -o [10,10,10,10] -n 100
    - $(run_parallel) ./op2p_check -o [100,100] -n 1000
    - $(run_parallel) ./op2p_check -o [201,201] -n 2000
- mpi_check:
   runs:
    - $(run_parallel) ./mpi_check
endif
- linalg_check:
   runs:
    - OMP_NUM_THREADS=1 ./linalg_check -o 75 -v 10000
    - ./linalg_check -o 75 -v 10000
    - ./linalg_check -o 75 -v 10000 -s Yes
    - ./linalg_check -o 75 -v 10000 -c 2
    - ./linalg_check -o 75 -v 10000 -c 2 -s Yes
endef
export F_REGTEST_INSTRUCTIONS
