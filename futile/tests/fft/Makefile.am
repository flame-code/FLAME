## Process this file with automake to produce Makefile.in

SHORT_CHECK=
LONG_CHECK= fft_check.report fft_cache.report

OTHERS = fft_check.f90 fft_cache.f90

EXTRA_DIST = fft_check.out.ref fft_cache.out.ref fft_cache.py

check_PROGRAMS = fft_check fft_cache fft_results

fft_check_SOURCES = fft_check.f90 $(top_srcdir)/wrappers/fft/fft3d.f90 $(mpi_source)
fft_check_LDFLAGS =
fft_check_LDADD = @LINALG_LIBS@

fft_cache_SOURCES = fft_cache.f90 $(top_srcdir)/wrappers/fft/fft3d.f90 $(mpi_source)
fft_cache_LDFLAGS =
fft_cache_LDADD = @LINALG_LIBS@

fft_results_SOURCES = fft_results.f90 $(top_srcdir)/wrappers/fft/fft3d.f90 $(mpi_source)
fft_results_LDFLAGS =
fft_results_LDADD = @LIB_FUTILE_LIBS@ @LINALG_LIBS@

CLEANFILES = fft_check.out fft_cache.out fft_cache.dat fft_columns.dat fft_perf.dat \
	fft_cache fft_cache.gnuplot *.report *.@MODULE_EXT@

AM_FCFLAGS = -I$(top_builddir)/includes @LIB_YAML_CFLAGS@ @MPI_INCLUDE@ -I. -I$(srcdir)

fft_cache.o: fft3d.o
fft_check.o: fft3d.o

fft_check.out: fft_check
	$(run_serial) ./fft_check > fft_check.out

fft_cache.dat: fft_cache
	$(run_serial) ./fft_cache 0   >  fft_cache.out 2>&1
	$(run_serial) ./fft_cache 6   >> fft_cache.out 2>&1
	$(run_serial) ./fft_cache 12  >> fft_cache.out 2>&1
	$(run_serial) ./fft_cache 24  >> fft_cache.out 2>&1
	$(run_serial) ./fft_cache 50  >> fft_cache.out 2>&1
	$(run_serial) ./fft_cache 75  >> fft_cache.out 2>&1
	$(run_serial) ./fft_cache 100 >> fft_cache.out 2>&1
	python $(srcdir)/fft_cache.py

fft_check.report: fft_check.out
	python $(pythondir)/fldiff.py --mode=pseudo --discrepancy=1.e-8 \
		fft_check.out $(srcdir)/fft_check.out.ref | tee fft_check.report

fft_cache.report: fft_cache.dat
	python $(pythondir)/fldiff.py --mode=pseudo \
		fft_cache.out $(srcdir)/fft_cache.out.ref | tee fft_cache.report

check: fft_check.report fft_cache.report report

report:
	@if test $(MAKELEVEL) = 0 ; then python $(pythondir)/report.py ; fi

regtest: fft_results
	python $(pythondir)/f_regtest.py --fldiff=$(pythondir)/fldiff_yaml.py --tols=$(srcdir)/tols.yaml --srcdir=$(srcdir)

define F_REGTEST_INSTRUCTIONS
- fft_results:
   runs:
    - ./fft_results > fft_results.out.yaml
endef
export F_REGTEST_INSTRUCTIONS
