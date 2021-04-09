## Process this file with automake to produce Makefile.in
EXTRA_DIST = \
	tols.yaml \
	yaml_test.ref.yaml \
	yaml_argparse.ref.yaml \
	f_allocations.ref.yaml \
	hooks.ref.yaml \
	enum_test.ref.yaml

dynlibdir = $(libdir)

check_PROGRAMS = yaml_test hooks yaml_argparse f_allocations test_plugin enum_test
noinst_LIBRARIES = libpong.a

AM_FCFLAGS = -I$(top_builddir)/dicts/ -I../../flib -I../../src @FC_PIC_FLAG@

AM_LDFLAGS = -L$(top_builddir)/src
#LIBS += -lflib -lrt @LIB_YAML_LIBS@

YAML_DATABASE= database.yaml


if BUILD_DYNAMIC_LIBS
flib_libraries = -lfutile-1
AM_LDFLAGS += @FC_RPATH@$(abs_top_builddir)/src
else
if USE_SIMGRID_SHARED_ALLOCS
AM_LDFLAGS += @FC_RPATH@$(abs_top_builddir)/src @FC_RPATH@$(DESTDIR)$(libdir)
endif
flib_libraries = $(top_builddir)/src/libfutile-1.a
LIBS += @LIBS_EXT_DEPS@
endif

yaml_test_SOURCES = yaml_test.f90 \
	errs.f90 \
	yamlout.f90 \
	dicts.f90 \
	dynmem.f90 \
	yaml_invoice.f90 \
	utls.f90
nodist_yaml_test_SOURCES = $(DATABASE_SOURCES)

yaml_test_LDADD = $(flib_libraries)

enum_test_SOURCES = enum.f90
enum_test_LDADD = $(flib_libraries)

yaml_argparse_SOURCES = yaml_argparse.f90
yaml_argparse_LDADD = $(flib_libraries)

f_allocations_SOURCES = f_allocations.f90
f_allocations_LDADD = $(flib_libraries)


hooks_SOURCES = hooks.f90
hooks_LDADD = $(flib_libraries)

libpong_a_SOURCES = pong.f90
test_plugin_SOURCES = test_plugin.f90
test_plugin_LDADD = $(flib_libraries) @FC_EXPORTS@

CLEANFILES = *.out.yaml *.report.yaml *.@MODULE_EXT@ memstatus.yaml report_remarks *.err libpong.so citations.bib libfmalloc-1.so

include $(srcdir)/import_database.mk

yaml_test.o: $(abs_top_builddir)/src/libfutile-1.a
hooks.o: $(abs_top_builddir)/src/libfutile-1.a

libpong.so: $(libpong_a_OBJECTS)
	$(FC) $(FCFLAGS) @FC_BUILD_SHARED@ $(LDFLAGS) @FC_RPATH@$(DESTDIR)$(dynlibdir) -o $@ $<
	@chmod a+x $@

yaml_test.out.yaml: yaml_test
	./yaml_test > yaml_test.out.yaml

yaml_test.report.yaml: tols.yaml yaml_test.ref.yaml yaml_test.out.yaml
	$(PYTHON) $(top_builddir)/tests/fldiff_yaml.py --label yaml_test -t $(srcdir)/tols.yaml 	\
		-r $(srcdir)/yaml_test.ref.yaml -d yaml_test.out.yaml --output $@

yaml_argparse.out.yaml: yaml_argparse
	./yaml_argparse > yaml_argparse.out.yaml
	./yaml_argparse -n 100 >> yaml_argparse.out.yaml
	./yaml_argparse -n 100 -d [89,90,90] >> yaml_argparse.out.yaml

yaml_argparse.report.yaml: tols.yaml yaml_argparse.ref.yaml yaml_argparse.out.yaml
	$(PYTHON) $(top_builddir)/tests/fldiff_yaml.py --label yaml_argparse -t $(srcdir)/tols.yaml 	\
		-r $(srcdir)/yaml_argparse.ref.yaml -d yaml_argparse.out.yaml --output $@

hooks.report.yaml: tols.yaml hooks.ref.yaml
	$(PYTHON) $(top_builddir)/tests/fldiff_yaml.py --label hooks -t $(srcdir)/tols.yaml 	\
		-r $(srcdir)/hooks.ref.yaml -d hooks.out.yaml --output $@

libfmalloc-1.so:
	ln -sf $(top_builddir)/flib/$@

check: yaml_test.report.yaml regtest report

report:
	@if test $(MAKELEVEL) = 0 ; then $(PYTHON) $(top_builddir)/tests/report.py ; fi

regtest: yaml_argparse test_plugin hooks libpong.so f_allocations libfmalloc-1.so
	$(PYTHON) $(top_srcdir)/tests/f_regtest.py --fldiff=$(top_builddir)/tests/fldiff_yaml.py \
	       --tols=$(srcdir)/tols.yaml --srcdir=$(srcdir) \
	       --exclusive=$(CHECK_ONLY)

define F_REGTEST_INSTRUCTIONS
- yaml_argparse:
   binary: yaml_argparse
   runs:
   - ./yaml_argparse
   - ./yaml_argparse -n 100
   - ./yaml_argparse -n 100 -d [89,90,90]
   reference: $(srcdir)/yaml_argparse.ref.yaml
   output: yaml_argparse.out.yaml
- f_allocations:
   runs:
   - ./f_allocations
   - $(top_srcdir)/src/f_memcheck ./f_allocations 2>f_allocations.err
- hooks:
   runs:
   - ./hooks
- enum_test:
   runs:
   - ./enum_test
endef
export F_REGTEST_INSTRUCTIONS
