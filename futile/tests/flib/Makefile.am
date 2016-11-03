## Process this file with automake to produce Makefile.in
EXTRA_DIST = \
	tols.yaml \
	yaml_test.ref.yaml \
	yaml_argparse.ref.yaml \
	f_allocations.ref.yaml \
	hooks.ref.yaml

dynlibdir = $(libdir)

check_PROGRAMS = yaml_test hooks yaml_argparse f_allocations
if BUILD_DYNAMIC_LIBS
check_PROGRAMS += test_plugin
noinst_LIBRARIES = libpong.a
endif

AM_FCFLAGS = -I../../flib -I../../src
AM_LDFLAGS = -L$(top_builddir)/src
#LIBS += -lflib -lrt @LIB_YAML_LIBS@


if BUILD_DYNAMIC_LIBS
flib_libraries = -lfutile-1
AM_LDFLAGS += @FC_RPATH@$(abs_top_builddir)/src
else
if USE_SIMGRID_SHARED_ALLOCS
AM_LDFLAGS += @FC_RPATH@$(abs_top_builddir)/src -Wl,-rpath=$(DESTDIR)$(libdir)
endif
flib_libraries = $(top_builddir)/src/libfutile-1.a
LIBS += @LIBS_EXT_DEPS@
endif

yaml_test_SOURCES = yaml_test.f90 errs.f90 yamlout.f90 dicts.f90 dynmem.f90 yaml_invoice.f90 utls.f90
yaml_test_LDADD = $(flib_libraries)

yaml_argparse_SOURCES = yaml_argparse.f90
yaml_argparse_LDADD = $(flib_libraries)

f_allocations_SOURCES = f_allocations.f90
f_allocations_LDADD = $(flib_libraries)


hooks_SOURCES = hooks.f90
hooks_LDADD = $(flib_libraries)

libpong_a_SOURCES = pong.f90
test_plugin_SOURCES = test_plugin.f90
test_plugin_LDADD = $(flib_libraries)

CLEANFILES = *.out.yaml *.report.yaml *.@MODULE_EXT@ memstatus.yaml report_remarks *.err

yaml_test.o: $(abs_top_builddir)/src/libfutile-1.a
hooks.o: $(abs_top_builddir)/src/libfutile-1.a

libpong.so: $(libpong_a_OBJECTS)
	$(FC) $(FCFLAGS) @FC_BUILD_SHARED@ $(LDFLAGS) @FC_RPATH@$(DESTDIR)$(dynlibdir) -o $@ $< $(flib_libraries)
	@chmod a+x $@

yaml_test.out.yaml: yaml_test
	./yaml_test > yaml_test.out.yaml

yaml_test.report.yaml: tols.yaml yaml_test.ref.yaml yaml_test.out.yaml
	python $(pythondir)/fldiff_yaml.py --label yaml_test -t $(srcdir)/tols.yaml 	\
		-r $(srcdir)/yaml_test.ref.yaml -d yaml_test.out.yaml --output $@

yaml_argparse.out.yaml: yaml_argparse
	./yaml_argparse > yaml_argparse.out.yaml
	./yaml_argparse -n 100 >> yaml_argparse.out.yaml
	./yaml_argparse -n 100 -d [89,90,90] >> yaml_argparse.out.yaml

yaml_argparse.report.yaml: tols.yaml yaml_argparse.ref.yaml yaml_argparse.out.yaml
	python $(pythondir)/fldiff_yaml.py --label yaml_argparse -t $(srcdir)/tols.yaml 	\
		-r $(srcdir)/yaml_argparse.ref.yaml -d yaml_argparse.out.yaml --output $@

hooks.out.yaml: hooks
	./hooks > hooks.out.yaml

hooks.report.yaml: tols.yaml hooks.ref.yaml hooks.out.yaml
	python $(pythondir)/fldiff_yaml.py --label hooks -t $(srcdir)/tols.yaml 	\
		-r $(srcdir)/hooks.ref.yaml -d hooks.out.yaml --output $@

check: yaml_test.report.yaml hooks.report.yaml regtest report

report:
	@if test $(MAKELEVEL) = 0 ; then python $(pythondir)/report.py ; fi

regtest: yaml_argparse
	python $(pythondir)/f_regtest.py --fldiff=$(pythondir)/fldiff_yaml.py \
	       --tols=$(srcdir)/tols.yaml --srcdir=$(srcdir) \
	       --exclusive=$(CHECK_ONLY)

define F_REGTEST_INSTRUCTIONS
- yaml_argparse:
   binary: yaml_argparse
   runs:
   - ./yaml_argparse > yaml_argparse.out.yaml
   - ./yaml_argparse -n 100 >> yaml_argparse.out.yaml
   - ./yaml_argparse -n 100 -d [89,90,90] >> yaml_argparse.out.yaml
   reference: $(srcdir)/yaml_argparse.ref.yaml
   output: yaml_argparse.out.yaml
- f_allocations:
   runs:
   - ./f_allocations > f_allocations.out.yaml
   - $(DESTDIR)$(bindir)/f_memcheck ./f_allocations >> f_allocations.out.yaml 2>f_allocations.err
endef
export F_REGTEST_INSTRUCTIONS
