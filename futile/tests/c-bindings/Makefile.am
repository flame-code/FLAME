if HAVE_PYTHON
pyscript = pyscript
pyscriptreport = pyscript.report.yaml
endif

check_PROGRAMS = dict $(pyscript)

AM_LDFLAGS = -L$(top_builddir)/src
if BUILD_DYNAMIC_LIBS
flib_libraries = -lfutile-1
AM_LDFLAGS += @FC_RPATH@$(abs_top_builddir)/src
else
flib_libraries = $(top_builddir)/src/libfutile-1.a @LIBS_EXT_DEPS@
endif

dict_SOURCES = dict.c
dict_LDADD = $(flib_libraries) @GLIB_LIBS@ $(FCLIBS)

pyscript_SOURCES = objs.f90
pyscript_LDADD = $(flib_libraries)

AM_CPPFLAGS = -I$(top_builddir)/c-bindings -I$(top_srcdir)/c-bindings @GLIB_CFLAGS@ @PYTHON_INCLUDES@
AM_FCFLAGS = -I$(top_builddir)/flib -I$(top_builddir)/src

CLEANFILES = *.@MODULE_EXT@

pyscript.out.yaml: pyscript
	./pyscript > pyscript.out.yaml

pyscript.report.yaml: tols.yaml pyscript.ref.yaml pyscript.out.yaml
	python $(pythondir)/fldiff_yaml.py --label pyscript -t $(srcdir)/tols.yaml	\
		-r $(srcdir)/pyscript.ref.yaml -d pyscript.out.yaml --output $@

check: $(pyscriptreport) report

report:
	@if test $(MAKELEVEL) = 0 ; then python $(pythondir)/report.py ; fi
