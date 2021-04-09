if HAVE_PYTHON
pyscript = pyscript
pyscriptreport = pyscript.report.yaml
endif

check_PROGRAMS = dict cobjs yaml_output $(pyscript)

AM_LDFLAGS = -L$(top_builddir)/src
if BUILD_DYNAMIC_LIBS
flib_libraries = -lfutile-1
AM_LDFLAGS += @FC_RPATH@$(abs_top_builddir)/src
else
flib_libraries = $(top_builddir)/src/libfutile-1.a @LIBS_EXT_DEPS@
endif

cobjs_SOURCES = cobjs.c cobjs_main.f90
cobjs_LDADD = $(flib_libraries) @GLIB_LIBS@ $(FCLIBS)

dict_SOURCES = dict.c
dict_LDADD =  $(flib_libraries) @GLIB_LIBS@ $(FCLIBS)

yaml_output_SOURCES = yaml_output.c
yaml_output_LDADD =  $(flib_libraries) @GLIB_LIBS@ $(FCLIBS)

pyscript_SOURCES = objs.f90
pyscript_LDADD = $(flib_libraries)

AM_CPPFLAGS = -I$(top_builddir)/dicts -I$(top_srcdir)/dicts -I$(top_builddir)/c-bindings -I$(top_srcdir)/c-bindings @GLIB_CFLAGS@ @PYTHON_INCLUDES@
AM_FCFLAGS = -I$(top_builddir)/dicts -I$(top_builddir)/flib -I$(top_builddir)/src
AM_CFLAGS = -I$(top_builddir)/dicts -I$(top_srcdir)/dicts

CLEANFILES = *.@MODULE_EXT@

pyscript.out.yaml: pyscript
	PYTHONPATH=$PYTHONPATH:$(abs_top_builddir)/src ./pyscript > pyscript.out.yaml

pyscript.report.yaml: tols.yaml pyscript.ref.yaml pyscript.out.yaml
	$(PYTHON) $(top_builddir)/tests/fldiff_yaml.py --label pyscript -t $(srcdir)/tols.yaml	\
		-r $(srcdir)/pyscript.ref.yaml -d pyscript.out.yaml --output $@

check: $(pyscriptreport) report

report:
	@if test $(MAKELEVEL) = 0 ; then $(PYTHON) $(top_builddir)/tests/report.py ; fi
