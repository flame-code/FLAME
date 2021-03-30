
# To be installed:
dynlibdir        = $(libdir)
lib_LIBRARIES    = libfutile-1.a
noinst_LIBRARIES = $(pylib)
dynlib_DATA      = $(futile_dynlib)
noinst_DATA      = $(futile_linklib)
module_DATA      = $(modules)
pyexec_DATA      = $(pymodule)
dist_bin_SCRIPTS = $(memcheck_script)

if CC_AUTOBUILD_DYNAMIC_LIBS
memcheck_script = f_memcheck
else
memcheck_script =
endif 

if HAVE_PYTHON
if BUILD_DYNAMIC_LIBS
pymodule = futile.so
pylib = libpyimp.a
endif
endif

# Static or dynamic library
futile_dynamic_library = libfutile-1.so.@FUTILE_MINOR_VERSION@.0.@FUTILE_MICRO_VERSION@
futile_library_soname  = libfutile-1.so.@FUTILE_MINOR_VERSION@

# We install futile library only if PSolver is deactivated.
if BUILD_DYNAMIC_LIBS
futile_dynlib = $(futile_dynamic_library)
futile_linklib = $(futile_library_soname)
else
futile_dynlib =
futile_linklib =
endif

libfutile_1_a_SOURCES = f_lib_highlev.f90

libpyimp_a_SOURCES = pyobjs.c

f_lib_highlev.o: f_lib_highlev.f90 configure.inc

if CAPITALIZE
modules = F_LIB_PACKAGE.@MODULE_EXT@ FUTILE.@MODULE_EXT@
else
modules = f_lib_package.@MODULE_EXT@ futile.@MODULE_EXT@
endif

AM_FCFLAGS = -I. -I$(top_builddir)/dicts/ -I$(top_builddir)/flib @MPI_INCLUDE@ @FC_PIC_FLAG@
AM_CPPFLAGS = -I. -I$(top_builddir) -I$(top_builddir)/c-bindings -I$(top_srcdir)/c-bindings -I$(top_builddir)/dicts/ -I$(top_srcdir)/dicts @GLIB_CFLAGS@ @PYTHON_INCLUDES@ -D'PYTHON_EXECDIR="$(pyexecdir)"'

if BUILD_DYNAMIC_LIBS
flib_libraries = -lfutile-1
else
flib_libraries = libfutile-1.a @LIBS_EXT_DEPS@
endif

tmp-libfutile/extract.stamp: $(top_builddir)/dicts/libdicts.a $(top_builddir)/flib/libflib.a $(top_builddir)/c-bindings/libbindings.a $(top_builddir)/wrappers/libwrappers.a $(top_builddir)/wrappers/mpi/libfmpi.a
	test -e "tmp-libfutile" || $(INSTALL) -d -m 755 tmp-libfutile
	cd tmp-libfutile ; \
	rm -f *; \
	for lib in $^; do \
	  $(AR) x ../$$lib ; \
	done ; \
	touch extract.stamp

libfutile-1.a: $(libfutile_1_a_OBJECTS) tmp-libfutile/extract.stamp insmod
	$(AR) $(ARFLAGS) $@ $(libfutile_1_a_OBJECTS) tmp-libfutile/*.o
	$(RANLIB) $@

# Dynamic library building
$(futile_dynamic_library): $(libfutile_1_a_OBJECTS) tmp-libfutile/extract.stamp
	$(FC) $(FCFLAGS) @FC_BUILD_SHARED@ $(LDFLAGS) @FC_RPATH@$(DESTDIR)$(dynlibdir) -Wl,-E -Wl,-soname=$(futile_library_soname) -o $@ $(libfutile_1_a_OBJECTS) tmp-libfutile/*.o @LIBS_EXT_DEPS@ $(LIBS)
	@chmod a+x $@

$(futile_library_soname): $(futile_dynamic_library)
	ln -fs $^ $@
	ln -fs $@ libfutile-1.so

futile.so: $(libpyimp_a_OBJECTS) libfutile-1.a $(futile_library_soname)
	$(CC) $(CFLAGS) $(CPPFLAGS) @CC_BUILD_SHARED@ -L$(builddir) $(LDFLAGS) $(AM_LDFLAGS) @CC_RPATH@$(DESTDIR)$(dynlibdir) -Wl,-soname=$@ -o $@ $(libpyimp_a_OBJECTS) $(flib_libraries) $(LIBS) @PYTHON_LIBS@
	@chmod a+x $@

install-data-hook: 
	if test -n "$(futile_dynlib)" ; then \
	  cd $(DESTDIR)$(dynlibdir) ; \
	  chmod a+x $(futile_dynamic_library) ; \
	  ln -fs $(futile_dynamic_library) $(futile_library_soname) ; \
	  ln -fs $(futile_library_soname) libfutile-1.so ; \
	fi
	if test -n "$(pymodule)" ; then \
	  cd $(DESTDIR)$(pyexecdir) ; \
	  chmod a+x $(pymodule) ; \
	fi

CLEANFILES = libfutile-1.so* *.@MODULE_EXT@ futile.so
clean-local:
	rm -rf tmp-libfutile

$(modules): $(libfutile_1_a_OBJECTS)

# Copy modules in a common include directory.
insmod: $(modules)
	test -e "$(top_builddir)/includes" || $(INSTALL) -d -m 755 $(top_builddir)/includes
	for mod in "$(modules)" ; do \
	  $(INSTALL) -m 644  $$mod $(top_builddir)/includes ; \
	done
