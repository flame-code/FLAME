noinst_LIBRARIES = libbindings.a

futileincdir = $(includedir)/futile
futileinc_HEADERS = tree.h misc.h enum.h f_string.h f_progress_bar.h
include_HEADERS = futile.h futileC_objects.h
#nodist_include_HEADERS  = futile_cst.h

futileinc_HEADERS += yaml_output.h

libbindings_a_SOURCES = tree.c tree.h misc.c misc.h enum.h futile.h python.c futileC_objects.c treef.f90

libbindings_a_SOURCES += yaml_output_c.c yaml_output.h yaml_output_bind.f90

AM_FCFLAGS = -I$(top_builddir)/dicts -I$(top_builddir)/flib
AM_CPPFLAGS = -I. -I$(top_builddir) -I$(top_builddir)/dicts -I$(top_srcdir)/dicts @GLIB_CFLAGS@ @PYTHON_INCLUDES@ -D'PYTHON_EXECDIR="$(pyexecdir)"'
AM_CFLAGS = -I$(top_builddir)/dicts
