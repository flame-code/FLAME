noinst_LIBRARIES = libbindings.a
include_HEADERS  = dict.h futile.h
nodist_include_HEADERS  = futile_cst.h

libbindings_a_SOURCES = misc.c futile.h dict.c dict.h dict-fapi.h dictf.f90 python.c

AM_FCFLAGS = -I$(top_builddir)/flib
AM_CPPFLAGS = -I. -I$(top_builddir) @GLIB_CFLAGS@ @PYTHON_INCLUDES@ -D'PYTHON_EXECDIR="$(pyexecdir)"'
