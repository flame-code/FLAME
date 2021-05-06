#aclocaldir = $(datadir)/aclocal
#dist_aclocal_DATA = config/m4/acc.m4 \
#	config/m4/defs.m4 \
#	config/m4/dynlib.m4 \
#	config/m4/fortran.m4 \
#	config/m4/linalg.m4 \
#	config/m4/mpi.m4 \
#	config/m4/pic.m4 \
#	config/m4/pkg.m4 \
#	config/m4/python.m4 \
#	config/m4/yaml.m4 \
#	config/m4/flib.m4
#
EXTRA_DIST= config/m4

ACLOCAL_AMFLAGS = -I config/m4

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = futile.pc libdicts.pc

SUBDIRS = dicts flib c-bindings wrappers src python

dist_bin_SCRIPTS = unienv.sh
nodist_bin_SCRIPTS = futilevars.sh futile_environment.sh

DISTCLEANFILES = config.inc

if HAVE_PYYAML
SUBDIRS += tests

check:
	@if test $(MAKELEVEL) = 0 ; then $(MAKE) report ; fi

report:
	python $(top_builddir)/tests/report.py
endif

SUBDIRS += doc

clean-local:
	rm -rf includes
