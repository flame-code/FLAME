# To be installed:
lib_LIBRARIES = libdicts.a
module_DATA   = $(modules)

futileincdir = $(includedir)/futile
futileinc_HEADERS = dict.h err.h
nodist_futileinc_HEADERS = futile_cst.h

libdicts_a_SOURCES = \
	yaml_strings.f90 dictionaries_base.f90 dictionaries.f90 \
	callbacks.f90 f_precisions.f90 \
	dict.c dict-fapi.h dictf.f90 \
	err.c err-fapi.h errf.f90 \
	addresses.h addresses.c

STRINGS_INCLUDES = yaml_toa-inc.f90 yaml_toa-arr-inc.f90 yaml_plus-inc.f90
ERR_INCLUDES = get_err-inc.f90
DICT_INCLUDES = error_handling.f90 $(ERR_INCLUDES) dict_add-inc.f90 dict_getvec-inc.f90 \
	set_arr-inc.f90 dict_cont-inc.f90 dict_cont_arr-inc.f90
NOOMP_INCLUDES = halt_omp-inc.f90

SUBDIRS = tests

EXTRA_DIST = $(STRINGS_INCLUDES) $(DICT_INCLUDES) $(NOOMP_INCLUDES)
DISTCLEANFILES = futile_cst.h

# Module files.
if CAPITALIZE
modules = YAML_STRINGS.@MODULE_EXT@ \
	DICTIONARIES.@MODULE_EXT@ \
	DICTIONARIES_BASE.@MODULE_EXT@ \
	F_PRECISIONS.@MODULE_EXT@ \
	EXCEPTION_CALLBACKS.@MODULE_EXT@
else
modules = yaml_strings.@MODULE_EXT@ \
	dictionaries.@MODULE_EXT@ \
	dictionaries_base.@MODULE_EXT@ \
	f_precisions.@MODULE_EXT@ \
	exception_callbacks.@MODULE_EXT@
endif

CLEANFILES =  *.@MODULE_EXT@
AM_CFLAGS = @CC_PIC_FLAG@
AM_FCFLAGS = -I. -I$(srcdir) @FC_PIC_FLAG@
AM_CPPFLAGS = -I. -I$(top_builddir) @GLIB_CFLAGS@

# Dependencies
$(modules): $(libdicts_a_OBJECTS)
yaml_strings.o: $(STRINGS_INCLUDES) f_precisions.o
dictionaries_base.o: yaml_strings.o
dictionaries.o: yaml_strings.o dictionaries_base.o $(DICT_INCLUDES) callbacks.o
callbacks.o: $(NOOMP_INCLUDES) f_precisions.o
dictf.o: dictionaries.o
errf.o: dictionaries.o
