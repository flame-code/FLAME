AC_PREREQ(2.59)

dnl Init basic data for package
dnl define following variables :
dnl  - PACKAGE_NAME, the name ;
dnl  - PACKAGE_TARNAME, the name used for the tarball ;
dnl  - PACKAGE_VERSION, the version ;
dnl  - PACKAGE_STRING, the name + the version ;
dnl  - PACKAGE_BUGREPORT, the mail.
AC_INIT([Futile],[1.8],[Damien.Caliste@cea.fr],[futile])

dnl All temporary building files are put in this directory.
dnl his must be called before AM_INIT_AUTOMAKE
AC_CONFIG_AUX_DIR(config)
dnl Use a config.h file to store macro definitions, required for CUDA.
AC_CONFIG_HEADERS(config.h)
AC_CONFIG_MACRO_DIR([config/m4])

dnl Init something, don't know exactly what...
dnl AM_INIT_AUTOMAKE([-Wno-portability, subdir-objects])
AM_INIT_AUTOMAKE(-Wno-portability)
dnl To keep the object files in the same structure as the source files (pb does not work with libABINIT)
dnl AM_INIT_AUTOMAKE([--warnings=no-portability, subdir-objects])

dnl Versioning
FUTILE_MAJOR_VERSION=$(echo $PACKAGE_VERSION | cut -d'.' -f1)
FUTILE_MINOR_VERSION=$(echo $PACKAGE_VERSION | cut -d'.' -f2 | cut -d'-' -f1)
FUTILE_MICRO_VERSION=$(echo $PACKAGE_VERSION | cut -d'.' -f3)
if test -z "$FUTILE_MICRO_VERSION" ; then
  FUTILE_MICRO_VERSION=0
fi
FUTILE_NANO_VERSION=$(echo $PACKAGE_VERSION | cut -d'.' -f4)
AC_SUBST(FUTILE_MAJOR_VERSION)
AC_SUBST(FUTILE_MINOR_VERSION)
AC_SUBST(FUTILE_MICRO_VERSION)
AC_SUBST(FUTILE_NANO_VERSION)

dnl Check for sed and put it in the SED variable
AC_PROG_SED()

dnl We defined CC for the bindings and CXX for FFTW.
AC_PROG_CC()
AM_PROG_CC_C_O()
AX_CC_RPATH()

dnl If FC environment variable is not already set, it
dnl looks for a modern Fortran compiler, prefering newest
dnl fortran available (e. g. Fortran95). Then FC variable
dnl is set.
AC_PROG_FC([], [Fortran 90])

dnl Fortran 77 (to be removed)
AC_PROG_F77([], [Fortran 77])

dnl Defining Fortran language for tests.
AC_LANG(Fortran)
AC_FC_SRCEXT(f90)

AC_FC_WRAPPERS()
AX_FC_FLUSH()
AX_FC_RECL_KIND()
AC_SUBST(RECL_INT_KIND, $ax_fc_recl_kind)
AX_FC_DEFS()
AX_FC_F2003()
AX_FC_OPENMP()
AX_FC_RPATH()

dnl Set if the module name are capitalized or not with the FC compiler.
AX_FC_MOD()
AM_CONDITIONAL(CAPITALIZE, test "$ax_fc_mod_capitalize" = "yes")
AC_SUBST(MODULE_EXT, $ax_fc_mod_ext)
dnl default installation directories
AX_FC_MODULEDIR()

dnl Set the FCFLAGS and FFLAGS variable
if test -z "$FFLAGS" ; then
  AC_SUBST(FFLAGS, "-O2")
fi
if test -n "$FCFLAGS_f90" ; then
  FCFLAGS="$FCFLAGS $FCFLAGS_f90"
fi
if test -z "$FCFLAGS" ; then
  AC_SUBST(FCFLAGS, "-O2")
fi

dnl Look for ranlib and canonicalize it with the $host variable
AC_CHECK_TOOL(RANLIB, ranlib, :)
if test "$ranlib" = ":" ; then
  AC_MSG_ERROR(["No 'ranlib' program found."])
fi
dnl Look for ar and canonicalize it with the $host variable
AC_CHECK_TOOL(AR, ar, :)
if test "$ar" = ":" ; then
  AC_MSG_ERROR(["No 'ar' program found."])
fi

dnl Enable or not MPI.
AX_MPI()

dnl C-level YAML support.
AX_YAML()

dnl Python YAML support.
AX_PYYAML()

dnl Produce dynamic excecutables.
AX_DYNAMIC_LIBRARIES()

dnl Test the availability of clock_gettime() and rt library.
dnl This will be added in libs, on top, therefore care should be taken for a linker which require ordering
AC_LANG_PUSH(C)
AC_CHECK_HEADERS([time.h])
AC_CHECK_LIB([rt], [clock_gettime])
AC_CHECK_LIB([dl], [dlopen])
AC_CHECK_FUNCS([clock_gettime])
AC_CHECK_FUNCS([strndup])
AC_CHECK_FUNCS([aligned_alloc])
AC_LANG_POP(C)

dnl Test Python availability
AX_PYTHON_DEV($ax_build_dynamic)

dnl Test Glib availability
AC_ARG_WITH(gobject, AS_HELP_STRING([--with-gobject], [Compile with GLib (disabled by default).]),
              ax_have_glib=$withval, ax_have_glib="no")
if test x"$ax_have_glib" = x"yes" ; then
  PKG_CHECK_MODULES(GLIB, glib-2.0 gobject-2.0 >= 2.22, [ax_have_glib=yes], [ax_have_glib=no])
fi
if test x"$ax_have_glib" = x"yes" ; then
  AC_DEFINE([HAVE_GLIB], [], [If set, we can call glib.h])
  AC_SUBST([GLIB_PKG], ["glib-2.0"])
  AC_SUBST([GLIB_TRUE], [""])
  AC_SUBST([GLIB_END_TRUE], [""])
  AC_SUBST([GLIB_FALSE], ["/*"])
  AC_SUBST([GLIB_END_FALSE], ["*/"])
else
  AC_SUBST([GLIB_PKG], [""])
  AC_SUBST([GLIB_TRUE], ["/*"])
  AC_SUBST([GLIB_END_TRUE], ["*/"])
  AC_SUBST([GLIB_FALSE], [""])
  AC_SUBST([GLIB_END_FALSE], [""])
fi

dnl Test accelerators
AX_ACC_CUDA()
AX_ACC_OCL()
AX_ACC_MIC()

dnl do we use simgrid shared allocators or fake version?
AX_SIMGRID_SHARED()

dnl Linalg libraries.
AX_LINALG()
if test $ax_have_linalg != "yes" ; then
  AC_MSG_ERROR([Linear algebra libraries are not available.])
fi
dnl Use the dgemmsy
AX_DGEMMSY()

dnl Agregate all the optional libs.
LIBS_EXT_DEPS="$LIBCUDA_LIBS $LINALG_LIBS $LIB_YAML_LIBS $MPI_LDFLAGS $MPI_LIBS $GLIB_LIBS $PYTHON_LIBS"
AC_SUBST(LIBS_EXT_DEPS, [$LIBS_EXT_DEPS])
LIB_FUTILE_LIBS="-lfutile-1"
if test $ax_build_dynamic != "yes" ; then
  LIB_FUTILE_LIBS="$LIB_FUTILE_LIBS $LIBS_EXT_DEPS $LIBS"
fi
AC_SUBST(LIB_FUTILE_LIBS)
LIB_DICTS_LIBS=-ldicts
AC_SUBST(LIB_DICTS_LIBS)

dnl Test for the developers' documentation.
AC_CHECK_PROG(ac_have_doxygen, [doxygen], [yes], [no])
AC_ARG_ENABLE(devel-doc,
              AS_HELP_STRING([--enable-devel-doc],
                             [compile documentation for developers using Doxygen.]),
              ac_devel_doc=$enableval, ac_devel_doc="no")
AM_CONDITIONAL(BUILD_DEVDOC, test x"$ac_have_doxygen" = x"yes" -a x"$ac_devel_doc" = x"yes")
dnl Test for the documentation.
AC_CHECK_PROG(ac_have_pdflatex, [pdflatex], [yes], [no])
AM_CONDITIONAL(BUILD_PDFDOC, test x"$ac_have_pdflatex" = x"yes" -a x"$ac_devel_doc" = x"yes")

dnl LIB_PSOLVER_DEPS="$LIBETSFIO_LIBS $LIBCUDA_LIBS $LINALG_LIBS $MPI_LDFLAGS $MPI_LIBS"
dnl LIB_PSOLVER_PKG=""
dnl if test $ax_build_dynamic != "yes" ; then
dnl   LIB_PSOLVER_PKG="$LIB_PSOLVER_DEPS"
dnl fi
dnl AC_SUBST(LIB_PSOLVER_PKG)
dnl LIB_PSOLVER_DEPS="$LIB_PSOLVER_DEPS $LIB_FUTILE_LIBS"
dnl AC_SUBST(LIB_PSOLVER_DEPS)


dnl default installation directories
docdir="${docdir}"
AC_SUBST(docdir)

AC_SUBST(CONFIG_ARGS, [$ac_configure_args])

dnl Give the name of file.in to transform to file
AC_CONFIG_FILES([
Makefile
futile.pc
libdicts.pc
futile_environment.sh
futilevars.sh
config/Doxyfile
dicts/Makefile
dicts/futile_cst.h
dicts/tests/Makefile
doc/Makefile
doc/Doxyfile
flib/Makefile
flib/f_utils.inc
wrappers/Makefile
wrappers/mpi/Makefile
wrappers/mpi/mpi_config.inc
c-bindings/Makefile
src/configure.inc
src/Makefile
python/Makefile
python/futile/Makefile
tests/Makefile
tests/fft/Makefile
tests/flib/Makefile
tests/wrappers/Makefile
tests/c-bindings/Makefile
])

AC_OUTPUT

echo "
Futile configure summary:
========================

Basics:
  Fortran90 compiler:        $FC
  Fortran90 flags:           $FCFLAGS
  Fortran77 compiler:        $F77
  Fortran77 flags:           $FFLAGS
  Linker flags:              $LDFLAGS
  Linked libraries:          $LIBS_EXT_DEPS $LIBS

Build:
  Dynamic library:           $ax_build_dynamic

Options:
  With OpenMP:               $ax_fc_openmp_msg
  With MPI:                  $ac_use_mpi
   | include dir.:           $MPI_INCLUDE
   | linker flags:           $MPI_LDFLAGS
   | linked libraries:       $MPI_LIBS
   | MPI_INIT_THREAD avail:  $ac_with_mpinithread
   | MPI2 support:           $ac_with_mpi2
  With dgemmsy support:      $ax_have_dgemmsy
  With libyaml:              $ax_have_yaml
   | C flags:                $LIB_YAML_CFLAGS
   | linked libraries:       $LIB_YAML_LIBS

Installation paths:
  Source code location:      ${srcdir}
  Prefix:                    $prefix
  Exec prefix:               $exec_prefix
  Binaries:                  $bindir
  Static libraries:          $libdir
  Fortran modules:           $moduledir
  Documentation:             $docdir
"

if test x"$ac_use_mpi" = x"warn" -o x"$withscalapack" = x"warn" ; then
  echo "Warnings:"
  if test x"$ac_use_mpi" = x"warn" ; then
    echo "  MPI:                       disabled by tests"
  fi
  if test x"$withscalapack" = x"warn" ; then
    echo "  Scalapack binding:         disabled by tests"
  fi

fi

dnl Do not touch this line. On the Cray machine in Basel, the preprocessor
dnl in Fortran sources does not allow C comments (/*  */), so we remove them.
$SED -n '1h;1!H;${;g;s%/\*[[^*]]*\*/%%g;p;}' config.h > config.inc
