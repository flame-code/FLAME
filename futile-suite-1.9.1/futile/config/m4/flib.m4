# -*- Autoconf -*-
#
# Copyright (c) 2015 BigDFT Group (Damien Caliste)
# All rights reserved.
#
# This file is part of the BigDFT software package. For license information,
# please see the COPYING file in the top-level directory of the BigDFT source
# distribution.
#
AC_DEFUN([AX_FLIB],
[dnl test for flib
AX_PACKAGE([FUTILE],[1.8],[-lfutile-1],[-lyaml -lrt],[],
 	     [program main
     use yaml_parse
     use yaml_output
     use f_utils
     use dynamic_memory
     use dictionaries
   
     call yaml_map("toto","titi")
   end program],
 	[call f_lib_initialize()
])
if test $ax_have_FUTILE != "yes" ; then
AX_PACKAGE([FUTILE],[1.8],[-lfutile-1],[-lyaml],[],
 	     [program main
     use yaml_parse
     use yaml_output
     use f_utils
     use dynamic_memory
     use dictionaries
   
     call yaml_map("toto","titi")
   end program],
 	[call f_lib_initialize()
])
if test $ax_have_FUTILE != "yes" ; then
  AC_MSG_ERROR([Futile library not found, cannot proceed.])
fi
fi

dnl Some installation directories.
AC_MSG_CHECKING([for installation python dir for FUTILE])
FUTILE_PYTHONDIR="/usr/local/lib/python2.7/dist-packages"
if test -n "$PKG_CONFIG"; then
    PKG_CHECK_EXISTS([futile],
                     [FUTILE_PYTHONDIR=`$PKG_CONFIG --variable=pythondir futile 2>/dev/null`],
                     [FUTILE_PYTHONDIR="none"])
fi
AC_MSG_RESULT([$FUTILE_PYTHONDIR])
AC_SUBST(FUTILE_PYTHONDIR)

AX_PACKAGE([LIBDICTS],[1.8],[-ldicts],[],[],
 	     [program main
     use dictionaries

     type(dictionary), pointer :: dict

     call dict_init(dict)
   end program],
 	[use f_precisions, only: f_address
	integer(f_address) :: iloc
	integer(f_address), external :: f_loc
	iloc=f_loc(iloc)
])
])
dnl AC_DEFUN([AX_FLIB],
dnl [dnl Test for FLib
dnl   AC_ARG_WITH(flib-libs, AS_HELP_STRING([--with-flib-libs], [Give the linker flags for an external FLib modules (default = None).]), ac_flib_libdir=$withval, ac_flib_libdir=)
dnl   AC_ARG_WITH(flib-incs, AS_HELP_STRING([--with-flib-incs], [Give the compiler include flags for an external FLib library (default = None).]), ac_flib_incdir=$withval, ac_flib_incdir=)
dnl   
dnl   dnl try first with pkg-config
dnl   AC_MSG_CHECKING([for FUTILE via PKG_CONFIG]) 
dnl   PKG_CHECK_MODULES([FUTILE],
dnl                     [futile >= 1.8],
dnl                     [ax_have_futile=yes],
dnl                     [ax_have_futile=no])
dnl   if test "$ax_have_futile" = "yes" ; then
dnl     if test -z "${FUTILE_CFLAGS// }" -a -n "$C_INCLUDE_PATH" ; then
dnl       for path in ${C_INCLUDE_PATH//:/ }; do
dnl         ax_futile_incdir="$ax_futile_incdir -I$path"
dnl       done
dnl       LIB_FUTILE_CFLAGS=$ax_futile_incdir
dnl     else
dnl       LIB_FUTILE_CFLAGS=$FUTILE_CFLAGS
dnl     fi
dnl     LIB_FUTILE_LIBS=$FUTILE_LIBS
dnl   fi
dnl   AC_MSG_RESULT("libs=" $LIB_FUTILE_LIBS " flags=" $LIB_FUTILE_CFLAGS "have=" $ax_have_futile)
dnl 
dnl   dnl try by hand search if failed
dnl   if test "$ax_have_futile" != "yes" ; then
dnl     dnl Test the modules for compilation
dnl     AC_LANG_PUSH(Fortran)
dnl     AC_REQUIRE([AC_PROG_FC])
dnl     
dnl     dnl Test the modules for compilation
dnl     AC_MSG_CHECKING([for FLib modules])
dnl     FCFLAGS_SVG=$FCFLAGS
dnl     if test -n "$ac_flib_incdir" ; then
dnl       FCFLAGS="$FCFLAGS $ac_flib_incdir"
dnl     elif test -n "$C_INCLUDE_PATH" ; then
dnl       for path in ${C_INCLUDE_PATH//:/ }; do
dnl         ac_flib_incdir="$ac_flib_incdir -I$path"
dnl       done
dnl       FCFLAGS="$FCFLAGS $ac_flib_incdir"
dnl     fi
dnl     AC_COMPILE_IFELSE([[program main
dnl     use yaml_parse
dnl     use yaml_output
dnl     use f_utils
dnl     use dynamic_memory
dnl     use dictionaries
dnl   
dnl     call yaml_map("toto", "titi")
dnl   end program]], withflibmod=yes, withflibmod=no)
dnl     AC_MSG_RESULT($withflibmod)
dnl     FCFLAGS=$FCFLAGS_SVG
dnl   
dnl     dnl Test the library of flib.
dnl     AC_MSG_CHECKING([for flib library])
dnl     LIBS_SVG=$LIBS
dnl     if test -z "$ac_flib_libdir" ; then
dnl       ac_flib_libdir="-lfutile-1"
dnl     fi
dnl     LIBS="$ac_flib_libdir $LIBS_SVG"
dnl     AC_LINK_IFELSE(
dnl       AC_LANG_PROGRAM([], [[
dnl   call f_lib_initialize()
dnl   ]]),
dnl       [ax_have_flib=yes],
dnl       [ax_have_flib=no])
dnl     if test $ax_have_flib != "yes" ; then
dnl       dnl Static case, need to link with additional libs.
dnl       ac_flib_libdir="$ac_flib_libdir -lyaml -lrt"
dnl       LIBS="$ac_flib_libdir $LIBS_SVG"
dnl       AC_LINK_IFELSE(
dnl         AC_LANG_PROGRAM([], [[
dnl   call f_lib_initialize()
dnl   ]]),
dnl         [ax_have_flib=yes],
dnl         [ax_have_flib=no])
dnl     fi
dnl     AC_MSG_RESULT($ax_have_flib)
dnl     LIBS=$LIBS_SVG
dnl   
dnl     if test "$ax_have_flib" = "yes" -a "$withflibmod" = "yes" ; then
dnl       LIB_FUTILE_CFLAGS=$ac_flib_incdir
dnl       LIB_FUTILE_LIBS=$ac_flib_libdir
dnl       ax_have_flib="yes"
dnl     else
dnl       ax_have_flib="no"
dnl     fi
dnl   fi
dnl   
dnl   dnl LIB_XC_CFLAGS="-I/usr/include"
dnl   dnl   PKG_CHECK_MODULES(LIB_XC, flib >= 2.0, ax_have_flib="yes", ax_have_flib="no")
dnl   
dnl   AC_SUBST(LIB_FUTILE_CFLAGS)
dnl   AC_SUBST(LIB_FUTILE_LIBS)
dnl 
dnl   dnl Try to find libflib-1.a for possible later inclusion.
dnl   
dnl   AC_LANG_POP(Fortran)
dnl ])
