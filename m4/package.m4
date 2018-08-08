# -*- Autoconf -*-
#
# Copyright (c) 2016 BigDFT Group (Bastian Schaefer)
# All rights reserved.
#
# This file is part of the BigDFT software package. For license information,
# please see the COPYING file in the top-level directory of the BigDFT source
# distribution.
#
AC_DEFUN([AX_EXTRA],
[dnl test for the presence of the package
  define([lcv],[translit([[$1]], [A-Z], [a-z])])
  AC_ARG_WITH(lcv, AS_HELP_STRING([--with-lcv],
              [Give the path of the lcv libraries (example = /home/<username>/lcv/). Do not use the -L before the path(es), just give the plain path.]),
              ac_lcv_dir=$withval, ac_lcv_dir=)
  ax_have_cv="no"
  if test -n "$ac_lcv_dir" ; then
    LIB_$1_LIBS="$ac_lcv_dir"
    ax_have_lcv="yes"
  fi
  AM_CONDITIONAL(HAVE_$1, test "$ax_have_lcv" = "yes")
])

AC_ARG_VAR([LAMMPS_ROOT], [lammps root path])
AC_ARG_VAR([SPGLIB_ROOT], [spglib root path])
AC_ARG_VAR([BDIR], [BigDFT root path])

