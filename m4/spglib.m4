# -*- Autoconf -*-
#
# Copyright (c) 2016 BigDFT Group (Bastian Schaefer)
# All rights reserved.
#
# This file is part of the BigDFT software package. For license information,
# please see the COPYING file in the top-level directory of the BigDFT source
# distribution.
#

AC_DEFUN([AX_SPGLIB],
[dnl test for the presence of spglib
  AC_ARG_WITH(spglib, AS_HELP_STRING([--with-spglib],
              [Give the path of the spglib libraries (example = /home/<username>/spglib/). Do not use the -L before the path(es), just give the plain path.]),
              ac_spglib_dir=$withval, ac_spglib_dir=)
  ax_have_spglib="no"
  if test -n "$ac_spglib_dir" ; then
    LIB_SPGLIB_LIBS="$ac_spglib_dir"
    ax_have_spglib="yes"
  fi
  AM_CONDITIONAL(HAVE_SPGLIB, test $ax_have_spglib = "yes")
])

AC_ARG_VAR([SPGLIB_ROOT], [spglib root path])

