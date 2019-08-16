# -*- Autoconf -*-
#
# Copyright (c) 2016 BigDFT Group (Bastian Schaefer)
# All rights reserved.
#
# This file is part of the BigDFT software package. For license information,
# please see the COPYING file in the top-level directory of the BigDFT source
# distribution.
#

AC_DEFUN([AX_MKL],
[dnl test for the presence of mkl
  AC_ARG_WITH(mkl, AS_HELP_STRING([--with-mkl],
              [Give the path of the mkl libraries (example = /home/<username>/mkl/). Do not use the -L before the path(es), just give the plain path.]),
              ac_mkl_dir=$withval, ac_mkl_dir=)
  ax_have_mkl="no"
  if test -n "$ac_mkl_dir" ; then
    LIB_MKL_LIBS="$ac_mkl_dir"
    ax_have_mkl="yes"
  fi
  AM_CONDITIONAL(HAVE_MKL, test $ax_have_mkl = "yes")
])


