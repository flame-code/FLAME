# -*- Autoconf -*-
#
# Copyright (c) 2016 BigDFT Group (Bastian Schaefer)
# All rights reserved.
#
# This file is part of the BigDFT software package. For license information,
# please see the COPYING file in the top-level directory of the BigDFT source
# distribution.
#

AC_DEFUN([AX_BPS],
[dnl test for the presence of bps
  AC_ARG_WITH(bps, AS_HELP_STRING([--with-bps],
              [Give the path of the bps libraries (example = /home/<username>/bps/). Do not use the -L before the path(es), just give the plain path.]),
              ac_bps_dir=$withval, ac_bps_dir=)
  ax_have_bps="no"
  if test -n "$ac_bps_dir" ; then
    LIB_BPS_LIBS="$ac_bps_dir"
    ax_have_bps="yes"
  fi
  AM_CONDITIONAL(HAVE_BPS, test $ax_have_bps = "yes")
])

AC_ARG_VAR([BDIR], [BigDFT root path])

