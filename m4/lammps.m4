# -*- Autoconf -*-
#
# Copyright (c) 2016 BigDFT Group (Bastian Schaefer)
# All rights reserved.
#
# This file is part of the BigDFT software package. For license information,
# please see the COPYING file in the top-level directory of the BigDFT source
# distribution.
#

AC_DEFUN([AX_LAMMPS],
[dnl test for the presence of lammps
  AC_ARG_WITH(lammps, AS_HELP_STRING([--with-lammps],
              [Give the path of the lammps libraries (example = /home/<username>/lammps/). Do not use the -L before the path(es), just give the plain path.]),
              ac_lammps_dir=$withval, ac_lammps_dir=)
  ax_have_lammps="no"
  if test -n "$ac_lammps_dir" ; then
    LIB_LAMMPS_LIBS="$ac_lammps_dir"
    ax_have_lammps="yes"
  fi
  AM_CONDITIONAL(HAVE_LAMMPS, test $ax_have_lammps = "yes")
])

AC_ARG_VAR([LAMMPS_ROOT], [lammps root path])

