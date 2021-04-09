# Define a macro to test flush() as intrinsic statement.
#
# Copyright (c) 2011-2011 BigDFT Group (Damien Caliste)
# All rights reserved.
#
# This file is part of the BigDFT software package. For license information,
# please see the COPYING file in the top-level directory of the BigDFT source
# distribution.

AC_DEFUN([AX_FC_FLUSH],
[
  AC_MSG_CHECKING([for flush(6) in Fortran.])

  AC_LANG_PUSH(Fortran)
  AC_REQUIRE([AC_PROG_FC])

  cat > flushtest.f90 <<EOF
program test_flush

  implicit none

  write(*,"(A)") "yes"
  flush(6)

end program test_flush
EOF
  cat > flushtest_sub.f90 <<EOF
program test_flush

  implicit none

  write(*,"(A)") "yes"
  call flush(6)

end program test_flush
EOF
  dnl Assume first that it should compile and run.
  ax_fc_flush="no"
  ac_try='$FC $FCFLAGS $LDFLAGS -o flushtest.x flushtest.f90 $LIBS 1>&AC_FD_CC'
  if AC_TRY_EVAL(ac_try); then
    ac_try=""
    ax_fc_flush=`./flushtest.x 2> /dev/null`;
    if test "$?" != 0 ; then
      ax_fc_flush="no"
    fi
   fi
  dnl Assume second that it should compile and run with Intel option.
  FCFLAGS_SVG="$FCFLAGS"
  if test x"$ax_fc_flush" == x"no" ; then
    FCFLAGS="$FCFLAGS -assume noold_unit_star"
    ac_try='$FC $FCFLAGS $LDFLAGS -o flushtest.x flushtest.f90 $LIBS 1>&AC_FD_CC'
    if AC_TRY_EVAL(ac_try); then
      ac_try=""
      ax_fc_flush=`./flushtest.x 2> /dev/null`;
      if test "$?" != 0 ; then
        ax_fc_flush="no"
      fi
    fi
  fi
  if test x"$ax_fc_flush" != x"yes" ; then
    FCFLAGS="$FCFLAGS_SVG"
  fi
  dnl Assume third that it should compile to have it.
  if test x"$ax_fc_flush" == x"no" ; then
    ac_try='$FC $FCFLAGS $LDFLAGS -o flushtest.x flushtest.f90 $LIBS 1>&AC_FD_CC'
    if AC_TRY_EVAL(ac_try); then
      ax_fc_flush="yes"
    fi
  fi
  if test x"$ax_fc_flush" == x"yes" ; then
  AC_DEFINE([HAVE_FC_FLUSH], [1], [Flush(6) can be used safely in fortran])
  fi
  AM_CONDITIONAL([HAVE_FC_FLUSH], [test x"$ax_fc_flush" == x"yes"])
  dnl then assume that flush is a subroutine
  ax_fc_flush_sub="no"
  if test x"$ax_fc_flush" == x"no" ; then
    ac_try='$FC $FCFLAGS $LDFLAGS -o flushtest.x flushtest_sub.f90 $LIBS 1>&AC_FD_CC'
    if AC_TRY_EVAL(ac_try); then
      ac_try=""
      ax_fc_flush=`./flushtest.x 2> /dev/null`;
      if test "$?" == 0 ; then
        ax_fc_flush_sub="yes"
      fi
    fi
  fi
  rm -f flushtest*
  if test x"$ax_fc_flush_sub" == x"yes" ; then
    AC_DEFINE([HAVE_FC_FLUSH_SUB], [1], [call flush(6) can be used safely in fortran])
  fi
  AM_CONDITIONAL([HAVE_FC_FLUSH_SUB], [test x"$ax_fc_flush_sub" == x"yes"])
  AC_LANG_POP(Fortran)

  AC_MSG_RESULT([$ax_fc_flush])
])

# Define a macro to test if the inquire intrinsic statement accepts long integers as recl parameters
AC_DEFUN([AX_FC_RECL_KIND],
[
  AC_MSG_CHECKING([for selected_int_kind accepted for recl parameter in inquire Fortran statement])

  AC_LANG_PUSH(Fortran)
  AC_REQUIRE([AC_PROG_FC])

  cat > recltest.f90 <<EOF
program inq
  implicit none
  integer, parameter :: unt=16 !test a closed unit
  integer, parameter :: k8=selected_int_kind(16)
  integer :: ierr
  integer(k8) :: recl_file

  recl_file=int(-1234567891,k8)

  inquire(unit=unt,recl=recl_file,iostat=ierr)

end program inq
EOF
  dnl Assume that it should compile to have it, then the f_utils_recl function should take care of it
  ax_fc_recl_kind="8"
  ac_try='$FC $FCFLAGS $LDFLAGS -o recltest.x recltest.f90 1>&AC_FD_CC'
  if AC_TRY_EVAL(ac_try); then
    ax_fc_recl_kind="16"
  fi
  rm -f recltest*

  AC_LANG_POP(Fortran)

  AC_MSG_RESULT([$ax_fc_recl_kind])
])
