# Define a macro to test how the Fortran compiler pass option to the preprocessor.
#
# Copyright (c) 2011-2011 BigDFT Group (Damien Caliste)
# All rights reserved.
#
# This file is part of the BigDFT software package. For license information,
# please see the COPYING file in the top-level directory of the BigDFT source
# distribution.

AC_DEFUN([AX_FC_DEFS],
[
  AC_MSG_CHECKING([for preprocessor option flag.])

  AC_LANG_PUSH(Fortran)
  AC_REQUIRE([AC_PROG_FC])

  cat > defstest.F90 <<EOF
#if defined TEST_FLAG
! choke me
#else
choke me
#endif

program test_defs

  write(*,*) "yes"
end program test_defs
EOF

  dnl Assume first -D.
  ax_fc_defs="unknown"
  ac_try='$FC -DTEST_FLAG $FCFLAGS -c defstest.F90 1>&AC_FD_CC'
  if AC_TRY_EVAL(ac_try); then
    ac_try=""
    ax_fc_defs="";
  fi
  dnl Assume second -WF,-D (xlf and fujitsu).
  if test "$ax_fc_defs" == "unknown" ; then
    ac_try='$FC -WF,-DTEST_FLAG $FCFLAGS -c defstest.F90 1>&AC_FD_CC'
    if AC_TRY_EVAL(ac_try); then
      ac_try=""
      ax_fc_defs="-WF,";
    fi
  fi
  rm -f defstest*

  AC_SUBST([FCDEFS], [$ax_fc_defs])
  AC_LANG_POP(Fortran)

  AC_MSG_RESULT([$ax_fc_defs"-D"])
])


AC_DEFUN([AX_SIMGRID_SHARED],
[dnl Use the shared allocators from simgrid
  AC_REQUIRE([AX_FC_RPATH])
  ax_have_simgrid_shared="no"
  AC_ARG_ENABLE(simgrid-shared, AS_HELP_STRING([--enable-simgrid-shared], [Use Simgrid shared allocators for explicitely marked buffers]),
                           ax_have_simgrid_shared=$enableval, ax_have_simgrid_shared="no")
  AC_MSG_CHECKING([for simgrid shared allocators])
  AC_MSG_RESULT([$ax_have_simgrid_shared])
  AM_CONDITIONAL(USE_SIMGRID_SHARED_ALLOCS, test "$ax_have_simgrid_shared" = "yes")
  if test "$ax_have_simgrid_shared" = x"yes" ; then
    AC_DEFINE([HAVE_SIMGRID_SHARED_ALLOCS], [1], [use SIMGRID allocators.])
  fi
])

