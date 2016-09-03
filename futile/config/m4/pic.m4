# Define a macro to test if the Fortran and C compiler support -fPIC
#
# Copyright (c) 2011-2013 BigDFT Group (Damien Caliste)
# All rights reserved.
#
# This file is part of the BigDFT software package. For license information,
# please see the COPYING file in the top-level directory of the BigDFT source
# distribution.

AC_DEFUN([AX_FLAG_PIC],
[
  f90_test_pic()
  {
    FCFLAGS_SVG=$FCFLAGS
    FCFLAGS="$FCFLAGS $[1]"
    AC_COMPILE_IFELSE([AC_LANG_SOURCE([
  subroutine api_func(a)
    integer, intent(out) :: a
  
    a = 1
  end subroutine api_func])],
      [ax_fc_pic=$[1]], [ax_fc_pic="no"])
    FCFLAGS=$FCFLAGS_SVG
  }
  f90_search_pic()
  {
    AC_LANG_PUSH(Fortran)
    for opt in $[@] ; do
      f90_test_pic $opt
      if test x"$ax_fc_pic" != x"no" ; then
        return
      fi
    done
    AC_LANG_POP(Fortran)
  }

  AC_MSG_CHECKING([for position-independant code option flag for $FC])
  dnl -qpic should be before -fPIC because -fPIC means something for xlf...
  f90_search_pic "-qpic" "-PIC" "-fPIC"
  AC_MSG_RESULT([$ax_fc_pic])

  cc_test_pic()
  {
    CFLAGS_SVG=$CFLAGS
    CFLAGS="$CFLAGS $[1]"
    AC_COMPILE_IFELSE([AC_LANG_SOURCE([
int api_func()
{
 return 1;
}])],
      [ax_cc_pic=$[1]], [ax_cc_pic="no"])
    CFLAGS=$CFLAGS_SVG
  }
  cc_search_pic()
  {
    AC_LANG_PUSH(C)
    for opt in $[@] ; do
      cc_test_pic $opt
      if test x"$ax_cc_pic" != x"no" ; then
        return
      fi
    done
    AC_LANG_POP(C)
  }

  AC_MSG_CHECKING([for position-independant code option flag for $CC])
  cc_search_pic "-qpic" "-PIC" "-fPIC"
  AC_MSG_RESULT([$ax_cc_pic])
])

  test_compiler_id()
  {
    compiler_basename=$($[1] --version | $SED 1q 2> /dev/null)
    if test -z "$compiler_basename" ; then
      compiler_basename=$($[1] -V | $SED 5q 2> /dev/null)
    fi
    echo $compiler_basename
  }

  test_compiler_options()
  {
     case $[1] in
      # GCC compiler.
      gcc*)
        ax_compiler_wl='-Wl,'
        ax_compiler_pic='-fPIC'
        ax_compiler_static='-static'
        ;;
      # old Intel for x86_64 which still supported -KPIC.
      ecc*)
        ax_compiler_wl='-Wl,'
        ax_compiler_pic='-KPIC'
        ax_compiler_static='-static'
        ;;
      # icc used to be incompatible with GCC.
      # ICC 10 doesn't accept -KPIC any more.
      icc* | ifort*)
        ax_compiler_wl='-Wl,'
        ax_compiler_pic='-fPIC'
        ax_compiler_static='-static'
        ;;
      # Lahey Fortran 8.1.
      lf95*)
        ax_compiler_wl='-Wl,'
        ax_compiler_pic='--shared'
        ax_compiler_static='--static'
    ;;
      nagfor*)
    # NAG Fortran compiler
        ax_compiler_wl='-Wl,-Wl,,'
        ax_compiler_pic='-PIC'
        ax_compiler_static='-Bstatic'
    ;;
      pgcc* | pgf77* | pgf90* | pgf95* | pgfortran*)
        # Portland Group compilers (*not* the Pentium gcc compiler,
    # which looks to be a dead project)
        ax_compiler_wl='-Wl,'
        ax_compiler_pic='-fpic'
        ax_compiler_static='-Bstatic'
        ;;
      ccc*)
        ax_compiler_wl='-Wl,'
        # All Alpha code is PIC.
        ax_compiler_static='-non_shared'
        ;;
      xl* | bgxl* | bgf* | mpixl*)
    # IBM XL C 8.0/Fortran 10.1, 11.1 on PPC and BlueGene
        ax_compiler_wl='-Wl,'
        ax_compiler_pic='-qpic'
        ax_compiler_static='-qstaticlink'
        ;;
     *Sun\ Ceres\ Fortran* | *Sun*Fortran*\ [1-7].* | *Sun*Fortran*\ 8.[0-3]*)
      # Sun Fortran 8.3 passes all unrecognized flags to the linker
        ax_compiler_pic='-KPIC'
        ax_compiler_static='-Bstatic'
        ax_compiler_wl=''
        ;;
     *Sun\ F* | *Sun*Fortran*)
        ax_compiler_pic='-KPIC'
        ax_compiler_static='-Bstatic'
        ax_compiler_wl='-Qoption ld '
        ;;
     *Sun\ C*)
      # Sun C 5.9
        ax_compiler_pic='-KPIC'
        ax_compiler_static='-Bstatic'
        ax_compiler_wl='-Wl,'
        ;;
     *Intel*\ [CF]*Compiler*)
        ax_compiler_wl='-Wl,'
        ax_compiler_pic='-fPIC'
        ax_compiler_static='-static'
        ;;
     *Portland\ Group*)
        ax_compiler_wl='-Wl,'
        ax_compiler_pic='-fpic'
        ax_compiler_static='-Bstatic'
        ;;
    esac
  }