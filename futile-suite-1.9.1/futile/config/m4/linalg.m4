# -*- Autoconf -*-
#
# Copyright (c) 2014 BigDFT Group (Damien Caliste)
# All rights reserved.
#
# This file is part of the BigDFT software package. For license information,
# please see the COPYING file in the top-level directory of the BigDFT source
# distribution.
#

AC_DEFUN([AX_LINALG],
[dnl Substitute Lapack and Blas by another specific library
  AC_REQUIRE([AX_MPI])
  AC_LANG_PUSH(Fortran)
  
  ax_have_linalg="no"

  withlinalg=no
  AC_ARG_WITH(ext-linalg, AS_HELP_STRING([--with-ext-linalg],
              [Give the name of the libraries replacing Blas and Lapack (default = none specified). Use the -l before the name(s).]),
              ac_linalg=$withval, ac_linalg=)
  AC_ARG_WITH(ext-linalg-path, AS_HELP_STRING([--with-ext-linalg-path],
              [Give the path of the other linear algebra libraries (default = -L/usr/lib). Use the -L before the path(es).]),
              ac_linalg_dir=$withval, ac_linalg_dir=)
  if test -n "$ac_linalg_dir" ; then
    LDFLAGS="$LDFLAGS $ac_linalg_dir"
  fi
  if test -n "$ac_linalg" ; then
  #  AC_LANG_CONFTEST([[program main
  #                     call dsysv
  #                     end program main]])
    AC_MSG_CHECKING([for dsysv in $ac_linalg])
    LIBS_OLD=$LIBS
    LIBS="$ac_linalg $LIBS"
    AC_LINK_IFELSE([[program main
    call dsysv
  end]], withlinalg=yes, withlinalg=no)
  #  AC_CHECK_LIB($ac_linalg, dsysv, withlinalg=yes, withlinalg=no)
    if test "$withlinalg" = "yes"; then
      LINALG_LIBS=$ac_linalg
      ax_have_linalg="yes"
    fi
    AC_MSG_RESULT([$withlinalg])
    dnl test for scalapack presence in the ext-linalg
    ax_have_scalapack=no
    if test "$withlinalg" = "yes"; then
     if test x"$withmpi" = x"yes" ; then
      AC_MSG_CHECKING([for pdsygvx and blacs_get in $ac_linalg])
      AC_LINK_IFELSE([[program main
      call pdsygvx
      call blacs_get
  end]], ax_have_scalapack=yes, ax_have_scalapack=warn)
      AC_MSG_RESULT([$ax_have_scalapack])
     fi
    fi
    LIBS=$LIBS_OLD
  fi
  
  if test "$withlinalg" = "no" ; then
  dnl Get the Blas library
  blas=yes
  AC_ARG_WITH(blas, AS_HELP_STRING([--with-blas], [Link with Blas library (default = yes).]), blas=$withval, blas=yes)
  AC_ARG_WITH(blas-path, AS_HELP_STRING([--with-blas-path], [Give the path of the Blas library (default = /usr/lib).]),
              ac_blas_dir=$withval, ac_blas_dir=)
  if test -n "$ac_blas_dir" ; then
    LDFLAGS="$LDFLAGS -L$ac_blas_dir"
  fi
  withblas=no
  if test "$blas" = "yes" ; then
    AC_CHECK_LIB(blas, dcopy, withblas=yes, withblas=no)
    if test "$withblas" = "yes"; then
      LINALG_LIBS="-lblas"
    fi
  fi
  
  dnl Get the lapack library
  lapack=yes
  AC_ARG_WITH(lapack, AS_HELP_STRING([--with-lapack], [Link with Lapack library (default = yes).]), lapack=$withval, lapack=yes)
  AC_ARG_WITH(lapack-path, AS_HELP_STRING([--with-lapack-path], [Give the path of the Lapack library (default = /usr/lib).]),
              ac_lapack_dir=$withval, ac_lapack_dir=)
  if test -n "$ac_lapack_dir" ; then
    LDFLAGS="$LDFLAGS -L$ac_lapack_dir"
  fi
  withlapack=no
  if test "$lapack" = "yes" ; then
    AC_CHECK_LIB(lapack, dsysv, withlapack=yes, withlapack=no,$LINALG_LIBS)
    if test "$withlapack" = "yes"; then
      LINALG_LIBS="-llapack $LINALG_LIBS"
    fi
  fi

  if test $withblas = "yes" -a $withlapack = "yes" ; then
    ax_have_linalg="yes"
  fi
  
  dnl Test for Scalapack
  if test x"$withmpi" = x"yes" ; then
  dnl Test for blacs
  ac_blacs=-lblacs
  AC_ARG_WITH(blacs, AS_HELP_STRING([--with-blacs], [Link with blacs library (default = no).]), ac_blacs=$withval, ac_blacs=-lblacs)
  AC_ARG_WITH(blacs-path, AS_HELP_STRING([--with-blacs-path], [Give the path of the Blacs library (default = /usr/lib).]),
              ac_blacs_dir=$withval, ac_blacs_dir=)
  if test -n "$ac_blacs_dir" ; then
    LDFLAGS="$LDFLAGS -L$ac_blacs_dir"
  fi
  withblacs=no
  if test "$ac_blacs" != "no" ; then
    case $ac_blacs in
      -l*)
      AC_MSG_CHECKING([for blacs libraries])
      AC_LANG_CONFTEST([AC_LANG_CALL([], [blacs_get])])
      LIBS_SVG=$LIBS
      for ac_lib in $ac_blacs '-lblacs -lblacsF77init' '-lblacs-openmpi' '-lblacs-openmpi -lblacsF77init-openmpi'; do
        LIBS="$ac_lib $LIBS_SVG"
        AC_LINK_IFELSE([], [ac_blacs=$ac_lib], [ac_blacs="no"])
        if test x"$ac_blacs" != x"no" ; then
          break
        fi
      done
      rm conftest*
      LIBS=$LIBS_SVG
      AC_MSG_RESULT([$ac_blacs])
      if test x"$ac_blacs" != x"no" ; then
        withblacs=yes
      else
        withblacs=warn
      fi
      ;;
      *)
      AC_CHECK_LIB($ac_blacs, blacs_get, withblacs=yes, withblacs=warn,$LINALG_LIBS)
      ac_blacs=-l$ac_blacs
      ;;
    esac
    if test "$withblacs" = "yes"; then
      LINALG_LIBS="$ac_blacs $LINALG_LIBS"
    fi
  fi
  dnl get the scalapack library
  ac_scalapack=-lscalapack
  AC_ARG_WITH(scalapack, AS_HELP_STRING([--with-scalapack], [Link with scalapack library (default = no).]), ac_scalapack=$withval, ac_scalapack=-lscalapack)
  AC_ARG_WITH(scalapack-path, AS_HELP_STRING([--with-scalapack-path], [Give the path of the Scalapack library (default = /usr/lib).]),
              ac_scalapack_dir=$withval, ac_scalapack_dir=)
  if test -n "$ac_scalapack_dir" ; then
    LDFLAGS="$LDFLAGS -L$ac_scalapack_dir"
  fi
  ax_have_scalapack=no
  if test "$ac_scalapack" != "no" -a "$ac_blacs" != "no" ; then
    case $ac_scalapack in
      -l*)
      AC_MSG_CHECKING([for scalapack libraries])
      AC_LANG_CONFTEST([AC_LANG_CALL([], [pdsygvx])])
      LIBS_SVG=$LIBS
      for ac_lib in $ac_scalapack '-lscalapack-openmpi'; do
        LIBS="$ac_lib $LIBS_SVG"
        AC_LINK_IFELSE([], [ac_scalapack=$ac_lib], [ac_scalapack="no"])
        if test x"$ac_scalapack" != x"no" ; then
          break
        fi
      done
      rm conftest*
      LIBS=$LIBS_SVG
      AC_MSG_RESULT([$ac_scalapack])
      if test x"$ac_scalapack" != x"no" ; then
        ax_have_scalapack=yes
      else
        ax_have_scalapack=warn
      fi
      ;;
      *)
      AC_CHECK_LIB($ac_scalapack, pdsygvx, ax_have_scalapack=yes, ax_have_scalapack=warn,$LINALG_LIBS)
      ac_scalapack=-l$ac_scalapack
      ;;
    esac
    if test "$ax_have_scalapack" = "yes"; then
      LINALG_LIBS="$ac_scalapack $LINALG_LIBS"
    fi
  fi
     
  fi #end $withmpi = yes
  
  fi #end $withlinalg = no
  AC_SUBST(LINALG_LIBS, $LINALG_LIBS)

  AC_LANG_POP([Fortran])
  dnl ---------------------------------------
  AM_CONDITIONAL(USE_BLACS, test x"$ax_have_scalapack" = x"yes")
])

AC_DEFUN([AX_DGEMMSY],
[dnl define the dgemmsy compilation or not
  DGEMMSY_CPPFLAGS=
  AC_ARG_ENABLE(dgemmsy, AS_HELP_STRING([--enable-dgemmsy], [Use dgemmsy (disabled by default).]),
                           ax_have_dgemmsy=$enableval, ax_have_dgemmsy="no")
  if test x"$ax_have_dgemmsy" = "xyes" ; then
    AC_LANG_PUSH(C)
  
    DGEMMSY_CPPFLAGS="-msse3"
    CPPFLAGS_SVG=$CPPFLAGS
    CPPFLAGS="$CPPFLAGS $DGEMMSY_CPPFLAGS"
    AC_CHECK_HEADER([pmmintrin.h],
                    [ax_have_dgemmsy="yes"],
                    [ax_have_dgemmsy="no"])
    CPPFLAGS=$CPPFLAGS_SVG
    AC_LANG_POP(C)
    if test "$ax_have_dgemmsy" = "no"; then
      DGEMMSY_CPPFLAGS=
      AC_MSG_WARN(["No SSE header found, dgemmsy disabled."])
    fi
  fi
  AC_SUBST(DGEMMSY_CPPFLAGS)
  AM_CONDITIONAL(USE_DGEMMSY, test "$ax_have_dgemmsy" = "yes")
])
