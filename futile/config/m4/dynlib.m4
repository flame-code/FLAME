# -*- Autoconf -*-
#
# Copyright (c) 2016 BigDFT Group (Damien Caliste)
# All rights reserved.
#
# This file is part of the BigDFT software package. For license information,
# please see the COPYING file in the top-level directory of the BigDFT source
# distribution.
#

AC_DEFUN([AX_FC_BUILD_SHARED],
[dnl Try to find the flag to build shared objects.
  AC_LANG_PUSH(Fortran)
  AC_REQUIRE([AC_PROG_FC])
  AC_REQUIRE([AX_FLAG_PIC])

  AC_LANG_CONFTEST([
subroutine test()
  implicit none

  write(*,*) "test"  
end subroutine test])

  AC_MSG_CHECKING([for option to create shared libraries with $FC])
  ac_try='$FC $FCFLAGS $ax_fc_pic -c conftest.f90 -o conftest.o 1>&AC_FD_CC'
  if AC_TRY_EVAL(ac_try); then
    ac_try=""
  else
    echo "configure: failed program was:" >&AC_FD_CC
    cat conftest.f90 >&AC_FD_CC
    AC_MSG_WARN(Fortran compiler cannot compile subroutine.)
  fi

  f90_test_shared()
  {
    ac_try='$FC $FCFLAGS $[1] ${ax_fc_linker_wl}-soname=libtest.so.0 -o libtest.so.0.0.0 conftest.o 1>&AC_FD_CC'
    if AC_TRY_EVAL(ac_try); then
      ac_try=""
      ax_fc_build_shared=$[1]
    else
      ax_fc_build_shared=""
    fi
  }

  f90_test_shared "-shared"
  if test -z "$ax_fc_build_shared" ; then
   f90_test_shared "-qmkshrobj"
  fi
  if test -z "$ax_fc_build_shared" ; then
   f90_test_shared "-dynamic"
  fi
  if test -z "$ax_fc_build_shared" ; then
    ax_fc_build_shared="no"
    AC_MSG_WARN(Fortran compiler cannot create shared library.)
  else
    AC_SUBST([FC_BUILD_SHARED], [$ax_fc_build_shared])
    AC_MSG_RESULT([$ax_fc_build_shared])
  fi
  rm -f conftest.o libtest.so.0.0.0
  AM_CONDITIONAL([FC_AUTOBUILD_DYNAMIC_LIBS], [test "x$ax_fc_build_shared" != "xno"])
  AC_LANG_POP(Fortran)
])

AC_DEFUN([AX_CC_BUILD_SHARED],
[dnl Try to find the flag to build shared objects.
  AC_LANG_PUSH(C)
  AC_REQUIRE([AC_PROG_CC])
  AC_REQUIRE([AX_FLAG_PIC])

  AC_LANG_CONFTEST([AC_LANG_SOURCE([
extern int testint=5;

void mtrace_init(int val)
{
testint=val;
}])])

  AC_MSG_CHECKING([for option to create shared libraries with $CC])
  ac_try='$CC $CFLAGS $ax_cc_pic -c conftest.c -o conftest.o 1>&AC_FD_CC'
  if AC_TRY_EVAL(ac_try); then
    ac_try=""
  else
    echo "configure: failed program was:" >&AC_FD_CC
    cat conftest.c >&AC_FD_CC
    AC_MSG_WARN(C compiler cannot compile function.)
  fi

  cc_test_shared()
  {
    ac_try='$CC $CFLAGS $[1] ${ax_cc_linker_wl}-soname=libtest.so.0 -o libtest.so.0.0.0 conftest.o 1>&AC_FD_CC'
    if AC_TRY_EVAL(ac_try); then
      ac_try=""
      ax_cc_build_shared=$[1]
    else
      ax_cc_build_shared=""
    fi
  }

  cc_test_shared "-shared"
  if test -z "$ax_cc_build_shared" ; then
    cc_test_shared "-qmkshrobj"
  fi
  if test -z "$ax_cc_build_shared" ; then
    cc_test_shared "-dynamic"
  fi
  if test -z "$ax_cc_build_shared" ; then
    ax_cc_build_shared="no"
    AC_MSG_WARN(C compiler cannot create shared library.)
  else
    AC_SUBST([CC_BUILD_SHARED], [$ax_cc_build_shared])
    AC_MSG_RESULT([$ax_cc_build_shared])
  fi
  rm -f conftest.o libtest.so.0.0.0
  AM_CONDITIONAL([CC_AUTOBUILD_DYNAMIC_LIBS], [test "x$ax_cc_build_shared" != "xno"])
  AC_LANG_POP(C)
])

#link a given program in the chosen LD for testing if a given option is supported
#usage: $1 option name
#usage: $2 actual option to be tested
#usage: $3 linker for be used (FC, CC)
#usage: $4 program to test
#usage: $5 additional options to be passed (but not stored in the variable)

AC_DEFUN([AX_LD_OPT],
[dnl try to find the option flag to specify the rpath to the fortran compiler when linking.
  AC_REQUIRE([AX_FLAG_PIC])

  AC_MSG_CHECKING([for linker option to specify $2 to $[$3]])

  LDFLAGS_SVG=$LDFLAGS
  LDFLAGS="$LDFLAGS_SVG ${ax_fc_linker_wl}$2$5"
  AC_LINK_IFELSE([AC_LANG_SOURCE($4)], [ax_$3_$1="${ax_fc_linker_wl}$2"], [ax_$3_$1="no"])
  LDFLAGS=$LDFLAGS_SVG

  if test x"$ax_$3_$1" = x"no" ; then
    AC_MSG_WARN($[$3] linker does not accept option $2)
  else
    AC_MSG_RESULT([$ax_$3_$1])
    AC_SUBST([$3_$1], [$ax_$3_$1])
  fi
])

AC_DEFUN([AX_LINKER_OPTS],
[dnl try to find the option flag to specify various options when linking.
   AC_REQUIRE([AX_FLAG_PIC])
   #Fortran first
   AC_LANG_PUSH(Fortran)	
   #test for rpath
   AX_LD_OPT([RPATH],[$ax_fc_linker_rpath],[FC],
   [program test
write(*,*) "test"
end program test],[$[PWD]])
   #test for export symbols
   AX_LD_OPT([EXPORTS],[$ax_fc_linker_export_symbols],[FC],
   [program test
write(*,*) "test"
end program test],[])
   AC_LANG_POP(Fortran)
   AC_LANG_PUSH(C)
   AX_LD_OPT([RPATH],[$ax_cc_linker_rpath],[CC],
   [int main(int argc, char **argv)
{
  return 0;
}],[$[PWD]])
   #test for export symbols
   AX_LD_OPT([EXPORTS],[$ax_cc_linker_export_symbols],[CC],
   [int main(int argc, char **argv)
{
  return 0;
}],[])
   
   AC_LANG_POP(C)
])


AC_DEFUN([AX_FC_RPATH],
[dnl try to find the option flag to specify the rpath to the fortran compiler when linking.
  AC_REQUIRE([AX_FLAG_PIC])
  AC_LANG_PUSH(Fortran)

  f90_test_rpath()
  {
    LDFLAGS_SVG=$LDFLAGS
    LDFLAGS="$LDFLAGS_SVG $[1]$[PWD]"
    AC_LINK_IFELSE(AC_LANG_SOURCE([
program test
write(*,*) "test"
end program test
]), [ax_fc_rpath="$[1]"], [ax_fc_rpath="no"])
    LDFLAGS=$LDFLAGS_SVG
  }

  AC_MSG_CHECKING([for option to specify rpath to $FC])
  f90_test_rpath "${ax_fc_linker_wl}${ax_fc_linker_rpath}"
  if test x"$ax_fc_rpath" = x"no" ; then
    AC_MSG_WARN(Fortran compiler expand rpath.)
  else
    AC_SUBST([FC_RPATH], [$ax_fc_rpath])
    AC_MSG_RESULT([$ax_fc_rpath])
  fi
  
  AC_LANG_POP(Fortran)
])

AC_DEFUN([AX_CC_RPATH],
[dnl try to find the option flag to specify the rpath to the c compiler when linking.
  AC_LANG_PUSH(C)

  c_test_rpath()
  {
    LDFLAGS_SVG=$LDFLAGS
    LDFLAGS="$LDFLAGS_SVG $[1]$[PWD]"
    AC_LINK_IFELSE([AC_LANG_SOURCE([
int main(int argc, char **argv)
{
  return 0;
}
])], [ax_cc_rpath="$[1]"], [ax_cc_rpath="no"])
    LDFLAGS=$LDFLAGS_SVG
  }

  AC_MSG_CHECKING([for option to specify rpath to $CC])
  c_test_rpath "${ax_cc_linker_wl}${ax_fc_linker_rpath}"
  if test x"$ax_cc_rpath" = x"no" ; then
    AC_MSG_WARN(C compiler expand rpath.)
  else
    AC_SUBST([CC_RPATH], [$ax_cc_rpath])
    AC_MSG_RESULT([$ax_cc_rpath])
  fi
  
  AC_LANG_POP(C)
])

# AX_DYNAMIC_LIBRARIES([DEFAULT = "no"])
AC_DEFUN([AX_DYNAMIC_LIBRARIES],
[dnl Produce dynamic libraries and executables.
  AC_ARG_ENABLE(dynamic-libraries, AS_HELP_STRING([--enable-dynamic-libraries],
                                                 [Build dynamical libraries (disabled by default).]),
                [ax_build_dynamic=$enableval; ax_user_input="yes"],
                [ax_build_dynamic=m4_default([$1], ["no"]); ax_user_input="no"])
  dnl Test for library building tools.
  if test x"$ax_build_dynamic" = x"yes" ; then
    AC_REQUIRE([AX_FLAG_PIC])
    AC_REQUIRE([AX_FC_BUILD_SHARED])
    AC_REQUIRE([AX_CC_BUILD_SHARED])
    dnl AC_REQUIRE([AX_FC_RPATH])
    dnl AC_REQUIRE([AX_CC_RPATH])
    AC_REQUIRE([AX_LINKER_OPTS])
    if test "$ax_fc_pic" = "no" -o "$ax_cc_pic" = "no" -o "$ax_fc_build_shared" = "no" -o "$ax_FC_RPATH" = "no"  -o "$ax_FC_EXPORTS" = "no" ; then
      AC_MSG_WARN(["Dynamic libraries disabled (see reason before)."])
      ax_build_dynamic=no
    else
      dnl try to link with dependencies to see if its working.
      dnl dependencies may not have been compiled with PIC.
      ax_dynamic_deps=yes
      ax_dyndeps_libs=m4_default([$3], [""])
      if test -n "$ax_dyndeps_libs" -a -n "$2" ; then
        AC_MSG_CHECKING([for dynamic linking of dependencies])
        AC_LANG_PUSH(Fortran)
        LDFLAGS_SVG=$LDFLAGS
        LDFLAGS="-shared $LDFLAGS"
        LIBS_SVG=$LIBS
        LIBS="$ax_dyndeps_libs $LIBS"
        FCFLAGS_SVG="$FCFLAGS"
        FCFLAGS="$ax_fc_pic $FCFLAGS"
        AC_LINK_IFELSE([AC_LANG_SOURCE([
subroutine testlib()
call $2
end subroutine testlib
        ])], ax_dynamic_deps=yes, ax_dynamic_deps=no)
        AC_LANG_POP([Fortran])
        FCFLAGS=$FCFLAGS_SVG
        LIBS=$LIBS_SVG
        LDFLAGS=$LDFLAGS_SVG
        AC_MSG_RESULT([$ax_dynamic_deps])
      fi

      if test "$ax_dynamic_deps" = "yes" ; then
        case $CFLAGS   in *"$ax_cc_pic"*) ;; *) CFLAGS="$CFLAGS $ax_cc_pic";; esac
        case $CXXFLAGS in *"$ax_cc_pic"*) ;; *) CXXFLAGS="$CXXFLAGS $ax_cc_pic";; esac
        case $FCFLAGS  in *"$ax_fc_pic"*) ;; *) FCFLAGS="$FCFLAGS $ax_fc_pic";; esac
    
dnl         eval "set x $ac_configure_args"
dnl         shift
dnl         ac_configure_args=
dnl         for ac_arg ; do
dnl           case $ac_arg in
dnl             CFLAGS=* | CXXFLAGS=* | FCFLAGS=*) ;;
dnl             *) ac_configure_args="$ac_configure_args '$ac_arg'"
dnl           esac
dnl         done
dnl         ac_configure_args="$ac_configure_args 'CFLAGS=$CFLAGS' 'CXXFLAGS=$CXXFLAGS' 'FCFLAGS=$FCFLAGS'"
      else
        ax_build_dynamic="no"
      fi
    fi

    dnl Raise error if asked for dynamic but cannot proceed.
    if test "$ax_user_input" == "yes" -a "$ax_build_dynamic" != "yes" ; then
      AC_MSG_ERROR(["Dynamic build is not possible."])
    fi
  else
    ax_build_dynamic="no"
  fi

  AM_CONDITIONAL([BUILD_DYNAMIC_LIBS], [test "x$ax_build_dynamic" = "xyes"])
])
