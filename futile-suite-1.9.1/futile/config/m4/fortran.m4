# -*- Autoconf -*-
#
# Copyright (c) 2005-2008 ABINIT Group (Yann Pouillon)
# All rights reserved.
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#

#
# Fortran compilers support
#



# _ABI_CHECK_FC_ABSOFT(COMPILER)
# ------------------------------
#
# Checks whether the specified Fortran compiler is the ABSoft Fortran compiler.
# If yes, tries to determine its version number and sets the fc_type
# and fc_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_FC_ABSOFT],
[dnl Do some sanity checking of the arguments
 m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

 dnl AC_MSG_CHECKING([if we are using the ABSoft Fortran compiler])
 fc_info_string=`$1 -V 2> /dev/null`
 abi_result=`echo "${fc_info_string}" | grep '^Pro Fortran'`
 if test "${abi_result}" = ""; then
  abi_result="no"
  fc_info_string=""
  fc_type="UNKNOWN"
  fc_version="UNKNOWN"
 else
  AC_DEFINE([FC_ABSOFT],1,[Define to 1 if you are using the ABSOFT Fortran compiler])
  fc_type="absoft"
  fc_version=`echo "${abi_result}" | sed -e 's/Pro Fortran //'`
  if test "${fc_version}" = "${abi_result}"; then
   fc_version="UNKNOWN"
  fi
  abi_result="yes"
 fi
 dnl AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_FC_ABSOFT



# _ABI_CHECK_FC_COMPAQ(COMPILER)
# ------------------------------
#
# Checks whether the specified Fortran compiler is the COMPAQ Fortran compiler.
# If yes, tries to determine its version number and sets the fc_type
# and fc_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_FC_COMPAQ],
[dnl Do some sanity checking of the arguments
 m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

 dnl AC_MSG_CHECKING([if we are using the Compaq Fortran compiler])
 fc_info_string=`$1 -version 2>&1 | sed -e 's/^	//' | grep '^Compaq Fortran Compiler'`
 abi_result="${fc_info_string}"
 if test "${abi_result}" = ""; then
  fc_info_string=`$1 -version 2>&1 | sed -e 's/^	//' | grep '^HP Fortran Compiler'`
  abi_result="${fc_info_string}"
 fi
 if test "${abi_result}" = ""; then
  abi_result="no"
  fc_info_string=""
  fc_type="UNKNOWN"
  fc_version="UNKNOWN"
 else
  AC_DEFINE([FC_COMPAQ],1,[Define to 1 if you are using the COMPAQ Fortran compiler])
  fc_type="compaq"
  fc_version=`echo "${abi_result}" | sed -e 's/.* V//;s/-.*//'`
  if test "${fc_version}" = "${abi_result}"; then
   fc_version="UNKNOWN"
  fi
  abi_result="yes"
 fi
 dnl AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_FC_COMPAQ



# _ABI_CHECK_FC_FUJITSU(COMPILER)
# -------------------------------
#
# Checks whether the specified Fortran compiler is the Fujitsu Fortran compiler.
# If yes, tries to determine its version number and sets the fc_type
# and fc_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_FC_FUJITSU],
[dnl Do some sanity checking of the arguments
 m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

 dnl AC_MSG_CHECKING([if we are using the Fujitsu Fortran compiler])
 fc_info_string=`$1 -V 2> /dev/null`
 abi_result=`echo "${fc_info_string}" | grep '^Fujitsu Fortran'`
 if test "${abi_result}" = ""; then
  abi_result="no"
  fc_info_string=""
  fc_type="UNKNOWN"
  fc_version="UNKNOWN"
 else
  AC_DEFINE([FC_FUJITSU],1,[Define to 1 if you are using the Fujitsu Fortran compiler])
  fc_type="fujitsu"
  fc_version=`echo "${abi_result}" | sed -e 's/.*Driver //;s/ .*//'`
  if test "${fc_version}" = "${abi_result}"; then
   fc_version="UNKNOWN"
  fi
  abi_result="yes"
 fi
 dnl AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_FC_FUJITSU



# _ABI_CHECK_FC_G95(COMPILER)
# ---------------------------
#
# Checks whether the specified Fortran compiler is the G95 Fortran compiler.
# If yes, tries to determine its version number and sets the fc_type
# and fc_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_FC_G95],
[dnl Do some sanity checking of the arguments
 m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

 dnl AC_MSG_CHECKING([if we are using the G95 Fortran compiler])
 fc_info_string=`$1 --version 2>&1`
 abi_result=`echo "${fc_info_string}" | grep '^G95'`
 if test "${abi_result}" = ""; then
  abi_result="no"
  fc_info_string=""
  fc_type="UNKNOWN"
  fc_version="UNKNOWN"
 else
  AC_DEFINE([FC_G95],1,[Define to 1 if you are using the G95 Fortran compiler])
  fc_type="g95"
  fc_version=`echo ${abi_result} | sed -e 's/.*GCC //; s/ .*//'`
  if test "${fc_version}" = "${abi_result}"; then
   fc_version="UNKNOWN"
  fi
  abi_result="yes"
 fi
 dnl AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_FC_G95



# _ABI_CHECK_FC_GCC(COMPILER)
# ---------------------------
#
# Checks whether the specified Fortran compiler is the GCC Fortran compiler.
# If yes, tries to determine its version number and sets the fc_type
# and fc_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_FC_GCC],
[dnl Do some sanity checking of the arguments
 m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

 dnl AC_MSG_CHECKING([if we are using the GCC Fortran compiler])
 fc_info_string=`$1 --version 2>&1`
 abi_result=`echo "${fc_info_string}" | grep '^GNU Fortran'`
 if test "${abi_result}" = ""; then
  abi_result="no"
  fc_info_string=""
  fc_type="UNKNOWN"
  fc_version="UNKNOWN"
 else
  AC_DEFINE([FC_GCC],1,[Define to 1 if you are using the GNU Fortran compiler])
  fc_type="gcc"
  fc_version=`echo ${abi_result} | sed -e 's/.*(GCC) //; s/.*GCC //; s/ .*//'`
  abi_result="yes"
 fi
 dnl AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_FC_GCC



# _ABI_CHECK_FC_HITACHI(COMPILER)
# -------------------------------
#
# Checks whether the specified Fortran compiler is the Hitachi Fortran compiler.
# If yes, tries to determine its version number and sets the fc_type
# and fc_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_FC_HITACHI],
[dnl Do some sanity checking of the arguments
 m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

 dnl AC_MSG_CHECKING([if we are using the Hitachi Fortran compiler])
 fc_info_string=`$1 -V 2> /dev/null`
 abi_result=`echo "${fc_info_string}" | grep '^Hitachi Fortran'`
 if test "${abi_result}" = ""; then
  abi_result="no"
  fc_info_string=""
  fc_type="UNKNOWN"
  fc_version="UNKNOWN"
 else
  AC_DEFINE([FC_HITACHI],1,[Define to 1 if you are using the Hitachi Fortran compiler])
  fc_type="hitachi"
  fc_version=`echo "${abi_result}" | sed -e 's/.*Driver //;s/ .*//'`
  if test "${fc_version}" = "${abi_result}"; then
   fc_version="UNKNOWN"
  fi
  abi_result="yes"
 fi
 dnl AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_FC_HITACHI



# _ABI_CHECK_FC_IBM(COMPILER)
# ---------------------------
#
# Checks whether the specified Fortran compiler is the IBM XL Fortran compiler.
# If yes, tries to determine its version number and sets the fc_type
# and fc_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_FC_IBM],
[dnl Do some sanity checking of the arguments
 m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

 dnl AC_MSG_CHECKING([if we are using the IBM XL Fortran compiler])
 fc_info_string=`$1 -qversion 2>&1`
 fc_garbage=`$1 -qversion 2>&1 | wc -l | sed -e 's/ //g'`
 abi_result=`echo "${fc_info_string}" | grep 'IBM(R) XL Fortran'`
 if test "${abi_result}" = ""; then
  abi_result=`echo "${fc_info_string}" | grep 'IBM XL Fortran'`
 fi
 if test "${abi_result}" = ""; then
  abi_result="no"
  fc_info_string=""
  fc_type="UNKNOWN"
  fc_version="UNKNOWN"
  if test "${fc_garbage}" -gt 50; then
   AC_DEFINE([FC_IBM],1,[Define to 1 if you are using the IBM XL Fortran compiler])
   fc_type="ibm"
   fc_version="UNKNOWN"
   abi_result="yes"
  fi
 else
  AC_DEFINE([FC_IBM],1,[Define to 1 if you are using the IBM XL Fortran compiler])
  fc_type="ibm"
  fc_version=`echo "${abi_result}" | sed -e 's/.* V\([[0-9\.]]*\) .*/\1/'`
  if test "${fc_version}" = "${abi_result}"; then
   fc_version="UNKNOWN"
  fi
  abi_result="yes"
 fi
 dnl AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_FC_IBM



# _ABI_CHECK_FC_INTEL(COMPILER)
# -----------------------------
#
# Checks whether the specified Fortran compiler is the Intel Fortran compiler.
# If yes, tries to determine its version number and sets the fc_type
# and fc_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_FC_INTEL],
[dnl Do some sanity checking of the arguments
 m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

 dnl AC_MSG_CHECKING([if we are using the Intel Fortran compiler])
 fc_info_string=`$1 -v -V 2>&1 | sed -e '/^ifc: warning/d'`
 abi_result=`echo "${fc_info_string}" | grep '^Intel(R) Fortran'`
 if test "${abi_result}" = ""; then
  abi_result="no"
  fc_info_string=""
  fc_type="UNKNOWN"
  fc_version="UNKNOWN"
 else
  AC_DEFINE([FC_INTEL],1,[Define to 1 if you are using the Intel Fortran compiler])
  fc_type="intel"
  fc_version=`echo "${fc_info_string}" | grep '^Version' | sed -e 's/Version //;s/ .*//;s/ //g' | head -n 1`
  if test "${fc_version}" = ""; then
   fc_version="UNKNOWN"
  fi
  abi_result="yes"
 fi
 dnl AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_FC_INTEL



# _ABI_CHECK_FC_MIPSPRO(COMPILER)
# -------------------------------
#
# Checks whether the specified Fortran compiler is the MIPSpro Fortran
# compiler.
# If yes, tries to determine its version number and sets the fc_type
# and fc_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_FC_MIPSPRO],
[dnl Do some sanity checking of the arguments
 m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

 dnl AC_MSG_CHECKING([if we are using the MIPSpro Fortran compiler])
 fc_info_string=`$1 -version 2>&1 | sed -e '/^$/d'`
 abi_result=`echo "${fc_info_string}" | grep '^MIPSpro'`
 if test "${abi_result}" = ""; then
  abi_result="no"
  fc_info_string=""
  fc_type="UNKNOWN"
  fc_version="UNKNOWN"
 else
  AC_DEFINE([FC_MIPSPRO],1,[Define to 1 if you are using the MIPSpro Fortran compiler])
  fc_type="mipspro"
  fc_version=`echo "${abi_result}" | sed -e 's/.*Version //'`
  if test "${fc_version}" = "${abi_result}"; then
   fc_version="UNKNOWN"
  fi
  abi_result="yes"
 fi
 dnl AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_FC_MIPSPRO



# _ABI_CHECK_FC_OPEN64(COMPILER)
# ------------------------------
#
# Checks whether the specified Fortran compiler is the Open64
# Fortran compiler.
# If yes, tries to determine its version number and sets the fc_type
# and fc_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_FC_OPEN64],
[dnl Do some sanity checking of the arguments
 m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

 dnl AC_MSG_CHECKING([if we are using the PathScale Fortran compiler])
 fc_info_string=`$1 --version 2>&1`
 abi_result=`echo "${fc_info_string}" | grep '^Open64'`
 if test "${abi_result}" = ""; then
  abi_result="no"
  fc_info_string=""
  fc_type="UNKNOWN"
  fc_version="UNKNOWN"
 else
  AC_DEFINE([FC_OPEN64],1,[Define to 1 if you are using the Open64 Fortran compiler])
  fc_type="open64"
  fc_version=`echo "${abi_result}" | sed -e 's/.* Version //; s/ .*//'`
  if test "${fc_version}" = "${abi_result}"; then
   fc_version="UNKNOWN"
  fi
  abi_result="yes"
 fi
 dnl AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_FC_OPEN64



# _ABI_CHECK_FC_PATHSCALE(COMPILER)
# ---------------------------------
#
# Checks whether the specified Fortran compiler is the PathScale
# Fortran compiler.
# If yes, tries to determine its version number and sets the fc_type
# and fc_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_FC_PATHSCALE],
[dnl Do some sanity checking of the arguments
 m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

 dnl AC_MSG_CHECKING([if we are using the PathScale Fortran compiler])
 fc_info_string=`$1 -version 2>&1`
 abi_result=`echo "${fc_info_string}" | grep '^PathScale'`
 if test "${abi_result}" = ""; then
  abi_result="no"
  fc_info_string=""
  fc_type="UNKNOWN"
  fc_version="UNKNOWN"
 else
  AC_DEFINE([FC_PATHSCALE],1,[Define to 1 if you are using the PathScale Fortran compiler])
  fc_type="pathscale"
  fc_version=`echo "${abi_result}" | sed -e 's/.* Version //; s/ .*//'`
  if test "${fc_version}" = "${abi_result}"; then
   fc_version="UNKNOWN"
  fi
  abi_result="yes"
 fi
 dnl AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_FC_PATHSCALE



# _ABI_CHECK_FC_PGI(COMPILER)
# ---------------------------
#
# Checks whether the specified Fortran compiler is the Portland Group
# Fortran compiler.
# If yes, tries to determine its version number and sets the fc_type
# and fc_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_FC_PGI],
[dnl Do some sanity checking of the arguments
 m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

 dnl AC_MSG_CHECKING([if we are using the Portland Group Fortran compiler])
 fc_info_string=`$1 -V 2>&1 | sed -e '/^$/d'`
 abi_result=`echo "${fc_info_string}" | grep '^pgf9[[05]]' | grep -v 'No files to process'`
 if test "${abi_result}" = ""; then
  abi_result="no"
  fc_info_string=""
  fc_type="UNKNOWN"
  fc_version="UNKNOWN"
 else
  AC_DEFINE([FC_PGI],1,[Define to 1 if you are using the Portland Group Fortran compiler])
  fc_type="pgi"
  fc_version=`echo "${abi_result}" | sed -e 's/^pgf9[[05]] //' | sed -e 's/-.*//'`
  if test "${fc_version}" = "${abi_result}"; then
   fc_version="UNKNOWN"
  else
   if test "${fc_version}" = "6.0"; then
        AC_DEFINE([FC_PGI6],1,[Define to 1 if you are using the Portland Group Fortran compiler version 6])
   fi
  fi
  abi_result="yes"
 fi
 dnl AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_FC_PGI



# _ABI_CHECK_FC_SUN(COMPILER)
# ---------------------------
#
# Checks whether the specified Fortran compiler is the Sun WorkShop Fortran compiler.
# If yes, tries to determine its version number and sets the fc_type
# and fc_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_FC_SUN],
[dnl Do some sanity checking of the arguments
 m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

 dnl AC_MSG_CHECKING([if we are using the Sun WorkShop Fortran compiler])
 fc_info_string=`$1 -V 2>&1 | head -n 1`
 abi_result=`echo "${fc_info_string}" | grep 'Sun' | grep 'Fortran 95'`
 if test "${abi_result}" = ""; then
  abi_result="no"
  fc_info_string=""
  fc_type="UNKNOWN"
  fc_version="UNKNOWN"
 else
  AC_DEFINE([FC_SUN],1,[Define to 1 if you are using the Sun WorkShop])
  fc_type="sun"
  fc_version=`echo "${abi_result}" | sed -e 's/.* Fortran 95 //;s/ .*//'`
  if test "${fc_version}" = "${abi_result}" -o "${fc_version}" = ""; then
   fc_version="UNKNOWN"
  fi
  abi_result="yes"
 fi
 dnl AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_FC_SUN



 ##############################################################################



# _ABI_CHECK_FC_EXIT()
# --------------------
#
# Checks whether the Fortran compiler supports the exit() subroutine.
#
AC_DEFUN([_ABI_CHECK_FC_EXIT],
[dnl Init
 fc_has_exit="no"

 AC_MSG_CHECKING([whether the Fortran compiler accepts exit()])

 dnl Try to compile a program calling exit()
 AC_LANG_PUSH([Fortran])
 AC_LINK_IFELSE([AC_LANG_PROGRAM([],
  [[
      call exit(1)
  ]])], [fc_has_exit="yes"])
 AC_LANG_POP()

 if test "${fc_has_exit}" = "yes"; then
  AC_DEFINE([HAVE_FC_EXIT],1,
   [Define to 1 if your Fortran compiler supports exit()])
 fi

 AC_MSG_RESULT(${fc_has_exit})
]) # _ABI_CHECK_FC_EXIT



# _ABI_CHECK_FC_FLUSH()
# ---------------------
#
# Checks whether the Fortran compiler supports the flush() subroutine.
#
AC_DEFUN([_ABI_CHECK_FC_FLUSH],
[dnl Init
 fc_has_flush="no"

 AC_MSG_CHECKING([whether the Fortran compiler accepts flush()])

 dnl Try to compile a program calling flush()
 AC_LANG_PUSH([Fortran])
 AC_LINK_IFELSE([AC_LANG_PROGRAM([],
  [[
      call flush()
  ]])], [fc_has_flush="yes"])
 AC_LANG_POP()

 if test "${fc_has_flush}" = "yes"; then
  AC_DEFINE([HAVE_FC_FLUSH],1,
   [Define to 1 if your Fortran compiler supports flush()])
 fi

 AC_MSG_RESULT(${fc_has_flush})
]) # _ABI_CHECK_FC_FLUSH



# ABI_PROG_FC()
# -------------
#
# Tries to determine which type of Fortran compiler is installed.
#
AC_DEFUN([ABI_PROG_FC],
[dnl Init
 if test "${fc_type}" = ""; then
  fc_type="UNKNOWN"
 fi
 if test "${fc_version}" = ""; then
  fc_version="UNKNOWN"
 fi
 fc_wrap="no"

 dnl Determine Fortran compiler type (the order is important)
 AC_MSG_CHECKING([which type of Fortran compiler we have])

 if test "${fc_type}" = "UNKNOWN"; then
  _ABI_CHECK_FC_G95(${FC})
 fi
 if test "${fc_type}" = "UNKNOWN"; then
  _ABI_CHECK_FC_GCC(${FC})
 fi
 if test "${fc_type}" = "UNKNOWN"; then
  _ABI_CHECK_FC_INTEL(${FC})
 fi
 if test "${fc_type}" = "UNKNOWN"; then
  _ABI_CHECK_FC_PATHSCALE(${FC})
 fi
 if test "${fc_type}" = "UNKNOWN"; then
  _ABI_CHECK_FC_PGI(${FC})
 fi
 if test "${fc_type}" = "UNKNOWN"; then
  _ABI_CHECK_FC_COMPAQ(${FC})
 fi
 if test "${fc_type}" = "UNKNOWN"; then
  _ABI_CHECK_FC_ABSOFT(${FC})
 fi
 if test "${fc_type}" = "UNKNOWN"; then
  _ABI_CHECK_FC_MIPSPRO(${FC})
 fi
 if test "${fc_type}" = "UNKNOWN"; then
  _ABI_CHECK_FC_OPEN64(${FC})
 fi
 if test "${fc_type}" = "UNKNOWN"; then
  _ABI_CHECK_FC_FUJITSU(${FC})
 fi
 if test "${fc_type}" = "UNKNOWN"; then
  _ABI_CHECK_FC_SUN(${FC})
 fi
 if test "${fc_type}" = "UNKNOWN"; then
  _ABI_CHECK_FC_HITACHI(${FC})
 fi
 dnl Always keep that one at the end
 if test "${fc_type}" = "UNKNOWN"; then
  _ABI_CHECK_FC_IBM(${FC})
 fi

 dnl Fall back to generic when detection fails
 if test "${fc_type}" = "UNKNOWN"; then
  fc_type="generic"
  fc_version="0.0"
 fi

 dnl Normalise Fortran compiler version
 fc_version=`echo ${fc_version} | cut -d. -f1-2`

 dnl Display final result
 AC_MSG_RESULT([${fc_type} ${fc_version}])

 dnl Schedule compiler info for substitution
 AC_SUBST(fc_type)
 AC_SUBST(fc_version)
 AC_SUBST(fc_wrap)

 dnl Further explore compiler peculiarities
 _ABI_CHECK_FC_EXIT
 _ABI_CHECK_FC_FLUSH
]) # ABI_PROG_FC


# Define a macro to test Fortran2003 implementation.
#
# Copyright (c) 2011-2011 BigDFT Group (Damien Caliste)
# All rights reserved.
#
# This file is part of the BigDFT software package. For license information,
# please see the COPYING file in the top-level directory of the BigDFT source
# distribution.
AC_DEFUN([AX_FC_F2003],
[
  AC_LANG_PUSH(Fortran)
  AC_REQUIRE([AC_PROG_FC])

  dnl We start with get_command_argument().
  AC_MSG_CHECKING([for get_command_argument() in Fortran.])

  AC_COMPILE_IFELSE([AC_LANG_SOURCE([
program test
  character(len = 128) :: arg

  call get_command_argument(1, arg)

end program test])],
  [ax_fc_get_command_argument="yes"], [ax_fc_get_command_argument="no"])
  if test x"$ax_fc_get_command_argument" == x"yes" ; then
    AC_DEFINE([HAVE_FC_GET_COMMAND_ARGUMENT], [1], [get_command_argument() can be used safely in Fortran])
  fi
  AM_CONDITIONAL([HAVE_FC_GET_COMMAND_ARGUMENT], [test x"$ax_fc_get_command_argument" == x"yes"])
  AC_MSG_RESULT([$ax_fc_get_command_argument])

  AC_LANG_POP(Fortran)
])

# Define a macro to test module output of the fortran compiler.
#
# Copyright (c) 2011-2011 BigDFT Group (Damien Caliste)
# All rights reserved.
#
# This file is part of the BigDFT software package. For license information,
# please see the COPYING file in the top-level directory of the BigDFT source
# distribution.

AC_DEFUN([AX_FC_MOD],
[
  AC_MSG_CHECKING([for module output in Fortran])

  AC_LANG_PUSH(Fortran)
  AC_REQUIRE([AC_PROG_FC])

  ax_fc_mod_compile=no
  AC_COMPILE_IFELSE([AC_LANG_SOURCE([
module modtest
  integer, public :: value
end module modtest
])], [ax_fc_mod_compile=yes],
   [AC_MSG_FAILURE(Fortran compiler cannot compile modules.)])
  if test $ax_fc_mod_compile = "yes" ; then
    ax_fc_mod_name="unknown"
    if test -s modtest.mod ; then
      ax_fc_mod_ext="mod"
      ax_fc_mod_capitalize="no"
      ax_fc_mod_name="module"
      rm -f modtest.mod
    fi
    if test -s modtest.MOD ; then
      ax_fc_mod_ext="MOD"
      ax_fc_mod_capitalize="no"
      ax_fc_mod_name="module"
      rm -f modtest.MOD
    fi
    if test -s MODTEST.MOD ; then
      ax_fc_mod_ext="MOD"
      ax_fc_mod_capitalize="yes"
      ax_fc_mod_name="MODULE"
      rm -f MODTEST.MOD
    fi
    if test -s MODTEST.mod ; then
      ax_fc_mod_ext="mod"
      ax_fc_mod_capitalize="yes"
      ax_fc_mod_name="MODULE"
      rm -f MODTEST.mod
    fi
    if test $ax_fc_mod_name = "unknown" ; then
       AC_MSG_WARN(Unknown module naming scheme for Fortran compiler.)
       ax_fc_mod_capitalize="no"
       ax_fc_mod_ext=""
    fi  
  fi
  AC_MSG_RESULT([$ax_fc_mod_name.$ax_fc_mod_ext])
])

# Define a macro to test the C binding of a Fortran pointer.
#
# Copyright (c) 2012-2012 BigDFT Group (Damien Caliste)
# All rights reserved.
#
# This file is part of the BigDFT software package. For license information,
# please see the COPYING file in the top-level directory of the BigDFT source
# distribution.
AC_DEFUN([AX_FC_POINTER],
[
  f90_pointer_test()
  {
  cat >pttest.f90 <<EOF
subroutine pt_test(shift, size)
  implicit none
  integer, intent(out) :: shift, size

  type pt_type
     double precision, dimension($[2]), pointer :: pt
  end type pt_type
  type(pt_type) :: pt(2)
  interface
     subroutine inqPt(pt1, pt2, start1, shift, size)
       double precision, dimension($[2]), pointer :: pt1, pt2
       double precision, intent(in) :: start1
       integer, intent(out) :: shift, size
     end subroutine inqPt
  end interface

  allocate(pt(1)%pt($[3]))
  call inqPt(pt(1)%pt, pt(2)%pt, pt(1)%pt($[4]), shift, size)
  deallocate(pt(1)%pt)
end subroutine pt_test
EOF
  }

  test_dim()
  {
  AC_MSG_CHECKING([for $[1] pointer structure in Fortran])

  AC_LANG_PUSH(Fortran)
  AC_REQUIRE([AC_PROG_FC])

  f90_pointer_test "$[1]" "$[2]" "$[3]" "$[4]"

  ac_try='$FC $FCFLAGS -c pttest.f90 1>&AC_FD_CC'
  if AC_TRY_EVAL(ac_try); then
    ac_try=""
  else
    echo "configure: failed program was:" >&AC_FD_CC
    cat intsizetest.f90 >&AC_FD_CC
    AC_MSG_FAILURE(Fortran compiler cannot compile subroutine.)
  fi

  AC_LANG_PUSH(C)
  AC_REQUIRE([AC_PROG_CC])
  
  LIBS_SVG="$LIBS"
  LIBS="pttest.o $LIBS $FCLIBS"
  AC_FC_FUNC([pt_test])
  AC_FC_FUNC([inqPt])
  AC_RUN_IFELSE([AC_LANG_PROGRAM([
#include <stdio.h>

void $inqPt(void **pt1, void **pt2, double *start1, int *shift, int *size)
{
  unsigned int i;

  for (i = 0; i < 20 && pt1[[i]] != (void*)start1; i++);
  *shift = (int)i;
  *size = (int)(((long)pt2 - (long)pt1) / sizeof(void*));
}],
[unsigned size, shift;

  $pt_test(&shift, &size);
  if (shift > size)
    return 1;
  fprintf(stdout, "%d %d\n", shift, size);])], [ax_fc_run=`./conftest$EXEEXT`],
 [AC_MSG_WARN([C compiler cannot link Fortran and C or cannot find the pointer shift value.])
  ax_fc_run="0 0"],
 [AC_MSG_WARN([Cross compiling, cannot test pointer length, using Gfortran values.])
  ax_fc_run="0 $[5]"])
  LIBS="$LIBS_SVG"
  rm -f pttest.*

  AC_LANG_POP(C)

  AC_LANG_POP(Fortran)
  }

  test_dim "1D" ":" "2" "1" "6"
  F90_1D_POINTER_SHIFT=`echo $ax_fc_run | cut -d' ' -f1`
  AC_SUBST(F90_1D_POINTER_SHIFT)
  F90_1D_POINTER_SIZE=`echo $ax_fc_run | cut -d' ' -f2`
  AC_SUBST(F90_1D_POINTER_SIZE)

  test_dim "2D" ":,:" "2,1" "1,1" "9"
  F90_2D_POINTER_SHIFT=`echo $ax_fc_run | cut -d' ' -f1`
  AC_SUBST(F90_2D_POINTER_SHIFT)
  F90_2D_POINTER_SIZE=`echo $ax_fc_run | cut -d' ' -f2`
  AC_SUBST(F90_2D_POINTER_SIZE)

  test_dim "3D" ":,:,:" "2,1,1" "1,1,1" "12"
  F90_3D_POINTER_SHIFT=`echo $ax_fc_run | cut -d' ' -f1`
  AC_SUBST(F90_3D_POINTER_SHIFT)
  F90_3D_POINTER_SIZE=`echo $ax_fc_run | cut -d' ' -f2`
  AC_SUBST(F90_3D_POINTER_SIZE)

  test_dim "4D" ":,:,:,:" "2,1,1,1" "1,1,1,1" "15"
  F90_4D_POINTER_SHIFT=`echo $ax_fc_run | cut -d' ' -f1`
  AC_SUBST(F90_4D_POINTER_SHIFT)
  F90_4D_POINTER_SIZE=`echo $ax_fc_run | cut -d' ' -f2`
  AC_SUBST(F90_4D_POINTER_SIZE)

  test_dim "5D" ":,:,:,:,:" "2,1,1,1,1" "1,1,1,1,1" "18"
  F90_5D_POINTER_SHIFT=`echo $ax_fc_run | cut -d' ' -f1`
  AC_SUBST(F90_5D_POINTER_SHIFT)
  F90_5D_POINTER_SIZE=`echo $ax_fc_run | cut -d' ' -f2`
  AC_SUBST(F90_5D_POINTER_SIZE)
])

# Define a macro to test Fortran compiler OpenMP flags.
#
# Copyright (c) 2014-2014 BigDFT Group (Damien Caliste)
# All rights reserved.
#
# This file is part of the BigDFT software package. For license information,
# please see the COPYING file in the top-level directory of the BigDFT source
# distribution.
AC_DEFUN([AX_FC_OPENMP],
[
  AC_LANG_PUSH(Fortran)
  AC_REQUIRE([AC_PROG_FC])
  dnl We start with get_command_argument().
  AC_MSG_CHECKING([for OpenMP flag in Fortran.])

  test_flag()
  {
  FCFLAGS_SVG=$FCFLAGS
  FCFLAGS=$FCFLAGS" $[1]"
  dnl First check if flag is valid for compiler.
  AC_COMPILE_IFELSE([AC_LANG_SOURCE([
program test
  write(*,*) "hello"
end program test
  ])],
  [ax_fc_openmp="yes"], [ax_fc_openmp="no"])
  if test x"$ax_fc_openmp" == x"yes" ; then
    dnl Now, test if the flag see the OpenMP error in the following code.
    AC_COMPILE_IFELSE([AC_LANG_SOURCE([
program test
    integer :: it, nt
!\$  integer :: ithread,omp_get_thread_num

    nt = 5
!\$omp parallel default (private) shared(nt)
!\$  ithread = omp_get_thread_num()
!\$omp CHOKE_ME
!\$omp do schedule(static,1)
    do it=1,nt,1
      write(*,*) ithread, it, nt
    end do
!\$omp enddo

end program test])],
    [ax_fc_openmp="no"], [ax_fc_openmp="$[1]"])
  fi
  FCFLAGS=$FCFLAGS_SVG
  }

  AC_ARG_WITH(openmp, AS_HELP_STRING([--with-openmp],
              [specify the flag to be used for OpenMP parts of the code (or autodetect if empty).]),
              ax_fc_openmp=$withval, ax_fc_openmp="auto")
  if test x"$ax_fc_openmp" == x"auto" -o x"$ax_fc_openmp" == x"yes" ; then
    test_flag "-openmp"
    if test x"$ax_fc_openmp" == x"no" ; then
      test_flag "-qsmp=omp"
      if test x"$ax_fc_openmp" == x"no" ; then
        test_flag "-fopenmp"
      fi
    fi
    if test x"$ax_fc_openmp" == x"no" ; then
      ax_fc_openmp=""
      ax_fc_openmp_msg="unknown"
    else
      ax_fc_openmp_msg="$ax_fc_openmp"
    fi
    AC_MSG_RESULT([$ax_fc_openmp_msg])
  elif test -n "$ax_fc_openmp" -a x"$ax_fc_openmp" != x"no" ; then
    test_flag "$ax_fc_openmp"
    ax_fc_openmp_msg="$ax_fc_openmp"
    AC_MSG_RESULT([$ax_fc_openmp_msg])
    if test x"$ax_fc_openmp" == x"no" ; then
      ax_fc_openmp=""
      AC_MSG_WARN([provided OpenMP flags are not working.])
    fi
  else
    ax_fc_openmp=""
    ax_fc_openmp_msg="not used"
    AC_MSG_RESULT([$ax_fc_openmp_msg])
  fi

  AC_SUBST(FCFLAGS_OPENMP, $ax_fc_openmp)

  AC_LANG_POP(Fortran)
])

# Define a macro to define the module installation directory.
#
# Copyright (c) 2016 BigDFT Group (Damien Caliste)
# All rights reserved.
#
# This file is part of the BigDFT software package. For license information,
# please see the COPYING file in the top-level directory of the BigDFT source
# distribution.
AC_DEFUN([AX_FC_MODULEDIR],
[
  AC_ARG_WITH(moduledir,
              AS_HELP_STRING([--with-moduledir],
                             [installation directory for module files [[INCLUDEDIR]]]),
                ac_moduledir=$withval, ac_moduledir="no")
  if test x"$ac_moduledir" != x"no" ; then
    moduledir=$ac_moduledir
  else
    moduledir=${includedir}
  fi
  AC_SUBST(moduledir)
])
