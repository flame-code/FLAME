# Define a macro to test if documantation utilities are is available 
#
# Copyright (c) 2011-2013 BigDFT Group (Damien Caliste)
# All rights reserved.
#
# This file is part of the BigDFT software package. For license information,
# please see the COPYING file in the top-level directory of the BigDFT source
# distribution.
AC_DEFUN([AX_DOC],
[
 AC_CHECK_PROG(ac_have_doxygen, [doxygen], [yes], [no])
 AM_CONDITIONAL(BUILD_DEVDOC, test x"$ac_have_doxygen" = x"yes" -a x"$ac_devel_doc" = x"yes")
 dnl default installation directories
 docdir="${docdir}"
 AC_SUBST(docdir)
])


