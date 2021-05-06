# -*- Autoconf -*-
#
# Copyright (c) 2014 BigDFT Group (Damien Caliste)
# All rights reserved.
#
# This file is part of the BigDFT software package. For license information,
# please see the COPYING file in the top-level directory of the BigDFT source
# distribution.
#

AC_DEFUN([AX_YAML],
[dnl C-level YAML support.
  AC_ARG_WITH([yaml-path],
              AS_HELP_STRING([--with-yaml-path], [give a path to find libyaml.]),
              [ax_yaml_path=$withval])
  
  CPPFLAGS_SVG="$CPPFLAGS"

  if test x"$ax_yaml_path" != x"" ; then
     CPPFLAGS="$CPPFLAGS -I$ax_yaml_path/include"
  fi

  AC_LANG_PUSH(C)
  AC_CHECK_HEADER([yaml.h],
                  [ax_have_yaml="yes"],
                  [ax_have_yaml="no"])
  if test x"$ax_have_yaml" = x"yes"; then
     if test x"$ax_yaml_path" != x"" ; then
        LIB_YAML_CFLAGS="-I$ax_yaml_path/include"
     else
        for path in ${C_INCLUDE_PATH//:/ }; do
           LIB_YAML_CFLAGS="$LIB_YAML_CFLAGS -I$path"
        done
     fi     
  else
     AC_MSG_ERROR([libyaml is not available, install YAML and provide path --with-yaml-path.])
  fi
  CPPFLAGS="$CPPFLAGS_SVG"

  LDFLAGS_SVG="$LDFLAGS"
  LDFLAGS="$LDFLAGS -L$ax_yaml_path/lib"
  AC_CHECK_LIB([yaml], [yaml_parser_parse],
               [ax_have_yaml=yes], [ax_have_yaml=no])
  if test x"$ax_have_yaml" = x"yes"; then
     if test x"$ax_yaml_path" != x"/usr" ; then
        LIB_YAML_LIBS="-L$ax_yaml_path/lib "
     fi
     LIB_YAML_LIBS=$LIB_YAML_LIBS"-lyaml"
  else
     AC_MSG_ERROR([libyaml is not available, install YAML and provide path --with-yaml-path.])
  fi
  AC_LANG_POP(C)
  LDFLAGS="$LDFLAGS_SVG"

  AC_SUBST(LIB_YAML_CFLAGS)
  AC_SUBST(LIB_YAML_LIBS)
])
  
AC_DEFUN([AX_PYYAML],
[dnl Test for libXC
  dnl Python itself (set to ":" if not present).
  AC_REQUIRE([AM_PATH_PYTHON])

  AC_ARG_WITH([pyyaml-path],
              AS_HELP_STRING([--with-pyyaml-path], [give a path to find YAML python module.]),
              [AX_PYYAML_PATH=$withval], [AX_PYYAML_PATH=$pyexecdir])
  ax_have_pyyaml="no"
  if test "$PYTHON" != ":" ; then
    AC_MSG_CHECKING([for PyYAML and CLoader from: "$AX_PYYAML_PATH"])
   topath=$(echo $AX_PYYAML_PATH | sed  s/\$\{exec_prefix\}//g)
   cat > pytest << EOF
try:
  import yaml

  try:
    l = yaml.CLoader
    print ("0")
  except:
    print ("2")
except:
  print ("1")
EOF
    eval=$(PYTHONPATH=$prefix$topath $PYTHON pytest)
    rm -f pytestiii
    case "$eval" in
      "0")
        AC_MSG_RESULT([$AX_PYYAML_PATH])
        ax_have_pyyaml="yes"
        ;;
      "2")
        AC_MSG_RESULT([no])
        AC_MSG_WARN([CLoader not available, check that the PyYAML python module was linked with the C-level YAML implementation.])
        ax_have_pyyaml="yes"
        ;;
      *)
        AC_MSG_RESULT([no])
        AC_MSG_WARN([PyYAML not available, use --with-pyyaml-path to specify the module location.])
        ax_have_pyyaml="no"
        ;;
    esac
  fi
  AM_CONDITIONAL(HAVE_PYYAML, test x"$ax_have_pyyaml" = x"yes")
  AC_SUBST(AX_PYYAML_PATH)
])
