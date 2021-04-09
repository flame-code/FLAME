## this one is commonly used with AM_PATH_PYTHONDIR ...
dnl AM_CHECK_PYMOD(MODNAME [,SYMBOL [,ACTION-IF-FOUND [,ACTION-IF-NOT-FOUND]]])
dnl Check if a module containing a given symbol is visible to machine's python.
AC_DEFUN([AM_CHECK_PYMOD],
[AC_REQUIRE([AM_PATH_PYTHON])
py_mod_var=`echo $1['_']$2 | sed 'y%./+-%__p_%'`
AC_MSG_CHECKING(for ifelse([$2],[],,[$2 in ])python module $1)
AC_CACHE_VAL(py_cv_mod_$py_mod_var, [
ifelse([$2],[], [prog="
import sys
try:
        import $1
except ImportError:
        sys.exit(1)
except:
        sys.exit(0)
sys.exit(0)"], [prog="
import $1
$1.$2"])
if $PYTHON -c "$prog" 1>&AC_FD_CC 2>&AC_FD_CC
  then
    eval "py_cv_mod_$py_mod_var=yes"
  else
    eval "py_cv_mod_$py_mod_var=no"
  fi
])
py_val=`eval "echo \`echo '$py_cv_mod_'$py_mod_var\`"`
if test "x$py_val" != xno; then
  AC_MSG_RESULT(yes)
  ifelse([$3], [],, [$3
])dnl
else
  AC_MSG_RESULT(no)
  ifelse([$4], [],, [$4
])dnl
fi
])

dnl a macro to check for ability to create python extensions
dnl  AM_CHECK_PYTHON_HEADERS([ACTION-IF-POSSIBLE], [ACTION-IF-NOT-POSSIBLE])
dnl function also defines PYTHON_INCLUDES
AC_DEFUN([AM_CHECK_PYTHON_HEADERS],
[AC_REQUIRE([AM_PATH_PYTHON])
AC_MSG_CHECKING(for headers required to compile python extensions)
dnl deduce PYTHON_INCLUDES
AC_ARG_WITH(python-prefix, AS_HELP_STRING([--with-python-prefix],
              [Force Python prefix (disable by default).]),
              [ax_python_prefix=$withval],
              [ax_python_prefix="no"])
if test x"$ax_python_prefix" = x"no" ; then
  py_prefix=`$PYTHON -c "import sys; print (sys.prefix)"`
  py_exec_prefix=`$PYTHON -c "import sys; print (sys.exec_prefix)"`
else
  py_prefix=$ax_python_prefix
  py_exec_prefix=$ax_python_prefix
fi
if test -x "$PYTHON-config"; then
PYTHON_INCLUDES=`$PYTHON-config --includes 2>/dev/null`
else
PYTHON_INCLUDES="-I${py_prefix}/include/python${PYTHON_VERSION}"
if test "$py_prefix" != "$py_exec_prefix"; then
  PYTHON_INCLUDES="$PYTHON_INCLUDES -I${py_exec_prefix}/include/python${PYTHON_VERSION}"
fi
fi
AC_SUBST(PYTHON_INCLUDES)
AC_LANG_PUSH(C)
dnl check if the headers exist:
save_CPPFLAGS="$CPPFLAGS"
CPPFLAGS="$CPPFLAGS $PYTHON_INCLUDES"
AC_TRY_CPP([#include <Python.h>],dnl
[AC_MSG_RESULT(found)
$1],dnl
[AC_MSG_RESULT(not found)
$2])
CPPFLAGS="$save_CPPFLAGS"
  LDFLAGS_SVG="$LDFLAGS"
  LDFLAGS="$LDFLAGS -L${py_prefix}/lib"
  if test "$py_prefix" != "$py_exec_prefix"; then
    LDFLAGS="$LDFLAGS -L${py_exec_prefix}/lib"
  fi
  AC_CHECK_LIB([python${PYTHON_VERSION}], [Py_Initialize],
               [PYTHON_LIBS="-L${py_prefix}/lib -lpython${PYTHON_VERSION}"
  AC_SUBST([PYTHON_LIBS])
$1], [ax_have_pythondev=no
$2])
  LDFLAGS="$LDFLAGS_SVG"
AC_LANG_POP(C)
])

AC_DEFUN([AX_PYTHON_DEV],
[dnl Test for Python availability and usage for bindings.
  AC_ARG_WITH(python, AS_HELP_STRING([--with-python], [Compile Python support (enabled by default).]),
              [ax_have_python=$withval],
              [ax_have_python=m4_default([$1], ["no"])])
  if test x"$ax_have_python" = x"yes" ; then
    AM_PATH_PYTHON([2.3.5], [ax_have_python="yes"], [ax_have_python="no"])
    if test x"$ax_have_python" = x"yes" ; then
      AM_CHECK_PYTHON_HEADERS(,[AC_MSG_WARN(could not find Python development files); ax_have_python="no"])
      if test x"$ax_have_python" = x"yes" ; then
        AC_DEFINE([HAVE_PYTHON], [], [if set, we can call Python.h])
      fi
    
      AC_ARG_WITH(numpy-incs, AS_HELP_STRING([--with-numpy-incs], [Include flag for numpy development files ($PYTHON_INCLUDES).]),
              [ax_python_numpy=$withval],
              [ax_python_numpy=""])
      AC_LANG_PUSH(C)
      CPPFLAGS_SVG=$CPPFLAGS
      CPPFLAGS=$CPPFLAGS" $PYTHON_INCLUDES $ax_python_numpy"
      AC_CHECK_HEADERS([numpy/ndarrayobject.h], [ax_have_numpy=yes], [ax_have_numpy=no],
        [#ifdef HAVE_PYTHON
    #include <Python.h>
    #endif
    ])
      CPPFLAGS=$CPPFLAGS_SVG
      AC_LANG_POP(C)
      if test x"$ax_have_numpy" = x"yes" ; then
        if test x"$ax_python_numpy" != x"$PYTHON_INCLUDES" ; then
          PYTHON_INCLUDES="$PYTHON_INCLUDES $ax_python_numpy"
        fi
        AC_DEFINE([HAVE_PYTHON_NUMPY], [], [if set, we can call numpy/ndarrayobject.h])
      fi
    fi
  fi
  AM_CONDITIONAL([HAVE_PYTHON], [test x"$ax_have_python" = x"yes"])
])
