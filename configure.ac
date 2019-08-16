AC_PREREQ(2.59)
dnl Init basic data for package
dnl define following variables :
dnl  - PACKAGE_NAME, the name ;
dnl  - PACKAGE_TARNAME, the name used for the tarball ;
dnl  - PACKAGE_VERSION, the version ;
dnl  - PACKAGE_STRING, the name + the version ;
dnl  - PACKAGE_BUGREPORT, the mail.
AC_INIT([FLAME - Fully Loaded library for Atomistic Modelling Environment],[1.0],[alirezagh76@gmail.com],[flame])
AC_PROG_RANLIB()
dnl for configuring the FCFUNC macro
AC_CONFIG_HEADERS(config.h)
dnl Init something, don't know exactly what...
AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE([--warnings=no-portability, subdir-objects])
prefix_pwd=$(pwd)
AC_PREFIX_DEFAULT([$prefix_pwd])

dnl If FC environment variable is not already set, it
dnl looks for a modern Fortran compiler, prefering newest
dnl fortran available (e. g. Fortran95). Then FC variable
dnl is set.
AC_PROG_FC([], [Fortran 90])
AC_PROG_F77([], [Fortran 77])
AC_FC_SRCEXT(f90)
AC_LANG(Fortran)
AC_FC_LINE_LENGTH([unlimited])
AC_FC_WRAPPERS()

AC_PROG_CC()
AC_PROG_CXX()

dnl Check for existence of other options
AX_LAMMPS()
AX_SPGLIB()
AX_BPS()
AX_MKL()

dnl pass the FUTILE variable to the Makefiles
AC_ARG_VAR([FUTILE], [FUTILE root path])

dnl Check for sed and put it in the SED variable
AC_PROG_SED()

dnl Give the name of file.in to transform to file
AC_CONFIG_FILES([
Makefile
src/Makefile
src/lenosky_tb/Makefile
tests/Makefile
])

AC_OUTPUT
