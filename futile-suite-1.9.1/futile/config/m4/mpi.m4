# -*- Autoconf -*-
#
# Copyright (c) 2014 BigDFT Group (Damien Caliste)
# All rights reserved.
#
# This file is part of the BigDFT software package. For license information,
# please see the COPYING file in the top-level directory of the BigDFT source
# distribution.
#

AC_DEFUN([AX_MPI],
[dnl Enable or not MPI.
  AC_LANG_PUSH(Fortran)
  AC_REQUIRE([AC_PROG_FC])
  
  ac_use_mpi="yes"
  AC_ARG_ENABLE(mpi, AS_HELP_STRING([--disable-mpi], [Don't use MPI (enabled by default).]),
                [ac_use_mpi=$enableval; poi_mpi_req="yes"],
                [ac_use_mpi="yes"; poi_mpi_req="no"])
  
  if test "$ac_use_mpi" = "yes"; then
    dnl We set all flags by hand
    AC_ARG_WITH(mpi-include, AS_HELP_STRING([--with-mpi-include],
                [specify INCLUDE for MPI in case of non-built-in library. Use the -I before the path(es).]),
                ac_mpi_include=$withval, ac_mpi_include=$MPI_INCLUDE)
    AC_ARG_WITH(mpi-ldflags, AS_HELP_STRING([--with-mpi-ldflags],
                [specify LDFLAGS for MPI in case of non-built-in library. Use the -L before the path(es).]),
                ac_mpi_ldflags=$withval, ac_mpi_ldflags=$MPI_LDFLAGS)
    AC_ARG_WITH(mpi-libs, AS_HELP_STRING([--with-mpi-libs],
                [specify the libraries for MPI in case of non-built-in library. Use the -l before the name(s).]),
                ac_mpi_libs=$withval, ac_mpi_libs=$MPI_LIBS)
  
    LIBS_OLD=$LIBS
    LIBS="$LIBS $ac_mpi_libs"
    LDFLAGS_OLD=$LDFLAGS
    LDFLAGS="$LDFLAGS $ac_mpi_ldflags"
    FCFLAGS_OLD=$FCFLAGS
    FCFLAGS="$FCFLAGS $ac_mpi_include"
  
    AC_MSG_CHECKING([for MPI capabilities])
    AC_LINK_IFELSE([[program main
    include "mpif.h"
  
    call MPI_INIT
  end]], withmpi=yes, withmpi=no)
    AC_MSG_RESULT([$withmpi])
  
    if test x"$withmpi" = x"yes" ; then
      dnl The compiler will pass everything needed
      MPI_LIBS=$ac_mpi_libs
      MPI_INCLUDE=$ac_mpi_include
      MPI_LDFLAGS=$ac_mpi_ldflags
  
      dnl Test for MPI_INIT_THREAD. Use link instead of run due to execution prob.
      AC_MSG_CHECKING([for MPI_INIT_THREAD availability])
      AC_LINK_IFELSE([[
  program main
  implicit none
  
  integer :: ierr,provided
  
  include "mpif.h"
  
  call MPI_INIT_THREAD(MPI_THREAD_FUNNELED,provided,ierr)
  write(*, "(I1)", advance = "NO") provided
  
  !test the availability of max active levels
  !$   call OMP_SET_NESTED(.false.) 
  !$   call OMP_SET_MAX_ACTIVE_LEVELS(1)
  
  end
  ]], withmpinithread=yes, withmpinithread=no, withmpinithread=yes)
      AC_MSG_RESULT([$withmpinithread])
  
      if test x"$withmpinithread" = x"yes" ; then
        ac_with_mpinithread="yes"
        AC_DEFINE([HAVE_MPI_INIT_THREAD], [1], [use MPI_INIT_THREAD])
      else
        ac_with_mpinithread="no"
      fi
  
      dnl Test for MPI2.
      AC_ARG_WITH(mpi2, AS_HELP_STRING([--with-mpi2],
                  [Use MPI2 capabilities like MPI_IN_PLACE.]),
                  ac_with_mpi2=$withval, ac_with_mpi2="auto")
      if test x"$ac_with_mpi2" == x"yes" -o x"$ac_with_mpi2" == x"auto" ; then
        AC_MSG_CHECKING([for MPI 2 capabilities])
        AC_LINK_IFELSE([[
  program main
  implicit none
  
  integer :: ierr
  real :: rho(3)
  
  include "mpif.h"
  
  call MPI_ALLREDUCE(MPI_IN_PLACE,rho,3,&
                MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_ACCUMULATE()

  end
  ]], ac_with_mpi2="yes", ac_with_mpi2="no")
        AC_MSG_RESULT([$ac_with_mpi2])
      fi
  
      if test x"$ac_with_mpi2" = x"yes" ; then
        AC_DEFINE([HAVE_MPI2], [1], [use MPI2 capabilities.])
      fi
      #AM_CONDITIONAL(HAVE_MPI2, test "$ac_with_mpi2" = "yes")
      dnl Test for MPI3.
      AC_ARG_WITH(mpi3, AS_HELP_STRING([--with-mpi3],
                [Use MPI3 capabilities like MPI_IALLREDUCE.]),
                ac_with_mpi3=$withval, ac_with_mpi3="auto")
      if test x"$ac_with_mpi3" == x"yes" -o x"$ac_with_mpi3" == x"auto" ; then
      AC_MSG_CHECKING([for MPI 3 capabilities])
      AC_LINK_IFELSE([[
program main
implicit none

integer :: ierr,request,nproc,iproc,jproc
real :: rho(3),rho2(3)
real,dimension(:),allocatable :: sendbuf,recvbuf
integer,dimension(:),allocatable :: sendcounts,sdispls,recvcounts,rdispls

include "mpif.h"

call mpi_comm_rank(MPI_COMM_WORLD,iproc,ierr)
call mpi_comm_size(MPI_COMM_WORLD,nproc,ierr)
allocate(sendbuf(nproc),recvbuf(nproc))
allocate(sendcounts(0:nproc-1),sdispls(0:nproc-1),recvcounts(0:nproc-1),rdispls(0:nproc-1))
do jproc=0,nproc-1
	sendcounts(jproc)=1
	sdispls(jproc)=jproc
	recvcounts(jproc)=1
	rdispls(jproc)=jproc
end do

call MPI_IALLREDUCE(rho2,rho,3,&
              MPI_REAL,MPI_SUM,MPI_COMM_WORLD,request,ierr)
call MPI_IALLTOALLV(sendbuf,sendcounts,sdispls,MPI_REAL,&
					recvbuf,recvcounts,rdispls,MPI_REAL,MPI_COMM_WORLD,request,ierr)
end
]], ac_with_mpi3="yes", ac_with_mpi3="no")
      AC_MSG_RESULT([$ac_with_mpi3])
     fi

     #if test x"$ac_with_mpi3" = x"yes" ; then
     # AC_DEFINE([HAVE_MPI3], [1], [use MPI3 capabilities (like MPI_IALLREDUCE and MPI_IALLTOALLV).])
     #fi
  
    else
      AC_MSG_WARN([MPI is not working and has been disabled.])
      AC_MSG_WARN([  If you don't want to see this warning,])
      AC_MSG_WARN([  set --disable-mpi in the configure line.])
      ac_use_mpi="warn"
      MPI_LIBS=
      MPI_INCLUDE=
      MPI_LDFLAGS=
    fi
  
    LIBS="$LIBS_OLD"
    LDFLAGS="$LDFLAGS_OLD"
    FCFLAGS="$FCFLAGS_OLD"
  else
    MPI_LIBS=
    MPI_INCLUDE=
    MPI_LDFLAGS=
  fi

  if test "${ac_use_mpi}" = "warn" -a "${poi_mpi_req}" = "yes"; then
    AC_MSG_FAILURE([MPI was requested by the user but is not available - aborting])
  fi
  
  AM_CONDITIONAL(USE_MPI, test "$ac_use_mpi" = "yes")
  AM_CONDITIONAL(HAVE_MPI3, test "$ac_with_mpi3" = "yes")
  AC_SUBST(MPI_LIBS)
  AC_SUBST(MPI_INCLUDE)
  AC_SUBST(MPI_LDFLAGS)
  
  AC_LANG_POP(Fortran)
]) # AX_MPI
