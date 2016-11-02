!> @file
!! Include fortran file for f_malloc routines
!! initialize the internal info
!! @author
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
  !local variables
  integer :: i 

  m%pin=.false.
  m%profile=f_malloc_default_profiling !mems(ictrl)%profile_routine !< here omp can be used to know whether to profile or not
  m%put_to_zero=.false.
  m%srcdata_add=int(0,f_address)
  m%rank=1
  m%shape=0
  m%lbounds=1
  m%ubounds=0
  do i=1,f_malloc_namelen
     m%array_id(i:i)=' '
     m%routine_id(i:i)=' '
  end do
  m%info=' '
