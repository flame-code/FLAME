!> @file
!! Include fortran file for allocation template
!! @author
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


  integer :: ierror
  ! integer :: sizeof
  integer(kind=8) :: iadd
  ! integer(kind=8) :: ilsize
  !$ logical :: not_omp
  !$ logical, external :: omp_in_parallel,omp_get_nested

  if (f_err_raise(ictrl == 0,&
       'ERROR (f_malloc): the routine f_malloc_initialize has not been called',&
       ERR_MALLOC_INTERNAL)) return

  !$ not_omp=.not. (omp_in_parallel() .or. omp_get_nested())

  !here we should add a control of the OMP behaviour of allocation
  !in particular for what concerns the OMP nesting procedure
  !the following action is the allocation
  !$ if(not_omp) then
  call f_timer_interrupt(TCAT_ARRAY_ALLOCATIONS)
  !$ end if
