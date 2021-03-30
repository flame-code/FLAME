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
  character(len=info_length) :: val
  type(dictionary), pointer :: info
  logical :: c_allocation
  integer(f_kind) :: alignment
  integer(f_kind) :: aligned_size
  integer :: padding
  ! integer(kind=8) :: ilsize

  if (f_err_raise(ictrl == 0,&
       'ERROR (f_malloc): the routine f_malloc_initialize has not been called',&
       ERR_MALLOC_INTERNAL)) return

  !no padding by default
  padding=ndebug !which is always zero so far

  !here we should add a control of the OMP behaviour of allocation
  !in particular for what concerns the OMP nesting procedure
  !the following action is the allocation
  call f_timer_interrupt(TCAT_ARRAY_ALLOCATIONS)
  if (f_nan_pad_size > 0) then
     if (product(int(m%shape(1:m%rank),f_long)) > int(0,f_long)) then
        padding=f_nan_pad_size
        call togglepadding(product(int(m%shape(1:m%rank-1),f_long))*&
             kind(array)*(m%shape(m%rank)+padding))
     end if
  end if
