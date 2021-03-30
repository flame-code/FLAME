!> @file
!! Include fortran file for f_malloc routines
!! initialize the internal info
!! @author
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
  type(array_bounds), dimension(:), intent(in) :: bounds
  !local variables
  integer :: i
  include 'f_malloc-base-inc.f90'
  m%rank=size(bounds)
  do i=1,m%rank
     m%lbounds(i)=bounds(i)%nlow
     m%ubounds(i)=bounds(i)%nhigh
     m%shape(i)=m%ubounds(i)-m%lbounds(i)+1
  end do
