!> @file
!! Include fortran file for f_malloc routines
!! initialize the internal info
!! @author
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
  integer, dimension(:), intent(in), optional :: sizes,lbounds,ubounds
  !local variables
  integer :: i
  include 'f_malloc-base-inc.f90'
  include 'f_malloc-extra-inc.f90'
