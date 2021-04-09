!> @file
!! Include fortran file for f_malloc routines
!! initialize the internal info
!! @author
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
  integer(f_integer), intent(in) :: size
  include 'f_malloc-base-inc.f90'
  m%rank=1
  m%shape(1)=int(size,f_kind)
  m%ubounds(1)=m%shape(1)
