!> @file
!! Other fortran file for f_malloc routines for f_buffers, declaration
!! @author
!!    Copyright (C) 2012-2015 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
      type(malloc_information_buf) :: m
      integer, intent(in), optional :: sizes
      !local variables
      integer :: i

