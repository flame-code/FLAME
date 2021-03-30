!> @file
!! Flush routine 
!! @author
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS

!> Routine used to flush a unit file
subroutine f_utils_flush(unit)
  implicit none
  integer, intent(in) :: unit !< Unit file

  ! Do nothing because the compiler does not support flush.
end subroutine f_utils_flush
