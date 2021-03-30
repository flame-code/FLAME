!> @file
!! Include file used in yaml_strings.f90
!! yaml plus values
!! @author
!!    Copyright (C) 2013-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
 ! character(len=*), intent(in) :: a
!  character(len=len(trim(adjustl(a)))+len(trim(yaml_toa(num)))) :: c
  c=a+yaml_toa(num)
