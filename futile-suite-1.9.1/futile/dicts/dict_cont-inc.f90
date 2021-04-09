!> @file
!!  template to be used to create a dict container from a value
!! @author
!!    Copyright (C) 2012-2015 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS

!! included in dictionaries.f90
  character(len = *), intent(in) :: key
  type(dictionary_container) :: cont
  
  cont%key(1:max_field_length) = key
  cont%value(1:max_field_length) = yaml_toa(val)
