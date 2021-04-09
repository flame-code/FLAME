!> @file
!!  template to be used to create a dict container from array
!! @author
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS

!! included in dictionaries.f90
  character(len = *), intent(in) :: key
  type(dictionary_container) :: cont
  cont=dict_cont_new_with_dict(key,list_new(.item. val))
