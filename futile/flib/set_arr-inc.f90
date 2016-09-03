!> @file
!!  template to be used to add different items to a dictionary
!! @author
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS

!! included in dictionaries.f90
!header of routine
!!$subroutine set(dict,list)
!!$  implicit none
!!$  <generic>, dimension(:), intent(in) :: list
  type(dictionary), pointer :: dict
  !local variables
  integer :: item,nitems,nitems_old

  nitems=size(list)
  nitems_old=dict_len(dict)
  do item=1,nitems
     call set(dict//(item-1),list(item))
  end do
  do item=nitems_old-1,nitems,-1
     call dict_remove(dict,item)
  end do
!!$end subroutine set
