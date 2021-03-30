!> @file
!! Include fortran file for metadata adresses of characters and other variables
!! file included in module metadata_interfaces of getadd.f90
!! @author
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
  integer(f_address) :: la
  !local variables
  integer(f_address), external :: f_loc

  if (size(array)==0) then
     la=int(0,f_address)
     return
  end if
