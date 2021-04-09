!> @file
!! Include fortran file for allocation template
!! 
!! @author
!!    Copyright (C) 2012-2015 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
  if (f_nan_pad_size > 0) call togglepadding(int(0,f_long))
  if (.not.  malloc_validate(ierror,size(shape(array)),m)) return
  call pad_array(array,m%put_to_zero,m%shape,padding)
  !also fill the array with the values of the source if the address is identified in the source
