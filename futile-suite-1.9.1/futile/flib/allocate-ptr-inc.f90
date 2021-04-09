!> @file
!! Include fortran file for allocation template
!! 
!! @author
!!    Copyright (C) 2012-2015 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
  if (m%srcdata_add == int(-1,f_address)) then
     nullify(array)
     if (f_nan_pad_size > 0) call togglepadding(0)
     call f_timer_resume()!TCAT_ARRAY_ALLOCATIONS
     !not possible: the pointer could have not been nullified call f_free_ptr(array) !to avoid memory leaks
     return
  end if

  
