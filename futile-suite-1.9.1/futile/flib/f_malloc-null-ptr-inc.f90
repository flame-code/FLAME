  !> @file
  !! Other fortran file for f_malloc routines
  !! @author
  !!    Copyright (C) 2012-2015 BigDFT group
  !!    This file is distributed under the terms of the
  !!    GNU General Public License, see ~/COPYING file
  !!    or http://www.gnu.org/copyleft/gpl.txt .
  !!    For the list of contributors, see ~/AUTHORS
  
  !if the src_ptr is nullified the f_malloc should provide nullification
  if (.not. associated(src_ptr)) then
     m%srcdata_add=int(-1,f_address)
     return
  end if
