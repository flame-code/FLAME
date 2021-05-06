!> @file
!! Include fortran file for memcpy interfaces
!! @author
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
  if (nd < ns) then
     call f_err_throw('Error in f_memcpy; the size of the source ('//trim(yaml_toa(ns))//&
          ') and of the destination buffer ('//trim(yaml_toa(nd))//&
          ') are not compatible',err_id=ERR_INVALID_COPY)
     return
  end if
  if (ns <=0) return
  if (f_loc(src) /= f_loc(dest)) call c_memcopy(dest,f_loc(src),ns)
