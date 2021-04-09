  !> @file
  !! Include fortran file for f_malloc routines
  !! @author
  !!    Copyright (C) 2012-2015 BigDFT group
  !!    This file is distributed under the terms of the
  !!    GNU General Public License, see ~/COPYING file
  !!    or http://www.gnu.org/copyleft/gpl.txt .
  !!    For the list of contributors, see ~/AUTHORS
  m%rank=size(shape(src))
  m%shape(1:m%rank)=shape(src)
  if (present(lbounds)) then
     if (size(lbounds) /= m%rank) then
        call f_err_throw('The "lbounds" array has not conformal shape',ERR_INVALID_MALLOC)
     end if
     m%lbounds(1:m%rank)=lbounds(1:m%rank)
  end if
  !then check the presence and validity of the ubounds otherwise fill it in a meaningful way
  m%ubounds(1:m%rank)=m%lbounds(1:m%rank)+m%shape(1:m%rank)-1
  if (present(ubounds)) then
     if (size(ubounds) /= m%rank) then
        call f_err_throw('The "ubounds" array has not conformal shape',err_id=ERR_INVALID_MALLOC)
     end if
     if (any(m%ubounds(1:m%rank)/=ubounds(1:m%rank))) then
        call f_err_throw('The "ubounds" and/or "lbounds arrays are not conformal'//&
             ' with the shape of "src"',ERR_INVALID_MALLOC)
     end if
  end if
  m%srcdata_add=f_loc(src)
