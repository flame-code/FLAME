!> @file
!! Other fortran file for f_malloc routines for f_buffers
!! @author
!!    Copyright (C) 2012-2015 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
  if (present(src) .eqv. present(sizes)) then
     call f_err_throw('at least sizes or src should be defined'//&
          ',array "'//trim(m%array_id)//'", routine "'//trim(m%routine_id)//'"',&
          ERR_INVALID_MALLOC)
     return
  end if
  m%rank=1
  if (present(src)) then
     m%shape(1:m%rank)=int(size(src),f_kind) !still the rank of a buffer is always one
     !then check the presence and validity of the ubounds otherwise fill it in a meaningful way
     m%srcdata_add=f_loc(src)
  else if (present(sizes)) then
     m%shape(1:m%rank)=int(sizes,f_kind)
  end if
  m%ubounds(1:m%rank)=int(m%lbounds(1:m%rank),f_kind)+m%shape(1:m%rank)-1
