!> @file
!! Other fortran file for f_malloc routines
!! @author
!!    Copyright (C) 2012-2015 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS

  !guess the rank
  m%rank=0

  if (present(lbounds)) then
     m%rank=size(lbounds)
     m%lbounds(1:m%rank)=int(lbounds,f_kind)
  end if

  if (present(sizes)) then
     if (m%rank == 0) then
        m%rank=size(sizes)
     else if (m%rank/=size(sizes)) then
        call f_err_throw('sizes not conformal with lbounds'//&
             ',array "'//trim(m%array_id)//'", routine "'//trim(m%routine_id)//'"',ERR_INVALID_MALLOC)
        return
     end if
     m%shape(1:m%rank)=int(sizes,f_kind)
     do i=1,m%rank
        m%ubounds(i)=int(m%lbounds(i),f_kind)+m%shape(i)-1
     end do
     if (present(ubounds)) then
        if (m%rank/=size(ubounds)) then
           call f_err_throw('sizes not conformal with ubounds'//&
             ',array "'//trim(m%array_id)//'", routine "'//trim(m%routine_id)//'"',ERR_INVALID_MALLOC)
           return
        end if
        do i=1,m%rank
           if (m%ubounds(i) /=int(ubounds(i),f_kind)) then
              call f_err_throw('ubounds not conformal with sizes and lbounds'//&
              ',array "'//trim(m%array_id)//'", routine "'//trim(m%routine_id)//'"',ERR_INVALID_MALLOC)
              return
           end if
        end do
     end if
  else
     if (present(ubounds)) then
        if (m%rank == 0) then
           m%rank=size(ubounds)
        else if (m%rank/=size(ubounds)) then
           call f_err_throw('ubounds not conformal with lbounds'//&
               ',array "'//trim(m%array_id)//'", routine "'//trim(m%routine_id)//'"',ERR_INVALID_MALLOC)
           return
        end if
        m%ubounds(1:m%rank)=int(ubounds,f_kind)
        do i=1,m%rank
           m%shape(i)=m%ubounds(i)-m%lbounds(i)+1
        end do
     else
        call f_err_throw('at least sizes or ubounds should be defined'//&
            ',array "'//trim(m%array_id)//'", routine "'//trim(m%routine_id)//'"',ERR_INVALID_MALLOC)
        return
     end if
  end if
