!> @file
!! Include fortran file for maxdiff interfaces
!! @author
!!    Copyright (C) 2012-2014 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
  if (ns ==-1 .and. nd >=0) then
     cnt=nd/kind(b)
  else if (nd == -1 .and. ns >= 0) then
     cnt=ns/kind(a)
  else if (nd==-1 .and. ns==-1) then
     cnt=1
     if (present(n)) cnt=n
  else
     cnt=min(ns/kind(a),nd/kind(b))
  end if
  if (present(n)) then
     if (n > cnt) then
        call f_err_throw('Error in f_maxdiff; sizes of the source ('//&
             trim(yaml_toa(ns))//&
             ') and of the destination buffer ('//trim(yaml_toa(nd))//&
             ') are not compatible with the given value of n ('//&
             trim(yaml_toa(n))//')',err_id=ERR_INVALID_COPY)
        return
     end if
     cnt=n
  end if
  call f_diff(cnt,a,b,maxdiff)
