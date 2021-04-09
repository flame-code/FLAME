!> @file
!! Include file used in yaml_strings.f90
!! Body of the yaml_toa template. To be used for arrays.
!! yaml: Yet Another Markup Language (ML for Human)
!! @author
!!    Copyright (C) 2013-2016 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
  character(len=max_value_length) :: vec_toa
  character(len=*), optional, intent(in) :: fmt
  !local variables
  character(len=max_value_length) :: tmp
  integer :: nl,nu,i,length,pos

  tmp=repeat(' ',max_value_length)
  vec_toa=tmp

  nl=lbound(vec,1)
  nu=ubound(vec,1)

  if (nl > nu) then
     !Special case for size 0 (nl is > nu!)
     vec_toa(1:2) = '[]'
  else
     vec_toa(1:2)='[ '
     pos=3
     do i=nl,nu
        if (present(fmt)) then
           tmp=yaml_toa(vec(i),fmt=fmt)
        else
           tmp=yaml_toa(vec(i))
        end if
        length=len(trim(tmp))-1

        !Truncate the chain if too long
        if (pos+length+2 > max_value_length) exit
        
        vec_toa(pos:pos+length)=tmp(1:length+1)
        if (i < nu) then
           vec_toa(pos+length+1:pos+length+2)=', '
        else
           vec_toa(pos+length+1:pos+length+2)=' ]'
        end if
        pos=pos+length+3
     end do
  end if

  vec_toa=yaml_adjust(vec_toa)
