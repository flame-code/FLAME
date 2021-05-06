!> @file
!! Include fortran file used by the error handling
!! @author
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Piece of routine to identify the errors, consider two optional arguments err_id and err_name
!! included in error_handling.f90
  !local variables
  logical :: isthere
  integer :: get_error 
  integer :: nerr,ierr,jerr
!  character(len=max_field_length) :: name

  get_error=-1 !no error specified
  nerr=dict_len(dict_present_error)

  if (present(err_name)) then
     get_error=0
     do ierr=0,nerr-1
        !this one can be substituted by the values of the dictionary
        jerr=dict_present_error//ierr//ERRID
        isthere=err_name .in. dict_errors//jerr
        if (isthere) then
           get_error=1 !name
           exit
        end if
     end do

  else if (present(err_id)) then
     get_error=0
     do ierr=0,nerr-1
        jerr=dict_present_error//ierr//ERRID
        if (jerr==err_id) then
           get_error=2
           exit
        end if
     end do
  end if
