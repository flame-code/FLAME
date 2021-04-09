!> @file
!! Include file used in yaml_output.f90.
!! Body of the yaml_map template for matrices.
!! yaml: Yet Another Markup Language (ML for Human)
!! @author
!!    Copyright (C) 2013-2016 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS

  character(len=*), intent(in) :: mapname
  character(len=*), optional, intent(in) :: label,advance,fmt
  integer, optional, intent(in) :: unit
  !Local variables
  integer :: irow,icol

  !open the sequence associated to the matrix
  call yaml_sequence_open(mapname,label=label,advance=advance,unit=unit)
  do irow=lbound(mapvalue,2),ubound(mapvalue,2)
     call yaml_newline(unit=unit)
     call yaml_sequence(advance='no',unit=unit)
     call yaml_sequence_open(flow=.true.,unit=unit)
     do icol=lbound(mapvalue,1),ubound(mapvalue,1)
        call yaml_sequence(trim(yaml_toa(mapvalue(icol,irow),fmt=fmt)),&
             advance=advance,unit=unit)
     end do
     call yaml_sequence_close(unit=unit)
  end do
  call yaml_sequence_close(advance=advance,unit=unit)
