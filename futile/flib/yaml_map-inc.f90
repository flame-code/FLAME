!> @file
!! Include file used in yaml_output.f90.
!! Body of the yaml_map template.
!! yaml: Yet Another Markeup Language (ML for Human)
!! @author
!!    Copyright (C) 2013-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


  character(len=*), intent(in) :: mapname
  character(len=*), optional, intent(in) :: label,advance,fmt
  integer, optional, intent(in) :: unit
  !local variables
  integer :: msg_lgt,strm,unt
  character(len=3) :: adv
  character(len=tot_max_record_length) :: towrite

  unt=0
  if (present(unit)) unt=unit
  call get_stream(unt,strm)

  adv='def' !default value
  if (present(advance)) adv=advance

  msg_lgt=0
  !put the message
  call buffer_string(towrite,len(towrite),trim(mapname),msg_lgt)
  !put the semicolon
  call buffer_string(towrite,len(towrite),': ',msg_lgt)
  !put the optional name
  if (present(label)) then
     call buffer_string(towrite,len(towrite),' &',msg_lgt)
     call buffer_string(towrite,len(towrite),trim(label)//' ',msg_lgt)
  end if
  !put the value
  if (present(fmt)) then
     call buffer_string(towrite,len(towrite),trim(yaml_toa(mapvalue,fmt=fmt)),msg_lgt)
  else
     call buffer_string(towrite,len(towrite),trim(yaml_toa(mapvalue)),msg_lgt)
  end if
  call dump(streams(strm),towrite(1:msg_lgt),advance=trim(adv),event=MAPPING)
