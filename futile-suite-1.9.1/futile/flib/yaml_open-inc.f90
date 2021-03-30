  !> @file
  !! Include file used in yaml_output.f90.
  !! Body of the yaml_mapping_open or sequence template for arrays.
  !! yaml: Yet Another Markup Language (ML for Human)
  !! @author
  !!    Copyright (C) 2013-2013 BigDFT group
  !!    This file is distributed under the terms of the
  !!    GNU General Public License, see ~/COPYING file
  !!    or http://www.gnu.org/copyleft/gpl.txt .
  !!    For the list of contributors, see ~/AUTHORS
  
  !template:
  !subroutine yaml_open_<template>(mapname,label,tag,flow,tabbing,advance,unit)
  !implicit none
  
  
  ! SM: doxygen does seem to get the descrption with @copydoc if this file is used with include, therefore the arguments
  !     are now again back in the main routines.
  !!  character(len=*), optional, intent(in) :: mapname !< Key of the sequence. @copydoc doc::mapname
  !!  character(len=*), optional, intent(in) :: label   !< @copydoc doc::label
  !!  character(len=*), optional, intent(in) :: tag     !< @copydoc doc::tag
  !!  logical, optional, intent(in) :: flow             !< @copydoc doc::flow
  !!  character(len=*), optional, intent(in) :: advance !< @copydoc doc::advance
  !!  integer, optional, intent(in) :: unit             !< @copydoc doc::unit
  !!  integer, optional, intent(in) :: tabbing          !< @copydoc doc::tabbing
  !local variables
  logical :: doflow
  integer :: msg_lgt,tb,ipos
  integer :: unt,strm
  character(len=3) :: adv
  character(len=tot_max_record_length) :: towrite

  unt=0
  if (present(unit)) unt=unit
  call get_stream(unt,strm)

  doflow=streams(strm)%flowrite
  !override if already active
  if (present(flow)) doflow=flow .or. doflow

  !Position of the cursor
  ipos=max(streams(strm)%icursor,streams(strm)%indent)

  msg_lgt=0
  !put the message
  if (present(mapname)) then
     if (len_trim(mapname)>0) then
        call buffer_string(towrite,len(towrite),trim(mapname),msg_lgt)
        !add some spaces if required
        if (present(tabbing)) then
           ipos=ipos+msg_lgt
           tb=max(tabbing-ipos-1,1)
           call buffer_string(towrite,len(towrite),repeat(' ',tb),msg_lgt)
           ipos=ipos+tb
        end if
        !put the semicolon
        call buffer_string(towrite,len(towrite),':',msg_lgt)
     end if
  end if
  !put the optional tag
  if (present(tag)) then
     if (len_trim(tag)>0) then
        call buffer_string(towrite,len(towrite),' !',msg_lgt)
        call buffer_string(towrite,len(towrite),trim(tag),msg_lgt)
     end if
  end if
  !put the optional name
  if (present(label)) then
     if(len_trim(label)>0) then
        call buffer_string(towrite,len(towrite),' &',msg_lgt)
        call buffer_string(towrite,len(towrite),trim(label),msg_lgt)
     end if
  end if

  call open_level(streams(strm),doflow)

  if (doflow .or. msg_lgt==0) then
     adv='no '
  else
     adv='yes'
     if (present(advance)) adv = advance
  end if

