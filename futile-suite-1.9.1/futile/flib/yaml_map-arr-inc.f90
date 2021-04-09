!> @file
!! Include file used in yaml_output.f90.
!! Body of the yaml_map template for arrays.
!! yaml: Yet Another Markup Language (ML for Human)
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
  integer :: msg_lgt,strm,nl,nu,tmp_lgt,i,il,icursor
  character(len=tot_max_record_length) :: towrite

  strm=stream_id(unit=unit)

  nl=lbound(mapvalue,1)
  nu=ubound(mapvalue,1)

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

  !check whether the final message will be too long or not
  if (nu-nl > 0) then
     tmp_lgt=0
     do il=nl,nu
        !change strategy if the remaining space is too low
        !template of an element
        tmp_lgt=tmp_lgt+len_trim(yaml_toa(mapvalue(il),fmt=fmt))
        tmp_lgt=tmp_lgt+3 !comma and spaces
     end do
     !     tmp_lgt=tmp_lgt*(nu-nl)
     !print *,'debug',max(streams(strm)%icursor+msg_lgt+1,streams(strm)%tabref)+tmp_lgt,&
     !     streams(strm)%max_record_length
     if (max(streams(strm)%icursor+msg_lgt+1,streams(strm)%tabref)+tmp_lgt > &
          streams(strm)%max_record_length) then
        !implement the writing explicitly per element
        call yaml_sequence_open(mapname,flow=.true.,unit=unit)
        do i=nl,nu
         call yaml_stream_attributes(icursor=icursor)
         !print *,'i,icursor',i,icursor
         call yaml_sequence(trim(yaml_toa(mapvalue(i),fmt=fmt)),unit=unit)
        end do
        call yaml_sequence_close(unit=unit)
        return
     end if
  end if
!print *,'debug'
  !put the value
  call buffer_string(towrite,len(towrite),trim(yaml_toa(mapvalue,fmt=fmt)),msg_lgt)
!print *,'debug2',towrite(1:msg_lgt),'msggt'
  call dump(streams(strm),towrite(1:msg_lgt),advance=advance,event=MAPPING)
