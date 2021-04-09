!> @file
!! Include fortran file for allreduce operations

!! @author
!!    Copyright (C) 2017-2017 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
  type(f_enumerator), intent(in), optional :: op
  integer, intent(in), optional :: comm
  integer(fmpi_integer), intent(out), optional :: request

  if (.not. present(op)) call f_err_throw('MPI_OP should be present in the multiple fmpi_allred',&
       err_id=ERR_MPI_WRAPPERS)

  tmpsend=0
  tmpsend(1)=val1
  tmpsend(2)=val2
  if (present(val3)) tmpsend(3)=val3
  if (present(val4)) tmpsend(4)=val4
  if (present(val5)) tmpsend(5)=val5
  if (present(val6)) tmpsend(6)=val6
  if (present(val7)) tmpsend(7)=val7
  if (present(val8)) tmpsend(8)=val8
  if (present(val9)) tmpsend(9)=val9
  if (present(val10)) tmpsend(10)=val10

  call fmpi_allreduce(sendbuf=tmpsend,recvbuf=tmprecv,op=op,comm=comm,request=request)

  val1=tmprecv(1)
  val2=tmprecv(2)
  if (present(val3)) val3=tmprecv(3)
  if (present(val4)) val4=tmprecv(4)
  if (present(val5)) val5=tmprecv(5)
  if (present(val6)) val6=tmprecv(6)
  if (present(val7)) val7=tmprecv(7)
  if (present(val8)) val8=tmprecv(8)
  if (present(val9)) val9=tmprecv(9)
  if (present(val10)) val10=tmprecv(10)
