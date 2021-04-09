!> @file
!! Include fortran file for allreduce operations for arrays
!! @author
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
  include 'allreduce-decl-inc.f90'

  integer(f_long) :: ntotrecv
  ntot=int(f_size(sendbuf),fmpi_integer)
  if (present(recvbuf)) then
     ntotrecv=f_size(recvbuf)
  else
     ntotrecv=int(ntot,fmpi_integer)
  end if
  if (ntot /= ntotrecv) then
     call f_err_throw('Error in allreduce; the size of send ('//ntot//&
          ') and of receive buffer ('//ntotrecv//&
          ') does not coincide',err_id=ERR_MPI_WRAPPERS)
     return
  end if

  in_place=have_mpi2 .and. .not. present(recvbuf)

  if (.not. present(recvbuf) .and. .not. in_place) then
     !case without mpi_in_place
     copybuf = f_malloc(src=sendbuf,id='copybuf')
     !call f_memcpy(src=sendbuf,dest=copybuf,n=ntot)
  end if

  include 'allreduce-core-inc.f90'
