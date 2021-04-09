!> @file
!! Include fortran file for allreduce operations

!! @author
!!    Copyright (C) 2012-2015 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
  include 'allreduce-decl-inc.f90'

  integer, intent(in) :: count
  ntot=count

  in_place=have_mpi2 .and. .not. present(recvbuf)

  if (.not. present(recvbuf) .and. .not. in_place) then
     !case without mpi_in_place
     copybuf = f_malloc(ntot,id='copybuf')
     call f_memcpy(src=sendbuf,dest=copybuf,n=ntot)
  end if

  include 'allreduce-core-inc.f90'
