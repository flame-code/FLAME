!> @file
!! Include fortran file for gather operations
!! @author
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
  integer, intent(in), optional :: root !< 0 if absent
  integer, intent(in), optional :: comm !< MPI_COMM_WORLD if absent
  !local variables
  integer :: iroot,mpi_comm,ntot,ntotrecv,ntasks,ierr

  ntot=size(sendbuf)
  ntotrecv=size(recvbuf)

  include 'gather-inner-inc.f90'
