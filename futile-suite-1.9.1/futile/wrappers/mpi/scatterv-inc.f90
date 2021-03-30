!> @file
!! Include fortran file for scatterv operations
!! @author
!!    Copyright (C) 2015-2015 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS

  mpi_comm=MPI_COMM_WORLD
  if (present(comm)) mpi_comm=comm
  root_=0
  if (present(root)) root_=root 

  call MPI_SCATTERV(sendbuf,sendcounts,displs,mpitype(sendbuf),&
       recvbuf,recvcount,mpitype(recvbuf),root_,mpi_comm,ierr)

  if (ierr/=0) then
     call f_err_throw('Error in MPI_SCATTERV',err_id=ERR_MPI_WRAPPERS)
  end if
