!> @file
!! Include fortran file for allreduce operations for arrays
!! @author
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
  integer, intent(in) :: op
  integer, intent(in), optional :: comm
  !local variables
  logical :: in_place
  integer :: tcat,ierr,ntot,ntotrecv,mpi_comm

  ntot=size(sendbuf)
  if (present(recvbuf)) then
     ntotrecv=size(recvbuf)
  else
     ntotrecv=ntot
  end if
  if (ntot /= ntotrecv) then
     call f_err_throw('Error in allreduce; the size of send ('//ntot//&
          ') and of receive buffer ('//ntotrecv//&
          ') does not coincide',err_id=ERR_MPI_WRAPPERS)
     return
  end if

  if (ntot <= 0) return
  in_place=have_mpi2 .and. .not. present(recvbuf)

  if (present(comm)) then
     mpi_comm=comm
  else
     mpi_comm=MPI_COMM_WORLD !or bigdft_mpi%mpi_comm?
  end if

  tcat=TCAT_ALLRED_LARGE
  if (ntot<=smallsize) tcat=TCAT_ALLRED_SMALL
  if (in_place) then
     call f_timer_interrupt(tcat)
     call MPI_ALLREDUCE(MPI_IN_PLACE,sendbuf,ntot,&
          mpitype(sendbuf),op,mpi_comm,ierr)
     call f_timer_resume()
  else if (present(recvbuf)) then
     call f_timer_interrupt(tcat)
     call MPI_ALLREDUCE(sendbuf,recvbuf,ntot,&
          mpitype(sendbuf),op,mpi_comm,ierr)
     call f_timer_resume()
  else
     !case without mpi_in_place
     copybuf = f_malloc(src=sendbuf,id='copybuf')
     call f_timer_interrupt(tcat)
     call MPI_ALLREDUCE(copybuf,sendbuf,ntot,&
          mpitype(sendbuf),op,mpi_comm,ierr)
     call f_timer_resume()
     call f_free(copybuf)
  end if

  if (ierr /=0) then
     call f_err_throw('An error in calling to MPI_ALLREDUCE occured',&
          err_id=ERR_MPI_WRAPPERS)
     return
  end if
