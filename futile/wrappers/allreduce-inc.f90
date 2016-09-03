!> @file
!! Include fortran file for allreduce operations

!! @author
!!    Copyright (C) 2012-2015 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
  integer, intent(in) :: op,count
  integer, intent(in), optional :: comm
  !local variables
  logical :: in_place
  integer :: tcat,ierr,ntot,mpi_comm

  external :: MPI_ALLREDUCE

  ntot=count

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
     copybuf = f_malloc(ntot,id='copybuf')
     call f_memcpy(src=sendbuf,dest=copybuf(1),n=ntot)
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
