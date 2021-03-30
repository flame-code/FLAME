!> @file
!! Include fortran file for allreduce operations

!! @author
!!    Copyright (C) 2017-2017 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
  if (ntot <= 0) return

  mpi_comm=fmpi_comm(comm)

  tcat=TCAT_ALLRED_LARGE
  if (ntot<=smallsize) tcat=TCAT_ALLRED_SMALL
  if (present(request)) request=FMPI_REQUEST_NULL
  if (in_place) then
     if (present(request)) then
        call f_timer_interrupt(tcat)
        ierr=FMPI_SUCCESS
        call MPI_IALLREDUCE(FMPI_IN_PLACE,sendbuf,ntot,&
             mpitype(sendbuf),int(toi(op),fmpi_integer),mpi_comm,request,ierr)
        call f_timer_resume()
     else
        call f_timer_interrupt(tcat)
        call MPI_ALLREDUCE(FMPI_IN_PLACE,sendbuf,ntot,&
             mpitype(sendbuf),int(toi(op),fmpi_integer),mpi_comm,ierr)
        call f_timer_resume()
     end if
  else if (present(recvbuf)) then
     if (present(request)) then
        call f_timer_interrupt(tcat)
        call MPI_IALLREDUCE(sendbuf,recvbuf,ntot,&
             mpitype(sendbuf),int(toi(op),fmpi_integer),mpi_comm,request,ierr)
        call f_timer_resume()
     else
        call f_timer_interrupt(tcat)
        call MPI_ALLREDUCE(sendbuf,recvbuf,ntot,&
             mpitype(sendbuf),int(toi(op),fmpi_integer),mpi_comm,ierr)
        call f_timer_resume()
     end if
  else
     if (present(request)) then
        call f_timer_interrupt(tcat)
        call MPI_IALLREDUCE(copybuf,sendbuf,ntot,&
             mpitype(sendbuf),int(toi(op),fmpi_integer),mpi_comm,request,ierr)
        call f_timer_resume()
     else
        call f_timer_interrupt(tcat)
        call MPI_ALLREDUCE(copybuf,sendbuf,ntot,&
             mpitype(sendbuf),int(toi(op),fmpi_integer),mpi_comm,ierr)
        call f_timer_resume()
     end if
     call f_free(copybuf)
  end if

  if (ierr /= FMPI_SUCCESS) then
     call f_err_throw('An error in calling to MPI_ALLREDUCE occured',&
          err_id=ERR_MPI_WRAPPERS)
     return
  end if
