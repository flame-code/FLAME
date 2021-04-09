!> @file
!! Include fortran file for gather operations, inner file
!! @author
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
  if (present(root)) then
     iroot=root
  else
     iroot=0
  end if
  mpi_comm=fmpi_comm(comm)
  !verify the size of the receive buffer
  ntasks=mpisize(mpi_comm)
  if (ntotrecv*kind(recvbuf) < ntot*ntasks*kind(sendbuf)) then
     call f_err_throw('Error in mpigather; the size of the receive buffer ('//&
          (ntotrecv*kind(recvbuf))//&
          ') is not large enough to contain '//(ntot*kind(sendbuf))//&
          ' * '//ntasks//' elements',err_id=ERR_MPI_WRAPPERS)
     return
  end if
  !then one can proceed with the MPI operation
  ntotrecv=int(int(ntot,f_long)*kind(sendbuf)/int(kind(recvbuf),f_long))
  call f_timer_interrupt(TCAT_GATHER)
  call MPI_GATHER(sendbuf,ntot,mpitype(sendbuf),&
       recvbuf,ntotrecv,mpitype(recvbuf),iroot,mpi_comm,ierr)
  call f_timer_resume()
  if (ierr /=0) then
     call f_err_throw('An error in calling to MPI_GATHER occured',&
          err_id=ERR_MPI_WRAPPERS)
     return
  end if
