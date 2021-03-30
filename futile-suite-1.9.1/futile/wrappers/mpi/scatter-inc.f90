!> @file
!! Include fortran file for scatter operations
!! @author
!!    Copyright (C) 2015-2015 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
  integer, intent(in), optional :: root !< 0 if absent
  integer, intent(in), optional :: comm !< MPI_COMM_WORLD if absent
  !local variables
  integer :: iroot,mpi_comm,ntot,ntotsend,ntasks,ierr

  ntot=size(recvbuf)
  ntotsend=size(sendbuf)

      

  if (present(root)) then
     iroot=root
  else
     iroot=0
  end if
  if (present(comm)) then
     mpi_comm=comm
  else
     mpi_comm=MPI_COMM_WORLD !or bigdft_mpi%mpi_comm?
  end if

  !verify the size of the send buffer
  ntasks=mpisize(mpi_comm)
  if (ntotsend*kind(sendbuf) < ntot*ntasks*kind(recvbuf)) then
     call f_err_throw('Error in mpiscatter; the size of the send buffer ('//&
          (ntotsend*kind(sendbuf))//&
          ') is not large enough to contain '//(ntot*kind(recvbuf))//&
          ' * '//ntasks//' elements',err_id=ERR_MPI_WRAPPERS)
     return
  end if
  !then one can proceed with the MPI operation
  ntotsend=int(int(ntot,kind=8)*kind(sendbuf)/int(kind(recvbuf),kind=8))
  call f_timer_interrupt(TCAT_SCATTER)
  call MPI_SCATTER(sendbuf,ntot,mpitype(sendbuf),&
       recvbuf,ntotsend,mpitype(recvbuf),iroot,mpi_comm,ierr)
  call f_timer_resume()
  if (ierr /=0) then
     call f_err_throw('An error in calling to MPI_SCATTER occured',&
          err_id=ERR_MPI_WRAPPERS)
     return
  end if
