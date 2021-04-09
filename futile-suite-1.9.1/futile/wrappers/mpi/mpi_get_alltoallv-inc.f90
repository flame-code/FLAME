!> @file
!! Include fortran file for mpi_get_alltoallv operations
!! @author
!!    Copyright (C) 2017-2017 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


  integer,intent(in) :: iproc, nproc, comm
  integer,dimension(0:nproc-1),intent(in) :: nsendcounts, nsenddspls, nrecvcounts, nrecvdspls
  
  integer :: ierr, info, window, jproc
  integer,dimension(:),allocatable :: nsenddspls_remote

  if (size(sendbuf)<sum(nsendcounts)) then
      call f_err_throw("insufficient size of array 'sendbuf': "//trim(yaml_toa(size(sendbuf)))//&
           &" instead of "//trim(yaml_toa(sum(nsendcounts))))
  end if
  if (size(recvbuf)<sum(nrecvcounts)) then
      call f_err_throw("insufficient size of array 'recvbuf': "//trim(yaml_toa(size(recvbuf)))//&
           &" instead of "//trim(yaml_toa(sum(nrecvcounts))))
  end if
  
  nsenddspls_remote = f_malloc(0.to.nproc-1,id='nsenddspls_remote')
  
  call mpi_alltoall(nsenddspls, 1, mpi_integer, &
                    nsenddspls_remote, 1, mpi_integer, &
                    comm, ierr)
  
  !info=mpiinfo("no_locks", "true")
  window = mpiwindow(size(sendbuf), sendbuf(1), comm)
  do jproc=0,nproc-1
      if (nrecvcounts(jproc)>0) then
          call mpi_get(recvbuf(nrecvdspls(jproc)+1), nrecvcounts(jproc), mpitype(recvbuf), &
                       jproc, int(nsenddspls_remote(jproc),kind=mpi_address_kind), &
                       nrecvcounts(jproc), mpitype(sendbuf), window, ierr)
      end if
  end do
  call mpi_fenceandfree(window)
  !call mpiinfofree(info)
  
  call f_free(nsenddspls_remote)
