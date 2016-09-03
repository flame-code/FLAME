
  integer,dimension(0:),intent(in) :: sendcounts, sdispls
  integer,dimension(0:),intent(in) :: recvcounts, rdispls
  integer,intent(in) :: comm
  ! Local variables
  integer :: ierr
  
  call mpi_alltoallv(sendbuf, sendcounts, sdispls, mpitype(sendbuf), &
       recvbuf, recvcounts, rdispls, mpitype(recvbuf), comm, ierr)
  
  if (ierr/=0) then
     call f_err_throw('An error in calling to MPI_ALLTOALLV occured',&
          err_id=ERR_MPI_WRAPPERS)
     return
  end if
