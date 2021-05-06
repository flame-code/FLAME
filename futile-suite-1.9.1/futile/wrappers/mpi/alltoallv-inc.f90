!> @file
!! Include fortran file for mpi_alltoall

!! @author
!!    Copyright (C) 2017-2017 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
  integer, intent(in), optional :: count
  integer, dimension(0:), intent(in), optional :: sendcounts,sdispls
  integer, dimension(0:), intent(in), optional :: recvcounts,rdispls
  integer(fmpi_integer), intent(out), optional :: request
  type(fmpi_win), intent(out), optional :: win
  integer,intent(in), optional :: comm
  type(f_enumerator), intent(in), optional :: algorithm
  ! Local variables
  integer :: algo,nproc,jproc
  integer(fmpi_integer) :: ierr,cnt
  integer(f_long) :: sizets,sizetr
  integer, dimension(:), allocatable :: nsenddspls_remote
  type(dictionary), pointer :: dict_info
  type(fmpi_win) :: window
  external :: MPI_ALLTOALL,MPI_ALLTOALLV

  nproc=mpisize(comm)
  sizets=f_size(sendbuf)
  sizetr=f_size(recvbuf)
  if (present(count)) then
     cnt=count
  else
     !division of the components per nproc
     cnt=int(min(sizets,sizetr)/nproc,fmpi_integer)
     if (cnt*nproc /= min(sizets,sizetr) .and. .not. present(sendcounts)) &
          call f_err_throw(&
          'The provided array has not a size which is a multiple of the number of processes',&
          err_id=ERR_MPI_WRAPPERS)
  end if
  !otherwise determine algorithm
  algo=alltoall_algorithm(nproc,sizets,sizetr,&
       cnt,sendcounts,sdispls,&
       cnt,recvcounts,rdispls,algorithm,comm)
  if (present(request)) request=FMPI_REQUEST_NULL
  if (present(win)) algo=VARIABLE_ONE_SIDED_GET_ALGO
  select case(algo)
  case(NOT_VARIABLE_ALGO)
     call MPI_ALLTOALL(sendbuf,cnt,mpitype(sendbuf), &
          recvbuf,cnt,mpitype(recvbuf),comm,ierr)
  case(VARIABLE_ALGO)
     if (present(request)) then
        ierr=FMPI_SUCCESS
        call MPI_IALLTOALLV(sendbuf,sendcounts,sdispls,mpitype(sendbuf), &
             recvbuf,recvcounts,rdispls,mpitype(recvbuf),comm,request,ierr)
     else
        call MPI_ALLTOALLV(sendbuf,sendcounts,sdispls,mpitype(sendbuf), &
             recvbuf,recvcounts,rdispls,mpitype(recvbuf),comm,ierr)
     end if
  case(VARIABLE_ONE_SIDED_GET_ALGO)
     nsenddspls_remote = f_malloc(0.to.nproc-1,id='nsenddspls_remote')
     call fmpi_alltoall(sendbuf=sdispls,recvbuf=nsenddspls_remote,comm=comm)
     !info=mpiinfo("no_locks", "true")
     dict_info=>dict_new('nolocks' .is. 'true')
     !window = mpiwindow(size(sendbuf), sendbuf, comm)
     call fmpi_win_create(window,sendbuf,size=sizets,dict_info=dict_info,comm=comm)
     call fmpi_win_fence(window,FMPI_WIN_OPEN)

     call dict_free(dict_info)
     do jproc=0,nproc-1
        if (recvcounts(jproc)>0) then
           call fmpi_get(origin_addr=recvbuf,origin_displ=rdispls(jproc),target_rank=jproc,&
                target_disp=int(nsenddspls_remote(jproc),fmpi_address),count=recvcounts(jproc),win=window)
        end if
     end do
     if (present(win)) then
        win=window
        !there should be no need to nullify window
     else
        call fmpi_win_shut(window)
     end if


     call f_free(nsenddspls_remote)
     ierr=FMPI_SUCCESS
  end select
  if (ierr/=FMPI_SUCCESS) then
     call f_err_throw('An error in calling to FMPI_ALLTOALL occured',&
          err_id=ERR_MPI_WRAPPERS)
     return
  end if



!old version
!!!  integer,dimension(0:),intent(in) :: sendcounts, sdispls
!!!  integer,dimension(0:),intent(in) :: recvcounts, rdispls
!!!  integer,intent(in) :: comm
!!!  character(len=*),intent(in),optional :: algorithm
!!!  ! Local variables
!!!  integer :: ierr, iproc, nproc, approach
!!!  logical :: large
!!!  integer,parameter :: USE_MPI_GET_ALLTOALLV = 0
!!!  integer,parameter :: USE_MPI_ALLTOALLV = 1
!!!  integer,parameter :: max_size = 100000000 !maximal number of elements sent and/or received on a proc using the standard MPI function
!!!
!!!  ! Check whether the sizes are correct
!!!  if (size(sendbuf)<sum(sendcounts)) then
!!!      call f_err_throw("insufficient size of array 'sendbuf': "//trim(yaml_toa(size(sendbuf)))//&
!!!           &" instead of "//trim(yaml_toa(sum(sendcounts))))
!!!  end if
!!!  if (size(recvbuf)<sum(recvcounts)) then
!!!      call f_err_throw("insufficient size of array 'recvbuf': "//trim(yaml_toa(size(recvbuf)))//&
!!!           &" instead of "//trim(yaml_toa(sum(recvcounts))))
!!!  end if
!!!
!!!  ! Check whether the approach to be used is imposed
!!!  if (present(algorithm)) then
!!!      select case (trim(algorithm))
!!!      case ('mpi_get')
!!!          approach = USE_MPI_GET_ALLTOALLV
!!!      case ('native')
!!!          approach = USE_MPI_ALLTOALLV
!!!      case default
!!!          call f_err_throw("wrong value of 'algorithm', use 'mpi_get' or 'native'")
!!!      end select
!!!  else
!!!      ! Check whether we are having a "big" case. If so, use the hand-made
!!!      ! version using mpi_get, otherwise the standard function.
!!!      large = .false.
!!!      if (sum(sendcounts)>max_size .or. sum(recvcounts)>max_size) then
!!!          large = .true.
!!!      end if
!!!      call mpiallred(large, 1, mpi_lor, comm=comm)
!!!      if (large) then
!!!          approach = USE_MPI_GET_ALLTOALLV
!!!          if (iproc==0) then
!!!              call yaml_warning('Large arrays in this call to mpi_alltoallv, &
!!!                   &using mpi_get_alltoallv instead')
!!!          end if
!!!      else
!!!          approach = USE_MPI_ALLTOALLV
!!!      end if
!!!  end if
!!!
!!!  ! Now call the corresponding MPI function
!!!  if (approach==USE_MPI_GET_ALLTOALLV) then
!!!
!!!      iproc = mpirank(comm)
!!!      nproc = mpisize(comm)
!!!      !call mpi_get_alltoallv(iproc, nproc, comm, sendcounts, sdispls, &
!!!      !     recvcounts, rdispls, sendbuf, recvbuf)
!!!      !!call c_f_pointer(c_loc(sendbuf), sendbuf_1d, [size(sendbuf)])
!!!      !!call c_f_pointer(c_loc(recvbuf), recvbuf_1d, [size(recvbuf)])
!!!      call mpi_get_alltoallv(iproc, nproc, comm, sendcounts, sdispls, &
!!!           recvcounts, rdispls, sendbuf, recvbuf)
!!!
!!!   else if (approach==USE_MPI_ALLTOALLV) then
!!!  
!!!      call mpi_alltoallv(sendbuf, sendcounts, sdispls, mpitype(sendbuf), &
!!!           recvbuf, recvcounts, rdispls, mpitype(recvbuf), comm, ierr)
!!!      if (ierr/=0) then
!!!         call f_err_throw('An error in calling to MPI_ALLTOALLV occured',&
!!!              err_id=ERR_MPI_WRAPPERS)
!!!         return
!!!      end if
!!!
!!!  else
!!!
!!!      call f_err_throw("wrong value of 'approach'")
!!!
!!!  end if
