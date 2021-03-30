!> @file
!! Wrapper for mpi_alltoall flavours
!! Use error handling
!! @author
!!    Copyright (C) 2012-2016 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
module f_alltoall
  use f_enums
  use f_precisions
  use fmpi_types
  use dictionaries, only: f_err_throw
  use f_allreduce
  use f_onesided
  implicit none

  private

  interface fmpi_alltoall
     module procedure mpialltoallw_d11,mpialltoallw_i11
     module procedure mpialltoallw_li11
     module procedure mpialltoallw_d22,mpialltoallw_d12
     !module procedure mpialltoallw_d00,mpialltoallw_i00
  end interface fmpi_alltoall

!!$  interface mpialltoallv
!!$     module procedure mpialltoallv_int11, mpialltoallv_long11, mpialltoallv_double11
!!$     !module procedure mpialltoallv_double61
!!$  end interface mpialltoallv
!!$
!!$  interface mpiialltoallv
!!$     module procedure mpiialltoallv_double
!!$  end interface mpiialltoallv
!!$
!!$  interface mpi_get_alltoallv
!!$     module procedure mpi_get_alltoallv_i, mpi_get_alltoallv_l, mpi_get_alltoallv_d
!!$  end interface mpi_get_alltoallv

  public :: fmpi_alltoall

  contains

    !>choose the correct algorithm to be applied as a function of the sendrecv arrays
    function alltoall_algorithm(nproc,sizets,sizetr,&
         sendcount,sendcounts,sdispls,&
         recvcount,recvcounts,rdispls,algorithm,comm) result(algo)
      use yaml_strings
      implicit none
      integer, intent(in) :: nproc
      integer(f_long), intent(in) :: sizets,sizetr
      integer, intent(in), optional :: sendcount,recvcount
      integer, dimension(0:), intent(in), optional :: sendcounts,sdispls
      integer, dimension(0:), intent(in), optional :: recvcounts,rdispls
      type(f_enumerator), intent(in), optional :: algorithm
      integer, intent(in), optional :: comm
      integer :: algo
      !local variables
      !>maximal number of elements sent and/or received on a proc using the standard MPI function
      integer, parameter :: max_size = 100000000
      logical :: large
      integer(f_long) :: sizecomms,sizecommr

      sizecomms=0
      !priority to variable version
      if (present(sendcounts)) then
         sizecomms=sum(int(sendcounts,f_long))
      else if (present(sendcount)) then
         sizecomms=int(sendcount,f_long)*nproc
      end if
      sizecommr=sizecomms

      if (present(recvcounts)) then
         sizecommr=sum(int(recvcounts,f_long))
      else if (present(recvcount)) then
         sizecommr=int(recvcount,f_long)*nproc
      end if

      !check the validity of the buffers
      if (sizets < sizecomms .and. sizets /=0) call f_err_throw(&
           'Size of send buffer and of data are not consistent ('&
           //trim(yaml_toa([sizets,sizecomms]))//')',&
           err_name='ERR_MPI_WRAPPERS')

      !check the validity of the buffers
      if (sizetr < sizecommr .and. sizetr /=0) call f_err_throw(&
           'Size of recv buffer and of data are not consistent ('&
           //trim(yaml_toa([sizetr,sizecommr]))//')',&
           err_name='ERR_MPI_WRAPPERS')

      algo=AUTOMATIC_ALGO
      if (present(sendcounts) .and. present(sdispls) .and. present(recvcounts) .and. present(rdispls)) then
         algo=toi(VARIABLE_ENUM)
      else if(present(sendcount) .and. present(recvcount)) then
         algo=toi(NOT_VARIABLE_ENUM)
      else
         call f_err_throw('Illegal arguments for alltoall',err_id=ERR_MPI_WRAPPERS)
      end if

      if (present(algorithm)) then
         !if the algorithm has been already provided check if it is possible
         select case(toi(algorithm))
         case(AUTOMATIC_ALGO)
            !nothing to do here
         case(NOT_VARIABLE_ALGO)
            !check if the proposition is compatible
            if (algo == toi(VARIABLE_ENUM)) &
                 call f_err_throw('Algorithm not compatible with arguments',err_id=ERR_MPI_WRAPPERS)
            return !ok
         case(VARIABLE_ALGO,VARIABLE_ONE_SIDED_GET_ALGO)
            !check if the proposition is compatible
            if (algo == toi(NOT_VARIABLE_ENUM)) &
                 call f_err_throw('Algorithm (variable) not compatible with arguments',err_id=ERR_MPI_WRAPPERS)
            if (algorithm == ONESIDED_ENUM) algo=VARIABLE_ONE_SIDED_GET_ALGO
            return !ok
         case default
            call f_err_throw('Illegal value for algorithm',err_id=ERR_MPI_WRAPPERS)
         end select
      end if

      if (algo==VARIABLE_ALGO) then !here we are in the automatic case
         ! Check whether we are having a "big" case. If so, use the hand-made
         ! version using mpi_get, otherwise the standard function.
         !large = (sizets>max_size .or. sizetr>max_size) .and. associated(sdispls_) .and. associated(rdispls_)
         large = (sizecommr>max_size .or. sizecomms>max_size)
         call fmpi_allreduce(large, 1,FMPI_LOR, comm=comm)
         if (large) algo=VARIABLE_ONE_SIDED_GET_ALGO
      end if

    end function alltoall_algorithm

!!$    recursive subroutine mpialltoallw_i00(sendbuf, recvbuf, &
!!$         count, sendcounts, sdispls, recvcounts, rdispls, comm, request, win, algorithm)
!!$      use dynamic_memory
!!$      use f_utils, only: f_size
!!$      use dictionaries
!!$      implicit none
!!$      integer(f_integer) ::  sendbuf
!!$      integer(f_integer) :: recvbuf
!!$      include 'alltoallv-inc.f90'
!!$    end subroutine mpialltoallw_i00


    recursive subroutine mpialltoallw_i11(sendbuf, recvbuf, &
         count, sendcounts, sdispls, recvcounts, rdispls, comm, request, win, algorithm)
      use dynamic_memory
      use f_utils, only: f_size
      use dictionaries
      implicit none
      integer(f_integer), dimension(:),intent(in)  :: sendbuf
      integer(f_integer), dimension(:),intent(out) :: recvbuf
      include 'alltoallv-inc.f90'
    end subroutine mpialltoallw_i11

    recursive subroutine mpialltoallw_li11(sendbuf, recvbuf, &
         count, sendcounts, sdispls, recvcounts, rdispls, comm, request, win, algorithm)
      use dynamic_memory
      use f_utils, only: f_size
      use dictionaries
      implicit none
      integer(f_long), dimension(:),intent(in)  :: sendbuf
      integer(f_long), dimension(:),intent(out) :: recvbuf
      include 'alltoallv-inc.f90'
    end subroutine mpialltoallw_li11

    subroutine mpialltoallw_d11(sendbuf, recvbuf, &
         count, sendcounts, sdispls, recvcounts, rdispls, comm, request, win, algorithm)
      use dynamic_memory
      use f_utils, only: f_size
      use dictionaries
      implicit none
      real(f_double), dimension(:),intent(in)  :: sendbuf
      real(f_double), dimension(:),intent(out) :: recvbuf
      include 'alltoallv-inc.f90'
    end subroutine mpialltoallw_d11

    subroutine mpialltoallw_d22(sendbuf, recvbuf, &
         count, sendcounts, sdispls, recvcounts, rdispls, comm, request, win, algorithm)
      use dynamic_memory
      use f_utils, only: f_size
      use dictionaries
      implicit none
      real(f_double), dimension(:,:),intent(in)  :: sendbuf
      real(f_double), dimension(:,:),intent(out) :: recvbuf
      include 'alltoallv-inc.f90'
    end subroutine mpialltoallw_d22

    subroutine mpialltoallw_d12(sendbuf, recvbuf, &
         count, sendcounts, sdispls, recvcounts, rdispls, comm, request, win, algorithm)
      use dynamic_memory
      use f_utils, only: f_size
      use dictionaries
      implicit none
      real(f_double), dimension(:),intent(in)  :: sendbuf
      real(f_double), dimension(:,:),intent(out) :: recvbuf
      include 'alltoallv-inc.f90'
    end subroutine mpialltoallw_d12

!!$    subroutine mpialltoallw_d00(sendbuf, recvbuf, &
!!$         count, sendcounts, sdispls, recvcounts, rdispls, comm, request, win, algorithm)
!!$      use dynamic_memory
!!$      use f_utils, only: f_size
!!$      use dictionaries
!!$      implicit none
!!$      real(f_double) :: sendbuf
!!$      real(f_double) :: recvbuf
!!$      include 'alltoallv-inc.f90'
!!$    end subroutine mpialltoallw_d00


!!$    subroutine mpialltoallv_long11(sendbuf, sendcounts, sdispls, recvbuf, recvcounts, rdispls, comm, algorithm)
!!$      use dictionaries, only: f_err_throw,f_err_define
!!$      use dynamic_memory
!!$      use yaml_output
!!$      !use iso_c_binding
!!$      implicit none
!!$      integer(f_long),dimension(:),intent(in),target :: sendbuf
!!$      integer(f_long),dimension(:),intent(out),target :: recvbuf
!!$      integer(f_long),dimension(:),pointer :: sendbuf_1d
!!$      integer(f_long),dimension(:),pointer :: recvbuf_1d
!!$      include 'alltoallv-inc.f90'
!!$    end subroutine mpialltoallv_long11
!!$
!!$    subroutine mpialltoallv_double11(sendbuf, sendcounts, sdispls, recvbuf, recvcounts, rdispls, comm, algorithm)
!!$      use dictionaries, only: f_err_throw,f_err_define
!!$      use dynamic_memory
!!$      use yaml_output
!!$      implicit none
!!$      double precision,dimension(:),intent(in),target :: sendbuf
!!$      double precision,dimension(:),intent(out),target :: recvbuf
!!$      double precision,dimension(:),pointer :: sendbuf_1d
!!$      double precision,dimension(:),pointer :: recvbuf_1d
!!$      include 'alltoallv-inc.f90'
!!$    end subroutine mpialltoallv_double11
!!$
!!$    !!subroutine mpialltoallv_double61(sendbuf, sendcounts, sdispls, recvbuf, recvcounts, rdispls, comm)
!!$    !!  use dictionaries, only: f_err_throw,f_err_define
!!$    !!  use dynamic_memory
!!$    !!  use yaml_output
!!$    !!  use iso_c_binding
!!$    !!  implicit none
!!$    !!  double precision,dimension(:,:,:,:,:,:),intent(in),target :: sendbuf
!!$    !!  double precision,dimension(:),intent(out),target :: recvbuf
!!$    !!  double precision,dimension(:),pointer :: sendbuf_1d
!!$    !!  double precision,dimension(:),pointer :: recvbuf_1d
!!$    !!  include 'alltoallv-inc.f90'
!!$    !!end subroutine mpialltoallv_double61
!!$
!!$    subroutine mpiialltoallv_double(sendbuf, sendcounts, senddspls, sendtype, &
!!$         recvbuf, recvcounts, recvdspls, recvtype, comm, request)
!!$      use dictionaries, only: f_err_throw,f_err_define
!!$      implicit none
!!$      ! Calling arguments
!!$      integer,intent(in) :: sendcounts, senddspls, sendtype, recvcounts, recvdspls, recvtype, comm
!!$      double precision,intent(in) :: sendbuf
!!$      double precision,intent(out) :: recvbuf
!!$      integer,intent(out) :: request
!!$      ! Local variables
!!$      integer :: ierr
!!$
!!$#ifdef HAVE_MPI3
!!$      call mpi_ialltoallv(sendbuf, sendcounts, senddspls, sendtype, &
!!$           recvbuf, recvcounts, recvdspls, recvtype, comm, request, ierr)
!!$      if (ierr/=0) then
!!$         call f_err_throw('An error in calling to MPI_IALLTOALLV occured',&
!!$              err_id=ERR_MPI_WRAPPERS)
!!$         return
!!$      end if
!!$#else
!!$      call mpi_alltoallv(sendbuf, sendcounts, senddspls, sendtype, &
!!$           recvbuf, recvcounts, recvdspls, recvtype, comm, ierr)
!!$      if (ierr/=0) then
!!$         call f_err_throw('An error in calling to MPI_IALLTOALLV occured',&
!!$              err_id=ERR_MPI_WRAPPERS)
!!$         return
!!$      end if
!!$      request = MPI_REQUEST_NULL
!!$#endif
!!$
!!$    end subroutine mpiialltoallv_double
!!$
!!$    subroutine mpi_get_alltoallv_i(iproc, nproc, comm, nsendcounts, nsenddspls, nrecvcounts, nrecvdspls, sendbuf, recvbuf)
!!$      use dynamic_memory
!!$      implicit none
!!$      integer(f_integer),dimension(:),intent(in) :: sendbuf
!!$      integer(f_integer),dimension(:),intent(in) :: recvbuf
!!$
!!$      include 'mpi_get_alltoallv-inc.f90'
!!$
!!$    end subroutine mpi_get_alltoallv_i
!!$
!!$    subroutine mpi_get_alltoallv_l(iproc, nproc, comm, nsendcounts, nsenddspls, nrecvcounts, nrecvdspls, sendbuf, recvbuf)
!!$      use dynamic_memory
!!$      implicit none
!!$      integer(f_long),dimension(:),intent(in) :: sendbuf
!!$      integer(f_long),dimension(:),intent(in) :: recvbuf
!!$
!!$      include 'mpi_get_alltoallv-inc.f90'
!!$
!!$    end subroutine mpi_get_alltoallv_l
!!$
!!$    subroutine mpi_get_alltoallv_d(iproc, nproc, comm, nsendcounts, nsenddspls, nrecvcounts, nrecvdspls, sendbuf, recvbuf)
!!$      use dynamic_memory
!!$      implicit none
!!$      double precision,dimension(:),intent(in) :: sendbuf
!!$      double precision,dimension(:),intent(in) :: recvbuf
!!$
!!$      include 'mpi_get_alltoallv-inc.f90'
!!$
!!$    end subroutine mpi_get_alltoallv_d
!!$
!!$
end module f_alltoall
