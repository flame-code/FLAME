module f_allgather
  use time_profiling, only: TIMING_UNINITIALIZED
  use f_precisions
  use fmpi_types
  use yaml_strings
  use f_enums
  use f_onesided
  use dictionaries, only: f_err_throw
  use f_utils, only: f_size
  use dynamic_memory
  implicit none

  private

  logical, parameter :: have_mpi2=.true.

  integer, public, save :: TCAT_ALLGATHERV   = TIMING_UNINITIALIZED
  integer, public, save :: TCAT_ALLGATHER    = TIMING_UNINITIALIZED

  type, public :: fmpi_counts
     !fmpi structure for the counts and displacements
     !to be used in variable/nonvariable routines
     integer(fmpi_integer) :: count = int(-1,fmpi_integer)
     logical :: counts_allocated = .false.
     logical :: displs_allocated = .false.
     integer(fmpi_integer), dimension(:), pointer :: counts => null()
     integer(fmpi_integer), dimension(:), pointer :: displs => null()
  end type fmpi_counts


  !> Interface for MPI_ALLGATHERV routine
  interface fmpi_allgather
     module procedure mpiallgatherv_d0,mpiallgatherv_d1,mpiallgatherv_d2d3,mpiallgatherv_i2
     module procedure mpiallgatherv_i0,mpiallgatherv_i1!,mpiallgatherv_i0i1
  end interface fmpi_allgather

  interface mpi_get_to_allgatherv
     module procedure mpi_get_to_allgatherv_double
  end interface mpi_get_to_allgatherv

  public :: fmpi_allgather,mpi_get_to_allgatherv !to be moved inside allgather

contains

  subroutine get_ntot_and_offset(iproc,nproc,crt,ntot,offset)
    implicit none
    integer, intent(in) :: iproc,nproc
    type(fmpi_counts), intent(in) :: crt
    integer(f_long), intent(out) :: ntot,offset

    if (crt%count == -1) then
       ntot=sum(int(crt%counts,f_long))
       offset=int(crt%displs(iproc),f_long)
    else
       ntot=nproc*int(crt%count,f_long)
       offset=iproc*int(crt%count,f_long)
    end if
    
  end subroutine get_ntot_and_offset

  subroutine get_srcounts(me,buf_size,sendcnt,recvcnt,sendcount,recvcount,recvcounts)
    implicit none
    integer, intent(in) :: me
    integer(f_long), intent(in) :: buf_size
    integer, intent(out) :: sendcnt,recvcnt
    integer, intent(in), optional :: sendcount
    integer, intent(in), optional :: recvcount
    integer, dimension(:), intent(in), optional :: recvcounts
    
    if (present(recvcounts)) then
       recvcnt=recvcounts(me+1)
    elseif (present(recvcount)) then
       recvcnt=recvcount
    else if (present(sendcount)) then
       recvcnt=sendcount
    else
       recvcnt=int(buf_size)
    end if

    !determine sendcnt
    if (present(sendcount)) then
       sendcnt=sendcount
    else if (present(recvcount)) then
       sendcnt=recvcount
    else if (present(recvcounts)) then
       sendcnt=recvcounts(me+1)
    else
       sendcnt=int(buf_size)
    end if

  end subroutine get_srcounts


  !to be inserted
!!$  function allgather_algorithm(sendbuf,sendcount,recvbuf,recvcount,&
!!$       recvcounts,displs,comm,win,algorithm

  !>performs gathering of array portions into a receive buffer
  !! the arguments can be provided such as to meet either allgather
  !! or allgatherv APIs. The wrapper chooses the better routine to call
  !! as a function of the arguments
  subroutine mpiallgatherv_d0(sendbuf,sendcount,recvbuf,recvcount,&
       recvcounts,displs,comm,win,algorithm)
    !    use yaml_strings, only: yaml_toa
    use dictionaries, only: f_err_throw
    use dynamic_memory
    implicit none
    real(f_double), intent(in) :: sendbuf
    real(f_double),  optional :: recvbuf
    real(f_double), dimension(:), allocatable :: copybuf
    include 'allgather-inc.f90'
  end subroutine mpiallgatherv_d0
  subroutine mpiallgatherv_d1(sendbuf,sendcount,recvbuf,recvcount,&
       recvcounts,displs,comm,win,algorithm)
    !    use yaml_strings, only: yaml_toa
    use dictionaries, only: f_err_throw
    use dynamic_memory
    implicit none
    real(f_double), dimension(:), intent(in) :: sendbuf
    real(f_double), dimension(:), intent(inout), optional :: recvbuf
    real(f_double), dimension(:), allocatable :: copybuf
    include 'allgather-inc.f90'
  end subroutine mpiallgatherv_d1
  subroutine mpiallgatherv_d2d3(sendbuf,sendcount,recvbuf,recvcount,&
       recvcounts,displs,comm,win,algorithm)
    !    use yaml_strings, only: yaml_toa
    use dictionaries, only: f_err_throw
    use dynamic_memory
    implicit none
    real(f_double), dimension(:,:), intent(in) :: sendbuf
    real(f_double), dimension(:,:,:), intent(inout), optional :: recvbuf
    real(f_double), dimension(:), allocatable :: copybuf
    include 'allgather-inc.f90'
  end subroutine mpiallgatherv_d2d3

  subroutine mpiallgatherv_i0(sendbuf,sendcount,recvbuf,recvcount,&
       recvcounts,displs,comm,win,algorithm)
    !    use yaml_strings, only: yaml_toa
    use dictionaries, only: f_err_throw
    use dynamic_memory
    implicit none
    integer(f_integer), intent(in) :: sendbuf
    integer(f_integer), intent(inout), optional :: recvbuf
    integer(f_integer), dimension(:), allocatable :: copybuf
    include 'allgather-inc.f90'
  end subroutine mpiallgatherv_i0

  recursive subroutine mpiallgatherv_i1(sendbuf,sendcount,recvbuf,recvcount,&
       recvcounts,displs,comm,win,algorithm)
    !    use yaml_strings, only: yaml_toa
    use dictionaries, only: f_err_throw
    use dynamic_memory
    implicit none
    integer(f_integer), dimension(:), intent(in) :: sendbuf
    integer(f_integer), dimension(:), intent(inout), optional :: recvbuf
    integer(f_integer), dimension(:), allocatable :: copybuf
    include 'allgather-inc.f90'
  end subroutine mpiallgatherv_i1


!!$  subroutine mpiallgatherv_i0i1(sendbuf,sendcount,recvbuf,recvcount,&
!!$       recvcounts,displs,comm,win,algorithm)
!!$    !    use yaml_strings, only: yaml_toa
!!$    use dictionaries, only: f_err_throw
!!$    use dynamic_memory
!!$    implicit none
!!$    integer(f_integer), intent(in) :: sendbuf
!!$    integer(f_integer), dimension(:), intent(inout) :: recvbuf
!!$    integer(f_integer), dimension(:), allocatable :: copybuf
!!$    include 'allgather-inc.f90'
!!$  end subroutine mpiallgatherv_i0i1

  subroutine mpiallgatherv_i2(sendbuf,sendcount,recvbuf,recvcount,&
       recvcounts,displs,comm,win,algorithm)
    !    use yaml_strings, only: yaml_toa
    use dictionaries, only: f_err_throw
    use dynamic_memory
    implicit none
    integer, dimension(:,:), intent(in) :: sendbuf
    integer, dimension(:,:), intent(inout), optional :: recvbuf
    integer, dimension(:), allocatable :: copybuf
    include 'allgather-inc.f90'
  end subroutine mpiallgatherv_i2

!!$  subroutine mpiallgatherv_c2(sendbuf,sendbuf_len,sendcount,&
!!$       recvbuf,recvbuf_len,recvcount,&
!!$       recvcounts,displs,comm,win,algorithm)
!!$    !    use yaml_strings, only: yaml_toa
!!$    use dictionaries, only: f_err_throw
!!$    use dynamic_memory
!!$    implicit none
!!$    integer, intent(in) :: sendbuf_len,recvbuf_len
!!$    character(len=sendbuf_len), dimension(:,:), intent(in) :: sendbuf
!!$    character(len=recvbuf_len), dimension(:,:), intent(inout), optional :: recvbuf
!!$    character(len=sendbuf_len), dimension(:,:), allocatable :: copybuf
!!$    include 'allgather-inc.f90'
!!$  end subroutine mpiallgatherv_c2


  subroutine mpi_get_to_allgatherv_double(sendbuf,sendcount,recvbuf,recvcounts,displs,comm,window)
    use dictionaries, only: f_err_throw,f_err_define
    !    use yaml_strings, only: yaml_toa
    use f_onesided
    implicit none
    real(f_double) :: sendbuf
    real(f_double) :: recvbuf
    integer,dimension(:),intent(in) :: recvcounts, displs
    integer,intent(in) :: comm, sendcount
    type(fmpi_win), intent(out), optional :: window
    !local variables
    integer :: nproc,nrecvbuf
    !external :: getall
    type(fmpi_win) :: window_

    nproc=mpisize(comm)
    nrecvbuf=sum(recvcounts)

    !check coherence
    if (any([size(recvcounts),size(displs)] /= nproc)) then
       call f_err_throw("Error in get_to_gatherv, sizes not coherent with communicator"//&
            trim(yaml_toa([size(recvcounts),size(displs), nproc])),&
            err_id=ERR_MPI_WRAPPERS)
       return
    end if

    call fmpi_win_create(window_,sendbuf,sendcount,comm=comm)
    call fmpi_win_fence(window_,FMPI_WIN_OPEN)
    call getall_d(nproc,recvcounts,displs,window_,recvbuf)

    if (.not. present(window)) then
       call fmpi_win_shut(window_)
    else
       window=window_
    end if

  end subroutine mpi_get_to_allgatherv_double

  include 'control-routines-inc.f90'

end module f_allgather

!> used by get_to_allgatherv to pass the good addresses to the mpiget wrapper
!! should be generalized to multiple communications
subroutine getall_d(nproc,recvcounts,displs,window,recvbuffer)
  use f_onesided, only: fmpi_get,fmpi_win !wrapper_MPI, only: mpiget, mpi_address_kind
  use fmpi_types, only: fmpi_address
  use f_precisions
  implicit none
  type(fmpi_win), intent(in) :: window
  integer,intent(in) :: nproc!,nrecvbuffer!,window
  integer,dimension(0:nproc-1),intent(in) :: recvcounts,displs
  real(f_double),dimension(*),intent(out) :: recvbuffer
  !real(f_double) :: recvbuffer
  ! Local variables
  integer :: jproc, jcount, jst

  do jproc=0,nproc-1
     jcount=recvcounts(jproc)
     jst=displs(jproc)
     if (jcount>0) then
        call fmpi_get(recvbuffer(jst+1),jproc,window,jcount)
        !call fmpi_get(recvbuffer,jproc,window,jcount,origin_displ=jst)
        !call mpiget(recvbuffer(jst+1), jcount, jproc, int(0,kind=mpi_address_kind), window)
     end if
  end do

end subroutine getall_d
