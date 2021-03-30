!> @file
!! Wrapper for the basic types
!! @author
!!    Copyright (C) 2012-2016 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


module mpif_module
  !do not put implicit none to avoid implicit declaration of
  !datatypes in some MPI implementations
  include 'mpif.h'      !< MPI definitions and datatypes
end module mpif_module

module fmpi_types
  !renaming of the public constants
  use mpif_module, FMPI_IN_PLACE => MPI_IN_PLACE, FMPI_SUCCESS=>MPI_SUCCESS,&
       FMPI_REQUEST_NULL=>MPI_REQUEST_NULL,FMPI_INFO_NULL=>MPI_INFO_NULL,FMPI_WIN_NULL=>MPI_WIN_NULL,&
       FMPI_COMM_NULL=>MPI_COMM_NULL,FMPI_GROUP_NULL=>MPI_GROUP_NULL,&
       FMPI_MODE_NOPRECEDE=>MPI_MODE_NOPRECEDE,&
       FMPI_MODE_NOSUCCEED=>MPI_MODE_NOSUCCEED,FMPI_MODE_NOSTORE=>MPI_MODE_NOSTORE,&
       FMPI_MODE_NOPUT=>MPI_MODE_NOPUT,FMPI_STATUSES_IGNORE=>MPI_STATUSES_IGNORE,&
       FMPI_STATUS_IGNORE=>MPI_STATUS_IGNORE,FMPI_STATUS_SIZE => MPI_STATUS_SIZE,&
       FMPI_ANY_SOURCE=>MPI_ANY_SOURCE,FMPI_ANY_TAG=>MPI_ANY_TAG
  use f_precisions
  use f_enums
  use dictionaries, only: f_err_throw
  implicit none

  private

  !> Interface for MPITYPE routine, to be used in all the wrappers
  interface mpitype
     module procedure mpitype_i,mpitype_d,mpitype_r,mpitype_l,mpitype_c,mpitype_li
     module procedure mpitype_b,mpitype_b1,mpitype_b2,mpitype_b3
     module procedure mpitype_i1,mpitype_i2,mpitype_i3
     module procedure mpitype_l1,mpitype_l2,mpitype_l3
     !module procedure mpitype_ll3
     module procedure mpitype_r1,mpitype_r2,mpitype_r3,mpitype_r4
     module procedure mpitype_d1,mpitype_d2,mpitype_d3,mpitype_d4,mpitype_d5,mpitype_d6
     module procedure mpitype_c1,mpitype_c2
     module procedure mpitype_li1,mpitype_li2,mpitype_li3

  end interface mpitype

  interface mpitypesize
    module procedure mpitypesize_d0, mpitypesize_d1, mpitypesize_d2, mpitypesize_i0, mpitypesize_long0, mpitypesize_l0
    module procedure mpitypesize_i1, mpitypesize_li1, mpitypesize_i2
  end interface mpitypesize

  !might be included in a config.inc fine
  integer, parameter, public :: fmpi_integer=kind(MPI_INTEGER) !as MPI_INTEGER_KIND is a spec of C MPI
  integer, parameter, public :: fmpi_address=MPI_ADDRESS_KIND
  integer, parameter, public :: fmpi_offset=MPI_OFFSET_KIND

  !> Error codes
  integer, public, save :: ERR_MPI_WRAPPERS

  !qualification of enumerators
  !type(f_enumerator), public, target :: FMPI_OP=f_enumerator('MPI_OP',-10,null())

  !>enumerators of mpi ops
  type(f_enumerator), public :: FMPI_LAND=f_enumerator('MPI_LAND',int(MPI_LAND),null())
  type(f_enumerator), public :: FMPI_BAND=f_enumerator('MPI_BAND',int(MPI_BAND),null())
  type(f_enumerator), public :: FMPI_BOR=f_enumerator('MPI_BOR',int(MPI_BOR),null())
  type(f_enumerator), public :: FMPI_LOR=f_enumerator('MPI_LOR',int(MPI_LOR),null())
  type(f_enumerator), public :: FMPI_MAX=f_enumerator('MPI_MAX',int(MPI_MAX),null())
  type(f_enumerator), public :: FMPI_MIN=f_enumerator('MPI_MIN',int(MPI_MIN),null())
  type(f_enumerator), public :: FMPI_SUM=f_enumerator('MPI_SUM',int(MPI_SUM),null())

  !algorithms
  integer, parameter, public :: AUTOMATIC_ALGO=0
  integer, parameter, public :: NOT_VARIABLE_ALGO=10
  integer, parameter, public :: VARIABLE_ALGO=11
  integer, parameter, public :: VARIABLE_ONE_SIDED_GET_ALGO=12

  type(f_enumerator), public :: AUTOMATIC_ENUM =f_enumerator('AUTOMATIC',AUTOMATIC_ALGO,null())
  type(f_enumerator), public :: NOT_VARIABLE_ENUM =f_enumerator('ALLTOALL',NOT_VARIABLE_ALGO,null())
  type(f_enumerator), public :: VARIABLE_ENUM =f_enumerator('ALLTOALLV',VARIABLE_ALGO,null())
  type(f_enumerator), public :: ONESIDED_ENUM =f_enumerator('ALLTOALL_GET',VARIABLE_ONE_SIDED_GET_ALGO,null())


  !>enumerator of objects
  !type(f_enumerator), public :: FMPI_SUCCESS=f_enumerator('MPI_SUCCESS',int(MPI_SUCCESS),null())
  !type(f_enumerator), public :: FMPI_REQUEST_NULL=f_enumerator('MPI_REQUEST_NULL',int(MPI_REQUEST_NULL),null())

  !> while deciding the optimal strategy for handles leave them as integers
  public :: FMPI_IN_PLACE,FMPI_REQUEST_NULL,FMPI_SUCCESS,FMPI_INFO_NULL,FMPI_WIN_NULL,FMPI_COMM_NULL
  public :: FMPI_MODE_NOPUT,FMPI_MODE_NOPRECEDE,FMPI_MODE_NOSTORE,FMPI_MODE_NOSUCCEED,FMPI_GROUP_NULL
  public :: FMPI_STATUSES_IGNORE,FMPI_STATUS_IGNORE,FMPI_STATUS_SIZE
  public :: FMPI_ANY_TAG,FMPI_ANY_SOURCE

  public :: mpitype,mpirank,mpisize,fmpi_comm,mpitypesize,fmpi_barrier
  public :: fmpi_maxtag,mpirank_null,mpicomm_null

contains

  !> returns true if the mpi has been initialized
  function mpiinitialized()
    implicit none
    logical :: mpiinitialized
    !local variables
    logical :: flag
    integer(fmpi_integer) :: ierr

    mpiinitialized=.false.
    call mpi_initialized(flag, ierr)
    if (ierr /=FMPI_SUCCESS) then
       flag=.false.
       call f_err_throw('An error in calling to MPI_INITIALIZED occured',&
            err_id=ERR_MPI_WRAPPERS)
    end if
    mpiinitialized=flag

  end function mpiinitialized

  pure function fmpi_comm(comm) result(mpi_comm)
    implicit none
    integer(fmpi_integer), intent(in), optional :: comm
    integer(fmpi_integer) :: mpi_comm

    if (present(comm)) then
       mpi_comm=comm
    else
       mpi_comm=MPI_COMM_WORLD
    end if

  end function fmpi_comm

  !> Performs the barrier of a given communicator, if present
  subroutine fmpi_barrier(comm)
    use dictionaries, only: f_err_throw
    implicit none
    integer, intent(in), optional :: comm !< the communicator
    !local variables
    integer :: mpi_comm,ierr

    mpi_comm=fmpi_comm(comm)
    !call the barrier
    call MPI_BARRIER(mpi_comm,ierr)
    if (ierr /=0) then
       call f_err_throw('An error in calling to MPI_BARRIER occured',&
            err_id=ERR_MPI_WRAPPERS)
    end if
  end subroutine fmpi_barrier

  pure function mpirank_null() result(iproc)
    implicit none
    integer :: iproc
    iproc=MPI_PROC_NULL
  end function mpirank_null

  pure function mpicomm_null() result(comm)
    implicit none
    integer :: comm
    comm=MPI_PROC_NULL
  end function mpicomm_null

!!$  pure function mpirequest_null() result(request)
!!$    implicit none
!!$    integer :: request
!!$    request=MPI_REQUEST_NULL
!!$  end function mpirequest_null


  !> Function giving the mpi rank id for a given communicator
  function mpirank(comm)
    use dictionaries, only: f_err_throw
    implicit none
    integer, intent(in), optional :: comm
    integer :: mpirank
    !local variables
    integer(fmpi_integer) :: iproc,ierr,mpi_comm

    if (mpiinitialized()) then
       mpi_comm=fmpi_comm(comm)

       call MPI_COMM_RANK(mpi_comm, iproc, ierr)
       if (ierr /=FMPI_SUCCESS) then
          iproc=-1
          mpirank=iproc
          call f_err_throw('An error in calling to MPI_COMM_RANK occurred',&
               err_id=ERR_MPI_WRAPPERS)
       end if
       mpirank=iproc
    else
       mpirank=0
    end if

  end function mpirank

  !> Returns the number of mpi_tasks associated to a given communicator
  function mpisize(comm)
    use dictionaries, only: f_err_throw
    implicit none
    integer, intent(in), optional :: comm
    integer :: mpisize
    !local variables
    integer(fmpi_integer) :: nproc,ierr,mpi_comm

    if (mpiinitialized()) then
       mpi_comm=fmpi_comm(comm)

       !verify the size of the receive buffer
       call MPI_COMM_SIZE(mpi_comm,nproc,ierr)
       if (ierr /=FMPI_SUCCESS) then
          nproc=0
          mpisize=nproc
          call f_err_throw('An error in calling to MPI_COMM_SIZE occured',&
               err_id=ERR_MPI_WRAPPERS)
       end if
       mpisize=nproc
    else
       mpisize=1
    end if

  end function mpisize

  function fmpi_maxtag(comm) result(mpimaxtag)
    implicit none
    integer, intent(in), optional :: comm
    integer(kind=MPI_ADDRESS_KIND) :: mpimaxtag
    !local variables
    logical :: flag
    integer :: comm_,ierr

    comm_=fmpi_comm(comm)

    call MPI_COMM_GET_ATTR(comm_,MPI_TAG_UB,mpimaxtag,flag,ierr)

    if (ierr /= FMPI_SUCCESS) then
       call f_err_throw('An error in the first call to mpimaxtag occured',&
            err_id=ERR_MPI_WRAPPERS)
    end if

    !in case the flag has not been found with this communicator it is perhaps
    ! a global flag
    if (.not. flag .and. comm_ /= fmpi_comm()) then
      call MPI_COMM_GET_ATTR(fmpi_comm(),MPI_TAG_UB,mpimaxtag,flag,ierr)
    end if

    !error check
    if (ierr /= FMPI_SUCCESS .or. .not. flag) then
       call f_err_throw('An error in calling to mpimaxtag occured',&
            err_id=ERR_MPI_WRAPPERS)
    end if

  end function fmpi_maxtag


  function mpitypesize_d0(foo) result(sizeof)
    use dictionaries, only: f_err_throw,f_err_define
    implicit none
    real(f_double), intent(in) :: foo
    integer(fmpi_integer) :: sizeof, ierr
    integer :: kindt
    kindt=kind(foo) !to remove compilation warning

    call mpi_type_size(mpi_double_precision, sizeof, ierr)
    if (ierr/=0) then
       call f_err_throw('Error in mpi_type_size',&
            err_id=ERR_MPI_WRAPPERS)
    end if
  end function mpitypesize_d0

  function mpitypesize_d1(foo) result(sizeof)
    implicit none
    real(f_double), dimension(:), intent(in) :: foo
    integer :: sizeof
    integer :: kindt
    kindt=kind(foo) !to remove compilation warning
    sizeof=mpitypesize(1.0_f_double)
  end function mpitypesize_d1

  function mpitypesize_d2(foo) result(sizeof)
    implicit none
    real(f_double), dimension(:,:), intent(in) :: foo
    integer :: sizeof
    integer :: kindt
    kindt=kind(foo) !to remove compilation warning
    sizeof=mpitypesize(1.0_f_double)
  end function mpitypesize_d2


  function mpitypesize_i0(foo) result(sizeof)
    use dictionaries, only: f_err_throw,f_err_define
    implicit none
    integer, intent(in) :: foo
    integer(fmpi_integer) :: sizeof, ierr
    integer :: kindt
    kindt=kind(foo) !to remove compilation warning

    call mpi_type_size(mpi_integer, sizeof, ierr)
    if (ierr/=0) then
       call f_err_throw('Error in mpi_type_size',&
            err_id=ERR_MPI_WRAPPERS)
    end if
  end function mpitypesize_i0

  function mpitypesize_i1(foo) result(sizeof)
    implicit none
    integer(f_integer), dimension(:), intent(in) :: foo
    integer :: sizeof
    integer :: kindt
    kindt=kind(foo) !to remove compilation warning
    sizeof=mpitypesize(int(1,f_integer))
  end function mpitypesize_i1

  function mpitypesize_i2(foo) result(sizeof)
    implicit none
    integer(f_integer), dimension(:,:), intent(in) :: foo
    integer :: sizeof
    integer :: kindt
    kindt=kind(foo) !to remove compilation warning
    sizeof=mpitypesize(int(1,f_integer))
  end function mpitypesize_i2


  function mpitypesize_long0(foo) result(sizeof)
    use dictionaries, only: f_err_throw,f_err_define
    implicit none
    integer(f_long), intent(in) :: foo
    integer(fmpi_integer) :: sizeof, ierr
    integer :: kindt
    kindt=kind(foo) !to remove compilation warning

    !call mpi_type_size(MPI_LONG, sizeof, ierr)
    call mpi_type_size(MPI_INTEGER8, sizeof, ierr)
    if (ierr/=0) then
       call f_err_throw('Error in mpi_type_size',&
            err_id=ERR_MPI_WRAPPERS)
    end if
  end function mpitypesize_long0

  function mpitypesize_li1(foo) result(sizeof)
    implicit none
    integer(f_long), dimension(:), intent(in) :: foo
    integer :: sizeof
    integer :: kindt
    kindt=kind(foo) !to remove compilation warning
    sizeof=mpitypesize(int(1,f_long))
  end function mpitypesize_li1

  function mpitypesize_l0(foo) result(sizeof)
    use dictionaries, only: f_err_throw,f_err_define
    implicit none
    logical, intent(in) :: foo
    integer(fmpi_integer) :: sizeof, ierr
    integer :: kindt
    kindt=kind(foo) !to remove compilation warning

    call mpi_type_size(mpi_logical, sizeof, ierr)
    if (ierr/=0) then
       call f_err_throw('Error in mpi_type_size',&
            err_id=ERR_MPI_WRAPPERS)
    end if
  end function mpitypesize_l0



  pure function mpitype_i(data) result(mt)
    implicit none
    integer(f_integer), intent(in) :: data
    integer(fmpi_integer) :: mt
    integer :: kindt
    kindt=kind(data) !to remove compilation warning
    mt=MPI_INTEGER
  end function mpitype_i
  pure function mpitype_i1(data) result(mt)
    implicit none
    integer(f_integer), dimension(:), intent(in) :: data
    integer(fmpi_integer) :: mt
    integer :: kindt
    kindt=kind(data) !to remove compilation warning
    mt=MPI_INTEGER
  end function mpitype_i1
  pure function mpitype_i2(data) result(mt)
    implicit none
    integer(f_integer), dimension(:,:), intent(in) :: data
    integer(fmpi_integer) :: mt
    integer :: kindt
    kindt=kind(data) !to remove compilation warning
    mt=MPI_INTEGER
  end function mpitype_i2
  pure function mpitype_i3(data) result(mt)
    implicit none
    integer(f_integer), dimension(:,:,:), intent(in) :: data
    integer(fmpi_integer) :: mt
    integer :: kindt
    kindt=kind(data) !to remove compilation warning
    mt=MPI_INTEGER
  end function mpitype_i3

  pure function mpitype_l1(data) result(mt)
    implicit none
    logical, dimension(:), intent(in) :: data
    integer(fmpi_integer) :: mt
    integer :: kindt
    kindt=kind(data) !to remove compilation warning
    mt=MPI_LOGICAL
  end function mpitype_l1

  pure function mpitype_l2(data) result(mt)
    implicit none
    logical, dimension(:,:), intent(in) :: data
    integer(fmpi_integer) :: mt
    integer :: kindt
    kindt=kind(data) !to remove compilation warning
    mt=MPI_LOGICAL
  end function mpitype_l2

  pure function mpitype_l3(data) result(mt)
    implicit none
    logical, dimension(:,:,:), intent(in) :: data
    integer(fmpi_integer) :: mt
    integer :: kindt
    kindt=kind(data) !to remove compilation warning
    mt=MPI_LOGICAL
  end function mpitype_l3

  !pure function mpitype_ll3(data) result(mt)
  !  implicit none
  !  logical(f_byte), dimension(:,:,:), intent(in) :: data
  !  integer(fmpi_integer) :: mt
  !  integer :: kindt
  !  kindt=kind(data) !to remove compilation warning
  !  mt=MPI_LOGICAL
  !end function mpitype_ll3

  pure function mpitype_b1(data) result(mt)
    implicit none
    logical(f_byte), dimension(:), intent(in) :: data
    integer(fmpi_integer) :: mt
    integer :: kindt
    kindt=kind(data) !to remove compilation warning
    mt=MPI_BYTE
  end function mpitype_b1

  pure function mpitype_b2(data) result(mt)
    implicit none
    logical(f_byte), dimension(:,:), intent(in) :: data
    integer(fmpi_integer) :: mt
    integer :: kindt
    kindt=kind(data) !to remove compilation warning
    mt=MPI_BYTE
  end function mpitype_b2

  pure function mpitype_b3(data) result(mt)
    implicit none
    logical(f_byte), dimension(:,:,:), intent(in) :: data
    integer(fmpi_integer) :: mt
    integer :: kindt
    kindt=kind(data) !to remove compilation warning
    mt=MPI_BYTE
  end function mpitype_b3


  pure function mpitype_li(data) result(mt)
    implicit none
    integer(f_long), intent(in) :: data
    integer(fmpi_integer) :: mt
    integer :: kindt
    kindt=kind(data) !to remove compilation warning
    mt=MPI_INTEGER8
  end function mpitype_li
  pure function mpitype_li1(data) result(mt)
    implicit none
    integer(f_long), dimension(:), intent(in) :: data
    integer(fmpi_integer) :: mt
    integer :: kindt
    kindt=kind(data) !to remove compilation warning
    mt=MPI_INTEGER8
  end function mpitype_li1
  pure function mpitype_li2(data) result(mt)
    implicit none
    integer(f_long), dimension(:,:), intent(in) :: data
    integer(fmpi_integer) :: mt
    integer :: kindt
    kindt=kind(data) !to remove compilation warning
    mt=MPI_INTEGER8
  end function mpitype_li2
  pure function mpitype_li3(data) result(mt)
    implicit none
    integer(f_long), dimension(:,:,:), intent(in) :: data
    integer(fmpi_integer) :: mt
    integer :: kindt
    kindt=kind(data) !to remove compilation warning
    mt=MPI_INTEGER8
  end function mpitype_li3


  pure function mpitype_r(data) result(mt)
    implicit none
    real, intent(in) :: data
    integer(fmpi_integer) :: mt
    integer :: kindt
    kindt=kind(data) !to remove compilation warning
    mt=MPI_REAL
  end function mpitype_r
  pure function mpitype_d(data) result(mt)
    implicit none
    real(f_double), intent(in) :: data
    integer(fmpi_integer) :: mt
    integer :: kindt
    kindt=kind(data) !to remove compilation warning
    mt=MPI_DOUBLE_PRECISION
  end function mpitype_d
  pure function mpitype_d1(data) result(mt)
    implicit none
    real(f_double), dimension(:), intent(in) :: data
    integer(fmpi_integer) :: mt
    integer :: kindt
    kindt=kind(data) !to remove compilation warning
    mt=MPI_DOUBLE_PRECISION
  end function mpitype_d1
  pure function mpitype_d2(data) result(mt)
    implicit none
    real(f_double), dimension(:,:), intent(in) :: data
    integer(fmpi_integer) :: mt
    integer :: kindt
    kindt=kind(data) !to remove compilation warning
    mt=MPI_DOUBLE_PRECISION
  end function mpitype_d2
  pure function mpitype_d3(data) result(mt)
    implicit none
    real(f_double), dimension(:,:,:), intent(in) :: data
    integer(fmpi_integer) :: mt
    integer :: kindt
    kindt=kind(data) !to remove compilation warning
    mt=MPI_DOUBLE_PRECISION
  end function mpitype_d3
  pure function mpitype_d4(data) result(mt)
    implicit none
    real(f_double), dimension(:,:,:,:), intent(in) :: data
    integer(fmpi_integer) :: mt
    integer :: kindt
    kindt=kind(data) !to remove compilation warning
    mt=MPI_DOUBLE_PRECISION
  end function mpitype_d4
  pure function mpitype_d5(data) result(mt)
    implicit none
    real(f_double), dimension(:,:,:,:,:), intent(in) :: data
    integer(fmpi_integer) :: mt
    integer :: kindt
    kindt=kind(data) !to remove compilation warning
    mt=MPI_DOUBLE_PRECISION
  end function mpitype_d5
  pure function mpitype_d6(data) result(mt)
    implicit none
    real(f_double), dimension(:,:,:,:,:,:), intent(in) :: data
    integer(fmpi_integer) :: mt
    integer :: kindt
    kindt=kind(data) !to remove compilation warning
    mt=MPI_DOUBLE_PRECISION
  end function mpitype_d6

  pure function mpitype_r1(data) result(mt)
    implicit none
    real, dimension(:), intent(in) :: data
    integer(fmpi_integer) :: mt
    integer :: kindt
    kindt=kind(data) !to remove compilation warning
    mt=MPI_REAL
  end function mpitype_r1
  pure function mpitype_r2(data) result(mt)
    implicit none
    real, dimension(:,:), intent(in) :: data
    integer(fmpi_integer) :: mt
    integer :: kindt
    kindt=kind(data) !to remove compilation warning
    mt=MPI_REAL
  end function mpitype_r2
  pure function mpitype_r3(data) result(mt)
    implicit none
    real, dimension(:,:,:), intent(in) :: data
    integer(fmpi_integer) :: mt
    integer :: kindt
    kindt=kind(data) !to remove compilation warning
    mt=MPI_REAL
  end function mpitype_r3
  pure function mpitype_r4(data) result(mt)
    implicit none
    real, dimension(:,:,:,:), intent(in) :: data
    integer(fmpi_integer) :: mt
    integer :: kindt
    kindt=kind(data) !to remove compilation warning
    mt=MPI_REAL
  end function mpitype_r4

  pure function mpitype_b(data) result(mt)
    implicit none
    logical(f_byte), intent(in) :: data
    integer(fmpi_integer) :: mt
    integer :: kindt
    kindt=kind(data) !to remove compilation warning
    mt=MPI_BYTE
  end function mpitype_b

  pure function mpitype_l(data) result(mt)
    implicit none
    logical, intent(in) :: data
    integer(fmpi_integer) :: mt
    integer :: kindt
    kindt=kind(data) !to remove compilation warning
    mt=MPI_LOGICAL
  end function mpitype_l
  pure function mpitype_c(data) result(mt)
    implicit none
    character, intent(in) :: data
    integer(fmpi_integer) :: mt
    integer :: kindt
    kindt=kind(data) !to remove compilation warning
    mt=MPI_CHARACTER
  end function mpitype_c
  pure function mpitype_c1(data) result(mt)
    implicit none
    character, dimension(:), intent(in) :: data
    integer(fmpi_integer) :: mt
    integer :: kindt
    kindt=kind(data) !to remove compilation warning
    mt=MPI_CHARACTER
  end function mpitype_c1
  pure function mpitype_c2(data) result(mt)
    implicit none
    character, dimension(:,:), intent(in) :: data
    integer(fmpi_integer) :: mt
    integer :: kindt
    kindt=kind(data) !to remove compilation warning
    mt=MPI_CHARACTER
  end function mpitype_c2


end module fmpi_types
