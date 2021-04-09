!> @file
!! Wrapper for mpi_gather
!! Use error handling
!! @author
!!    Copyright (C) 2018-2018 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
module f_gather
  use time_profiling, only: TIMING_UNINITIALIZED
  use f_precisions
  use fmpi_types
  use yaml_strings
  use dictionaries, only: f_err_throw
  implicit none

  integer, public, save :: TCAT_GATHER       = TIMING_UNINITIALIZED
    
  private

  interface fmpi_gather
     module procedure mpigather_d0d2,mpigather_d1d1,mpigather_d1d2,mpigather_d2,mpigather_d2d1
     module procedure mpigather_i0i2,mpigather_i1,mpigather_i1i2,mpigather_i2
     module procedure mpigather_li0li2,mpigather_li1,mpigather_li1li2,mpigather_li2
     module procedure mpigather_c1i2,mpigather_c0c1
     module procedure mpigather_c1li2
  end interface fmpi_gather

  interface fmpi_gather_ptr
     module procedure mpigathered_d2
  end interface fmpi_gather_ptr

  public :: fmpi_gather,fmpi_gather_ptr

contains
  
  !> Gather the results of a given array into the root proc
  subroutine mpigather_d1d1(sendbuf,recvbuf,root,comm)
    use dictionaries, only: f_err_throw,f_err_define

    implicit none
    double precision, dimension(:), intent(in) :: sendbuf
    double precision, dimension(:), intent(inout) :: recvbuf
    include 'gather-inc.f90'
  end subroutine mpigather_d1d1

  subroutine mpigather_d1d2(sendbuf,recvbuf,root,comm)
    use dictionaries, only: f_err_throw,f_err_define

    implicit none
    double precision, dimension(:), intent(in) :: sendbuf
    double precision, dimension(:,:), intent(inout) :: recvbuf
    include 'gather-inc.f90'
  end subroutine mpigather_d1d2

  subroutine mpigather_i1i2(sendbuf,recvbuf,root,comm)
    use dictionaries, only: f_err_throw,f_err_define

    implicit none
    integer(f_integer), dimension(:), intent(in) :: sendbuf
    integer(f_integer), dimension(:,:), intent(inout) :: recvbuf
    include 'gather-inc.f90'
  end subroutine mpigather_i1i2

  subroutine mpigather_i2(sendbuf,recvbuf,root,comm)
    use dictionaries, only: f_err_throw,f_err_define

    implicit none
    integer(f_integer), dimension(:,:), intent(in) :: sendbuf
    integer(f_integer), dimension(:,:), intent(inout) :: recvbuf
    include 'gather-inc.f90'
  end subroutine mpigather_i2

  subroutine mpigather_i1(sendbuf,recvbuf,root,comm)
    use dictionaries, only: f_err_throw,f_err_define

    implicit none
    integer(f_integer), dimension(:), intent(in) :: sendbuf
    integer(f_integer), dimension(:), intent(inout) :: recvbuf
    include 'gather-inc.f90'
  end subroutine mpigather_i1

  subroutine mpigather_li1(sendbuf,recvbuf,root,comm)
    use dictionaries, only: f_err_throw,f_err_define

    implicit none
    integer(f_long), dimension(:), intent(in) :: sendbuf
    integer(f_long), dimension(:), intent(inout) :: recvbuf
    include 'gather-inc.f90'
  end subroutine mpigather_li1

  subroutine mpigather_li1li2(sendbuf,recvbuf,root,comm)
    use dictionaries, only: f_err_throw,f_err_define

    implicit none
    integer(f_long), dimension(:), intent(in) :: sendbuf
    integer(f_long), dimension(:,:), intent(inout) :: recvbuf
    include 'gather-inc.f90'
  end subroutine mpigather_li1li2

  subroutine mpigather_li2(sendbuf,recvbuf,root,comm)
    use dictionaries, only: f_err_throw,f_err_define

    implicit none
    integer(f_long), dimension(:,:), intent(in) :: sendbuf
    integer(f_long), dimension(:,:), intent(inout) :: recvbuf
    include 'gather-inc.f90'
  end subroutine mpigather_li2


  subroutine mpigather_d2d1(sendbuf,recvbuf,root,comm)
    use dictionaries, only: f_err_throw,f_err_define

    implicit none
    double precision, dimension(:,:), intent(in) :: sendbuf
    double precision, dimension(:), intent(inout) :: recvbuf
    include 'gather-inc.f90'
  end subroutine mpigather_d2d1

  subroutine mpigather_d2(sendbuf,recvbuf,root,comm)
    use dictionaries, only: f_err_throw,f_err_define

    implicit none
    double precision, dimension(:,:), intent(in) :: sendbuf
    double precision, dimension(:,:), intent(inout) :: recvbuf
    include 'gather-inc.f90'
  end subroutine mpigather_d2

  !> Gather the results of a given array into the root proc, version
  !! working with adresses
  subroutine mpigather_i0i2(sendbuf,sendcount,recvbuf,root,comm)
    use dictionaries, only: f_err_throw,f_err_define

    implicit none
    integer(f_integer), intent(inout) :: sendbuf
    integer, intent(in) :: sendcount
    integer(f_integer), dimension(:,:), intent(inout) :: recvbuf
    !---like gather-inc
    integer, intent(in), optional :: root !< 0 if absent
    integer, intent(in), optional :: comm !< MPI_COMM_WORLD if absent
    !local variables
    integer :: iroot,mpi_comm,ntot,ntotrecv,ntasks,ierr

    ntot=sendcount
    ntotrecv=size(recvbuf)

    include 'gather-inner-inc.f90'
    !-end gather-inc
  end subroutine mpigather_i0i2

  !> Gather the results of a given array into the root proc, version
  !! working with adresses
  subroutine mpigather_c1i2(sendbuf_c,recvbuf,root,comm)
    use dynamic_memory
    use dictionaries, only: f_err_throw,f_err_define

    implicit none
    character, dimension(:), intent(inout) :: sendbuf_c
    integer(f_integer), dimension(:,:), intent(inout) :: recvbuf
    !---like gather-inc
    integer, intent(in), optional :: root !< 0 if absent
    integer, intent(in), optional :: comm !< MPI_COMM_WORLD if absent
    !local variables
    integer :: iroot,mpi_comm,ntot,ntotrecv,ntasks,ierr
    integer(f_integer), dimension(:), allocatable :: sendbuf
    ntot=size(sendbuf_c)
    ntotrecv=size(recvbuf)
    sendbuf=f_malloc(ntot,id='sendbuf')
    call f_memcpy(src=sendbuf_c,dest=sendbuf)
    include 'gather-inner-inc.f90'
    call f_free(sendbuf)
    !-end gather-inc
  end subroutine mpigather_c1i2

  subroutine mpigather_c0c1(length,sendbuf,recvbuf,root,comm)
    use dynamic_memory
    use dictionaries, only: f_err_throw,f_err_define

    implicit none
    integer, intent(in) :: length
    character(len=length), intent(in) :: sendbuf
    character(len=length), dimension(:), intent(inout) :: recvbuf
    !---like gather-inc
    integer, intent(in), optional :: root !< 0 if absent
    integer, intent(in), optional :: comm !< MPI_COMM_WORLD if absent
    !local variables
    integer :: iroot,mpi_comm,ntot,ntotrecv,ntasks,ierr
    ntot=length
    ntotrecv=size(recvbuf)
    !include 'gather-inner-inc.f90'

    if (present(root)) then
       iroot=root
    else
       iroot=0
    end if
    mpi_comm=fmpi_comm(comm)

    !verify the size of the receive buffer
    ntasks=mpisize(mpi_comm)
    if (ntotrecv < ntasks) then
       call f_err_throw('Error in mpigather; the size of the receive buffer ('//&
            (ntotrecv)//&
            ') is not large enough to contain '//ntasks//&
            ' elements',err_id=ERR_MPI_WRAPPERS)
       return
    end if
    !then one can proceed with the MPI operation
    ntotrecv=int(int(ntot,f_long)*kind(sendbuf)/int(kind(recvbuf),f_long))
    call f_timer_interrupt(TCAT_GATHER)
    call MPI_GATHER(sendbuf,ntot,mpitype(sendbuf),&
         recvbuf,ntotrecv,mpitype(sendbuf),iroot,mpi_comm,ierr)
    call f_timer_resume()
    if (ierr /=0) then
       call f_err_throw('An error in calling to MPI_GATHER occured',&
            err_id=ERR_MPI_WRAPPERS)
       return
    end if

  end subroutine mpigather_c0c1


  !> Gather the results of a given array into the root proc, version
  !! working with adresses
  subroutine mpigather_c1li2(sendbuf_c,recvbuf,root,comm)
    use dynamic_memory
    use dictionaries, only: f_err_throw,f_err_define

    implicit none
    character, dimension(:), intent(inout) :: sendbuf_c
    integer(f_long), dimension(:,:), intent(inout) :: recvbuf
    !---like gather-inc
    integer, intent(in), optional :: root !< 0 if absent
    integer, intent(in), optional :: comm !< MPI_COMM_WORLD if absent
    !local variables
    integer :: iroot,mpi_comm,ntot,ntotrecv,ntasks,ierr
    integer(f_long), dimension(:), allocatable :: sendbuf
    ntot=size(sendbuf_c)
    ntotrecv=size(recvbuf)
    sendbuf=f_malloc(ntot,id='sendbuf')
    call f_memcpy(src=sendbuf_c,dest=sendbuf)
    include 'gather-inner-inc.f90'
    call f_free(sendbuf)
  end subroutine mpigather_c1li2

  !> Gather the results of a given array into the root proc, version
  !! working with adresses
  subroutine mpigather_li0li2(sendbuf,sendcount,recvbuf,root,comm)
    use dynamic_memory
    use dictionaries, only: f_err_throw,f_err_define

    implicit none
    integer(f_long), intent(inout) :: sendbuf
    integer, intent(in) :: sendcount
    integer(f_long), dimension(:,:), intent(inout) :: recvbuf
    !---like gather-inc
    integer, intent(in), optional :: root !< 0 if absent
    integer, intent(in), optional :: comm !< MPI_COMM_WORLD if absent
    !local variables
    integer :: iroot,mpi_comm,ntot,ntotrecv,ntasks,ierr

    ntot=sendcount
    ntotrecv=size(recvbuf)

    include 'gather-inner-inc.f90'
    !-end gather-inc
  end subroutine mpigather_li0li2


  subroutine mpigather_d0d2(sendbuf,sendcount,recvbuf,root,comm)
    use dynamic_memory
    use dictionaries, only: f_err_throw,f_err_define

    implicit none
    double precision, intent(inout) :: sendbuf
    integer, intent(in) :: sendcount
    double precision, dimension(:,:), intent(inout) :: recvbuf
    !---like gather-inc
    integer, intent(in), optional :: root !< 0 if absent
    integer, intent(in), optional :: comm !< MPI_COMM_WORLD if absent
    !local variables
    integer :: iroot,mpi_comm,ntot,ntotrecv,ntasks,ierr

    ntot=sendcount
    ntotrecv=size(recvbuf)

    include 'gather-inner-inc.f90'
    !-end gather-inc
  end subroutine mpigather_d0d2

  !> retrieve the buffers of each of the tasks, allocate a memory space in root
  !! as provide the pointer to this memory space
  function mpigathered_d2(sendbuf,root,comm) result(ptr)
    use f_utils
    use dynamic_memory
    implicit none
    real(f_double), dimension(:,:), intent(in) :: sendbuf
    integer, intent(in), optional :: comm,root
    real(f_double), dimension(:), pointer :: ptr
    !local variables
    integer :: root_,jproc,iproc,nproc,comm_,ierr
    integer, dimension(1) :: ncount_
    integer, dimension(:), allocatable :: ncounts,ndispls
    real(f_double), dimension(:,:), pointer :: recvbuf

    ncount_(1)=size(sendbuf)   
    iproc=mpirank(comm)
    nproc=mpisize(comm)

    ncounts=f_malloc0(nproc,id='ncounts')
    ndispls=f_malloc(nproc,id='ndispls')

    call fmpi_gather(sendbuf=ncount_,recvbuf=ncounts,root=root,comm=comm)

    ndispls(1)=0
    do jproc=2,nproc
       ndispls(jproc)=ndispls(jproc-1)+ncounts(jproc-1)
    end do
    root_=0
    if (present(root)) root_=root
    comm_=fmpi_comm(comm)

    !allocate the pointer with the good size
    if (iproc == root_) then
       ptr=f_malloc_ptr(sum(ncounts),id='ptr')
    else
       ptr=f_malloc_ptr(1,id='ptr')
    end if

    if (nproc==1) then
       call f_memcpy(src=sendbuf,dest=ptr)
    else
       !then perform the call to the gatherv routine
       call MPI_GATHERV(sendbuf,ncount_(1),mpitype(sendbuf),&
            ptr,ncounts,ndispls,mpitype(recvbuf),&
            root_,comm_,ierr)

       if (ierr /=0) then
          call f_err_throw('An error in calling to MPI_GATHERV occured',&
               err_id=ERR_MPI_WRAPPERS)
          return
       end if
    end if

    if (iproc /= root_) call f_free_ptr(ptr)

    call f_free(ncounts,ndispls)

  end function mpigathered_d2

end module f_gather
