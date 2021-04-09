!> @file
!! Control operations for counts and displacements
!! @author
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
!!$module f_control
!!$  use fmpi_types
!!$  use f_precisions
!!$  use dictionaries, only: f_err_throw
!!$  implicit none
!!$
!!$  private
!!$
!!$  type, public :: fmpi_counts
!!$     !fmpi structure for the counts and displacements
!!$     !to be used in variable/nonvariable routines
!!$     integer(fmpi_integer) :: count = int(-1,fmpi_integer)
!!$     logical :: counts_allocated = .false.
!!$     logical :: displs_allocated = .false.
!!$     integer(fmpi_integer), dimension(:), pointer :: counts => null()
!!$     integer(fmpi_integer), dimension(:), pointer :: displs => null()
!!$  end type fmpi_counts
!!$
!!$  public :: get_control_buffers
!!$
!!$contains

  pure function fmpi_counts_null() result(cnt)
    implicit none
    type(fmpi_counts) :: cnt
  end function fmpi_counts_null

  subroutine free_fmpi_counts(cnt)
    implicit none
    type(fmpi_counts), intent(inout) :: cnt

    if (cnt%counts_allocated) nullify(cnt%counts)
    if (cnt%displs_allocated) nullify(cnt%displs)
    call f_free_ptr(cnt%counts)
    call f_free_ptr(cnt%displs)
    cnt%counts_allocated=.false.
    cnt%displs_allocated=.false.

  end subroutine free_fmpi_counts

  function get_control_buffers(buf_size,count,counts,displs,comm,control) result(ctr)
    implicit none
    integer(f_long), intent(in) :: buf_size
    logical, intent(in), optional :: control
    integer, intent(in), optional :: count,comm
    integer, dimension(0:), intent(in), optional, target :: counts,displs
    type(fmpi_counts) :: ctr
    !local variables
    integer :: ourselves,me
    logical :: variable

    ctr=fmpi_counts_null()
    if (present(count) .or. present(counts)) then
       !extract from the passed arguments the counts and displacements to deal with
       ctr%counts=>get_counts(count,counts,comm,control)
    else
       !if nothing is specified as counts use the buffer size and gather it
       ourselves=mpisize(comm=comm)
       me=mpirank(comm=comm)
       ctr%counts=f_malloc_ptr(0.to.ourselves-1,id='count_')
       ctr%counts(me)=int(buf_size)
       call fmpi_allgather(ctr%counts,recvcount=1,comm=comm)
    end if
    ctr%counts_allocated=present(counts)
    ctr%displs=>get_displacements(ctr%counts,displs)
    ctr%displs_allocated=present(displs)

    !verify if we are in a variable regime or not
    !now determine if counts are variable or not
    ctr%count=ctr%counts(0)
    variable=any(ctr%counts /= ctr%count)
    if (.not. variable) variable=get_variable_displacements(ctr%count,ctr%displs)
    !in the case of not variable just nullify the pointers and keep count
    if (.not. variable) then
       call free_fmpi_counts(ctr)
    else
       ctr%count=-1
    end if

  end function get_control_buffers

  function get_counts(count,counts,comm,control) result(counts_)
    use f_bcast, only: fmpi_maxdiff
    use dynamic_memory
    implicit none
    logical, intent(in), optional :: control
    integer, optional :: count !intent in implicit
    integer, intent(in), optional :: comm
    integer, dimension(0:), optional, target :: counts !intent in implicit
    integer, dimension(:), pointer :: counts_
    !local variables
    logical :: variable,ctrl
    integer :: ourselves,mdiff

    nullify(counts_)
    ourselves=mpisize(comm=comm)
    ctrl=.false.
    if (present(control)) ctrl=control
    if(present(counts)) then
       if (size(counts) /= ourselves) then
          call f_err_throw('Size of counts is not of length group size',&
               err_name='ERR_MPI_WRAPPERS')
          return
       end if
       if (ctrl) then
          mdiff=fmpi_maxdiff(counts,source=0,comm=comm)
          if (mdiff /= 0) then
             call f_err_throw('The counts array is not identical for all the processes, cannot continue',&
                  err_name='ERR_MPI_WRAPPERS')
             return
          end if
       end if
       counts_=>counts
    else if (present(count)) then
       counts_=f_malloc_ptr(0.to.ourselves-1,id='count_')
       if (ctrl) then
          !todo: control to be activated
          !call fmpi_allgather(sendbuf=count,sendcount=1,recvbuf=counts_,comm=comm)
       else
          counts_=count
       end if
    else
       call f_err_throw('At least count or counts should be present',&
            err_name='ERR_MPI_WRAPPERS')
    end if

  end function get_counts

  function get_displacements(counts,displs) result(displs_)
    use yaml_strings
    use dynamic_memory
    implicit none
    integer, dimension(0:), intent(in) :: counts
    integer, dimension(0:), intent(in), optional, target :: displs
    integer, dimension(:), pointer :: displs_
    !local variables
    logical :: variable
    integer :: dspl,jproc,chunk,nproc

    nproc=0
    variable=.false.
    nullify(displs_)
    nproc=size(counts)
    if (present(displs)) then
       if (size(displs) /= nproc) &
            call f_err_throw('Illegal sizes of displs ('//&
            trim(yaml_toa(nproc))//' /= '//&
            trim(yaml_toa(size(displs)))//')',&
            err_name='ERR_MPI_WRAPPERS')
       displs_=>displs
    else
       !in the variable case, determine the displacements in the case they are not given
       displs_=f_malloc_ptr(0.to.nproc-1,id='displs')
       dspl=0
       do jproc=0,nproc-1
          displs_(jproc)=dspl
          dspl=dspl+counts(jproc)
       end do
    end if

  end function get_displacements

  function get_variable_displacements(chunk,displs) result(variable)
    implicit none
    integer, intent(in) :: chunk
    integer, dimension(0:), intent(in) :: displs
    logical :: variable
    !local variables
    integer :: dspl,nproc,jproc

    variable=.false.
    !control if the provided displacements are all trivial
    dspl=0
    nproc=size(displs)
    do jproc=0,nproc-1
       if (dspl /= displs(jproc)) then
          variable=.true.
          exit
       end if
       dspl=dspl+chunk
    end do

  end function get_variable_displacements


!!$  subroutine verify_recvbuf_size(buf_size,counts,count,nproc)
!!$    use yaml_strings
!!$    implicit none
!!$    integer(f_long), intent(in) :: buf_size
!!$    integer, intent(in), optional :: count,nproc
!!$    integer, dimension(:), intent(in), optional :: counts
!!$    !local variables
!!$    integer(f_long) :: sizet
!!$
!!$    if (buf_size ==0) return
!!$    if(present(counts)) then
!!$       sizet=sum(int(counts,f_long))
!!$    else if (present(count)) then
!!$       sizet=nproc*int(count,f_long)
!!$    end if
!!$    !check the validity of the buffers
!!$    if (sizet < buf_size) call f_err_throw(&
!!$         'Size of buffer and of data are not consistent ('&
!!$         //trim(yaml_toa([sizet,buf_size]))//')',&
!!$         err_name='ERR_MPI_WRAPPERS')
!!$
!!$  end subroutine verify_recvbuf_size
!!$
!!$  function determine_displs(buf_size,count,counts,displs,comm) result(displs_)
!!$    implicit none
!!$    integer(f_long), intent(in) :: buf_size
!!$    integer, intent(in), optional :: count,comm
!!$    integer, dimension(0:), intent(in), optional :: counts
!!$    integer, dimension(0:), intent(in), optional, target :: displs
!!$    integer, dimension(:), pointer :: displs_
!!$    !local variables
!!$    integer :: ourselves
!!$
!!$    ourselves=mpisize(comm=comm)
!!$    displs_=>get_displacements(counts,displs)
!!$    call verify_recvbuf_size(buf_size,counts,count,ourselves)
!!$
!!$  end function determine_displs
!!$

!!$  subroutine counts_and_displs(buf_size,count,count_alt,counts,displs,displs_,null_displ,comm)
!!$    implicit none
!!$
!!$    me=mpirank(comm=comm)
!!$    ourselves=mpisize(comm=comm)
!!$    !calculate the receive displacements
!!$    if (present(count)) then
!!$       cnt=count
!!$       sizet=ourselves*int(cnt,f_long)
!!$    else if(present(count_alt)) then
!!$       cnt=count
!!$       sizet=ourselves*int(cnt,f_long)
!!$    else if(present(counts)) then
!!$       cnt=counts(me)
!!$       sizet=sum(int(counts,f_long))
!!$       if (size(counts) /= ourselves) &
!!$            call f_err_throw('Size of counts is not of length group size',&
!!$            err_name='ERR_MPI_WRAPPERS')
!!$    end if
!!$
!!$    displs_=>get_displacements(count,counts,displs)
!!$    null_displ=present(rdispls)
!!$    if (null_displ) null_displ=associated(displs_,target=displs)
!!$    call verify_recvbuf_size(buf_size,counts,count,ourselves)
!!$
!!$
!!$  end subroutine counts_and_displs
!!$
!!$  subroutine sendrecv_sizes(sendbuf_size,recvbuf_size,sendcount,sendcounts,sdispls,&
!!$       recvcount,recvcounts,rdispls,sdispls_,rdispls_,null_sdispl,null_rdispl,comm)
!!$    implicit none
!!$    integer, intent(in) :: sendbuf_size
!!$    integer, intent(in) :: recvbuf_size
!!$    integer, intent(in), optional :: comm
!!$    integer, dimension(0:), intent(in), optional :: recvcounts
!!$    !local variables
!!$    logical :: allv
!!$    integer :: me,ourselves
!!$    integer, dimension(:), pointer :: displs_
!!$
!!$    call counts_and_displs(recvbuf_size,recvcount,sendcount,recvcounts,rdispls,rdispls_,null_rdispl,comm)
!!$    call counts_and_displs(sendbuf_size,sendcount,recvcount,sendcounts,sdispls,sdispls_,null_sdispl,comm)
!!$
!!$  end subroutine sendrecv_sizes

!!$end module f_control
