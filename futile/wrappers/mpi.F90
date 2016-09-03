!> @file
!! Wrapper for the MPI call (this file is preprocessed.)
!! Use error handling
!! @author
!!    Copyright (C) 2012-2015 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


#if defined HAVE_CONFIG_H
#include <config.inc>
#endif


!> Module defining the routines which wrap the MPI calls
module wrapper_MPI
  use time_profiling, only: TIMING_UNINITIALIZED
  use yaml_strings, only: operator(//)
  use f_precisions
  use f_refcnts
  use dictionaries, only: f_err_throw
  implicit none

  ! MPI handling
#ifdef HAVE_MPI2
  logical, parameter :: have_mpi2 = .true.  !< Flag to use in the code to switch between MPI1 and MPI2
#else
  integer :: MPI_IN_PLACE = 0               !< Fake MPI_IN_PLACE variable to allow compilation in sumrho.
  logical, parameter :: have_mpi2 = .false. !< Flag to use in the code to switch between MPI1 and MPI2
#endif

  include 'mpif.h'      !< MPI definitions and datatypes

  logical :: mpi_thread_funneled_is_supported=.false. !< Control the OMP_NESTED based overlap, checked by bigdft_mpi_init below

  !timing categories for MPI wrapper
  integer, parameter :: smallsize=5 !< limit for a communication with small size
  character(len=*), parameter, public :: tgrp_mpi_name='Communications'
  !timing categories
  integer, public, save :: TCAT_ALLRED_SMALL = TIMING_UNINITIALIZED
  integer, public, save :: TCAT_ALLRED_LARGE = TIMING_UNINITIALIZED
  integer, public, save :: TCAT_ALLGATHERV   = TIMING_UNINITIALIZED
  integer, public, save :: TCAT_ALLGATHER    = TIMING_UNINITIALIZED
  integer, public, save :: TCAT_GATHER       = TIMING_UNINITIALIZED
  integer, public, save :: TCAT_SCATTER      = TIMING_UNINITIALIZED
  integer, public, save :: TCAT_FENCE        = TIMING_UNINITIALIZED
  integer, public, save :: TCAT_SEND         = TIMING_UNINITIALIZED
  integer, public, save :: TCAT_RECV         = TIMING_UNINITIALIZED
  integer, public, save :: TCAT_WAIT         = TIMING_UNINITIALIZED

  !> Error codes
  integer, public, save :: ERR_MPI_WRAPPERS

  !> Interface for MPITYPE routine
  interface mpitype
     module procedure mpitype_i,mpitype_d,mpitype_r,mpitype_l,mpitype_c,mpitype_li
     module procedure mpitype_i1,mpitype_i2,mpitype_i3
     module procedure mpitype_l3
     module procedure mpitype_r1,mpitype_r2,mpitype_r3,mpitype_r4
     module procedure mpitype_d1,mpitype_d2,mpitype_d3,mpitype_d4,mpitype_d5
     module procedure mpitype_c1
     module procedure mpitype_li1,mpitype_li2,mpitype_li3
  end interface mpitype

  interface mpimaxdiff
     module procedure mpimaxdiff_i0,mpimaxdiff_li0,mpimaxdiff_d0
     module procedure mpimaxdiff_i1,mpimaxdiff_li1,mpimaxdiff_d1
     module procedure mpimaxdiff_i2,mpimaxdiff_d2
  end interface mpimaxdiff

  !> Interface for MPI_ALLREDUCE routine, to be updated little by little
  interface mpiallred
     module procedure mpiallred_int,mpiallred_real, &
          & mpiallred_double,&!,mpiallred_double_1,mpiallred_double_2,&
          & mpiallred_log
     module procedure mpiallred_long
     module procedure mpiallred_r1,mpiallred_r2,mpiallred_r3,mpiallred_r4
     module procedure mpiallred_d1,mpiallred_d2,mpiallred_d3,mpiallred_d4,mpiallred_d5
     module procedure mpiallred_i1,mpiallred_i2,mpiallred_i3
     module procedure mpiallred_l3
  end interface mpiallred

  interface mpigather
     module procedure mpigather_d0d2,mpigather_d1d1,mpigather_d1d2,mpigather_d2,mpigather_d2d1
     module procedure mpigather_i0i2,mpigather_i1,mpigather_i1i2,mpigather_i2
     module procedure mpigather_li0li2,mpigather_li1,mpigather_li1li2,mpigather_li2
     module procedure mpigather_c1i2
     module procedure mpigather_c1li2
  end interface mpigather

  interface mpibcast
     module procedure mpibcast_i0,mpibcast_li0,mpibcast_d0,mpibcast_c0
     module procedure mpibcast_c1,mpibcast_d1,mpibcast_d2,mpibcast_i1,mpibcast_i2
  end interface mpibcast

  interface mpiscatter
      module procedure mpiscatter_i1i1
  end interface mpiscatter

  interface mpiscatterv
     module procedure mpiscatterv_d0
     module procedure mpiscatterv_d2d3,mpiscatterv_d3d2
  end interface mpiscatterv

  interface mpi_get_to_allgatherv
     module procedure mpi_get_to_allgatherv_double
  end interface mpi_get_to_allgatherv

  interface mpiget
    module procedure mpiget_d0
  end interface mpiget

  interface mpisend
     module procedure mpisend_d0, mpisend_gpu
  end interface mpisend

  interface mpirecv
     module procedure mpirecv_d0,mpirecv_gpu
  end interface mpirecv

  interface mpiput
     module procedure mpiput_d0
  end interface mpiput

  interface mpiaccumulate
     module procedure mpiaccumulate_d0
  end interface mpiaccumulate


  interface mpitypesize
    module procedure mpitypesize_d0, mpitypesize_d1, mpitypesize_i0, mpitypesize_l0
  end interface mpitypesize

  interface mpiwindow
    module procedure mpiwindow_d0, mpiwindow_i0, mpiwindow_l0
  end interface mpiwindow

  !> Interface for MPI_ALLGATHERV routine
  interface mpiallgather
     module procedure mpiallgatherv_d0,mpiallgatherv_d1,mpiallgatherv_d2d3
  end interface mpiallgather

  interface mpiiallred
      module procedure mpiiallred_double
  end interface mpiiallred

  interface mpialltoallv
      module procedure mpialltoallv_int, mpialltoallv_long, mpialltoallv_double
  end interface mpialltoallv

  interface mpiialltoallv
      module procedure mpiialltoallv_double
  end interface mpiialltoallv

!!$  interface mpiaccumulate
!!$      module procedure mpiaccumulate_double
!!$  end interface mpiaccumulate

  !> Global MPI communicator which contains all information related to the MPI process
  type, public :: mpi_environment
     !>reference counter of the communicator.
     !!used to understand whether the communicator has to be destroyed
     type(f_reference_counter) :: refcnt
     integer :: mpi_comm !< MPI communicator
     integer :: iproc    !< Process Id
                         !! @ingroup RESERVED
     integer :: nproc    !< Number of MPI processes (in the given communicator)
                         !! @ingroup RESERVED
     integer :: igroup   !< MPI Group Id
     integer :: ngroup   !< Number of MPI groups
  end type mpi_environment

  public :: mpi_environment_null
  public :: release_mpi_environment
  public :: mpi_environment_set
  public :: mpi_environment_set1 !to be removed

  !>fake type to enhance documentation
  type, private :: doc
     !>number of entries in buffer (integer). Useful for buffer passed by reference
     integer :: count
     !> rank of mpitask executing the operation (default value is root=0)
     integer :: root
     !> communicator of the communication
     integer :: comm
  end type doc

  private :: operator(//),f_err_throw

contains

  pure subroutine nullify_mpi_environment(mpi)
    implicit none
    type(mpi_environment), intent(out) :: mpi
    call nullify_f_ref(mpi%refcnt)
    mpi%mpi_comm=MPI_COMM_NULL !better to put an invalid comm?
    mpi%igroup=-1
    mpi%ngroup=-1
    mpi%iproc=-1
    mpi%nproc=-1
  end subroutine nullify_mpi_environment

  pure function mpi_environment_null() result(mpi)
    implicit none
    type(mpi_environment) :: mpi
    call nullify_mpi_environment(mpi)
  end function mpi_environment_null

  subroutine release_mpi_environment(mpi_env)
    use yaml_strings!, only: yaml_toa,operator(//),f_string
!    use dictionaries, only: f_err_throw
    implicit none
    type(mpi_environment), intent(inout) :: mpi_env
    !local variables
    integer :: ierr,count

    !first check if we are in a nullified status. If so, do nothing
    if (mpi_env%mpi_comm /= MPI_COMM_NULL .and. f_associated(mpi_env%refcnt)) then
       call f_unref(mpi_env%refcnt,count=count)
       !if the communicator is still active, just destroy the structure
       if (count==0) then
          call f_ref_free(mpi_env%refcnt)
          !free the MPI communicator if it is not WORLD
          !also there is no need to
          if (mpi_env%mpi_comm /= MPI_COMM_WORLD) then
             call MPI_COMM_FREE(mpi_env%mpi_comm,ierr)
             if (ierr /=0) then
                call f_err_throw('Problem in MPI_COMM_FREE, ierr='//ierr,&
                     err_name='BIGDFT_MPI_ERROR')
                return
             end if
          end if
       end if
    end if
    !in any case nullify the status of the mpi_env
    mpi_env=mpi_environment_null()
  end subroutine release_mpi_environment

  subroutine deepcopy_mpi_environment(dest,src)
    implicit none
    ! Calling arguments
    type(mpi_environment),intent(in) :: src
    type(mpi_environment),intent(out) :: dest
    ! Local variables
    integer :: ierr

    dest=mpi_environment_null()

    if (src%mpi_comm/=MPI_COMM_NULL) then
       call mpi_comm_dup(src%mpi_comm, dest%mpi_comm, ierr)
       if (ierr /=0) then
          call f_err_throw('Problem in MPI_COMM_DUP, ierr='//ierr,&
               err_name='BIGDFT_MPI_ERROR')
          return
       end if
       dest%nproc=mpisize(dest%mpi_comm)
       dest%iproc=mpirank(dest%mpi_comm)
       !call mpi_comm_size(dest%mpi_comm, dest%nproc, ierr)
       !LG: here there was a BIGGG mistake! (nproc instead of iproc)
       !call mpi_comm_rank(dest%mpi_comm, dest%nproc, ierr)
       dest%igroup = src%igroup
       dest%ngroup = src%ngroup
       !a new reference counter has to be activated, if the
       !source has been
       dest%refcnt=f_ref_new('mpi_copied')
    end if

  end subroutine deepcopy_mpi_environment

  !>shallow copy of the mpi_environment.
  !! it has no effect if the src has a null communicator
  subroutine copy_mpi_environment(dest,src)
    implicit none
    ! Calling arguments
    type(mpi_environment),intent(in) :: src
    type(mpi_environment),intent(out) :: dest

    dest=mpi_environment_null()

    if (src%mpi_comm/=MPI_COMM_NULL) then
       !if meaningful copy the reference counter
       if (f_associated(src%refcnt)) &
            call f_ref_associate(src=src%refcnt,dest=dest%refcnt)
       dest%nproc=src%nproc
       dest%iproc=src%iproc
       dest%igroup = src%igroup
       dest%ngroup = src%ngroup
       dest%mpi_comm=src%mpi_comm
    end if

  end subroutine copy_mpi_environment

  !> Set the MPI environment (i.e. taskgroup or MPI communicator)
  subroutine mpi_environment_set(mpi_env,iproc,nproc,mpi_comm,groupsize)
    use dynamic_memory
    use yaml_output
    implicit none
    integer, intent(in) :: iproc                   !< Proc id
    integer, intent(in) :: nproc                   !< Total number of MPI processes
    integer, intent(in) :: mpi_comm                !< Global MPI_communicator
    integer, intent(in) :: groupsize               !< Number of MPI processes by (task)group
                                                   !! if 0 one taskgroup (MPI_COMM_WORLD)
    type(mpi_environment), intent(out) :: mpi_env  !< MPI environment (out)
    !local variables
    integer :: j,base_grp
    integer, dimension(:), allocatable :: group_list

    call f_routine(id='mpi_environment_set')
    mpi_env=mpi_environment_null()
    mpi_env%mpi_comm=mpi_comm
    !these assignments would depend on the strategy adopted for
    !taskgroup creations
    mpi_env%igroup=iproc/groupsize
    mpi_env%ngroup=nproc/groupsize
    mpi_env%iproc=mod(iproc,groupsize)
    mpi_env%nproc=groupsize
    if (groupsize /= nproc) then
       !define the strategy for the taskgroups
       group_list=f_malloc(groupsize,id='group_list')
       !iproc in the same group are close to each other
       do j=0,groupsize-1
          group_list(j+1)=mpi_env%igroup*groupsize+j
       enddo
       base_grp=mpigroup(mpi_comm)
       call mpi_env_create_group(iproc/groupsize,nproc/groupsize,mpi_comm,&
            base_grp,groupsize,group_list,mpi_env)
       call mpigroup_free(base_grp)
       !call create_group_comm(mpi_comm,groupsize,group_list,mpi_env%mpi_comm)
       if (iproc == 0) then
          call yaml_map('Total No. of Taskgroups created',nproc/mpi_env%nproc)
       end if
       call f_free(group_list)

    end if
    call f_release_routine()
  end subroutine mpi_environment_set


!!! PSolver n1-n2 plane mpi partitioning !!!
  !> This is exactly like mpi_environment_set but it always creates groups
  !! the routine above should be modified accordingly
!!$  subroutine mpi_environment_set2(mpi_env,iproc,nproc,mpi_comm,groupsize)
!!$    use yaml_output
!!$    implicit none
!!$    integer, intent(in) :: iproc,nproc,mpi_comm,groupsize
!!$    type(mpi_environment), intent(out) :: mpi_env
!!$    !local variables
!!$    integer :: j
!!$    integer, dimension(:), allocatable :: group_list
!!$
!!$    call f_routine(id='mpi_environment_set2')
!!$    mpi_env=mpi_environment_null()
!!$
!!$    mpi_env%mpi_comm=mpi_comm
!!$
!!$    mpi_env%igroup=iproc/groupsize
!!$    mpi_env%ngroup=nproc/groupsize
!!$    mpi_env%iproc=mod(iproc,groupsize)
!!$    mpi_env%nproc=groupsize
!!$
!!$    !define the strategy for the taskgroups
!!$    group_list=f_malloc(groupsize,id='group_list')
!!$    !iproc in the same group are close to each other
!!$    do j=0,groupsize-1
!!$       group_list(j+1)=mpi_env%igroup*groupsize+j
!!$    enddo
!!$
!!$    call create_group_comm(mpi_comm,nproc,mpi_env%igroup,mpi_env%nproc,group_list,mpi_env%mpi_comm)
!!$!    if (iproc == 0) then
!!$!       call yaml_map('Total No. of Taskgroups created',nproc/mpi_env%nproc)
!!$!    end if
!!$    call f_free(group_list)
!!$    call f_release_routine()
!!$  end subroutine mpi_environment_set2


  !> This is a different procedure to assign the iproc according to the groups.
  subroutine mpi_environment_set1(mpi_env,iproc,mpi_comm,groupsize,ngroup)
    use yaml_output
    use dynamic_memory
    implicit none
    integer, intent(in) :: iproc,mpi_comm,groupsize,ngroup
    type(mpi_environment), intent(out) :: mpi_env
    !local variables
    integer :: j
    integer, dimension(:), allocatable :: group_list

    call f_routine(id='mpi_environment_set1')

    mpi_env=mpi_environment_null()

    mpi_env%igroup=-1

    mpi_env%ngroup=ngroup
    if (iproc < groupsize*ngroup) mpi_env%igroup=mod(iproc,ngroup)
    mpi_env%iproc=iproc/ngroup
    mpi_env%nproc=groupsize
    mpi_env%mpi_comm=mpi_comm

    !define the strategy for the taskgroups
    group_list=f_malloc(groupsize,id='group_list')
    !round-robin strategy
    if (mpi_env%igroup >0) then
       do j=0,groupsize-1
          group_list(j+1)=mpi_env%igroup+j*mpi_env%ngroup
       enddo
    else
       !these processes have MPI_COMM_NULL
       group_list=-1
       mpi_env%mpi_comm=MPI_COMM_NULL
    end if

    !call create_group_comm1(mpi_comm,nproc,mpi_env%igroup,ngroup,mpi_env%nproc,mpi_env%mpi_comm)
    call create_group_comm(mpi_comm,mpi_env%nproc,group_list,mpi_env%mpi_comm)
    !    if (iproc == 0) then
    !       call yaml_map('Total No. of Taskgroups created',ngroup)
    !    end if
    call f_free(group_list)
    mpi_env%refcnt=f_ref_new('MPI_env1')
    call f_release_routine()
  end subroutine mpi_environment_set1

  !> create a mpi_environment from a group list in a base group
  subroutine mpi_env_create_group(igrp,ngrp,base_comm,base_grp,group_size,group_list,&
       mpi_env)
    implicit none
    integer, intent(in) :: igrp,ngrp !<id and no of the groups
    integer, intent(in) :: group_size,base_comm,base_grp
    integer, dimension(group_size), intent(in) :: group_list
    type(mpi_environment), intent(out) :: mpi_env
    !local variables
    integer :: grp,ierr

    !create the groups with the list
    grp=mpigroupincl(base_grp,group_size,group_list)
    !create the communicator (the communicator can be also null)
    call MPI_COMM_CREATE(base_comm,grp,mpi_env%mpi_comm,ierr)
    if (ierr /= 0) call f_err_throw('Problem in comm_create, ierr:'//&
         ierr,err_name='BIGDFT_MPI_ERROR')
    !free temporary group
    call mpigroup_free(grp)
    !then fill iproc and nproc
    if (mpi_env%mpi_comm /= MPI_COMM_NULL) then
       mpi_env%iproc=mpirank(mpi_env%mpi_comm)
       mpi_env%nproc=mpisize(mpi_env%mpi_comm) !this should be group_size
       mpi_env%igroup=igrp
       mpi_env%ngroup=ngrp
       mpi_env%refcnt=f_ref_new('MPI_env_from_grp')
    end if
  end subroutine mpi_env_create_group

  !> function that returns the handle of the group of different
  !! processes that will belong to a rma communication
  function mpigroupincl(base_grp,group_size,group_list) result(grp)
    implicit none
    integer, intent(in) :: base_grp,group_size
    !>ranks of the groups in the rma access pattern
    integer, dimension(group_size), intent(in) :: group_list
    integer :: grp
    !local variables
    integer :: ierr
    !create the groups with the list
    call MPI_GROUP_INCL(base_grp,group_size,group_list,grp,ierr)
    if (ierr /= 0) call f_err_throw('Problem in group inclusion, ierr:'//&
         ierr,err_name='BIGDFT_MPI_ERROR')
  end function mpigroupincl

  function mpigroup(comm)
    implicit none
    integer, intent(in), optional :: comm
    integer :: mpigroup
    !local variables
    integer :: ierr,mpi_comm

    if (present(comm)) then
       mpi_comm=comm
    else
       mpi_comm=MPI_COMM_WORLD
    end if
    call MPI_COMM_GROUP(mpi_comm,mpigroup,ierr)
    if (ierr /= 0) then
       mpigroup=mpigroup_null()
       call f_err_throw('Problem in group identification, ierr:'//&
            ierr,err_name='BIGDFT_MPI_ERROR')
    end if
  end function mpigroup

  pure function mpigroup_null() result(grp)
    implicit none
    integer :: grp
    grp=MPI_GROUP_NULL
  end function mpigroup_null

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

  pure function mpirequest_null() result(request)
    implicit none
    integer :: request
    request=MPI_REQUEST_NULL
  end function mpirequest_null


  subroutine mpigroup_free(grp)
    implicit none
    integer, intent(inout) :: grp
    !local variables
    integer :: ierr
    ierr=0
    if (grp /= MPI_GROUP_NULL) call MPI_GROUP_FREE(grp,ierr)
    if (ierr /= 0) then
       call f_err_throw('Problem in group free, ierr:'//&
            ierr,err_name='BIGDFT_MPI_ERROR')
    end if
    grp=MPI_GROUP_NULL
  end subroutine mpigroup_free

  !> Create communicators associated to the groups of size group_size
  subroutine create_group_comm(base_comm,group_size,group_list,group_comm)
    implicit none
    integer, intent(in) :: base_comm,group_size
    integer, dimension(group_size), intent(in) :: group_list !< list of id of the group identified by group_id in units of base_comm
    integer, intent(out) :: group_comm
    !local variables
    integer :: grp,ierr,base_grp

    !take the base group
    call MPI_COMM_GROUP(base_comm,base_grp,ierr)
    if (ierr /= 0) then
       call check_ierr(ierr,'group creation')
       return
    end if
    !create the groups with the list
    call MPI_GROUP_INCL(base_grp,group_size,group_list,grp,ierr)
    if (ierr /= 0) then
       call check_ierr(ierr,'group inclusion')
       return
    end if
    !create the communicator (the communicator can be also null)
    call MPI_COMM_CREATE(base_comm,grp,group_comm,ierr)
    if (ierr /= 0) then
       call check_ierr(ierr,'communicator creator')
       return
    end if
    !free temporary group
    call MPI_GROUP_FREE(grp,ierr)
    if (ierr /= 0) then
       call check_ierr(ierr,'new_group free')
       return
    end if

    !free base group
    call MPI_GROUP_FREE(base_grp,ierr)
    if (ierr /= 0) then
       call check_ierr(ierr,'base_group free')
       return
    end if


  contains

    subroutine check_ierr(ierr,message)
      use yaml_strings, only: yaml_toa
      use dictionaries, only: f_err_throw
      implicit none
      integer, intent(in) :: ierr
      character(len=*), intent(in) :: message
      if (ierr /= 0) then
         call f_err_throw('Problem in '//trim(message)//&
              ', ierr:'//yaml_toa(ierr),err_name='BIGDFT_MPI_ERROR')
      end if
    end subroutine check_ierr

  end subroutine create_group_comm


!!! PSolver n1-n2 plane mpi partitioning !!!
  !> This routine is like create_group_comm with a different group_list
  subroutine create_group_comm1(base_comm,group_id,ngroup,group_size,group_comm)
    use dynamic_memory
    use yaml_output
    implicit none
    integer, intent(in) :: base_comm,group_size,group_id,ngroup
    integer, intent(out) :: group_comm
    !local variables
    character(len=*), parameter :: subname='create_group_comm'
    integer :: grp,ierr,i,j,base_grp,temp_comm!,i_stat,i_all
    integer, dimension(:), allocatable :: group_list

    ! allocate(group_list(group_size+ndebug),stat=i_stat)
    group_list = f_malloc(group_size,id='group_list')

    !take the base group
    call MPI_COMM_GROUP(base_comm,base_grp,ierr)
    if (ierr /=0) then
       call yaml_warning('Problem in group creation, ierr='//ierr)
       call MPI_ABORT(base_comm,1,ierr)
    end if
    do i=0,ngroup-1
       !define the new groups and thread_id
       do j=0,group_size-1
          group_list(j+1)=i+j*ngroup
       enddo
       call MPI_GROUP_INCL(base_grp,group_size,group_list,grp,ierr)
       if (ierr /=0) then
          call yaml_warning('Problem in group inclusion, ierr='//ierr)
          call MPI_ABORT(base_comm,1,ierr)
       end if
       call MPI_COMM_CREATE(base_comm,grp,temp_comm,ierr)
       if (ierr /=0) then
          call yaml_warning('Problem in communicator creator, ierr='//ierr)
          call MPI_ABORT(base_comm,1,ierr)
       end if
       !print *,'i,group_id,temp_comm',i,group_id,temp_comm
       if (i.eq. group_id) group_comm=temp_comm
    enddo

    !i_all=-product(shape(group_list ))*kind(group_list )
    ! deallocate(group_list,stat=i_stat)
    call f_free(group_list)
  end subroutine create_group_comm1


  !> Create a communicator between proc of same rank between the taskgroups.
  subroutine create_rank_comm(group_comm, rank_comm)
    use dynamic_memory
    use yaml_output
    implicit none
    integer, intent(in) :: group_comm
    integer, intent(out) :: rank_comm
    !local variables
    character(len=*), parameter :: subname='create_group_master'
    integer :: iproc_group, nproc, nproc_group, ngroups
    integer :: ierr, i, j
    integer, dimension(:), allocatable :: lrank, ids

    call MPI_COMM_RANK(group_comm, iproc_group, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)
    call MPI_COMM_SIZE(group_comm, nproc_group, ierr)
    ngroups = nproc / nproc_group

    ! Put in lrank the group rank of each process, indexed by global iproc.
    !   allocate(lrank(nproc+ndebug), stat = i_stat)
    lrank = f_malloc(nproc,id='lrank')
    call mpi_allgather(iproc_group, 1, MPI_INTEGER, lrank, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

    ! Put in ids, the global iproc of each process that share the same group iproc.
    !   allocate(ids(ngroups+ndebug), stat = i_stat)
    ids = f_malloc(ngroups,id='ids')
    j = 1
    do i = 1, nproc
       if (lrank(i) == iproc_group) then
          ids(j) = i - 1
          j = j + 1
       end if
    end do
    !  i_all=-product(shape(lrank ))*kind(lrank )
    !   deallocate(lrank,stat=i_stat)
    call f_free(lrank)

!!$    call mpi_comm_rank(MPI_COMM_WORLD, iproc_group, ierr)
!!$    write(*,*) iproc_group, "->", ids

    ! Create a new comminucator for the list of ids.
    call create_group_comm(MPI_COMM_WORLD, ngroups, ids, rank_comm)
    !  i_all=-product(shape(ids ))*kind(ids )
    !   deallocate(ids,stat=i_stat)
    call f_free(ids)
  END SUBROUTINE create_rank_comm


  subroutine wmpi_init_thread(ierr)
    use dictionaries, only: f_err_throw
    implicit none
    integer, intent(out) :: ierr
#ifdef HAVE_MPI_INIT_THREAD
    integer :: provided
    call MPI_INIT_THREAD(MPI_THREAD_FUNNELED,provided,ierr)
    if (ierr /= MPI_SUCCESS) then
       write(*,*)'BigDFT_mpi_INIT: Error in MPI_INIT_THREAD',ierr
    else if (provided < MPI_THREAD_FUNNELED) then
       !write(*,*)'WARNING: MPI_THREAD_FUNNELED not supported!',provided,ierr
       !call MPI_INIT(ierr)
    else
       mpi_thread_funneled_is_supported=.true.
    endif
#else
    call MPI_INIT(ierr)
    if (ierr /= MPI_SUCCESS) then
       call f_err_throw('An error in calling to MPI_INIT (THREAD) occured',&
            err_id=ERR_MPI_WRAPPERS)
    end if
#endif
  end subroutine wmpi_init_thread

  function mpihostname()
    implicit none
    character(len=MPI_MAX_PROCESSOR_NAME) :: mpihostname
    !local variables
    integer :: ierr,namelen,ipos,i

    call MPI_GET_PROCESSOR_NAME(mpihostname,namelen,ierr)
    if (ierr /= MPI_SUCCESS) then
       call f_err_throw('An error in calling to MPI_GET_PROCESSOR_NAME occured',&
            err_id=ERR_MPI_WRAPPERS)
    end if

    !clean the hostname such as to include only the last word
    !this solves a problem in ibm machines
    ipos=index(mpihostname,' ',back=.true.)
    if (ipos > 0) then
       do i=1,len(mpihostname)
          if (i+ipos+1 <= len(mpihostname)) then
             mpihostname(i:i)=mpihostname(i+ipos+1:i+ipos+1)
          else
             mpihostname(i:i)=' '
          end if
       end do
    end if

  end function mpihostname

  !>initialization of the mpi library
  subroutine mpiinit(inithread)
    use dictionaries, only: f_err_throw
    implicit none
    !>if present, set the initialization to the
    !!mpi_init_thread case (mpi_thread_funneled is supported)
    !! default is false, traditional mpi_init
    logical, intent(in), optional :: inithread
    !local variables
    logical :: thd
    integer :: ierr
    external :: MPI_INIT

    thd=.false.
    if (present(inithread)) thd=inithread

    if (thd) then
       call wmpi_init_thread(ierr)
    else
       call MPI_INIT(ierr)
    end if
    if (ierr /=0) call f_err_throw('An error in calling to MPI_INIT (THREAD) occured',&
         err_id=ERR_MPI_WRAPPERS)

  end subroutine mpiinit

  pure function mpiworld()
    implicit none
    integer :: mpiworld
    mpiworld=MPI_COMM_WORLD
  end function mpiworld

  !> Finalization of the mpi
  subroutine mpifinalize()
    use dictionaries, only: f_err_throw
    implicit none
    !local variables
    integer :: ierr

    call MPI_FINALIZE(ierr)
    if (ierr /= MPI_SUCCESS) then
       call f_err_throw('An error in calling to MPI_FINALIZE occured',&
            err_id=ERR_MPI_WRAPPERS)
    end if
  end subroutine mpifinalize


  !> Initialize timings and also mpi errors
  subroutine mpi_initialize_timing_categories()
    use time_profiling, only: f_timing_category_group,f_timing_category
    use dictionaries, only: f_err_throw,f_err_define
    use yaml_strings, only: yaml_toa
    implicit none

    call f_timing_category_group(tgrp_mpi_name,&
         'Operations between MPI tasks')

    call f_timing_category('Allreduce, Small Size',tgrp_mpi_name,&
         'Allreduce operations for less than'//&
         trim(yaml_toa(smallsize))//' elements',&
         TCAT_ALLRED_SMALL)
    call f_timing_category('Allreduce, Large Size',tgrp_mpi_name,&
         'Allreduce operations for more than'//&
         trim(yaml_toa(smallsize))//' elements',&
         TCAT_ALLRED_LARGE)
    call f_timing_category('Allgatherv',tgrp_mpi_name,&
         'Variable allgather operations',&
         TCAT_ALLGATHERV)
    call f_timing_category('Allgather',tgrp_mpi_name,&
         'Allgather operations',&
         TCAT_ALLGATHER)
    call f_timing_category('Gather',tgrp_mpi_name,&
         'Gather operations, in general moderate size arrays',&
         TCAT_GATHER)
    call f_timing_category('Scatter',tgrp_mpi_name,&
         'Scatter operations, in general moderate size arrays',&
         TCAT_SCATTER)
    call f_timing_category('Fence',tgrp_mpi_name,&
         'Fence, waiting for a RMA operation to end',&
         TCAT_FENCE)
    call f_timing_category('Send',tgrp_mpi_name,&
         'Time spent in MPI_Send or MPI_Isend',&
         TCAT_SEND)
    call f_timing_category('Recv',tgrp_mpi_name,&
         'Time spent in MPI_Recv or MPI_Irecv',&
         TCAT_RECV)
    call f_timing_category('Wait',tgrp_mpi_name,&
         'Time spent in MPI_Wait or MPI_Waitall',&
         TCAT_WAIT)
    call f_err_define(err_name='ERR_MPI_WRAPPERS',err_msg='Error of MPI library',&
         err_id=ERR_MPI_WRAPPERS,&
         err_action='Some MPI library returned an error code, inspect runtime behaviour')

  end subroutine mpi_initialize_timing_categories

  pure function mpitype_i(data) result(mt)
    implicit none
    integer(f_integer), intent(in) :: data
    integer :: mt
    mt=MPI_INTEGER
  end function mpitype_i
  pure function mpitype_i1(data) result(mt)
    implicit none
    integer(f_integer), dimension(:), intent(in) :: data
    integer :: mt
    mt=MPI_INTEGER
  end function mpitype_i1
  pure function mpitype_i2(data) result(mt)
    implicit none
    integer(f_integer), dimension(:,:), intent(in) :: data
    integer :: mt
    mt=MPI_INTEGER
  end function mpitype_i2
  pure function mpitype_i3(data) result(mt)
    implicit none
    integer(f_integer), dimension(:,:,:), intent(in) :: data
    integer :: mt
    mt=MPI_INTEGER
  end function mpitype_i3

  pure function mpitype_l3(data) result(mt)
    implicit none
    logical, dimension(:,:,:), intent(in) :: data
    integer :: mt
    mt=MPI_LOGICAL
  end function mpitype_l3

  pure function mpitype_li(data) result(mt)
    implicit none
    integer(f_long), intent(in) :: data
    integer :: mt
    mt=MPI_INTEGER8
  end function mpitype_li
  pure function mpitype_li1(data) result(mt)
    implicit none
    integer(f_long), dimension(:), intent(in) :: data
    integer :: mt
    mt=MPI_INTEGER8
  end function mpitype_li1
  pure function mpitype_li2(data) result(mt)
    implicit none
    integer(f_long), dimension(:,:), intent(in) :: data
    integer :: mt
    mt=MPI_INTEGER8
  end function mpitype_li2
  pure function mpitype_li3(data) result(mt)
    implicit none
    integer(f_long), dimension(:,:,:), intent(in) :: data
    integer :: mt
    mt=MPI_INTEGER8
  end function mpitype_li3


  pure function mpitype_r(data) result(mt)
    implicit none
    real, intent(in) :: data
    integer :: mt
    mt=MPI_REAL
  end function mpitype_r
  pure function mpitype_d(data) result(mt)
    implicit none
    double precision, intent(in) :: data
    integer :: mt
    mt=MPI_DOUBLE_PRECISION
  end function mpitype_d
  pure function mpitype_d1(data) result(mt)
    implicit none
    double precision, dimension(:), intent(in) :: data
    integer :: mt
    mt=MPI_DOUBLE_PRECISION
  end function mpitype_d1
  pure function mpitype_d2(data) result(mt)
    implicit none
    double precision, dimension(:,:), intent(in) :: data
    integer :: mt
    mt=MPI_DOUBLE_PRECISION
  end function mpitype_d2
  pure function mpitype_d3(data) result(mt)
    implicit none
    double precision, dimension(:,:,:), intent(in) :: data
    integer :: mt
    mt=MPI_DOUBLE_PRECISION
  end function mpitype_d3
  pure function mpitype_d4(data) result(mt)
    implicit none
    double precision, dimension(:,:,:,:), intent(in) :: data
    integer :: mt
    mt=MPI_DOUBLE_PRECISION
  end function mpitype_d4
  pure function mpitype_d5(data) result(mt)
    implicit none
    double precision, dimension(:,:,:,:,:), intent(in) :: data
    integer :: mt
    mt=MPI_DOUBLE_PRECISION
  end function mpitype_d5

  pure function mpitype_r1(data) result(mt)
    implicit none
    real, dimension(:), intent(in) :: data
    integer :: mt
    mt=MPI_REAL
  end function mpitype_r1
  pure function mpitype_r2(data) result(mt)
    implicit none
    real, dimension(:,:), intent(in) :: data
    integer :: mt
    mt=MPI_REAL
  end function mpitype_r2
  pure function mpitype_r3(data) result(mt)
    implicit none
    real, dimension(:,:,:), intent(in) :: data
    integer :: mt
    mt=MPI_REAL
  end function mpitype_r3
  pure function mpitype_r4(data) result(mt)
    implicit none
    real, dimension(:,:,:,:), intent(in) :: data
    integer :: mt
    mt=MPI_REAL
  end function mpitype_r4

  pure function mpitype_l(data) result(mt)
    implicit none
    logical, intent(in) :: data
    integer :: mt
    mt=MPI_LOGICAL
  end function mpitype_l
  pure function mpitype_c(data) result(mt)
    implicit none
    character, intent(in) :: data
    integer :: mt
    mt=MPI_CHARACTER
  end function mpitype_c
  pure function mpitype_c1(data) result(mt)
    implicit none
    character, dimension(:), intent(in) :: data
    integer :: mt
    mt=MPI_CHARACTER
  end function mpitype_c1

  !> Function giving the mpi rank id for a given communicator
  function mpirank(comm)
    use dictionaries, only: f_err_throw
    implicit none
    integer, intent(in), optional :: comm
    integer :: mpirank
    !local variables
    integer :: iproc,ierr,mpi_comm

    if (mpiinitialized()) then
       if (present(comm)) then
          mpi_comm=comm
       else
          mpi_comm=MPI_COMM_WORLD
       end if

       call MPI_COMM_RANK(mpi_comm, iproc, ierr)
       if (ierr /=0) then
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
    integer :: nproc,ierr,mpi_comm

    if (mpiinitialized()) then
       if (present(comm)) then
          mpi_comm=comm
       else
          mpi_comm=MPI_COMM_WORLD
       end if

       !verify the size of the receive buffer
       call MPI_COMM_SIZE(mpi_comm,nproc,ierr)
       if (ierr /=0) then
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

  !> Give the number of MPI processes per node (nproc_node) and before iproc (iproc_node)
  subroutine mpinoderanks(iproc,nproc,mpi_comm,iproc_node,nproc_node)
    use dynamic_memory
    implicit none
    integer, intent(in) :: iproc,nproc,mpi_comm
    integer, intent(out) :: iproc_node, nproc_node
    !local variables
    character(len=*), parameter :: subname='processor_id_per_node'
    integer :: ierr,jproc
    character(len=MPI_MAX_PROCESSOR_NAME) :: nodename_local
    character(len=MPI_MAX_PROCESSOR_NAME), dimension(:), allocatable :: nodename

    call f_routine(id=subname)

    if (nproc == 1) then
       iproc_node=0
       nproc_node=1
    else
       nodename=f_malloc0_str(MPI_MAX_PROCESSOR_NAME,0 .to. nproc-1,id='nodename')

       nodename_local=mpihostname()

       !gather the result between all the processes
       call MPI_ALLGATHER(nodename_local,MPI_MAX_PROCESSOR_NAME,MPI_CHARACTER,&
            nodename(0),MPI_MAX_PROCESSOR_NAME,MPI_CHARACTER,&
            mpi_comm,ierr)
       if (ierr /=0) call f_err_throw('An error in calling to MPI_ALLGATHER occured',&
            err_id=ERR_MPI_WRAPPERS)


       !found the processors which belong to the same node
       !before the processor iproc
       iproc_node=0
       do jproc=0,iproc-1
          if (trim(nodename(jproc)) == trim(nodename(iproc))) then
             iproc_node=iproc_node+1
          end if
       end do
       nproc_node=iproc_node
       do jproc=iproc,nproc-1
          if (trim(nodename(jproc)) == trim(nodename(iproc))) then
             nproc_node=nproc_node+1
          end if
       end do

       call f_free_str(MPI_MAX_PROCESSOR_NAME,nodename)
    end if
    call f_release_routine()
  END SUBROUTINE mpinoderanks


  subroutine mpihostnames_list(comm,dict)
    use dictionaries
    use dynamic_memory
    implicit none
    integer, intent(in) :: comm
    type(dictionary), pointer :: dict
    !local variables
    integer :: nproc,ierr
    character(len=MPI_MAX_PROCESSOR_NAME) :: nodename_local
    character(len=MPI_MAX_PROCESSOR_NAME), dimension(:), allocatable :: nodename

    nproc=mpisize(comm)

    nodename=f_malloc0_str(MPI_MAX_PROCESSOR_NAME,0 .to. nproc-1,id='nodename')

    nodename_local=mpihostname()

    !gather the result between all the processes
    call MPI_ALLGATHER(nodename_local,MPI_MAX_PROCESSOR_NAME,MPI_CHARACTER,&
         nodename(0),MPI_MAX_PROCESSOR_NAME,MPI_CHARACTER,&
         comm,ierr)
    if (ierr /=0) call f_err_throw('An error in calling to MPI_ALLGATHER occured',&
         err_id=ERR_MPI_WRAPPERS)
    
    dict=>list_new(.item. nodename)

    call f_free_str(MPI_MAX_PROCESSOR_NAME,nodename)

  end subroutine mpihostnames_list

  !> returns true if the mpi has been initialized
  function mpiinitialized()
    use dictionaries, only: f_err_throw
    implicit none
    logical :: mpiinitialized
    !local variables
    logical :: flag
    integer :: ierr

    mpiinitialized=.false.
    call mpi_initialized(flag, ierr)
    if (ierr /=0) then
       flag=.false.
       call f_err_throw('An error in calling to MPI_INITIALIZED occured',&
            err_id=ERR_MPI_WRAPPERS)
    end if
    mpiinitialized=flag

  end function mpiinitialized


  !> Performs the barrier of a given communicator, if present
  subroutine mpibarrier(comm)
    use dictionaries, only: f_err_throw
    implicit none
    integer, intent(in), optional :: comm !< the communicator
    !local variables
    integer :: mpi_comm,ierr

    if (present(comm)) then
       mpi_comm=comm
    else
       mpi_comm=MPI_COMM_WORLD
    end if
    !call the barrier
    call MPI_BARRIER(mpi_comm,ierr)
    if (ierr /=0) then
       call f_err_throw('An error in calling to MPI_BARRIER occured',&
            err_id=ERR_MPI_WRAPPERS)
    end if
  end subroutine mpibarrier

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
    !-end gather-inc
  end subroutine mpigather_c1i2

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
    !-end gather-inc
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

  !>performs gathering of array portions into a receive buffer
  !! the arguments can be provided such as to meet either allgather
  !! or allgatherv APIs. The wrapper chooses the better routine to call
  !! as a function of the arguments
  subroutine mpiallgatherv_d0(sendbuf,sendcount,recvbuf,recvcount,&
       recvcounts,displs,comm)
    use yaml_strings, only: yaml_toa
    use dictionaries, only: f_err_throw
    use dynamic_memory
    implicit none
    double precision, intent(in) :: sendbuf
    double precision, intent(inout), optional :: recvbuf
    double precision, dimension(:), allocatable :: copybuf
    include 'allgather-inc.f90'
  end subroutine mpiallgatherv_d0
  subroutine mpiallgatherv_d1(sendbuf,sendcount,recvbuf,recvcount,&
       recvcounts,displs,comm)
    use yaml_strings, only: yaml_toa
    use dictionaries, only: f_err_throw
    use dynamic_memory
    implicit none
    double precision, dimension(:), intent(in) :: sendbuf
    double precision, dimension(:), intent(inout), optional :: recvbuf
    double precision, dimension(:), allocatable :: copybuf
    include 'allgather-inc.f90'
  end subroutine mpiallgatherv_d1
  subroutine mpiallgatherv_d2d3(sendbuf,sendcount,recvbuf,recvcount,&
       recvcounts,displs,comm)
    use yaml_strings, only: yaml_toa
    use dictionaries, only: f_err_throw
    use dynamic_memory
    implicit none
    double precision, dimension(:,:), intent(in) :: sendbuf
    double precision, dimension(:,:,:), intent(inout), optional :: recvbuf
    double precision, dimension(:), allocatable :: copybuf
    include 'allgather-inc.f90'
  end subroutine mpiallgatherv_d2d3

  subroutine mpiallgatherv_i2(sendbuf,sendcount,recvbuf,recvcount,&
       recvcounts,displs,comm)
    use yaml_strings, only: yaml_toa
    use dictionaries, only: f_err_throw
    use dynamic_memory
    implicit none
    integer, dimension(:,:), intent(in) :: sendbuf
    integer, dimension(:,:), intent(inout), optional :: recvbuf
    integer, dimension(:), allocatable :: copybuf
    include 'allgather-inc.f90'
  end subroutine mpiallgatherv_i2




  subroutine mpialltoallv_int(sendbuf, sendcounts, sdispls, recvbuf, recvcounts, rdispls, comm)
    use dictionaries, only: f_err_throw
    use dynamic_memory
    implicit none
    integer(f_integer),intent(in) :: sendbuf
    integer(f_integer),intent(out) :: recvbuf
    include 'alltoallv-inc.f90'
  end subroutine mpialltoallv_int

  subroutine mpialltoallv_long(sendbuf, sendcounts, sdispls, recvbuf, recvcounts, rdispls, comm)
    use dictionaries, only: f_err_throw,f_err_define
    use dynamic_memory
    implicit none
    integer(f_long),intent(in) :: sendbuf
    integer(f_long),intent(out) :: recvbuf
    include 'alltoallv-inc.f90'
  end subroutine mpialltoallv_long

  subroutine mpialltoallv_double(sendbuf, sendcounts, sdispls, recvbuf, recvcounts, rdispls, comm)
    use dictionaries, only: f_err_throw,f_err_define
    use dynamic_memory
    implicit none
    double precision,intent(in) :: sendbuf
    double precision,intent(out) :: recvbuf
    include 'alltoallv-inc.f90'
  end subroutine mpialltoallv_double


  !> Interface for MPI_ALLREDUCE operations
  subroutine mpiallred_int(sendbuf,count,op,comm,recvbuf)
    use dictionaries, only: f_err_throw,f_err_define
    use dynamic_memory
    implicit none
    integer(f_integer), intent(inout) :: sendbuf
    integer(f_integer), intent(inout), optional :: recvbuf
    integer(f_integer), dimension(:), allocatable :: copybuf
    include 'allreduce-inc.f90'
  end subroutine mpiallred_int

  subroutine mpiallred_long(sendbuf,count,op,comm,recvbuf)
    use dictionaries, only: f_err_throw,f_err_define
    use dynamic_memory
    implicit none
    integer(f_long), intent(inout) :: sendbuf
    integer(f_long), intent(inout), optional :: recvbuf
    integer(f_long), dimension(:), allocatable :: copybuf
    include 'allreduce-inc.f90'
  end subroutine mpiallred_long

  !> Interface for MPI_ALLREDUCE operations
  subroutine mpiallred_real(sendbuf,count,op,comm,recvbuf)
    use dynamic_memory
    use dictionaries, only: f_err_throw,f_err_define
    implicit none
    real, intent(inout) :: sendbuf
    real, intent(inout), optional :: recvbuf
    real, dimension(:), allocatable :: copybuf
    include 'allreduce-inc.f90'
  end subroutine mpiallred_real

  subroutine mpiallred_double(sendbuf,count,op,recvbuf,comm)
    use dynamic_memory
    use dictionaries, only: f_err_throw,f_err_define
    implicit none
    double precision :: sendbuf
    double precision, intent(inout), optional :: recvbuf
    double precision, dimension(:), allocatable :: copybuf
    include 'allreduce-inc.f90'
  end subroutine mpiallred_double

  subroutine mpiallred_log(sendbuf,count,op,comm,recvbuf)
    use dynamic_memory
    use dictionaries, only: f_err_throw,f_err_define
    implicit none
    logical, intent(inout) :: sendbuf
    logical, intent(inout), optional :: recvbuf
    logical, dimension(:), allocatable :: copybuf
    include 'allreduce-inc.f90'
  end subroutine mpiallred_log

  subroutine mpiallred_i1(sendbuf,op,comm,recvbuf)
    use dynamic_memory
    use dictionaries, only: f_err_throw!,f_err_define
    use yaml_strings, only: yaml_toa
    implicit none
    integer, dimension(:), intent(inout) :: sendbuf
    integer, dimension(:), intent(inout), optional :: recvbuf
    integer, dimension(:), allocatable :: copybuf
    include 'allreduce-arr-inc.f90'
  end subroutine mpiallred_i1

  subroutine mpiallred_i2(sendbuf,op,comm,recvbuf)
    use dynamic_memory
    use dictionaries, only: f_err_throw!,f_err_define
    use yaml_strings, only: yaml_toa
    implicit none
    integer, dimension(:,:), intent(inout) :: sendbuf
    integer, dimension(:,:), intent(inout), optional :: recvbuf
    integer, dimension(:,:), allocatable :: copybuf
    include 'allreduce-arr-inc.f90'
  end subroutine mpiallred_i2

  subroutine mpiallred_i3(sendbuf,op,comm,recvbuf)
    use dynamic_memory
    use dictionaries, only: f_err_throw!,f_err_define
    use yaml_strings, only: yaml_toa
    implicit none
    integer, dimension(:,:,:), intent(inout) :: sendbuf
    integer, dimension(:,:,:), intent(inout), optional :: recvbuf
    integer, dimension(:,:,:), allocatable :: copybuf
    include 'allreduce-arr-inc.f90'
  end subroutine mpiallred_i3

  subroutine mpiallred_l3(sendbuf,op,comm,recvbuf)
    use dynamic_memory
    use dictionaries, only: f_err_throw!,f_err_define
    use yaml_strings, only: yaml_toa
    implicit none
    logical, dimension(:,:,:), intent(inout) :: sendbuf
    logical, dimension(:,:,:), intent(inout), optional :: recvbuf
    logical, dimension(:,:,:), allocatable :: copybuf
    include 'allreduce-arr-inc.f90'
  end subroutine mpiallred_l3


  subroutine mpiallred_r1(sendbuf,op,comm,recvbuf)
    use dynamic_memory
    use dictionaries, only: f_err_throw!,f_err_refine
    use yaml_strings, only: yaml_toa
    implicit none
    real, dimension(:), intent(inout) :: sendbuf
    real, dimension(:), intent(inout), optional :: recvbuf
    real, dimension(:), allocatable :: copybuf
    include 'allreduce-arr-inc.f90'
  end subroutine mpiallred_r1

  subroutine mpiallred_r2(sendbuf,op,comm,recvbuf)
    use dynamic_memory
    use dictionaries, only: f_err_throw!,f_err_define
    use yaml_strings, only: yaml_toa
    implicit none
    real, dimension(:,:), intent(inout) :: sendbuf
    real, dimension(:,:), intent(inout), optional :: recvbuf
    real, dimension(:,:), allocatable :: copybuf
    include 'allreduce-arr-inc.f90'
  end subroutine mpiallred_r2

  subroutine mpiallred_r3(sendbuf,op,comm,recvbuf)
    use dynamic_memory
    use dictionaries, only: f_err_throw!,f_err_define
    use yaml_strings, only: yaml_toa
    implicit none
    real, dimension(:,:,:), intent(inout) :: sendbuf
    real, dimension(:,:,:), intent(inout), optional :: recvbuf
    real, dimension(:,:,:), allocatable :: copybuf
    include 'allreduce-arr-inc.f90'
  end subroutine mpiallred_r3

  subroutine mpiallred_r4(sendbuf,op,comm,recvbuf)
    use dynamic_memory
    use dictionaries, only: f_err_throw!,f_err_define
    use yaml_strings, only: yaml_toa
    implicit none
    real, dimension(:,:,:,:), intent(inout) :: sendbuf
    real, dimension(:,:,:,:), intent(inout), optional :: recvbuf
    real, dimension(:,:,:,:), allocatable :: copybuf
    include 'allreduce-arr-inc.f90'
  end subroutine mpiallred_r4

  subroutine mpiallred_d1(sendbuf,op,comm,recvbuf)
    use dynamic_memory
    use dictionaries, only: f_err_throw!,f_err_define
    use yaml_strings, only: yaml_toa
    implicit none
    double precision, dimension(:), intent(inout) :: sendbuf
    double precision, dimension(:), intent(inout), optional :: recvbuf
    double precision, dimension(:), allocatable :: copybuf
    include 'allreduce-arr-inc.f90'
  end subroutine mpiallred_d1

  subroutine mpiallred_d2(sendbuf,op,comm,recvbuf)
    use dynamic_memory
    use dictionaries, only: f_err_throw!,f_err_define
    use yaml_strings, only: yaml_toa
    implicit none
    double precision, dimension(:,:), intent(inout) :: sendbuf
    double precision, dimension(:,:), intent(inout), optional :: recvbuf
    double precision, dimension(:,:), allocatable :: copybuf
    include 'allreduce-arr-inc.f90'
  end subroutine mpiallred_d2

  subroutine mpiallred_d3(sendbuf,op,comm,recvbuf)
    use dynamic_memory
    use dictionaries, only: f_err_throw!,f_err_define
    use yaml_strings, only: yaml_toa
    implicit none
    double precision, dimension(:,:,:), intent(inout) :: sendbuf
    double precision, dimension(:,:,:), intent(inout), optional :: recvbuf
    double precision, dimension(:,:,:), allocatable :: copybuf
    include 'allreduce-arr-inc.f90'
  end subroutine mpiallred_d3

  subroutine mpiallred_d4(sendbuf,op,comm,recvbuf)
    use dynamic_memory
    use dictionaries, only: f_err_throw!,f_err_define
    use yaml_strings, only: yaml_toa
    implicit none
    double precision, dimension(:,:,:,:), intent(inout) :: sendbuf
    double precision, dimension(:,:,:,:), intent(inout), optional :: recvbuf
    double precision, dimension(:,:,:,:), allocatable :: copybuf
    include 'allreduce-arr-inc.f90'
  end subroutine mpiallred_d4

  subroutine mpiallred_d5(sendbuf,op,comm,recvbuf)
    use dynamic_memory
    use dictionaries, only: f_err_throw!,f_err_define
    use yaml_strings, only: yaml_toa
    implicit none
    double precision, dimension(:,:,:,:,:), intent(inout) :: sendbuf
    double precision, dimension(:,:,:,:,:), intent(inout), optional :: recvbuf
    double precision, dimension(:,:,:,:,:), allocatable :: copybuf
    include 'allreduce-arr-inc.f90'
  end subroutine mpiallred_d5


  recursive subroutine mpibcast_i0(buffer,count,root,comm,check,maxdiff)
    use dynamic_memory
    use dictionaries, only: f_err_throw
    use yaml_output !for check=.true.
    use f_utils, only: f_zero
    implicit none
    integer(f_integer) ::  buffer
    integer(f_integer), intent(out), optional :: maxdiff
    integer(f_integer), dimension(:), allocatable :: array_diff
    include 'bcast-decl-inc.f90'
    include 'bcast-inc.f90'
  end subroutine mpibcast_i0

  subroutine mpibcast_li0(buffer,count,root,comm,check,maxdiff)
    use dynamic_memory
    use dictionaries, only: f_err_throw
    use yaml_output !for check=.true.
    use f_utils, only: f_zero
    implicit none
    integer(f_long) ::  buffer
    integer(f_long), intent(out), optional :: maxdiff
    integer(f_long), dimension(:), allocatable :: array_diff
    include 'bcast-decl-inc.f90'
    include 'bcast-inc.f90'
  end subroutine mpibcast_li0

  recursive subroutine mpibcast_d0(buffer,count,root,comm,check,maxdiff)
    use dynamic_memory
    use dictionaries, only: f_err_throw
    use yaml_output !for check=.true.
    use f_utils, only: f_zero
    implicit none
    double precision ::  buffer
    double precision, intent(out), optional :: maxdiff
    double precision, dimension(:), allocatable :: array_diff
    include 'bcast-decl-inc.f90'
    include 'bcast-inc.f90'
  end subroutine mpibcast_d0

  subroutine mpibcast_c0(buffer,root,comm,check,maxdiff)
    use dynamic_memory
    use dictionaries, only: f_err_throw
    use yaml_output !for check=.true.
    use f_utils, only: f_zero
    implicit none
    character(len=*) ::  buffer
    integer(f_integer), intent(out), optional :: maxdiff
    integer(f_integer), dimension(:), allocatable :: array_diff !<the difference is performed with ascii value
    ! 'bcast-decl-arr-inc.f90'
    integer, intent(in), optional :: root  !< @copydoc doc::root
    integer, intent(in), optional :: comm  !< @copydoc doc::comm
    logical, intent(in), optional :: check !< performs the check of the arguments
    !local variables
    logical chk
    integer :: n,iroot,mpi_comm,ierr
    integer, dimension(3) :: iarg_check
    external :: MPI_BCAST

    chk=.false.
    n=len(buffer)
    if (present(maxdiff)) then
       call f_zero(maxdiff)
       array_diff=f_malloc(n,id='array_diff')
       call f_memcpy(src=buffer,dest=array_diff)
    end if
    ! end bcast_decl
    include 'bcast-inc.f90'
  end subroutine mpibcast_c0

  subroutine mpibcast_c1(buffer,root,comm,check,maxdiff)
    use dynamic_memory
    use dictionaries, only: f_err_throw
    use yaml_output !for check=.true.
    use f_utils, only: f_zero
    implicit none
    character, dimension(:), intent(inout) ::  buffer
    integer, intent(out), optional :: maxdiff
    integer, dimension(:), allocatable :: array_diff !<the difference is performed with ascii value
    include 'bcast-decl-arr-inc.f90'
    include 'bcast-inc.f90'
  end subroutine mpibcast_c1

  subroutine mpibcast_i1(buffer,root,comm,check,maxdiff)
    use dynamic_memory
    use dictionaries, only: f_err_throw
    use yaml_output !for check=.true.
    use f_utils, only: f_zero
    implicit none
    integer, dimension(:), intent(inout) ::  buffer
    integer, intent(out), optional :: maxdiff
    integer, dimension(:), allocatable :: array_diff
    include 'bcast-decl-arr-inc.f90'
    include 'bcast-inc.f90'
  end subroutine mpibcast_i1

  subroutine mpibcast_i2(buffer,root,comm,check,maxdiff)
    use dynamic_memory
    use dictionaries, only: f_err_throw
    use yaml_output !for check=.true.
    use f_utils, only: f_zero
    implicit none
    integer, dimension(:,:), intent(inout) ::  buffer
    integer, intent(out), optional :: maxdiff
    integer, dimension(:), allocatable :: array_diff
    include 'bcast-decl-arr-inc.f90'
    include 'bcast-inc.f90'
  end subroutine mpibcast_i2

  subroutine mpibcast_d1(buffer,root,comm,check,maxdiff)
    use dynamic_memory
    use dictionaries, only: f_err_throw
    use yaml_output !for check=.true.
    use f_utils, only: f_zero
    implicit none
    double precision, dimension(:), intent(inout) ::  buffer
    double precision, intent(out), optional :: maxdiff
    double precision, dimension(:), allocatable :: array_diff
    include 'bcast-decl-arr-inc.f90'
    include 'bcast-inc.f90'
  end subroutine mpibcast_d1

  subroutine mpibcast_d2(buffer,root,comm,check,maxdiff)
    use dynamic_memory
    use dictionaries, only: f_err_throw
    use yaml_output !for check=.true.
    use f_utils, only: f_zero
    implicit none
    double precision, dimension(:,:), intent(inout) ::  buffer
    double precision, intent(out), optional :: maxdiff
    double precision, dimension(:), allocatable :: array_diff
    include 'bcast-decl-arr-inc.f90'
    include 'bcast-inc.f90'
  end subroutine mpibcast_d2

  subroutine mpiscatter_i1i1(sendbuf, recvbuf, root, comm)
    use dictionaries, only: f_err_throw,f_err_define
    implicit none
    integer,dimension(:),intent(in) :: sendbuf
    integer,dimension(:),intent(inout) :: recvbuf
    include 'scatter-inc.f90'
  end subroutine mpiscatter_i1i1

  subroutine mpiscatterv_d0(sendbuf, sendcounts, displs, recvbuf, recvcount, root, comm)
    use dictionaries, only: f_err_throw
    implicit none
    real(f_double) :: sendbuf
    real(f_double), intent(inout) :: recvbuf
    include 'scatterv-decl-inc.f90'
    include 'scatterv-inc.f90'
  end subroutine mpiscatterv_d0

  subroutine mpiscatterv_d2d3(sendbuf, sendcounts, displs, recvbuf, root, comm)
    use dictionaries, only: f_err_throw
    implicit none
    real(f_double), dimension(:,:), intent(in) :: sendbuf
    real(f_double), dimension(:,:,:), intent(out) :: recvbuf
    include 'scatterv-decl-inc.f90'
    recvcount=size(sendbuf)
    include 'scatterv-inc.f90'
  end subroutine mpiscatterv_d2d3

  subroutine mpiscatterv_d3d2(sendbuf, sendcounts, displs, recvbuf, root, comm)
    use dictionaries, only: f_err_throw
    implicit none
    real(f_double), dimension(:,:,:), intent(in) :: sendbuf
    real(f_double), dimension(:,:), intent(out) :: recvbuf
    include 'scatterv-decl-inc.f90'
    recvcount=size(sendbuf)
    include 'scatterv-inc.f90'
  end subroutine mpiscatterv_d3d2


  !> Detect the maximum difference between arrays all over a given communicator
  function mpimaxdiff_i0(n,array,root,source,comm,bcast) result(maxdiff)
    use dynamic_memory
    implicit none
    integer, intent(in) :: n !<number of elements to be controlled
    integer(f_integer), intent(inout) :: array !< starting point of the array
    integer(f_integer), dimension(:,:), allocatable :: array_glob
    integer(f_integer) :: maxdiff
    include 'maxdiff-decl-inc.f90'

    ndims = n
    maxdiff=0

    include 'maxdiff-inc.f90'
  end function mpimaxdiff_i0

  !> Detect the maximum difference between arrays all over a given communicator
  function mpimaxdiff_li0(n,array,root,source,comm,bcast) result(maxdiff)
    use dynamic_memory
    implicit none
    integer(f_long), intent(in) :: n !<number of elements to be controlled
    integer(f_long), intent(inout) :: array !< starting point of the array
    integer(f_long), dimension(:,:), allocatable :: array_glob
    integer(f_long) :: maxdiff
    include 'maxdiff-decl-inc.f90'
    ndims = int(n,kind=4)
    maxdiff=int(0,f_long)
    include 'maxdiff-inc.f90'
  end function mpimaxdiff_li0

  function mpimaxdiff_c1(array,root,source,comm,bcast) result(maxdiff)
    use dynamic_memory
    implicit none
    !> array to be checked
    character, dimension(:), intent(inout) :: array
    integer, dimension(:,:), allocatable :: array_glob
    integer :: maxdiff
    include 'maxdiff-decl-inc.f90'
    ndims = size(array)
    maxdiff=0
    include 'maxdiff-arr-inc.f90'
  end function mpimaxdiff_c1

  function mpimaxdiff_d0(n,array,root,source,comm,bcast) result(maxdiff)
    use dynamic_memory
    implicit none
    integer, intent(in) :: n !<number of elements to be controlled
    double precision, intent(inout) :: array !< starting point of the array
    double precision, dimension(:,:), allocatable :: array_glob
    double precision :: maxdiff
    include 'maxdiff-decl-inc.f90'

    ndims = n
    maxdiff=0.d0

    include 'maxdiff-inc.f90'
  end function mpimaxdiff_d0

  function mpimaxdiff_d1(array,root,source,comm,bcast) result(maxdiff)
    use dynamic_memory
    implicit none
    !> array to be checked
    double precision, dimension(:), intent(inout) :: array
    double precision, dimension(:,:), allocatable :: array_glob
    double precision :: maxdiff
    include 'maxdiff-decl-inc.f90'

    ndims = size(array)

    maxdiff=0.d0

    include 'maxdiff-arr-inc.f90'
  end function mpimaxdiff_d1

  function mpimaxdiff_i1(array,root,source,comm,bcast) result(maxdiff)
    use dynamic_memory
    implicit none
    !> array to be checked
    integer(f_integer), dimension(:), intent(inout) :: array
    integer(f_integer), dimension(:,:), allocatable :: array_glob
    integer(f_integer) :: maxdiff
    include 'maxdiff-decl-inc.f90'

    ndims = size(array)

    maxdiff=0

    include 'maxdiff-arr-inc.f90'
  end function mpimaxdiff_i1

  function mpimaxdiff_li1(array,root,source,comm,bcast) result(maxdiff)
    use dynamic_memory
    implicit none
    !> array to be checked
    integer(f_long), dimension(:), intent(inout) :: array
    integer(f_long), dimension(:,:), allocatable :: array_glob
    integer(f_long) :: maxdiff
    include 'maxdiff-decl-inc.f90'

    ndims = size(array)

    maxdiff=0

    include 'maxdiff-arr-inc.f90'
  end function mpimaxdiff_li1

  function mpimaxdiff_i2(array,root,source,comm,bcast) result(maxdiff)

    use dynamic_memory
    implicit none
    !> array to be checked
    integer(f_integer), dimension(:,:), intent(inout) :: array
    integer(f_integer), dimension(:,:), allocatable :: array_glob
    integer :: maxdiff
    include 'maxdiff-decl-inc.f90'

    ndims = size(array)

    maxdiff=0

    include 'maxdiff-arr-inc.f90'
  end function mpimaxdiff_i2


  function mpimaxdiff_d2(array,root,source,comm,bcast) result(maxdiff)
    use dynamic_memory
    implicit none
    !> array to be checked
    double precision, dimension(:,:), intent(inout) :: array
    double precision, dimension(:,:), allocatable :: array_glob
    double precision :: maxdiff
    include 'maxdiff-decl-inc.f90'

    ndims = size(array)

    maxdiff=0.d0

    include 'maxdiff-arr-inc.f90'
  end function mpimaxdiff_d2

  !!function mpitypesize_d(foo) result(sizeof)
  !!  use dictionaries, only: f_err_throw,f_err_define
  !!  implicit none
  !!  double precision, intent(in) :: foo
  !!  integer :: sizeof, ierr

  !!  call mpi_type_size(mpi_double_precision, sizeof, ierr)
  !!  if (ierr/=0) then
  !!      call f_err_throw('Error in mpi_type_size',&
  !!           err_id=ERR_MPI_WRAPPERS)
  !!  end if
  !!end function mpitypesize_d

  function mpitypesize_d0(foo) result(sizeof)
    use dictionaries, only: f_err_throw,f_err_define
    implicit none
    double precision, intent(in) :: foo
    integer :: sizeof, ierr

    call mpi_type_size(mpi_double_precision, sizeof, ierr)
    if (ierr/=0) then
       call f_err_throw('Error in mpi_type_size',&
            err_id=ERR_MPI_WRAPPERS)
    end if
  end function mpitypesize_d0

  function mpitypesize_d1(foo) result(sizeof)
    implicit none
    double precision, dimension(:), intent(in) :: foo
    integer :: sizeof
    sizeof=mpitypesize(1.d0)
  end function mpitypesize_d1

  function mpitypesize_i0(foo) result(sizeof)
    use dictionaries, only: f_err_throw,f_err_define
    implicit none
    integer, intent(in) :: foo
    integer :: sizeof, ierr

    call mpi_type_size(mpi_integer, sizeof, ierr)
    if (ierr/=0) then
       call f_err_throw('Error in mpi_type_size',&
            err_id=ERR_MPI_WRAPPERS)
    end if
  end function mpitypesize_i0

  function mpitypesize_l0(foo) result(sizeof)
    use dictionaries, only: f_err_throw,f_err_define
    implicit none
    logical, intent(in) :: foo
    integer :: sizeof, ierr

    call mpi_type_size(mpi_logical, sizeof, ierr)
    if (ierr/=0) then
       call f_err_throw('Error in mpi_type_size',&
            err_id=ERR_MPI_WRAPPERS)
    end if
  end function mpitypesize_l0

  function mpiinfo(key,val) result(info)
    use dictionaries, only: f_err_throw,f_err_define
    implicit none
    character(len=*), intent(in) :: key
    character(len=*), intent(in) :: val
    integer :: info, ierr

    call mpi_info_create(info, ierr)
    if (ierr/=0) then
       call f_err_throw('Error in mpi_info_create',&
            err_id=ERR_MPI_WRAPPERS)
       return
    end if
    call mpi_info_set(info, "no_locks", "true", ierr)
    if (ierr/=0) then
       !!call f_err_throw('Error in mpi_info_set, key='//trim(key)//&
       !!     ', value=',trim(val),err_id=ERR_MPI_WRAPPERS)
       call f_err_throw('Error in mpi_info_set, key='//trim(key)//&
            ', value='//trim(val),err_id=ERR_MPI_WRAPPERS)
    end if

  end function mpiinfo

  subroutine mpiinfofree(info)
    use dictionaries, only: f_err_throw,f_err_define
    implicit none
    integer, intent(inout) :: info
    ! Local variables
    integer :: ierr
    call mpi_info_free(info, ierr)
    if (ierr/=0) then
       call f_err_throw('Error in mpi_info_free',&
            err_id=ERR_MPI_WRAPPERS)
    end if
  end subroutine mpiinfofree

  function mpiwindow_d0(size,base,comm) result(window)
    use dictionaries, only: f_err_throw,f_err_define
    implicit none
    integer,intent(in) :: size
    double precision,intent(in) :: base
    integer,intent(in) :: comm
    !local variables
    integer :: sizeof,info,ierr
    integer :: window

    sizeof=mpitypesize(base)
    info=mpiinfo("no_locks", "true")

    call mpi_win_create(base, int(size,kind=mpi_address_kind)*int(sizeof,kind=mpi_address_kind), &
         sizeof, info,comm, window, ierr)

    if (ierr/=0) then
       call f_err_throw('Error in mpi_win_create',&
            err_id=ERR_MPI_WRAPPERS)
    end if

    call mpiinfofree(info)

    call mpi_win_fence(MPI_MODE_NOPRECEDE, window, ierr)
    if (ierr/=0) then
       call f_err_throw('Error in mpi_win_fence',&
            err_id=ERR_MPI_WRAPPERS)
    end if


  end function mpiwindow_d0

  function mpiwindow_i0(size,base,comm) result(window)
    use dictionaries, only: f_err_throw,f_err_define
    implicit none
    integer,intent(in) :: size
    integer,intent(in) :: base
    integer,intent(in) :: comm
    !local variables
    integer :: sizeof,info,ierr
    integer :: window

    sizeof=mpitypesize(base)
    info=mpiinfo("no_locks", "true")

    call mpi_win_create(base, int(size,kind=mpi_address_kind)*int(sizeof,kind=mpi_address_kind), &
         sizeof, info,comm, window, ierr)

    if (ierr/=0) then
       call f_err_throw('Error in mpi_win_create',&
            err_id=ERR_MPI_WRAPPERS)
    end if

    call mpiinfofree(info)

    call mpi_win_fence(MPI_MODE_NOPRECEDE, window, ierr)
    if (ierr/=0) then
       call f_err_throw('Error in mpi_win_fence',&
            err_id=ERR_MPI_WRAPPERS)
    end if


  end function mpiwindow_i0

  function mpiwindow_l0(size,base,comm) result(window)
    use dictionaries, only: f_err_throw,f_err_define
    implicit none
    integer,intent(in) :: size
    logical,intent(in) :: base
    integer,intent(in) :: comm
    !local variables
    integer :: sizeof,info,ierr
    integer :: window

    sizeof=mpitypesize(base)
    info=mpiinfo("no_locks", "true")

    call mpi_win_create(base, int(size,kind=mpi_address_kind)*int(sizeof,kind=mpi_address_kind), &
         sizeof, info,comm, window, ierr)

    if (ierr/=0) then
       call f_err_throw('Error in mpi_win_create',&
            err_id=ERR_MPI_WRAPPERS)
    end if

    call mpiinfofree(info)

    call mpi_win_fence(MPI_MODE_NOPRECEDE, window, ierr)
    if (ierr/=0) then
       call f_err_throw('Error in mpi_win_fence',&
            err_id=ERR_MPI_WRAPPERS)
    end if


  end function mpiwindow_l0

  !> create a peer_to_peer group, to use RMA calls instead of send-receive
  function p2p_group(base_grp,p1,p2,p3) result(grp)
    use yaml_strings, only: yaml_toa
    implicit none
    integer, intent(in) :: base_grp
    integer, intent(in) :: p1,p2
    integer, intent(in), optional :: p3
    integer :: grp
    !local variables
    integer :: i,nlist
    integer, dimension(3) :: list,ipiv,list2

    if (present(p3)) then
       list(1)=p1
       list(2)=p2
       list(3)=p3
       call sort_positions(3,real(list,kind=8),ipiv)
       nlist=3
       do i=1,3
          if (i > 1) then
             if (list(ipiv(i))==list2(i-1)) then
                nlist=nlist-1
             else
                list2(i)=list(ipiv(i))
             end if
          else
             list2(i)=list(ipiv(i))
          end if
       end do
       grp=mpigroupincl(base_grp,nlist,list2)
    else
       list(1)=min(p1,p2)
       list(2)=max(p1,p2)
       grp=mpigroupincl(base_grp,2,list)
    end if
    if (grp==MPI_GROUP_NULL) then
       call f_err_throw('Error in the group creation for list='//trim(yaml_toa(list)),&
            err_id=ERR_MPI_WRAPPERS)
    end if

  end function p2p_group


  subroutine mpi_fence(window, assert)
    use dictionaries, only: f_err_throw,f_err_define
    ! Calling arguments
    integer,intent(inout) :: window !<window to be synchronized
    integer,intent(in),optional :: assert

    ! Local variables
    integer :: ierr, assert_, tcat

    if (present(assert)) then
       assert_ = assert
    else
       assert_ = 0
    end if
    tcat=TCAT_FENCE
    ! Synchronize the communication
    call f_timer_interrupt(tcat)
    call mpi_win_fence(assert_, window, ierr)
    call f_timer_resume()
    if (ierr/=0) then
       call f_err_throw('Error in mpi_win_fence',&
            err_id=ERR_MPI_WRAPPERS)
    end if
  end subroutine mpi_fence

  subroutine mpiwinstart(grp,win,assert)
    implicit none
    integer, intent(in) :: grp
    integer, intent(in) :: win
    integer, intent(in), optional :: assert
    !local variables
    integer :: assert_,ierr
    assert_=0
    if (present(assert)) assert_=assert

    if (grp==mpigroup_null()) then
       call f_err_throw('Error in mpi_win_start, passed a null group',&
            err_id=ERR_MPI_WRAPPERS)
    end if


    call MPI_WIN_START(grp,assert_,win,ierr)
    if (ierr /=0) then
       call f_err_throw('Error in mpi_win_start',&
            err_id=ERR_MPI_WRAPPERS)
    end if

  end subroutine mpiwinstart

  subroutine mpiwinpost(grp,win,assert)
    implicit none
    integer, intent(in) :: grp
    integer, intent(in) :: win
    integer, intent(in), optional :: assert
    !local variables
    integer :: assert_,ierr
    assert_=0
    if (present(assert)) assert_=assert

    if (grp==mpigroup_null()) then
       call f_err_throw('Error in mpi_win_post, passed a null group',&
            err_id=ERR_MPI_WRAPPERS)
    end if

    call MPI_WIN_POST(grp,assert_,win,ierr)
    if (ierr /=0) then
       call f_err_throw('Error in mpi_win_post',&
            err_id=ERR_MPI_WRAPPERS)
    end if

  end subroutine mpiwinpost

  subroutine mpiwincomplete(win)
    implicit none
    integer, intent(in) :: win
    !local variables
    integer :: ierr

    call MPI_WIN_COMPLETE(win,ierr)
    if (ierr /=0) then
       call f_err_throw('Error in mpi_win_complete',&
            err_id=ERR_MPI_WRAPPERS)
    end if

  end subroutine mpiwincomplete

  subroutine mpiwinwait(win)
    implicit none
    integer, intent(in) :: win
    !local variables
    integer :: ierr

    call MPI_WIN_WAIT(win,ierr)
    if (ierr /=0) then
       call f_err_throw('Error in mpi_win_wait',&
            err_id=ERR_MPI_WRAPPERS)
    end if

  end subroutine mpiwinwait


  subroutine mpi_fenceandfree(window, assert)
    use dictionaries, only: f_err_throw,f_err_define
    ! Calling arguments
    integer,intent(inout) :: window !<window to be synchronized and freed
    integer,intent(in),optional :: assert

    ! Local variables
    integer :: ierr, assert_

    if (present(assert)) then
       assert_ = assert
    else
       assert_ = 0
    end if

    ! Synchronize the communication
    call mpi_win_fence(assert_, window, ierr)
    if (ierr/=0) then
       call f_err_throw('Error in mpi_win_fence',&
            err_id=ERR_MPI_WRAPPERS)
    end if
    call mpi_win_free(window, ierr)
    if (ierr/=0) then
       call f_err_throw('Error in mpi_win_fence',&
            err_id=ERR_MPI_WRAPPERS)
    end if
  end subroutine mpi_fenceandfree

  subroutine mpiget_d0(origin,count,target_rank,target_disp,window)
    use dictionaries, only: f_err_throw,f_err_define
    implicit none
    double precision,intent(inout) :: origin !<fake intent(in)
    integer,intent(in) :: count, target_rank,window
    integer(kind=mpi_address_kind),intent(in) :: target_disp

    ! Local variables
    integer :: ierr

    call mpi_get(origin,count,mpitype(origin),target_rank, &
         target_disp,count,mpitype(origin), window, ierr)
    if (ierr/=0) then
       call f_err_throw('Error in mpi_get',&
            err_id=ERR_MPI_WRAPPERS)
    end if
  end subroutine mpiget_d0

  subroutine mpiput_d0(origin,count,target_rank,target_disp,window)
    use dictionaries, only: f_err_throw,f_err_define
    implicit none
    double precision,intent(inout) :: origin !<fake intent(in)
    integer,intent(in) :: count, target_rank,window
    integer(kind=mpi_address_kind),intent(in) :: target_disp

    ! Local variables
    integer :: ierr

    call mpi_put(origin,count,mpitype(origin),target_rank, &
         target_disp,count,mpitype(origin), window, ierr)
    if (ierr/=0) then
       call f_err_throw('Error in mpi_put',&
            err_id=ERR_MPI_WRAPPERS)
    end if
  end subroutine mpiput_d0

  subroutine mpiaccumulate_d0(origin,count,target_rank,target_disp,op,window)
    use dictionaries, only: f_err_throw,f_err_define
    implicit none
    double precision,intent(inout) :: origin !<fake intent(in)
    integer,intent(in) :: count, target_rank,window,op
    integer(kind=mpi_address_kind),intent(in) :: target_disp

    ! Local variables
    integer :: ierr

    call mpi_accumulate(origin,count,mpitype(origin),target_rank, &
         target_disp,count,mpitype(origin), op, window, ierr)
    if (ierr/=0) then
       call f_err_throw('Error in mpi_accumulate',&
            err_id=ERR_MPI_WRAPPERS)
    end if
  end subroutine mpiaccumulate_d0


  subroutine mpi_get_to_allgatherv_double(sendbuf,sendcount,recvbuf,recvcounts,displs,comm,check_,window_)
    use dictionaries, only: f_err_throw,f_err_define
    use yaml_strings, only: yaml_toa
    implicit none
    double precision,intent(in) :: sendbuf
    double precision,intent(inout) :: recvbuf
    integer,dimension(:),intent(in) :: recvcounts, displs
    integer,intent(in) :: comm, sendcount
    logical,intent(in),optional :: check_
    integer,intent(out),pointer,optional :: window_
    !local variables
    integer :: nproc,nrecvbuf
    !external :: getall
    logical :: check
    integer,target:: window

    nproc=mpisize(comm)
    nrecvbuf=sum(recvcounts)

    if (present(check_)) then
       check = check_
    else
       check = .false.
    end if

    if (check) then
       !check coherence
       if (any([size(recvcounts),size(displs)] /= nproc)) then
          call f_err_throw("Error in get_to_gatherv, sizes not coherent with communicator"//&
               trim(yaml_toa([size(recvcounts),size(displs), nproc])),&
               err_id=ERR_MPI_WRAPPERS)
          return
       end if
    end if

    if (present(window_)) then
       window_ => window
    end if
    window = mpiwindow(sendcount,sendbuf,comm)


    call getall_d(nproc,recvcounts,displs,window,nrecvbuf,recvbuf)

    if (.not. present(window_)) then
       call mpi_fenceandfree(window)
    end if

  end subroutine mpi_get_to_allgatherv_double


!!$  subroutine mpiaccumulate_double(origin_addr, origin_count, target_rank, target_disp, target_count, op, wind)
!!$    use dictionaries, only: f_err_throw,f_err_define
!!$    use yaml_strings, only: yaml_toa
!!$    implicit none
!!$    double precision,intent(in) :: origin_addr
!!$    integer,intent(in) :: origin_count, target_rank, target_count, op
!!$    integer(kind=mpi_address_kind),intent(in) :: target_disp
!!$    integer,intent(inout) :: wind
!!$    !local variables
!!$    integer :: nproc,jproc,nrecvbuf,ierr
!!$    external :: getall
!!$    logical :: check
!!$    integer,target:: window
!!$
!!$
!!$    call mpi_accumulate(origin_addr, origin_count, mpitype(origin_addr), &
!!$         target_rank, target_disp, target_count, mpitype(origin_addr), op, wind, ierr)
!!$    if (ierr/=0) then
!!$       call f_err_throw('An error in calling to MPI_ACCUMULATE occured',&
!!$            err_id=ERR_MPI_WRAPPERS)
!!$       return
!!$    end if
!!$
!!$  end subroutine mpiaccumulate_double
!!$

  subroutine mpiiallred_double(sendbuf, recvbuf, ncount, op, comm, request)
    use dictionaries, only: f_err_throw,f_err_define
    implicit none
    ! Calling arguments
    integer,intent(in) :: ncount, op, comm
    double precision,intent(in) :: sendbuf
    double precision,intent(out) :: recvbuf
    integer,intent(out) :: request
    ! Local variables
    integer :: ierr

#ifdef HAVE_MPI3
    call mpi_iallreduce(sendbuf, recvbuf, ncount, mpitype(sendbuf), op, comm, request, ierr)
    if (ierr/=0) then
       call f_err_throw('An error in calling to MPI_IALLREDUCE occured',&
            err_id=ERR_MPI_WRAPPERS)
       return
    end if
#else
    call mpi_allreduce(sendbuf, recvbuf, ncount, mpitype(sendbuf), op, comm, ierr)
    if (ierr/=0) then
       call f_err_throw('An error in calling to MPI_ALLREDUCE occured',&
            err_id=ERR_MPI_WRAPPERS)
       return
    end if
    request = MPI_REQUEST_NULL
#endif

  end subroutine mpiiallred_double


  subroutine mpiialltoallv_double(sendbuf, sendcounts, senddspls, sendtype, &
       recvbuf, recvcounts, recvdspls, recvtype, comm, request)
    use dictionaries, only: f_err_throw,f_err_define
    implicit none
    ! Calling arguments
    integer,intent(in) :: sendcounts, senddspls, sendtype, recvcounts, recvdspls, recvtype, comm
    double precision,intent(in) :: sendbuf
    double precision,intent(out) :: recvbuf
    integer,intent(out) :: request
    ! Local variables
    integer :: ierr

#ifdef HAVE_MPI3
    call mpi_ialltoallv(sendbuf, sendcounts, senddspls, sendtype, &
         recvbuf, recvcounts, recvdspls, recvtype, comm, request, ierr)
    if (ierr/=0) then
       call f_err_throw('An error in calling to MPI_IALLTOALLV occured',&
            err_id=ERR_MPI_WRAPPERS)
       return
    end if
#else
    call mpi_alltoallv(sendbuf, sendcounts, senddspls, sendtype, &
         recvbuf, recvcounts, recvdspls, recvtype, comm, ierr)
    if (ierr/=0) then
       call f_err_throw('An error in calling to MPI_IALLTOALLV occured',&
            err_id=ERR_MPI_WRAPPERS)
       return
    end if
    request = MPI_REQUEST_NULL
#endif

  end subroutine mpiialltoallv_double

  subroutine mpisend_d0(buf,count,dest,tag,comm,request,simulate,verbose,type)
    use yaml_output
    implicit none
    real(f_double) :: buf !fake intent(in)
    integer, intent(in) :: count
    integer, intent(in) :: dest
    integer, intent(in), optional :: tag
    integer, intent(in), optional :: comm
    integer, intent(out), optional :: request !<toggle the isend operation
    logical, intent(in), optional :: simulate,verbose
    integer, intent(in), optional :: type
    !local variables
    logical :: verb,sim
    integer :: mpi_comm,ierr,tag_,tcat

    mpi_comm=MPI_COMM_WORLD
    if (present(comm)) mpi_comm=comm
    if (present(tag)) then
       tag_=tag
    else
       tag_=mpirank(mpi_comm)
    end if

    verb=.false.
    if (present(verbose)) verb=verbose .and. dest /=mpirank_null()

    if (verb) then
       call yaml_mapping_open('MPI_(I)SEND')
       call yaml_map('Elements',count)
       call yaml_map('Source',mpirank(mpi_comm))
       call yaml_map('Dest',dest)
       call yaml_map('Tag',tag_)
       call yaml_mapping_close()
    end if

    sim=.false.
    if (present(simulate)) sim=simulate
    if (sim) return

    tcat=TCAT_SEND
    ! Synchronize the communication
    call f_timer_interrupt(tcat)

    if (present(type)) then
      if (present(request)) then
         call MPI_ISEND(buf,count,type,dest,tag,mpi_comm,request,ierr)
      else
         call MPI_SEND(buf,count,type,dest,tag,mpi_comm,ierr)
      end if
    else
      if (present(request)) then
         call MPI_ISEND(buf,count,mpitype(buf),dest,tag,mpi_comm,request,ierr)
      else
         call MPI_SEND(buf,count,mpitype(buf),dest,tag,mpi_comm,ierr)
      end if
    end if
    call f_timer_resume()
    if (ierr/=0) call f_err_throw('An error in calling to MPI_(I)SEND occured',&
            err_id=ERR_MPI_WRAPPERS)

  end subroutine mpisend_d0

  subroutine mpisend_gpu(buf,count,dest,tag,comm,request,simulate,verbose,type,offset)
    use yaml_output
    use iso_c_binding
    use f_precisions, only: f_address
    implicit none
    type(c_ptr) :: buf !fake intent(in)
    integer, intent(in) :: count
    integer, intent(in) :: dest
    integer, intent(in), optional :: tag
    integer, intent(in), optional :: comm
    integer, intent(out), optional :: request !<toggle the isend operation
    logical, intent(in), optional :: simulate,verbose
    integer, intent(in) :: type
    integer, intent(in), optional :: offset
    real(f_double),pointer :: a !fake intent(in)
    !local variables
    logical :: verb,sim
    integer :: mpi_comm,ierr,tag_,tmpsize,tcat
    integer(f_address) tmpint
    type(c_ptr) :: tmpaddr

    mpi_comm=MPI_COMM_WORLD
    if (present(comm)) mpi_comm=comm
    if (present(tag)) then
       tag_=tag
    else
       tag_=mpirank(mpi_comm)
    end if

    verb=.false.
    if (present(verbose)) verb=verbose .and. dest /=mpirank_null()

    if (verb) then
       call yaml_mapping_open('MPI_(I)SEND')
       call yaml_map('Elements',count)
       call yaml_map('Source',mpirank(mpi_comm))
       call yaml_map('Dest',dest)
       call yaml_map('Tag',tag_)
       call yaml_mapping_close()
    end if

    sim=.false.
    if (present(simulate)) sim=simulate
    if (sim) return

    tcat=TCAT_SEND
    ! Synchronize the communication
    call f_timer_interrupt(tcat)
    
    !LG: this cannot be written like that (segfault on some compilers, see fortran spec)
    !if(present(offset) .and. offset/=0)then
    if (present(offset)) then
       if (offset /=0) then
          tmpint = TRANSFER(buf, tmpint)
          call mpi_type_size(type, tmpsize, ierr)
          tmpint = tmpint + offset*tmpsize
          tmpaddr= TRANSFER(tmpint, tmpaddr)
          call c_f_pointer(tmpaddr, a)
       end if
    else
      call c_f_pointer(buf, a)
    end if
    if (present(request)) then
       call MPI_ISEND(a,count,type,dest,tag,mpi_comm,request,ierr)
    else
       call MPI_SEND(a,count,type,dest,tag,mpi_comm,ierr)
    end if
    call f_timer_resume()

    if (ierr/=0) call f_err_throw('An error in calling to MPI_(I)SEND occured',&
            err_id=ERR_MPI_WRAPPERS)

  end subroutine mpisend_gpu

  subroutine mpirecv_d0(buf,count,source,tag,comm,status,request,simulate,verbose,type)
    use yaml_output
    implicit none
    real(f_double), intent(inout) :: buf !fake intent(out)
    integer, intent(in) :: count
    integer, intent(in), optional :: source
    integer, intent(in), optional :: tag
    integer, intent(in), optional :: comm
    integer, intent(out), optional :: request !<toggle the isend operation
    integer, dimension(MPI_STATUS_SIZE), intent(out), optional :: status !<for the blocking operation
    logical, intent(in), optional :: simulate,verbose
    integer, intent(in), optional :: type
    !local variables
    logical :: verb,sim
    integer :: mpi_comm,ierr,mpi_source,mpi_tag,mpi_type,tcat

    mpi_comm=MPI_COMM_WORLD
    if (present(comm)) mpi_comm=comm
    mpi_source=MPI_ANY_SOURCE
    mpi_tag=MPI_ANY_TAG
    if (present(source)) then
       mpi_source=source
       mpi_tag=source
    end if
    if (present(tag)) mpi_tag=tag
    verb=.false.
    if (present(verbose)) verb=verbose .and. source /= mpirank_null()
    if (verb) call yaml_comment('Receiving'//count//'elements from'//source//'in'//mpirank(mpi_comm))

    if (verb) then
       call yaml_mapping_open('MPI_(I)RECV')
       call yaml_map('Elements',count)
       call yaml_map('Source',source)
       call yaml_map('Dest',mpirank(mpi_comm))
       call yaml_map('Tag',mpi_tag)
       call yaml_mapping_close()
    end if

    sim=.false.
    if (present(simulate)) sim=simulate
    if (sim) return

    tcat=TCAT_RECV
    ! Synchronize the communication
    call f_timer_interrupt(tcat)

    if (present(type)) then
      mpi_type=type
    else
      mpi_type=mpitype(buf)
    end if
    if (present(request)) then
       call MPI_IRECV(buf,count,mpi_type,mpi_source,mpi_tag,mpi_comm,request,ierr)
    else
       if (present(status)) then
          call MPI_RECV(buf,count,mpi_type,mpi_source,mpi_tag,mpi_comm,status,ierr)
       else
          call MPI_RECV(buf,count,mpi_type,mpi_source,mpi_tag,mpi_comm,MPI_STATUS_IGNORE,ierr)
       end if
    end if
    call f_timer_resume()
    if (ierr/=0) call f_err_throw('An error in calling to MPI_(I)RECV occured',&
         err_id=ERR_MPI_WRAPPERS)

  end subroutine mpirecv_d0

  subroutine mpirecv_gpu(buf,count,source,tag,comm,status,request,simulate,verbose,type)
    use yaml_output
    use iso_c_binding
    implicit none
    type(c_ptr) :: buf !fake intent(in)
    real(f_double),pointer:: a !fake intent(out)
    integer, intent(in) :: count
    integer, intent(in), optional :: source
    integer, intent(in), optional :: tag
    integer, intent(in), optional :: comm
    integer, intent(out), optional :: request !<toggle the isend operation
    integer, dimension(MPI_STATUS_SIZE), intent(out), optional :: status !<for the blocking operation
    logical, intent(in), optional :: simulate,verbose
    integer, intent(in) :: type
    !local variables
    logical :: verb,sim
    integer :: mpi_comm,ierr,mpi_source,mpi_tag,mpi_type,tcat

    mpi_comm=MPI_COMM_WORLD
    if (present(comm)) mpi_comm=comm
    mpi_source=MPI_ANY_SOURCE
    mpi_tag=MPI_ANY_TAG
    if (present(source)) then
       mpi_source=source
       mpi_tag=source
    end if
    if (present(tag)) mpi_tag=tag
    verb=.false.
    if (present(verbose)) verb=verbose .and. source /= mpirank_null()
    if (verb) call yaml_comment('Receiving'//count//'elements from'//source//'in'//mpirank(mpi_comm))

    if (verb) then
       call yaml_mapping_open('MPI_(I)RECV')
       call yaml_map('Elements',count)
       call yaml_map('Source',source)
       call yaml_map('Dest',mpirank(mpi_comm))
       call yaml_map('Tag',mpi_tag)
       call yaml_mapping_close()
    end if

    sim=.false.
    if (present(simulate)) sim=simulate
    if (sim) return

    tcat=TCAT_RECV
    ! Synchronize the communication
    call f_timer_interrupt(tcat)

    mpi_type=type
    call c_f_pointer(buf, a)
    if (present(request)) then
       call MPI_IRECV(a,count,mpi_type,mpi_source,mpi_tag,mpi_comm,request,ierr)
    else
       if (present(status)) then
          call MPI_RECV(a,count,mpi_type,mpi_source,mpi_tag,mpi_comm,status,ierr)
       else
          call MPI_RECV(a,count,mpi_type,mpi_source,mpi_tag,mpi_comm,MPI_STATUS_IGNORE,ierr)
       end if
    end if
    call f_timer_resume()
    if (ierr/=0) call f_err_throw('An error in calling to MPI_(I)RECV occured',&
         err_id=ERR_MPI_WRAPPERS)

  end subroutine mpirecv_gpu

  subroutine mpiwaitall(ncount, array_of_requests,array_of_statuses,simulate)
    use dictionaries, only: f_err_throw,f_err_define
    implicit none
    ! Local variables
    integer, intent(in) :: ncount
    integer, dimension(ncount),intent(in) :: array_of_requests
    logical, intent(in), optional :: simulate
    integer, dimension(MPI_STATUS_SIZE,ncount), intent(out), optional :: array_of_statuses
    ! Local variables
    logical :: sim
    integer :: ierr,tcat

    !no wait if no requests
    if (ncount==0) return

    sim=.false.
    if (present(simulate)) sim=simulate
    if (sim) return

    tcat=TCAT_WAIT
    ! Synchronize the communication
    call f_timer_interrupt(tcat)
    if (present(array_of_statuses)) then
       call mpi_waitall(ncount, array_of_requests,array_of_statuses, ierr)
    else
       call mpi_waitall(ncount, array_of_requests, MPI_STATUSES_IGNORE, ierr)
    end if
    call f_timer_resume()
    if (ierr/=0) then
       call f_err_throw('An error in calling to MPI_WAITALL occured',&
            err_id=ERR_MPI_WRAPPERS)
       return
    end if

  end subroutine mpiwaitall


  subroutine mpiwait(request)
    use dictionaries, only: f_err_throw,f_err_define
    implicit none
    ! Local variables
    integer,intent(in) :: request
    ! Local variables
    integer :: ierr,tcat

    if (request /= MPI_REQUEST_NULL) then
       tcat=TCAT_WAIT
       ! Synchronize the communication
       call f_timer_interrupt(tcat)
       call mpi_wait(request, MPI_STATUSES_IGNORE, ierr)
       call f_timer_resume()
       if (ierr/=0) then
          call f_err_throw('An error in calling to MPI_WAIT occured',&
               err_id=ERR_MPI_WRAPPERS)
       end if
    end if
  end subroutine mpiwait

end module wrapper_MPI


!> Routine to gather the clocks of all the instances of flib time module
subroutine gather_timings(ndata,nproc,mpi_comm,src,dest)
  use wrapper_MPI
  implicit none
  integer, intent(in) :: ndata !< number of categories of the array
  integer, intent(in) :: nproc,mpi_comm !< number of MPI tasks and communicator
  real(kind=8), dimension(ndata), intent(in) :: src !< total timings of the instance
  real(kind=8), dimension(ndata,nproc), intent(inout) :: dest !< gathered timings
  call mpigather(sendbuf=src,recvbuf=dest,root=0,comm=mpi_comm)

end subroutine gather_timings


!> used by get_to_allgatherv to pass the good addresses to the mpiget wrapper
subroutine getall_d(nproc,recvcounts,displs,window,nrecvbuffer,recvbuffer)
  use wrapper_MPI, only: mpiget, mpi_address_kind
  implicit none
  integer,intent(in) :: nproc,nrecvbuffer,window
  integer,dimension(0:nproc-1),intent(in) :: recvcounts,displs
  double precision,dimension(nrecvbuffer),intent(out) :: recvbuffer
  ! Local variables
  integer :: jproc, jcount, jst

  do jproc=0,nproc-1
     jcount=recvcounts(jproc)
     jst=displs(jproc)
     if (jcount>0) then
         call mpiget(recvbuffer(jst+1), jcount, jproc, int(0,kind=mpi_address_kind), window)
     end if
  end do

end subroutine getall_d

!> Activates the nesting for UNBLOCK_COMMS performance case
subroutine bigdft_open_nesting(num_threads)
  use wrapper_mpi
  implicit none
  integer, intent(in) :: num_threads
#ifdef HAVE_MPI_INIT_THREAD
  !$ call OMP_SET_NESTED(.true.)
  !$ call OMP_SET_MAX_ACTIVE_LEVELS(2)
  !$ call OMP_SET_NUM_THREADS(num_threads)
#else
  integer :: idummy
  write(*,*)'BigDFT_open_nesting is not active!'
  !call MPI_ABORT(bigdft_mpi%mpi_comm,ierr)
  stop
  idummy=num_threads
#endif
end subroutine bigdft_open_nesting


!> Activates the nesting for UNBLOCK_COMMS performance case
subroutine bigdft_close_nesting(num_threads)
  use wrapper_mpi
  implicit none
  integer, intent(in) :: num_threads
#ifdef HAVE_MPI_INIT_THREAD
  !$ call OMP_SET_MAX_ACTIVE_LEVELS(1) !redundant
  !$ call OMP_SET_NESTED(.false.)
  !$ call OMP_SET_NUM_THREADS(num_threads)
#else
  integer :: idummy
  write(*,*)'BigDFT_close_nesting is not active!'
  stop
  !call MPI_ABORT(bigdft_mpi%mpi_comm,ierr)
  idummy=num_threads
#endif
end subroutine bigdft_close_nesting
