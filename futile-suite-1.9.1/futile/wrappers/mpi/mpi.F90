!> @file
!! Wrapper for the MPI call (this file is preprocessed.)
!! Use error handling
!! @author
!!    Copyright (C) 2012-2017 BigDFT group
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
  use yaml_strings !, only: operator(//)
  use f_precisions
  use f_refcnts
  use dictionaries, only: f_err_throw
  use mpif_module
  use fmpi_types!, only: ERR_MPI_WRAPPERS
  use f_sendrecv
  use f_allreduce
  use f_onesided
  use f_alltoall
  use f_allgather
  use f_bcast
  use f_gather
  implicit none

  ! MPI handling
#ifdef HAVE_MPI2
  logical, parameter :: have_mpi2 = .true.  !< Flag to use in the code to switch between MPI1 and MPI2
#else
  ! Provided now by mpif_module
  !integer :: MPI_IN_PLACE               !< Fake MPI_IN_PLACE variable to allow compilation in sumrho.
  logical, parameter :: have_mpi2 = .false. !< Flag to use in the code to switch between MPI1 and MPI2
#endif

!  include 'mpif.h'      !< MPI definitions and datatypes, now within mpif_module

  logical :: mpi_thread_funneled_is_supported=.false. !< Control the OMP_NESTED based overlap, checked by bigdft_mpi_init below

  character(len=*), parameter, public :: tgrp_mpi_name='Communications'

  !> Timing categories
  integer, public, save :: TCAT_SCATTER      = TIMING_UNINITIALIZED

  interface mpiscatter
      module procedure mpiscatter_i1i1
  end interface mpiscatter

  interface mpiscatterv
     module procedure mpiscatterv_d0
     module procedure mpiscatterv_d2d3,mpiscatterv_d3d2
  end interface mpiscatterv


!!$  interface mpiaccumulate
!!$      module procedure mpiaccumulate_double
!!$  end interface mpiaccumulate

  !> Global MPI communicator which contains all information related to the MPI process
  type, public :: mpi_environment
     !> Reference counter of the communicator.
     !! used to understand whether the communicator has to be destroyed
     type(f_reference_counter) :: refcnt
     integer :: mpi_comm=-1 !< MPI communicator
     integer :: iproc=0    !< Process Id
                         !! @ingroup RESERVED
     integer :: nproc=1    !< Number of MPI processes (in the given communicator)
                         !! @ingroup RESERVED
     integer :: igroup=0   !< MPI Group Id
     integer :: ngroup=1   !< Number of MPI groups
  end type mpi_environment

  public :: mpi_environment_null
  public :: release_mpi_environment
  public :: mpi_environment_set
  public :: mpi_environment_set1 !to be removed
  public :: mpi_environment_dict
  public :: mpi_environment_comm

  !> Fake type to enhance documentation
  type, private :: doc
     !> number of entries in buffer (integer). Useful for buffer passed by reference
     integer :: count
     !> rank of mpitask executing the operation (default value is root=0)
     integer :: root
     !> communicator of the communication
     integer :: comm
  end type doc

  private :: operator(//),f_err_throw,operator(+) ! To avoid an export from yaml_strings module


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
!    use yaml_strings!, only: yaml_toa,operator(//),f_string
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


  function mpi_environment_comm(comm) result(mpi_env)
    !create the mpi environment associated to mpi_comm_worl
    integer(fmpi_integer), intent(in), optional :: comm
    type(mpi_environment) :: mpi_env

    mpi_env=mpi_environment_null()
    mpi_env%mpi_comm=fmpi_comm(comm)
    mpi_env%nproc=mpisize(mpi_env%mpi_comm)
    mpi_env%iproc=mpirank(mpi_env%mpi_comm)

  end function mpi_environment_comm


  !> Deep copy of the mpi_environment.
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


  !> Shallow copy of the mpi_environment.
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

  subroutine mpi_environment_dict(mpi_env,dict_info)
    !build the information of the mpi environment associated and write that in a dictionary
    !also integrate the OMP information
    use dictionaries
    use dynamic_memory
    implicit none
    type(mpi_environment), intent(in) :: mpi_env
    type(dictionary), pointer :: dict_info !dictionary (valid) in which to store the information
    !local variables
    integer :: nthreads
    character(len=MPI_MAX_PROCESSOR_NAME) :: nodename_local
    character(len=MPI_MAX_PROCESSOR_NAME), dimension(:), allocatable :: nodename
    !$ integer :: omp_get_max_threads

    nthreads = 0
    !$  nthreads=omp_get_max_threads()
    call set(dict_info//'CPU parallelism'//'MPI tasks',mpi_env%nproc)
    if (nthreads /= 0) call set(dict_info//'CPU parallelism'//'OMP threads',&
         nthreads)
    nodename=f_malloc0_str(MPI_MAX_PROCESSOR_NAME,0.to.mpi_env%nproc-1,id='nodename')
    if (mpi_env%nproc>1) then
       nodename_local=mpihostname()
       call fmpi_gather(MPI_MAX_PROCESSOR_NAME,nodename_local,nodename,&
            comm=mpi_env%mpi_comm)
       if (mpi_env%iproc==0) call set(dict_info//'Hostnames',&
            list_new(.item. nodename))
    end if
    call f_free_str(MPI_MAX_PROCESSOR_NAME,nodename)
   
  end subroutine mpi_environment_dict

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


  !> Create a mpi_environment from a group list in a base group
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

!!$    if (present(comm)) then
!!$       mpi_comm=comm
!!$    else
!!$       mpi_comm=MPI_COMM_WORLD
!!$    end if
    mpi_comm=fmpi_comm(comm)
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
!      use yaml_strings, only: yaml_toa
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
    integer :: ierr,namelen,ipos

    call MPI_GET_PROCESSOR_NAME(mpihostname,namelen,ierr)
    if (ierr /= MPI_SUCCESS) then
       call f_err_throw('An error in calling to MPI_GET_PROCESSOR_NAME occured',&
            err_id=ERR_MPI_WRAPPERS)
    end if
    !clean the hostname such as to include only the last word
    !this solves a problem in ibm machines
    ipos=index(mpihostname(1:namelen),' ',back=.true.)
    if (ipos > 0) mpihostname=mpihostname(1:ipos)

  end function mpihostname

  !> Initialization of the mpi library
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


  !> Function to give MPI_COMM_WORLD
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
    use f_onesided, only: TCAT_FENCE
    use f_allreduce, only: smallsize
!    use yaml_strings, only: yaml_toa
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


  function mpireduce_i0(sendbuf,op,root,comm) result(recv)
    implicit none
    integer(f_integer), intent(in) :: sendbuf
    integer, intent(in) :: op
    integer, intent(in), optional :: comm,root
    integer(f_integer) :: recv
    !local variables
    integer :: root_,count,ierr,comm_
    count=1

    root_=0
    if (present(root)) root_=root
    comm_=mpiworld()
    if(present(comm))comm_=comm
    recv=0
    call MPI_REDUCE(sendbuf,recv,count,mpitype(sendbuf),op,comm_,ierr)
    if (ierr /=0) then
       call f_err_throw('An error in calling to MPI_REDUCE occured',&
            err_id=ERR_MPI_WRAPPERS)
       return
    end if

  end function mpireduce_i0

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

  !> create a peer_to_peer group, to use RMA calls instead of send-receive
  function p2p_group(base_grp,p1,p2,p3) result(grp)
!    use yaml_strings, only: yaml_toa
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

end module wrapper_MPI

!> Routine to gather the clocks of all the instances of flib time module
subroutine gather_timings(ndata,nproc,mpi_comm,src,dest)
  use wrapper_MPI
  implicit none
  integer, intent(in) :: ndata !< number of categories of the array
  integer, intent(in) :: nproc,mpi_comm !< number of MPI tasks and communicator
  real(kind=8), dimension(ndata), intent(in) :: src !< total timings of the instance
  real(kind=8), dimension(ndata,nproc), intent(inout) :: dest !< gathered timings
  call fmpi_gather(sendbuf=src,recvbuf=dest,root=0,comm=mpi_comm)

end subroutine gather_timings

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
