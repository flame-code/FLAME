!> @file
!! Manage dynamic memory allocation
!! @author
!!    Copyright (C) 2012-2015 BigDFT group <br>
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Module used to manage memory allocations and de-allocations
module dynamic_memory_base
  use memory_profiling
  use dictionaries, info_length => max_field_length
  use yaml_strings!, only: yaml_toa,yaml_date_and_time_toa,operator(//),f_string
  use module_f_malloc 
  use f_precisions
  use yaml_parse, only: yaml_load
  use f_utils, only: f_time,f_zero
  implicit none

  private 

  logical :: track_origins=.true.      !< When true keeps track of all the allocation statuses using dictionaries
  logical :: bigdebug=.false.      !< Experimental parameter to explore the usage of f_routine as a debugger
  integer, parameter :: namelen=f_malloc_namelen  !< Length of the character variables
  integer, parameter :: error_string_len=80       !< Length of error string
  integer, parameter :: ndebug=0                  !< Size of debug parameters
!!$  integer, parameter :: max_rank=7             !< Maximum rank in fortran
  !> Maximum size of f_lib control variables
  integer, parameter :: max_ctrl = 5 !< Maximum number of nested levels
  integer :: ictrl=0                 !< Id of active control structure (<=max_ctrl)


  !> Parameters for defitions of internal dictionary
  character(len=*), parameter :: arrayid='Array Id'
  character(len=*), parameter :: routineid='Allocating Routine Id'
  character(len=*), parameter :: sizeid='Size (Bytes)'
  character(len=*), parameter :: metadatadd='Address of metadata'
  character(len=*), parameter :: firstadd='Address of first element'
  character(len=*), parameter :: processid='Process Id'
  character(len=*), parameter :: subprograms='Subroutines'
  character(len=*), parameter :: no_of_calls='No. of calls'
  character(len=*), parameter :: t0_time='Time of last opening'
  character(len=*), parameter :: tot_time='Total time (s)'
  character(len=*), parameter :: prof_enabled='Profiling Enabled'
  character(len=*), parameter :: main='Main_program'

  !> Error codes
  integer, save :: ERR_ALLOCATE
  integer, save :: ERR_DEALLOCATE
  integer, save :: ERR_MEMLIMIT
  integer, save :: ERR_INVALID_COPY
  integer, save :: ERR_MALLOC_INTERNAL

  !> Timing categories
  integer, public, save :: TCAT_ARRAY_ALLOCATIONS
  integer, public, save :: TCAT_ROUTINE_PROFILING

  !> Control structure of flib library. 
  !! Contains all global variables of interest in a separate instance of f_lib
  type :: mem_ctrl 
     logical :: profile_initialized  !< global variables for initialization
     logical :: routine_opened       !< global variable (can be stored in dictionaries)
     logical :: profile_routine      !< decide whether the routine has to be profiled
     character(len=namelen) :: present_routine !< name of the active routine 
     character(len=256) :: logfile !<file in which reports are written
     integer :: logfile_unit !< unit of the logfile stream
     integer :: output_level !< decide the level of reporting
     integer :: depth !< present level of depth of the routine
     integer :: profiling_depth !< maximum value of the profiling level
     !> Dictionaries needed for profiling storage
     type(dictionary), pointer :: dict_global    !<status of the memory at higher level
     type(dictionary), pointer :: dict_routine   !<status of the memory inside the routine
     type(dictionary), pointer :: dict_calling_sequence !<profiling of the routines
     type(dictionary), pointer :: dict_codepoint !<points to where we are in the previous dictionary
  end type mem_ctrl
  
  !> Global variable controlling the different instances of the calls
  !! the 0 component is supposed to be unused, it is allocated to avoid segfaults
  !! if the library routines are called without initialization
  type(mem_ctrl), dimension(0:max_ctrl), save :: mems
  !> Global status of the memory, used in memory_profiling module
  !! it considers the overall status of the memory regardless of the
  !! number of active instances in mems
  type(memory_state), save :: memstate 


  interface assignment(=)
     module procedure i1_all,i2_all,i3_all,i4_all
     module procedure l1_all,l2_all,l3_all
     module procedure ll1_all
     module procedure d1_all,d2_all,d3_all,d4_all,d5_all,d6_all,d7_all
     module procedure r1_all,r2_all,r3_all,r4_all
     module procedure z2_all,z3_all
     module procedure li1_all,li2_all,li3_all,li4_all
     module procedure d1_ptr,d2_ptr,d3_ptr,d4_ptr,d5_ptr,d6_ptr
     module procedure i1_ptr,i2_ptr,i3_ptr,i4_ptr
     module procedure l1_ptr, l2_ptr, l3_ptr
     module procedure li1_ptr
     module procedure z1_ptr
     !strings and pointers for characters
     module procedure c1_all
     module procedure c1_ptr
  end interface

  interface f_free
     module procedure i1_all_free,i2_all_free,i3_all_free,i4_all_free
     !     module procedure il1_all_free, il2_all_free
     module procedure i1_all_free_multi
     module procedure l1_all_free,l2_all_free,l3_all_free
     module procedure ll1_all_free
     module procedure d1_all_free,d2_all_free,d1_all_free_multi,d3_all_free,d4_all_free,d5_all_free,d6_all_free,d7_all_free
     module procedure r1_all_free,r2_all_free,r3_all_free,r4_all_free
     module procedure z2_all_free,z3_all_free
     module procedure li1_all_free,li2_all_free,li3_all_free,li4_all_free
  end interface f_free

  interface f_free_ptr
     module procedure i1_ptr_free,i2_ptr_free,i3_ptr_free,i4_ptr_free
     module procedure i1_ptr_free_multi
     module procedure d1_ptr_free,d2_ptr_free,d3_ptr_free,d4_ptr_free,d5_ptr_free,d6_ptr_free
     module procedure l1_ptr_free, l2_ptr_free, l3_ptr_free
     module procedure li1_ptr_free
     module procedure z1_ptr_free
  end interface f_free_ptr

  interface f_memcpy
     module procedure f_memcpy_i0,f_memcpy_i1,f_memcpy_i2
     module procedure f_memcpy_i0i1,f_memcpy_i1i2,f_memcpy_i2i1,f_memcpy_i2i0
     module procedure f_memcpy_li0,f_memcpy_li1
     module procedure f_memcpy_li0li1,f_memcpy_li1li2,f_memcpy_li2li1,f_memcpy_li2li0
     module procedure f_memcpy_l1
     module procedure f_memcpy_r0
     module procedure f_memcpy_d0,f_memcpy_d1,f_memcpy_d2,f_memcpy_d0d1
     module procedure f_memcpy_d1d2,f_memcpy_d2d1,f_memcpy_d2d3,f_memcpy_d3,f_memcpy_d4,f_memcpy_d1d0
     module procedure f_memcpy_d0d3,f_memcpy_d0d2,f_memcpy_d3d0,f_memcpy_d2d0,f_memcpy_d3d2
     module procedure f_memcpy_l0
     module procedure f_memcpy_c1i1,f_memcpy_i1c1,f_memcpy_c0i1
     module procedure f_memcpy_c1li1,f_memcpy_li1c1,f_memcpy_c0li1
  end interface f_memcpy

  interface f_maxdiff
     module procedure f_maxdiff_i0,f_maxdiff_i1
     module procedure f_maxdiff_li0,f_maxdiff_li1
     module procedure f_maxdiff_i0i1,f_maxdiff_i1i2,f_maxdiff_i2i1
     module procedure f_maxdiff_li0li1,f_maxdiff_li1li2,f_maxdiff_li2li1
     module procedure f_maxdiff_r0
     module procedure f_maxdiff_d0,f_maxdiff_d1,f_maxdiff_d2
     module procedure f_maxdiff_d0d1,f_maxdiff_d1d2,f_maxdiff_d2d1,f_maxdiff_d2d3
     module procedure f_maxdiff_l0
     module procedure f_maxdiff_c1i1,f_maxdiff_c0i1
     module procedure f_maxdiff_c1li1,f_maxdiff_c0li1
  end interface f_maxdiff

  interface f_subptr
     module procedure f_subptr_d0
  end interface f_subptr

  public :: f_free,f_free_ptr,f_free_str,f_free_str_ptr,f_malloc_dump_status
  public :: f_routine,f_release_routine,f_malloc_set_status,f_malloc_initialize,f_malloc_finalize
  public :: f_memcpy,f_maxdiff,f_update_database,f_purge_database,f_subptr
  public :: assignment(=),operator(.to.)

  !for internal f_lib usage
  public :: dynamic_memory_errors

contains

  pure function mem_ctrl_null() result(mem)
    type(mem_ctrl) :: mem
    call nullify_mem_ctrl(mem)
  end function mem_ctrl_null
  pure subroutine nullify_mem_ctrl(mem)
    implicit none
    type(mem_ctrl), intent(out) :: mem
    mem%profile_initialized=.false. 
    mem%routine_opened=.false.      
    mem%profile_routine=.true.
    call f_zero(mem%present_routine)!=repeat(' ',namelen)
    mem%logfile=repeat(' ',len(mem%logfile))
    mem%logfile_unit=-1 !< not initialized
    mem%output_level=0
    mem%depth=0
    mem%profiling_depth=-1 !<disabled
    !> Dictionaries needed for profiling storage
    nullify(mem%dict_global)
    nullify(mem%dict_routine)
    nullify(mem%dict_calling_sequence)
    nullify(mem%dict_codepoint)
  end subroutine nullify_mem_ctrl

  !pure 
  function mem_ctrl_init() result(mem)
    type(mem_ctrl) :: mem
    call nullify_mem_ctrl(mem)
    call initialize_mem_ctrl(mem)
  end function mem_ctrl_init
  !pure 
  subroutine initialize_mem_ctrl(mem)
    implicit none
    type(mem_ctrl), intent(inout) :: mem
    mem%profile_initialized=.true.
    !initalize the dictionary with the allocation information
    nullify(mem%dict_routine)
    call dict_init(mem%dict_global)
    call dict_init(mem%dict_calling_sequence)
    !in principle the calling sequence starts from the main
    mem%dict_codepoint => mem%dict_calling_sequence
    call set_routine_info(mem%present_routine,mem%profile_routine)
  end subroutine initialize_mem_ctrl

  !> Transfer to the f_malloc_module the information of the routine
  subroutine set_routine_info(name,profile)
    implicit none
    logical, intent(in) :: profile
    character(len=*), intent(in) :: name

    f_malloc_routine_name(1:len(f_malloc_routine_name))=name
    f_malloc_default_profiling=profile
  end subroutine set_routine_info


  !> routine to associate the rank-1 array to a workspace 
  subroutine map_workspace_d(pos,lb,lu,work,wsz,ptr)
    implicit none
    real(f_double), dimension(:), pointer :: work,wtmp,ptr

!!$    !include file
!!$    integer(f_long), intent(in) :: lb,lu
!!$    integer(f_long), intent(inout) :: pos,wsz
!!$    !local variables
!!$    integer(f_long) :: szm1
!!$    szm1=lu-lb
!!$    !check if the position and the sizes are compatible with the 
!!$    !allocation of the pointer, otherwise resize the pointer
!!$    if (wsz > pos+szm1) then
!!$       wsz=2*wsz
!!$       wtmp=>work
!!$       nullify(work)
!!$       work = f_malloc_ptr(wsz,id='work_d')
!!$       call f_memcpy(src=wtmp,dest=work)
!!$       call f_free_ptr(wtmp)
!!$    end if
!!$    call f_map_ptr(lb,lu,work(pos:pos+szm1),ptr)
!!$    !increment the position
!!$    pos=pos+szm1+1
!!$    !end of include file
    include 'f_map-inc.f90'
  end subroutine map_workspace_d

  !> routine to associate the rank-1 array to a workspace 
  subroutine map_workspace_r(pos,lb,lu,work,wsz,ptr)
    implicit none
    real(f_double), dimension(:), pointer :: work,wtmp,ptr
    include 'f_map-inc.f90'
  end subroutine map_workspace_r

  !> routine to associate the rank-1 array to a workspace 
  subroutine map_workspace_i(pos,lb,lu,work,wsz,ptr)
    implicit none
    integer(f_integer), dimension(:), pointer :: work,wtmp,ptr
    include 'f_map-inc.f90'
  end subroutine map_workspace_i

  !> routine to associate the rank-1 array to a workspace 
  subroutine map_workspace_li(pos,lb,lu,work,wsz,ptr)
    implicit none
    integer(f_long), dimension(:), pointer :: work,wtmp,ptr
    include 'f_map-inc.f90'
  end subroutine map_workspace_li

  !> routine to associate the rank-1 array to a workspace 
  subroutine map_workspace_l(pos,lb,lu,work,wsz,ptr)
    implicit none
    logical, dimension(:), pointer :: work,wtmp,ptr
    include 'f_map-inc.f90'
  end subroutine map_workspace_l

  !> Copy the contents of an array into another one
  include 'f_memcpy-inc.f90'

  subroutine set_depth(depth)
    implicit none
    integer, intent(in) :: depth
    mems(ictrl)%depth=mems(ictrl)%depth+depth
    track_origins = &
         mems(ictrl)%depth <= mems(ictrl)%profiling_depth .or. &
         mems(ictrl)%profiling_depth ==-1
  end subroutine set_depth

  !> This routine adds the corresponding subprogram name to the dictionary
  !! and prepend the dictionary to the global info dictionary
  !! if it is called more than once for the same name it has no effect
  subroutine f_routine(id,profile)
    use yaml_output, only: yaml_map,yaml_flush_document,yaml_mapping_open
    use yaml_strings, only: yaml_time_toa
    implicit none
    logical, intent(in), optional :: profile     !< ???
    character(len=*), intent(in), optional :: id !< name of the subprogram

    !local variables
    integer :: lgt,ncalls,unit_dbg
    integer(kind=8) :: itime


    if (f_err_raise(ictrl == 0,&
         'ERROR (f_routine): the routine f_malloc_initialize has not been called',&
         ERR_MALLOC_INTERNAL)) return

    if (.not. present(id)) return !no effect

    !take the time
    itime=f_time()

    !profile the profiling
    call f_timer_interrupt(TCAT_ROUTINE_PROFILING)

    call set_depth(1)

    !desactivate profile_routine if the mother routine has desactivated it
    if (present(profile)) mems(ictrl)%profile_routine=mems(ictrl)%profile_routine .and. profile

    !print *,'test',trim(mems(ictrl)%present_routine),'test2',trim(id),'ccc',mems(ictrl)%routine_opened
    !if (trim(mems(ictrl)%present_routine) /= trim(id) .or. &
    !         (trim(mems(ictrl)%present_routine) == trim(id) .and. .not. mems(ictrl)%routine_opened) ) then
    !debug
!!$    call yaml_mapping_open('Status before the opening of the routine')
!!$      call yaml_map('Level',ictrl)
!!$      call yaml_map('Routine opened',mems(ictrl)%routine_opened)
!!$      call yaml_map('Codepoint',trim(dict_key(mems(ictrl)%dict_codepoint)))
!!$      call yaml_mapping_open('Codepoint dictionary')
!!$      call yaml_dict_dump(mems(ictrl)%dict_codepoint)
!!$      call yaml_mapping_close()
!!$    call yaml_mapping_close()
    !end debug  
    !test
    if (track_origins) then !.true.) then
       if(associated(mems(ictrl)%dict_routine)) then
          !          call yaml_map('adding routine; '//trim(id),&
          !               (/dict_size(mems(ictrl)%dict_global),dict_size(mems(ictrl)%dict_routine)/))
          !          call yaml_map('Dict to add',mems(ictrl)%dict_routine)
          call prepend(mems(ictrl)%dict_global,mems(ictrl)%dict_routine)
          !          call yaml_map('verify; '//trim(id),dict_size(mems(ictrl)%dict_global))

          nullify(mems(ictrl)%dict_routine)
       end if
       !this means that the previous routine has not been closed yet
       if (mems(ictrl)%routine_opened) then
          !call open_routine(dict_codepoint)
          mems(ictrl)%dict_codepoint=>mems(ictrl)%dict_codepoint//subprograms
       end if
       mems(ictrl)%routine_opened=.true.
       !call add(dict_codepoint,trim(id))
       !see if the key existed in the codepoint
       if (has_key(mems(ictrl)%dict_codepoint,trim(id))) then
          !retrieve number of calls and increase it
          ncalls=mems(ictrl)%dict_codepoint//trim(id)//no_of_calls
          call set(mems(ictrl)%dict_codepoint//trim(id)//no_of_calls,ncalls+1)
          !write the starting point for the time
          call set(mems(ictrl)%dict_codepoint//trim(id)//t0_time,itime)
          call set(mems(ictrl)%dict_codepoint//trim(id)//prof_enabled,mems(ictrl)%profile_routine)
       else
          !create a new dictionary
          call set(mems(ictrl)%dict_codepoint//trim(id),&
               dict_new((/no_of_calls .is. yaml_toa(1), t0_time .is. yaml_toa(itime),&
               tot_time .is. yaml_toa(0.d0,fmt='(f4.1)'), &
               prof_enabled .is. yaml_toa(mems(ictrl)%profile_routine)/)))
       end if
       !then fix the new codepoint from this one
       mems(ictrl)%dict_codepoint=>mems(ictrl)%dict_codepoint//trim(id)

       lgt=min(len(id),namelen)
       call f_zero(mems(ictrl)%present_routine)
       mems(ictrl)%present_routine(1:lgt)=id(1:lgt)

    end if
    call set_routine_info(mems(ictrl)%present_routine,mems(ictrl)%profile_routine)
    if (bigdebug) then
       unit_dbg=bigdebug_stream()
       call yaml_mapping_open(mems(ictrl)%present_routine,unit=unit_dbg)
       call yaml_map('Entering',yaml_time_toa(),unit=unit_dbg)
       call yaml_flush_document(unit=unit_dbg)
    end if
    call f_timer_resume()
  end subroutine f_routine

  function bigdebug_stream() result(unit_dbg)
    use f_utils, only: f_get_free_unit
    use yaml_strings, only: operator(+)
    use yaml_output, only: yaml_stream_connected,yaml_set_stream
    implicit none
    integer :: unit_dbg
    !local variables
    integer :: istat, iproc
    character(len=32) :: strm
    !check if the stream for debugging exist
    iproc=0
    iproc= mems(ictrl)%dict_global .get. processid
    call f_strcpy(dest=strm,src='Routines-'+iproc)
    call yaml_stream_connected(strm, unit_dbg, istat)
    !otherwise create it
    if (istat /=0) then
       unit_dbg=f_get_free_unit(117)
       call yaml_set_stream(unit_dbg,filename=strm,position='rewind',setdefault=.false.)
    end if
  end function bigdebug_stream

  !> Close a previously opened routine
  subroutine f_release_routine()
    use yaml_output, only: yaml_dict_dump,yaml_map,yaml_flush_document,yaml_mapping_close,&
         yaml_comment
    use f_utils, only: f_rewind
    use yaml_strings, only: yaml_time_toa
    implicit none
    integer :: jproc,unit_dbg

    if (f_err_raise(ictrl == 0,&
         '(f_release_routine): the routine f_malloc_initialize has not been called',&
         ERR_MALLOC_INTERNAL)) return

    !profile the profiling
    call f_timer_interrupt(TCAT_ROUTINE_PROFILING)


    if (associated(mems(ictrl)%dict_routine)) then
       call prepend(mems(ictrl)%dict_global,mems(ictrl)%dict_routine)
       nullify(mems(ictrl)%dict_routine)
    end if
    !call yaml_map('Closing routine',trim(dict_key(dict_codepoint)))
    if (bigdebug) then
       unit_dbg=bigdebug_stream()
       call yaml_map('Exiting',yaml_time_toa(),unit=unit_dbg,advance='no')
       call yaml_comment(mems(ictrl)%present_routine,unit=unit_dbg)
       call yaml_mapping_close(unit=unit_dbg)
       call yaml_flush_document(unit=unit_dbg)
    end if
    !test
    if (.not. track_origins) then
       call set_depth(-1)
       call f_timer_resume()
       return
    end if
    call close_routine(mems(ictrl)%dict_codepoint,.not. mems(ictrl)%routine_opened)!trim(dict_key(dict_codepoint)))

!!$    if (f_err_check()) then
!!$       call yaml_warning('ERROR found!')
!!$       call f_dump_last_error()
!!$       call yaml_comment('End of ERROR')
!!$       call f_timer_resume()
!!$       return
!!$    end if

    !last_opened_routine=trim(dict_key(dict_codepoint))!repeat(' ',namelen)
    !the main program is opened until there is a subprograms keyword
    if (f_err_raise(.not. associated(mems(ictrl)%dict_codepoint%parent),&
         'parent not associated(A)',&
         ERR_MALLOC_INTERNAL)) then
       call f_timer_resume()
       return
    end if
    if (dict_key(mems(ictrl)%dict_codepoint%parent) == subprograms) then    
       mems(ictrl)%dict_codepoint=>mems(ictrl)%dict_codepoint%parent
       if (f_err_raise(.not. associated(mems(ictrl)%dict_codepoint%parent),&
            'parent not associated(B)',&
            ERR_MALLOC_INTERNAL)) then
          call f_timer_resume()
          return
       end if
       mems(ictrl)%dict_codepoint=>mems(ictrl)%dict_codepoint%parent
    else !back in the main program
       mems(ictrl)%routine_opened=.false.
    end if

    mems(ictrl)%present_routine(1:len(mems(ictrl)%present_routine))=&
         trim(dict_key(mems(ictrl)%dict_codepoint))
    if (.not. has_key(mems(ictrl)%dict_codepoint,prof_enabled)) then
       call yaml_dict_dump(mems(ictrl)%dict_codepoint)
       call f_err_throw('The key '//prof_enabled//' is not present in the codepoint',&
            err_id=ERR_MALLOC_INTERNAL)
       call f_timer_resume()
       return
    end if

    mems(ictrl)%profile_routine=mems(ictrl)%dict_codepoint//prof_enabled! 

    call set_routine_info(mems(ictrl)%present_routine,mems(ictrl)%profile_routine)

    !debug
!!$    call yaml_mapping_open('Codepoint after closing')
!!$    call yaml_map('Potential Reference Routine',trim(dict_key(mems(ictrl)%dict_codepoint)))
!!$    call yaml_dict_dump(mems(ictrl)%dict_codepoint)
!!$    call yaml_mapping_close()
!!$    call yaml_comment('End of release routine',hfill='=')
    !end debug

    !write the report for the output_level
    if (mems(ictrl)%output_level == 1) then
       jproc = mems(ictrl)%dict_global//processid
       if (jproc ==0) then
          !rewind unit
          call f_rewind(mems(ictrl)%logfile_unit)
          !write present value
          call dump_status_line(memstate,mems(ictrl)%logfile_unit,mems(ictrl)%present_routine)
       end if
    end if
    call set_depth(-1)
    call f_timer_resume()
  end subroutine f_release_routine

  !> Update the memory database with the data provided
  !! Use when allocating Fortran structures
  subroutine f_update_database(size,kind,rank,address,id,routine)
    use metadata_interfaces, only: long_toa
    use yaml_output, only: yaml_flush_document
    implicit none
    !> Number of elements of the buffer
    integer(kind=8), intent(in) :: size
    !> Size in bytes of one buffer element
    integer, intent(in) :: kind
    !> Rank of the array
    integer, intent(in) :: rank
    !> Address of the first buffer element.
    !! Used to store the address in the dictionary.
    !! If this argument is zero, only the memory is updated
    !! Otherwise the dictionary is created and stored
    integer(kind=8), intent(in) :: address
    !> Id of the array
    character(len=*), intent(in) :: id
    !> Id of the allocating routine
    character(len=*), intent(in) :: routine
    !local variables
    integer(kind=8) :: ilsize,jproc
    !$ include 'remove_omp-inc.f90' 

    ilsize=max(int(kind,kind=8)*size,int(0,kind=8))
    !store information only for array of size /=0
    if (track_origins .and. address /= int(0,kind=8) .and. ilsize /= int(0,kind=8)) then
       !create the dictionary array
       if (.not. associated(mems(ictrl)%dict_routine)) then
          call dict_init(mems(ictrl)%dict_routine)
       end if
       call set(mems(ictrl)%dict_routine//long_toa(address),&
            '[ '//trim(id)//', '//trim(routine)//', '//&
            trim(yaml_toa(ilsize))//', '//trim(yaml_toa(rank))//']')
    end if
    call memstate_update(memstate,ilsize,id,routine)
    if (mems(ictrl)%output_level==2) then
       jproc = mems(ictrl)%dict_global//processid
       if (jproc ==0) then
          call dump_status_line(memstate,mems(ictrl)%logfile_unit,trim(routine),trim(id))
          call yaml_flush_document(unit=mems(ictrl)%logfile_unit)
       end if
    end if
  end subroutine f_update_database


  !> Clean the database with the information of the array
  !! Use when allocating Fortran structures
  subroutine f_purge_database(size,kind,address,id,routine)
    use metadata_interfaces, only: long_toa
    use yaml_output, only: yaml_flush_document
    use yaml_strings, only: f_strcpy
    implicit none
    !> Number of elements of the buffer
    integer(kind=8), intent(in) :: size
    !> Size in bytes of one buffer element
    integer, intent(in) :: kind
    !> Address of the first buffer element.
    !! Used to store the address in the dictionary.
    !! If this argument is zero, only the memory is updated
    !! Otherwise the dictionary is created and stored
    integer(kind=8), intent(in), optional :: address
    !> Id of the array
    character(len=*), intent(in), optional :: id
    !> Id of the allocating routine
    character(len=*), intent(in), optional :: routine
    !local variables
    logical :: use_global
    integer :: jproc
    integer(kind=8) :: ilsize,jlsize,iadd
    character(len=namelen) :: array_id,routine_id
    character(len=info_length) :: array_info
    type(dictionary), pointer :: dict_add
    !$ include 'remove_omp-inc.f90' 

    iadd=int(0,kind=8)
    if (present(address)) iadd=address
    ilsize=max(int(kind,kind=8)*size,int(0,kind=8))
    !address of first element (not needed for deallocation)
    nullify(dict_add)
    if (track_origins .and. iadd/=int(0,f_address) .and. &
         ilsize/=int(0,kind=8)) then
       !hopefully only address is necessary for the deallocation

       !search in the dictionaries the address
       dict_add=>find_key(mems(ictrl)%dict_routine,long_toa(iadd))
       if (.not. associated(dict_add)) then
          dict_add=>find_key(mems(ictrl)%dict_global,long_toa(iadd))
          if (.not. associated(dict_add) .and. mems(ictrl)%profiling_depth == -1) then
             call f_err_throw('Address '//trim(long_toa(iadd))//&
                  ' not present in dictionary',ERR_INVALID_MALLOC)
             return
          else
             use_global=.true.
          end if
       else
          use_global=.false.
       end if

       if (associated(dict_add)) then
          !transform the dict_add in a list
          !retrieve the string associated to the database
          array_info=dict_add
          dict_add => yaml_load(array_info)
          !then retrieve the array information
          array_id=dict_add//0
          routine_id=dict_add//1
          jlsize=dict_add//2

          call dict_free(dict_add)

          if (ilsize /= jlsize) then
             call f_err_throw('Size of array '//trim(array_id)//&
                  ' ('//trim(yaml_toa(ilsize))//') not coherent with dictionary, found='//&
                  trim(yaml_toa(jlsize)),ERR_MALLOC_INTERNAL)
             return
          end if
          if (use_global) then
             call dict_remove(mems(ictrl)%dict_global,long_toa(iadd))
          else
             call dict_remove(mems(ictrl)%dict_routine,long_toa(iadd))
          end if
       end if
    end if
    if (.not. associated(dict_add)) then
       if (present(id)) then
          call f_strcpy(dest=array_id,src=id)
       else
          call f_strcpy(dest=array_id,src='Unknown')
       end if
       if (present(routine)) then
          call f_strcpy(dest=routine_id,src=routine)
       else
          call f_strcpy(dest=routine_id,src='Unknown')
       end if
    end if

    call memstate_update(memstate,-ilsize,trim(array_id),trim(routine_id))
    !here in the case of output_level == 2 the data can be extracted  
    if (mems(ictrl)%output_level==2) then
       jproc = mems(ictrl)%dict_global//processid
       if (jproc ==0) then
          if (len_trim(array_id) == 0) then
             if (len_trim(routine_id) == 0) then
                call dump_status_line(memstate,mems(ictrl)%logfile_unit,&
                     'Unknown','Unknown')
             else
                call dump_status_line(memstate,mems(ictrl)%logfile_unit,trim(routine_id),'Unknown')
             end if
          else if (len_trim(routine_id) == 0) then
             call dump_status_line(memstate,mems(ictrl)%logfile_unit,'Unknown',trim(array_id))
          else
             call dump_status_line(memstate,mems(ictrl)%logfile_unit,trim(routine_id),trim(array_id))
          end if
          call yaml_flush_document(unit=mems(ictrl)%logfile_unit)
       end if
    end if
  end subroutine f_purge_database


  !> Create the id of a new routine in the codepoint and points to it.
  !! works for sequences
!!$ subroutine open_routine(dict)
!!$   implicit none
!!$   type(dictionary), pointer :: dict
!!$   !local variables
!!$   integer :: ival
!!$   character(len=info_length) :: routinename
!!$   type(dictionary), pointer :: dict_tmp
!!$
!!$   !now imagine that a new routine is created
!!$   ival=dict_len(dict)-1
!!$   routinename=dict//ival
!!$
!!$   !call yaml_map('The routine which has to be converted is',trim(routinename))
!!$
!!$   call dict_remove(dict,ival)
!!$
!!$   dict_tmp=>dict//ival//trim(routinename)
!!$
!!$   dict => dict_tmp
!!$   nullify(dict_tmp)
!!$
!!$ end subroutine open_routine


  subroutine close_routine(dict,jump_up)
    !use yaml_output !debug
    implicit none
    type(dictionary), pointer :: dict
    logical, intent(in) :: jump_up
    !character(len=*), intent(in) :: name
    !local variables
    integer(kind=8) :: itime,jtime
    real(kind=8) :: rtime
    type(dictionary), pointer :: dict_tmp

    if (f_err_raise(.not. associated(dict),'routine not associated',ERR_MALLOC_INTERNAL)) return

    itime=f_time()
    !debug
!!$    call yaml_mapping_open('codepoint'//trim(dict_key(dict)))
!!$    call yaml_dict_dump(dict)
!!$    call yaml_mapping_close()
!!$    call yaml_comment('We should jump up '//trim(yaml_toa(jump_up)),hfill='}')
    !end debug

    !update the total time, if the starting point is present
    if (has_key(dict,t0_time)) then
       jtime=dict//t0_time
       jtime=itime-jtime
       rtime=dict//tot_time
       call set(dict//tot_time,rtime+real(jtime,kind=8)*1.d-9,fmt='(1pe15.7)')
       call dict_remove(dict,t0_time)
    else
       call f_err_throw('Key '//t0_time//&
            ' not found, most likely f_release_routine has been called too many times',&
            err_id=ERR_INVALID_MALLOC)
    end if

    !we should go up of three levels
    if (jump_up) then
       dict_tmp=>dict%parent
       if (f_err_raise(.not. associated(dict_tmp),'parent not associated(1)',&
         ERR_MALLOC_INTERNAL)) return
!       call yaml_map('Present Key 1',dict_key(dict_tmp))
       dict_tmp=>dict_tmp%parent
       if (f_err_raise(.not. associated(dict_tmp),'parent not associated(2)',&
            ERR_MALLOC_INTERNAL)) return
!       call yaml_map('Present Key 2',dict_key(dict_tmp))
       if (f_err_raise(.not. associated(dict_tmp%parent),'parent not associated(3)',&
            ERR_MALLOC_INTERNAL)) return
       dict_tmp=>dict_tmp%parent
       if (f_err_raise(.not. associated(dict_tmp%parent),'parent not associated(4)',&
            ERR_MALLOC_INTERNAL)) return
       dict=>dict_tmp%parent
    end if

  end subroutine close_routine

  !> Routine which is called for most of the errors of the module
  subroutine f_malloc_callback()
    use yaml_output, only: yaml_warning,yaml_flush_document
    use exception_callbacks, only: severe_callback_add 
    implicit none
    call yaml_warning('An error occured in dynamic memory module. Printing info')
    !if f_err_severe is not overridden, dump memory
    !status in the default stream
    if (severe_callback_add == 0) call f_malloc_dump_status()
    call yaml_flush_document()
    call f_err_severe()
  end subroutine f_malloc_callback

  !> Decide the error messages associated to the dynamic memory
  subroutine dynamic_memory_errors()
    use dictionaries, only: f_err_define
    implicit none
    
    call f_err_define(err_name='ERR_ALLOCATE',err_msg='Allocation error',err_id=ERR_ALLOCATE,&
         err_action='Control the order of the allocation or if the memory limit has been reached',&
         callback=f_malloc_callback)
    call f_err_define(err_name='ERR_DEALLOCATE',err_msg='Deallocation error',err_id=ERR_DEALLOCATE,&
         err_action='Control the order of the allocation or if the memory limit has been reached',&
         callback=f_malloc_callback)
    call f_err_define(err_name='ERR_MEMLIMIT',err_msg='Memory limit reached',err_id=ERR_MEMLIMIT,&
         err_action='Control the size of the arrays needed for this run with bigdft-tool program',&
         callback=f_malloc_callback)
    call f_err_define(err_name='ERR_INVALID_COPY',err_msg='Copy not allowed',&
         err_id=ERR_INVALID_COPY,&
         err_action=&
         'A f_memcpy command failed, probably invalid sizes: check sizes of arrays at runtime',&
         callback=f_malloc_callback)
    call f_err_define(err_name='ERR_INVALID_MALLOC',err_msg='Invalid specification of f_malloc',&
         err_id=ERR_INVALID_MALLOC,&
         err_action='Put coherent data for the memory space allocation',&
         callback=f_malloc_callback)
    call f_err_define(err_name='ERR_MALLOC_INTERNAL',err_msg='Internal error of memory profiler',&
         err_id=ERR_MALLOC_INTERNAL,&
         err_action='An invalid operation occurs, submit bug report to developers',&
         callback=f_malloc_callback)
  end subroutine dynamic_memory_errors

  !> Opens a new instance of the dynamic memory handling
  subroutine f_malloc_initialize()
    implicit none
    !local variables
    integer :: istat
    character(len=1) :: val
    !increase the number of active instances
    ictrl=ictrl+1
    if (f_err_raise(ictrl > max_ctrl,&
         'The number of active instances cannot exceed'//trim(yaml_toa(max_ctrl)),&
         ERR_MALLOC_INTERNAL)) return

    !extra options can be passed at the initialization
    mems(ictrl)=mem_ctrl_init()

    !in the first instance initialize the global memory info
    if (ictrl==1) then
       call memstate_init(memstate)
       !check if we are in the bigdebug mode or not
       call get_environment_variable('FUTILE_DEBUG_MODE', val,status=istat)
       bigdebug = val == '1' .and. istat == 0
    end if

    !initialize the memprofiling counters
    call set(mems(ictrl)%dict_global//'Timestamp of Profile initialization',&
         trim(yaml_date_and_time_toa()))
    !Process Id (used to dump)
    call set(mems(ictrl)%dict_global//processid,0)
    !start the profiling of the main program
    call f_routine(id=main)

    !set status of library to the initial case
    call f_malloc_set_status(memory_limit=0.e0)

  end subroutine f_malloc_initialize

  !> Initialize the library
  subroutine f_malloc_set_status(memory_limit,output_level,logfile_name,iproc,profiling_depth)
    use yaml_output!, only: yaml_date_and_time_toa
    use f_utils
    use yaml_strings
    implicit none
    !Arguments
    character(len=*), intent(in), optional :: logfile_name   !< Name of the logfile
    real(kind=4), intent(in), optional :: memory_limit       !< Memory limit
    integer, intent(in), optional :: output_level            !< Level of output for memocc
                                                             !! 0 no file, 1 light, 2 full
    integer, intent(in), optional :: iproc                   !< Process Id (used to dump, by default 0)
    !> innermost profiling level that will be applied to the code
    integer, intent(in), optional :: profiling_depth
    !local variables
    integer :: unt,jctrl,jproc

    if (f_err_raise(ictrl == 0,&
         'ERROR (f_malloc_set_status): the routine f_malloc_initialize has not been called',&
         ERR_MALLOC_INTERNAL)) return

!!$    if (.not. mems(ictrl)%profile_initialized) then
!!$       profile_initialized=.true.
!!$       !call malloc_errors()
!!$       !initalize the dictionary with the allocation information
!!$       nullify(dict_routine)
!!$       call dict_init(dict_global)
!!$       call set(dict_global//'Timestamp of Profile initialization',trim(yaml_date_and_time_toa()))
!!$       !Process Id (used to dump)
!!$       call set(dict_global//processid,0)
!!$       call dict_init(dict_calling_sequence)
!!$       !in principle the calling sequence starts from the main
!!$       dict_codepoint => dict_calling_sequence
!!$       call f_routine(id='Main program')
!!$    end if

    if (present(output_level)) then
       if (output_level > 0) then
          !first, check if we already know which proc we are
          jproc=0
          jproc=mems(ictrl)%dict_global .get. processid
          !if iproc is present, overrides
          if (present(iproc)) jproc=iproc

          if (.not. present(logfile_name)) then
               call f_err_throw('Error, f_malloc_set_status needs logfile_name for nontrivial output level',&
               err_id=ERR_INVALID_MALLOC)
               return
          end if
          !first, close the previously opened stream
          if (mems(ictrl)%logfile_unit > 0 .and. jproc==0) then
             call yaml_close_stream(unit=mems(ictrl)%logfile_unit)
          end if
          !check if it is assigned to 
          !a previous instance of malloc_set_status, and raise and exception if it is so
          do jctrl=ictrl-1,1,-1
             if (trim(logfile_name)==mems(jctrl)%logfile) &
!!$                  call f_err_throw('Logfile name "'//trim(logfile_name)//&
!!$                  '" in f_malloc_set_status invalid, aleady in use for instance No.'//&
!!$                  trim(yaml_toa(jctrl)),err_id=ERR_INVALID_MALLOC)
             call f_err_throw('Logfile name "'//trim(logfile_name)//&
                  '" in f_malloc_set_status invalid, already in use for instance No.'//jctrl&
                  ,err_id=ERR_INVALID_MALLOC)

             exit
          end do
          unt=-1 !SM: unt otherwise not defined for jproc/=0
          if (jproc == 0) then
             !eliminate the previous existing file if it has the same name
             !check if the file is opened
             call f_file_unit(trim(logfile_name),unt)
             !after this check an opened filename may now be closed
             call f_close(unt)
             !now the file can be opened
             !get a free unit, starting from 98
             unt=f_get_free_unit(98)
             call yaml_set_stream(unit=unt,filename=trim(logfile_name),position='rewind',setdefault=.false.,&
                  record_length=131)
             if (output_level==2) then
                call yaml_comment(&
                     'Present Array,Present Routine, Present Memory, Peak Memory, Peak Array, Peak Routine',&
                     unit=unt)
                call yaml_sequence_open('List of allocations',unit=unt)
             end if
          end if
          !store the found unit in the structure
          mems(ictrl)%logfile_unit=unt
          call f_strcpy(dest=mems(ictrl)%logfile,src=logfile_name)
       end if
       mems(ictrl)%output_level=output_level
    end if

    if (present(memory_limit)) call f_set_memory_limit(memory_limit)
       
    if (present(iproc)) call set(mems(ictrl)%dict_global//processid,iproc)

    if (present(profiling_depth)) then
       !set the limit of the profiling to the maximum
       !it cannot be lower than the present depth
       if (profiling_depth < mems(ictrl)%profiling_depth) then
          call f_err_throw('The profiling_depth level cannot be lowered',&
               err_id=ERR_INVALID_MALLOC)
       end if
       mems(ictrl)%profiling_depth=max(profiling_depth,mems(ictrl)%depth)
       if (profiling_depth == -1) mems(ictrl)%profiling_depth=-1 !to disab
    end if

  end subroutine f_malloc_set_status

  !> Finalize f_malloc (Display status)
  subroutine f_malloc_finalize(dump,process_id)
    use yaml_output
    use f_utils
    implicit none
    !Arguments
    logical, intent(in), optional :: dump !< Dump always information, 
                                          !! otherwise only for Process Id == 0 and errors
    integer, intent(out), optional :: process_id !< retrieve the process_id  
                                                 !! useful to dump further information after finalization
    !local variables
    integer :: pid
    logical :: dump_status
    !integer :: unt
    
    if (f_err_raise(ictrl == 0,&
         'ERROR (f_malloc_finalize): the routine f_malloc_initialize has not been called',&
         ERR_MALLOC_INTERNAL)) return

    if (present(process_id)) process_id=-1
    pid=0
    !quick return if variables not associated
    if (associated(mems(ictrl)%dict_global)) then
       !put the last values in the dictionary if not freed
       if (associated(mems(ictrl)%dict_routine)) then
          call prepend(mems(ictrl)%dict_global,mems(ictrl)%dict_routine)
          nullify(mems(ictrl)%dict_routine)
       end if
       if (present(process_id)) process_id = mems(ictrl)%dict_global//processid

       !retrieve nonetheless
       pid = mems(ictrl)%dict_global//processid

       if (present(dump)) then
          dump_status=dump
       else 
          if (pid == 0) then
             dump_status=.true.
          else
             dump_status=.false.
          end if
          !Print if error
          !if (dict_size(mems(ictrl)%dict_global) == 2) dump_status=.false.
       end if

       if (dump_status .and. dict_size(mems(ictrl)%dict_global) /= 2) then
          call yaml_warning('Heap memory has not been completely freed! See status at finalization to find where.')
          call yaml_map('Size of the global database',dict_size(mems(ictrl)%dict_global))
          !call yaml_map('Raw version',mems(ictrl)%dict_global)
          call yaml_mapping_open('Status of the memory at finalization')
          !call yaml_dict_dump(dict_global)
          call dump_leaked_memory(mems(ictrl)%dict_global)
          call yaml_mapping_close()
       end if
       call f_release_routine() !release main
       call dict_free(mems(ictrl)%dict_global)
       !    call yaml_mapping_open('Calling sequence')
       !    call yaml_dict_dump(dict_calling_sequence)
       !    call yaml_mapping_close()
       call dict_free(mems(ictrl)%dict_calling_sequence)
    end if

    if (mems(ictrl)%profile_initialized) call memstate_report(memstate,dump=dump_status)

    !close or delete report file
    if (mems(ictrl)%output_level >= 1 .and. pid==0) then
       if (mems(ictrl)%output_level == 2) call yaml_sequence_close(unit=mems(ictrl)%logfile_unit)
       call yaml_close_stream(unit=mems(ictrl)%logfile_unit)

       !in this case the file has to be removed
       !as it is used only to clarify the active codepoint
       if (mems(ictrl)%output_level == 1) call f_delete_file(mems(ictrl)%logfile)
    end if

    !nullify control structure
    mems(ictrl)=mem_ctrl_null()
    !lower the level
    ictrl=ictrl-1
    !clean the memory report (reentrant)
    if (ictrl == 0)  call memstate_init(memstate)
  end subroutine f_malloc_finalize

  !> Dump all allocations
  subroutine dump_leaked_memory(dict,unit)
    use metadata_interfaces, only: address_toi
     use yaml_output
     implicit none
     type(dictionary), pointer, intent(in) :: dict  !< dictionary containing the memory status???
     integer, intent(in), optional :: unit          !< unit to which the status should be dumped
     !Local variables
     type(dictionary), pointer :: dict_ptr
!!$     type(dictionary), pointer :: dict_list
!!$     character(len=namelen) :: array_id
!!$     character(len=info_length) :: array_info
     integer :: iunt

     if (present(unit)) then
        iunt=unit
     else
        call yaml_get_default_stream(iunt)
     end if
     dict_ptr => dict_next(dict)
     do while(associated(dict_ptr))
        !can be used if one wants more verbose information
!!$        array_info=dict_ptr
!!$        dict_list => yaml_load(array_info)
!!$        !then retrieve the array information
!!$        array_id=dict_list//0
!!$        routine_id=dict_list//1
!!$        jlsize=dict_list//2
!!$        if (has_key(dict_ptr,trim(arrayid))) then
!!$           array_id = dict_ptr//arrayid
!!$           call yaml_mapping_open(trim(array_id),unit=iunt)
!!$           call yaml_dict_dump(dict_ptr,unit=iunt)
!!$           call yaml_map(metadatadd,trim(dict_key(dict_ptr)),unit=iunt)
!!$           call yaml_mapping_close(unit=iunt)
!!$        else
        call yaml_map(trim(dict_key(dict_ptr)),dict_ptr,unit=iunt)
!!$        end if
        dict_ptr=>dict_next(dict_ptr)
     end do
  end subroutine dump_leaked_memory


  !> Dump the status of the allocated memory (and all allocations)
  subroutine f_malloc_dump_status(filename,dict_summary)
    use yaml_output
    implicit none
    character(len=*), intent(in), optional :: filename  !< file to which the memory should be dumped
    !> If present, this dictionary is filled with the summary of the 
    !! dumped dictionary. Its presence disables the normal dumping
    type(dictionary), pointer, optional, intent(out) :: dict_summary 
    !local variables
    integer, parameter :: iunit=97 !<if used switch to default
    integer :: iunt,iunit_def,istat
    type(dictionary), pointer :: dict_compact

    if (f_err_raise(ictrl == 0,&
         'ERROR (f_malloc_dump_status): the routine f_malloc_initialize has not been called',&
         ERR_MALLOC_INTERNAL)) return
    if (present(dict_summary)) then
       call dict_init(dict_summary)
       call postreatment_of_calling_sequence(-1.d0,&
            mems(ictrl)%dict_calling_sequence,dict_summary)
       !call yaml_map('Codepoint',trim(dict_key(mems(ictrl)%dict_codepoint)))
       return
    end if

    !retrieve current unit
    call yaml_get_default_stream(iunit_def)
    iunt=iunit_def
    !inquire for presence of unit iunit in the case of filename
    if (present(filename)) then
       call yaml_set_stream(unit=iunit,filename=filename,position='rewind',setdefault=.false.,&
            istat=istat)
       if (istat /= 0) then
          call yaml_newline()
          call yaml_warning('Memory allocation status filename '//trim(filename)//&
               ' not created, dumping in default stream')
          !in the case of a filename, we might investigate if an error is produced
       else
          iunt=iunit
          call f_dump_all_errors(iunit)
       end if

    end if

    call yaml_newline(unit=iunt)
!    call yaml_mapping_open('Calling sequence of Main program',unit=iunt)
!      call yaml_dict_dump(mems(ictrl)%dict_calling_sequence,unit=iunt)
!    call yaml_mapping_close(unit=iunt)
    !use the new compact version for the calling sequence
    call dict_init(dict_compact)
    call postreatment_of_calling_sequence(-1.d0,&
         mems(ictrl)%dict_calling_sequence,dict_compact)
    call yaml_map('Calling sequence of Main program (routines with * are not closed yet)',&
         dict_compact,unit=iunt)
    call dict_free(dict_compact)
    if (associated(mems(ictrl)%dict_routine)) then
       call yaml_mapping_open('Routine dictionary',unit=iunt)
       call dump_leaked_memory(mems(ictrl)%dict_routine,unit=iunt)
       call yaml_mapping_close(unit=iunt)
    end if
    call yaml_mapping_open('Global dictionary (size'//&
         trim(yaml_toa(dict_size(mems(ictrl)%dict_global)))//')',unit=iunt)
    call dump_leaked_memory(mems(ictrl)%dict_global,unit=iunt)
    call yaml_mapping_close(unit=iunt)

    !then close the file
    if (iunt /= iunit_def) then
       call yaml_close_stream(unit=iunt)
    end if

  end subroutine f_malloc_dump_status

  !>points toward a region of a given pointer
  function f_subptr_d0(ptr,region,from,size,lbound) result(win)
    implicit none
    real(f_double), dimension(:), target :: ptr
    real(f_double), dimension(:), pointer :: win
    type(array_bounds), intent(in), optional :: region
    integer, intent(in), optional :: from
    integer, intent(in), optional :: size
    integer, intent(in), optional :: lbound !<in the case of different bounds for the pointer
    !local variables
    integer(f_kind) :: lb,ub,is,ie
    
    if (present(region) .eqv. present(size)) then
       call f_err_throw('Error in f_subptr, size of the window unknown or redundant',&
            ERR_INVALID_MALLOC)
       return
    end if

    nullify(win)
    if (present(region)) then
       is=region%nlow
       ie=region%nhigh
    else if (present(from)) then
       is=from
       ie=from+size-1
    else
       is=1
       ie=size
    end if

    lb=1
    if (present(lbound)) lb=lbound
    ub=ie-is+lb
    if (ub < lb) return

    !perform the association on the window
    if (lb==1) then
       win => ptr(is:ie)
    else
       call f_map_ptr(lb,ub,ptr(is:ie),win)
    end if

    !then perform the check for the subpointer region
    if (get_lbnd(win) /= lb) call f_err_throw(&
         'ERROR (f_subptr): expected shape does not match, '//&
         trim(yaml_toa(get_lbnd(win)))//' vs. '//trim(yaml_toa(lb)),&
         ERR_MALLOC_INTERNAL)

    if (f_loc(win(lb)) /= f_loc(ptr(is)) .or. &
         f_loc(win(ub)) /= f_loc(ptr(ie))) call f_err_throw(&
         'ERROR (f_subptr): addresses do not match, the allocating system has performed a copy',&
         ERR_MALLOC_INTERNAL)

  end function f_subptr_d0

  pure function get_lbnd(win)
    implicit none
    real(f_double), dimension(:), pointer :: win
    integer :: get_lbnd

    get_lbnd=lbound(win,1)
  end function get_lbnd

  !> This routine identifies for each of the routines the most time consuming parts and print it in the logfile
  recursive subroutine postreatment_of_calling_sequence(base_time,&
       dict_cs,dict_pt)
    implicit none
    !> Time on which percentages has to be given
    double precision, intent(in) :: base_time !< needs to be explained
    type(dictionary), pointer :: dict_pt      !< needs to be explained
    type(dictionary), pointer :: dict_cs      !< needs to be explained
    !local variables
    logical :: found
    integer :: ikey,jkey,nkey,icalls,ikeystar
    integer(kind=8) :: itime,jtime
    double precision :: bt
    character(len=1) :: extra
    character(len=info_length) :: keyval,percent
    type(dictionary), pointer :: dict_tmp
    integer, dimension(:), allocatable :: ipiv
    double precision, dimension(:), allocatable :: time
    character(len=info_length), dimension(:), allocatable :: keys

    !build the dictionary of timings
    nkey=dict_size(dict_cs)
    !here f_malloc cannot be used as this routine
    !is also called in case of errors 
    if (f_err_check()) then
       allocate(time(nkey),keys(nkey),ipiv(nkey))
    else
       time=f_malloc(nkey,id='time')
       keys=f_malloc_str(info_length,nkey,id='keys')
       ipiv=f_malloc(nkey,id='ipiv')
    end if

    !now fill the timings with the values
    dict_tmp=>dict_iter(dict_cs)
    ikey=0
    ikeystar=-1
    do while (associated(dict_tmp))
       ikey=ikey+1       
       keys(ikey)=dict_key(dict_tmp)
       if (t0_time .in. dict_tmp) then
          !measure the time from last opening
          itime=f_time()
          jtime=dict_tmp//t0_time
          jtime=itime-jtime
          time(ikey)=real(jtime,kind=8)*1.d-9
          !add an asterisk to the key
          ikeystar=ikey
       else if (tot_time .in. dict_tmp) then
          time(ikey)=dict_tmp//tot_time
       else
          time(ikey)=0.d0
       end if
       dict_tmp=>dict_next(dict_tmp)
    end do

    !now order the arrays from the most to the less expensive
    call sort_positions(nkey,time,ipiv)
    do ikey=1,nkey
       jkey=ipiv(ikey)
       icalls=dict_cs//trim(keys(jkey))//no_of_calls
       !add to the dictionary the information associated to this routine
       if (jkey == ikeystar) then
          extra='*'
       else
          extra=' '
       end if

       if (base_time > 0.d0) then
!!$          call f_strcpy(src=trim(yaml_toa(time(jkey)/base_time*100.d0,fmt='(f6.2)'))//'%'//extra,&
!!$               dest=percent)
          percent(1:len(percent))=&
               trim(yaml_toa(time(jkey)/base_time*100.d0,fmt='(f6.2)'))//'%'//extra
       else
!!$          call f_strcpy(src='~'//extra,dest=percent)
          percent(1:len(percent))='~'//extra
       end if
       call add(dict_pt,dict_new(trim(keys(jkey)) .is. &
            list_new(.item. yaml_toa(time(jkey),fmt='(1pg12.3)'),&
            .item. yaml_toa(icalls), .item. percent)&
            ))
    end do
    if (f_err_check()) then
       deallocate(ipiv,keys,time)
    else
       call f_free(ipiv)
       call f_free_str(info_length,keys)
       call f_free(time)
    end if
    !now that the routine level has been ordered, inspect lower levels,

    do ikey=1,nkey
       !the first key is the name of the routine
       dict_tmp=>dict_iter(dict_pt//(ikey-1))
       keyval=dict_key(dict_tmp)
       bt=dict_tmp//0
       found = subprograms .in. dict_cs//trim(keyval)  
       if (found) then
          call postreatment_of_calling_sequence(bt,dict_cs//trim(keyval)//subprograms,&
               dict_pt//(ikey-1)//subprograms)
       end if
    end do

  end subroutine postreatment_of_calling_sequence

  !---Templates start here
  include 'malloc_templates-inc.f90'

end module dynamic_memory_base


module dynamic_memory
  use module_f_malloc
  use dynamic_memory_base
  implicit none

  public 

  private :: ERR_INVALID_MALLOC,f_malloc_namelen
end module dynamic_memory
