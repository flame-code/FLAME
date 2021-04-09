!> @file
!! Memory profiling module (flib)
!! @author
!!    Copyright (C) 2010-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Module which contains routines to profile the code.
!! Currently only memory occupation are provided here.
!! Ideally, timab should be incorporated here.
module memory_profiling
  implicit none

  public !<low-level module

  !Memory profiling
  type, public :: memstat
     character(len=36) :: routine,array
     integer(kind=8) :: memory,peak
  end type memstat

  type, public :: memory_state
     integer :: memalloc     !< Number of allocations recorded
     integer :: memdealloc   !< Number of deallocations recorded
     type(memstat) :: memloc !< State of the memory in the local routine
     type(memstat) :: memtot !< Global state of the memory in profiler instance
  end type memory_state

  real :: memorylimit = 0.e0 !< Limit of the memory allowed, in Gb

  private :: transform_to_MB

contains

  !> Set a memory limit
  subroutine f_set_memory_limit()
    use f_environment, only: f_memorylimit
    implicit none
    !real, intent(in) :: limit
    memorylimit = real(f_memorylimit)
  end subroutine f_set_memory_limit

  !> Retrieve the information of the present memory state
  subroutine memstate_report(memstate,dump,peak,memory,array,routine,&
       array_peak,routine_peak,&
       true_proc_peak,true_proc_memory)
    use yaml_output
    use dictionaries
    use yaml_strings, only: f_strcpy
    implicit none
    type(memory_state), intent(in) :: memstate
    !> if true, dump on default yaml stream the information about memory occupation
    logical, intent(in), optional :: dump
    integer, intent(out), optional :: peak,memory
    character(len=*), intent(out), optional :: array_peak, routine_peak,true_proc_peak,true_proc_memory
    character(len=*), intent(out), optional :: array, routine
    !local variables
    logical :: dmp
    character(len=max_field_length) :: proc_peak,proc_memory
    type(dictionary), pointer :: procstatus
    dmp=.true.
    if (present(dump)) dmp=dump
    
    !determine if procstatus has to be filled
    if (dmp .or. present(true_proc_memory) .or. &
         present(true_proc_memory)) then
       call get_proc_status_dict(procstatus)
       call f_strcpy(src='unknown',dest=proc_peak)
       call f_strcpy(src='unknown',dest=proc_memory)
       !retrieve values if they exists
       proc_peak = procstatus .get. 'VmHWM'   !Peak of memory
       proc_memory = procstatus .get. 'VmRSS' !Instant occupation memory
       call transform_to_MB(proc_peak)
       call transform_to_MB(proc_memory)
       call dict_free(procstatus)
    end if

    if (dmp) then
       call yaml_mapping_open('Memory Consumption Report')
         call yaml_map('Tot. No. of Allocations',memstate%memalloc)
         call yaml_map('Tot. No. of Deallocations',memstate%memdealloc)
         call yaml_map('Remaining Memory (B)',memstate%memtot%memory)
         call yaml_mapping_open('Memory occupation')
          call yaml_map('Peak Value (MB)',real(memstate%memtot%peak,kind=8)/1048576d0,fmt="(f15.3)")
          call yaml_map('for the array',trim(memstate%memtot%array))
          call yaml_map('in the routine',trim(memstate%memtot%routine))
          call yaml_map('Memory Peak of process',proc_peak)
         call yaml_mapping_close()
       call yaml_mapping_close()
    end if

    !retrieve separate variables, units are in MB here
    if (present(peak)) peak=int(memstate%memtot%peak/int(1048576,kind=8))
    !and kB here (the limit for local memory reporesentation is therefore 200 GB)
    if (present(memory)) memory=int(memstate%memtot%memory/int(1024,kind=8))
    if (present(array_peak)) then
       if (len_trim(memstate%memtot%array) > 0) then
          call f_strcpy(src=memstate%memtot%array,dest=array_peak)
       else
          call f_strcpy(src='Unknown',dest=array_peak)
       end if
    end if
    if (present(routine_peak)) then
       if (len_trim(memstate%memtot%routine) ==0) then
          call f_strcpy(src='Unknown',dest=routine_peak)
       else
          call f_strcpy(src=memstate%memtot%routine,dest=routine_peak)
       end if
    end if
    if (present(true_proc_peak)) call f_strcpy(src=proc_peak,dest=true_proc_peak)
    if (present(true_proc_memory)) call f_strcpy(src=proc_memory,dest=true_proc_memory)
    if (present(array)) call f_strcpy(src=memstate%memloc%array,dest=array)
    if (present(routine)) call f_strcpy(src=memstate%memloc%routine,dest=routine)

  end subroutine memstate_report

  pure subroutine transform_to_MB(data)
    use yaml_strings, only: f_strcpy
    implicit none
    character(len=*), intent(inout) :: data
    !local variables
    integer :: ipos,jpos
    character(len=256) :: tmp !< should be enough

    !search the kB unit
    ipos=index(data,'kB')
    if (ipos==0) ipos=index(data,'KB')
    if (ipos==0) return !nothing to convert in this case
    if (ipos > 0) then
       call f_strcpy(src=adjustl(data(1:ipos-1)),dest=tmp)
       !add the points before the last digits
       ipos=len_trim(tmp)
       if (ipos > 3 .and. ipos < len(tmp)) then
          !make way for the point
          do jpos=ipos,ipos-2,-1
             tmp(jpos+1:jpos+1)=tmp(jpos:jpos)
          end do
          !put the point
          tmp(ipos-2:ipos-2)='.'
          !then rewrite the data in data
          call f_strcpy(src=tmp(1:ipos+1)//' MB',dest=data)
       end if
    end if

  end subroutine transform_to_MB


  !> Dump in the unit yaml stream the status of the memory
  subroutine dump_status_line(memstate,unit,routine,array)
    use yaml_output
    use dictionaries
    use yaml_strings, only: f_strcpy,yaml_toa
    implicit none
    integer, intent(in) :: unit
    type(memory_state), intent(in) :: memstate
    character(len=*), intent(in) :: routine
    character(len=*), intent(in), optional :: array
    !local variables
    integer :: peak,memory
    character(len=40) :: localpoint
    character(len=24) :: peakpoint
    character(len=max_field_length) :: peakstr,memstr,arr,rout

    !retrieve values
    call memstate_report(memstate,dump=.false.,peak=peak,memory=memory,&
         true_proc_peak=peakstr,true_proc_memory=memstr,array_peak=arr,routine_peak=rout)

    !describe localpoint
    if (present(array)) then
       call f_strcpy(src=trim(array)//'('//trim(routine)//')',&
            dest=localpoint)
    else
       call f_strcpy(src=trim(routine),dest=localpoint)
    end if

    !describe peakpoint
    call f_strcpy(src=trim(arr)//'('//trim(rout)//')',&
         dest=peakpoint)

    call yaml_sequence(advance='no',unit=unit)
    call yaml_sequence_open(flow=.true.,unit=unit)
      call yaml_sequence(trim(localpoint),advance='no',unit=unit,padding=len(localpoint))
      call yaml_sequence(trim(yaml_toa(memory)),advance='no',unit=unit,padding=10)
      call yaml_sequence(trim(yaml_toa(peak)),advance='no',unit=unit,padding=5)
      call yaml_sequence(trim(peakpoint),advance='no',unit=unit,padding=len(peakpoint))
      call yaml_sequence(trim(memstr),advance='no',unit=unit,padding=12)
      call yaml_sequence(trim(peakstr),advance='no',unit=unit,padding=12)
    call yaml_sequence_close(unit=unit)
  end subroutine dump_status_line


  !> Test the existence of /proc/<pid>/status and read it (yaml format)
  !! for linux architecture
  subroutine get_proc_status_dict(dict)
    use yaml_parse
    use f_utils
    use yaml_strings, only: f_strcpy,yaml_toa
    use dictionaries
    implicit none
    type(dictionary), pointer :: dict
    !local variables
    logical :: exists
    integer :: pid
    character(len=64) :: filename
    type(dictionary), pointer :: dict_loaded

    nullify(dict)
    !inquire for the existence of /proc/<pid>/status
    pid=f_getpid()
    call f_strcpy(src='/proc/'//trim(adjustl(yaml_toa(pid)))//'/status',dest=filename)
    call f_file_exists(trim(filename),exists)

    if (exists) then
       !This file has a yaml structure!
       call yaml_parse_from_file(dict_loaded,trim(filename))
       dict => dict_loaded .pop. 0
       call dict_free(dict_loaded)
    end if

  end subroutine get_proc_status_dict


  !> Put to zero memocc counters
  pure subroutine memstate_init(memstate)
    use f_precisions
    implicit none
    type(memory_state), intent(out) :: memstate

    memstate%memtot%memory=f_0
    memstate%memtot%peak=f_0
    memstate%memtot%routine=f_0
    memstate%memtot%array=f_0

    memstate%memalloc=f_0
    memstate%memdealloc=f_0

    memstate%memloc%routine='routine'
    memstate%memloc%array='array'
    memstate%memloc%memory=f_0 !fake initialisation to print the first routine
    memstate%memloc%peak=f_0
  end subroutine memstate_init


  subroutine memstate_update(memstate,isize,array,routine)
    use yaml_output
    use dictionaries!error_handling
    implicit none
    ! Arguments
    type(memory_state), intent(inout) :: memstate
    integer(kind=8), intent(in) :: isize
    character(len=*), intent(in) :: array,routine
    ! Local variables
    !logical :: lmpinit
    !logical :: linq
    !integer :: ierr,istat_del,impinit
    character(len=256) :: message

    !Total counter, for all the processes
    memstate%memtot%memory=memstate%memtot%memory+isize
    if (memstate%memtot%memory > memstate%memtot%peak) then
       memstate%memtot%peak=memstate%memtot%memory
       memstate%memtot%routine=routine
       memstate%memtot%array=array
    end if
    if (isize > int(0,kind=8)) then
       memstate%memalloc=memstate%memalloc+1
    else if (isize < int(0,kind=8)) then
       memstate%memdealloc=memstate%memdealloc+1
    end if

    if (memorylimit /= 0.e0 .and. &
         memstate%memtot%memory > int(real(memorylimit,kind=8)*1073741824.d0,kind=8)) then !memory limit is in GB
       write(message,'(a,f7.3,a,i0,a)')&
            'Limit of ',memorylimit,' GB reached, total memory is ',memstate%memtot%memory,' B. '

       call f_err_throw(trim(message)//' Array '//trim(memstate%memtot%array)//&
            ', routine '//trim(memstate%memtot%routine),err_name='ERR_MEMLIMIT')
       return
    end if
    if (trim(memstate%memloc%routine) /= routine) then
       memstate%memloc%routine=routine
       memstate%memloc%array=array
       memstate%memloc%memory=isize
       memstate%memloc%peak=isize
    else
       memstate%memloc%memory=memstate%memloc%memory+isize
       if (memstate%memloc%memory > memstate%memloc%peak) then
          memstate%memloc%peak=memstate%memloc%memory
          memstate%memloc%array=array
       end if
    end if

  END SUBROUTINE memstate_update

end module memory_profiling
