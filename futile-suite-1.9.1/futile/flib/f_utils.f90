!> @file
!! Manage different low-level operations
!! like operations on external files and basic operations in memory
!! @author
!!    Copyright (C) 2012-2015 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Manage low-level operations (external files and basic memory operations)
module f_utils
  use dictionaries, only: f_err_throw,f_err_define, &
       & dictionary, dict_len, dict_iter, dict_next, dict_value, max_field_length
  use yaml_strings, only: yaml_toa
  use f_precisions
  !use module_razero
  implicit none

  private

  ! This error is public
  integer, public, save :: INPUT_OUTPUT_ERROR

  integer, public, save :: TCAT_INIT_TO_ZERO

  character(len=*), parameter :: NULL_='nonechar'

  !preprocessed include file with processor-specific values
  !defines recl_kind
  include 'f_utils.inc'

  !> This type can be used to dump strings at bunches in a file
  type, public :: f_dump_buffer
     integer :: ipos
     character(len=1), dimension(:), pointer :: buf
  end type f_dump_buffer

  !>none structure, to be used for nullification
  !!together with assignment overload
  type, public :: f_none_object
     character(len=len(NULL_)) :: none
  end type f_none_object

  type, public :: f_progress_bar
     integer :: nstep !< number of steps for the progress
     integer(f_long) :: t0 !< creation time of the progress bar
     character(len=90) :: message !< Message of the progress bar to be updated
     integer :: ncall !< number of times the progress bar is called
  end type f_progress_bar

  !> Interface for difference between two intrinsic types
  interface f_diff
     module procedure f_diff_i,f_diff_r,f_diff_d,f_diff_li,f_diff_l
     module procedure f_diff_d2d3,f_diff_d2d1,f_diff_d3d1,f_diff_d1d2,f_diff_d2,f_diff_d1
     module procedure f_diff_d3
     module procedure f_diff_i2i1,f_diff_i3i1,f_diff_i1,f_diff_i2,f_diff_i1i2
     module procedure f_diff_li2li1,f_diff_li1,f_diff_li2,f_diff_li1li2
     module procedure f_diff_d0d1,f_diff_i0i1, f_diff_li0li1
     module procedure f_diff_c1i1,f_diff_c0i1
     module procedure f_diff_c1li1,f_diff_c0li1
  end interface f_diff

  interface f_sizeof
     module procedure f_sizeof_i1,f_sizeof_i2,f_sizeof_i3,f_sizeof_i4,f_sizeof_i5
     module procedure f_sizeof_li1,f_sizeof_li2,f_sizeof_li3,f_sizeof_li4,f_sizeof_li5
     module procedure f_sizeof_b1,f_sizeof_b2,f_sizeof_b3,f_sizeof_b4,f_sizeof_b5
     module procedure f_sizeof_l1,f_sizeof_l2,f_sizeof_l3,f_sizeof_l4,f_sizeof_l5
     module procedure f_sizeof_r1,f_sizeof_r2,f_sizeof_r3,f_sizeof_r4,f_sizeof_r5
     module procedure f_sizeof_d1,f_sizeof_d2,f_sizeof_d3,f_sizeof_d4,f_sizeof_d5,f_sizeof_d6,f_sizeof_d7
     module procedure f_sizeof_z1,f_sizeof_z2,f_sizeof_z3,f_sizeof_z4,f_sizeof_z5
     module procedure f_sizeof_c0,f_sizeof_c1,f_sizeof_c2
  end interface f_sizeof

  interface f_size
     module procedure f_size_i0,f_size_i1,f_size_i2,f_size_i3,f_size_i4,f_size_i5
     module procedure f_size_li1,f_size_li2,f_size_li3,f_size_li4,f_size_li5
     module procedure f_size_b1,f_size_b2,f_size_b3
     module procedure f_size_l1,f_size_l2,f_size_l3,f_size_l4,f_size_l5
     module procedure f_size_r1,f_size_r2,f_size_r3,f_size_r4,f_size_r5
     module procedure f_size_d0,f_size_d1,f_size_d2,f_size_d3,f_size_d4,f_size_d5,f_size_d6,f_size_d7
     module procedure f_size_z1,f_size_z2,f_size_z3,f_size_z4,f_size_z5
     module procedure f_size_c0,f_size_c1
  end interface f_size

  !> Initialize to zero an array (should be called f_memset)
  interface f_zero
     module procedure zero_string
     module procedure zero_li,zero_i,zero_r,zero_d,zero_l,zero_ll
     !module procedure put_to_zero_simple, put_to_zero_long
     module procedure put_to_zero_double, put_to_zero_double_1, put_to_zero_double_2
     module procedure put_to_zero_double_3, put_to_zero_double_4, put_to_zero_double_5
     module procedure put_to_zero_double_6, put_to_zero_double_7
     module procedure put_to_zero_integer,put_to_zero_integer1,put_to_zero_integer2
     module procedure put_to_zero_integer3
     module procedure put_to_zero_long,put_to_zero_long1,put_to_zero_long2
     module procedure put_to_zero_long3
     module procedure put_to_zero_r1
  end interface f_zero

  interface f_increment
     module procedure f_inc_i0
  end interface f_increment

  interface f_humantime
     module procedure f_humantime,f_ht_long
  end interface f_humantime

  interface f_assert
     module procedure f_assert, f_assert_double,f_assert_str,f_assert_byte
  end interface f_assert

  interface f_savetxt
     module procedure f_savetxt_d2
  end interface f_savetxt

  !to be verified if clock_gettime is without side-effect, otherwise the routine cannot be pure
  interface
     pure subroutine nanosec(itime)
       use f_precisions, only: f_long
       implicit none
       integer(f_long), intent(out) :: itime
     end subroutine nanosec
  end interface

  interface f_get_option
     module procedure f_get_option_l
  end interface f_get_option

  interface assignment(=)
     module procedure f_null_i0,f_null_d0,f_null_r0
     module procedure f_null_d1_ptr
     module procedure f_null_i1_ptr,f_null_i2_ptr
     module procedure f_null_l0
  end interface assignment(=)

  public :: f_diff,f_file_unit,f_mkdir,f_savetxt
  public :: f_utils_errors,f_utils_recl,f_file_exists,f_close,f_zero
  public :: f_get_free_unit,f_delete_file,f_getpid,f_rewind,f_open_file
  public :: f_increment
  public :: f_time,f_pause,f_move_file
  public :: f_progress_bar_new,update_progress_bar,f_tty,f_humantime,f_system
  public :: assignment(=),f_none,f_assert,f_sizeof,f_size,f_get_option

contains

  subroutine f_utils_errors()

    call f_err_define('INPUT_OUTPUT_ERROR',&
         'Some of intrinsic I/O fortran routines returned an error code.',&
         INPUT_OUTPUT_ERROR,&
         err_action='Check if you have correct file system permission in I/O library or check the fortran runtime library.')

  end subroutine f_utils_errors

  pure function f_none()
    implicit none
    type(f_none_object) :: f_none
    f_none%none=NULL_
  end function f_none

  pure function f_time()
    integer(f_long) :: f_time
    !local variables
    integer(f_long) :: itime
    call nanosec(itime)
    f_time=itime
  end function f_time

  subroutine f_assert(condition,id,err_id,err_name)
    use module_f_malloc, only: f_malloc_routine_name
    use yaml_strings
    use dictionaries
    implicit none
    logical, intent(in) :: condition
    character(len=*), intent(in) :: id
    integer, intent(in), optional :: err_id
    character(len=*), intent(in), optional :: err_name
    if (condition) return
    call f_err_throw('Assertion id="'+id+'" in routine="'+&
         f_malloc_routine_name+'" not satisfied. Raising error...',&
         err_id=err_id,err_name=err_name)
  end subroutine f_assert

  subroutine f_assert_byte(condition,id,err_id,err_name)
    use yaml_strings
    use dictionaries
    implicit none
    logical(f_byte), intent(in) :: condition
    character(len=*), intent(in) :: id
    integer, intent(in), optional :: err_id
    character(len=*), intent(in), optional :: err_name
    !local variables
    logical :: cond
    cond = condition
    if (cond) return
    call f_assert(cond,id,err_id,err_name)
  end subroutine f_assert_byte

  subroutine f_assert_str(condition,id,err_id,err_name)
    use yaml_strings
    use dictionaries
    implicit none
    logical, intent(in) :: condition
    type(f_string), intent(in) :: id
    integer, intent(in), optional :: err_id
    character(len=*), intent(in), optional :: err_name
    if (condition) return
    call f_assert(condition,id%msg,err_id,err_name)
  end subroutine f_assert_str

  subroutine f_assert_double(condition,id,err_id,err_name,tol)
    use module_f_malloc, only: f_malloc_namelen
    use yaml_strings
    use dictionaries
    implicit none
    real(f_double), intent(in) :: condition
    character(len=*), intent(in) :: id
    integer, intent(in), optional :: err_id
    character(len=*), intent(in), optional :: err_name
    real(f_double), intent(in), optional :: tol
    !local variables
    real(f_double) :: tol_

    tol_=1.e-12_f_double
    if (present(tol)) tol_=tol

    call f_assert(abs(condition)< tol_,id=id,err_id=err_id,err_name=err_name)

  end subroutine f_assert_double

  pure function f_get_option_l(default,opt) result(val)
    implicit none
    logical, intent(in) :: default
    logical, intent(in), optional :: opt
    logical :: val
    val=default
    if (present(opt)) val=opt
  end function f_get_option_l


  pure function f_progress_bar_new(nstep) result(bar)
    implicit none
    integer, intent(in), optional :: nstep
    type(f_progress_bar) :: bar

    bar%nstep=-1
    if (present(nstep)) bar%nstep=nstep
    bar%t0=f_time()
    call f_zero(bar%message)
    bar%ncall=0
  end function f_progress_bar_new

  pure function ticker(ncall) result(t)
    implicit none
    integer, intent(in) :: ncall
    character :: t
    select case(modulo(ncall,4))
    case(0)
       t='|'
    case(1)
       t='/'
    case(2)
       t='-'
    case(3)
       t=achar(92) !the backslash
    end select
  end function ticker

  !> Routine to build the message to be dump
  subroutine update_progress_bar(bar,istep)
    use yaml_strings
    implicit none
    integer, intent(in) :: istep
    type(f_progress_bar), intent(inout) :: bar
    !local variables
    integer, parameter :: nstars=25
    integer :: j
    real(f_double) :: percent
    real(f_double) :: time_elapsed, it_s !< in seconds
    real(f_double) :: time_remaining !< seconds, estimation

    character(len=3) :: prc
    character(len=32) ::endtime
    character(len=nstars) :: stars

    percent=real(istep,f_double)/real(bar%nstep,f_double)*100.0_f_double
    j=int(nstars*percent)/100
    write(unit=prc,fmt="(i3)")int(percent)
    if (j>0) then
       call f_strcpy(src=repeat('=',j-1)//'>',dest=stars)
    else
       call f_zero(stars)
    end if
    bar%ncall=bar%ncall+1
    time_elapsed=&
         (f_time()-bar%t0)*real(1.e-9,f_double)
    it_s=real(istep,f_double)/time_elapsed
    if (percent==0.0_f_double) then
       time_remaining=0.d0
    else
       time_remaining=time_elapsed*&
            (100.0_f_double/percent-1.0_f_double)
    end if

    if (istep <  bar%nstep) then
       call f_strcpy(src='ETA '//trim(f_humantime(time_remaining*1.e9_f_double,.true.)),&
            dest=endtime)
    else
       call f_strcpy(src='Tot '//trim(f_humantime(time_elapsed*1.e9_f_double,.true.)),dest=endtime)
    end if

    !compose the message
    call f_strcpy(src='('+yaml_time_toa()+')'//prc//&
         '% '//ticker(bar%ncall)//&
         ' ['//stars//'] ('//trim(yaml_toa(istep))//'/'//&
         trim(adjustl(yaml_toa(bar%nstep)))//&
         ', '//&
         trim(yaml_toa(it_s,fmt='(1pg12.2)'))//&
         ' it/s), '//trim(endtime),&
         dest=bar%message)

  end subroutine update_progress_bar

  pure function f_ht_long(ns,short) result(time)
    implicit none
    integer(f_long), intent(in) :: ns !<nanoseconds
    logical, intent(in), optional :: short !<if .true. only shows one units after the leading one
    character(len=95) :: time

    time=f_humantime(real(ns,f_double),short)

  end function f_ht_long

  !> Convert a time in seconds into a string of the format e.g 3.5s,10m3s,12h10m,350d12h,1y120d
  pure function f_humantime(ns,short) result(time)
    use yaml_strings
    implicit none
    real(f_double), intent(in) :: ns !<nanoseconds
    logical, intent(in), optional :: short !<if .true. only shows one units after the leading one
    character(len=95) :: time
    !local variables
    logical :: sht
    character(len=*), parameter :: fmt='(i2.2)'
    integer(f_long), parameter :: billion=int(1000000000,f_long)
    integer(f_long), parameter :: sixty=int(60,f_long)
    integer(f_long), parameter :: tsf=int(365,f_long),tf=int(24,f_long),zr=int(0,f_long)
    integer(f_long) :: nsn,m,h,d,y,si
    real(f_double) :: s

    sht=.false.
    if (present(short)) sht=short

    !get the seconds
    !s=int(ns/billion,kind=f_long)
    s=ns/real(billion,f_double)
    !then get nanosecs
    nsn=int(ns,kind=f_long)-int(s,f_long)*billion
    si=int(s,f_long)
    !then take minutes from seconds
    m=si/sixty; si=si-m*sixty
    !and hours from minutes
    h=m/sixty; m=m-h*sixty
    !days
    d=h/tf; h=h-d*tf
    !years
    y=d/tsf; d=d-y*tsf

    if (sht) then
       !find the first unit which is not zero
       if (y > zr) then
          call f_strcpy(dest=time,src=(yaml_toa(y)+'y')+(yaml_toa(d)+'d'))
       else if (d > zr) then
          call f_strcpy(dest=time,src=(yaml_toa(d)+'d')+(yaml_toa(h,fmt)+'h'))
       else if (h > zr) then
          call f_strcpy(dest=time,src=(yaml_toa(h,fmt)+'h')+(yaml_toa(m,fmt)+'m'))
       else if (m > zr) then
          call f_strcpy(dest=time,src=(yaml_toa(m,fmt)+'m')+(yaml_toa(si,fmt)+'s'))
       else
          call f_strcpy(dest=time,src=yaml_toa(real(s,f_double),'(f5.1)')+'s')
       end if
    else
       !test with new API to deal with strings
       !that would be the best solution
       call f_strcpy(dest=time,src=&
            h**fmt+':'+m**fmt+':'+si**fmt+'.'+nsn**'(i9.9)')
!!$
       !split the treatment in the case of multiple days
       if (d >0.0_f_double .or. y > 0.0_f_double ) call f_strcpy(&
            dest=time,src=(yaml_toa(y)+'y')+(yaml_toa(d)+'d')+time)
    end if

  end function f_humantime


  !>returns true if a unit is tty.
  !! essentially only check if the unit is related to stdout and check if stdout is
  !! tty
  function f_tty(unit)
    implicit none
    integer, intent(in) :: unit
    logical :: f_tty
    !local variables
    integer :: itis
    itis=0
    if (unit == 6) call stdoutistty(itis)
    f_tty=itis==1
  end function f_tty

  !>call posix sleep function for sec seconds
  subroutine f_pause(sec)
    implicit none
    integer, intent(in) :: sec !< seconds to be waited
    if (sec <=0) return
    call csleep(sec)
  end subroutine f_pause

  !> gives the maximum record length allowed for a given unit
  subroutine f_utils_recl(unt,recl_max,recl)
    implicit none
    integer, intent(in) :: unt !< unit to be checked for record length
    integer, intent(in) :: recl_max !< maximum value for record length
    !> Value for the record length. This corresponds to the minimum between recl_max and the processor-dependent value
    !! provided by inquire statement
    integer, intent(out) :: recl
    !local variables
    logical :: unit_is_open
    integer :: ierr,ierr_recl
    integer(kind=recl_kind) :: recl_file

    !in case of any error, the value is set to recl_max
    recl=recl_max
    ierr_recl=-1
    !initialize the value of recl_file
    recl_file=int(-1234567891,kind=recl_kind)
    inquire(unit=unt,opened=unit_is_open,iostat=ierr)
    if (ierr == 0 .and. .not. unit_is_open) then
       !inquire the record length for the unit
       inquire(unit=unt,recl=recl_file,iostat=ierr_recl)
    end if
    if (ierr_recl == 0) then
       recl=int(min(int(recl_max,kind=recl_kind),recl_file))
    end if
    if (recl <=0) recl=recl_max
  end subroutine f_utils_recl

  !> inquire for the existence of a file
  subroutine f_file_exists(file,exists)
    implicit none
    character(len=*), intent(in) :: file
    logical, intent(out) :: exists
    !local variables
    integer :: ierr

    exists=.false.
    inquire(file=trim(file),exist=exists,iostat=ierr)
    if (ierr /=0) then
       call f_err_throw('Error in inquiring file='//&
         trim(file)//', iostat='//trim(yaml_toa(ierr)),&
         err_id=INPUT_OUTPUT_ERROR)
    end if
    exists = exists .and. ierr==0

  end subroutine f_file_exists

  !> call the close statement and retrieve the error
  !! do not close the file if the unit is not connected
  subroutine f_close(unit)
    implicit none
    integer, intent(in) :: unit
    !local variables
    integer :: ierr

    if (unit > 0) then
       close(unit,iostat=ierr)
       if (ierr /= 0) call f_err_throw('Error in closing unit='//&
               trim(yaml_toa(unit))//', iostat='//trim(yaml_toa(ierr)),&
               err_id=INPUT_OUTPUT_ERROR)
    end if
  end subroutine f_close

  !> Search the unit associated to a filename.
  !! the unit is -1 if the file does not exists or if the file is
  !! not connected
  subroutine f_file_unit(file,unit)
    implicit none
    character(len=*), intent(in) :: file
    integer, intent(out) :: unit
    !local variables
    logical ::  exists
    integer :: unt,ierr

    unit=-1
    call f_file_exists(file,exists)
    if (exists) then
       inquire(file=trim(file),number=unt,iostat=ierr)
       if (ierr /= 0) then
          call f_err_throw('Error in inquiring file='//&
               trim(file)//' for number, iostat='//trim(yaml_toa(ierr)),&
               err_id=INPUT_OUTPUT_ERROR)
       else
          unit=unt
       end if
    end if
  end subroutine f_file_unit

  !> Get a unit which is not opened at present
  !! start the search from the unit
  function f_get_free_unit(unit) result(unt2)
    implicit none
    !> putative free unit. Starts to search from this value
    integer, intent(in), optional :: unit
    integer :: unt2
    !local variables
    logical :: unit_is_open
    integer :: unt,ierr

    unit_is_open=.true.
    unt=7
    if (present(unit)) unt=unit
    inquire(unit=unt,opened=unit_is_open,iostat=ierr)
    do while(unit_is_open .and. ierr==0)
       unt=unt+1
       inquire(unit=unt,opened=unit_is_open,iostat=ierr)
    end do
    if (ierr /=0) then
       call f_err_throw('Error in inquiring unit='//&
            trim(yaml_toa(unt))//', iostat='//trim(yaml_toa(ierr)),&
            err_id=INPUT_OUTPUT_ERROR)
    end if
    unt2=unt
  end function f_get_free_unit


  !> Create a directory from CWD path
  subroutine f_mkdir(dir,path)
    use f_precisions, only: f_integer
    implicit none
    character(len=*), intent(in) :: dir   !< Directory to be created
    character(len=*), intent(out) :: path !< Path of the created directory (trailing slash added)
    !local variables
    integer(f_integer) :: ierr
    integer(f_integer) :: lin,lout

    call f_zero(path)
    lin=int(len_trim(dir),f_integer)
    lout=int(len(path),f_integer)

    call getdir(dir,lin,path,lout,ierr)
    if (ierr /= 0 .and. ierr /= 1) then
       call f_err_throw('Error in creating directory ='//&
            trim(dir)//', iostat='//trim(yaml_toa(ierr)),&
            err_id=INPUT_OUTPUT_ERROR)
    end if

  end subroutine f_mkdir


  subroutine f_delete_file(file)
    implicit none
    character(len=*), intent(in) :: file
    !local variables
    logical :: exists
    integer :: ierr,unit
    external :: delete

    call f_file_exists(trim(file),exists)
    if (exists) then
       !close the corresponding fortran unit if the file is connected to it
       call f_file_unit(trim(file),unit)
       call f_close(unit)
       !c-function in utils.c
       call delete(trim(file),len_trim(file),ierr)
       if (ierr /=0) call f_err_throw('Error in deleting file='//&
            trim(file)//'iostat='//trim(yaml_toa(ierr)),&
            err_id=INPUT_OUTPUT_ERROR)
    end if

  end subroutine f_delete_file

  subroutine f_move_file(src,dest)
    implicit none
    character(len=*), intent(in) :: src,dest
    !local variables
    integer(f_integer) :: ierr

    call movefile(trim(src),int(len_trim(src),f_integer),&
         trim(dest),int(len_trim(dest),f_integer),ierr)
    if (ierr /= 0) call f_err_throw('Error in moving file='//&
         trim(src)//' into='//trim(dest)//', iostat='//trim(yaml_toa(ierr)),&
         err_id=INPUT_OUTPUT_ERROR)
  end subroutine f_move_file

  !> get process id
  function f_getpid()
    implicit none
    integer :: f_getpid
    !local variables
    integer :: pid
    external :: getprocid

    call getprocid(pid)
    f_getpid=pid

  end function f_getpid

  !> run the system command "command" and raise an
  !! error if needed
  subroutine f_system(command)
    implicit none
    character(len=*), intent(in) :: command
    !local variables
    integer :: ierr

    call callsystem(trim(command),len_trim(command),ierr)

    if (ierr ==-1) call f_err_throw('Error in system call "'//&
         trim(command),err_id=INPUT_OUTPUT_ERROR)

  end subroutine f_system

  !> rewind a unit
  subroutine f_rewind(unit)
    implicit none
    integer, intent(in) :: unit
    !local variables
    integer :: ierr

    rewind(unit,iostat=ierr)
    if (ierr /=0) call f_err_throw('Error in rewind unit='//&
         trim(yaml_toa(unit))//'iostat='//trim(yaml_toa(ierr)),&
         err_id=INPUT_OUTPUT_ERROR)

  end subroutine f_rewind

  !>tentative example of writing the data in a buffer
  subroutine f_write(unit,msg,advance,buffer)
    use f_precisions, only: cr => f_cr
    use yaml_strings
    !use dynamic_memory, only: f_memcpy
    implicit none
    integer, intent(in) :: unit
    character(len=*), intent(in) :: msg
    character(len=*), intent(in), optional :: advance
    type(f_dump_buffer), optional, intent(inout) :: buffer
    !local variables
    integer :: lpos
    character(len=3) :: adv
    character(len=len(cr)) :: crtmp

    adv='yes'
    if (present(advance)) call f_strcpy(src=advance,dest=adv)

    if (present(buffer)) then
       !determine the size of the input
       lpos=len(msg)
       call f_zero(crtmp)
       if (adv .eqv. 'yes') crtmp=cr
       lpos=lpos+len_trim(cr)
       !copy the values we would like to add in the buffer
       !check if the total length is bigger than buffer size
       if (lpos+buffer%ipos > size(buffer%buf)) then
          write(unit=unit,fmt='(a)') buffer%buf(:buffer%ipos)
          buffer%ipos=1
       end if
       !copy the data
       !call f_memcpy(n=lpos,src=msg+crtmp,dest=buffer%buf(buffer%ipos))
       buffer%ipos=buffer%ipos+lpos
    else
       !we should inquire if the unit is formatted or not
       write(unit=unit,fmt='(a)',advance=adv) msg
    end if
  end subroutine f_write

  !> Open a filename and retrieve the integer for the unit
  subroutine f_open_file(unit,file,status,position,action,binary)
    use yaml_strings, only: f_strcpy
    implicit none
    !> integer of the unit. On entry, it indicates the
    !! suggested unit number. On exit, it indicates the free unit
    !! which has been used for the file opening
    integer, intent(inout) :: unit
    !> filename
    character(len=*), intent(in) :: file
    !> status (unknown, old)
    character(len=*), intent(in), optional :: status
    !> position
    character(len=*), intent(in), optional :: position
    !> action (readwrite by default)
    character(len=*), intent(in), optional :: action
    !> if true, the file will be opened in the unformatted i/o
    !! if false or absent, the file will be opened for formatted i/o
    logical, intent(in), optional :: binary
    !local variables
    integer :: unt,ierror
    character(len=7) :: f_status
    character(len=11) :: f_form
    character(len=6) :: f_position
    character(len=9) :: f_action

    !first, determine if the file is already opened.
    call f_file_unit(file,unt)
    if (unt /= -1) then
       unit=unt
    else
       !find the first free unit
       unt=f_get_free_unit(unit)

       !useful open specifiers
       call f_strcpy(src='unknown',dest=f_status)
       if (present(status)) call f_strcpy(src=status,dest=f_status)

       call f_strcpy(src='formatted',dest=f_form)
       if (present(binary)) then
          if (binary) call f_strcpy(src='unformatted',dest=f_form)
       end if

       call f_strcpy(src='asis',dest=f_position)
       if (present(position)) call f_strcpy(src=position,dest=f_position)

       call f_strcpy(src='readwrite',dest=f_action)
       if (present(action)) call f_strcpy(src=action,dest=f_action)

       !then open the file with the given unit
       open(unit=unt,file=trim(file),status=f_status,form=f_form,&
            position=f_position,action=f_action,iostat=ierror)
       if (ierror /= 0) then
          call f_err_throw('Error in opening file='//&
               trim(file)//' with unit='//trim(yaml_toa(unt,fmt='(i0)'))//&
               ', iostat='//trim(yaml_toa(ierror)),&
               err_id=INPUT_OUTPUT_ERROR)
       else
          !when everything succeded, assign the unit
          unit=unt
       end if
    end if

  end subroutine f_open_file

  subroutine f_savetxt_d2(file,data)
    implicit none
    real(f_double), dimension(:,:), intent(in) :: data
    character(len=*), intent(in) :: file
    !local variables
    integer :: unt,nlines,ncols,icol,iline
    unt=90
    call f_open_file(unt,file)

    nlines=size(data,dim=2)
    ncols=size(data,dim=1)

    do iline=1,nlines
       do icol=1,ncols-1
          write(unt,'(1pg26.16e3)',advance='no')data(icol,iline)
       end do
       write(unt,'(1pg26.16e3)')data(ncols,iline)
    end do

    call f_close(unt)
  end subroutine f_savetxt_d2

  !>nullification information
  pure subroutine f_null_i0(val,nl)
    implicit none
    type(f_none_object), intent(in) :: nl
    integer, intent(out) :: val
    if (nl%none==NULL_) val=-123456789
  end subroutine f_null_i0

  !>nullification information
  pure subroutine f_null_r0(val,nl)
    implicit none
    type(f_none_object), intent(in) :: nl
    real(f_simple), intent(out) :: val
    if (nl%none==NULL_) val=-123456789.e0_f_simple
  end subroutine f_null_r0

  !>nullification information
  pure subroutine f_null_d0(val,nl)
    implicit none
    type(f_none_object), intent(in) :: nl
    real(f_double), intent(out) :: val
    if (nl%none==NULL_) val=-123456789.e0_f_double
  end subroutine f_null_d0

  pure subroutine f_null_d1_ptr(val,nl)
    implicit none
    type(f_none_object), intent(in) :: nl
    real(f_double), dimension(:), intent(out), pointer :: val
    if (nl%none==NULL_) nullify(val)
  end subroutine f_null_d1_ptr

  pure subroutine f_null_i1_ptr(val,nl)
    implicit none
    type(f_none_object), intent(in) :: nl
    integer, dimension(:), intent(out), pointer :: val
    if (nl%none==NULL_) nullify(val)
  end subroutine f_null_i1_ptr

  pure subroutine f_null_i2_ptr(val,nl)
    implicit none
    type(f_none_object), intent(in) :: nl
    integer, dimension(:,:), intent(out), pointer :: val
    if (nl%none==NULL_) nullify(val)
  end subroutine f_null_i2_ptr

  !>nullification information
  pure subroutine f_null_l0(val,nl)
    implicit none
    type(f_none_object), intent(in) :: nl
    logical, intent(out) :: val
    if (nl%none==NULL_) val=.false.
  end subroutine f_null_l0


  !>increment a integer, to be used in low-performance routines
  !to improve readability
  pure elemental subroutine f_inc_i0(i,inc)
    implicit none
    integer, intent(inout) :: i
    integer, intent(in), optional :: inc
    !local variables
    integer :: inc_

    inc_=1
    if (present(inc))inc_=inc
    i=i+inc_
  end subroutine f_inc_i0

  !>perform a difference of two objects (of similar kind)
  subroutine f_diff_i(n,a_add,b_add,diff)
    implicit none
    integer(f_long), intent(in) :: n
    integer(f_integer) :: a_add
    integer(f_integer) :: b_add
    integer(f_integer), intent(out) :: diff
    !local variables
    integer(f_long) :: idiff
    external :: diff_i
    call diff_i(n,a_add,b_add,diff,idiff)
  end subroutine f_diff_i
  subroutine f_diff_i2i1(n,a,b,diff)
    implicit none
    integer(f_long), intent(in) :: n
    integer(f_integer), dimension(:,:),   intent(in) :: a
    integer(f_integer), dimension(:), intent(in) :: b
    integer(f_integer), intent(out) :: diff
    !local variables
    integer(f_long) :: idiff
    external :: diff_i
    call diff_i(n,a(1,1),b(1),diff,idiff)
  end subroutine f_diff_i2i1
  subroutine f_diff_i3i1(n,a,b,diff)
    implicit none
    integer(f_long), intent(in) :: n
    integer(f_integer), dimension(:,:,:),   intent(in) :: a
    integer(f_integer), dimension(:), intent(in) :: b
    integer(f_integer), intent(out) :: diff
    !local variables
    integer(f_long) :: idiff
    external :: diff_i
    call diff_i(n,a(1,1,1),b(1),diff,idiff)
  end subroutine f_diff_i3i1
  subroutine f_diff_i2(n,a,b,diff)
    implicit none
    integer(f_long), intent(in) :: n
    integer(kind=4), dimension(:,:),   intent(in) :: a
    integer(kind=4), dimension(:,:), intent(in) :: b
    integer(kind=4), intent(out) :: diff
    !local variables
    integer(f_long) :: idiff
    external ::  diff_i
    call diff_i(n,a(1,1),b(1,1),diff,idiff)
  end subroutine f_diff_i2
  subroutine f_diff_i1(n,a,b,diff)
    implicit none
    integer(f_long), intent(in) :: n
    integer(kind=4), dimension(:),   intent(in) :: a
    integer(kind=4), dimension(:), intent(in) :: b
    integer(kind=4), intent(out) :: diff
    !local variables
    integer(f_long) :: idiff
    external ::  diff_i
    call diff_i(n,a(1),b(1),diff,idiff)
  end subroutine f_diff_i1
  subroutine f_diff_i1i2(n,a,b,diff)
    implicit none
    integer(f_long), intent(in) :: n
    integer(kind=4), dimension(:),   intent(in) :: a
    integer(kind=4), dimension(:,:), intent(in) :: b
    integer(kind=4), intent(out) :: diff
    !local variables
    integer(f_long) :: idiff
    external ::  diff_i
    call diff_i(n,a(1),b(1,1),diff,idiff)
  end subroutine f_diff_i1i2

  subroutine f_diff_li(n,a_add,b_add,diff)
    implicit none
    integer(f_long), intent(in) :: n
    integer(kind=8), intent(inout) :: a_add
    integer(kind=8), intent(inout) :: b_add
    integer(kind=8), intent(out) :: diff
    !local variables
    integer(f_long) :: idiff
    external ::  diff_li
    call diff_li(n,a_add,b_add,diff,idiff)
  end subroutine f_diff_li
  subroutine f_diff_li2li1(n,a,b,diff)
    implicit none
    integer(f_long), intent(in) :: n
    integer(kind=8), dimension(:,:),   intent(in) :: a
    integer(kind=8), dimension(:), intent(in) :: b
    integer(kind=8), intent(out) :: diff
    !local variables
    integer(f_long) :: idiff
    external ::  diff_li
    call diff_li(n,a(1,1),b(1),diff,idiff)
  end subroutine f_diff_li2li1
  subroutine f_diff_li2(n,a,b,diff)
    implicit none
    integer(f_long), intent(in) :: n
    integer(kind=8), dimension(:,:),   intent(in) :: a
    integer(kind=8), dimension(:,:), intent(in) :: b
    integer(kind=8), intent(out) :: diff
    !local variables
    integer(f_long) :: idiff
    external ::  diff_li
    call diff_li(n,a(1,1),b(1,1),diff,idiff)
  end subroutine f_diff_li2
  subroutine f_diff_li1(n,a,b,diff)
    implicit none
    integer(f_long), intent(in) :: n
    integer(kind=8), dimension(:),   intent(in) :: a
    integer(kind=8), dimension(:), intent(in) :: b
    integer(kind=8), intent(out) :: diff
    !local variables
    integer(f_long) :: idiff
    external ::  diff_li
    call diff_li(n,a(1),b(1),diff,idiff)
  end subroutine f_diff_li1
  subroutine f_diff_li1li2(n,a,b,diff)
    implicit none
    integer(f_long), intent(in) :: n
    integer(kind=8), dimension(:),   intent(in) :: a
    integer(kind=8), dimension(:,:), intent(in) :: b
    integer(kind=8), intent(out) :: diff
    !local variables
    integer(f_long) :: idiff
    external ::  diff_li
    call diff_li(n,a(1),b(1,1),diff,idiff)
  end subroutine f_diff_li1li2


  subroutine f_diff_r(n,a_add,b_add,diff)
    implicit none
    integer(f_long), intent(in) :: n
    real, intent(inout) :: a_add
    real, intent(inout) :: b_add
    real, intent(out) :: diff
    !local variables
    integer(f_long) :: idiff
    external ::  diff_r
    call diff_r(n,a_add,b_add,diff,idiff)
  end subroutine f_diff_r

  subroutine f_diff_d(n,a_add,b_add,diff)
    implicit none
    integer(f_long), intent(in) :: n
    double precision, intent(inout) :: a_add
    double precision, intent(inout) :: b_add
    double precision, intent(out) :: diff
    !local variables
    integer(f_long) :: idiff
    external ::  diff_d
    call diff_d(n,a_add,b_add,diff,idiff)
  end subroutine f_diff_d
  subroutine f_diff_d1(n,a,b,diff,ind)
    implicit none
    integer(f_long), intent(in) :: n
    double precision, dimension(:),   intent(in) :: a
    double precision, dimension(:), intent(in) :: b
    double precision, intent(out) :: diff
    integer(f_long), intent(out), optional :: ind
    !local variables
    integer(f_long) :: idiff
    external ::  diff_d
    call diff_d(n,a(1),b(1),diff,idiff)
    if (present(ind)) ind=idiff
  end subroutine f_diff_d1
  subroutine f_diff_d2d3(n,a,b,diff)
    implicit none
    integer(f_long), intent(in) :: n
    double precision, dimension(:,:),   intent(in) :: a
    double precision, dimension(:,:,:), intent(in) :: b
    double precision, intent(out) :: diff
    !local variables
    integer(f_long) :: idiff
    external ::  diff_d
    call diff_d(n,a(1,1),b(1,1,1),diff,idiff)
  end subroutine f_diff_d2d3
  subroutine f_diff_d2d1(n,a,b,diff)
    implicit none
    integer(f_long), intent(in) :: n
    double precision, dimension(:,:),   intent(in) :: a
    double precision, dimension(:), intent(in) :: b
    double precision, intent(out) :: diff
    !local variables
    integer(f_long) :: idiff
    external ::  diff_d
    call diff_d(n,a(1,1),b(1),diff,idiff)
  end subroutine f_diff_d2d1
  subroutine f_diff_d3d1(n,a,b,diff)
    implicit none
    integer(f_long), intent(in) :: n
    double precision, dimension(:,:,:),   intent(in) :: a
    double precision, dimension(:), intent(in) :: b
    double precision, intent(out) :: diff
    !local variables
    integer(f_long) :: idiff
    external ::  diff_d
    call diff_d(n,a(1,1,1),b(1),diff,idiff)
  end subroutine f_diff_d3d1
  subroutine f_diff_d0d1(n,a,b,diff)
    implicit none
    integer(f_long), intent(in) :: n
    double precision, intent(inout) :: a
    double precision, dimension(:), intent(in) :: b
    double precision, intent(out) :: diff
    !local variables
    integer(f_long) :: idiff
    external ::  diff_d
    call diff_d(n,a,b(1),diff,idiff)
  end subroutine f_diff_d0d1

  subroutine f_diff_d2(n,a,b,diff)
    implicit none
    integer(f_long), intent(in) :: n
    double precision, dimension(:,:),   intent(in) :: a
    double precision, dimension(:,:), intent(in) :: b
    double precision, intent(out) :: diff
    !local variables
    integer(f_long) :: idiff
    external ::  diff_d
    call diff_d(n,a(1,1),b(1,1),diff,idiff)
  end subroutine f_diff_d2

  subroutine f_diff_d3(n,a,b,diff)
    implicit none
    integer(f_long), intent(in) :: n
    real(f_double), dimension(:,:,:),   intent(in) :: a
    real(f_double), dimension(:,:,:), intent(in) :: b
    real(f_double), intent(out) :: diff
    !local variables
    integer(f_long) :: idiff
    external ::  diff_d
    call diff_d(n,a(1,1,1),b(1,1,1),diff,idiff)
  end subroutine f_diff_d3


  subroutine f_diff_d1d2(n,a,b,diff)
    implicit none
    integer(f_long), intent(in) :: n
    double precision, dimension(:),   intent(in) :: a
    double precision, dimension(:,:), intent(in) :: b
    double precision, intent(out) :: diff
    !local variables
    integer(f_long) :: idiff
    external ::  diff_d
    call diff_d(n,a(1),b(1,1),diff,idiff)
  end subroutine f_diff_d1d2



  subroutine f_diff_l(n,a_add,b_add,diff)
    implicit none
    integer(f_long), intent(in) :: n
    logical, intent(inout) :: a_add
    logical, intent(inout) :: b_add
    logical, intent(out) :: diff
    external :: diff_l
    call diff_l(n,a_add,b_add,diff)
  end subroutine f_diff_l

  subroutine f_diff_c1i1(n,a,b,diff)
    implicit none
    integer(f_long), intent(in) :: n
    character, dimension(:),   intent(in) :: a
    integer(f_integer), dimension(:), intent(in) :: b
    integer(f_integer), intent(out) :: diff
    !local variables
    integer(f_long) :: idiff
    external ::  diff_ci
    call diff_ci(n,a(1),b(1),diff,idiff)
  end subroutine f_diff_c1i1

  subroutine f_diff_c1li1(n,a,b,diff)
    implicit none
    integer(f_long), intent(in) :: n
    character, dimension(:),   intent(in) :: a
    integer(f_long), dimension(:), intent(in) :: b
    integer(f_long), intent(out) :: diff
    !local variables
    integer(f_long) :: idiff
    external ::  diff_ci
    call diff_ci(n,a(1),b(1),diff,idiff)
  end subroutine f_diff_c1li1

  subroutine f_diff_c0i1(n,a,b,diff)
    implicit none
    integer(f_long), intent(in) :: n
    character(len=*),   intent(in) :: a
    integer(f_integer), dimension(:), intent(in) :: b
    integer(f_integer), intent(out) :: diff
    !local variables
    integer(f_long) :: idiff
    external ::  diff_ci
    call diff_ci(n,a,b(1),diff,idiff)
  end subroutine f_diff_c0i1

  subroutine f_diff_c0li1(n,a,b,diff)
    implicit none
    integer(f_long), intent(in) :: n
    character(len=*),   intent(in) :: a
    integer(f_long), dimension(:), intent(in) :: b
    integer(f_long), intent(out) :: diff
    !local variables
    integer(f_long) :: idiff
    external ::  diff_ci
    call diff_ci(n,a,b(1),diff,idiff)
  end subroutine f_diff_c0li1

  subroutine f_diff_li0li1(n,a,b,diff)
    implicit none
    integer(f_long), intent(in) :: n
    integer(f_long), intent(inout) :: a
    integer(f_long), dimension(:), intent(in) :: b
    integer(f_long), intent(out) :: diff
    !local variables
    integer(f_long) :: idiff
    external ::  diff_li
    call diff_li(n,a,b(1),diff,idiff)
  end subroutine f_diff_li0li1

  subroutine f_diff_i0i1(n,a,b,diff)
    implicit none
    integer(f_long), intent(in) :: n
    integer(f_integer), intent(inout) :: a
    integer(f_integer), dimension(:), intent(in) :: b
    integer(f_integer), intent(out) :: diff
    !local variables
    integer(f_long) :: idiff
    external ::  diff_i
    call diff_i(n,a,b(1),diff,idiff)
  end subroutine f_diff_i0i1

  pure subroutine zero_string(str)
    use yaml_strings, only: f_strcpy
    implicit none
    character(len=*), intent(out) :: str
    call f_strcpy(src=' ',dest=str)
  end subroutine zero_string

  pure subroutine zero_li(val)
    implicit none
    integer(f_long), intent(out) :: val
    val=int(0,f_long)
  end subroutine zero_li

  pure subroutine zero_i(val)
    implicit none
    integer(f_integer), intent(out) :: val
    val=0
  end subroutine zero_i

  pure subroutine zero_r(val)
    implicit none
    real, intent(out) :: val
    val=0.e0
  end subroutine zero_r

  pure subroutine zero_d(val)
    implicit none
    double precision, intent(out) :: val
    val=0.d0
  end subroutine zero_d

  pure subroutine zero_l(val)
    implicit none
    logical, intent(out) :: val
    val=.false.
  end subroutine zero_l

  pure subroutine zero_ll(val)
    implicit none
    logical(f_byte), intent(out) :: val
    val=f_F
  end subroutine zero_ll

!!$  subroutine put_to_zero_byte_1(da)
!!$    implicit none
!!$    logical(f_byte), dimension(:), intent(out) :: da
!!$    call f_timer_interrupt(TCAT_INIT_TO_ZERO)
!!$    !call razero(size(da),da(lbound(da,1)))
!!$    call setzero(int(size(da),f_long)*kind(da),da)
!!$    call f_timer_resume()
!!$  end subroutine put_to_zero_byte_1


  subroutine put_to_zero_simple(n,da)
    implicit none
    integer, intent(in) :: n
    real :: da

    call f_timer_interrupt(TCAT_INIT_TO_ZERO)
    !call razero_simple(n,da)
    call setzero(int(n,f_long)*kind(da),da)
    call f_timer_resume()
  end subroutine put_to_zero_simple

  subroutine put_to_zero_r1(da)
    implicit none
    real(f_simple), dimension(:), intent(out) :: da
    call f_timer_interrupt(TCAT_INIT_TO_ZERO)
    !call razero(size(da),da(lbound(da,1)))
    call setzero(int(size(da),f_long)*kind(da),da)
    call f_timer_resume()
  end subroutine put_to_zero_r1


  subroutine put_to_zero_double(n,da)
    implicit none
    integer, intent(in) :: n
    double precision, intent(out) :: da
    call f_timer_interrupt(TCAT_INIT_TO_ZERO)
    !call razero(n,da)
    call setzero(int(n,f_long)*kind(da),da)
    call f_timer_resume()
  end subroutine put_to_zero_double

  subroutine put_to_zero_double_1(da)
    implicit none
    double precision, dimension(:), intent(out) :: da
    call f_timer_interrupt(TCAT_INIT_TO_ZERO)
    !call razero(size(da),da(lbound(da,1)))
    call setzero(int(size(da),f_long)*kind(da),da)
    call f_timer_resume()
  end subroutine put_to_zero_double_1

  subroutine put_to_zero_double_2(da)
    implicit none
    double precision, dimension(:,:), intent(out) :: da
    call f_timer_interrupt(TCAT_INIT_TO_ZERO)
    !call razero(size(da),da(lbound(da,1),lbound(da,2)))
    call setzero(int(size(da),f_long)*kind(da),da)
    call f_timer_resume()
  end subroutine put_to_zero_double_2

  subroutine put_to_zero_double_3(da)
    implicit none
    double precision, dimension(:,:,:), intent(out) :: da
    call f_timer_interrupt(TCAT_INIT_TO_ZERO)
    !call razero(size(da),da(lbound(da,1),lbound(da,2),lbound(da,3)))
    call setzero(int(size(da),f_long)*kind(da),da)
    call f_timer_resume()
  end subroutine put_to_zero_double_3

  subroutine put_to_zero_double_4(da)
    implicit none
    double precision, dimension(:,:,:,:), intent(out) :: da
    call f_timer_interrupt(TCAT_INIT_TO_ZERO)
    !call razero(size(da),da(lbound(da,1),lbound(da,2),lbound(da,3),lbound(da,4)))
    call setzero(int(size(da),f_long)*kind(da),da)
    call f_timer_resume()
  end subroutine put_to_zero_double_4

  subroutine put_to_zero_double_5(da)
    implicit none
    double precision, dimension(:,:,:,:,:), intent(out) :: da
    call f_timer_interrupt(TCAT_INIT_TO_ZERO)
    !call razero(size(da),da(lbound(da,1),lbound(da,2),lbound(da,3),lbound(da,4),lbound(da,5)))
    call setzero(int(size(da),f_long)*kind(da),da)
    call f_timer_resume()
  end subroutine put_to_zero_double_5

  subroutine put_to_zero_double_6(da)
    implicit none
    double precision, dimension(:,:,:,:,:,:), intent(out) :: da
    call f_timer_interrupt(TCAT_INIT_TO_ZERO)
    !call razero(size(da),da(lbound(da,1),lbound(da,2),lbound(da,3),lbound(da,4),lbound(da,5),lbound(da,6)))
    call setzero(int(size(da),f_long)*kind(da),da)
    call f_timer_resume()
  end subroutine put_to_zero_double_6

  subroutine put_to_zero_double_7(da)
    implicit none
    double precision, dimension(:,:,:,:,:,:,:), intent(out) :: da
    call f_timer_interrupt(TCAT_INIT_TO_ZERO)
    !call razero(size(da),da(lbound(da,1),lbound(da,2),lbound(da,3),lbound(da,4),lbound(da,5),lbound(da,6),lbound(da,7)))
    call setzero(int(size(da),f_long)*kind(da),da)
    call f_timer_resume()
  end subroutine put_to_zero_double_7

  subroutine put_to_zero_integer(n,da)
    implicit none
    integer, intent(in) :: n
    integer(f_integer) :: da
    call f_timer_interrupt(TCAT_INIT_TO_ZERO)
    !call razero_integer(n,da)
    call setzero(int(n,f_long)*kind(da),da)
    call f_timer_resume()
  end subroutine put_to_zero_integer

  subroutine put_to_zero_integer1(da)
    implicit none
    integer(f_integer), dimension(:), intent(out) :: da
    call f_timer_interrupt(TCAT_INIT_TO_ZERO)
    !call razero_integer(size(da),da(lbound(da,1)))
    call setzero(int(size(da),f_long)*kind(da),da)
    call f_timer_resume()
  end subroutine put_to_zero_integer1

  subroutine put_to_zero_integer2(da)
    implicit none
    integer(f_integer), dimension(:,:), intent(out) :: da
    call f_timer_interrupt(TCAT_INIT_TO_ZERO)
    !call razero_integer(size(da),da(lbound(da,1),lbound(da,2)))
    call setzero(int(size(da),f_long)*kind(da),da)
    call f_timer_resume()
  end subroutine put_to_zero_integer2

  subroutine put_to_zero_integer3(da)
    implicit none
    integer(f_integer), dimension(:,:,:), intent(out) :: da
    call f_timer_interrupt(TCAT_INIT_TO_ZERO)
    !call razero_integer(size(da),da(lbound(da,1),lbound(da,2),lbound(da,3)))
    call setzero(int(size(da),f_long)*kind(da),da)
    call f_timer_resume()
  end subroutine put_to_zero_integer3


  subroutine put_to_zero_long(n,da)
    implicit none
    integer, intent(in) :: n
    integer(f_long) :: da
    call f_timer_interrupt(TCAT_INIT_TO_ZERO)
    !call razero_integerlong(n,da)
    call setzero(int(n,f_long)*kind(da),da)
    call f_timer_resume()
  end subroutine put_to_zero_long

  subroutine put_to_zero_long1(da)
    implicit none
    integer(f_long), dimension(:), intent(out) :: da
    call f_timer_interrupt(TCAT_INIT_TO_ZERO)
    !call razero_integerlong(size(da),da(lbound(da,1)))
    call setzero(int(size(da),f_long)*kind(da),da)
    call f_timer_resume()
  end subroutine put_to_zero_long1

  subroutine put_to_zero_long2(da)
    implicit none
    integer(f_long), dimension(:,:), intent(out) :: da
    call f_timer_interrupt(TCAT_INIT_TO_ZERO)
    !call razero_integerlong(size(da),da(lbound(da,1),lbound(da,2)))
    call setzero(int(size(da),f_long)*kind(da),da)
    call f_timer_resume()
  end subroutine put_to_zero_long2

  subroutine put_to_zero_long3(da)
    implicit none
    integer(f_long), dimension(:,:,:), intent(out) :: da
    call f_timer_interrupt(TCAT_INIT_TO_ZERO)
    !call razero_integerlong(size(da),da(lbound(da,1),lbound(da,2),lbound(da,3)))
    call setzero(int(size(da),f_long)*kind(da),da)
    call f_timer_resume()
  end subroutine put_to_zero_long3

  pure function f_sizeof_i1(datatype) result(s)
    integer(f_integer), dimension(:), intent(in) :: datatype
    integer(f_long) :: s; s=product(int(shape(datatype),f_long))*int(kind(datatype),f_long)
  end function f_sizeof_i1
  pure function f_sizeof_i2(datatype) result(s)
    integer(f_integer), dimension(:,:), intent(in) :: datatype
    integer(f_long) :: s; s=product(int(shape(datatype),f_long))*int(kind(datatype),f_long)
  end function f_sizeof_i2
  pure function f_sizeof_i3(datatype) result(s)
    integer(f_integer), dimension(:,:,:), intent(in) :: datatype
    integer(f_long) :: s; s=product(int(shape(datatype),f_long))*int(kind(datatype),f_long)
  end function f_sizeof_i3
  pure function f_sizeof_i4(datatype) result(s)
    integer(f_integer), dimension(:,:,:,:), intent(in) :: datatype
    integer(f_long) :: s; s=product(int(shape(datatype),f_long))*int(kind(datatype),f_long)
  end function f_sizeof_i4
  pure function f_sizeof_i5(datatype) result(s)
    integer(f_integer), dimension(:,:,:,:,:), intent(in) :: datatype
    integer(f_long) :: s; s=product(int(shape(datatype),f_long))*int(kind(datatype),f_long)
  end function f_sizeof_i5


  pure function f_sizeof_li1(datatype) result(s)
    integer(f_long), dimension(:), intent(in) :: datatype
    integer(f_long) :: s; s=product(int(shape(datatype),f_long))*int(kind(datatype),f_long)
  end function f_sizeof_li1
  pure function f_sizeof_li2(datatype) result(s)
    integer(f_long), dimension(:,:), intent(in) :: datatype
    integer(f_long) :: s; s=product(int(shape(datatype),f_long))*int(kind(datatype),f_long)
  end function f_sizeof_li2
  pure function f_sizeof_li3(datatype) result(s)
    integer(f_long), dimension(:,:,:), intent(in) :: datatype
    integer(f_long) :: s; s=product(int(shape(datatype),f_long))*int(kind(datatype),f_long)
  end function f_sizeof_li3
  pure function f_sizeof_li4(datatype) result(s)
    integer(f_long), dimension(:,:,:,:), intent(in) :: datatype
    integer(f_long) :: s; s=product(int(shape(datatype),f_long))*int(kind(datatype),f_long)
  end function f_sizeof_li4
  pure function f_sizeof_li5(datatype) result(s)
    integer(f_long), dimension(:,:,:,:,:), intent(in) :: datatype
    integer(f_long) :: s; s=product(int(shape(datatype),f_long))*int(kind(datatype),f_long)
  end function f_sizeof_li5


  pure function f_sizeof_l1(datatype) result(s)
    logical, dimension(:), intent(in) :: datatype
    integer(f_long) :: s; s=product(int(shape(datatype),f_long))*int(kind(datatype),f_long)
  end function f_sizeof_l1
  pure function f_sizeof_l2(datatype) result(s)
    logical, dimension(:,:), intent(in) :: datatype
    integer(f_long) :: s; s=product(int(shape(datatype),f_long))*int(kind(datatype),f_long)
  end function f_sizeof_l2
  pure function f_sizeof_l3(datatype) result(s)
    logical, dimension(:,:,:), intent(in) :: datatype
    integer(f_long) :: s; s=product(int(shape(datatype),f_long))*int(kind(datatype),f_long)
  end function f_sizeof_l3
  pure function f_sizeof_l4(datatype) result(s)
    logical, dimension(:,:,:,:), intent(in) :: datatype
    integer(f_long) :: s; s=product(int(shape(datatype),f_long))*int(kind(datatype),f_long)
  end function f_sizeof_l4
  pure function f_sizeof_l5(datatype) result(s)
    logical, dimension(:,:,:,:,:), intent(in) :: datatype
    integer(f_long) :: s; s=product(int(shape(datatype),f_long))*int(kind(datatype),f_long)
  end function f_sizeof_l5

  pure function f_sizeof_b1(datatype) result(s)
    logical(f_byte), dimension(:), intent(in) :: datatype
    integer(f_long) :: s; s=product(int(shape(datatype),f_long))*int(kind(datatype),f_long)
  end function f_sizeof_b1
  pure function f_sizeof_b2(datatype) result(s)
    logical(f_byte), dimension(:,:), intent(in) :: datatype
    integer(f_long) :: s; s=product(int(shape(datatype),f_long))*int(kind(datatype),f_long)
  end function f_sizeof_b2
  pure function f_sizeof_b3(datatype) result(s)
    logical(f_byte), dimension(:,:,:), intent(in) :: datatype
    integer(f_long) :: s; s=product(int(shape(datatype),f_long))*int(kind(datatype),f_long)
  end function f_sizeof_b3
  pure function f_sizeof_b4(datatype) result(s)
    logical(f_byte), dimension(:,:,:,:), intent(in) :: datatype
    integer(f_long) :: s; s=product(int(shape(datatype),f_long))*int(kind(datatype),f_long)
  end function f_sizeof_b4
  pure function f_sizeof_b5(datatype) result(s)
    logical(f_byte), dimension(:,:,:,:,:), intent(in) :: datatype
    integer(f_long) :: s; s=product(int(shape(datatype),f_long))*int(kind(datatype),f_long)
  end function f_sizeof_b5

  pure function f_sizeof_d1(datatype) result(s)
    real(f_double), dimension(:), intent(in) :: datatype
    integer(f_long) :: s; s=product(int(shape(datatype),f_long))*int(kind(datatype),f_long)
  end function f_sizeof_d1

  pure function f_sizeof_d2(datatype) result(s)
    real(f_double), dimension(:,:), intent(in) :: datatype
    integer(f_long) :: s; s=product(int(shape(datatype),f_long))*int(kind(datatype),f_long)
  end function f_sizeof_d2

  pure function f_sizeof_d3(datatype) result(s)
    real(f_double), dimension(:,:,:), intent(in) :: datatype
    integer(f_long) :: s; s=product(int(shape(datatype),f_long))*int(kind(datatype),f_long)
  end function f_sizeof_d3

  pure function f_sizeof_d4(datatype) result(s)
    real(f_double), dimension(:,:,:,:), intent(in) :: datatype
    integer(f_long) :: s; s=product(int(shape(datatype),f_long))*int(kind(datatype),f_long)
  end function f_sizeof_d4

  pure function f_sizeof_d5(datatype) result(s)
    real(f_double), dimension(:,:,:,:,:), intent(in) :: datatype
    integer(f_long) :: s; s=product(int(shape(datatype),f_long))*int(kind(datatype),f_long)
  end function f_sizeof_d5

  pure function f_sizeof_d6(datatype) result(s)
    real(f_double), dimension(:,:,:,:,:,:), intent(in) :: datatype
    integer(f_long) :: s; s=product(int(shape(datatype),f_long))*int(kind(datatype),f_long)
  end function f_sizeof_d6

  pure function f_sizeof_d7(datatype) result(s)
    real(f_double), dimension(:,:,:,:,:,:,:), intent(in) :: datatype
    integer(f_long) :: s; s=product(int(shape(datatype),f_long))*int(kind(datatype),f_long)
  end function f_sizeof_d7


  pure function f_sizeof_r1(datatype) result(s)
    real(f_simple), dimension(:), intent(in) :: datatype
    integer(f_long) :: s; s=product(int(shape(datatype),f_long))*int(kind(datatype),f_long)
  end function f_sizeof_r1

  pure function f_sizeof_r2(datatype) result(s)
    real(f_simple), dimension(:,:), intent(in) :: datatype
    integer(f_long) :: s; s=product(int(shape(datatype),f_long))*int(kind(datatype),f_long)
  end function f_sizeof_r2

  pure function f_sizeof_r3(datatype) result(s)
    real(f_simple), dimension(:,:,:), intent(in) :: datatype
    integer(f_long) :: s; s=product(int(shape(datatype),f_long))*int(kind(datatype),f_long)
  end function f_sizeof_r3

  pure function f_sizeof_r4(datatype) result(s)
    real(f_simple), dimension(:,:,:,:), intent(in) :: datatype
    integer(f_long) :: s; s=product(int(shape(datatype),f_long))*int(kind(datatype),f_long)
  end function f_sizeof_r4

  pure function f_sizeof_r5(datatype) result(s)
    real(f_simple), dimension(:,:,:,:,:), intent(in) :: datatype
    integer(f_long) :: s; s=product(int(shape(datatype),f_long))*int(kind(datatype),f_long)
  end function f_sizeof_r5

  pure function f_sizeof_c0(datatype) result(s)
    character(len=*), intent(in) :: datatype
    integer(f_long) :: s; s=int(len(datatype),f_long)
  end function f_sizeof_c0
  pure function f_sizeof_c1(ln,datatype) result(s)
    integer, intent(in) :: ln
    character(len=ln), dimension(:), intent(in) :: datatype
    integer(f_long) :: s; s=product(int(shape(datatype),f_long))*int(ln,f_long)
  end function f_sizeof_c1
  pure function f_sizeof_c2(ln,datatype) result(s)
    integer, intent(in) :: ln
    character(len=ln), dimension(:,:), intent(in) :: datatype
    integer(f_long) :: s; s=product(int(shape(datatype),f_long))*int(ln,f_long)
  end function f_sizeof_c2


  pure function f_sizeof_z1(datatype) result(s)
    complex(f_double), dimension(:), intent(in) :: datatype
    integer(f_long) :: s; s=product(int(shape(datatype),f_long))*int(2*kind(datatype),f_long)
  end function f_sizeof_z1
  pure function f_sizeof_z2(datatype) result(s)
    complex(f_double), dimension(:,:), intent(in) :: datatype
    integer(f_long) :: s; s=product(int(shape(datatype),f_long))*int(2*kind(datatype),f_long)
  end function f_sizeof_z2
  pure function f_sizeof_z3(datatype) result(s)
    complex(f_double), dimension(:,:,:), intent(in) :: datatype
    integer(f_long) :: s; s=product(int(shape(datatype),f_long))*int(2*kind(datatype),f_long)
  end function f_sizeof_z3
  pure function f_sizeof_z4(datatype) result(s)
    complex(f_double), dimension(:,:,:,:), intent(in) :: datatype
    integer(f_long) :: s; s=product(int(shape(datatype),f_long))*int(2*kind(datatype),f_long)
  end function f_sizeof_z4
  pure function f_sizeof_z5(datatype) result(s)
    complex(f_double), dimension(:,:,:,:,:), intent(in) :: datatype
    integer(f_long) :: s; s=product(int(shape(datatype),f_long))*int(2*kind(datatype),f_long)
  end function f_sizeof_z5

  pure function f_size_i0(datatype) result(s)
    integer(f_integer), intent(in) :: datatype
    integer(f_long) :: s
    !local variable
    integer :: mt; mt=kind(datatype)
    s=int(1,f_long)
  end function f_size_i0
  pure function f_size_i1(datatype) result(s)
    integer(f_integer), dimension(:), intent(in) :: datatype
    integer(f_long) :: s; s=product(int(shape(datatype),f_long))!*int(kind(datatype),f_long)
  end function f_size_i1
  pure function f_size_i2(datatype) result(s)
    integer(f_integer), dimension(:,:), intent(in) :: datatype
    integer(f_long) :: s; s=product(int(shape(datatype),f_long))!*int(kind(datatype),f_long)
  end function f_size_i2
  pure function f_size_i3(datatype) result(s)
    integer(f_integer), dimension(:,:,:), intent(in) :: datatype
    integer(f_long) :: s; s=product(int(shape(datatype),f_long))!*int(kind(datatype),f_long)
  end function f_size_i3
  pure function f_size_i4(datatype) result(s)
    integer(f_integer), dimension(:,:,:,:), intent(in) :: datatype
    integer(f_long) :: s; s=product(int(shape(datatype),f_long))!*int(kind(datatype),f_long)
  end function f_size_i4
  pure function f_size_i5(datatype) result(s)
    integer(f_integer), dimension(:,:,:,:,:), intent(in) :: datatype
    integer(f_long) :: s; s=product(int(shape(datatype),f_long))!*int(kind(datatype),f_long)
  end function f_size_i5


  pure function f_size_li1(datatype) result(s)
    integer(f_long), dimension(:), intent(in) :: datatype
    integer(f_long) :: s; s=product(int(shape(datatype),f_long))!*int(kind(datatype),f_long)
  end function f_size_li1
  pure function f_size_li2(datatype) result(s)
    integer(f_long), dimension(:,:), intent(in) :: datatype
    integer(f_long) :: s; s=product(int(shape(datatype),f_long))!*int(kind(datatype),f_long)
  end function f_size_li2
  pure function f_size_li3(datatype) result(s)
    integer(f_long), dimension(:,:,:), intent(in) :: datatype
    integer(f_long) :: s; s=product(int(shape(datatype),f_long))!*int(kind(datatype),f_long)
  end function f_size_li3
  pure function f_size_li4(datatype) result(s)
    integer(f_long), dimension(:,:,:,:), intent(in) :: datatype
    integer(f_long) :: s; s=product(int(shape(datatype),f_long))!*int(kind(datatype),f_long)
  end function f_size_li4
  pure function f_size_li5(datatype) result(s)
    integer(f_long), dimension(:,:,:,:,:), intent(in) :: datatype
    integer(f_long) :: s; s=product(int(shape(datatype),f_long))!*int(kind(datatype),f_long)
  end function f_size_li5


  pure function f_size_l1(datatype) result(s)
    logical, dimension(:), intent(in) :: datatype
    integer(f_long) :: s; s=product(int(shape(datatype),f_long))!*int(kind(datatype),f_long)
  end function f_size_l1
  pure function f_size_l2(datatype) result(s)
    logical, dimension(:,:), intent(in) :: datatype
    integer(f_long) :: s; s=product(int(shape(datatype),f_long))!*int(kind(datatype),f_long)
  end function f_size_l2
  pure function f_size_l3(datatype) result(s)
    logical, dimension(:,:,:), intent(in) :: datatype
    integer(f_long) :: s; s=product(int(shape(datatype),f_long))!*int(kind(datatype),f_long)
  end function f_size_l3
  pure function f_size_l4(datatype) result(s)
    logical, dimension(:,:,:,:), intent(in) :: datatype
    integer(f_long) :: s; s=product(int(shape(datatype),f_long))!*int(kind(datatype),f_long)
  end function f_size_l4
  pure function f_size_l5(datatype) result(s)
    logical, dimension(:,:,:,:,:), intent(in) :: datatype
    integer(f_long) :: s; s=product(int(shape(datatype),f_long))!*int(kind(datatype),f_long)
  end function f_size_l5

  pure function f_size_b1(datatype) result(s)
    logical(f_byte), dimension(:), intent(in) :: datatype
    integer(f_long) :: s; s=product(int(shape(datatype),f_long))!*int(kind(datatype),f_long)
  end function f_size_b1
  pure function f_size_b2(datatype) result(s)
    logical(f_byte), dimension(:,:), intent(in) :: datatype
    integer(f_long) :: s; s=product(int(shape(datatype),f_long))!*int(kind(datatype),f_long)
  end function f_size_b2
  pure function f_size_b3(datatype) result(s)
    logical(f_byte), dimension(:,:,:), intent(in) :: datatype
    integer(f_long) :: s; s=product(int(shape(datatype),f_long))!*int(kind(datatype),f_long)
  end function f_size_b3


  pure function f_size_d0(datatype) result(s)
    real(f_double), intent(in) :: datatype
    integer(f_long) :: s
    !local variable
    integer :: mt; mt=kind(datatype)
    s=int(1,f_long)
  end function f_size_d0
  pure function f_size_d1(datatype) result(s)
    real(f_double), dimension(:), intent(in) :: datatype
    integer(f_long) :: s; s=product(int(shape(datatype),f_long))!*int(kind(datatype),f_long)
  end function f_size_d1

  pure function f_size_d2(datatype) result(s)
    real(f_double), dimension(:,:), intent(in) :: datatype
    integer(f_long) :: s; s=product(int(shape(datatype),f_long))!*int(kind(datatype),f_long)
  end function f_size_d2

  pure function f_size_d3(datatype) result(s)
    real(f_double), dimension(:,:,:), intent(in) :: datatype
    integer(f_long) :: s; s=product(int(shape(datatype),f_long))!*int(kind(datatype),f_long)
  end function f_size_d3

  pure function f_size_d4(datatype) result(s)
    real(f_double), dimension(:,:,:,:), intent(in) :: datatype
    integer(f_long) :: s; s=product(int(shape(datatype),f_long))!*int(kind(datatype),f_long)
  end function f_size_d4

  pure function f_size_d5(datatype) result(s)
    real(f_double), dimension(:,:,:,:,:), intent(in) :: datatype
    integer(f_long) :: s; s=product(int(shape(datatype),f_long))!*int(kind(datatype),f_long)
  end function f_size_d5

  pure function f_size_d6(datatype) result(s)
    real(f_double), dimension(:,:,:,:,:,:), intent(in) :: datatype
    integer(f_long) :: s; s=product(int(shape(datatype),f_long))!*int(kind(datatype),f_long)
  end function f_size_d6

  pure function f_size_d7(datatype) result(s)
    real(f_double), dimension(:,:,:,:,:,:,:), intent(in) :: datatype
    integer(f_long) :: s; s=product(int(shape(datatype),f_long))!*int(kind(datatype),f_long)
  end function f_size_d7


  pure function f_size_r1(datatype) result(s)
    real(f_simple), dimension(:), intent(in) :: datatype
    integer(f_long) :: s; s=product(int(shape(datatype),f_long))!*int(kind(datatype),f_long)
  end function f_size_r1

  pure function f_size_r2(datatype) result(s)
    real(f_simple), dimension(:,:), intent(in) :: datatype
    integer(f_long) :: s; s=product(int(shape(datatype),f_long))!*int(kind(datatype),f_long)
  end function f_size_r2

  pure function f_size_r3(datatype) result(s)
    real(f_simple), dimension(:,:,:), intent(in) :: datatype
    integer(f_long) :: s; s=product(int(shape(datatype),f_long))!*int(kind(datatype),f_long)
  end function f_size_r3

  pure function f_size_r4(datatype) result(s)
    real(f_simple), dimension(:,:,:,:), intent(in) :: datatype
    integer(f_long) :: s; s=product(int(shape(datatype),f_long))!*int(kind(datatype),f_long)
  end function f_size_r4

  pure function f_size_r5(datatype) result(s)
    real(f_simple), dimension(:,:,:,:,:), intent(in) :: datatype
    integer(f_long) :: s; s=product(int(shape(datatype),f_long))!*int(kind(datatype),f_long)
  end function f_size_r5

  pure function f_size_c0(datatype) result(s)
    character(len=*), intent(in) :: datatype
    integer(f_long) :: s; s=1
  end function f_size_c0
  pure function f_size_c1(ln,datatype) result(s)
    integer, intent(in) :: ln
    character(len=ln), dimension(:), intent(in) :: datatype
    integer(f_long) :: s; s=product(int(shape(datatype),f_long))!*int(ln,f_long)
  end function f_size_c1

  pure function f_size_z1(datatype) result(s)
    complex(f_double), dimension(:), intent(in) :: datatype
    integer(f_long) :: s; s=product(int(shape(datatype),f_long))!*int(2*kind(datatype),f_long)
  end function f_size_z1
  pure function f_size_z2(datatype) result(s)
    complex(f_double), dimension(:,:), intent(in) :: datatype
    integer(f_long) :: s; s=product(int(shape(datatype),f_long))!*int(2*kind(datatype),f_long)
  end function f_size_z2
  pure function f_size_z3(datatype) result(s)
    complex(f_double), dimension(:,:,:), intent(in) :: datatype
    integer(f_long) :: s; s=product(int(shape(datatype),f_long))!*int(2*kind(datatype),f_long)
  end function f_size_z3
  pure function f_size_z4(datatype) result(s)
    complex(f_double), dimension(:,:,:,:), intent(in) :: datatype
    integer(f_long) :: s; s=product(int(shape(datatype),f_long))!*int(2*kind(datatype),f_long)
  end function f_size_z4
  pure function f_size_z5(datatype) result(s)
    complex(f_double), dimension(:,:,:,:,:), intent(in) :: datatype
    integer(f_long) :: s; s=product(int(shape(datatype),f_long))!*int(2*kind(datatype),f_long)
  end function f_size_z5


end module f_utils
