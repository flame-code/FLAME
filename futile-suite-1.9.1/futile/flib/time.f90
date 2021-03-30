!> @file
!!  Define routines for timing
!! @author
!!    Copyright (C) 2010-2013 BigDFT group (LG)
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Module containing variables used for the timing for BigDFT
module time_profiling
  use f_precisions
  use dictionaries
  use f_utils, only: f_time
  implicit none

  private 

  !>maximum number of allowed categories
  integer, parameter :: ncat_max=300
  !>maximum number of partial counters active
  integer, parameter :: nctr_max=20
  !>integer indicating the file unit
  integer, parameter :: timing_unit=60
  !>integer indicating the tabbing of the file unit
  integer, parameter :: tabfile=25
  !>maximum number of nested profiling levels
  integer, parameter :: max_ctrl = 5 
  !>id of active control structure (<=max_ctrl)
  integer :: ictrl=0                 
  !entries of the dictionaries 
  character(len=*), parameter :: catname='Category'
  character(len=*), parameter :: grpname='Group'
  character(len=*), parameter :: catinfo='Info'
  !uninitalized category, used for distinguish if a cat_id is active of not
  integer, parameter, public :: TIMING_UNINITIALIZED=-2

  !unspecified timings
  integer, save :: TCAT_UNSPECIFIED=TIMING_UNINITIALIZED !<after initialization this should be 0
  !Error codes
  integer, public, save :: TIMING_INVALID
 
  !> Contains all global variables associated to time profiling
  type :: time_ctrl
     logical :: master !<flag to store whether the instance can write on a file
     logical :: debugmode !<flag to store how to process the information
     integer :: cat_on !<id of the active category
     integer :: cat_paused !<id of paused category when interrupt action
     integer :: timing_ncat !<number of categories
     integer :: timing_nctr !<number of partial counters
     integer(f_long) :: epoch !<time of the creation of the routine
     double precision :: time0 !<reference time since last checkpoint
     double precision :: t0 !<reference time since last opening action
     double precision, dimension(ncat_max+1) :: clocks       !< timings of different categories
     double precision, dimension(nctr_max) :: counter_clocks !< times of the partial counters
     character(len=10), dimension(nctr_max) :: counter_names !< names of the partial counters, to be assigned
     character(len=128) :: report_file                    !< name of the file to write the report on
     type(dictionary), pointer :: dict_timing_categories  !< categories definitions
     type(dictionary), pointer :: dict_timing_groups      !< group definitions
  end type time_ctrl

  !>global variable controlling the different instances of the calls
  type(time_ctrl), dimension(max_ctrl) :: times

  public :: f_timing_reset,f_timing,f_timing_checkpoint,f_timing_stop,timing_errors
  public :: f_timing_category,f_timing_category_group,f_timing_finalize,f_timing_initialize
  public :: get_category_name,f_clock,f_profile,f_profile_end

  contains

    pure function time_ctrl_null() result(time)
      implicit none
      type(time_ctrl) :: time
      call nullify_time_ctrl(time)
    end function time_ctrl_null
    pure subroutine nullify_time_ctrl(time)
      implicit none
      type(time_ctrl), intent(out) :: time
      
      time%master=.false.
      time%debugmode=.false. 
      time%cat_on=0 
      time%cat_paused=0 
      time%timing_ncat=0 
      time%timing_nctr=0
      time%epoch=f_time() !take the initial time
      time%time0=0.d0
      time%t0=0.d0
      time%clocks=0.d0
      time%counter_clocks=0.d0
      time%counter_names=repeat(' ',len(time%counter_names))
      time%report_file=repeat(' ',len(time%report_file))
      nullify(time%dict_timing_categories)
      nullify(time%dict_timing_groups)
    end subroutine nullify_time_ctrl


    !> Check if the module has been initialized
    subroutine check_initialization()
      implicit none
      if (ictrl==0) then 
         call f_lib_err_severe_external(&
              'Timing library not initialized, f_lib_initialized should be called')
      end if
    end subroutine check_initialization


    !> For the moment the timing callback is a severe error.
    !! we should decide what to do to override this
    subroutine f_timing_callback()
      use yaml_output
      use dictionaries, only: f_err_severe
      implicit none
      !local variables
      integer, parameter :: iunit=96
      integer :: iunit_def,iunt,istat
      call yaml_warning('An error occured in timing module!')
      !retrieve current unit
      call yaml_get_default_stream(iunit_def)
      iunt=iunit_def
      !inquire for presence of unit iunit
      call yaml_set_stream(unit=iunit,filename='timing-categories.yaml',&
           position='rewind',setdefault=.false.,istat=istat)
      if (istat == 0) then
         call yaml_comment('Dumping active categories in file timing-categories.yaml')
         iunt=iunit
      end if
      call yaml_map('Dictionary of category groups',&
           times(ictrl)%dict_timing_groups,unit=iunt)
      call yaml_map('Dictionary of active categories',&
           times(ictrl)%dict_timing_categories,unit=iunt)

      !then close the file
      if (iunt /= iunit_def) then
         call yaml_close_stream(unit=iunt)
      end if

      call f_err_severe()
    end subroutine f_timing_callback

    !> initialize error codes of timing module
    subroutine timing_errors()
      use dictionaries, only: f_err_define
      implicit none
       call f_err_define(err_name='TIMING_INVALID',err_msg='Error in timing routines',&
            err_id=TIMING_INVALID,&
            err_action='Control the running conditions of f_timing routines called',&
            callback=f_timing_callback)
    end subroutine timing_errors

    !> define a class, which is a group of categories
    subroutine f_timing_category_group(grp_name,grp_info)
      implicit none
      character(len=*), intent(in) :: grp_name !< name of the class
      character(len=*), intent(in) :: grp_info !< description of it
      !local variables

      call check_initialization()

      if (grp_name .in. times(ictrl)%dict_timing_groups) then
         call f_err_throw('The timing category group '//grp_name//' has already been defined',&
              err_id=TIMING_INVALID)
      end if
      !in case of dry run override the commentary nonetheless  
      call set(times(ictrl)%dict_timing_groups//grp_name,grp_info)

    end subroutine f_timing_category_group

    !> define a new timing category with its description
    subroutine f_timing_category(cat_name,grp_name,cat_info,cat_id)
      use yaml_strings, only: yaml_toa
      implicit none
      character(len=*), intent(in) :: cat_name !< name of the category
      character(len=*), intent(in) :: grp_name !<class to which category belongs (see f_timing_class)
      character(len=*), intent(in) :: cat_info !< description of it
      integer, intent(out) :: cat_id !< id of the defined class, to be used for reference
      !local variables
      type(dictionary), pointer :: dict_cat
     
      call check_initialization()

      !check that the time has not started yet 
      if (times(ictrl)%time0/=0.d0) then
         call f_err_throw('Categories cannot be initialized when time counting started',&
              err_id=TIMING_INVALID)
         return
      end if

      if (.not. (grp_name .in. times(ictrl)%dict_timing_groups)) then
         call f_err_throw('The timing category group '//grp_name//' has not been defined',&
              err_id=TIMING_INVALID)
         return
      end if

      !then proceed to the definition of the category
      cat_id=dict_len(times(ictrl)%dict_timing_categories)

      call dict_init(dict_cat)
      call set(dict_cat//catname,cat_name)
      call set(dict_cat//grpname,grp_name)
      call set(dict_cat//catinfo,cat_info)

      call add(times(ictrl)%dict_timing_categories,dict_cat)
      times(ictrl)%timing_ncat=times(ictrl)%timing_ncat+1
      if (times(ictrl)%timing_ncat > ncat_max) then
         call f_err_throw('The number of initialized categories cannot exceed'//&
              trim(yaml_toa(ncat_max))//'. Change ncat_max in profile_time module',&
              err_id=TIMING_INVALID)
      end if
    end subroutine f_timing_category


    !> Initialize the timing by putting to zero all the chronometers
    subroutine f_timing_initialize()
      use yaml_strings, only: yaml_toa
      implicit none
      !create the general category for unspecified timings
      ictrl=ictrl+1
      if (f_err_raise(ictrl > max_ctrl,&
           'Timing: the number of active instances cannot exceed'//&
           trim(yaml_toa(max_ctrl)),TIMING_INVALID)) return
      call nullify_time_ctrl(times(ictrl))

      if (ictrl==1) then
         call dict_init(times(ictrl)%dict_timing_groups)
         call dict_init(times(ictrl)%dict_timing_categories)
      else
         !pre-existing categories can be modified
         times(ictrl)%dict_timing_groups=>times(ictrl-1)%dict_timing_groups
         
      end if

      !define the main groups and categories
      call f_timing_category_group('NULL','Nullified group to contain unspecifed category')
      call f_timing_category('UNSPECIFIED','NULL',&
           'Unspecified category collecting garbage timings',TCAT_UNSPECIFIED)
      times(ictrl)%timing_ncat=0
    end subroutine f_timing_initialize


    !> Get the walltime since most recent call of the f_timing initialize
    function f_clock()
      implicit none
      integer(f_long) :: f_clock !< elapsed walltime since last call of the initialize
      f_clock=f_time()-times(ictrl)%epoch
    end function f_clock


    !> Finalize the timing by putting to zero all the chronometers
    subroutine f_timing_finalize(walltime)
      use yaml_output
      implicit none
      integer(f_long), intent(out), optional :: walltime !< elapsed walltime since last call of the initialize
      !create the general category for unspecified timings
      call dict_free(times(ictrl)%dict_timing_categories)
      call dict_free(times(ictrl)%dict_timing_groups)
      !put to zero the number of categories
      if (present(walltime)) walltime=f_clock()
      times(ictrl)=time_ctrl_null()
      ictrl=ictrl-1
    end subroutine f_timing_finalize


    !> Re-initialize the timing by putting to zero all the chronometers (old action IN)
    subroutine f_timing_reset(filename,master,verbose_mode)
      use yaml_output, only: yaml_new_document
      implicit none
      logical, intent(in) :: master !<true if the task is the one responsible for file writing
      character(len=*), intent(in) :: filename !<name of the file where summary have to be written
      !>Toggle verbose mode in the module. In case of parallel execution all the processors will
      !! print out their information on the counters.
      logical, intent(in), optional :: verbose_mode 
      !local variables
      integer :: ictr,i,iunit_def
      integer(kind=8) :: itns

      !global timer
      itns=f_time()

      call check_initialization()

      !check if some categories have been initialized
      if (f_err_raise(times(ictrl)%timing_ncat==0,'No timing categories have been initialized, no counters to reset .'//&
           'Use f_timing_category(_group) routine(s)',err_id=TIMING_INVALID)) return

      times(ictrl)%time0=real(itns,kind=8)*1.d-9
      !reset partial counters and categories
      do i=1,times(ictrl)%timing_ncat
         times(ictrl)%clocks(i)=0.d0
      end do
      do ictr=1,times(ictrl)%timing_nctr
         times(ictrl)%counter_clocks(ictr)=0.d0
      enddo
      !store filename where report have to be written
      !default stream can be used when the filename is empty
      times(ictrl)%report_file(1:len(times(ictrl)%report_file))=filename
      !store debug mode
      if (present(verbose_mode)) then
         times(ictrl)%debugmode=verbose_mode
      else
         times(ictrl)%debugmode=.false.
      end if

      !no category has been used so far
      !init=.false.
      times(ictrl)%master=master
      times(ictrl)%cat_on=0
      times(ictrl)%cat_paused=0 !no stopped category
      times(ictrl)%timing_nctr=0 !no partial counters activated
      !initialize the document
      if (times(ictrl)%master) then
         call timing_open_stream(iunit_def)
         call yaml_new_document() !in principle is active only when the document is released
         call timing_close_stream(iunit_def)
      end if
    end subroutine f_timing_reset


    !> Perform a checkpoint of the chronometer with a partial counter
    !! the last active category is halted and a summary of the timing 
    !! is printed out
    subroutine f_timing_checkpoint(ctr_name,mpi_comm,nproc,gather_routine)
      use yaml_output, only: yaml_map
      use yaml_strings, only: yaml_toa
      use dynamic_memory
      implicit none
      !> name of the partial counter for checkpoint identification
      character(len=*), intent(in) :: ctr_name 
      !> handle of the mpi_communicator associated with the checkpoint 
      integer, intent(in), optional :: mpi_comm,nproc 
!!$      interface
!!$         subroutine gather_routine(ndata,nproc,mpi_comm,src,dest)
!!$           implicit none
!!$           integer, intent(in) :: ndata !< number of categories of the array
!!$           integer, intent(in) :: nproc !< number of MPI tasks
!!$           real(kind=8), dimension(ndata), intent(in) :: src !< total timings of the instance
!!$           real(kind=8), dimension(ndata,nproc), intent(inout) :: dest 
!!$        end interface
      external :: gather_routine
      optional :: gather_routine !< routine to perform the gathering of counters
      !local variables
      integer :: i,nnodes
      integer(kind=8) :: itns
      logical, dimension(3) :: test

      !global timer
      itns=f_time()

      call check_initialization()

      !stop partial counters and restart from the beginning
      if (times(ictrl)%cat_on/=0) then
         call f_err_throw('TIMING IS INITIALIZED BEFORE PARTIAL RESULTS'//&
              trim(yaml_toa(times(ictrl)%cat_on)),&
              err_id=TIMING_INVALID)
      end if
      test=(/present(mpi_comm),present(nproc),present(gather_routine)/)
      if (any(test) .and. .not. all(test)) then
         call f_err_throw('mpi_comm and nproc should be present together '//&
             'and consistent with each other. Also gather_routine should appear',&
              err_id=TIMING_INVALID)
      end if
      nnodes=1
      if (present(nproc)) nnodes=nproc

      times(ictrl)%timing_nctr=times(ictrl)%timing_nctr+1
      if (f_err_raise(times(ictrl)%timing_nctr > nctr_max,&
           'Max No. of partial counters reached',err_id=TIMING_INVALID)) return
      !name of the partial counter
      times(ictrl)%counter_names(times(ictrl)%timing_nctr)(1:len(times(ictrl)%counter_names))=&
           trim(ctr_name)
      !total time elapsed in it
      times(ictrl)%clocks(times(ictrl)%timing_ncat+1)&
           =real(itns,kind=8)*1.d-9-times(ictrl)%time0
      times(ictrl)%counter_clocks(times(ictrl)%timing_nctr)=&
           times(ictrl)%clocks(times(ictrl)%timing_ncat+1)

      if (present(mpi_comm)) then
         call gather_and_dump_results(times(ictrl)%master,&
              times(ictrl)%timing_ncat,nnodes,&
              times(ictrl)%counter_names(times(ictrl)%timing_nctr),&
              times(ictrl)%clocks,mpi_comm,gather_routine)
      else
         call gather_and_dump_results(times(ictrl)%master,&
              times(ictrl)%timing_ncat,nnodes,&
              times(ictrl)%counter_names(times(ictrl)%timing_nctr),&
              times(ictrl)%clocks)
      end if

!!$      call sum_results(times(ictrl)%timing_ncat,mpi_comm,&
!!$           times(ictrl)%counter_names(times(ictrl)%timing_nctr),&
!!$           times(ictrl)%clocks)

      !reset all timings
      times(ictrl)%time0=real(itns,kind=8)*1.d-9
      do i=1,times(ictrl)%timing_ncat
         times(ictrl)%clocks(i)=0.d0
      enddo

      !temporary: just dump the present categories
      !call f_timing_callback()

    end subroutine f_timing_checkpoint

    subroutine gather_and_dump_results(master,ncat,nnodes,message,clocks,mpi_comm,gather_routine,dict_info)
      use dynamic_memory
      implicit none
      logical, intent(in) :: master
      integer, intent(in) :: ncat,nnodes
      character(len=*), intent(in) :: message
      double precision, dimension(ncat+1), intent(in) :: clocks
      integer, intent(in), optional :: mpi_comm
      type(dictionary), pointer, optional :: dict_info
      external :: gather_routine
      optional :: gather_routine
      real(kind=8), dimension(:,:), allocatable :: timeall

      !allocate total timings
      timeall=f_malloc((/1.to.ncat+1,0.to.nnodes/),&
           id='timeall')
      !gather the results
      if (nnodes > 1) then
         if (present(mpi_comm)) call gather_routine(ncat+1,nnodes,mpi_comm,&
              clocks,timeall)
      else
         call f_memcpy(src=clocks,dest=timeall)
      end if
      if (master) then
         call timing_dump_results(ncat,nnodes,trim(message),timeall,dict_info)
      endif
      call f_free(timeall)
    end subroutine gather_and_dump_results

    subroutine gather_and_dump_counters(master,ncnt,nnodes,pcnames,pctimes,dict_info,&
         mpi_comm,gather_routine)
      use dynamic_memory
      implicit none
      logical, intent(in) :: master
      integer, intent(in) :: ncnt,nnodes
      double precision, dimension(ncnt), intent(in) :: pctimes
      character(len=10), dimension(ncnt), intent(in) :: pcnames
      type(dictionary), pointer :: dict_info
      integer, intent(in), optional :: mpi_comm
      external :: gather_routine
      optional :: gather_routine
      !local variables
      double precision, dimension(:,:), allocatable :: timecnt 

      !allocate total timings
      timecnt=f_malloc((/1.to.ncnt,0.to.nnodes/),&
           id='timecnt')
      !gather the results
      if (nnodes > 1) then
         if (present(mpi_comm)) call gather_routine(ncnt,nnodes,mpi_comm,&
              pctimes,timecnt)
      else
         call f_memcpy(src=pctimes,dest=timecnt)
      end if
      if (master) then
         call timing_dump_counters(ncnt,nnodes,pcnames,timecnt,dict_info)
      endif
      call f_free(timecnt)
    end subroutine gather_and_dump_counters

    !>stop the timing and dump information of the partial counters
    subroutine f_timing_stop(mpi_comm,nproc,gather_routine,dict_info)
      use dynamic_memory
      use dictionaries
      implicit none
      !> handle of the mpi_communicator associated with the results
      integer, intent(in), optional :: mpi_comm,nproc 
      type(dictionary), pointer, optional :: dict_info
!!$      interface
!!$         subroutine gather_routine(ndata,nproc,mpi_comm,src,dest)
!!$           implicit none
!!$           integer, intent(in) :: ndata !< number of categories of the array
!!$           integer, intent(in) :: nproc !< number of MPI tasks
!!$           real(kind=8), dimension(ndata), intent(in) :: src !< total timings of the instance
!!$           real(kind=8), dimension(ndata,nproc), intent(inout) :: dest 
!!$        end interface
      external :: gather_routine
      optional :: gather_routine !< routine to perform the gathering of counters

      !local variables
      logical, dimension(3) :: test
      integer :: nnodes
      integer(kind=8) :: itns
      type(dictionary), pointer :: dict_tmp
!!$      !$ integer :: omp_get_max_threads

      !global timer
      itns=f_time()

      !stop partial counters and restart from the beginning
      if (times(ictrl)%cat_on/=0) then
         call f_err_throw('TIMING IS INITIALIZED BEFORE RESULTS',&
              err_id=TIMING_INVALID)
      end if

      test=(/present(mpi_comm),present(nproc),present(gather_routine)/)
      if (any(test) .and. .not. all(test)) then
         call f_err_throw('Timing Stop: mpi_comm and nproc should be present together '//&
              'and consistent with each other. Also gather_routine should appear',&
              err_id=TIMING_INVALID)
      end if
      nnodes=1
      if (present(nproc)) nnodes=nproc

      if (times(ictrl)%timing_nctr == 0) then !no partial counters selected
         times(ictrl)%clocks(times(ictrl)%timing_ncat+1)&
              =real(itns,kind=8)*1.d-9-times(ictrl)%time0

         if (present(mpi_comm)) then
            call gather_and_dump_results(times(ictrl)%master,&
                 times(ictrl)%timing_ncat,nnodes,'ALL',times(ictrl)%clocks,mpi_comm,gather_routine,dict_info=dict_info)
         else
            call gather_and_dump_results(times(ictrl)%master,&
                 times(ictrl)%timing_ncat,nnodes,'ALL',times(ictrl)%clocks,dict_info=dict_info)
         end if

!!$         call sum_results(times(ictrl)%timing_ncat,mpi_comm,'ALL',&
!!$              times(ictrl)%clocks)
      else !consider only the results of the partial counters
         !creation of the dict_info        
         if (present(dict_info)) then
            dict_tmp=>dict_info
         else
            nullify(dict_tmp)
         end if
         if (present(mpi_comm)) then
            call gather_and_dump_counters(times(ictrl)%master,&
                 times(ictrl)%timing_nctr,nnodes,times(ictrl)%counter_names,&
                 times(ictrl)%counter_clocks,dict_tmp,&
                 mpi_comm,gather_routine)
         else
            call gather_and_dump_counters(times(ictrl)%master,&
                 times(ictrl)%timing_nctr,nnodes,times(ictrl)%counter_names,&
                 times(ictrl)%counter_clocks,dict_tmp)
         end if
!!$         call sum_counters(times(ictrl)%counter_clocks,&
!!$              times(ictrl)%counter_names,times(ictrl)%timing_nctr,mpi_comm,&
!!$              times(ictrl)%debugmode)
      end if

      !restore timing, categories can be manipulated now
      times(ictrl)%time0=0.d0
      times(ictrl)%timing_nctr=0 !no partial counters activated anymore
    end subroutine f_timing_stop


    !> The same timing routine but with system_clock (in case of a supported specs)
    subroutine f_timing(cat_id,action)
      use dictionaries, only: f_err_raise,f_err_throw
      use f_utils, only: f_time
      use yaml_strings, only: yaml_toa
      use nvtx !for nvidia profiler
      implicit none
      !Variables
      integer, intent(in) :: cat_id
      character(len=2), intent(in) :: action      ! possibilities: INitialize, ON, OFf, REsults
      !Local variables
      integer(kind=8) :: itns
      real(kind=8) :: t1
      character(len=max_field_length):: catname
      !$ include 'remove_omp-inc.f90'
      
      !first of all, read the time
      itns=f_time()

      !routine has no effect if the module has not been initialized
      !this is important to guarantee the possibility of
      !putting timing calls to lower-level routines which might have been called
      !before initialization
      if (ictrl==0 .or. cat_id==TIMING_UNINITIALIZED) return
      if (cat_id <= 0) then
         call f_err_throw('Timing category id must be a valid positive number',&
              err_id=TIMING_INVALID)
      end if
      select case(action)
      case('ON')
         if (times(ictrl)%cat_paused /=0) then
            !no action except for misuse of interrupts
            call f_err_throw('INTERRUPTS SHOULD NOT BE altered by ON',&
                 err_id=TIMING_INVALID)
            return
         end if
         !if some other category was initalized before, return (no action)
         if (times(ictrl)%cat_on /= 0) return
         times(ictrl)%t0=real(itns,kind=8)*1.d-9
         times(ictrl)%cat_on=cat_id !category which has been activated
        !catname=f_malloc((/1.to.128/), id="catname")
        call get_category_name(cat_id, catname)
         call nvtxrangepusha(catname//CHAR(0));
        !call f_free(catname)

!call Extrae_event(6000019,cat_id)
      case('OF')
         if (times(ictrl)%cat_paused /=0) then
            !no action except for misuse of interrupts
            call f_err_throw('INTERRUPTS SHOULD NOT BE ALTERED by OF',&
                 err_id=TIMING_INVALID)
            return
         else if (cat_id==times(ictrl)%cat_on) then 
            !switching off the good category
            t1=real(itns,kind=8)*1.d-9
            times(ictrl)%clocks(cat_id)=times(ictrl)%clocks(cat_id)+&
                 t1-times(ictrl)%t0
            times(ictrl)%cat_on=0
         call nvtxrangepop();
            !call Extrae_event(6000019,0)
         end if
         !otherwise no action as the off mismatches
      case('IR') !interrupt category
         if (times(ictrl)%cat_paused /=0) then
            call f_err_throw('Category No. '//&
                 trim(yaml_toa(times(ictrl)%cat_paused))//&
                 ' already interrupted ('//trim(yaml_toa(times(ictrl)%cat_on))//&
                 ' is active),  cannot interrupt again with'//&
                 trim(yaml_toa(cat_id)),&
                 err_id=TIMING_INVALID)
            return
         end if
         !time
         t1=real(itns,kind=8)*1.d-9
         if (times(ictrl)%cat_on /=0) then !there is already something active
            !stop the active counter
            times(ictrl)%clocks(times(ictrl)%cat_on)=&
                 times(ictrl)%clocks(times(ictrl)%cat_on)+t1-times(ictrl)%t0
            times(ictrl)%cat_paused=times(ictrl)%cat_on
         else
            times(ictrl)%cat_paused=-1 !start by pausing
         end if
         times(ictrl)%cat_on=cat_id
         times(ictrl)%t0=t1
         call get_category_name(cat_id, catname)
         call nvtxrangepusha(catname//CHAR(0));

      case('RS') !resume the category by supposing it has been activated by IR
         if (f_err_raise(times(ictrl)%cat_paused==0,&
              'It appears no category has to be resumed',&
              err_id=TIMING_INVALID)) return
         if (cat_id /= times(ictrl)%cat_on) then
            call f_err_throw('The category id '//trim(yaml_toa(cat_id))//&
                 ' is not active',&
                 err_id=TIMING_INVALID)
            return
         end if
         !time
         t1=real(itns,kind=8)*1.d-9
         times(ictrl)%clocks(cat_id)=times(ictrl)%clocks(cat_id)+&
              t1-times(ictrl)%t0
         !restore normal counter
         if (times(ictrl)%cat_paused/=-1) then
            times(ictrl)%cat_on=times(ictrl)%cat_paused
            times(ictrl)%t0=t1
         else
            times(ictrl)%cat_on=0
         end if
         times(ictrl)%cat_paused=0
         call nvtxrangepop();

      case('RX') !resume the interrupted category, cat_id is ignored here
         !dry run if expert mode active and not initialized categories
         if (times(ictrl)%cat_paused==0 .or. times(ictrl)%cat_on==0) return
!!$         if (f_err_raise(times(ictrl)%cat_paused==0 .or. &
!!$              times(ictrl)%cat_on==0 ,&
!!$              'It appears no category has to be resumed (RX case),'//&
!!$              ' control whether the interrupted category has been initialized',&
!!$              err_id=TIMING_INVALID)) return
         !time
         t1=real(itns,kind=8)*1.d-9
         times(ictrl)%clocks(times(ictrl)%cat_on)=&
              times(ictrl)%clocks(times(ictrl)%cat_on)+&
              t1-times(ictrl)%t0
         !restore normal counter
         if (times(ictrl)%cat_paused/=-1) then
            times(ictrl)%cat_on=times(ictrl)%cat_paused
            times(ictrl)%t0=t1
         else
            times(ictrl)%cat_on=0
         end if
         times(ictrl)%cat_paused=0
         call nvtxrangepop();

      case default
         call f_err_throw('TIMING ACTION UNDEFINED',err_id=TIMING_INVALID)
      end select

    END SUBROUTINE f_timing


    !> Opens the file of the timing unit
    subroutine timing_open_stream(iunit_def)
      use yaml_output, only: yaml_get_default_stream,yaml_set_stream
      implicit none
      integer, intent(out) :: iunit_def !< previous default unit
      !first get the default stream
      call yaml_get_default_stream(iunit_def)
      if (iunit_def /= timing_unit) then
         call yaml_set_stream(unit=timing_unit,&
              filename=trim(times(ictrl)%report_file),&
              record_length=120,tabbing=tabfile)
      end if

    end subroutine timing_open_stream


    !> Close the stream and restore old default unit
    subroutine timing_close_stream(iunit_def)
      use yaml_output, only: yaml_set_default_stream,yaml_close_stream
      implicit none
      integer, intent(in) :: iunit_def !< previous default unit
      !local variables
      integer :: ierr

      call yaml_set_default_stream(iunit_def,ierr)
      !close the previous one
      if (iunit_def /= timing_unit) call yaml_close_stream(unit=timing_unit)
    end subroutine timing_close_stream


    !> Dump the line of the timings for time.yaml form
    subroutine timing_dump_line(name,tabbing,pc,secs,unit,loads)
      use yaml_output
      use yaml_strings
      implicit none
      integer, intent(in) :: tabbing !<vlue of the tabbing for pretty printing
      double precision, intent(in) :: pc !< percent of the time for line id
      double precision, intent(in) :: secs !< seconds spent for line id
      character(len=*), intent(in) :: name !< id of the line printed
      integer, intent(in), optional :: unit !< @copydoc yaml_output::doc::unit
      !> extra info containing the load for each task for the line id, 
      !! calculated with respect to the average value given by secs
      double precision, dimension(0:), intent(in), optional :: loads 
      !local variables
      character(len=*), parameter :: fmt_pc='(f5.1)'
      character(len=*), parameter :: fmt_secs='(1pg9.2)'
      character(len=*), parameter :: fmt_extra='(f5.2)'
      integer :: iextra,nextra,unt

      unt=0
      if (present(unit)) unt=unit
      !determine the presence of extra information
      nextra=0
      if (present(loads)) nextra=size(loads)

      call yaml_sequence_open(name,flow=.true.,tabbing=tabbing,unit=unt)
      call yaml_sequence(yaml_toa(pc,fmt=fmt_pc),unit=unt)
      call yaml_sequence(yaml_toa(secs,fmt=fmt_secs),unit=unt)
      do iextra=0,nextra-1
         call yaml_sequence(yaml_toa(loads(iextra),fmt=fmt_extra),unit=unt)
      end do
      call yaml_sequence_close(unit=unt)


    end subroutine timing_dump_line


    !> Put the average value of timeall in the timesum array
    !! then rewrite each element with the deviation from it (in debug mode)
    !! in normal mode write only the max and min deviations (only in parallel)
    subroutine timing_data_synthesis(nproc,ncats,timeall,timesum_tot)
      implicit none
      integer, intent(in) :: nproc,ncats
      real(kind=8), dimension(ncats,0:nproc-1), intent(inout) :: timeall
      real(kind=8), dimension(ncats), intent(out) :: timesum_tot
      !local variables
      integer :: icat,jproc
      real(kind=8) :: tmin,tmax

      do icat=1,ncats
         timesum_tot(icat)=0.d0
         do jproc=0,nproc-1
            timesum_tot(icat)=timesum_tot(icat)+timeall(icat,jproc)
         end do
         timesum_tot(icat)=timesum_tot(icat)/real(nproc,kind=8)
         if (timesum_tot(icat)>0.d0) then
            if (times(ictrl)%debugmode) then
               do jproc=0,nproc-1
                  timeall(icat,jproc)=timeall(icat,jproc)/timesum_tot(icat)
               end do
            else if (nproc >1) then
               tmax=0.0d0
               tmin=1.0d300
               do jproc=0,nproc-1
                  tmax=max(timeall(icat,jproc),tmax)
                  tmin=min(timeall(icat,jproc),tmin)
               end do
               timeall(icat,0)=tmax/timesum_tot(icat)
               timeall(icat,1)=tmin/timesum_tot(icat)
            end if
         end if
      end do
    end subroutine timing_data_synthesis


    !> Dump the final information of the partial counters
    subroutine timing_dump_counters(ncounters,nproc,pcnames,timecnt,dict_info)
      use yaml_output
      use yaml_strings, only: yaml_date_and_time_toa
      implicit none
      integer, intent(in) :: ncounters,nproc
      character(len=10), dimension(ncounters), intent(in) :: pcnames
      double precision, dimension(ncounters,0:nproc), intent(inout) :: timecnt
      type(dictionary), pointer :: dict_info
      !local variables
      integer :: iunit_def,i
      double precision :: pc
      !synthesis of the counters
      call timing_data_synthesis(nproc,ncounters,timecnt,timecnt(1,nproc))

      call timing_open_stream(iunit_def)
      call yaml_mapping_open('SUMMARY',advance='no')
      call yaml_comment('     % ,  Time (s)',tabbing=tabfile)

      !sum all the information by counters
      do i=1,ncounters
         pc=100.d0*timecnt(i,nproc)/sum(timecnt(1:ncounters,nproc))
         call timing_dump_line(trim(pcnames(i)),tabfile,pc,timecnt(i,nproc))
      end do
      call timing_dump_line('Total',tabfile,100.d0,sum(timecnt(1:ncounters,nproc)))
      call yaml_mapping_close() !summary

      !dump extra info dictionary
      call dump_extra_info_dict(dict_info)
      !restore the default stream
      call timing_close_stream(iunit_def)

    end subroutine timing_dump_counters


    subroutine dump_extra_info_dict(dict_info)
      use yaml_strings, only: yaml_date_and_time_toa
      use yaml_output
      implicit none
      type(dictionary), pointer, optional :: dict_info
      
      if (.not. present(dict_info)) return
      if (associated(dict_info)) call yaml_dict_dump(dict_info)
      call yaml_map('Report timestamp',trim(yaml_date_and_time_toa()))
      
    end subroutine dump_extra_info_dict

    
    !> Dump the results of the nonzero timings of the categories in the file indicated by filename_time
    !! the array timesum should contain the timings for each processor (from 0 to nproc-1)
    !! and will also contain the average value (in position nproc)
    subroutine timing_dump_results(ncat,nproc,message,timeall,dict_info)
      use dynamic_memory
      use yaml_output
      use yaml_strings
      implicit none
      integer, intent(in) :: ncat,nproc
      character(len=*), intent(in) :: message
      real(kind=8), dimension(ncat+1,0:nproc), intent(inout) :: timeall
      type(dictionary), pointer, optional :: dict_info
      !local variables
      integer :: ncls,i,j,icls,icat,jproc,iunit_def,nextra
      real(kind=8) :: total_pc,pc
      type(dictionary), pointer :: dict_cat
      character(len=max_field_length) :: name
      integer, dimension(ncat) :: isort !< automatic array should be enough
      real(kind=8), dimension(:,:), allocatable :: timecls
      character(len=max_field_length), dimension(:), allocatable :: group_names

      call f_routine(id='timing_dump_results')

      ncls=dict_size(times(ictrl)%dict_timing_groups)-1 !the first is always null group
      !regroup the data for each category in any processor
      timecls=f_malloc0((/1.to.ncls,0.to.nproc/),id='timecls')

      !this has to be done via the dictionary of the categories
      !store the keys of the valid groups (eliminate the first)
      group_names=f_malloc_str(max_field_length,ncls+1,id='group_names')
      group_names=dict_keys(times(ictrl)%dict_timing_groups)

      dict_cat=>dict_iter(times(ictrl)%dict_timing_categories)
      !neglect the first one by calling dict_next immediately.
      do icat=1,ncat
         !categories are always in order
         dict_cat=>dict_next(dict_cat)
         
         if (.not. associated(dict_cat)) then
            call f_err_throw('Dictionary of categories not compatible with total number, icat='//icat,&
                 err_id=TIMING_INVALID)
            exit
         end if
         name=dict_cat//grpname

         !then for each processor adds the timing category to the corresponding group
         find_group: do icls=1,ncls
            if (trim(name)==trim(group_names(icls+1))) then
               do jproc=0,nproc-1
                  timecls(icls,jproc)=timecls(icls,jproc)+timeall(icat,jproc)
               end do
               exit find_group
            end if
         end do find_group
      end do

      !synthesis of the categories
      call timing_data_synthesis(nproc,ncat+1,timeall,timeall(1,nproc))
      !synthesis of the classes
      call timing_data_synthesis(nproc,ncls,timecls,timecls(1,nproc))
      
      if (nproc >1) then
         nextra=2
         if (times(ictrl)%debugmode) nextra=nproc
      else
         nextra=0
      end if

      !calculate the summary of the category
      call sort_positions(ncat,timeall(1,nproc),isort)

      !use yaml to write time.yaml
      call timing_open_stream(iunit_def)

      call yaml_mapping_open(trim(message),advance='no')
      if (nproc==1) then
         call yaml_comment('     % ,  Time (s)',tabbing=tabfile)
      else if (times(ictrl)%debugmode) then
         call yaml_comment('     % ,  Time (s), Load per MPI proc (relative) ',tabbing=tabfile)
      else
         call yaml_comment('     % ,  Time (s), Max, Min Load (relative) ',tabbing=tabfile)
      end if
      call yaml_mapping_open('Classes')
      total_pc=0.d0
      do icls=1,ncls
         pc=0.0d0
         if (timeall(ncat+1,nproc)/=0.d0) &
              pc=100.d0*timecls(icls,nproc)/timeall(ncat+1,nproc)!times(ictrl)%clocks(ncat+1)
         total_pc=total_pc+pc
         !only nonzero classes are printed out
         if (timecls(icls,nproc) /= 0.d0) then
            call timing_dump_line(trim(group_names(icls+1)),tabfile,pc,timecls(icls,nproc),&
                 loads=timecls(icls,0:nextra-1))
         end if
      end do
      call timing_dump_line('Total',tabfile,total_pc,times(ictrl)%clocks(ncat+1),&
           loads=timeall(ncat+1,0:nextra-1))
      call yaml_mapping_close() !Classes
      call yaml_mapping_open('Categories',advance='no')
      call yaml_comment('Ordered by time consumption')
      do j=1,ncat
         i=isort(j)
         pc=0.d0
         !only nonzero categories are printed out
         if (timeall(i,nproc) /= 0.d0) then
            dict_cat=>times(ictrl)%dict_timing_categories//i
            if (timeall(ncat+1,nproc)/=0.d0)&
                 pc=100.d0*timeall(i,nproc)/timeall(ncat+1,nproc)
            name=dict_cat//catname
            call yaml_mapping_open(trim(name))
            call timing_dump_line('Data',tabfile,pc,timeall(i,nproc),&
                 loads=timeall(i,0:nextra-1))
            name=dict_cat//grpname
            call yaml_map('Class',trim(name))
            name=dict_cat//catinfo
            call yaml_map('Info',trim(name))
            call yaml_mapping_close()
         end if
      enddo

      call yaml_mapping_close() !categories
      call yaml_mapping_close() !counter

      !dump extra info dictionary
      call dump_extra_info_dict(dict_info)
      !restore the default stream
      call timing_close_stream(iunit_def)

      call f_free_str(max_field_length,group_names)
      call f_free(timecls)
      call f_release_routine()
    end subroutine timing_dump_results

    !>extract the category name 
    subroutine get_category_name(cat_id,getname)
      use yaml_output, only: yaml_map
      implicit none
      integer, intent(in) :: cat_id
      character(len=*), intent(inout) :: getname
      !local variables
      character(len=max_field_length) :: name

      getname(1:len(getname))=' '
      if (cat_id > 0 .and. cat_id < dict_len(times(ictrl)%dict_timing_categories)) then
         name=times(ictrl)%dict_timing_categories//cat_id//catname
         getname(1:len(getname))=name
      end if
    end subroutine get_category_name

    subroutine f_profile(entry_point,id,repeat,jmpbuf,dump_results,unit)
      use f_jmp
      use yaml_output
      use f_utils, only: f_humantime
      implicit none
      logical, intent(in) :: dump_results
      character(len=*), intent(in) :: id
      integer, intent(in) :: repeat
      type(f_jmpbuf), intent(inout) :: jmpbuf
      integer, intent(in), optional :: unit !< stream id
      external :: entry_point
      !local variables
      integer(f_long) :: t1

      !first of all determine if setjmp has already been called
      jmpbuf%id=id
      jmpbuf%destroy_signal=repeat
      jmpbuf%t0 = f_time()
      jmpbuf%callback=f_loc(entry_point)
      !print *,'start'
      call f_jmpbuf_set(jmpbuf)
      !otherwise we have already set the jump point, nothing to do 

      !here we also end the calculation
      if (dump_results) then
         t1=f_time()
         call yaml_mapping_open('Ending profiling section',unit=unit,&
              flow=.true.)
         call yaml_map('Id',jmpbuf%id,unit=unit)
         call yaml_map('Times repeated',jmpbuf%destroy_signal,unit=unit)
         call yaml_map('Elapsed time','"'//trim(f_humantime(t1-jmpbuf%t0))//'"',unit=unit)
         call yaml_mapping_close(unit=unit)
      end if
      call f_jmpbuf_free(jmpbuf)

    end subroutine f_profile

    subroutine f_profile_end(entry_point,jmpbuf)
      use f_jmp
      implicit none
      type(f_jmpbuf), intent(inout) :: jmpbuf

      external :: entry_point

      !print *,'there',jmpbuf%signal,jmpbuf%destroy_signal,&
      !     associated(jmpbuf%jmp_buf)
      if (.not. associated(jmpbuf%jmp_buf)) return
      !print *,'there2',jmpbuf%signal,jmpbuf%destroy_signal
      if (f_loc(entry_point) /= jmpbuf%callback) &
           call f_err_throw('Inconsistent entry points, exiting',&
              err_id=TIMING_INVALID)

      !here we should jmp back if the jmpbuffer need it
      if (jmpbuf%signal /= jmpbuf%destroy_signal) then
         call f_longjmp(jmpbuf) !increase signal
      end if

    end subroutine f_profile_end

  end module time_profiling

!!!  subroutine sum_counters(pctimes,pcnames,ncounters,mpi_comm,debugmode)
!!!    use yaml_output
!!!    use dynamic_memory
!!!    use time_profiling, only: timing_unit,timing_dump_line,timing_data_synthesis,&
!!!         timing_open_stream,timing_close_stream,tabfile
!!!    use dictionaries
!!!    
!!!  implicit none
!!!  include 'mpif.h'
!!!  logical, intent(in) :: debugmode
!!!  integer, intent(in) :: mpi_comm,ncounters
!!!  real(kind=8), dimension(ncounters), intent(in) :: pctimes
!!!  character(len=10), dimension(ncounters), intent(in) :: pcnames
!!!  !local variables
!!!  logical :: parallel
!!!  integer :: i,ierr,iproc,jproc,icat,nthreads,namelen,iunit_def,nproc
!!!  real(kind=8) :: pc
!!!  
!!!  character(len=MPI_MAX_PROCESSOR_NAME) :: nodename_local
!!!  double precision, dimension(:,:), allocatable :: timecnt 
!!!  character(len=MPI_MAX_PROCESSOR_NAME), dimension(:), allocatable :: nodename
!!!  type(dictionary), pointer :: dict_info
!!!  !$ integer :: omp_cannot intget_max_threads
!!!
!!!  ! Not initialised case.
!!!  if (mpi_comm==MPI_COMM_NULL) return
!!!  call f_routine(id='sum_counters')
!!!  
!!!  call MPI_COMM_SIZE(mpi_comm,nproc,ierr)
!!!  parallel=nproc>1
!!!
!!!  nodename=f_malloc_str(MPI_MAX_PROCESSOR_NAME,0.to.nproc-1,id='nodename')
!!!  timecnt=f_malloc((/1.to.ncounters,0.to.nproc/),id='timecnt')
!!!
!!!  if (parallel) then 
!!!     call MPI_COMM_RANK(mpi_comm,iproc,ierr)
!!!     call MPI_GATHER(pctimes,ncounters,MPI_DOUBLE_PRECISION,&
!!!          timecnt,ncounters,MPI_DOUBLE_PRECISION,0,mpi_comm,ierr)
!!!     if (debugmode) then
!!!        !initalise nodenames
!!!        do jproc=0,nproc-1
!!!           nodename(jproc)=repeat(' ',MPI_MAX_PROCESSOR_NAME)
!!!        end do
!!!
!!!        call MPI_GET_PROCESSOR_NAME(nodename_local,namelen,ierr)
!!!
!!!        !gather the result between all the process
!!!        call MPI_GATHER(nodename_local,MPI_MAX_PROCESSOR_NAME,MPI_CHARACTER,&
!!!             nodename(0),MPI_MAX_PROCESSOR_NAME,MPI_CHARACTER,0,&
!!!             mpi_comm,ierr)
!!!     end if
!!!
!!!  else
!!!     do i=1,ncounters
!!!        timecnt(i,0)=pctimes(i)
!!!     end do
!!!     iproc=0
!!!  endif
!!!
!!!  if (iproc == 0) then
!!!!!$     !synthesis of the counters
!!!!!$     call timing_data_synthesis(nproc,ncounters,timecnt,timecnt(1,nproc))
!!!!!$
!!!!!$     call timing_open_stream(iunit_def)
!!!!!$     call yaml_mapping_open('SUMMARY',advance='no')
!!!!!$     call yaml_comment('     % ,  Time (s)',tabbing=tabfile)
!!!!!$     
!!!!!$     !sum all the information by counters
!!!!!$     do i=1,ncounters
!!!!!$        pc=100.d0*timecnt(i,nproc)/sum(timecnt(1:ncounters,nproc))
!!!!!$        call timing_dump_line(trim(pcnames(i)),tabfile,pc,timecnt(i,nproc))
!!!!!$     end do
!!!!!$     call timing_dump_line('Total',tabfile,100.d0,sum(timecnt(1:ncounters,nproc)))
!!!!!$     call yaml_mapping_close() !summary
!!!
!!!     call dict_init(dict_info)
!!!     !here this information can be dumped by adding an extra dictionary to the routine arguments
!!!     !call yaml_mapping_open('CPU parallelism')
!!!     !call yaml_map('MPI_tasks',nproc)
!!!     nthreads = 0
!!!     !$  nthreads=omp_get_max_threads()
!!!     !if (nthreads /= 0) call yaml_map('OMP threads',nthreads)
!!!     !call yaml_mapping_close()
!!!     call set(dict_info//'CPU parallelism'//'MPI tasks',nproc)
!!!     if (nthreads /= 0) call set(dict_info//'CPU parallelism'//'OMP threads',&
!!!          nthreads)
!!!
!!!     if (debugmode .and. parallel) then
!!!        !call yaml_sequence_open('Hostnames')
!!!        !do jproc=0,nproc-1
!!!        !   call yaml_sequence(trim(nodename(jproc)))
!!!        !end do
!!!        !call yaml_sequence_close()
!!!        call set(dict_info//'Hostnames',&
!!!             list_new(.item. nodename))
!!!     end if
!!!
!!!     call timing_dump_counters(ncounters,nproc,pcnames,timecnt,dict_info)
!!!
!!!!!$     !dump extra info dictionary
!!!!!$     if (associated(dict_info)) call yaml_dict_dump(dict_info)
!!!!!$     call yaml_map('Report timestamp',trim(yaml_date_and_time_toa()))
!!!!!$     !restore the default stream
!!!!!$     call timing_close_stream(iunit_def)
!!!     call dict_free(dict_info)
!!!  end if
!!!  call f_free(timecnt)
!!!  call f_free_str(MPI_MAX_PROCESSOR_NAME,nodename)
!!!  
!!!  call f_release_routine()
!!!end subroutine sum_counters
!!!
!!!
!!!subroutine sum_results(ncat,mpi_comm,message,timesum)
!!!  use dynamic_memory
!!!  use yaml_output
!!!  use time_profiling, only: timing_dump_results
!!!  implicit none
!!!  include 'mpif.h'
!!!  integer, intent(in) :: mpi_comm,ncat
!!!  character(len=*), intent(in) :: message
!!!  real(kind=8), dimension(ncat+1), intent(inout) :: timesum
!!!   !local variables
!!!  integer :: i,ierr,j,icls,icat,jproc,iextra,iproc,iunit_def,nproc
!!!  integer, dimension(ncat) :: isort
!!!  real(kind=8), dimension(:,:), allocatable :: timeall
!!!
!!!  ! Not initialised case.
!!!  if (mpi_comm==MPI_COMM_NULL) return
!!!  call f_routine(id='sum_results')
!!!
!!!  call MPI_COMM_SIZE(mpi_comm,nproc,ierr)
!!!  !allocate total timings
!!!  timeall=f_malloc((/1.to.ncat+1,0.to.nproc/),id='timeall')
!!!
!!!  if (nproc>1) then
!!!     call MPI_COMM_RANK(mpi_comm,iproc,ierr)
!!!     call MPI_GATHER(timesum,ncat+1,MPI_DOUBLE_PRECISION,&
!!!          timeall,ncat+1,MPI_DOUBLE_PRECISION,0,mpi_comm,ierr)
!!!  else
!!!     do i=1,ncat+1
!!!        timeall(i,0)=timesum(i)
!!!     end do
!!!     iproc=0
!!!  endif
!!!
!!!  if (iproc == 0) then
!!!     call timing_dump_results(ncat,nproc,message,timeall)
!!!  endif
!!!  call f_free(timeall)
!!!  call f_release_routine()
!!!
!!!END SUBROUTINE sum_results
!!!

!> interrupts all timing activities to profile the category indicated by cat_id
!! see e.g. http://en.wikipedia.org/wiki/Interrupt 
!! The action is then finished by calling the routine f_timer_resume
subroutine f_timer_interrupt(cat_id)
  use time_profiling, only: f_timing
  implicit none
  integer, intent(in) :: cat_id

  call f_timing(cat_id,'IR')
end subroutine f_timer_interrupt

!>restore the previous status of the timer and stop the counter of the interruption
subroutine f_timer_resume()
  use time_profiling, only: f_timing
  implicit none
  
  call f_timing(1,'RX')
end subroutine f_timer_resume

subroutine sort_positions(n,a,ipiv)
  implicit none
  integer, intent(in) :: n
  real(kind=8), dimension(n), intent(in) :: a
  integer, dimension(n), intent(out) :: ipiv
  !local variables
  integer :: i,j,jmax,imax
  real(kind=8) :: locmax

  !neutral permutation
  do i=1,n
     ipiv(i)=i
  end do
  !find the order for all the arrays
  do j=1,n
  !search the maximum
     locmax=-1.d300
     do i=j,n
        if (locmax < a(ipiv(i))) then
           locmax=a(ipiv(i))
           jmax=ipiv(i)
           imax=i
        end if
     end do
     !swap the position with j
     ipiv(imax)=ipiv(j) !throw in the present element
     ipiv(j)=jmax       !take out the present maximum
  end do
  !do i=1,n
  !   print *,'a',i,a(i),ipiv(i),a(ipiv(i))
  !end do
  !stop
end subroutine sort_positions

