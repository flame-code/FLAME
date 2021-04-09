!> @file
!! Module to load environment variables and store tem in global variables
!! It can be used bu developers for tune futile behaviour
!! as well as the verbosity of the output
!! like operations on external files and basic operations in memory
!! @author
!!    Copyright (C) 2016-2016 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
module f_environment
  implicit none

  !these variables have to be accessed by all the associating routines
  public 

  !> simple flag to see if we are in debug mode or not
  logical :: f_debug=.false.
  !> integer flag specifying the debug level
  integer :: f_debug_level=0
  !> flag to sontrol if we are in simulation mode
  logical :: f_simulation_mode=.false.
  !> ndebug integer, associated to the extra elements which are padded with NaN (whenever possible)
  !! in the last dimension of the array
  integer :: f_nan_pad_size=0
  !> memorylimit of the allocation, raise an error whenever the allocated buffers pass above a given threshold
  integer :: f_memorylimit=0
  !> maximum depth of the profiling routines (unlimited by default)
  integer :: f_maximum_profiling_depth=-1

  interface f_getenv
     module procedure getenv_i0,getenv_l0
  end interface f_getenv

  private :: getenv_i0,getenv_l0

  contains

    subroutine f_environment_acquire()
      implicit none

      call f_getenv('FUTILE_DEBUG_MODE',f_debug_level)
      f_debug = f_debug_level /= 0

      call f_getenv('FUTILE_SIMULATION_MODE',f_simulation_mode)

      call f_getenv('FUTILE_MALLOC_NAN_PADDING',f_nan_pad_size)

      call f_getenv('FUTILE_MEMORYLIMIT_GB',f_memorylimit)

      call f_getenv('FUTILE_PROFILING_DEPTH',f_maximum_profiling_depth)

    end subroutine f_environment_acquire

    subroutine getenv_i0(envvar,data)
      implicit none
      character(len=*), intent(in) :: envvar
      integer, intent(inout) :: data
      !local variables
      integer :: istat,idum
      character(len=8) :: val

      call get_environment_variable(envvar,val,status=istat)
      if (istat==0) then
         read(val,*)idum !let it crash if the reading is unsafe
         data=idum
      end if

    end subroutine getenv_i0

    subroutine getenv_l0(envvar,data)
      implicit none
      character(len=*), intent(in) :: envvar
      logical, intent(inout) :: data
      !local variables
      integer :: istat,idum
      character(len=8) :: val

      call get_environment_variable(envvar,val,status=istat)
      if (istat==0) then
         read(val,*)idum !let it crash if the reading is unsafe
         data=idum==1
      end if

    end subroutine getenv_l0



end module f_environment
