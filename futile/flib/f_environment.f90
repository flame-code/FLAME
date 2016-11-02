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

  contains

    subroutine f_environment_acquire()
      implicit none
      !local variables
      integer :: istat,idum
      character(len=8) :: val
          
      !check if we are in the bigdebug mode or not
      call get_environment_variable('FUTILE_DEBUG_MODE',&
           val,status=istat)
      !take the debug level
      if (istat==0) then
         read(val,*)f_debug_level
         f_debug=f_debug_level /=0
      end if

      call get_environment_variable('FUTILE_SIMULATION_MODE',&
           val,status=istat)
      if (istat==0) then
         read(val,*)idum
         f_simulation_mode=idum ==1
      end if

      call get_environment_variable('FUTILE_MALLOC_NAN_PADDING',&
           val,status=istat)
      if (istat==0) then
         read(val,*)idum
         f_nan_pad_size=idum
      end if

      call get_environment_variable('FUTILE_MALLOC_MEMORYLIMIT',&
           val,status=istat)
      if (istat==0) then
         read(val,*)idum
         f_memorylimit=idum
      end if

    end subroutine f_environment_acquire

end module f_environment
