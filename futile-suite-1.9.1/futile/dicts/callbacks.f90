!> @file
!! Contains definition of callback mechanisms for BigDFT
!! @author Luigi Genovese
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> @brief Module interfacing the callback to different functions
!! @details Should be used in the error handling module
!! The user of this module is able to define the stopping routine, the new error and the error codes
module exception_callbacks
  use f_precisions
  implicit none
  private

  !> Address of the generic callback functions, valid for errors with non-specific callbacks
  integer(f_address) :: callback_add=0
  integer(f_address) :: callback_data_add=0
  !> Address of the overrided severe error
  integer(f_address) :: severe_callback_add=0
  !> Addresses of the error dump callbacks functions
  integer(f_address) :: last_error_callback_add=0
  !> Addresses of the all errors error dump callbacks functions
  integer(f_address) :: all_errors_callback_add=0


  interface f_err_set_callback
     module procedure err_set_callback_simple,err_set_callback_advanced
  end interface


  public :: f_err_set_callback
  public :: f_err_unset_callback
  public :: f_err_severe, f_err_severe_override, f_err_severe_restore
  public :: f_err_ignore,f_dump_last_error,f_dump_possible_errors
  public :: f_err_set_last_error_callback
  public :: f_err_unset_last_error_callback
  public :: f_err_set_all_errors_callback
  public :: f_err_unset_all_errors_callback

  !> Internal variables for f_lib usage
  public :: callback_add
  public :: callback_data_add
  public :: severe_callback_add
  public :: err_abort


contains


  !> Subroutine which defines the way the system stops
  subroutine err_abort(callback,callback_data)
    !use metadata_interfaces
    implicit none
    integer(f_address), intent(in) :: callback,callback_data

    if (callback_data /=0 .and. callback /=0) then
       call callable_void(callback) !for the moment data are ignored
    else if (callback /=0) then
       call callable_void(callback)
    else
       call f_err_severe()
    end if
  end subroutine err_abort


  !> Defines the error routine which have to be used
  subroutine err_set_callback_simple(callback)
    implicit none
    external :: callback !< Error routine which will be called
    !$ include 'halt_omp-inc.f90'

    callback_add=f_loc(callback)
    callback_data_add=0

  end subroutine err_set_callback_simple


  subroutine err_set_callback_advanced(callback,callback_data_address)
    implicit none
    integer(f_address), intent(in) :: callback_data_address
    external :: callback
    !$ include 'halt_omp-inc.f90'

    callback_add=f_loc(callback)
    callback_data_add=callback_data_address

  end subroutine err_set_callback_advanced


  subroutine f_err_unset_callback()
    implicit none
    !!$ include 'halt_omp-inc.f90'
    callback_add=0
    callback_data_add=0
  end subroutine f_err_unset_callback


  subroutine f_err_severe_override(callback)
    implicit none
    external :: callback
    !$ include 'halt_omp-inc.f90'

    severe_callback_add=f_loc(callback)
  end subroutine f_err_severe_override


  subroutine f_err_severe_restore()
    implicit none
    !!$ include 'halt_omp-inc.f90'
    severe_callback_add=0
  end subroutine f_err_severe_restore


  !> Wrapper to ignore errors, do not dump
  subroutine f_err_ignore()
    implicit none
  end subroutine f_err_ignore


  !> Wrapper for severe errors, the can be desactivated
  subroutine f_err_severe()
    implicit none
    if (severe_callback_add == 0) then
       call f_err_severe_internal()
    else
       call callable_void(severe_callback_add)
    end if
  end subroutine f_err_severe

  !> Callback routine for severe errors
  subroutine f_err_severe_internal()
    implicit none
    call f_dump_last_error()
    stop 'Severe error, cannot proceed'
  end subroutine f_err_severe_internal

  !> Print error information about last error
  subroutine f_dump_last_error()
    implicit none
    if (last_error_callback_add /= 0) then
       call callable_void(last_error_callback_add)
    end if
  end subroutine f_dump_last_error

  subroutine f_dump_possible_errors(msg)
    implicit none
    character, dimension(*), intent(in) :: msg
    if (all_errors_callback_add /= 0) then
       call callable_arg(all_errors_callback_add,f_loc(msg))
    end if
  end subroutine f_dump_possible_errors

  !> Defines the error routine which have to be used
  subroutine f_err_set_last_error_callback(callback)
    implicit none
    external :: callback !< Error routine which will be called
    !$ include 'halt_omp-inc.f90'

    last_error_callback_add=f_loc(callback)

  end subroutine f_err_set_last_error_callback

  subroutine f_err_unset_last_error_callback()
    implicit none
    last_error_callback_add=0
  end subroutine f_err_unset_last_error_callback

  subroutine f_err_unset_all_errors_callback()
    implicit none
    all_errors_callback_add=0
  end subroutine f_err_unset_all_errors_callback

  !> Defines the error routine which have to be used
  subroutine f_err_set_all_errors_callback(callback)
    implicit none
    interface 
       subroutine callback(msg)
         implicit none
         character, dimension(*), intent(in) :: msg
       end subroutine callback
    end interface
    !$ include 'halt_omp-inc.f90'

    all_errors_callback_add=f_loc(callback)

  end subroutine f_err_set_all_errors_callback

end module exception_callbacks
