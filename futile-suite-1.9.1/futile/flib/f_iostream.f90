!> @file
!! Manage different low-level operations on io streams
!! like operations on external files and basic operations in memory
!! @author
!!    Copyright (C) 2012-2015 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
module f_iostream
  use dictionaries
  implicit none

  private

  !> This type can be used to get strings from a file or a dictionary long string.
  !! it should be amended to also contain the same information from a character buffer
  type, public :: io_stream
     integer :: iunit = 0
     type(dictionary), pointer :: lstring => null()
  end type io_stream

  public :: f_iostream_from_file,f_iostream_from_lstring
  public :: f_iostream_get_line,f_iostream_release


contains

  subroutine f_iostream_from_file(ios, filename)
    use yaml_strings
    use f_utils
    implicit none
    type(io_stream), intent(out) :: ios
    character(len = *), intent(in) :: filename
    !Local variables
    integer :: ierror

    ios%iunit=f_get_free_unit(742)
    open(unit=ios%iunit,file=trim(filename),status='old',iostat=ierror)
    !Check the open statement
    if (ierror /= 0) call f_err_throw('Error in opening file='//&
         trim(filename)//' iostat='//trim(yaml_toa(ierror)),&
         err_name='INPUT_OUTPUT_ERROR')
    nullify(ios%lstring)
  end subroutine f_iostream_from_file

  subroutine f_iostream_from_lstring(ios, dict)
    implicit none
    type(io_stream), intent(out) :: ios
    type(dictionary), pointer :: dict

    ios%iunit = 0
    if (dict_len(dict) < 0) call f_err_throw('Error dict is not a long string',&
         err_name='INPUT_OUTPUT_ERROR')
    ios%lstring => dict
  end subroutine f_iostream_from_lstring

  subroutine f_iostream_get_line(ios, line, eof)
    implicit none
    !Arguments
    type(io_stream), intent(inout) :: ios
    character(len=max_field_length), intent(out) :: line
    logical, optional, intent(out) :: eof
    !Local variables
    integer :: i_stat
    character(len=8) :: fmt

    if (ios%iunit > 0) then
       write(fmt, "(A,I0,A)") "(A", max_field_length, ")"
       if (present(eof)) then
          read(ios%iunit, trim(fmt), iostat = i_stat) line
       else
          read(ios%iunit, trim(fmt)) line
          i_stat = 0
       end if
       if (i_stat /= 0) then
          close(ios%iunit)
          ios%iunit = 0
       end if
       if (present(eof)) eof = (i_stat /= 0)
    else if (associated(ios%lstring)) then
       if (dict_len(ios%lstring) > 0) ios%lstring => dict_iter(ios%lstring)
       line = dict_value(ios%lstring)
       ios%lstring => dict_next(ios%lstring)
       if (present(eof)) eof = .not. associated(ios%lstring)
    else if (present(eof)) then
       eof = .true.
    end if
  end subroutine f_iostream_get_line

  subroutine f_iostream_release(ios)
    implicit none
    type(io_stream), intent(inout) :: ios

    if (ios%iunit > 0) close(ios%iunit)
    ios%iunit = 0
    nullify(ios%lstring)
  end subroutine f_iostream_release
  
!!$  !does not increment the line size
!!$  function f_iostream_get_line_len(ios)
!!$  end function f_iostream_get_line_len
!!$
!!$  function f_iostream_readline(ios,line,eor) result(eof)
!!$    type(io_stream), intent(inout) :: ios
!!$    character(len=*), intent(out) :: line
!!$    logical :: eof
!!$    !if provided, determine if the line character string is too short to
!!$    !contain the line in the length of the buffer line
!!$    !if false, it implies that the line is large enough and padded
!!$    !with blanks. 
!!$    !When true, the line is filled with blanks and cannot be used for reading.
!!$    !the function f_iostream_get_line_len should be used to determine the
!!$    !size of the new length
!!$    !if not provided, the line is filled with the available
!!$    !characters and truncated in the case the record is longer
!!$    logical, intent(out), optional :: eor
!!$    
!!$  end function f_iostream_readline
  
end module f_iostream
