!> @file
!!  Define the get_command_argument
!! @author Luigi Genovese
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Define the get_command_argument routine
subroutine get_command_argument(count, value, status)
  implicit none
  integer, intent(in) :: count
  integer, intent(out) :: status
  character(len = *), intent(out) :: value

  write(value, "(A)") ""
  status = 1
end subroutine get_command_argument
