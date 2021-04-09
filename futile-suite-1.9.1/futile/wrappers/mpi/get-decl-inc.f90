!> @file
!! Include fortran file for mpi_get

!! @author
!!    Copyright (C) 2017-2017 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
  integer, intent(in) :: count
  integer(fmpi_integer), intent(in) :: target_rank
  type(fmpi_win), intent(in) :: win
  integer(fmpi_address), intent(in), optional :: target_disp
  integer, intent(in), optional :: origin_displ
  ! Local variables
  integer :: from
  integer(fmpi_integer) :: ierr
  integer(fmpi_address) :: tdispl
  external :: MPI_GET
  from=1
  if(present(origin_displ)) from=origin_displ+1
