!> @file
!! Include fortran file for allreduce operations

!! @author
!!    Copyright (C) 2017-2017 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
  type(f_enumerator), intent(in) :: op
  integer, intent(in), optional :: comm
  integer(fmpi_integer), intent(out), optional :: request
  !local variables
  logical :: in_place
  integer :: tcat
  integer(fmpi_integer) ::mpi_comm,ntot,ierr

  external :: MPI_ALLREDUCE,MPI_IALLREDUCE
