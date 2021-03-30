!> @file
!! Include fortran file for declaration scatterv operations
!! @author
!!    Copyright (C) 2015-2015 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
  integer :: recvcount  !intent in for scalars    
  integer, dimension(:), intent(in) :: displs,sendcounts
  integer, intent(in), optional :: root,comm
  !local variables
  integer :: mpi_comm,root_,ierr
  external :: MPI_SCATTERV
