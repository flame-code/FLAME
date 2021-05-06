!> @file
!! Include fortran file for maxdiff declarations
!! @author
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS

  !> if .true. all the proc will have the same maxdiff (default .false.)
  !! if source is present, the array is overwritten with the value coming
  !! from rank source
  logical, intent(in), optional :: bcast 
  !> determine the rank of the processor whose data are used for comparison.
  !! Implies bcast=.true. as at the end all the differences are broadcasted
  integer, intent(in), optional :: source 
  integer, intent(in), optional :: root !<rank of the process retrieving the diff
  integer, intent(in), optional :: comm !<communicator
  !local variables
  logical :: bcst
  integer :: ndims,nproc,mpi_comm,iroot,i,jproc,srce,irank

  bcst=.false.
  if (present(bcast)) bcst=bcast
  srce=-1
  if (present(source)) then
     srce=source
  end if
  mpi_comm=fmpi_comm(comm)
  if (present(root)) then
     iroot=root
  else
     iroot=0
  end if
  nproc=mpisize(mpi_comm)
  irank=mpirank(mpi_comm)
