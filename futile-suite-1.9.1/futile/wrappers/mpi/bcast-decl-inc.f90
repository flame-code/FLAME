!> @file
!! Include fortran file for broadcast operations
!! declarations for scalars
!! @author
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
  integer, intent(in), optional :: count 
  integer, intent(in), optional :: root  
  integer, intent(in), optional :: comm  
  logical, intent(in), optional :: check
  !local variables
  logical chk
  integer :: n,iroot,mpi_comm,ierr
  integer, dimension(3) :: iarg_check
  external :: MPI_BCAST

  chk=.false.
  n=1
  if (present(count)) n=count
  if (present(maxdiff)) then
     call f_zero(maxdiff)
     array_diff=f_malloc(n,id='array_diff')
     call f_memcpy(n=n,src=buffer,dest=array_diff)
  end if
