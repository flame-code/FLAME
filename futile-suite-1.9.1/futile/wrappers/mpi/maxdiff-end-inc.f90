!> @file
!! Include fortran file for end maxdiff operations
!! @author
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS

  if (srce == -1) then
     if ( irank == iroot) then
        do jproc=2,nproc
           do i=1,ndims
              maxdiff=max(maxdiff,&
                   abs(array_glob(i,jproc)-array_glob(i,1)))
           end do
        end do
     end if
  else
     do i=1,ndims
        maxdiff=max(maxdiff,&
             abs(array_glob(i,2)-array_glob(i,1)))
     end do
  end if

  call f_free(array_glob)

  !in case of broadcasting the difference should be known by everyone
  if (bcst .and. srce==-1) then
     !never put check =.true. here, otherwise stack overflow as recursion occurs
     call fmpi_bcast(maxdiff,1,root=iroot,comm=mpi_comm,check=.false.)
     call fmpi_barrier(mpi_comm) !redundant, but just in case
  else if(srce >=0) then
     call fmpi_allreduce(maxdiff,1,FMPI_MAX,comm=mpi_comm)
  end if
