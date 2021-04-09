!> @file
!! Include fortran file for maxdiff operations with scalars
!! @author
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
  if (nproc == 1 .or. ndims == 0) return

  !choose the communication strategy according to the input variables
  if (srce ==-1) then
     !check that the positions are identical for all the processes
     array_glob=f_malloc((/ndims,nproc/),id='array_glob')
     call fmpi_gather(sendbuf=array,sendcount=ndims,recvbuf=array_glob,&
          root=iroot,comm=mpi_comm)
  else
     !check that the positions are identical for all the processes
     array_glob=f_malloc((/ndims,2/),id='array_glob')
     call f_memcpy(n=ndims,src=array,dest=array_glob(1,2))
     if (irank == srce) call f_memcpy(n=ndims,src=array,dest=array_glob(1,1))
     call fmpi_bcast(array_glob(1,1),ndims,root=srce,comm=mpi_comm)
     if (bcst .and. irank /= srce) &
          call f_memcpy(n=ndims,dest=array,src=array_glob(1,1))
  end if
  include 'maxdiff-end-inc.f90'
