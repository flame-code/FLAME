!> @file
!! Include fortran file for mpi_get

!! @author
!!    Copyright (C) 2017-2017 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
  tdispl=int(0,fmpi_address)
  if (present(target_disp)) tdispl=target_disp
  call MPI_GET(origin_ptr,int(count,fmpi_integer),mpitype(origin_ptr),target_rank, &
       tdispl,int(count,fmpi_integer),mpitype(origin_ptr), win%handle, ierr)
  if (ierr/=FMPI_SUCCESS) call f_err_throw('Error in mpi_get',err_id=ERR_MPI_WRAPPERS)
