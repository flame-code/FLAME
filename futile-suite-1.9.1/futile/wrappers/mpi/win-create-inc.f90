!> @file
!! Include fortran file for window creation

!! @author
!!    Copyright (C) 2017-2017 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
  integer(fmpi_integer),intent(in), optional :: comm
  type(dictionary), pointer, optional :: dict_info
  type(fmpi_info), intent(in), optional :: info
  type(fmpi_win) :: win
  !local variables
  integer(fmpi_integer) :: ierr,infohandle
  type(fmpi_info) :: info_
  external :: MPI_WIN_CREATE

  infohandle=FMPI_INFO_NULL
  if (present(dict_info)) then
     call fmpi_info_create(info_)
     call fmpi_info_set(info_,dict_info)
     infohandle=info_%handle
  end if

  if (present(info)) then
     infohandle=info%handle
  end if

  !info=mpiinfo("no_locks", "true")

  win%disp_unit=mpitypesize(base)
  win%size=int(size,fmpi_address)
  win%comm=fmpi_comm(comm)
  call MPI_WIN_CREATE(base, win%size*win%disp_unit, &
       win%disp_unit, infohandle,win%comm,win%handle, ierr)

  if (ierr/=FMPI_SUCCESS) then
     call f_err_throw('Error in mpi_win_create',&
          err_id=ERR_MPI_WRAPPERS)
  end if

  if (present(dict_info)) then
     call fmpi_info_free(info_)
  end if
  !call mpiinfofree(info)

!!$      call mpi_win_fence(MPI_MODE_NOPRECEDE, window, ierr)
!!$      if (ierr/=0) then
!!$         call f_err_throw('Error in mpi_win_fence',&
!!$              err_id=ERR_MPI_WRAPPERS)
!!$      end if
