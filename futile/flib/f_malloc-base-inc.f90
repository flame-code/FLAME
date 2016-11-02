!> @file
!! Include fortran file for f_malloc routines
!! initialize the structure
!! @author
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
  logical, intent(in), optional :: profile
  character(len=*), intent(in), optional :: id,routine_id,info
  !local variables
  !integer(f_integer) :: lgt
!!$  logical :: within_openmp
!!$  !$ logical :: omp_in_parallel,omp_get_nested
!!$  within_openmp=.false.
!!$  !$ within_openmp=omp_in_parallel() .or. omp_get_nested()
!!$
!!$  if (within_openmp) then
!!$     call f_err_throw('It appears that f_malloc set of instructions is being used within OMP',&
!!$          err_id=ERR_INVALID_MALLOC)
!!$  end if

  call nullify_malloc_information(m)

  if (present(id)) then
     !lgt=min(len(id),f_malloc_namelen)
     !m%array_id(1:lgt)=id(1:lgt)
     call f_strcpy(src=id,dest=m%array_id)
  else
     !m%array_id(1:len(m%array_id))='UNKNOWN_ID'
     call f_strcpy(src='UNKNOWN_ID',dest=m%array_id)
  end if
  if (present(routine_id)) then
     !lgt=min(len(routine_id),f_malloc_namelen)
     !m%routine_id(1:lgt)=routine_id(1:lgt)
     call f_strcpy(src=routine_id,dest=m%routine_id)
  else
     !m%routine_id(1:len(m%routine_id))=f_malloc_routine_name!mems(ictrl)%present_routine
     call f_strcpy(src=f_malloc_routine_name,dest=m%routine_id)
  end if

  if(present(profile)) m%profile=profile

  if(present(info)) m%info=info

