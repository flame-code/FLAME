!> @file
!! Include fortran file for memory remapping
!! @author
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
  integer(f_long), intent(in) :: lb,lu
  integer(f_long), intent(inout) :: pos,wsz
  !local variables
  integer(f_long) :: szm1
  szm1=lu-lb
  !check if the position and the sizes are compatible with the 
  !allocation of the pointer, otherwise resize the pointer
  if (wsz > pos+szm1) then
     wsz=2*wsz
     wtmp=>work
     nullify(work)
     work = f_malloc_ptr(wsz,id='work_d')
     call f_memcpy(src=wtmp,dest=work)
     call f_free_ptr(wtmp)
  end if
  call f_map_ptr(lb,lu,work(pos:pos+szm1),ptr)
  !increment the position
  pos=pos+szm1+1
