!> @file
!! Include fortran file for deallocation template
!! @author
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
  if(present(shared)) then
     if(shared) then 
        !$ if (not_omp) then
        call f_timer_interrupt(TCAT_ARRAY_ALLOCATIONS)
        !$ end if

        !here the size should be corrected with ndebug (or maybe not)
        ilsize=product(int(shape(array),kind=8))
        !retrieve the address of the first element if the size is not zero
        !iadd=int(0,kind=8)
        !if (ilsize /= int(0,kind=8)) 
        iadd=loc_arr(array)!call getlongaddress(array,iadd)
        !fortran deallocation
        !LG: this in not standard and creates memory leak as
        !the fortran pointer descriptor is never deallocated
        !I am forced to comment it out as iso_c_binding does not provide 
        !a good alternative. We should find a workaround
        call bindfree(iadd)
        !call smpi_shared_free(c_loc(c_add))

        call f_purge_database(ilsize,kind(array),iadd)

        !$ if (not_omp) then
        call f_timer_resume()!TCAT_ARRAY_ALLOCATIONS
        !$ end if
        nullify(array)
        return
     end if
  end if
