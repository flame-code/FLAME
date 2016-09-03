!> @file
!! Include fortran file for allocation template
!! 
!! @author
!!    Copyright (C) 2012-2015 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
  if (ierror/=0) then
     !$ if(not_omp) then
     call f_timer_resume()!TCAT_ARRAY_ALLOCATIONS
     !$ end if
     call f_err_throw('array ' // trim(m%array_id) // &
          & '(' // trim(yaml_toa(product(m%shape))) // &
          & '), error code '//trim(yaml_toa(ierror)),ERR_ALLOCATE)
     return
  end if
  if (size(shape(array))==m%rank) then
     call pad_array(array,m%put_to_zero,m%shape,ndebug)
     !also fill the array with the values of the source if the address is identified in the source
     if (m%srcdata_add > int(0,kind=8)) call c_memcopy(array,m%srcdata_add,product(shape(array))*kind(array))
     !profile the array allocation
     iadd=int(0,kind=8)
     !write the address of the first element in the address string
     if (m%profile .and. track_origins) iadd=loc_arr(array)!call getlongaddress(array,iadd)

     call f_update_database(product(int(m%shape(1:m%rank),kind=8)),kind(array),m%rank,&
          iadd,m%array_id,m%routine_id)

  else
     !$ if(not_omp) then
     call f_timer_resume()!TCAT_ARRAY_ALLOCATIONS
     !$ end if
     call f_err_throw('Rank specified by f_malloc ('+yaml_toa(m%rank)// ',routine_id=' // trim(m%routine_id) // &
          & ',id=' // trim(m%array_id) // ') is not coherent with the one of the array ('+yaml_toa(size(shape(array)))//')',&
          & ERR_INVALID_MALLOC)
     return
  end if
  !$ if(not_omp) then
  call f_timer_resume()!TCAT_ARRAY_ALLOCATIONS
  !$ end if
