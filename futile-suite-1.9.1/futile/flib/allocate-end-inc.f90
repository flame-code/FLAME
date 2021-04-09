!> @file
!! Include fortran file for allocation template
!! 
!! @author
!!    Copyright (C) 2012-2015 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
  !profile the array allocation
  iadd=int(0,f_address)
  !write the address of the first element in the address string
  if (m%profile .and. track_origins) iadd=loc_arr(array)!call getlongaddress(array,iadd)

  call f_update_database(product(int(m%shape(1:m%rank),kind=8)),kind(array),m%rank,&
       iadd,m%array_id,m%routine_id,m%info)

  call f_timer_resume()!TCAT_ARRAY_ALLOCATIONS
