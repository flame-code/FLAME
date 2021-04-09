!> @file
!! Include fortran file for allocation template
!! 
!! @author
!!    Copyright (C) 2012-2015 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
  logical :: ok
  ok=.false.
  if (ierror/=0) then
     call f_timer_resume()!TCAT_ARRAY_ALLOCATIONS
     call f_err_throw('array ' // trim(m%array_id) // &
          & '(' // trim(yaml_toa(product(m%shape(1:m%rank)))) // &
          & '), error code '//trim(yaml_toa(ierror)),ERR_ALLOCATE)
     return
  end if
  if (bigdebug .and. any(m%shape < 0)) then
     call f_err_throw('array has suspect shape (size < 0 in at least one dimension) ' // trim(m%array_id) // &
          & '(' // trim(yaml_toa(product(m%shape))) // &
          & ')',ERR_ALLOCATE)
  end if
  if (rank/=m%rank) then
     call f_timer_resume()!TCAT_ARRAY_ALLOCATIONS
     call f_err_throw('Rank specified by f_malloc ('+yaml_toa(m%rank)// ',routine_id=' // trim(m%routine_id) // &
          & ',id=' // trim(m%array_id) // ') is not coherent with the one of the array ('+yaml_toa(rank)//')',&
          & ERR_INVALID_MALLOC)
     return
  end if
  ok=.true.
