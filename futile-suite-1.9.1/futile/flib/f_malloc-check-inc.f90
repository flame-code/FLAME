!> @file
!! Other fortran file for f_malloc routines
!! @author
!!    Copyright (C) 2012-2015 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS

  if (present(src) .and. present(src_ptr)) then
     call f_err_throw('Source arrays "src" and "src_ptr" cannot be present'//&
          ' at the same time in a heap allocation instance',ERR_INVALID_MALLOC)
  end if
  if (present(src_ptr)) then
     if ((present(lbounds) .or. present(ubounds) .or. present(sizes))) &
          call f_err_throw(&
          'The presence of lbounds, ubounds or sizes arrays is forbidden when src_ptr is present',&
          ERR_INVALID_MALLOC)
     include 'f_malloc-null-ptr-inc.f90'
     include 'f_malloc-ptr-inc.f90'
     !when src_ptr is given there is no need anymore to continue the routine
     return
  else if (present(src)) then
     if (present(sizes)) call f_err_throw(&
          'The presence of "sizes" array is forbidden when "src" is present',&
          ERR_INVALID_MALLOC)
     include 'f_malloc-inc.f90'
     !when src is given there is no need anymore to continue the routine
     return
  end if
