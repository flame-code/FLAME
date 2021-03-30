!> @file
!! Define ternary operators for fortran intrinsics
!! @author
!!    Copyright (C) 2016-2016 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
module f_ternary
  use f_precisions
  implicit none

  private
  
  integer, parameter :: cl=128

  type, public :: ternary_condition
     logical :: conditional
     integer(f_integer) :: i4=0
     integer(f_long) :: i8=0
     real(f_simple) :: r4=0
     real(f_double) :: r8=0
     character(len=cl) :: c0=repeat(' ',cl)
  end type ternary_condition

  interface operator(.if.)
     module procedure ternary_init
  end interface operator(.if.)

  interface operator(.then.)
     module procedure ternary_then_i0,ternary_then_d0,ternary_then_c0
  end interface operator(.then.)

  interface operator(.else.)
     module procedure ternary_else_i0,ternary_else_d0,ternary_else_c0
  end interface operator(.else.)

  public :: operator(.if.),operator(.then.),operator(.else.)

  contains

    pure function ternary_init(condition) result(cond)
      implicit none
      logical, intent(in) :: condition
      type(ternary_condition) :: cond
      cond%conditional=condition
    end function ternary_init

    pure function ternary_then_i0(cond,i4) result(th)
      implicit none
      integer(f_integer), intent(in) :: i4
      type(ternary_condition), intent(in) :: cond
      type(ternary_condition) :: th
      th%conditional=cond%conditional
      th%i4=i4
    end function ternary_then_i0

    pure function ternary_else_i0(cond,i4) result(i0)
      implicit none
      integer(f_integer), intent(in) :: i4
      type(ternary_condition), intent(in) :: cond
      integer(f_integer) :: i0
      i0=i4
      if (cond%conditional) i0=cond%i4
    end function ternary_else_i0

    pure function ternary_then_d0(cond,r8) result(th)
      implicit none
      real(f_double), intent(in) :: r8
      type(ternary_condition), intent(in) :: cond
      type(ternary_condition) :: th
      th%conditional=cond%conditional
      th%r8=r8
    end function ternary_then_d0

    pure function ternary_else_d0(cond,r8) result(d0)
      implicit none
      real(f_double), intent(in) :: r8
      type(ternary_condition), intent(in) :: cond
      real(f_double) :: d0
      d0=r8
      if (cond%conditional) d0=cond%r8
    end function ternary_else_d0

    pure function ternary_then_c0(cond,c0) result(th)
      implicit none
      character(len=*), intent(in) :: c0
      type(ternary_condition), intent(in) :: cond
      type(ternary_condition) :: th
      th%conditional=cond%conditional
      th%c0=c0
    end function ternary_then_c0

    !>this only works if the length of the results are similar
    pure function ternary_else_c0(cond,c) result(c0)
      implicit none
      character(len=*), intent(in) :: c
      type(ternary_condition), intent(in) :: cond
      character(len=len(c)) :: c0
      c0=c
      if (cond%conditional) c0=cond%c0
    end function ternary_else_c0

end module f_ternary
