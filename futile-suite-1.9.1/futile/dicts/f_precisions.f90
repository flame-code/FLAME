!> @file
!! Define precisions in portable format for future usage
!! @author
!!    Copyright (C) 2014-2015 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
module f_precisions
  !This module enhances the portability of various kind of variables
  !and defines other objects that might be used in the host code
  implicit none

  public

  !for reals and complex, to be verified if supported
  integer, parameter :: f_simple = selected_real_kind(6, 37) !simple precision specification, to be used for real variables
  integer, parameter :: f_double = selected_real_kind(15, 307) !double precision for real 
  integer, parameter :: f_quadruple = selected_real_kind(33, 4931) !parameter indicating the quadrupole precision, if supported by the fortran processor

  !for integers to be verified if they are supported, especially long
  integer, parameter :: f_short=selected_int_kind(4) !16-bit integer
  integer, parameter :: f_integer=selected_int_kind(8) !32-bit integer, usual precision for default integer
  integer, parameter :: f_long=selected_int_kind(16) !64-bit integer, usually needed for pointer adresses on modern machines

  !logicals to be done also, and tested against bits and bytes with f_loc
  !integer, parameter :: bit=0 !not supported
  integer, parameter :: f_byte=1 !precision for 1-byte logical, to be declared as logical(f_byte)

  integer, parameter :: f_address = f_long ! kind of a pointer address, long integer on most machines. Might be detected at configure time

  character(len=*), parameter :: f_cr=char(13)//char(10) !portable carriage return, contains both CR for unix and DOS

  character(len=*), parameter :: f_backslash=char(92)  ! portable backslash character  as some compilers do not appreciate it even in strings

  !>characters for the equality and the association of the parameters
  integer, parameter, private :: num_size=4
  character(len=num_size), parameter, private :: c_0='zero'
  character(len=num_size), parameter, private :: c_1='one '

  type, public :: f_parameter
     !type for the definition of a f_parameter
     character(len=num_size) :: val
  end type f_parameter

  type(f_parameter), parameter :: f_0=f_parameter(c_0)
  type(f_parameter), parameter :: f_1=f_parameter(c_1)

  !bytes representation of true and false
  logical(f_byte), parameter :: f_T=.true._f_byte
  logical(f_byte), parameter :: f_F=.false._f_byte

  !> function to localize the address of anything
  integer(f_address), external :: f_loc

  interface assignment(=)
     module procedure set_d,set_i,set_li,set_c
  end interface assignment(=)

  interface f_sizeof
     module procedure sizeof_r,sizeof_d,sizeof_i,sizeof_l,sizeof_li
  end interface f_sizeof

  private :: set_d,set_i,set_li,sizeof_r,sizeof_d,sizeof_i,sizeof_l,sizeof_li
  
  contains

    !> safe conversion from double to simple, avoid overflow and underflows
    elemental pure function f_simplify(d) result(r)
      implicit none
      real(f_double), intent(in) :: d
      real(f_simple) :: r
      !local variables
      real(f_double), parameter :: hg=3.4028234663852886E+38_f_double!real(huge(1.0_f_simple),f_double)
      real(f_double), parameter :: tn=1.17549435082228751E-38_f_double!real(tiny(1.0_f_simple),f_double)

      if (d > hg) then
         r=huge(1.0_f_simple)
      else if (d < tn) then
         r=tiny(1.0_f_simple)
      else
         r=real(d,f_simple)
      end if

    end function f_simplify


    function sizeof_r(av) result(k)
      !function documented here
      implicit none
      real(f_simple), intent(in) :: av
      real(f_simple), dimension(2) :: a
      integer :: k
      a(1)=av
      k=int(f_loc(a(2))-f_loc(a(1)))
    end function sizeof_r

    function sizeof_d(av) result(k)
      implicit none
      real(f_double), intent(in) :: av
      real(f_double), dimension(2) :: a
      integer :: k
      a(1)=av
      k=int(f_loc(a(2))-f_loc(a(1)))
    end function sizeof_d

    function sizeof_i(av) result(k)
      implicit none
      integer(f_integer), intent(in) :: av
      integer(f_integer), dimension(2) :: a
      integer :: k
      a(1)=av
      k=int(f_loc(a(2))-f_loc(a(1)))
    end function sizeof_i

    function sizeof_li(av) result(k)
      implicit none
      integer(f_long), intent(in) :: av
      integer(f_long), dimension(2) :: a
      integer :: k
      a(1)=av
      k=int(f_loc(a(2))-f_loc(a(1)))
    end function sizeof_li
    
    function sizeof_l(av) result(k)
      implicit none
      logical, intent(in) :: av
      logical, dimension(2) :: a
      integer :: k
      a(1)=av
      k=int(f_loc(a(2))-f_loc(a(1)))
    end function sizeof_l

    elemental subroutine set_d(val,par)
      implicit none
      real(f_double), intent(out) :: val
      type(f_parameter), intent(in) :: par
      select case(par%val)
      case(c_0)
         val=0.0_f_double
      case(c_1)
         val=1.0_f_double
      end select
    end subroutine set_d

    elemental subroutine set_c(val,par)
      implicit none
      character(len=*), intent(out) :: val
      type(f_parameter), intent(in) :: par
      select case(par%val)
      case(c_0)
         val=repeat(' ',len(val))
      end select
    end subroutine set_c


    elemental subroutine set_i(val,par)
      implicit none
      integer(f_integer), intent(out) :: val
      type(f_parameter), intent(in) :: par
      select case(par%val)
      case(c_0)
         val=int(0,f_integer)
      case(c_1)
         val=int(1,f_integer)
      end select
    end subroutine set_i

    elemental subroutine set_li(val,par)
      implicit none
      integer(f_long), intent(out) :: val
      type(f_parameter), intent(in) :: par
      select case(par%val)
      case(c_0)
         val=int(0,f_long)
      case(c_1)
         val=int(1,f_long)
      end select
    end subroutine set_li


end module f_precisions

!> Function which identifies the address of the scalar object
!! associated to a unknown quantity
function f_loc(routine)
  use f_precisions, only: f_address
  implicit none
  external :: routine       !< Object
  integer(f_address) :: f_loc  !< Address of the object routine

  call getlongaddress(routine,f_loc)

end function f_loc
