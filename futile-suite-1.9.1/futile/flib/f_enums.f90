!> @file
!! Define operations on enumerators
!! @author
!!    Copyright (C) 2014-2015 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
!! @usage
!!type(f_enumerator) :: primates=f_enumerator('PRIMATE',PRIMATE_INT,null()) !define the concept of primates
!!type(f_enumerator) :: monkey=f_enumerator('MONKEY',MONKEY_INT,null()) !define the monkeys
!!type(f_enumerator) :: man=f_enumerator('MAN',MAN_INT,null()) !define the man
!!
!!!the man is also a primate
!!call f_enum_attr(man,primates)
!!!as well as the monkey
!!call f_enum_attr(monkey,primates)
!!
!!!we now have as a variable the following enumerator:
!!type(f_enumerator) :: specimen
!!!which have been inutuialized somehow
!!
!!!usual usage of a flag
!!select case(str(speciment))
!!case('MAN')
!!!etc etc
!!case('MONKEY')
!!!etc !etc
!!end select
!!
!!!otherwise we might define an action that is executed either for monkey and man
!!if (specimen .hasattr. primates) then
!!!this source section is executed either for man and monkey case
!!end if
module f_enums
  implicit none

  private

  integer, parameter, private :: NULL_INT=-1024
  character(len=*), parameter, private :: null_name='nullified enumerator'

  !> enumerator type, useful to define different modes
  type, public :: f_enumerator
     character(len=64) :: name=null_name
     integer :: id=null_int
     type(f_enumerator), pointer :: family => null()
  end type f_enumerator

  type(f_enumerator), parameter, private :: &
       f_enum_null=f_enumerator(null_name,NULL_INT,null())

  interface operator(==)
     module procedure enum_is_enum,enum_is_char,enum_is_int, int_is_enum
  end interface operator(==)

  interface operator(/=)
     module procedure enum_is_not_enum,enum_is_not_char,enum_is_not_int
  end interface operator(/=)

  interface operator(.hasattr.)
     module procedure enum_has_attribute,enum_has_char,enum_has_int
  end interface operator(.hasattr.)

  interface operator(.getattr.)
     module procedure enum_get_from_int,enum_get_from_enum,enum_get_from_char
  end interface operator(.getattr.)


  interface toi
     module procedure int_enum
  end interface toi

  interface toa
     module procedure char_enum
  end interface toa

  public :: f_enum_attr,operator(.hasattr.),nullify_f_enum,f_enum_update
  public :: toi,toa,f_enumerator_null,operator(==),operator(/=),operator(.getattr.)

contains

  pure subroutine nullify_f_enum(en)
    implicit none
    type(f_enumerator), intent(out) :: en
    en=f_enumerator_null()
  end subroutine nullify_f_enum

  pure function f_enumerator_null() result(en)
    use yaml_strings, only: f_strcpy
    implicit none
    type(f_enumerator) :: en
    call f_strcpy(src=null_name,dest=en%name)
    en%id=null_int
    nullify(en%family)
  end function f_enumerator_null

  !>associate the enumerator to a family
  subroutine f_enum_attr(dest,attr)
    implicit none
    type(f_enumerator), intent(inout), target :: attr !< to avoid the parameter attribute
    type(f_enumerator), intent(inout) :: dest
    !local variables
    type(f_enumerator), pointer :: iter
    !print *,'ASSOCIATING',trim(char(dest)),trim(char(attr))
    !search the first iterator that has no family
    if (.not. associated(dest%family)) then
       !print *,'here'
       dest%family=>attr
    else
       !print *,'there'
       !if the enumerator already exists do not associate
       iter => dest%family
       do while(associated(iter%family) .and. (iter /= attr))
          iter => iter%family
       end do
       if (iter /= attr) iter%family=>attr
    end if
  end subroutine f_enum_attr

  function enum_get_from_int(en,family) result(iter)
    implicit none
    integer, intent(in) :: family

    type(f_enumerator), intent(in) :: en
    type(f_enumerator), pointer :: iter
    !local variables
    logical :: ok
    ok=.false.
    iter=>en%family
    do while(associated(iter))
       ok= iter == family
       if (ok) exit
       iter => iter%family
    end do

  end function enum_get_from_int

  function enum_get_from_char(en,family) result(iter)
    implicit none
    character(len=*), intent(in) :: family

    type(f_enumerator), intent(in) :: en
    type(f_enumerator), pointer :: iter
    !local variables
    logical :: ok
    ok=.false.
    iter=>en%family
    do while(associated(iter))
       ok= iter == family
       if (ok) exit
       iter => iter%family
    end do

  end function enum_get_from_char

  function enum_get_from_enum(en,family) result(iter)
    implicit none
    type(f_enumerator), intent(in) :: family

    type(f_enumerator), intent(in) :: en
    type(f_enumerator), pointer :: iter
    !local variables
    logical :: ok
    ok=.false.
    iter=>en%family
    do while(associated(iter))
       ok= iter == family
       if (ok) exit
       iter => iter%family
    end do

  end function enum_get_from_enum

  function enum_has_attribute(en,family) result(ok)
    implicit none
    type(f_enumerator), intent(in) :: family
    type(f_enumerator), intent(in) :: en
    logical :: ok

    ok =associated(en .getattr. family)
!!$    !local variables
!!$    type(f_enumerator), pointer :: iter
!!$
!!$    ok=.false.
!!$    iter=>en%family
!!$    do while(associated(iter) .and. .not. ok)
!!$       ok= iter == family
!!$       iter => iter%family
!!$    end do
  end function enum_has_attribute

  !> used in the assignment. Change the basic enumerator but combine hte attributes of the previous one 
  !! with the attributes of the present one
  subroutine f_enum_update(dest,src)
    implicit none
    type(f_enumerator), intent(in) :: src
    type(f_enumerator), intent(inout) :: dest
    !local variables
    type(f_enumerator), pointer :: iter

    !start first by copying the attributes of the first enumerator
    iter=>src%family
    do while(associated(iter))
       call f_enum_attr(dest,iter)
       iter => iter%family
    end do
    !then scratch the root of the enumerator with src
    dest%name=src%name
    dest%id=src%id
  end subroutine f_enum_update

  function enum_has_char(en,family) result(ok)
    implicit none
    character(len=*), intent(in) :: family
    type(f_enumerator), intent(in) :: en
    logical :: ok

    ok =associated(en .getattr. family)
!!$    !local variables
!!$    type(f_enumerator), pointer :: iter
!!$
!!$    ok=.false.
!!$    iter=>en%family
!!$    do while(associated(iter) .and. .not. ok)
!!$       ok= iter == family
!!$       iter => iter%family
!!$    end do
  end function enum_has_char

  function enum_has_int(en,family) result(ok)
    implicit none
    integer, intent(in) :: family
    type(f_enumerator), intent(in) :: en
    logical :: ok

    ok =associated(en .getattr. family)

!!$    !local variables
!!$    type(f_enumerator), pointer :: iter
!!$
!!$    ok=.false.
!!$    iter=>en%family
!!$    do while(associated(iter) .and. .not. ok)
!!$       ok= iter == family
!!$       iter => iter%family
!!$    end do
  end function enum_has_int

  elemental pure function enum_is_enum(en,en1) result(ok)
    implicit none
    type(f_enumerator), intent(in) :: en
    type(f_enumerator), intent(in) :: en1
    logical :: ok
    ok = (en == en1%id) .and. (en == en1%name)
  end function enum_is_enum

  elemental pure function enum_is_int(en,int) result(ok)
    implicit none
    type(f_enumerator), intent(in) :: en
    integer, intent(in) :: int
    logical :: ok
    ok = en%id == int
  end function enum_is_int

  elemental pure function int_is_enum(int,en) result(ok)
    implicit none
    type(f_enumerator), intent(in) :: en
    integer, intent(in) :: int
    logical :: ok
    ok = en%id == int
  end function int_is_enum
  
  elemental pure function enum_is_char(en,char) result(ok)
    use yaml_strings, only: operator(.eqv.)
    implicit none
    type(f_enumerator), intent(in) :: en
    character(len=*), intent(in) :: char
    logical :: ok
    ok = trim(en%name) .eqv. trim(char)
  end function enum_is_char

  elemental pure function enum_is_not_enum(en,en1) result(ok)
    implicit none
    type(f_enumerator), intent(in) :: en
    type(f_enumerator), intent(in) :: en1
    logical :: ok
    ok = .not. (en == en1)
  end function enum_is_not_enum

  elemental pure function enum_is_not_int(en,int) result(ok)
    implicit none
    type(f_enumerator), intent(in) :: en
    integer, intent(in) :: int
    logical :: ok
    !ok = .not. (en == int)
    ok = .not. enum_is_int(en,int)
  end function enum_is_not_int

  elemental pure function enum_is_not_char(en,char) result(ok)
    implicit none
    type(f_enumerator), intent(in) :: en
    character(len=*), intent(in) :: char
    logical :: ok
    ok = .not. (en == char)
  end function enum_is_not_char


  !>integer of f_enumerator type.
  elemental pure function int_enum(en)
    type(f_enumerator), intent(in) :: en
    integer :: int_enum
    int_enum=en%id
  end function int_enum

  !>char of f_enumerator type.
  elemental pure function char_enum(en)
    type(f_enumerator), intent(in) :: en
    character(len=len(en%name)) :: char_enum
    char_enum=en%name
  end function char_enum


end module f_enums
