  !> @file
  !! Define the reference counter and associated methods
  !! @author
  !!    Copyright (C) 2012-2015 BigDFT group
  !!    This file is distributed under the terms of the
  !!    GNU General Public License, see ~/COPYING file
  !!    or http://www.gnu.org/copyleft/gpl.txt .
  !!    For the list of contributors, see ~/AUTHORS
module f_refcnts
  use dictionaries
  implicit none

  private

  !> Reference counter. Can be used to control the pointer
  !! referencing to a derived datatype
  type, public :: f_reference_counter
     !> Counter of references. When nullified or zero, 
     !! the associated object is ready to be destroyed
     integer, pointer :: iref =>null()
     !> Information about the associated object
     type(dictionary), pointer :: info =>null()
  end type f_reference_counter

  integer, save :: ERR_REFERENCE_COUNTERS

  !reference counters
  public :: f_ref_new,f_ref_null,f_unref,f_ref_free,f_ref_associate
  public :: nullify_f_ref,f_ref,f_ref_count,f_associated

  !for internal f_lib usage
  public :: refcnts_errors

contains

 subroutine refcnts_errors()
   implicit none
   call f_err_define(err_name='ERR_REFERENCE_COUNTERS',&
        err_msg='Error in the usage of the reference counter',&
        err_id=ERR_REFERENCE_COUNTERS,&
        err_action='When a reference counter is present each pointer association should be tracked, check for it')
 end subroutine refcnts_errors

 pure function f_associated(f_ref)
   implicit none
   type(f_reference_counter), intent(in) :: f_ref
   logical :: f_associated
   f_associated=associated(f_ref%iref)
 end function f_associated
!!!reference counter objects
  pure function f_ref_null() result(f_ref)
    implicit none
    type(f_reference_counter) :: f_ref
    call nullify_f_ref(f_ref)
  end function f_ref_null
  pure subroutine nullify_f_ref(f_ref)
    implicit none
    type(f_reference_counter), intent(out) :: f_ref
    nullify(f_ref%iref)
    nullify(f_ref%info)
  end subroutine nullify_f_ref

  subroutine dump_ref_cnt(f_ref)
    use yaml_output, only: yaml_dict_dump,yaml_map
    implicit none
    type(f_reference_counter), intent(in) :: f_ref

    call yaml_dict_dump(f_ref%info)
    call yaml_map('References',f_ref%iref)
  end subroutine dump_ref_cnt

  !> Allocate a reference counter
  !! this function should be called whe the associated object starts
  !! to be non-trivial
  function f_ref_new(id,address) result(f_ref)
    use f_precisions, only: f_address
    implicit none
    character(len=*), intent(in) :: id
    integer(f_address), intent(in), optional :: address
    type(f_reference_counter) :: f_ref
    !local variables
    type(dictionary), pointer :: dict

    call nullify_f_ref(f_ref)
    allocate(f_ref%iref)
    call dict_init(f_ref%info)
    dict=>f_ref%info//'Reference counter'
    call set(dict//'Id',id)
    if (present(address)) &
         call set(dict//'Address of referenced object',address)
    f_ref%iref=1
  end function f_ref_new

  !> Dereferencing counter
  subroutine f_unref(f_ref,count)
    implicit none
    type(f_reference_counter), intent(inout) :: f_ref
    !> Reference counter. Gives the user the possiblity to free after unref
    !! if present, it returns the number of counters associated to the object
    !! if absent f_unref raise an exception in the case the object is orphan
    integer, intent(out), optional :: count

    if (.not. associated(f_ref%iref)) then
       call f_err_throw('Illegal dereference: nullified object',&
            err_id=ERR_REFERENCE_COUNTERS)
    else
       if (f_ref%iref <= 1 .and. .not. present(count) .or. f_ref%iref==0 .and. present(count)) then
          call dump_ref_cnt(f_ref)
          call f_err_throw('Illegal dereference:'//&
               ' object is orphan, it should be have been freed',&
               err_id=ERR_REFERENCE_COUNTERS)
       else
          f_ref%iref=f_ref%iref-1
          if (present(count)) count=f_ref%iref
       end if
    end if

  end subroutine f_unref

  !> Returns the number of reference to an object.
  !! it returns a negative number if the object is nullified
  pure function f_ref_count(f_ref) result(count)
    implicit none
    type(f_reference_counter), intent(in) :: f_ref
    integer :: count

    if (associated(f_ref%iref)) then
       count=f_ref%iref
    else
       count=-1
    end if
  end function f_ref_count

  !> Free and check a reference counter
  !!this should be called when calling the destructor of the
  !!associated object
  subroutine f_ref_free(f_ref)
    implicit none
    type(f_reference_counter), intent(inout) :: f_ref

    if (.not. associated(f_ref%iref)) then
       call f_err_throw('Illegal free of reference counter, nullified object',&
            err_id=ERR_REFERENCE_COUNTERS)
    else
       if (f_ref%iref > 1) then
          call dump_ref_cnt(f_ref)
          call f_err_throw('Illegal free : object is still referenced',&
               err_id=ERR_REFERENCE_COUNTERS)
       else
          deallocate(f_ref%iref)
          call dict_free(f_ref%info)
          f_ref=f_ref_null()
       end if
    end if

  end subroutine f_ref_free

  !> Increase the reference counter of a associated source
  subroutine f_ref(src)
    implicit none
    type(f_reference_counter), intent(inout) :: src
    !check source validity
    if (.not. associated(src%iref)) then
       call f_err_throw('Illegal reference: nullified source',&
            err_id=ERR_REFERENCE_COUNTERS)
    else if (src%iref <= 0) then
       call dump_ref_cnt(src)
       call f_err_throw('Illegal reference: the source'//&
            ' is already dereferenced, it must be destroyed',&
            err_id=ERR_REFERENCE_COUNTERS)
    else
       src%iref=src%iref+1
    end if
  end subroutine f_ref

  !> Associate two reference objects.
  !! the destination is supposed to be in a nullified status,
  !! and the second one is supposed to be valid
  subroutine f_ref_associate(src,dest)
    use yaml_output, only: yaml_dict_dump
    implicit none
    !> Source reference. Should be in a valid state, which means
    !! that iref should be at least one.
    type(f_reference_counter), intent(in) :: src
    type(f_reference_counter), intent(inout) :: dest

    !check source validity
    if (.not. associated(src%iref)) then
       call f_err_throw('Illegal association: nullified source',&
            err_id=ERR_REFERENCE_COUNTERS)
    else if (src%iref <= 0) then
       call dump_ref_cnt(src)
       call f_err_throw('Illegal association: the source '//&
            'is already dereferenced, it must be destroyed',&
            err_id=ERR_REFERENCE_COUNTERS)
    end if

    !check destination suitablity
    if (associated(dest%iref)) then
       if (dest%iref == 1) then
          call dump_ref_cnt(dest)
          call f_err_throw('Illegal association: destination '//&
               ' is orphan, it must be destroyed to avoid memory leak',&
               err_id=ERR_REFERENCE_COUNTERS)
       end if
       !otherwise decrement it
       dest%iref=dest%iref-1
    end if

    !then reassociate the pointers
    dest%iref=>src%iref
    dest%info=>src%info

    dest%iref=dest%iref+1

  end subroutine f_ref_associate

  !Not used
!!$  subroutine f_ref_identify(dest,src)
!!$    implicit none
!!$    type(f_reference_counter), intent(out) :: dest
!!$    type(f_reference_counter), intent(in) :: src
!!$    
!!$    !the destination is completely scratched
!!$    call nullify_f_ref(dest)
!!$    call f_ref_associate(src,dest)
!!$  end subroutine f_ref_identify

end module f_refcnts
