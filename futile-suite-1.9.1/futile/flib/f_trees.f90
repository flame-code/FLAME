  !> @file
  !! Define the trees objects and operators on them
  !! It inherits (manually...) the methods from the dictionaries
  !! like operations on external files and basic operations in memory
  !! @author
  !!    Copyright (C) 2012-2015 BigDFT group
  !!    This file is distributed under the terms of the
  !!    GNU General Public License, see ~/COPYING file
  !!    or http://www.gnu.org/copyleft/gpl.txt .
  !!    For the list of contributors, see ~/AUTHORS
module f_trees
  use dictionaries
  use f_refcnts
  implicit none

  private

  !> structure version for the dictionary
  !! does not require the pointer attribute and can therefore be more easily binded
  !! a reference counter is added to preserve the association status
  !! in a modern (sic!) version of fortran this structure might have its methods as type
  !! procedures
  !! @warning: this structure is for small usage only. It is supposed to be (even) slower 
  !! that the usual implementation of dictionaries. So use with care
  type, public :: f_tree
     type(f_reference_counter) :: rc !<reference counter of the tree
     type(dictionary), pointer :: d=>null()  !< the actual dictionary implementation
  end type f_tree

  interface operator(//)
     module procedure f_tree_getitem,f_tree_getkey
  end interface operator(//)
  interface assignment(=)
     module procedure f_tree_scpy
  end interface assignment(=)
     

  interface f_tree_new
     module procedure f_tree_new_empty!,f_tree_dict_new,f_tree_list_new,f_tree_dict_new
  end interface f_tree_new
  interface f_tree_push
     module procedure f_tree_push_c0
  end interface f_tree_push


  public :: assignment(=),operator(//)
  public :: f_tree_new,f_tree_cpy,f_tree_push,f_tree_free,f_tree_dump,f_tree_load

  contains

    !new, !push !copy !equal (shallowcopy) !free !dump !getitem
    !iterator

    function f_tree_new_empty() result(ft)
      implicit none
      type(f_tree) :: ft
      call dict_init(ft%d)
      call nullify_f_ref(ft%rc)
    end function f_tree_new_empty

    function f_tree_getitem(src,item) result(ft)
      implicit none
      type(f_tree), intent(in) :: src
      integer, intent(in) :: item
      type(f_tree) :: ft
      
      ft%d=>src%d//item
      !nullify the reference counter, to be considered when doing a copy
      call nullify_f_ref(ft%rc)
      
    end function f_tree_getitem
    function f_tree_getkey(src,key) result(ft)
      implicit none
      type(f_tree), intent(in) :: src
      character(len=*), intent(in) :: key
      type(f_tree) :: ft

      ft%d=>src%d//key
      !nullify the reference counter, to be considered when doing a copy
      call nullify_f_ref(ft%rc)

    end function f_tree_getkey

    !> initialize the character variable
    subroutine f_tree_push_c0(ft,val)
      implicit none
      type(f_tree), intent(in) :: ft
      character(len=*), intent(in) :: val
      call set(ft%d,val)
    end subroutine f_tree_push_c0

    function f_tree_load(stream) result(ft)
      use yaml_parse
      implicit none
      character(len=*), intent(in) :: stream
      type(f_tree) :: ft
      
      !ft%rc=f_ref_new('loaded')
      ft%d => yaml_load(stream)
    end function f_tree_load

    subroutine f_tree_free(ft)
      implicit none
      type(f_tree), intent(inout) :: ft
      !local variables
      integer :: count

      !verify if someone has already referenced the same object
      call f_unref(ft%rc,count=count)
      if (count==0) then
         call dict_free(ft%d)
         call f_ref_free(ft%rc)
      else
         nullify(ft%d)
      end if
    end subroutine f_tree_free

    subroutine f_tree_scpy(dest,src)
      implicit none
      type(f_tree), intent(inout) :: dest
      type(f_tree), intent(in) :: src

      !erase the pre-exisiting dictionary if needed
      if (associated(dest%d)) call f_tree_free(dest)
      !shallow copy
      !if the dictionary is stack-created, then create the reference
      if (f_ref_count(src%rc) == -1) then
         dest%rc=f_ref_new('newtree')
      else
         call nullify_f_ref(dest%rc)
         call f_ref_associate(src=src%rc,dest=dest%rc)
      end if
      dest%d=>src%d
      
    end subroutine f_tree_scpy

    !which one is more elegant in fortran?
    !dict1 = .copy. dict2
    !call f_tree_cpy(src=dict1,dest=dict2)

    subroutine f_tree_cpy(dest,src)
      implicit none
      type(f_tree), intent(inout) :: dest
      type(f_tree), intent(in) :: src
     
      dest%rc=f_ref_new('tree')
      nullify(dest%d)
      call dict_copy(src=src%d,dest=dest%d)

    end subroutine f_tree_cpy

    subroutine f_tree_dump(ft)
      use yaml_output
      implicit none
      type(f_tree), intent(in) :: ft
      
      if (f_ref_count(ft%rc) == -1) then
         call f_err_throw('The tree cannot be dumped without affecting it,'//&
              ' this operation is illegal as it will create a memory leak')
      end if

      if (associated(ft%d)) then
         call yaml_dict_dump(ft%d)
      else
         call yaml_mapping_open(advance='no')
         call yaml_mapping_close()
      end if
    end subroutine f_tree_dump

end module f_trees
