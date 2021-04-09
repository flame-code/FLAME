!> @file
!!  Module defining a dictionary
!! @author Luigi Genovese
!!    Copyright (C) 2012-2015 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!> Module which defines a dictionary and the pure functions for its basic usage rules (no dependency)
module dictionaries_base
  use f_precisions
  implicit none

  integer, parameter, public :: max_field_length = 256    !< Maximum length of a field
  character(len=max_field_length), parameter :: TYPE_DICT='__dict__'
  character(len=max_field_length), parameter :: TYPE_LIST='__list__'

  character(len = max_field_length), parameter :: NOT_A_VALUE = "__not_a_value__"

  !> Global variables associated to the number of dictionaries allocated
  integer, private :: ndicts=0         !< Number of dictionaries allocated simultaneously
  integer, private :: ndicts_max=0     !< Maximum number of dictionaries allocated
  integer, private :: nfolders=0       !< Number of libraries
  integer, private :: nfolders_max=0   !< Maximum number of libraries allocated
  
  integer, parameter, private :: nfolder_size=10000 !< Size of the folder of pre-allocated dictionaries. Should be about 7 MB.
  
  
  type, public :: storage
     sequence
     integer :: item   !< Id of the item associated to the list
     integer :: nitems !< No. of items in the list
     integer :: nelems !< No. of items in the dictionary
     character(len=max_field_length) :: key
     character(len=max_field_length) :: value
  end type storage

  !> structure of the dictionary element (this internal structure is in principle private)
  type, public :: dictionary
     type(storage) :: data
     type(dictionary), pointer :: parent => null()
     type(dictionary), pointer :: next => null()
     type(dictionary), pointer :: child => null()
     type(dictionary), pointer :: previous => null()
  end type dictionary

!!$  type :: values
!!$     integer  :: inte 
!!$     double precision :: rel
!!$     integer(kind=8) :: value_type
!!$     character, dimension(:), pointer :: char
!!$     type(dictionary), pointer :: dict => null()
!!$     type(list), pointer :: list => null()
!!$  end type values
!!$
!!$  type :: link_elem
!!$     character(len=max_field_length) :: key
!!$     type(values), pointer :: elem
!!$     type(link_elem), pointer :: next => null()
!!$  end type link_elem
!!$
!!$  type, public :: dictionary
!!$     integer :: nelems !< No. of items in the dictionary
!!$     type(link_elem), dimension(:), pointer :: elements
!!$  end type dictionary
!!$
!!$  type, public :: list
!!$     integer :: nitems !< No. of items in the list
!!$     type(values), dimension(:), pointer :: elements
!!$  end type list
!!$
!!$  type, public :: list_iterator
!!$     integer :: item
!!$     type(list), pointer :: list
!!$  end type list_iterator
!!$
!!$  type, public :: dict_iterator
!!$     integer :: hash_index
!!$     type(dictionary), pointer :: dict
!!$     type(link_elem), pointer :: current
!!$  end type dict_iterator


  !> the database book is a workspace of pre-allocated dictionaries, which is used to manage dictionary creation
  type, private :: database_book
     !> this is the registry of the dictionary. It keeps track of the address of the associated dictionaries
     !It is a stack of free ifolder
     integer, dimension(:), pointer :: free_folders => null()
     integer :: nfree
     !>this is the place occupied by a dictionary, which is pre-allocated by
     !! the allocation in database_book. The routine dict_init associates the dictionary to a database item
     type(dictionary), dimension(:), pointer :: files => null()
     !> this is useful to create new books of the library in case the book is full
     type(database_book), pointer :: next => null()
     type(database_book), pointer :: previous => null()
  end type database_book


  !> Database of objects, associated contiguously and used for storage
  type(database_book), pointer, private :: library => null()

  !> Operator to access and create a key in the dictionary
  interface operator(//)
     module procedure get_child_ptr,get_list_ptr
  end interface

  interface dict_free
     module procedure dict_free, dict_free_multi
  end interface 

!!$  interface operator(.val.)
!!$     module procedure list_iter_to_value,dict_iter_to_value
!!$  end interface

  private :: allocate_library, allocate_file, deallocate_file, destroy_library
  public :: dict_free

contains
  
!!$  function list_iter_to_value(iter) return(val)
!!$    implicit none
!!$    type(list_iterator), intent(in) :: iter
!!$    type(value) :: val
!!$
!!$    val=val_null()
!!$    if (associated(iter%list)) then
!!$       val = get_value(iter%list,iter%item)
!!$    end if
!!$  end function list_iter_to_value
!!$
!!$  function get_value(list,item) return(val)
!!$    type(list) :: list
!!$    type(value) :: val
!!$    
!!$    val=val_null()
!!$    if (item < nitems) then
!!$       val=list%elements(item+1)
!!$    end if
!!$  end function get_value
!!$
!!$  function dict_iter_to_value(iter) return(val)
!!$    implicit none
!!$    type(dict_iterator), intent(in) :: iter
!!$    type(value) :: val
!!$
!!$    val=val_null()
!!$    if (associated(iter%current)) then
!!$       val = get_value(iter%current)
!!$    end if
!!$  end function list_iter_to_value
!!$
!!$  function get_value_le(link_elem) return(val)
!!$    type(link_elem) :: link_
!!$    type(value) :: val
!!$
!!$    val=link_elem%elem
!!$  end function get_value_le

  !> Test if keys are present
  pure function no_key(dict)
    implicit none
    type(dictionary), intent(in) :: dict
    logical :: no_key
    no_key=(len_trim(dict%data%key) == 0 .and. dict%data%item == -1) .and. &
         associated(dict%parent)
  end function no_key

  pure function no_value(dict)
    implicit none
    type(dictionary), intent(in) :: dict
    logical :: no_value

    no_value=trim(dict%data%value) == NOT_A_VALUE .and. .not. associated(dict%child)
  end function no_value


  pure function storage_null() result(st)
    type(storage) :: st
    st%key=repeat(' ',max_field_length)
    st%value(1:max_field_length)=NOT_A_VALUE
    st%item=-1
    st%nitems=0
    st%nelems=0
  end function storage_null


  pure subroutine dictionary_nullify(dict)
    implicit none
    type(dictionary), intent(inout) :: dict
    dict%data=storage_null()
    nullify(dict%child,dict%next,dict%parent,dict%previous)
  end subroutine dictionary_nullify


  pure subroutine allocate_library(library)
    implicit none
    type(database_book), intent(inout) :: library 
    !local variables
    integer :: ifolder

    allocate(library%free_folders(nfolder_size))
    allocate(library%files(nfolder_size))
    nullify(library%next)
    !all free places in the library
    library%nfree = nfolder_size
    do ifolder=1,nfolder_size
       library%free_folders(ifolder)=ifolder
    end do
  end subroutine allocate_library


  !> Assign a place in the library to the newcoming dictionary
  function allocate_file() result(dict)
    implicit none
    type(dictionary), pointer :: dict
    !local variables
    integer :: ifolder
    integer(kind=8), external :: f_loc
    type(database_book), pointer :: lib_tmp
        
    !first create the place for a folder if there is none
    if (.not. associated(library%free_folders)) then
       call allocate_library(library)
       nfolders=nfolders+1
       nfolders_max=max(nfolders_max,nfolders)
    end if

    ! We use the saved free place.
    ifolder = library%free_folders(nfolder_size - library%nfree + 1)
    call dictionary_nullify(library%files(ifolder))
    dict=>library%files(ifolder)
    !and immediately put it unavailable
    library%nfree = library%nfree - 1
    library%free_folders(nfolder_size - library%nfree)=0

    do while (library%nfree == 0)
       !then the library which is used is the following one
       if (.not. associated(library%next)) then
          allocate(library%next) !the pointer
          lib_tmp => library%next
          lib_tmp%previous => library
          call allocate_library(library%next)
          nfolders=nfolders+1
          nfolders_max=max(nfolders_max,nfolders)
       end if
       library => library%next
    end do
  end function allocate_file

  subroutine deallocate_file(file)
    implicit none
    integer(kind=8), intent(in) :: file !< the target is already inside

    type(database_book), pointer :: lib_tmp
    integer(kind=8), external :: f_loc

    !first, find the containing library
    lib_tmp => library
    do while (associated(lib_tmp%next))
       lib_tmp => lib_tmp%next
    end do
    do while (associated(lib_tmp))
       if (file >= f_loc(lib_tmp%files(1)) .and. &
            & file <= f_loc(lib_tmp%files(nfolder_size))) exit
       lib_tmp => lib_tmp%previous
    end do
    if (associated(lib_tmp)) then
       call deallocate_file_(file, lib_tmp)
    else
       !this might also mean that a variable has not been found
       stop 'dictionary has not been allocated by pool'
    end if

    ! Point on first library with free space.
    lib_tmp => library
    do while (associated(lib_tmp%next))
       lib_tmp => lib_tmp%next
    end do
    do while (associated(lib_tmp))
       if (lib_tmp%nfree > 0) library => lib_tmp
       lib_tmp => lib_tmp%previous
    end do

  contains

    !> Free the space from the library
    subroutine deallocate_file_(file, lib)
      implicit none
      integer(kind=8), intent(in) :: file !< the target is already inside
      type(database_book), pointer :: lib
      !local variables
      integer(kind=8) :: ifolder

      !find the file in the registry
      !then search in the present library the first free place

      ifolder = (file - f_loc(lib%files(1))) /&
           (f_loc(lib%files(2)) - f_loc(lib%files(1))) + 1

      !if the place has been found the dictionary can be identified
      !the registry is freed and the dictionary nullified
      lib%free_folders(nfolder_size - lib%nfree) = int(ifolder)
      call dictionary_nullify(lib%files(ifolder))
      lib%nfree = lib%nfree + 1
    end subroutine deallocate_file_

  end subroutine deallocate_file


  !> Terminate the library and free all memory space
  recursive subroutine destroy_library()
    implicit none
    !first, find last library
    if (associated(library%next)) then
       library=>library%next
       call destroy_library()
    else
       !once we are at the last tree free memory spaces
       deallocate(library%free_folders)
       nullify(library%free_folders)
       deallocate(library%files)
       nullify(library%files)
    end if
    !then come back of one step
    if (associated(library%previous)) then
       library => library%previous
       !and free the memory space of the library
       deallocate(library%next)
       nullify(library%next)
       nfolders=nfolders-1
       call destroy_library()
    end if
  end subroutine destroy_library


  subroutine dict_init(dict)
    !Initialize a dictionary, ready for usage.
    !The dictionary may be undefined on input, thus care must be taken in freeing previously
    !initialized instances before calling :f:subr:`dict_init` again.
    !:Example:
    ! .. code-block:: fortran
    !
    !       program dictinit
    !         use dictionaries
    !         implicit none
    !         type(dictionary), pointer :: dict
    !         call dict_init(dict)
    !         ![...] perform actions with the dict
    !         call dict_free(dict)
    !         !and then dict (which is now nullified) might be reinitialized again
    !       end program dictinit
    implicit none
    type(dictionary), pointer :: dict !the dictionary to be iniitalized
    !local variables
    !integer :: nres

    if (nfolder_size == 0) then
       !this is normal allocation, let the system do the allocation
       allocate(dict)
       call dictionary_nullify(dict)
    else
       !this strategy retrieves the allocation from the library
       if (.not. associated(library)) allocate(library) ! the pointer
       dict=>allocate_file()
    end if
    !global variables counting the number of dictionaries simultaneously used
    ndicts=ndicts+1
    ndicts_max=max(ndicts_max,ndicts)
    !!debug, to check that the allocation has been done properly
    !call find_residual_dicts(nres)
    !if (nres /= ndicts) then
    !   print *,'internal error',nres,ndicts,nfolder_size-library%nfree
    !   call dictionary_check_leak()
    !   stop
    !end if
  end subroutine dict_init


  !> Destroy only one level
  subroutine dict_destroy(dict)
    implicit none
    type(dictionary), intent(in) :: dict !to copy the target
    !local variables
    !integer :: nres
    integer(kind=8), external :: f_loc
    !free a space in the library and let the allocation live
    call deallocate_file(f_loc(dict))
    ndicts=ndicts-1
    !!debug, to check that the allocation has been done properly
    !call find_residual_dicts(nres)
    !if (nres /= ndicts) then
      !print *,'internal error in free',nres,ndicts,nfolder_size-library%nfree
    !print *,'deall',trim(dict%data%key),'key',trim(dict%data%value),'val'
    !   print *,'floc',f_loc(dict)
    !   call dictionary_check_leak()
    !   stop
    !end if
    !if there are no dictionaries anymore the library can be destroyed
    if (ndicts==0 .and. associated(library)) then
       call destroy_library()
       !otherwise free the last tree 
       deallocate(library)
       nullify(library)
    end if
  end subroutine dict_destroy

  !> free different dictionaries at the same time
  subroutine dict_free_multi(dict0,dict1,dict2,dict3,dict4,dict5,dict6,dict7,dict8,dict9)
    type(dictionary), pointer :: dict0,dict1
    type(dictionary), pointer, optional :: dict2,dict3,dict4,dict5,dict6,dict7,dict8,dict9
    
    call dict_free(dict0)
    call dict_free(dict1)
    if (present(dict2)) call dict_free(dict2)
    if (present(dict3)) call dict_free(dict3)
    if (present(dict4)) call dict_free(dict4)
    if (present(dict5)) call dict_free(dict5)
    if (present(dict6)) call dict_free(dict6)
    if (present(dict7)) call dict_free(dict7)
    if (present(dict8)) call dict_free(dict8)
    if (present(dict9)) call dict_free(dict9)

  end subroutine dict_free_multi

  subroutine dict_free(dict)
    type(dictionary), pointer :: dict
    if (associated(dict)) then
       call dict_free_(dict)
       nullify(dict)
    end if

  contains

    !> Implement the recursivity only on children 
    recursive subroutine dict_free_(dict)
      implicit none
      type(dictionary), pointer :: dict
      !local variables
      type(dictionary), pointer :: dict_iter,child,current
      dict_iter=>dict
      do while(associated(dict_iter))
         child=>dict_iter%child
         current=>dict_iter
         dict_iter=>dict_iter%next
         !destroy current
         if (nfolder_size==0) then
            deallocate(current)
         else
            call dict_destroy(current)
         end if
         nullify(current)
         !destroy the children
         if (associated(child)) then
            call dict_free_(child)
         end if
      end do
    end subroutine dict_free_

  end subroutine dict_free

  pure &!do not remove this continuation line, needed for autodoc
       function dict_len(dict)
    !Provides the length of the dictionary in the case it is used as a *list* of objects.
    !Returns 0 if the dictionary is a mapping or a scalar. The return value is -1 is the 
    !dictionary is nullified. For mapping the :f:func:`dict_size` function will have to be used.
    !
    !:Attributes: pure
    !
    !:Example:
    ! .. code-block:: fortran
    !
    !       program dictlen
    !         use dictionaries
    !         implicit none
    !         integer :: ilen
    !         type(dictionary), pointer :: dict
    !         call dict_init(dict)
    !         call dict_add(dict,'item1')
    !         call dict_add(dict,'item2')
    !         call dict_add(dict,'item3')
    !         ilen = dict_len(dict) !ilen has value 3
    !         call dict_free(dict)
    !         ilen = dict_len(dict) !ilen has value -1
    !       end program dictlen
    implicit none
    type(dictionary), intent(in), pointer :: dict 
    integer :: dict_len !Length of the dict as list container
    
    if (associated(dict)) then
       dict_len=dict%data%nitems
    else
       dict_len=-1
    end if
  end function dict_len


  !> Return the size of the dictionary
  pure &!do not remove this continuation line, needed for autodoc
       function dict_size(dict)
    !Provides the length of the dictionary in the case it is used as a *mapping* of objects.
    !Returns 0 if the dictionary is a list or a scalar. The return value is -1 is the 
    !dictionary is nullified. For lists the :f:func:`dict_len` function will have to be used.
    !
    !:Attributes: pure
    !
    !:Example:
  ! .. code-block:: fortran
    !
    !       program dictlen
    !         use dictionaries
    !         implicit none
    !         integer :: isize
    !         type(dictionary), pointer :: dict
    !         call dict_init(dict)
    !         call dict_set(dict//'key1','value1')
    !         call dict_set(dict//'key2','value2')
    !         call dict_set(dict//'key3','value3')
    !         isize = dict_size(dict) !ilen has value 3
    !         call dict_free(dict)
    !         isize = dict_size(dict) !ilen has value -1
    !       end program dictlen
    implicit none
    type(dictionary), intent(in), pointer :: dict
    integer :: dict_size  !Length of the dict as mapping container
    
    if (associated(dict)) then
       dict_size=dict%data%nelems
    else
       dict_size=-1
    end if
  end function dict_size


  !> This function returns the key if present otherwise the value of the element if in a list
  !! this is a function indended for internal flib usage which 
  !! can be used for lists of hash table, as dictionary type is "polymorph"
  !pure 
  function name_is(dict,name)
    implicit none
    type(dictionary), pointer, intent(in) :: dict
    character(len=*), intent(in) :: name
    logical :: name_is

!!$    print *,'name is',trim(name),no_key(dict%child),no_value(dict%child),&
!!$         'value',trim(dict%data%value),'key',trim(dict%data%key),&
!!$         'item',dict%data%item,&
!!$         'child',associated(dict%child),'valueC',trim(dict%child%data%value),&
!!$         'keyC',trim(dict%child%data%key)

    name_is=.false.
    if (trim(name) == trim(dict%data%key)) then
       name_is=.true.
    else if (associated(dict%child)) then
       !most likely dict is a item list
       name_is=(trim(name) == trim(dict%child%data%key))
       !print *,'here',name_is,trim(dict%child%data%value),'ag',trim(name)
    else
       name_is=(trim(name) == trim(dict%data%value))
    end if
  end function name_is

  !> Returns the value of the key of the dictionary
  pure &!do not remove, for autodoc
       function dict_key(dict)
  !Returns the value of the key of the dictionary. Useful especially when using 
  !the dictionary as an iterator among a mapping
  !
  !:Attributes: pure
  !
  !:Example:
  !
  !   .. literalinclude:: ../tests/flib/dicts.f90
  !      :language: fortran
  !      :start-after: perform an iterator on dictA
  !      :end-before: end of perform an iterator on dictA
  !
  implicit none
    type(dictionary), pointer, intent(in) :: dict 
    character(len=max_field_length) :: dict_key !key of `dict`, empty string if nullified
    
    if (associated(dict)) then 
       !call check_key(dict)
       dict_key=dict%data%key
    else
       dict_key=repeat(' ',len(dict_key))
    end if
  end function dict_key


  !> Returns the value of the key of the dictionary
  pure &!for autodoc
       function dict_item(dict)
  !Returns the value of the item of the dictionary. Useful especially when using 
  !the dictionary as an iterator among a mapping
  !
  !:Attributes: pure
  !
  !:Example:
  !
  !   .. literalinclude:: ../tests/flib/dicts.f90
  !      :language: fortran
  !      :start-after: perform an iterator with item on dict
  !      :end-before: end of perform an iterator with item
  !
  implicit none
    type(dictionary), pointer, intent(in) :: dict
    integer :: dict_item !item of `dict`, -1 if nullified

    if (associated(dict)) then 
       dict_item=dict%data%item
    else
       dict_item=-1
    end if
  end function dict_item

  
  !> Value of the dictionary, if present, otherwise empty
  !! if the value is a dictionary, it returns __dict__ in the character
  pure &!
       function dict_value(dict)
  !Returns the value dictionary. Useful especially when using 
  !the dictionary as an iterator among a mapping. Useful for scalar values.
  !If the value is not a scalar, it returns either :f:var`TYPE_DICT` or :f:var:`TYPE_LIST`.
  !
  !:Attributes: pure
  !
  !:Example:
  !
  !   .. literalinclude:: ../tests/flib/dicts.f90
  !      :language: fortran
  !      :start-after: perform an iterator with item on dict
  !      :end-before: end of perform an iterator with item
  !
    implicit none
    type(dictionary), pointer, intent(in) :: dict
    character(len=max_field_length) :: dict_value !value of the dictionary, meaningful if scalar.

    if (associated(dict)) then 
       !call check_key(dict)
       if (associated(dict%child)) then
          if (dict%data%nitems > 0) then
             dict_value=TYPE_LIST
          else if (dict%data%nelems > 0) then
             dict_value=TYPE_DICT
          else
             !illegal condition
             dict_value= NOT_A_VALUE
          end if
       else if (trim(dict%data%value) == NOT_A_VALUE) then
          dict_value=repeat(' ',len(dict_value))
       else
          dict_value=dict%data%value
       end if
    else
       dict_value=repeat(' ',len(dict_value))
    end if

  end function dict_value


  !non-pure subroutines, due to pointer assignments


  !> Define the same parent(dict) for any of the elements of the linked chain (child)
  recursive subroutine define_parent(dict,child)
    implicit none
    !Arguments
    type(dictionary), target :: dict
    type(dictionary) :: child
!!$    type(dictionary), pointer :: dict
!!$    type(dictionary), pointer :: child

    !local variables
!!$    type(dictionary), pointer :: iter
!!$    !eliminate recursion
!!$    iter => child
!!$    do while(associated(iter))
!!$       iter%parent=>dict
!!$       iter => iter%next
!!$    end do

    child%parent=>dict
    if (associated(child%next)) call define_parent(dict,child%next)
  end subroutine define_parent


  !> Set brother as the previous element of dict
  subroutine define_brother(brother,dict)
    implicit none
    type(dictionary), target :: brother
    type(dictionary) :: dict

    dict%previous=>brother
  end subroutine define_brother


  !> This routine creates a key for the dictionary in case it is absent
  !! the it adds one to the number of elements of the parent dictionary
  pure subroutine set_elem(dict,key)
    use yaml_strings, only: f_strcpy
    implicit none
    type(dictionary), pointer :: dict !!TO BE VERIFIED
    character(len=*), intent(in) :: key

    !print *,'set_elem in ',trim(key),dict%data%nelems,dict%parent%data%nelems
    call f_strcpy(src=trim(key),dest=dict%data%key)
    if (associated(dict%parent)) then
       dict%parent%data%nelems=dict%parent%data%nelems+1
    else
       dict%data%nelems=dict%data%nelems+1
    end if
    !print *,'set_elem out ',trim(key),dict%data%nelems,dict%parent%data%nelems

  end subroutine set_elem


  !> Associates an extra item to the dictionary and takes care of the 
  !! increment of the number of items
  !! this routine adds the check that the number of items is preserved
  !! for the moment this check is removed
  pure subroutine set_item(dict,item)
    implicit none
    type(dictionary), pointer :: dict !< recently modified, to be checked
    integer, intent(in) :: item

    dict%data%item=item
    if (associated(dict%parent)) then
       if (len_trim(dict%parent%data%value) > 0) dict%parent%data%value=repeat(' ',max_field_length)
       dict%parent%data%nitems=dict%parent%data%nitems+1
    end if

  end subroutine set_item


  !> Retrieve the pointer to the dictionary which has this key.
  !! If the key does not exists, create it in the child chain
  function get_child_ptr(dict,key) result(subd_ptr)
    implicit none
    type(dictionary), intent(in), pointer :: dict 
    character(len=*), intent(in) :: key
    type(dictionary), pointer :: subd_ptr

    if (.not. associated(dict)) then
       nullify(subd_ptr)
       return
    end if
    !!commented out, the key is checked only when retrieving
    !call check_key(dict)
    if (associated(dict%child)) then
       !subd_ptr => get_dict_ptr(dict%child,key)
       subd_ptr => get_dict_from_key(dict%child,key,create=.true.)
    else
       call dict_init(dict%child)
       call define_parent(dict,dict%child)
       call set_elem(dict%child,key)
       subd_ptr => dict%child
    end if
  end function get_child_ptr


  !>points to the dictionary which has the key.
  function get_dict_from_key(dict,key,create) result (dict_ptr)
    implicit none
    !> root of the dictionary to start the search from
    type(dictionary), intent(in), pointer :: dict 
    !> key that has to be matched, trailing blanks excluded
    character(len=*), intent(in) :: key
    !> default .false. if present with value .true., a item is created in the dictionary with the key
    logical, intent(in), optional :: create 
    type(dictionary), pointer :: dict_ptr
    !local variables
    logical :: crt
    type(dictionary), pointer :: iter

    crt=.false.
    if (present(create)) crt=create
    !iterate until key found
    nullify(dict_ptr)
    iter => dict
    seek: do 
       if (iter%data%key == trim(key)) then
          dict_ptr=> iter
          exit seek
       else if (associated(iter%next)) then
          iter => iter%next
          cycle seek
       else 
          exit seek
       end if
    end do seek

    !this is useful for the first assignation, might be moved at the beginning
    if (crt) then
       if (no_key(iter)) then 
          call set_elem(iter,key)
          dict_ptr => iter
       end if
       !if we did not find the key, decide to create it
       if (.not. associated(dict_ptr)) then
          call dict_init(iter%next)
          call define_brother(iter,iter%next) !chain the list in both directions
          if (associated(iter%parent)) call define_parent(iter%parent,iter%next)
          call set_elem(iter%next,key)
          dict_ptr => iter%next
       end if
!!$       !might be refactored into
!!$       if (.not. associated(dict_ptr)) then
!!$          call init_next(iter)
!!$          call set_elem(iter,key)
!!$          dict_ptr=>iter
!!$       end if
    end if
  end function get_dict_from_key
!!$
!!$  !> Retrieve the pointer to the dictionary which has this key.
!!$  !! If the key does not exists, create it in the next chain 
!!$  !! Key Must be already present 
!!$  recursive function get_dict_ptr(dict,key) result (dict_ptr)
!!$    implicit none
!!$    type(dictionary), intent(in), pointer :: dict !hidden inout
!!$    character(len=*), intent(in) :: key
!!$    type(dictionary), pointer :: dict_ptr
!!$
!!$!    print *,'here',trim(key)
!!$    !follow the chain, stop at the first occurence
!!$    if (trim(dict%data%key) == trim(key)) then
!!$       dict_ptr => dict
!!$    else if (associated(dict%next)) then
!!$       dict_ptr => get_dict_ptr(dict%next,key)
!!$    else if (no_key(dict)) then !this is useful for the first assignation
!!$       call set_elem(dict,key)
!!$       dict_ptr => dict
!!$    else
!!$       call dict_init(dict%next)
!!$       call define_brother(dict,dict%next) !chain the list in both directions
!!$       if (associated(dict%parent)) call define_parent(dict%parent,dict%next)
!!$       call set_elem(dict%next,key)
!!$       dict_ptr => dict%next
!!$    end if
!!$
!!$  end function get_dict_ptr


!!$  !> Retrieve the pointer to the item of the list.
!!$  !! If the list does not exists, create it in the child chain.
!!$  !! If the list is too short, create it in the next chain
!!$  recursive function get_item_ptr(dict,item) result (item_ptr)
!!$    implicit none
!!$    type(dictionary), intent(in), pointer :: dict !hidden inout
!!$    integer, intent(in) :: item
!!$    type(dictionary), pointer :: item_ptr
!!$
!!$    !follow the chain, stop at  first occurence
!!$    if (dict%data%item == item) then
!!$       item_ptr => dict
!!$    else if (associated(dict%next)) then
!!$       item_ptr => get_item_ptr(dict%next,item)
!!$    else if (no_key(dict)) then
!!$       call set_item(dict,item)
!!$       item_ptr => dict
!!$    else
!!$       call dict_init(dict%next)
!!$       call define_brother(dict,dict%next) !chain the list in both directions
!!$       if (associated(dict%parent)) call define_parent(dict%parent,dict%next)
!!$       call set_item(dict%next,item)
!!$       item_ptr => dict%next
!!$    end if
!!$  end function get_item_ptr


  !> Retrieve the pointer to the item of the list.
  !! If the list does not exists, create it in the child chain.
  !! If the list is too short, create it in the next chain
  subroutine item_ptr_find(dict,item,item_ptr)
    implicit none
    type(dictionary), intent(in), pointer :: dict 
    integer, intent(in) :: item
    type(dictionary), pointer :: item_ptr

    item_ptr=>dict
    find_item: do 
       if (item_ptr%data%item == item) exit find_item
       if (associated(item_ptr%next)) then
          item_ptr=>item_ptr%next
          cycle find_item
       end if
       if (.not. no_key(item_ptr)) call init_next(item_ptr)
       call set_item(item_ptr,item)
       exit find_item
    end do find_item
  end subroutine item_ptr_find

  subroutine init_next(dict)
    implicit none
    type(dictionary), pointer :: dict

    call dict_init(dict%next)
    call define_brother(dict,dict%next) !chain the list in both directions
    if (associated(dict%parent)) &
         call define_parent(dict%parent,dict%next)
    dict=>dict%next          
  end subroutine init_next


  !> Retrieve the pointer to the item of the list.
  !! If the list does not exists, create it in the child chain.
  !! If the list is too short, create it in the next chain
  function get_list_ptr(dict,item) result (subd_ptr)
    implicit none
    type(dictionary), intent(in), pointer :: dict !hidden inout
    integer, intent(in) :: item
    type(dictionary), pointer :: subd_ptr

    !commented out for the moment
    !call check_key(dict)
    if (.not. associated(dict)) then
       nullify(subd_ptr)
       return
    end if
    
    !print *,'nelems',dict%data%nelems,dict%data%nitems
    !free previously existing dictionaries
    call clean_subdict(dict)
    
    if (associated(dict%child)) then         
       !subd_ptr => get_item_ptr(dict%child,item)
       call item_ptr_find(dict%child,item,subd_ptr)
    else
       call dict_init(dict%child)
       call define_parent(dict,dict%child)
       call set_item(dict%child,item)
       subd_ptr => dict%child
    end if

  end function get_list_ptr

  !> Defines a storage structure with a key-value couple
  elemental pure function storage_data(key,val)
    use yaml_strings, only: f_strcpy
    character(len=*), intent(in) :: key,val
    type(storage) :: storage_data

    storage_data=storage_null()

    call f_strcpy(src=key,dest=storage_data%key)
    call f_strcpy(src=val,dest=storage_data%value)

  end function storage_data

  !> Test to see the g95 behaviour
  pure function stored_key(st) result(key)
    use yaml_strings, only: f_strcpy
    implicit none
    type(storage), intent(in) :: st
    character(len=max_field_length) :: key
    call f_strcpy(src=st%key,dest=key)
  end function stored_key

  !> test to see the g95 behaviour
  pure function stored_value(st) result(val)
    use yaml_strings, only: f_strcpy
    implicit none
    type(storage), intent(in) :: st
    character(len=max_field_length) :: val
    call f_strcpy(src=st%value,dest=val)
  end function stored_value


  !> Clean the child and put to zero the number of elements in the case of a dictionary
  !! pure 
  subroutine clean_subdict(dict)
    implicit none
    type(dictionary), pointer :: dict
    
    if (associated(dict)) then
       if (dict%data%nelems > 0) then
          call dict_free(dict%child)
          dict%data%nelems=0
       end if
    end if
    
  end subroutine clean_subdict


  !> Export the number of dictionaries
  !! this routine is useful to understand the total usage of the 
  !! dictionaries in the f_lib module
  subroutine dict_get_num(ndict,ndict_max,nlibs,nlibs_max)
    implicit none
    integer, intent(out) :: ndict     !< actual number of dictionaries active
    integer, intent(out) :: ndict_max !< maximum number of dictionaries active during the session
    integer, intent(out) :: nlibs     !< actual number of libraries active at present
    integer, intent(out) :: nlibs_max !< maximum number of libraries active at present (should be nlibs)

    ndict=ndicts
    ndict_max=ndicts_max
    nlibs=nfolders
    nlibs_max=nfolders_max
  end subroutine dict_get_num

  !> Check for leaked dictionaries and prints out in stdout 
  !! the dictionaries which have not been freed, without mutual structure
  subroutine dictionary_check_leak()
    implicit none
    !local variables
    integer :: ndict,ndict_max,nlibs,nlibs_max

    call dict_get_num(ndict,ndict_max,nlibs,nlibs_max)
    if (ndict /= 0 ) then
       write(*,'(a,i0)')&
            '#Error, found unfreed dictionaries after finalization: ',&
            ndict
       call find_residual_dicts()
    end if

  end subroutine dictionary_check_leak


  !> Debug subroutine
  subroutine find_residual_dicts(nres)
    implicit none
    !>if present, counts the number of dictionaries
    !! and not dump the output
    integer, intent(out), optional :: nres
    !local variables
    integer :: i,j
    type(database_book), pointer :: lib
    logical, dimension(:), allocatable :: freed

    nullify(lib)
    if (associated(library)) then
       lib=>library
       do while(associated(lib%next))
          lib=>lib%next
       end do
    end if
    !then start by printing the residual keys
    allocate(freed(nfolder_size))
    freed(:) = .false.
    if (present(nres)) nres=0
    do while(associated(lib))
       !print *,'nowcheck',lib%nfree
       do i=1,lib%nfree
          j=lib%free_folders(nfolder_size - i + 1)
          if (j/=0) then
             if (.not. freed(j)) then
                freed(j) = .true.
             else
                print *,'internal, same place freed twice',j
                !stop
             end if
          end if
       end do
       do i = 1, nfolder_size
          if (.not. freed(i)) then
             if (present(nres)) then
                nres=nres+1
             else
                if (associated(lib%files(i)%child)) then
                   write(*,'(a,a,a,a,a)') &
                        & '#Unfreed key: "',trim(lib%files(i)%data%key),&
                        & '" value: "',trim(TYPE_DICT),'"'
                else
                   write(*,'(a,a,a,a,a)') &
                        & '#Unfreed key: "',trim(lib%files(i)%data%key),&
                        & '" value: "',trim(lib%files(i)%data%value),'"'
                end if
             end if
          end if
       end do
       lib=>lib%previous
    end do
    deallocate(freed)

  end subroutine find_residual_dicts
end module dictionaries_base
