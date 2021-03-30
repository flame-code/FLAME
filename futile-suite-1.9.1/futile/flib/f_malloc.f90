!> @file
!! Manage dynamic memory allocation control structures
!! @author
!!    Copyright (C) 2012-2015 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Module used to manage memory allocations control structures.
!! This module has to be intended as a submodule of dynamic_memory module
module module_f_malloc

  use f_precisions
  use dictionaries, only: f_err_throw,f_err_raise,dictionary
  use yaml_strings, only: f_strcpy,operator(+)

  implicit none

  private

  integer, parameter, public :: f_kind=f_long
                                      !! this parameter can be modified only by dynamic memory module
  integer(f_integer), parameter, public :: f_malloc_namelen=32          !< length of the character variables
  integer, parameter :: f_malloc_info_maxlength=128 !< for the extra information
  integer, parameter :: max_rank=7          !< maximum rank in fortran

  !to be initialized in the dynamic_memory module
  integer, save, public :: ERR_INVALID_MALLOC

  logical, save, public :: f_malloc_default_profiling=.true.
  character(len=f_malloc_namelen), save, public :: f_malloc_routine_name=repeat(' ',f_malloc_namelen)

  character(len=*), parameter, public :: INFO_TYPE_KEY='Type'
  character(len=*), parameter, public :: INFO_ALIGNMENT_KEY='alignment'
  character(len=*), parameter, public :: INFO_SHARED_TYPE='SHARED'


  type, public:: f_workspace
     integer(f_long) :: pos_r,pos_d,pos_i,pos_li,pos_l
     integer(f_integer), dimension(:), pointer :: ptr_i
     integer(f_long), dimension(:), pointer :: ptr_li
     real(f_simple), dimension(:), pointer :: ptr_r
     real(f_double), dimension(:), pointer :: ptr_d
     logical, dimension(:), pointer :: ptr_l
  end type f_workspace


  !> Structure needed to allocate an allocatable array
  type, public :: malloc_information_all
     logical :: pin                          !< flag to control the pinning of the address
     logical :: profile                      !< activate profiling for this allocation
     logical :: put_to_zero                  !< initialize to zero after allocation
     integer :: rank                         !< rank of the array
     integer(f_kind), dimension(max_rank) :: shape   !< shape of the structure
     integer(f_kind), dimension(max_rank) :: lbounds !< lower bounds
     integer(f_kind), dimension(max_rank) :: ubounds !< upper bounds
     integer(f_address) :: srcdata_add          !< physical address of source data
     character(len=f_malloc_namelen) :: array_id      !< label the array
     character(len=f_malloc_namelen) :: routine_id    !< label the routine
     character(len=f_malloc_info_maxlength) :: info !<information to the database for handling various approaches for allocation
  end type malloc_information_all

  !> Structure needed to allocate an allocatable array of string of implicit length (for non-2003 compilers)
  type, public :: malloc_information_str_all
     logical :: pin                          !< flag to control the pinning of the address
     logical :: profile                      !< activate profiling for this allocation
     logical :: put_to_zero                  !< initialize to zero after allocation
     integer :: rank                         !< rank of the array
     integer :: len                  !< length of the character
     integer(f_kind), dimension(max_rank) :: shape   !< shape of the structure
     integer(f_kind), dimension(max_rank) :: lbounds !< lower bounds
     integer(f_kind), dimension(max_rank) :: ubounds !< upper bounds
     integer(f_address) :: srcdata_add          !< physical address of source data
     character(len=f_malloc_namelen) :: array_id      !< label the array
     character(len=f_malloc_namelen) :: routine_id    !< label the routine
     character(len=f_malloc_info_maxlength) :: info !<information to the database for handling various approaches for allocation
  end type malloc_information_str_all

  !> Structure needed to allocate a pointer
  type, public :: malloc_information_ptr
     logical :: ptr                          !< just to make the structures different, to see if needed
     logical :: pin                          !< flag to control the pinning of the address
     logical :: profile                      !< activate profiling for this allocation
     logical :: put_to_zero                  !< initialize to zero after allocation
     integer :: rank                         !< rank of the pointer
     integer(f_kind), dimension(max_rank) :: shape   !< shape of the structure
     integer(f_kind), dimension(max_rank) :: lbounds !< lower bounds
     integer(f_kind), dimension(max_rank) :: ubounds !< upper bounds
     integer(f_address) :: srcdata_add          !< physical address of source data
     character(len=f_malloc_namelen) :: array_id      !< label the array
     character(len=f_malloc_namelen) :: routine_id    !< label the routine
     type(f_workspace), pointer :: w !<workspace structure if the array is a work array
     character(len=f_malloc_info_maxlength) :: info !<information to the database for handling various approaches for allocation
  end type malloc_information_ptr

  !> Structure needed to allocate a f_buffer object
  type, public :: malloc_information_buf
     logical :: pin                          !< flag to control the pinning of the address
     logical :: profile                      !< activate profiling for this allocation
     logical :: put_to_zero                  !< initialize to zero after allocation
     integer :: rank                         !< rank of the pointer
     integer(f_kind), dimension(max_rank) :: shape   !< shape of the structure
     integer(f_kind), dimension(max_rank) :: lbounds !< lower bounds
     integer(f_kind), dimension(max_rank) :: ubounds !< upper bounds
     integer(f_address) :: srcdata_add          !< physical address of source data
     character(len=f_malloc_namelen) :: array_id      !< label the array
     character(len=f_malloc_namelen) :: routine_id    !< label the routine
     character(len=f_malloc_info_maxlength) :: info !<information to the database for handling various approaches for allocation
     !all the components which are added from this line have to be initialized separately
  end type malloc_information_buf


  !> Structure needed to allocate a pointer of string of implicit length (for non-2003 complilers)
  type, public :: malloc_information_str_ptr
     logical :: ptr                          !< just to make the structures different, to see if needed
     logical :: pin                          !< flag to control the pinning of the address
     logical :: profile                      !< activate profiling for this allocation
     logical :: put_to_zero                  !< initialize to zero after allocation
     integer :: rank                         !< rank of the pointer
     integer :: len                  !< length of the character
     integer(f_kind), dimension(max_rank) :: shape   !< shape of the structure
     integer(f_kind), dimension(max_rank) :: lbounds !< lower bounds
     integer(f_kind), dimension(max_rank) :: ubounds !< upper bounds
     integer(f_address) :: srcdata_add          !< physical address of source data
     character(len=f_malloc_namelen) :: array_id      !< label the array
     character(len=f_malloc_namelen) :: routine_id    !< label the routine
     character(len=f_malloc_info_maxlength) :: info !<information to the database for handling various approaches for allocation
  end type malloc_information_str_ptr

  type, public :: array_bounds
     integer(f_kind) :: nlow  !<lower bounds
     integer(f_kind) :: nhigh !<higher bounds
  end type array_bounds

  interface operator(.to.)
     module procedure f_array_bounds_intint, f_array_bounds_longlong
     module procedure f_array_bounds_intlong, f_array_bounds_longint
  end interface

  interface operator(.plus.)
     module procedure f_array_concatenate_ii,f_array_concatenate_ii0
  end interface

  interface nullify_malloc_information
     module procedure nullify_malloc_information_all
     module procedure nullify_malloc_information_ptr
     module procedure nullify_malloc_information_buf
     module procedure nullify_malloc_information_str_all
     module procedure nullify_malloc_information_str_ptr
  end interface

  !> Fake structure needed to document common arguments of the module
  type, private :: doc
     !> integer indicating the size of the array. Can be specified only if the target
     !! array has to be of rank one.
     integer :: size
     !> identification of the allocation. Usually called as the name of the allocatable variable
     character(len=1) :: id
  end type doc

  interface f_map_ptr
     module procedure remap_bounds_d,remap_bounds_r,remap_bounds_l,remap_bounds_i,remap_bounds_li
  end interface

  !> Structure to perform heap allocation. Can be used to allocate allocatable arrays of any intrinsic
  !! type, kind and rank.
  !! The assignment of an allocatable array to this structure will allocate the array according to the
  !! specifications. It will also store the information about the allocation in an iternal database
  !! so that memory leaks and allocation problems can be more easily detected.
  !! @param mapname  @copydoc doc::size
  !! @param id       @copydoc doc::id
  interface f_malloc
     module procedure f_malloc,f_malloci_simple,f_mallocli_simple!,f_malloc_li
     module procedure f_malloc_bounds,f_malloc_bound
     !here also the procedures for the copying of arrays have to be defined
     module procedure f_malloc_b1,f_malloc_b2,f_malloc_b3
     module procedure f_malloc_l1,f_malloc_l2,f_malloc_l3
     module procedure f_malloc_i2,f_malloc_d2
     module procedure f_malloc_d1,f_malloc_i3
     module procedure f_malloc_d3,f_malloc_d4,f_malloc_d5
     module procedure f_malloc_r1,f_malloc_r2,f_malloc_r3,f_malloc_r4
  end interface

  interface f_malloc0
     module procedure f_malloc0,f_malloci0_simple,f_mallocli0_simple
     !!module procedure f_mallocli0
     module procedure f_malloc0_bounds,f_malloc0_bound
  end interface

  interface f_malloc_ptr
     module procedure f_malloc_ptr,f_malloci_ptr_simple,f_mallocli_ptr_simple
     module procedure f_malloc_ptr_bounds,f_malloc_ptr_bound
     !module procedure f_malloc_ptr_i2,f_malloc_ptr_i3
     module procedure f_malloc_ptr_i2_sp,f_malloc_ptr_i3_sp
     module procedure f_malloc_ptr_i4_sp
     !module procedure f_malloc_ptr_d1,f_malloc_ptr_d2
     !module procedure f_malloc_ptr_d3,f_malloc_ptr_d4
     module procedure f_malloc_ptr_d1_sp,f_malloc_ptr_d2_sp
     module procedure f_malloc_ptr_d3_sp,f_malloc_ptr_d4_sp,f_malloc_ptr_d5_sp
  end interface

  interface f_malloc0_ptr
     module procedure f_malloc0_ptr,f_malloci0_ptr_simple,f_mallocli0_ptr_simple
     !!module procedure f_mallocli0_ptr
     module procedure f_malloc0_ptr_bounds,f_malloc0_ptr_bound
  end interface

  interface f_malloc_str
     module procedure f_malloci_str,f_malloci_str_simple,f_mallocli_str_simple
     module procedure f_malloc_str_bounds,f_malloc_str_bound
     !!module procedure f_mallocli_str
     !!module procedure f_malloc_istr_simple
     !!module procedure f_malloc_listr_bounds
     !!module procedure f_malloc_listr_bound
  end interface

  interface f_malloc0_str
     module procedure f_malloci0_str,f_malloci0_str_simple,f_mallocli0_str_simple
     module procedure f_malloc0_str_bounds,f_malloc0_str_bound
     !!module procedure f_mallocli0_str
     !!module procedure f_malloc0_listr_bounds,f_malloc0_listr_bound
  end interface

  interface f_malloc_str_ptr
     module procedure f_malloci_str_ptr
     !!module procedure f_mallocli_str_ptr
     module procedure f_malloci_str_ptr_simple
     module procedure f_mallocli_str_ptr_simple
     module procedure f_malloc_str_ptr_bound
     module procedure f_malloc_str_ptr_bounds
  end interface

  interface f_malloc0_str_ptr
     module procedure f_malloci0_str_ptr
     !!module procedure f_mallocli0_str_ptr
     module procedure f_malloci0_str_ptr_simple
     module procedure f_malloc0_str_ptr_bounds,f_malloc0_str_ptr_bound
     !!module procedure f_malloc0_listr_ptr
     module procedure f_mallocli0_str_ptr_simple
  end interface

  !> Public routines
  public :: f_malloc,f_malloc0,f_malloc_ptr,f_malloc0_ptr,operator(.to.),operator(.plus.)
  public :: f_malloc_str,f_malloc0_str,f_malloc_str_ptr,f_malloc0_str_ptr
  public :: f_map_ptr,f_malloc_buf,f_malloc0_buf
  public :: remap_bounds_i2 !to verify if it works


contains

  pure function workspace_null() result(w)
    implicit none
    type(f_workspace) :: w
    call nullify_workspace(w)
  end function workspace_null
  pure subroutine nullify_workspace(w)
    implicit none
    type(f_workspace), intent(out) :: w
    w%pos_r=int(0,f_long)
    w%pos_d=int(0,f_long)
    w%pos_i=int(0,f_long)
    w%pos_li=int(0,f_long)
    w%pos_l=int(0,f_long)
    nullify(w%ptr_r)
    nullify(w%ptr_d)
    nullify(w%ptr_i)
    nullify(w%ptr_li)
    nullify(w%ptr_l)
  end subroutine nullify_workspace

  elemental pure function f_array_bounds_intint(nlow,nhigh) result(f_array_bounds)
    implicit none
    integer(f_integer), intent(in) :: nlow
    integer(f_integer), intent(in) :: nhigh
    type(array_bounds) :: f_array_bounds

    f_array_bounds%nlow=int(nlow,f_kind)
    f_array_bounds%nhigh=int(nhigh,f_kind)
  end function f_array_bounds_intint

  elemental pure function f_array_bounds_longlong(nlow,nhigh) result(f_array_bounds)
    implicit none
    integer(f_long), intent(in) :: nlow
    integer(f_long), intent(in) :: nhigh
    type(array_bounds) :: f_array_bounds

    f_array_bounds%nlow=int(nlow,f_kind)
    f_array_bounds%nhigh=int(nhigh,f_kind)
  end function f_array_bounds_longlong

  elemental pure function f_array_bounds_intlong(nlow,nhigh) result(f_array_bounds)
    implicit none
    integer(f_integer), intent(in) :: nlow
    integer(f_long), intent(in) :: nhigh
    type(array_bounds) :: f_array_bounds

    f_array_bounds%nlow=int(nlow,f_kind)
    f_array_bounds%nhigh=int(nhigh,f_kind)
  end function f_array_bounds_intlong

  elemental pure function f_array_bounds_longint(nlow,nhigh) result(f_array_bounds)
    implicit none
    integer(f_long), intent(in) :: nlow
    integer(f_integer), intent(in) :: nhigh
    type(array_bounds) :: f_array_bounds

    f_array_bounds%nlow=int(nlow,f_kind)
    f_array_bounds%nhigh=int(nhigh,f_kind)
  end function f_array_bounds_longint

  pure function f_array_concatenate_ii(a,b) result(c)
    implicit none
    integer(f_integer), dimension(:), intent(in) :: a,b
    integer(f_integer), dimension(size(a)+size(b)) :: c
    c(1:size(a))=a
    c(size(a)+1:size(a)+size(b))=b
  end function f_array_concatenate_ii

  pure function f_array_concatenate_ii0(a,b) result(c)
    implicit none
    integer(f_integer), intent(in) :: b
    integer(f_integer), dimension(:), intent(in) :: a
    integer(f_integer), dimension(size(a)+1) :: c
    c(1:size(a))=a
    c(size(a)+1)=b
  end function f_array_concatenate_ii0



  pure subroutine nullify_malloc_information_all(m)
    implicit none
    type(malloc_information_all), intent(out) :: m
    include 'f_malloc-null-inc.f90'
  end subroutine nullify_malloc_information_all

  pure subroutine nullify_malloc_information_ptr(m)
    implicit none
    type(malloc_information_ptr), intent(out) :: m
    include 'f_malloc-null-inc.f90'
    m%ptr=.true.
    nullify(m%w)
  end subroutine nullify_malloc_information_ptr

  pure subroutine nullify_malloc_information_buf(m)
    implicit none
    type(malloc_information_buf), intent(out) :: m
    include 'f_malloc-null-inc.f90'
  end subroutine nullify_malloc_information_buf

  pure subroutine nullify_malloc_information_str_all(m)
    implicit none
    type(malloc_information_str_all), intent(out) :: m
    include 'f_malloc-null-inc.f90'
    m%len=0
  end subroutine nullify_malloc_information_str_all

  pure subroutine nullify_malloc_information_str_ptr(m)
    implicit none
    type(malloc_information_str_ptr), intent(out) :: m
    include 'f_malloc-null-inc.f90'
    m%len=0
    m%ptr=.true.
  end subroutine nullify_malloc_information_str_ptr

  !> f95-compliant routine to remap pointer bounds, as suggested from (as of Sep. 2015)
  !! https://en.wikipedia.org/wiki/Fortran_95_language_features#Pointers_as_dynamic_aliases
  !! What has to be verified if compiler perform workarrays constructions
  subroutine remap_bounds_d(lb,lu,heap,ptr)
    implicit none
    integer(f_kind), intent(in) :: lb,lu
    real(f_double), dimension(lb:lu), intent(in), target :: heap
    real(f_double), dimension(:), pointer, intent(out) :: ptr
    ptr => heap
  end subroutine remap_bounds_d

  !> the same routine but with shape
  subroutine remap_bounds_d2(lb,lu,heap,ptr)
    implicit none
    integer(f_kind), dimension(2), intent(in) :: lb,lu
    real(f_double), dimension(lb(1):lu(1),lb(2):lu(2)), intent(in), target :: heap
    real(f_double), dimension(:,:), pointer, intent(out) :: ptr
    ptr => heap
  end subroutine remap_bounds_d2


  subroutine remap_bounds_r(lb,lu,heap,ptr)
    implicit none
    integer(f_kind), intent(in) :: lb,lu
    real(f_simple), dimension(lb:lu), intent(in), target :: heap
    real(f_simple), dimension(:), pointer, intent(out) :: ptr
    ptr => heap
  end subroutine remap_bounds_r

  subroutine remap_bounds_i(lb,lu,heap,ptr)
    implicit none
    integer(f_kind), intent(in) :: lb,lu
    integer(f_integer), dimension(lb:lu), intent(in), target :: heap
    integer(f_integer), dimension(:), pointer, intent(out) :: ptr
    ptr => heap
  end subroutine remap_bounds_i

  subroutine remap_bounds_li(lb,lu,heap,ptr)
    implicit none
    integer(f_kind), intent(in) :: lb,lu
    integer(f_long), dimension(lb:lu), intent(in), target :: heap
    integer(f_long), dimension(:), pointer, intent(out) :: ptr
    ptr => heap
  end subroutine remap_bounds_li

  subroutine remap_bounds_l(lb,lu,heap,ptr)
    implicit none
    integer(f_kind), intent(in) :: lb,lu
    logical, dimension(lb:lu), intent(in), target :: heap
    logical, dimension(:), pointer, intent(out) :: ptr
    ptr => heap
  end subroutine remap_bounds_l


  subroutine remap_bounds_i2(lb,lu,heap,ptr)
    implicit none
    integer(f_kind), dimension(2), intent(in) :: lb,lu
    integer(f_integer), dimension(lb(1):lu(1),lb(2):lu(2)), intent(in), target :: heap
    integer(f_integer), dimension(:,:), pointer, intent(out) :: ptr
    ptr => heap
  end subroutine remap_bounds_i2



  !---routines for low-level dynamic memory handling

  !> For rank-1 arrays
  pure function f_malloci_simple(size,id,routine_id,profile,info) result(m)
    implicit none
    type(malloc_information_all) :: m
    include 'f_malloci-simple-inc.f90'
  end function f_malloci_simple
  pure function f_mallocli_simple(size,id,routine_id,profile,info) result(m)
    implicit none
    type(malloc_information_all) :: m
    include 'f_mallocli-simple-inc.f90'
  end function f_mallocli_simple
!!$  function f_malloc_li(sizes,id,routine_id,lbounds,ubounds,profile,info,src,src_ptr) result(m)
!!$    implicit none
!!$    !the integer array src is here added to avoid problems in resolving the ambiguity with f_malloc_src
!!$    integer(f_long), dimension(:), intent(in), optional :: src
!!$    integer(f_long), dimension(:), pointer, intent(in), optional :: src_ptr
!!$    type(malloc_information_all) :: m
!!$    integer(f_long), dimension(:), intent(in), optional :: sizes,lbounds,ubounds
!!$    !local variables
!!$    integer :: i
!!$    include 'f_malloc-base-inc.f90'
!!$    include 'f_malloc-check-inc.f90'
!!$    include 'f_malloc-extra-inc.f90'
!!$  end function f_malloc_li

  !> For rank-1 arrays
  pure function f_malloci0_simple(size,id,routine_id,profile,info) result(m)
    implicit none
    type(malloc_information_all) :: m
    include 'f_malloci-simple-inc.f90'
    m%put_to_zero=.true.
  end function f_malloci0_simple
  pure function f_mallocli0_simple(size,id,routine_id,profile,info) result(m)
    implicit none
    type(malloc_information_all) :: m
    include 'f_mallocli-simple-inc.f90'
    m%put_to_zero=.true.
  end function f_mallocli0_simple



  pure function f_malloci_ptr_simple(size,id,routine_id,profile,info) result(m)
    implicit none
    type(malloc_information_ptr) :: m
    include 'f_malloci-simple-inc.f90'
  end function f_malloci_ptr_simple
  pure function f_mallocli_ptr_simple(size,id,routine_id,profile,info) result(m)
    implicit none
    type(malloc_information_ptr) :: m
    include 'f_mallocli-simple-inc.f90'
  end function f_mallocli_ptr_simple
  !> For rank-1 arrays
  pure function f_malloci0_ptr_simple(size,id,routine_id,profile,info) result(m)
    implicit none
    type(malloc_information_ptr) :: m
    include 'f_malloci-simple-inc.f90'
    m%put_to_zero=.true.
  end function f_malloci0_ptr_simple
  pure function f_mallocli0_ptr_simple(size,id,routine_id,profile,info) result(m)
    implicit none
    type(malloc_information_ptr) :: m
    include 'f_mallocli-simple-inc.f90'
    m%put_to_zero=.true.
  end function f_mallocli0_ptr_simple
  !> For rank-1 arrays
  pure function f_malloci_str_simple(length,size,id,routine_id,profile,info) result(m)
    implicit none
    type(malloc_information_str_all) :: m
    integer, intent(in) :: length
    include 'f_malloci-simple-inc.f90'
    m%len=length
  end function f_malloci_str_simple
  pure function f_mallocli_str_simple(length,size,id,routine_id,profile,info) result(m)
    implicit none
    type(malloc_information_str_all) :: m
    integer, intent(in) :: length
    include 'f_mallocli-simple-inc.f90'
    m%len=length
  end function f_mallocli_str_simple
  !> For rank-1 arrays
  pure function f_malloci0_str_simple(length,size,id,routine_id,profile,info) result(m)
    implicit none
    type(malloc_information_str_all) :: m
    integer, intent(in) :: length
    include 'f_malloci-simple-inc.f90'
    m%len=length
    m%put_to_zero=.true.
  end function f_malloci0_str_simple
  pure function f_mallocli0_str_simple(length,size,id,routine_id,profile,info) result(m)
    implicit none
    type(malloc_information_str_all) :: m
    integer, intent(in) :: length
    include 'f_mallocli-simple-inc.f90'
    m%len=length
    m%put_to_zero=.true.
  end function f_mallocli0_str_simple
  !> For rank-1 arrays
  pure function f_malloci_str_ptr_simple(length,size,id,routine_id,profile,info) result(m)
    implicit none
    type(malloc_information_str_ptr) :: m
    integer, intent(in) :: length
    include 'f_malloci-simple-inc.f90'
    m%len=length
  end function f_malloci_str_ptr_simple
  pure function f_mallocli_str_ptr_simple(length,size,id,routine_id,profile,info) result(m)
    implicit none
    type(malloc_information_str_ptr) :: m
    integer, intent(in) :: length
    include 'f_mallocli-simple-inc.f90'
    m%len=length
  end function f_mallocli_str_ptr_simple
  !> For rank-1 arrays
  pure function f_malloci0_str_ptr_simple(length,size,id,routine_id,profile,info) result(m)
    implicit none
    type(malloc_information_str_ptr) :: m
    integer, intent(in) :: length
    include 'f_malloci-simple-inc.f90'
    m%len=length
    m%put_to_zero=.true.
  end function f_malloci0_str_ptr_simple
  pure function f_mallocli0_str_ptr_simple(length,size,id,routine_id,profile,info) result(m)
    implicit none
    type(malloc_information_str_ptr) :: m
    integer, intent(in) :: length
    include 'f_mallocli-simple-inc.f90'
    m%len=length
    m%put_to_zero=.true.
  end function f_mallocli0_str_ptr_simple

  !> For rank-1 arrays, with bounds
  pure function f_malloc_bound(bounds,id,routine_id,profile,info) result(m)
    implicit none
    type(malloc_information_all) :: m
    include 'f_malloc-bound-inc.f90'
  end function f_malloc_bound
  !>for rank-1 arrays, with boundaries
  pure function f_malloc0_bound(bounds,id,routine_id,profile,info) result(m)
    implicit none
    type(malloc_information_all) :: m
    include 'f_malloc-bound-inc.f90'
    m%put_to_zero=.true.
  end function f_malloc0_bound
  !> For rank-1 arrays
  pure function f_malloc_ptr_bound(bounds,id,routine_id,profile,info) result(m)
    implicit none
    type(malloc_information_ptr) :: m
    include 'f_malloc-bound-inc.f90'
  end function f_malloc_ptr_bound
  !> For rank-1 arrays
  pure function f_malloc0_ptr_bound(bounds,id,routine_id,profile,info) result(m)
    implicit none
    type(malloc_information_ptr) :: m
    include 'f_malloc-bound-inc.f90'
    m%put_to_zero=.true.
  end function f_malloc0_ptr_bound
  !> For rank-1 arrays, with bounds
  pure function f_malloc_str_bound(length,bounds,id,routine_id,profile,info) result(m)
    implicit none
    type(malloc_information_str_all) :: m
    integer, intent(in) :: length
    include 'f_malloc-bound-inc.f90'
    m%len=length
  end function f_malloc_str_bound
  !!pure function f_malloc_listr_bound(length,bounds,id,routine_id,profile,info) result(m)
  !!  implicit none
  !!  type(malloc_information_str_all) :: m
  !!  integer(kind=4), intent(in) :: length
  !!  include 'f_malloc-bound-inc.f90'
  !!  m%len=int(length,kind=4)
  !!end function f_malloc_listr_bound
  !>for rank-1 arrays, with boundaries
  pure function f_malloc0_str_bound(length,bounds,id,routine_id,profile,info) result(m)
    implicit none
    type(malloc_information_str_all) :: m
    integer, intent(in) :: length
    include 'f_malloc-bound-inc.f90'
    m%len=length
    m%put_to_zero=.true.
  end function f_malloc0_str_bound
  !!pure function f_malloc0_listr_bound(length,bounds,id,routine_id,profile,info) result(m)
  !!  implicit none
  !!  type(malloc_information_str_all) :: m
  !!  integer(kind=4), intent(in) :: length
  !!  include 'f_malloc-bound-inc.f90'
  !!  m%len=int(length,kind=4)
  !!  m%put_to_zero=.true.
  !!end function f_malloc0_listr_bound
  !> For rank-1 arrays
  pure function f_malloc_str_ptr_bound(length,bounds,id,routine_id,profile,info) result(m)
    implicit none
    type(malloc_information_str_ptr) :: m
    integer, intent(in) :: length
    include 'f_malloc-bound-inc.f90'
    m%len=length
  end function f_malloc_str_ptr_bound
  !> For rank-1 arrays
  pure function f_malloc0_str_ptr_bound(length,bounds,id,routine_id,profile,info) result(m)
    implicit none
    type(malloc_information_str_ptr) :: m
    integer, intent(in) :: length
    include 'f_malloc-bound-inc.f90'
    m%len=length
    m%put_to_zero=.true.
  end function f_malloc0_str_ptr_bound
  pure function f_malloc0_listr_ptr_bound(length,bounds,id,routine_id,profile,info) result(m)
    implicit none
    type(malloc_information_str_ptr) :: m
    integer, intent(in) :: length
    include 'f_malloc-bound-inc.f90'
    m%len=length
    m%put_to_zero=.true.
  end function f_malloc0_listr_ptr_bound


  !> Define the allocation information for  arrays of different rank
  pure function f_malloc_bounds(bounds,id,routine_id,profile,info) result(m)
    implicit none
    type(malloc_information_all) :: m
    include 'f_malloc-bounds-inc.f90'
  end function f_malloc_bounds
  pure function f_malloc0_bounds(bounds,id,routine_id,profile,info) result(m)
    implicit none
    type(malloc_information_all) :: m
    include 'f_malloc-bounds-inc.f90'
    m%put_to_zero=.true.
  end function f_malloc0_bounds
  !> Define the allocation information for  arrays of different rank
  pure function f_malloc_ptr_bounds(bounds,id,routine_id,profile,info) result(m)
    implicit none
    type(malloc_information_ptr) :: m
    include 'f_malloc-bounds-inc.f90'
  end function f_malloc_ptr_bounds
  !> Define the allocation information for  arrays of different rank
  pure function f_malloc0_ptr_bounds(bounds,id,routine_id,profile,info) result(m)
    implicit none
    type(malloc_information_ptr) :: m
    include 'f_malloc-bounds-inc.f90'
    m%put_to_zero=.true.
  end function f_malloc0_ptr_bounds
  !> Define the allocation information for  arrays of different rank
  pure function f_malloc_str_bounds(length,bounds,id,routine_id,profile,info) result(m)
    implicit none
    type(malloc_information_str_all) :: m
    integer, intent(in) :: length
    include 'f_malloc-bounds-inc.f90'
    m%len=length
  end function f_malloc_str_bounds
  !!pure function f_malloc_listr_bounds(length,bounds,id,routine_id,profile,info) result(m)
  !!  implicit none
  !!  type(malloc_information_str_all) :: m
  !!  integer(kind=4), intent(in) :: length
  !!  include 'f_malloc-bounds-inc.f90'
  !!  m%len=int(length,kind=4)
  !!end function f_malloc_listr_bounds
  pure function f_malloc0_str_bounds(length,bounds,id,routine_id,profile,info) result(m)
    implicit none
    type(malloc_information_str_all) :: m
    integer, intent(in) :: length
    include 'f_malloc-bounds-inc.f90'
    m%len=length
    m%put_to_zero=.true.
  end function f_malloc0_str_bounds
  !!pure function f_malloc0_listr_bounds(length,bounds,id,routine_id,profile,info) result(m)
  !!  implicit none
  !!  type(malloc_information_str_all) :: m
  !!  integer(kind=4), intent(in) :: length
  !!  include 'f_malloc-bounds-inc.f90'
  !!  m%len=int(length,kind=4)
  !!  m%put_to_zero=.true.
  !!end function f_malloc0_listr_bounds
  !> Define the allocation information for  arrays of different rank
  pure function f_malloc_str_ptr_bounds(length,bounds,id,routine_id,profile,info) result(m)
    implicit none
    type(malloc_information_str_ptr) :: m
    integer, intent(in) :: length
    include 'f_malloc-bounds-inc.f90'
    m%len=length
  end function f_malloc_str_ptr_bounds
  !> Define the allocation information for  arrays of different rank
  pure function f_malloc0_str_ptr_bounds(length,bounds,id,routine_id,profile,info) result(m)
    implicit none
    type(malloc_information_str_ptr) :: m
    integer, intent(in) :: length
    include 'f_malloc-bounds-inc.f90'
    m%len=length
    m%put_to_zero=.true.
  end function f_malloc0_str_ptr_bounds

  !> Define the allocation information for arrays of different rank
  function f_malloc(sizes,id,routine_id,lbounds,ubounds,profile,info,src,src_ptr) result(m)
    implicit none
    !the integer array src is here added to avoid problems in resolving the ambiguity with f_malloc_src
    integer(f_integer), dimension(:), intent(in), optional :: src
    integer(f_integer), dimension(:), pointer, intent(in), optional :: src_ptr
    type(malloc_information_all) :: m
    integer(f_integer), dimension(:), intent(in), optional :: sizes,lbounds,ubounds
    !local variables
    integer :: i
    include 'f_malloc-base-inc.f90'
    include 'f_malloc-check-inc.f90'
    include 'f_malloc-extra-inc.f90'
  end function f_malloc
  !> define the allocation information for  arrays of different rank
  function f_malloc0(sizes,id,routine_id,lbounds,ubounds,profile,info) result(m)
    implicit none
    type(malloc_information_all) :: m
    include 'f_malloc-total-inc.f90'
    m%put_to_zero=.true.
  end function f_malloc0
  !!function f_mallocli0(sizes,id,routine_id,lbounds,ubounds,profile,info) result(m)
  !!  implicit none
  !!  type(malloc_information_all) :: m
  !!  include 'f_mallocli-total-inc.f90'
  !!  m%put_to_zero=.true.
  !!end function f_mallocli0
  !> Define the allocation information for  arrays of different rank
  function f_malloc_ptr(sizes,id,routine_id,lbounds,ubounds,profile,info,src,src_ptr) result(m)
    implicit none
    !the integer array src is here added to avoid problems in resolving the ambiguity
    integer, dimension(:), intent(in), optional :: src
    integer, dimension(:), pointer, intent(in), optional :: src_ptr
    type(malloc_information_ptr) :: m
    integer, dimension(:), intent(in), optional :: sizes,lbounds,ubounds
    !local variables
    integer :: i
    include 'f_malloc-base-inc.f90'
    include 'f_malloc-check-inc.f90'
    include 'f_malloc-extra-inc.f90'
  end function f_malloc_ptr
  !> Define the allocation information for  arrays of different rank
  function f_malloc0_ptr(sizes,id,routine_id,lbounds,ubounds,profile,info) result(m)
    implicit none
    type(malloc_information_ptr) :: m
    include 'f_malloc-total-inc.f90'
    m%put_to_zero=.true.
  end function f_malloc0_ptr
  !!function f_mallocli0_ptr(sizes,id,routine_id,lbounds,ubounds,profile,info) result(m)
  !!  implicit none
  !!  type(malloc_information_ptr) :: m
  !!  include 'f_mallocli-total-inc.f90'
  !!  m%put_to_zero=.true.
  !!end function f_mallocli0_ptr
  !> Define the allocation information for  arrays of different rank
  function f_malloci_str(length,sizes,id,routine_id,lbounds,ubounds,profile,info) result(m)
    implicit none
    type(malloc_information_str_all) :: m
    integer, intent(in) :: length
    include 'f_malloc-total-inc.f90'
    m%len=length
  end function f_malloci_str

  !> define the allocation information for  arrays of different rank
  function f_malloci0_str(length,sizes,id,routine_id,lbounds,ubounds,profile,info) result(m)
    implicit none
    type(malloc_information_str_all) :: m
    integer, intent(in) :: length
    include 'f_malloc-total-inc.f90'
    m%len=length
    m%put_to_zero=.true.
  end function f_malloci0_str

  !> Define the allocation information for  arrays of different rank
  function f_malloci_str_ptr(length,sizes,id,routine_id,lbounds,ubounds,profile,info) result(m)
    implicit none
    type(malloc_information_str_ptr) :: m
    integer, intent(in) :: length
    include 'f_malloc-total-inc.f90'
    m%len=length
  end function f_malloci_str_ptr

  !> Define the allocation information for  arrays of different rank
  function f_malloci0_str_ptr(length,sizes,id,routine_id,lbounds,ubounds,profile,info) result(m)
    implicit none
    type(malloc_information_str_ptr) :: m
    integer, intent(in) :: length
    include 'f_malloc-total-inc.f90'
    m%len=length
    m%put_to_zero=.true.
  end function f_malloci0_str_ptr

  function f_malloc_d1(src,lbounds,ubounds,id,routine_id,profile,info) result(m)
    implicit none
    double precision, dimension(:), intent(in) :: src
    integer, dimension(:), intent(in), optional :: lbounds,ubounds
    type(malloc_information_all) :: m
    include 'f_malloc-base-inc.f90'
    include 'f_malloc-inc.f90'
  end function f_malloc_d1

  function f_malloc_d2(src,id,routine_id,lbounds,ubounds,profile,info) result(m)
    implicit none
    double precision, dimension(:,:), intent(in) :: src
    integer, dimension(:), intent(in), optional :: lbounds,ubounds
    type(malloc_information_all) :: m
    include 'f_malloc-base-inc.f90'
    include 'f_malloc-inc.f90'
  end function f_malloc_d2

  function f_malloc_d3(src,id,routine_id,lbounds,ubounds,profile,info) result(m)
    implicit none
    double precision, dimension(:,:,:), intent(in) :: src
    integer, dimension(:), intent(in), optional :: lbounds,ubounds
    type(malloc_information_all) :: m
    include 'f_malloc-base-inc.f90'
    include 'f_malloc-inc.f90'
  end function f_malloc_d3

  function f_malloc_d4(src,id,routine_id,lbounds,ubounds,profile,info) result(m)
    implicit none
    double precision, dimension(:,:,:,:), intent(in) :: src
    integer, dimension(:), intent(in), optional :: lbounds,ubounds
    type(malloc_information_all) :: m
    include 'f_malloc-base-inc.f90'
    include 'f_malloc-inc.f90'
  end function f_malloc_d4

  function f_malloc_d5(src,id,routine_id,lbounds,ubounds,profile,info) result(m)
    implicit none
    double precision, dimension(:,:,:,:,:), intent(in) :: src
    integer, dimension(:), intent(in), optional :: lbounds,ubounds
    type(malloc_information_all) :: m
    include 'f_malloc-base-inc.f90'
    include 'f_malloc-inc.f90'
  end function f_malloc_d5

  function f_malloc_r1(src,lbounds,ubounds,id,routine_id,profile,info) result(m)
    implicit none
    real, dimension(:), intent(in) :: src
    integer, dimension(:), intent(in), optional :: lbounds,ubounds
    type(malloc_information_all) :: m
    include 'f_malloc-base-inc.f90'
    include 'f_malloc-inc.f90'
  end function f_malloc_r1

  function f_malloc_b1(src,lbounds,ubounds,id,routine_id,profile,info) result(m)
    implicit none
    logical(f_byte), dimension(:), intent(in) :: src
    integer, dimension(:), intent(in), optional :: lbounds,ubounds
    type(malloc_information_all) :: m
    include 'f_malloc-base-inc.f90'
    include 'f_malloc-inc.f90'
  end function f_malloc_b1

  function f_malloc_b2(src,lbounds,ubounds,id,routine_id,profile,info) result(m)
    implicit none
    logical(f_byte), dimension(:,:), intent(in) :: src
    integer, dimension(:), intent(in), optional :: lbounds,ubounds
    type(malloc_information_all) :: m
    include 'f_malloc-base-inc.f90'
    include 'f_malloc-inc.f90'
  end function f_malloc_b2

  function f_malloc_b3(src,lbounds,ubounds,id,routine_id,profile,info) result(m)
    implicit none
    logical(f_byte), dimension(:,:,:), intent(in) :: src
    integer, dimension(:), intent(in), optional :: lbounds,ubounds
    type(malloc_information_all) :: m
    include 'f_malloc-base-inc.f90'
    include 'f_malloc-inc.f90'
  end function f_malloc_b3

  function f_malloc_r2(src,id,routine_id,lbounds,ubounds,profile,info) result(m)
    implicit none
    real, dimension(:,:), intent(in) :: src
    integer, dimension(:), intent(in), optional :: lbounds,ubounds
    type(malloc_information_all) :: m
    include 'f_malloc-base-inc.f90'
    include 'f_malloc-inc.f90'
  end function f_malloc_r2

  function f_malloc_r3(src,id,routine_id,lbounds,ubounds,profile,info) result(m)
    implicit none
    real, dimension(:,:,:), intent(in) :: src
    integer, dimension(:), intent(in), optional :: lbounds,ubounds
    type(malloc_information_all) :: m
    include 'f_malloc-base-inc.f90'
    include 'f_malloc-inc.f90'
  end function f_malloc_r3

  function f_malloc_r4(src,id,routine_id,lbounds,ubounds,profile,info) result(m)
    implicit none
    real, dimension(:,:,:,:), intent(in) :: src
    integer, dimension(:), intent(in), optional :: lbounds,ubounds
    type(malloc_information_all) :: m
    include 'f_malloc-base-inc.f90'
    include 'f_malloc-inc.f90'
  end function f_malloc_r4

  function f_malloc_i2(src,id,routine_id,lbounds,ubounds,profile,info) result(m)
    implicit none
    integer, dimension(:,:), intent(in) :: src
    integer, dimension(:), intent(in), optional :: lbounds,ubounds
    type(malloc_information_all) :: m
    include 'f_malloc-base-inc.f90'
    include 'f_malloc-inc.f90'
  end function f_malloc_i2

  function f_malloc_i3(src,id,routine_id,lbounds,ubounds,profile,info) result(m)
    implicit none
    integer, dimension(:,:,:), intent(in) :: src
    integer, dimension(:), intent(in), optional :: lbounds,ubounds
    type(malloc_information_all) :: m
    include 'f_malloc-base-inc.f90'
    include 'f_malloc-inc.f90'
  end function f_malloc_i3

  function f_malloc_l1(src,lbounds,ubounds,id,routine_id,profile,info) result(m)
    implicit none
    logical, dimension(:), intent(in) :: src
    integer, dimension(:), intent(in), optional :: lbounds,ubounds
    type(malloc_information_all) :: m
    include 'f_malloc-base-inc.f90'
    include 'f_malloc-inc.f90'
  end function f_malloc_l1

  function f_malloc_l2(src,lbounds,ubounds,id,routine_id,profile,info) result(m)
    implicit none
    logical, dimension(:,:), intent(in) :: src
    integer, dimension(:), intent(in), optional :: lbounds,ubounds
    type(malloc_information_all) :: m
    include 'f_malloc-base-inc.f90'
    include 'f_malloc-inc.f90'
  end function f_malloc_l2

  function f_malloc_l3(src,id,routine_id,lbounds,ubounds,profile,info) result(m)
    implicit none
    logical, dimension(:,:,:), intent(in) :: src
    integer, dimension(:), intent(in), optional :: lbounds,ubounds
    type(malloc_information_all) :: m
    include 'f_malloc-base-inc.f90'
    include 'f_malloc-inc.f90'
  end function f_malloc_l3

!!$  function f_malloc_ptr_i2(src,id,routine_id,sizes,lbounds,ubounds,profile,info) result(m)
!!$    implicit none
!!$    integer, dimension(:,:), intent(in) :: src
!!$    integer, dimension(:), intent(in), optional :: sizes,lbounds,ubounds
!!$    type(malloc_information_ptr) :: m
!!$    include 'f_malloc-base-inc.f90'
!!$    include 'f_malloc-inc.f90'
!!$  end function f_malloc_ptr_i2

  !this template is here waiting for a new unambiguous module procedure
  function f_malloc_ptr_i2_sp(src_ptr,id,routine_id,profile,info) result(m)
    implicit none
    integer, dimension(:,:), pointer, intent(in) :: src_ptr
    type(malloc_information_ptr) :: m
    include 'f_malloc-base-inc.f90'
    include 'f_malloc-null-ptr-inc.f90'
    include 'f_malloc-ptr-inc.f90'
  end function f_malloc_ptr_i2_sp

!!$  function f_malloc_ptr_i3(src,id,routine_id,sizes,lbounds,ubounds,profile,info) result(m)
!!$    implicit none
!!$    integer, dimension(:,:,:), intent(in) :: src
!!$    integer, dimension(:), intent(in), optional :: sizes,lbounds,ubounds
!!$    type(malloc_information_ptr) :: m
!!$    include 'f_malloc-base-inc.f90'
!!$    include 'f_malloc-inc.f90'
!!$  end function f_malloc_ptr_i3

  !this template is here waiting for a new unambiguous module procedure
  function f_malloc_ptr_i3_sp(src_ptr,id,routine_id,profile,info) result(m)
    implicit none
    integer, dimension(:,:,:), pointer, intent(in) :: src_ptr
    type(malloc_information_ptr) :: m
    include 'f_malloc-base-inc.f90'
    include 'f_malloc-null-ptr-inc.f90'
    include 'f_malloc-ptr-inc.f90'
  end function f_malloc_ptr_i3_sp

  !this template is here waiting for a new unambiguous module procedure
  function f_malloc_ptr_i4_sp(src_ptr,id,routine_id,profile,info) result(m)
    implicit none
    integer, dimension(:,:,:,:), pointer, intent(in) :: src_ptr
    type(malloc_information_ptr) :: m
    include 'f_malloc-base-inc.f90'
    include 'f_malloc-null-ptr-inc.f90'
    include 'f_malloc-ptr-inc.f90'
  end function f_malloc_ptr_i4_sp

  !this template is here waiting for a new unambiguous module procedure
  function f_malloc_ptr_d1_sp(src_ptr,id,routine_id,profile,info) result(m)
    implicit none
    real(f_double), dimension(:), pointer, intent(in) :: src_ptr
    type(malloc_information_ptr) :: m
    include 'f_malloc-base-inc.f90'
    include 'f_malloc-null-ptr-inc.f90'
    include 'f_malloc-ptr-inc.f90'
  end function f_malloc_ptr_d1_sp

  !this template is here waiting for a new unambiguous module procedure
  function f_malloc_ptr_d2_sp(src_ptr,id,routine_id,profile,info) result(m)
    implicit none
    real(f_double), dimension(:,:), pointer, intent(in) :: src_ptr
    type(malloc_information_ptr) :: m
    include 'f_malloc-base-inc.f90'
    include 'f_malloc-null-ptr-inc.f90'
    include 'f_malloc-ptr-inc.f90'
  end function f_malloc_ptr_d2_sp

  !this template is here waiting for a new unambiguous module procedure
  function f_malloc_ptr_d3_sp(src_ptr,id,routine_id,profile,info) result(m)
    implicit none
    double precision, dimension(:,:,:), pointer, intent(in) :: src_ptr
    type(malloc_information_ptr) :: m
    include 'f_malloc-base-inc.f90'
    include 'f_malloc-null-ptr-inc.f90'
    include 'f_malloc-ptr-inc.f90'
  end function f_malloc_ptr_d3_sp

  !this template is here waiting for a new unambiguous module procedure
  function f_malloc_ptr_d4_sp(src_ptr,id,routine_id,profile,info) result(m)
    implicit none
    double precision, dimension(:,:,:,:), pointer, intent(in) :: src_ptr
    type(malloc_information_ptr) :: m
    include 'f_malloc-base-inc.f90'
    include 'f_malloc-null-ptr-inc.f90'
    include 'f_malloc-ptr-inc.f90'
  end function f_malloc_ptr_d4_sp

  function f_malloc_ptr_d5_sp(src_ptr,id,routine_id,profile,info) result(m)
    implicit none
    double precision, dimension(:,:,:,:,:), pointer, intent(in) :: src_ptr
    type(malloc_information_ptr) :: m
    include 'f_malloc-base-inc.f90'
    include 'f_malloc-null-ptr-inc.f90'
    include 'f_malloc-ptr-inc.f90'
  end function f_malloc_ptr_d5_sp

  function f_malloc_buf(sizes,id,routine_id,profile,info,src) result(m)
    implicit none
    !the integer array src is here added to avoid problems in resolving the ambiguity
    integer, dimension(:), intent(in), optional :: src
    include 'f_malloc-buf-base-inc.f90'
    include 'f_malloc-base-inc.f90'
    include 'f_malloc-buf-inc.f90'
  end function f_malloc_buf
  !> Define the allocation information for  arrays of different rank
  function f_malloc0_buf(sizes,id,routine_id,profile,info) result(m)
    implicit none
    include 'f_malloc-buf-base-inc.f90'
    include 'f_malloc-base-inc.f90'
    m%rank=1
    m%shape(1:m%rank)=int(sizes,f_kind)
    m%ubounds(1:m%rank)=int(m%lbounds(1:m%rank),f_kind)+m%shape(1:m%rank)-1
    m%put_to_zero=.true.
  end function f_malloc0_buf


!!$  function f_malloc_ptr_d1(id,src,routine_id,sizes,lbounds,ubounds,profile,info) result(m)
!!$    implicit none
!!$    double precision, dimension(:), optional, intent(in) :: src
!!$    integer, dimension(:), intent(in), optional :: sizes,lbounds,ubounds
!!$    type(malloc_information_ptr) :: m
!!$    include 'f_malloc-base-inc.f90'
!!$    include 'f_malloc-inc.f90'
!!$  end function f_malloc_ptr_d1

!!$  function f_malloc_ptr_d2(id,src,routine_id,sizes,lbounds,ubounds,profile,info) result(m)
!!$    implicit none
!!$    double precision, dimension(:,:), optional, intent(in) :: src
!!$    integer, dimension(:), intent(in), optional :: sizes,lbounds,ubounds
!!$    type(malloc_information_ptr) :: m
!!$    include 'f_malloc-base-inc.f90'
!!$    include 'f_malloc-inc.f90'
!!$  end function f_malloc_ptr_d2
!!$
!!$  function f_malloc_ptr_d3(id,src,routine_id,sizes,lbounds,ubounds,profile,info) result(m)
!!$    implicit none
!!$    double precision, dimension(:,:,:), optional, intent(in) :: src
!!$    integer, dimension(:), intent(in), optional :: sizes,lbounds,ubounds
!!$    type(malloc_information_ptr) :: m
!!$    include 'f_malloc-base-inc.f90'
!!$    include 'f_malloc-inc.f90'
!!$  end function f_malloc_ptr_d3
!!$
!!$  function f_malloc_ptr_d4(id,src,routine_id,sizes,lbounds,ubounds,profile,info) result(m)
!!$    implicit none
!!$    double precision, dimension(:,:,:,:), optional, intent(in) :: src
!!$    integer, dimension(:), intent(in), optional :: sizes,lbounds,ubounds
!!$    type(malloc_information_ptr) :: m
!!$    include 'f_malloc-base-inc.f90'
!!$    include 'f_malloc-inc.f90'
!!$  end function f_malloc_ptr_d4

end module module_f_malloc

!> trick to use the element of a pointer as an address for subptr
subroutine f_map_ptr_addr_d0(lb,ub,is,ie,heap,ptr)
  use module_f_malloc, only: f_map_ptr,f_kind
  use f_precisions, only: f_double,f_loc
  use dictionaries, only: f_err_throw
  implicit none
  integer(f_kind) :: lb,ub,is,ie
  real(f_double), dimension(*) :: heap
  real(f_double), dimension(:), pointer :: ptr

  call f_map_ptr(lb,ub,heap(is:ie),ptr)
  if (f_loc(ptr(lb)) /= f_loc(heap(is)) .or. &
       f_loc(ptr(ub)) /= f_loc(heap(ie))) call f_err_throw(&
       'ERROR (f_subptr): addresses do not match, the allocating system has performed a copy',&
       err_name='ERR_MALLOC_INTERNAL')

end subroutine f_map_ptr_addr_d0


subroutine f_map_ptr_addr_i0(lb,ub,is,ie,heap,ptr)
  use module_f_malloc, only: f_map_ptr,f_kind
  use f_precisions, only: f_integer,f_loc
  use dictionaries, only: f_err_throw
  implicit none
  integer(f_kind) :: lb,ub,is,ie
  integer(f_integer), dimension(*) :: heap
  integer(f_integer), dimension(:), pointer :: ptr

  call f_map_ptr(lb,ub,heap(is:ie),ptr)
  if (f_loc(ptr(lb)) /= f_loc(heap(is)) .or. &
       f_loc(ptr(ub)) /= f_loc(heap(ie))) call f_err_throw(&
       'ERROR (f_subptr): addresses do not match, the allocating system has performed a copy',&
       err_name='ERR_MALLOC_INTERNAL')

end subroutine f_map_ptr_addr_i0

subroutine f_map_ptr_addr_li0(lb,ub,is,ie,heap,ptr)
  use module_f_malloc, only: f_map_ptr,f_kind
  use f_precisions, only: f_long,f_loc
  use dictionaries, only: f_err_throw
  implicit none
  integer(f_kind) :: lb,ub,is,ie
  integer(f_long), dimension(*) :: heap
  integer(f_long), dimension(:), pointer :: ptr

  call f_map_ptr(lb,ub,heap(is:ie),ptr)
  if (f_loc(ptr(lb)) /= f_loc(heap(is)) .or. &
       f_loc(ptr(ub)) /= f_loc(heap(ie))) call f_err_throw(&
       'ERROR (f_subptr): addresses do not match, the allocating system has performed a copy',&
       err_name='ERR_MALLOC_INTERNAL')

end subroutine f_map_ptr_addr_li0
