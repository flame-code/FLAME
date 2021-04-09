module f_python
  use f_precisions, only: f_address

  implicit none

  !> Equivalent type than the numpy one, to be able
  !  to export a Fortran array into Python space.
  type ndarray
     integer(f_address) :: data
     integer :: ndims
     integer, dimension(7) :: shapes
     character(len = 2) :: kind
  end type ndarray

  interface toNdArray
     module procedure i0_to_ndarray, i1_to_ndarray, i2_to_ndarray
     module procedure f0_to_ndarray, f1_to_ndarray, f2_to_ndarray
  end interface toNdArray

  interface toNdArray_ptr
     module procedure pi0_to_ndarray, pi1_to_ndarray, pi2_to_ndarray
     module procedure pf0_to_ndarray, pf1_to_ndarray, pf2_to_ndarray
  end interface toNdArray_ptr

  interface
     subroutine f_python_initialize(iproc, nproc, igroup, ngroup)
       implicit none
       integer, intent(in) :: iproc, nproc, igroup, ngroup
     end subroutine f_python_initialize
  end interface

  interface
     subroutine f_python_finalize()
       implicit none
     end subroutine f_python_finalize
  end interface

  interface
     subroutine f_python_execute_dict(dict, status)
       use dictionaries_base, only : dictionary
       implicit none
       type(dictionary), pointer :: dict
       integer, intent(out) :: status
     end subroutine f_python_execute_dict
  end interface

  interface
     subroutine f_python_execute(script, status)
       implicit none
       character(len = *), intent(in) :: script
       integer, intent(out) :: status
     end subroutine f_python_execute
  end interface

contains

  subroutine ndarray_new(ptr, obj)
    use f_precisions, only: f_loc
    implicit none
    type(ndarray), pointer :: ptr
    integer(f_address), intent(out) :: obj

    allocate(ptr)
    ptr = ndarray_null()
    obj = f_loc(ptr)
  end subroutine ndarray_new

  subroutine ndarray_free(ptr)
    implicit none
    type(ndarray), pointer :: ptr

    deallocate(ptr)
  end subroutine ndarray_free

  function ndarray_null() result(arr)
    type(ndarray) :: arr

    arr%data = int(0, f_address)
    arr%ndims = 0
  end function ndarray_null

  function i0_to_ndarray(data) result(arr)
    use f_precisions, only: f_loc
    implicit none
    integer, intent(in) :: data
    type(ndarray) :: arr
    
    arr%data = f_loc(data)
    arr%ndims = 0
    !arr%shapes = shape(data)
    arr%kind = "i4"
  end function i0_to_ndarray

  function i1_to_ndarray(data) result(arr)
    use f_precisions, only: f_loc
    implicit none
    integer, dimension(:), intent(in) :: data
    type(ndarray) :: arr
    
    arr%data = f_loc(data(1))
    arr%ndims = 1
    arr%shapes(1:1) = shape(data)
    arr%kind = "i4"
  end function i1_to_ndarray

  function i2_to_ndarray(data) result(arr)
    use f_precisions, only: f_loc
    implicit none
    integer, dimension(:, :), intent(in) :: data
    type(ndarray) :: arr
    
    arr%data = f_loc(data(1,1))
    arr%ndims = 2
    arr%shapes(1:2) = shape(data)
    arr%kind = "i4"
  end function i2_to_ndarray

  function f0_to_ndarray(data) result(arr)
    use f_precisions, only: f_loc
    implicit none
    double precision, intent(in) :: data
    type(ndarray) :: arr
    
    arr%data = f_loc(data)
    arr%ndims = 0
    !arr%shapes = shape(data)
    arr%kind = "f8"
  end function f0_to_ndarray

  function f1_to_ndarray(data) result(arr)
    use f_precisions, only: f_loc
    implicit none
    double precision, dimension(:), intent(in) :: data
    type(ndarray) :: arr
    
    arr%data = f_loc(data(1))
    arr%ndims = 1
    arr%shapes(1:1) = shape(data)
    arr%kind = "f8"
  end function f1_to_ndarray

  function f2_to_ndarray(data) result(arr)
    use f_precisions, only: f_loc
    implicit none
    double precision, dimension(:, :), intent(in) :: data
    type(ndarray) :: arr
    
    arr%data = f_loc(data(1,1))
    arr%ndims = 2
    arr%shapes(1:2) = shape(data)
    arr%kind = "f8"
  end function f2_to_ndarray


  function pi0_to_ndarray(data) result(arr)
    implicit none
    integer, pointer :: data
    type(ndarray) :: arr

    if (associated(data)) then
       arr = i0_to_ndarray(data)
    else
       arr = ndarray_null()
    end if
  end function pi0_to_ndarray

  function pi1_to_ndarray(data) result(arr)
    implicit none
    integer, dimension(:), pointer :: data
    type(ndarray) :: arr

    if (associated(data)) then
       arr = i1_to_ndarray(data)
    else
       arr = ndarray_null()
    end if
  end function pi1_to_ndarray

  function pi2_to_ndarray(data) result(arr)
    implicit none
    integer, dimension(:, :), pointer :: data
    type(ndarray) :: arr

    if (associated(data)) then
       arr = i2_to_ndarray(data)
    else
       arr = ndarray_null()
    end if
  end function pi2_to_ndarray

  function pf0_to_ndarray(data) result(arr)
    implicit none
    double precision, pointer :: data
    type(ndarray) :: arr

    if (associated(data)) then
       arr = f0_to_ndarray(data)
    else
       arr = ndarray_null()
    end if
  end function pf0_to_ndarray

  function pf1_to_ndarray(data) result(arr)
    implicit none
    double precision, dimension(:), pointer :: data
    type(ndarray) :: arr

    if (associated(data)) then
       arr = f1_to_ndarray(data)
    else
       arr = ndarray_null()
    end if
  end function pf1_to_ndarray

  function pf2_to_ndarray(data) result(arr)
    implicit none
    double precision, dimension(:, :), pointer :: data
    type(ndarray) :: arr

    if (associated(data)) then
       arr = f2_to_ndarray(data)
    else
       arr = ndarray_null()
    end if
  end function pf2_to_ndarray

end module f_python

subroutine f_python_ndarray_init()
  use f_python
  
  call f_object_add_method("class", "ndarray_new", ndarray_new, 1)
  call f_object_add_method("class", "ndarray_free", ndarray_free, 0)
end subroutine f_python_ndarray_init

subroutine f_python_ndarray_get(arr, data, ndims, shapes, kind)
  use f_python, only: ndarray
  use f_precisions, only: f_address
  implicit none
  type(ndarray), intent(in) :: arr
  integer(f_address), intent(out) :: data
  integer, intent(out) :: ndims
  integer, dimension(7), intent(out) :: shapes
  character(len = 2), intent(out) :: kind
  
  data = arr%data
  ndims = arr%ndims
  shapes = arr%shapes
  kind = arr%kind
end subroutine f_python_ndarray_get
