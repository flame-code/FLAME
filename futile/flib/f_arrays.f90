!> @file
!!  Introduce the concept of referenced arrays of different dimension
!! @author
!!    Copyright (C) 2013-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
module f_arrays
  use f_precisions
  use f_refcnts
  use dictionaries, only: f_err_throw
  use dynamic_memory
  
  private

  !> matrix to be used to store refcounted information
  type, public :: f_vector
     type(f_reference_counter) :: rc
     real(f_double), dimension(:), pointer :: ptr=>null()
  end type f_vector

  !> matrix to be used to store refcounted information
  type, public :: f_matrix
     type(f_reference_counter) :: rc
     real(f_double), dimension(:,:), pointer :: ptr=>null()
  end type f_matrix

  interface assignment(=)
     module procedure f_matrix_allocate,f_vector_allocate
     module procedure f_matrix_shallow_copy,f_vector_shallow_copy
     module procedure f_vector_copy,f_matrix_copy
     module procedure f_vector_allocate_1,f_matrix_allocate_3
  end interface assignment(=)

  interface f_array_deallocate
     module procedure f_vector_deallocate,f_matrix_deallocate
  end interface f_array_deallocate

  interface f_array_free
     module procedure f_vector_free,f_matrix_free
  end interface f_array_free

  interface f_array_ptr_free
     module procedure f_free_vector_1,f_free_matrix_ptr
  end interface f_array_ptr_free
  

  public :: assignment(=)
  public :: dump_f_matrix_ptr,f_array_ptr_free

  contains

    pure subroutine nullify_f_vector(arr)
      implicit none
      type(f_vector), intent(inout) :: arr
      call nullify_f_ref(arr%rc)
      nullify(arr%ptr)
    end subroutine nullify_f_vector
    pure subroutine nullify_f_matrix(arr)
      implicit none
      type(f_matrix), intent(inout) :: arr
      call nullify_f_ref(arr%rc)
      nullify(arr%ptr)
    end subroutine nullify_f_matrix

    subroutine f_vector_deallocate(mat)
      implicit none
      type(f_vector), intent(inout) :: mat
      call f_ref_free(mat%rc)
      call f_free_ptr(mat%ptr)
    end subroutine f_vector_deallocate
    subroutine f_matrix_deallocate(mat)
      implicit none
      type(f_matrix), intent(inout) :: mat
      call f_ref_free(mat%rc)
      call f_free_ptr(mat%ptr)
    end subroutine f_matrix_deallocate

    subroutine f_release_matrix(arr)
      implicit none
      type(f_matrix), intent(inout) :: arr
      !local variables
      integer :: count
      call f_unref(arr%rc,count=count)
      if (count == 0) call f_array_deallocate(arr)
    end subroutine f_release_matrix
    subroutine f_release_vector(arr)
      implicit none
      type(f_vector), intent(inout) :: arr
      !local variables
      integer :: count
      call f_unref(arr%rc,count=count)
      if (count == 0) call f_array_deallocate(arr)
    end subroutine f_release_vector

    subroutine f_vector_free(mat)
      implicit none
      type(f_vector), intent(inout) :: mat
      if (f_ref_count(mat%rc) > 0) call f_release_vector(mat)
      call nullify_f_vector(mat)
    end subroutine f_vector_free
    subroutine f_matrix_free(mat)
      implicit none
      type(f_matrix), intent(inout) :: mat
      if (f_ref_count(mat%rc) > 0) call f_release_matrix(mat)
      call nullify_f_matrix(mat)
    end subroutine f_matrix_free


    subroutine f_vector_shallow_copy(dest,src)
      implicit none
      type(f_vector), intent(inout) :: dest
      type(f_vector), intent(in) :: src

      call f_array_free(dest)
      !then perform the association
      dest%ptr=>src%ptr
      call f_ref_associate(dest=dest%rc,src=src%rc)

    end subroutine f_vector_shallow_copy
    subroutine f_matrix_shallow_copy(dest,src)
      implicit none
      type(f_matrix), intent(inout) :: dest
      type(f_matrix), intent(in) :: src

      call f_array_free(dest)
      !then perform the association
      dest%ptr=>src%ptr
      call f_ref_associate(dest=dest%rc,src=src%rc)

    end subroutine f_matrix_shallow_copy

    subroutine f_matrix_copy_fmat(dest,src)
      implicit none
      type(f_matrix), intent(inout) :: dest
      type(f_matrix), intent(in) :: src

      call f_array_free(dest)
      dest%ptr=f_malloc_ptr(src_ptr=src%ptr,id='ptr')
      dest%rc=f_ref_new('f_matrix')
    end subroutine f_matrix_copy_fmat

    !> overload the allocation rule
    subroutine f_vector_allocate(dest,m)
      implicit none
      type(malloc_information_ptr), intent(in) :: m
      type(f_vector), intent(inout) :: dest

      call f_array_free(dest)
      dest%ptr=m
      dest%rc=f_ref_new('f_vector')
    end subroutine f_vector_allocate
    subroutine f_matrix_allocate(dest,m)
      implicit none
      type(malloc_information_ptr), intent(in) :: m
      type(f_matrix), intent(inout) :: dest

      call f_array_free(dest)
      dest%ptr=m
      dest%rc=f_ref_new('f_matrix')
    end subroutine f_matrix_allocate

    !> copy an array in the f_vector
    subroutine f_vector_copy(dest,src)
      implicit none
      type(f_vector), intent(inout) :: dest
      real(f_double), dimension(:), intent(in) :: src

      dest=f_malloc_ptr(size(src),id='copyvec')
      call f_memcpy(src=src,dest=dest%ptr)

    end subroutine f_vector_copy
    subroutine f_matrix_copy(dest,src)
      implicit none
      type(f_matrix), intent(inout) :: dest
      real(f_double), dimension(:,:), intent(in) :: src

      dest=f_malloc_ptr([size(src,dim=1),size(src,dim=2)],id='copymat')
      call f_memcpy(src=src,dest=dest%ptr)

    end subroutine f_matrix_copy

    subroutine f_free_vector_1(arr)
      implicit none
      type(f_vector), dimension(:), pointer :: arr
      integer :: i

      if (.not. associated(arr)) return
      do i=lbound(arr,1),ubound(arr,1)
         call f_array_free(arr(i))
      end do
      deallocate(arr)
      nullify(arr)
    end subroutine f_free_vector_1

    subroutine f_free_matrix_ptr(arr)
      implicit none
      type(f_matrix), dimension(:,:,:), pointer :: arr
      !local variables
      integer :: i1,i2,i3
      if (.not. associated(arr)) return
      do i3=lbound(arr,3),ubound(arr,3)
         do i2=lbound(arr,2),ubound(arr,2)
            do i1=lbound(arr,1),ubound(arr,1)
               call f_array_free(arr(i1,i2,i3))
            end do
         end do
      end do
      deallocate(arr)
      nullify(arr)
    end subroutine f_free_matrix_ptr

    
    subroutine f_vector_allocate_1(array,m)
      use yaml_output
      use yaml_strings
      implicit none
      type(f_vector), dimension(:), pointer, intent(inout) :: array
      type(malloc_information_ptr), intent(in) :: m
      !local variables
      integer :: ierror

      allocate(array(m%lbounds(1):m%ubounds(1)),stat=ierror)

      if (ierror/=0) then
         call f_err_throw('array ' // trim(m%array_id) // &
              '(' // trim(yaml_toa(product(m%shape))) // &
              '), error code '//trim(yaml_toa(ierror)),&
              err_name='ERR_ALLOCATE')
         return
      end if
      if (size(shape(array)) /= m%rank) then
         call f_err_throw('Rank specified by f_malloc ('+&
              yaml_toa(m%rank)// ',routine_id=' // trim(m%routine_id) // &
              ',id=' // trim(m%array_id) //&
              ') is not coherent with the one of the array ('+&
              yaml_toa(size(shape(array)))//')',&
              err_name='ERR_INVALID_MALLOC')
              return
      end if

      !here the database for the allocation might be updated

    end subroutine f_vector_allocate_1

    subroutine f_matrix_allocate_3(array,m)
      use yaml_output
      use yaml_strings
      implicit none
      type(f_matrix), dimension(:,:,:), pointer, intent(inout) :: array
      type(malloc_information_ptr), intent(in) :: m
      !local variables
      integer :: ierror

      allocate(array(m%lbounds(1):m%ubounds(1),m%lbounds(2):m%ubounds(2),&
           m%lbounds(3):m%ubounds(3)),stat=ierror)

      if (ierror/=0) then
         call f_err_throw('array ' // trim(m%array_id) // &
              '(' // trim(yaml_toa(product(m%shape))) // &
              '), error code '//trim(yaml_toa(ierror)),&
              err_name='ERR_ALLOCATE')
         return
      end if
      if (size(shape(array)) /= m%rank) then
         call f_err_throw('Rank specified by f_malloc ('+&
              yaml_toa(m%rank)// ',routine_id=' // trim(m%routine_id) // &
              ',id=' // trim(m%array_id) //&
              ') is not coherent with the one of the array ('+&
              yaml_toa(size(shape(array)))//')',&
              err_name='ERR_INVALID_MALLOC')
              return
      end if

      !here the database for the allocation might be updated

    end subroutine f_matrix_allocate_3

    subroutine dump_f_matrix_ptr(key,arr)
      use yaml_output
      use yaml_strings
      character(len=*), intent(in) :: key
      type(f_matrix), dimension(:,:,:), pointer :: arr
      !local variables
      integer :: i1,i2,i3

      if (.not. associated(arr)) then
         call yaml_map(trim(key),'Nullified pointer')
         return
      end if
      do i3=lbound(arr,3),ubound(arr,3)
         do i2=lbound(arr,2),ubound(arr,2)
            do i1=lbound(arr,1),ubound(arr,1)
               if (associated(arr(i1,i2,i3)%ptr)) then
                  call yaml_map('"'+yaml_toa([i1,i2,i3])+'"',arr(i1,i2,i3)%ptr)
               else
                  call yaml_map('"'+yaml_toa([i1,i2,i3])+'"','Nullified pointer')
               end if
            end do
         end do
      end do
    end subroutine dump_f_matrix_ptr


end module f_arrays
