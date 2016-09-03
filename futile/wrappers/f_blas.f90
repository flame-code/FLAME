!> @file
!!  Define the operation which are associated to blas
!! @author
!!    Copyright (C) 2013-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
module f_blas
  use f_precisions
  use f_refcnts
  use wrapper_linalg
  use dictionaries, only: f_err_throw
  use dynamic_memory
  implicit none

  type, public :: f_eye
     integer :: n !< size of the identity matrix
  end type f_eye

  !> matrix to be used to store refcounted information
  type, public :: f_matrix
     type(f_reference_counter) :: rc
     real(f_double), dimension(:,:), pointer :: dmat
  end type f_matrix

  interface f_gemv
     module procedure f_gemv_md0,f_gemv_id0
  end interface f_gemv

  interface assignment(=)
     module procedure f_matrix_shallow_copy
  end interface assignment(=)

  public :: f_gemv
  public :: f_free_matrix_ptr,f_matrix_allocate,assignment(=),f_matrix_allocate_ptr
  public :: dump_f_matrix_ptr

  contains

    pure subroutine nullify_f_matrix(mat)
      implicit none
      type(f_matrix), intent(out) :: mat
      call nullify_f_ref(mat%rc)
      nullify(mat%dmat)
    end subroutine nullify_f_matrix
    pure function f_matrix_null() result(mat)
      type(f_matrix) :: mat
      call nullify_f_matrix(mat)
    end function f_matrix_null

    subroutine f_release_matrix(mat)
      implicit none
      type(f_matrix), intent(inout) :: mat
      !local variables
      integer :: count
      call f_unref(mat%rc,count=count)
      if (count == 0) then
         call f_ref_free(mat%rc)
         call f_free_ptr(mat%dmat)
      end if
    end subroutine f_release_matrix
   
    subroutine f_free_matrix(mat)
      implicit none
      type(f_matrix), intent(inout) :: mat
      if (f_ref_count(mat%rc) > 0) call f_release_matrix(mat)
      call nullify_f_matrix(mat)
    end subroutine f_free_matrix

    subroutine f_free_matrix_ptr(arr)
      implicit none
      type(f_matrix), dimension(:,:,:), pointer :: arr
      !local variables
      integer :: i1,i2,i3
      if (.not. associated(arr)) return
      do i3=lbound(arr,3),ubound(arr,3)
         do i2=lbound(arr,2),ubound(arr,2)
            do i1=lbound(arr,1),ubound(arr,1)
               call f_free_matrix(arr(i1,i2,i3))
            end do
         end do
      end do
      deallocate(arr)
      nullify(arr)
    end subroutine f_free_matrix_ptr

    subroutine f_matrix_shallow_copy(dest,src)
      implicit none
      type(f_matrix), intent(inout) :: dest
      type(f_matrix), intent(in) :: src
      
      call f_free_matrix(dest)
     
      !then perform the association
      dest%dmat=>src%dmat
      call f_ref_associate(dest=dest%rc,src=src%rc)

    end subroutine f_matrix_shallow_copy

    subroutine f_matrix_copy_fmat(dest,src)
      implicit none
      type(f_matrix), intent(inout) :: dest
      type(f_matrix), intent(in) :: src

      call f_free_matrix(dest)
      dest%dmat=f_malloc_ptr(src_ptr=src%dmat,id='dmat')
      dest%rc=f_ref_new('f_matrix')
    end subroutine f_matrix_copy_fmat

    subroutine f_matrix_allocate(dest,n)
      implicit none
      type(f_matrix), intent(inout) :: dest
      integer, intent(in) :: n

      call f_free_matrix(dest)
      dest%dmat=f_malloc_ptr([n,n],id='dmat')
      dest%rc=f_ref_new('f_matrix')
    end subroutine f_matrix_allocate
   
    subroutine f_matrix_allocate_ptr(arr,bounds)
      implicit none
      type(array_bounds), dimension(3), intent(in) :: bounds
      type(f_matrix), dimension(:,:,:), pointer :: arr
      !local variables
      integer :: i1,i2,i3

      allocate(arr(bounds(1)%nlow:bounds(1)%nhigh,&
           bounds(2)%nlow:bounds(2)%nhigh,bounds(3)%nlow:bounds(3)%nhigh))
      do i3=lbound(arr,3),ubound(arr,3)
         do i2=lbound(arr,2),ubound(arr,2)
            do i1=lbound(arr,1),ubound(arr,1)
               call nullify_f_matrix(arr(i1,i2,i3))
            end do
         end do
      end do
    end subroutine f_matrix_allocate_ptr

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
               if (associated(arr(i1,i2,i3)%dmat)) then
                  call yaml_map('"'+yaml_toa([i1,i2,i3])+'"',arr(i1,i2,i3)%dmat)
               else
                  call yaml_map('"'+yaml_toa([i1,i2,i3])+'"','Nullified pointer')
               end if
            end do
         end do
      end do
    end subroutine dump_f_matrix_ptr

    subroutine f_gemv_md0(alpha,A,x,beta,y)
      implicit none
      real(f_double), intent(in) :: alpha,beta
      real(f_double), dimension(:,:), intent(in) :: A
      real(f_double) :: x,y
      !local variables
      integer,dimension(2) :: nm
      
      nm=shape(A)
      if (alpha==0.0_f_double .or. nm(2) ==0) return
      call dgemv('N',nm(1),nm(2),alpha,A,nm(1),x,1,beta,y,1)

    end subroutine f_gemv_md0


    subroutine f_gemv_id0(alpha,A,x,beta,y)
      implicit none
      real(f_double), intent(in) :: alpha,beta
      type(f_eye), intent(in) :: A
      real(f_double) :: x,y

      if (beta==0.0_f_double) then
         stop 'to be implemented'
      else if (beta==1.0_f_double) then
         call f_axpy_d0(A%n,alpha,x,y)         
      else
         stop 'again to be implemented'
      end if
    end subroutine f_gemv_id0

    subroutine f_axpy_d1(a,x,y,stride_x,stride_y)
      implicit none
      real(f_double), intent(in) :: a
      real(f_double), dimension(:), intent(in) :: x
      real(f_double), dimension(:), intent(inout) :: y
      integer, intent(in), optional :: stride_x,stride_y
      !local variables
      integer :: n,incx,incy

      n=size(x)
      if (n/=size(y)) &
           call f_err_throw('Error in axpy, the size of the x and y array do not coincide')!,&
      !err_name=ERROR_LINALG) to be defined for BLAS and LAPACK routines
      
      if (n==0) return
      if (a==0.0_f_double) return
      incx=1
      if (present(stride_x)) incx=stride_x
      incy=1
      if (present(stride_y)) incy=stride_y

      !here we also have to insert the broker for the GPU acceleration
      call axpy(n,a,x(1),incx,y(1),incy)

    end subroutine f_axpy_d1

    subroutine f_axpy_d0(n,a,x,y,stride_x,stride_y)
      implicit none
      integer, intent(in) :: n
      real(f_double), intent(in) :: a
      real(f_double) :: x
      real(f_double) :: y
      integer, intent(in), optional :: stride_x,stride_y
      !local variables
      integer :: incx,incy

      if (n==0) return
      incx=1
      if (present(stride_x)) incx=stride_x
      incy=1
      if (present(stride_y)) incy=stride_y

      !here we also have to insert the broker for the GPU acceleration
      call axpy(n,a,x,incx,y,incy)

    end subroutine f_axpy_d0

end module f_blas
