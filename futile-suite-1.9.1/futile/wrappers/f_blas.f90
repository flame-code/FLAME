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
  use wrapper_linalg
  implicit none

  private

  type, public :: f_eye
     integer :: n !< size of the identity matrix
  end type f_eye

  interface f_gemv
     module procedure f_gemv_md0,f_gemv_id0
  end interface f_gemv

  interface f_axpy
     module procedure f_axpy_d0,f_axpy_d1
  end interface f_axpy

  interface f_dot
     module procedure f_dot_d1,f_dot_d3,f_dot_d4d3
  end interface f_dot


  public :: f_gemv,f_axpy,f_dot

  contains

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
      use dictionaries, only: f_err_throw
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

    function f_dot_d1(x,y,stride_x,stride_y) result(c)
      use dictionaries, only: f_err_throw
      implicit none
      real(f_double), dimension(:), intent(in) :: x,y
      real(f_double) :: c
      integer, intent(in), optional :: stride_x,stride_y
      !local variables
      integer :: n,incx,incy

      n=size(x)
      if (n/=size(y)) then
         call f_err_throw('Error in dot, the size of the x and y array do not coincide')!,&
      !err_name=ERROR_LINALG) to be defined for BLAS and LAPACK routines
         return
      end if

      if (n==0) return
      incx=1
      if (present(stride_x)) incx=stride_x
      incy=1
      if (present(stride_y)) incy=stride_y

      !here we also have to insert the broker for the GPU acceleration
      c=dot(n,x(1),incx,y(1),incy)
      
    end function f_dot_d1

    function f_dot_d3(x,y,stride_x,stride_y) result(c)
      use dictionaries, only: f_err_throw
      implicit none
      real(f_double), dimension(:,:,:), intent(in) :: x,y
      real(f_double) :: c
      integer, intent(in), optional :: stride_x,stride_y
      !local variables
      integer :: n,incx,incy

      n=size(x)
      if (n/=size(y)) then
         call f_err_throw('Error in dot, the size of the x and y array do not coincide')!,&
         !err_name=ERROR_LINALG) to be defined for BLAS and LAPACK routines
         return
      end if

      if (n==0) return
      incx=1
      if (present(stride_x)) incx=stride_x
      incy=1
      if (present(stride_y)) incy=stride_y

      !here we also have to insert the broker for the GPU acceleration
      c=dot(n,x(1,1,1),incx,y(1,1,1),incy)

    end function f_dot_d3

    function f_dot_d4d3(x,y,stride_x,stride_y) result(c)
      use dictionaries, only: f_err_throw
      implicit none
      real(f_double), dimension(:,:,:,:), intent(in) :: x
      real(f_double), dimension(:,:,:), intent(in) :: y
      real(f_double) :: c
      integer, intent(in), optional :: stride_x,stride_y
      !local variables
      integer :: n,incx,incy

      n=size(x)
      if (n/=size(y)) then
         call f_err_throw('Error in dot, the size of the x and y array do not coincide')!,&
         !err_name=ERROR_LINALG) to be defined for BLAS and LAPACK routines
         return
      end if

      if (n==0) return
      incx=1
      if (present(stride_x)) incx=stride_x
      incy=1
      if (present(stride_y)) incy=stride_y

      !here we also have to insert the broker for the GPU acceleration
      c=dot(n,x(1,1,1,1),incx,y(1,1,1),incy)

    end function f_dot_d4d3

    

end module f_blas
