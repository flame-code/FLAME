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

  !> scalar type, to hnalde real or complex
  type, public :: f_scalar
     real(f_double) :: r = 0._f_double, i = 0._f_double
  end type f_scalar

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
     module procedure f_scalar_assign_r, f_scalar_assign_z, f_scalar_assign_r1
  end interface assignment(=)

  interface f_array_deallocate
     module procedure f_vector_deallocate,f_matrix_deallocate
  end interface f_array_deallocate

  interface f_array_free
     module procedure f_vector_free,f_matrix_free
     module procedure f_vector_v1_free,f_matrix_v1_free
  end interface f_array_free 

  interface f_array_ptr_free
     module procedure f_free_vector_1,f_free_matrix_ptr
  end interface f_array_ptr_free

  interface operator(+)
     module procedure f_scalar_add_sr, f_scalar_add_rs, f_scalar_add_zs, f_scalar_add_sz
     module procedure f_scalar_add_ss
  end interface operator(+)

  interface operator(*)
     module procedure f_scalar_prod_si, f_scalar_prod_is
     module procedure f_scalar_prod_sr, f_scalar_prod_rs, f_scalar_prod_zs, f_scalar_prod_sz
     module procedure f_scalar_prod_ss
  end interface operator(*)

  interface operator(-)
     module procedure f_scalar_opp
  end interface operator(-)

  interface abs
     module procedure f_scalar_abs
  end interface abs

  interface real
     module procedure f_scalar_real
  end interface real

  interface imag
     module procedure f_scalar_imag
  end interface imag

  public :: assignment(=), operator(+), operator(*), operator(-)
  public :: dump_f_matrix_ptr,f_array_ptr_free,f_array_free,abs,real,imag

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

    subroutine f_vector_v1_free(arr)
      implicit none
      type(f_vector), dimension(:), intent(inout) :: arr
      !local variables
      integer :: i
      do i=1,size(arr)
         call f_array_free(arr(i))
      end do
    end subroutine f_vector_v1_free

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
    subroutine f_matrix_v1_free(arr)
      implicit none
      type(f_matrix), dimension(:), intent(inout) :: arr
      !local variables
      integer :: i
      do i=1,size(arr)
         call f_array_free(arr(i))
      end do
    end subroutine f_matrix_v1_free

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
      !this section might be provided as an include file in the 
      !directory of installation to profile the allocation
      !of other data structures
      !local variables
      integer :: ierror

      call f_timer_interrupt(TCAT_ARRAY_ALLOCATIONS)

      allocate(array(m%lbounds(1):m%ubounds(1)),stat=ierror)

      !here the times is already resumed
      if (.not. malloc_validate(ierror,size(shape(array)),m)) return

      !here the database for the allocation might be updated

      call f_timer_resume()!TCAT_ARRAY_ALLOCATIONS
    end subroutine f_vector_allocate_1

    subroutine f_matrix_allocate_3(array,m)
      use yaml_output
      use yaml_strings
      implicit none
      type(f_matrix), dimension(:,:,:), pointer, intent(inout) :: array
      type(malloc_information_ptr), intent(in) :: m
      !local variables
      integer :: ierror

      call f_timer_interrupt(TCAT_ARRAY_ALLOCATIONS)
      
      allocate(array(m%lbounds(1):m%ubounds(1),m%lbounds(2):m%ubounds(2),&
           m%lbounds(3):m%ubounds(3)),stat=ierror)

      if (.not. malloc_validate(ierror,size(shape(array)),m)) return

      !here the database for the allocation might be updated
      call f_timer_resume()!TCAT_ARRAY_ALLOCATIONS
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

    pure elemental function cmplx_imag(z) result(i)
      complex(f_double), intent(in) :: z
      real(f_double) :: i

      complex(f_double) :: zintern
      real(f_double), dimension(2) :: intern

      equivalence(zintern, intern)
      zintern = z
      i = intern(2)
    end function cmplx_imag

    pure elemental function f_scalar_add_sr(sa, b) result(s)
      implicit none
      type(f_scalar), intent(in) :: sa
      real(f_double), intent(in) :: b
      type(f_scalar) :: s

      s = f_scalar(sa%r + b, sa%i)
    end function f_scalar_add_sr

    pure elemental function f_scalar_add_rs(b, sa) result(s)
      implicit none
      type(f_scalar), intent(in) :: sa
      real(f_double), intent(in) :: b
      type(f_scalar) :: s

      s = f_scalar(sa%r + b, sa%i)
    end function f_scalar_add_rs

    pure elemental function f_scalar_add_sz(sa, b) result(s)
      implicit none
      type(f_scalar), intent(in) :: sa
      complex(f_double), intent(in) :: b
      type(f_scalar) :: s

      s = f_scalar(sa%r + real(b, f_double), sa%i + cmplx_imag(b))
    end function f_scalar_add_sz

    pure elemental function f_scalar_add_zs(b, sa) result(s)
      implicit none
      type(f_scalar), intent(in) :: sa
      complex(f_double), intent(in) :: b
      type(f_scalar) :: s

      s = f_scalar(sa%r + real(b, f_double), sa%i + cmplx_imag(b))
    end function f_scalar_add_zs

    pure elemental function f_scalar_add_ss(sa, sb) result(s)
      implicit none
      type(f_scalar), intent(in) :: sa, sb
      type(f_scalar) :: s

      s = f_scalar(sa%r + sb%r, sa%i + sb%i)
    end function f_scalar_add_ss

    pure elemental function f_scalar_prod_si(sa, b) result(s)
      implicit none
      type(f_scalar), intent(in) :: sa
      integer, intent(in) :: b
      type(f_scalar) :: s

      s = f_scalar(sa%r * b, sa%i * b)
    end function f_scalar_prod_si

    pure elemental function f_scalar_prod_is(b, sa) result(s)
      implicit none
      type(f_scalar), intent(in) :: sa
      integer, intent(in) :: b
      type(f_scalar) :: s

      s = f_scalar(sa%r * b, sa%i * b)
    end function f_scalar_prod_is

    pure elemental function f_scalar_prod_sr(sa, b) result(s)
      implicit none
      type(f_scalar), intent(in) :: sa
      real(f_double), intent(in) :: b
      type(f_scalar) :: s

      s = f_scalar(sa%r * b, sa%i * b)
    end function f_scalar_prod_sr

    pure elemental function f_scalar_prod_rs(b, sa) result(s)
      implicit none
      type(f_scalar), intent(in) :: sa
      real(f_double), intent(in) :: b
      type(f_scalar) :: s

      s = f_scalar(sa%r * b, sa%i * b)
    end function f_scalar_prod_rs

    pure elemental function f_scalar_prod_sz(sa, b) result(s)
      implicit none
      type(f_scalar), intent(in) :: sa
      complex(f_double), intent(in) :: b
      type(f_scalar) :: s

      s = f_scalar(sa%r * real(b, f_double) - sa%i * cmplx_imag(b), &
           & sa%i * real(b, f_double) + sa%r * cmplx_imag(b))
    end function f_scalar_prod_sz

    pure elemental function f_scalar_prod_zs(b, sa) result(s)
      implicit none
      type(f_scalar), intent(in) :: sa
      complex(f_double), intent(in) :: b
      type(f_scalar) :: s

      s = f_scalar(sa%r * real(b, f_double) - sa%i * cmplx_imag(b), &
           & sa%i * real(b, f_double) + sa%r * cmplx_imag(b))
    end function f_scalar_prod_zs

    pure elemental function f_scalar_prod_ss(sa, sb) result(s)
      implicit none
      type(f_scalar), intent(in) :: sa, sb
      type(f_scalar) :: s

      s = f_scalar(sa%r * sb%r - sa%i * sb%i, sa%i * sb%r + sa%r * sb%i)
    end function f_scalar_prod_ss

    pure elemental function f_scalar_opp(sa) result(s)
      implicit none
      type(f_scalar), intent(in) :: sa
      type(f_scalar) :: s

      s = f_scalar(-sa%r, -sa%i)
    end function f_scalar_opp

    pure elemental subroutine f_scalar_assign_r(s, r)
      implicit none
      type(f_scalar), intent(out) :: s
      real(f_double), intent(in) :: r

      s = f_scalar(r, 0._f_double)
    end subroutine f_scalar_assign_r

    pure subroutine f_scalar_assign_r1(s, r)
      implicit none
      type(f_scalar), intent(out) :: s
      real(f_double), dimension(:), intent(in) :: r

      s = f_scalar(r(1), 0._f_double)
      if (size(r) > 1) s%i = r(2)
    end subroutine f_scalar_assign_r1

    pure elemental subroutine f_scalar_assign_z(s, z)
      implicit none
      type(f_scalar), intent(out) :: s
      complex(f_double), intent(in) :: z

      s = f_scalar(real(z, f_double), cmplx_imag(z))
    end subroutine f_scalar_assign_z

    pure elemental function f_scalar_abs(s) result(r)
      implicit none
      type(f_scalar), intent(in) :: s
      real(f_double) :: r

      r = sqrt(s%r ** 2 + s%i **2)
    end function f_scalar_abs

    pure elemental function f_scalar_real(s) result(r)
      implicit none
      type(f_scalar), intent(in) :: s
      real(f_double) :: r

      r = s%r
    end function f_scalar_real

    pure elemental function f_scalar_imag(s) result(r)
      implicit none
      type(f_scalar), intent(in) :: s
      real(f_double) :: r

      r = s%i
    end function f_scalar_imag

end module f_arrays
