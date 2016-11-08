!> @file
!! Define the concept of memory buffer, useful for DMA and device memory handling
!! @author
!!    Copyright (C) 2014-2015 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
module f_buffers
  use iso_c_binding
  use dictionaries

  implicit none

  character(len=*), parameter :: DATATYPE_KEY='datatype'
  character(len=*), parameter :: METHOD_KEY='method'

  character(len=*), parameter :: F_PTR_UC='F_PTR'
  character(len=*), parameter :: F_PTR_LC='f_ptr'
  character(len=*), parameter :: C_PTR_UC='C_PTR'
  character(len=*), parameter :: C_PTR_LC='c_ptr'
  character(len=*), parameter :: CUDA_PTR_UC='CUDA_PTR'
  character(len=*), parameter :: CUDA_PTR_LC='cuda_ptr'
  character(len=*), parameter :: OCL_PTR_UC='OCL_PTR'
  character(len=*), parameter :: OCL_PTR_LC='ocl_ptr'
  character(len=*), parameter :: DATATYPE_INT_LC='integer'
  character(len=*), parameter :: DATATYPE_INT_UC='INTEGER'
  character(len=*), parameter :: DATATYPE_DOUBLE_LC='double'
  character(len=*), parameter :: DATATYPE_DOUBLE_UC='DOUBLE'


  integer, parameter :: HOST=0
  integer, parameter :: CUDA=1
  integer, parameter :: OCL=2

  integer, parameter :: GPU_INTEGER=1
  integer, parameter :: GPU_SIMPLE=-1
  integer, parameter :: GPU_DOUBLE=2


  type, public :: f_buffer
     !> information about the buffer, extensible metadata
     type(dictionary), pointer :: info
     !intrinsic pointers
     integer(f_integer), dimension(:), pointer :: f_ptr_i
     real(f_double), dimension(:), pointer :: f_ptr_d
     type(c_ptr) :: ptr !< points to GPU data of C_PTR
  end type f_buffer

contains

  pure subroutine nullify_f_buffer(g)
    implicit none
    type(f_buffer), intent(out) :: g
    nullify(g%info)
    nullify(g%f_ptr)
    g%ptr=C_NULL_PTR
  end subroutine nullify_GPU_ptr

  pure function sizeof(g) result(b)
    implicit none
    type(f_buffer), intent(in) :: g
    !local variables
    character(len=20) :: datatype
    
    call f_zero(datatype)  
    datatype=g%info .get. DATATYPE_KEY
    select case(trim(g%datatype))
    case(DATATYPE_INT_LC,DATATYPE_INT_UC)
       b=4
    case(DATATYPE_DOUBLE_LC,DATATYPE_DOUBLE_UC)
       b=8
    case default
       b=0
    end select
    
  end function sizeof

  !>low level allocation of the buffer pointers
  subroutine allocate_f_buf(buf,size,datatype,method,ierror)
    implicit none
    type(f_buffer), intent(in) :: buf
    integer, intent(in) :: size
    character(len=*), intent(in) :: datatype,method
    integer, intent(out) :: ierror
    integer(f_long) :: size_t

    ierror=-1
    !initialize the buffer metadata
    if (associated(buf%info)) &
         call f_err_throw('Attempting to reallocate buffer, memory leak')
    size_t=size*sizeof(buf)
    select case(trim(method))
    case(F_PTR_UC,F_PTR_LC)
       select case(trim(datatype))
       case(DATATYPE_INT_LC,DATATYPE_INT_UC)
          allocate(buf%f_ptr_i(size),stat=ierror)
       case(DATATYPE_DOUBLE_LC,DATATYPE_DOUBLE_UC)
          allocate(buf%f_ptr_d(size),stat=ierror)
       case default
          call err()
       end select
    case(C_PTR_UC,C_PTR_LC)
       select case(trim(datatype))
          !allocation for c_pointer, might be wrapped in C
       case default
          call err()
       end select
    case(CUDA_PTR_UC,CUDA_PTR_LC)
       select case(trim(datatype))
       case(DATATYPE_DOUBLE_LC,DATATYPE_DOUBLE_UC)
          call cudamalloc(size,buf%ptr,ierror)
       case default
          call err()
       end select
    case(OCL_PTR_UC,OCL_PTR_LC)
       select case(trim(datatype))
       !case(DATATYPE_DOUBLE_LC,DATATYPE_DOUBLE_UC)
       case default
          call err()
       end select
    case default
       call f_err_throw('Method type "'//trim(method)//&
            '"unknown for f_buffer')
    end select 

   buf%info => dict_new(DATATYPE_KEY .is. datatype,&
        METHOD_KEY .is. method,&
        SIZEOF_KEY .is. size_t)
    
   contains
     subroutine err()
          call f_err_throw('Datatype '//trim(datatype)//&
               'not supported for method '//trim(method))
     end subroutine err

  end subroutine allocate_f_buf

  !>low level deallocation of the buffer pointers
  subroutine deallocate_f_buf(buf,ierror)
    implicit none
    type(f_buffer), intent(in) :: buf
    integer, intent(out) :: ierror
    !local variables
    character(len=32) :: datatype,method

    ierror=-1
    call f_zero(datatype)
    call f_zero(method)

    !initialize the buffer metadata
    if (.not. associated(buf%info)) &
         call f_err_throw('Attempting to deallocate free buffer, double free')
    datatype=buf%info .get. DATATYPE_KEY
    method=buf%info .get. METHOD_KEY
    select case(trim(method))
    case(F_PTR_UC,F_PTR_LC)
       select case(trim(datatype))
       case(DATATYPE_INT_LC,DATATYPE_INT_UC)
          deallocate(buf%f_ptr_i,stat=ierror)
       case(DATATYPE_DOUBLE_LC,DATATYPE_DOUBLE_UC)
          deallocate(buf%f_ptr_d,stat=ierror)
       case default
          call err()
       end select
    case(C_PTR_UC,C_PTR_LC)
       select case(trim(datatype))
          !deallocation for c_pointer, might be wrapped in C
       case default
          call err()
       end select
    case(CUDA_PTR_UC,CUDA_PTR_LC)
       select case(trim(datatype))
       case(DATATYPE_DOUBLE_LC,DATATYPE_DOUBLE_UC)
          call cudafree(buf%ptr)
          ierror=0 !we should wrap the cuda errors
       case default
          call err()
       end select
    case(OCL_PTR_UC,OCL_PTR_LC)
       select case(trim(datatype))
          !case(DATATYPE_DOUBLE_LC,DATATYPE_DOUBLE_UC)
       case default
          call err()
       end select
    case default
       call f_err_throw('Method type "'//trim(method)//&
            '"unknown for f_buffer')
    end select

    if (ierror==0) then
       call dict_free(buf%info)
       call nullify_f_buffer(buf)
    end if

  contains
    subroutine err()
      call f_err_throw('Datatype '//trim(datatype)//&
           'not supported for method '//trim(method))
    end subroutine err

  end subroutine deallocate_f_buf

  subroutine zero_f_buf(buf)
    implicit none
    type(f_buffer), intent(in) :: buf
    !local variables
    integer :: ierror,size
    integer(f_long) :: size_t
    character(len=32) :: datatype,method

    call f_zero(datatype)
    call f_zero(method)

    !initialize the buffer metadata
    if (.not. associated(buf%info)) &
         call f_err_throw('Attempting to deallocate free buffer, double free')
    datatype=buf%info .get. DATATYPE_KEY
    method=buf%info .get. METHOD_KEY
    select case(trim(method))
    case(F_PTR_UC,F_PTR_LC)
       select case(trim(datatype))
       case(DATATYPE_INT_LC,DATATYPE_INT_UC)
          call f_zero(buf%f_ptr_i)
       case(DATATYPE_DOUBLE_LC,DATATYPE_DOUBLE_UC)
          call f_zero(buf%f_ptr_d)
       case default
          call err()
       end select
    case(C_PTR_UC,C_PTR_LC)
       select case(trim(datatype))
          !deallocation for c_pointer, might be wrapped in C
       case default
          call err()
       end select
    case(CUDA_PTR_UC,CUDA_PTR_LC)
       select case(trim(datatype))
       case(DATATYPE_DOUBLE_LC,DATATYPE_DOUBLE_UC)
          !memset to be simplified
          size_t=buf%info//SIZEOF_KEY
          size=int(size_t/sizeof(buf))
          call cudamemset(buf%ptr, 0, size,ierror)
          if (ierror /=0) call f_err_throw('Error in zero_f_buf')
       case default
          call err()
       end select
    case(OCL_PTR_UC,OCL_PTR_LC)
       select case(trim(datatype))
          !case(DATATYPE_DOUBLE_LC,DATATYPE_DOUBLE_UC)
       case default
          call err()
       end select
    case default
       call f_err_throw('Method type "'//trim(method)//&
            '" unknown for f_buffer')
    end select

  contains
    subroutine err()
      call f_err_throw('Datatype '//trim(datatype)//&
           'not supported for method '//trim(method))
    end subroutine err

  end subroutine zero_f_buf

  subroutine deppcopy_f_buf(buf,src_add)
    implicit none
    type(f_buffer), intent(in) :: buf
    integer(f_long), intent(in) :: src_add
    !local variables
    integer :: ierror,size
    integer(f_long) :: size_t
    character(len=32) :: datatype,method

    call f_zero(datatype)
    call f_zero(method)

    !initialize the buffer metadata
    if (.not. associated(buf%info)) &
         call f_err_throw('Attempting to copy into free buffer')
    datatype=buf%info .get. DATATYPE_KEY
    method=buf%info .get. METHOD_KEY
    select case(trim(method))
    case(F_PTR_UC,F_PTR_LC)
       select case(trim(datatype))
       case default
          call err()
       end select
    case(C_PTR_UC,C_PTR_LC)
       select case(trim(datatype))
          !deallocation for c_pointer, might be wrapped in C
       case default
          call err()
       end select
    case(CUDA_PTR_UC,CUDA_PTR_LC)
       select case(trim(datatype))
       case(DATATYPE_DOUBLE_LC,DATATYPE_DOUBLE_UC)
          !memset to be simplified
          size_t=buf%info//SIZEOF_KEY
          call copy_buffer_to_gpu_fromadd(array%sizeof,m%srcdata_add,array%ptr,ierror)
          if (ierror /=0) call f_err_throw('Error in deepcopy_f_buf')
       case default
          call err()
       end select
    case(OCL_PTR_UC,OCL_PTR_LC)
       select case(trim(datatype))
          !case(DATATYPE_DOUBLE_LC,DATATYPE_DOUBLE_UC)
       case default
          call err()
       end select
    case default
       call f_err_throw('Method type "'//trim(method)//&
            '" unknown for f_buffer')
    end select

  contains
    subroutine err()
      call f_err_throw('Datatype '//trim(datatype)//&
           'not supported for method '//trim(method))
    end subroutine err

  end subroutine deppcopy_f_buf


end module f_buffers
