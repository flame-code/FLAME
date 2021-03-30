!> @file
!! Include fortran file for allocation templates
!! file included in module dynamic_memory.f90
!! @author
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
!subroutine xx_all(array,m)
!  use metadata_interfaces
!  implicit none
!  type(malloc_information_all), intent(in) :: m
!  integer, dimension(:), allocatable, intent(inout) :: array
!  !--- allocate_profile-inc.f90
!  integer :: ierror
!  integer(kind=8) :: iadd
!  character(len=info_length) :: val
!  type(dictionary), pointer :: info
!  logical :: c_allocation
!  integer :: padding
!  !$ logical :: not_omp
!  !$ logical, external :: omp_in_parallel,omp_get_nested
!
!  if (f_err_raise(ictrl == 0,&
!       'ERROR (f_malloc): the routine f_malloc_initialize has not been called',&
!       ERR_MALLOC_INTERNAL)) return
!
!  !$ not_omp=.not. (omp_in_parallel() .or. omp_get_nested())
!
!  padding=ndebug !which is always zero so far
!
!  !here we should add a control of the OMP behaviour of allocation
!  !in particular for what concerns the OMP nesting procedure
!  !the following action is the allocation
!  !$ if(not_omp) then
!  call f_timer_interrupt(TCAT_ARRAY_ALLOCATIONS)
!  !$ end if
!  if (f_nan_pad_size > 0) then
!     padding=f_nan_pad_size
!     call togglepadding(product(int(shape(1:m%rank-1),f_long))*&
!          kind(array)*padding)
!  end if
!  !END--- allocate_profile-inc.f90
!  alignment=0
!  if(len_trim(m%info) > 0) then
!     info => yaml_load(m%info)
!     alignment=info//'Alignment'
!     call dict_free(info)
!  end if
!  if (alignment /=0)
!    array=aligned_alloc(align_size,size)
!  else

!  !allocate the array
!  allocate(array(m%lbounds(1):m%ubounds(1)+ndebug),stat=ierror)
!end
!  !--- allocate-inc.f90
!  if (f_nan_pad_size > 0) call togglepadding(0)
!  if (ierror/=0) then
!     !$ if(not_omp) then
!     call f_timer_resume()!TCAT_ARRAY_ALLOCATIONS
!     !$ end if
!     call f_err_throw('Allocation problem, error code '//trim(yaml_toa(ierror)),ERR_ALLOCATE)
!     return
!  end if
!  if (size(shape(array))==m%rank) then
!     call pad_array(array,m%put_to_zero,m%shape,ndebug)
!     !also fill the array with the values of the source if the address is identified in the source
!     if (m%srcdata_add > int(0,kind=8)) call c_memcopy(array,m%srcdata_add,f_sizeof(array))
!     !profile the array allocation
!     iadd=int(0,kind=8)
!        !write the address of the first element in the address string
!     if (m%profile .and. track_origins) iadd=loc_arr(array)!call getlongaddress(array,iadd)
!
!     call f_update_database(product(int(m%shape(1:m%rank),kind=8)),kind(array),m%rank,&
!          iadd,m%array_id,m%routine_id)
!
!  else
!     !$ if(not_omp) then
!     call f_timer_resume()!TCAT_ARRAY_ALLOCATIONS
!     !$ end if
!     call f_err_throw('Rank specified by f_malloc ('//trim(yaml_toa(m%rank))//&
!          ') is not coherent with the one of the array ('//trim(yaml_toa(size(shape(array))))//')',&
!          ERR_INVALID_MALLOC)
!     return
!  end if
!  !$ if(not_omp) then
!  call f_timer_resume()!TCAT_ARRAY_ALLOCATIONS
!  !$ end if
!  !END--- allocate-inc.f90
!end subroutine xx_all
!
!
!subroutine xx_all_free(array)
!  use metadata_interfaces
!  implicit none
!  integer, dimension(:), allocatable, intent(inout) :: array
!  !--'deallocate-profile-inc.f90'
!  !local variables
!  integer :: ierror
!  !$ logical :: not_omp
!  !$ logical, external :: omp_in_parallel,omp_get_nested
!  integer(kind=8) :: ilsize,iadd
!
!  if (f_err_raise(ictrl == 0,&
!       'ERROR (f_free): the routine f_malloc_initialize has not been called',&
!       ERR_MALLOC_INTERNAL)) return
!
!  !$ not_omp=.not. (omp_in_parallel() .or. omp_get_nested())
!
!  !here we should add a control of the OMP behaviour of allocation
!  !in particular for what concerns the OMP nesting procedure
!
!  !END--'deallocate-profile-inc.f90'
!  !-- 'deallocate-inc.f90'
!  !$ if (not_omp) then
!  call f_timer_interrupt(TCAT_ARRAY_ALLOCATIONS)
!  !$ end if
!
!  !here the size should be corrected with ndebug (or maybe not)
!  ilsize=product(int(shape(array),kind=8))
!  !retrieve the address of the first element if the size is not zero
!  !iadd=int(0,kind=8)
!  !if (ilsize /= int(0,kind=8))
!  iadd=loc_arr(array)!call getlongaddress(array,iadd)
!  !fortran deallocation
!  deallocate(array,stat=ierror)
!
!  if (ierror/=0) then
!     !$ if (not_omp) then
!     call f_timer_resume()!TCAT_ARRAY_ALLOCATIONS
!     !$ end if
!     call f_err_throw('Deallocation problem, error code '//trim(yaml_toa(ierror)),&
!          ERR_DEALLOCATE)
!     return
!  end if
!
!  call f_purge_database(ilsize,kind(array),iadd)
!
!  !$ if (not_omp) then
!  call f_timer_resume()!TCAT_ARRAY_ALLOCATIONS
!  !$ end if
!  !END-- 'deallocate-inc.f90'
!end subroutine xx_all_free


subroutine i1_all(array,m)
  use metadata_interfaces, metadata_address => geti1
  implicit none
  type(malloc_information_all), intent(in) :: m
  integer(kind=4), dimension(:), allocatable, intent(inout) :: array
  include 'allocate-profile-inc.f90'
  !allocate the array
  allocate(array(m%lbounds(1):m%ubounds(1)+ndebug),stat=ierror)
  include 'allocate-inc.f90'
end subroutine i1_all

subroutine i1_all_free(array)
  use metadata_interfaces, metadata_address => geti1
  implicit none
  integer(kind=4), dimension(:), allocatable, intent(inout) :: array
  include 'deallocate-profile-inc.f90'
  include 'deallocate-inc.f90'
end subroutine i1_all_free

subroutine i2_all(array,m)
  use metadata_interfaces, metadata_address => geti2
  implicit none
  type(malloc_information_all), intent(in) :: m
  integer(kind=4), dimension(:,:), allocatable, intent(inout) :: array
  include 'allocate-profile-inc.f90'
  !allocate the array
  allocate(array(m%lbounds(1):m%ubounds(1),m%lbounds(2):m%ubounds(2)+ndebug),stat=ierror)
  include 'allocate-inc.f90'
end subroutine i2_all

subroutine i2_all_free(array)
  use metadata_interfaces, metadata_address => geti2
  implicit none
  integer(kind=4), dimension(:,:), allocatable, intent(inout) :: array
  include 'deallocate-profile-inc.f90'
  include 'deallocate-inc.f90'
end subroutine i2_all_free

subroutine i3_all(array,m)
  use metadata_interfaces, metadata_address => geti3
  implicit none
  type(malloc_information_all), intent(in) :: m
  integer(kind=4), dimension(:,:,:), allocatable, intent(inout) :: array
  include 'allocate-profile-inc.f90'
  !allocate the array
  allocate(array(m%lbounds(1):m%ubounds(1),m%lbounds(2):m%ubounds(2),m%lbounds(3):m%ubounds(3)+ndebug),stat=ierror)
  include 'allocate-inc.f90'
end subroutine i3_all

subroutine i3_all_free(array)
  use metadata_interfaces, metadata_address => geti3
  implicit none
  integer(kind=4), dimension(:,:,:), allocatable, intent(inout) :: array
  include 'deallocate-profile-inc.f90'
  include 'deallocate-inc.f90'
end subroutine i3_all_free

subroutine i4_all(array,m)
  use metadata_interfaces, metadata_address => geti4
  implicit none
  type(malloc_information_all), intent(in) :: m
  integer(kind=4), dimension(:,:,:,:), allocatable, intent(inout) :: array
  include 'allocate-profile-inc.f90'
  !allocate the array
  allocate(array(m%lbounds(1):m%ubounds(1),m%lbounds(2):m%ubounds(2),&
       m%lbounds(3):m%ubounds(3),m%lbounds(4):m%ubounds(4)+ndebug),stat=ierror)
  include 'allocate-inc.f90'
end subroutine i4_all

subroutine i4_all_free(array)
  use metadata_interfaces, metadata_address => geti4
  implicit none
  integer(kind=4), dimension(:,:,:,:), allocatable, intent(inout) :: array
  include 'deallocate-profile-inc.f90'
  include 'deallocate-inc.f90'
end subroutine i4_all_free

!!$! long integer arrays
!!$subroutine il1_all(array,m)
!!$  use metadata_interfaces, metadata_address => getil1
!!$  implicit none
!!$  type(malloc_information_all), intent(in) :: m
!!$  integer(kind=8), dimension(:), allocatable, intent(inout) :: array
!!$  include 'allocate-profile-inc.f90'
!!$  !allocate the array
!!$  allocate(array(m%lbounds(1):m%ubounds(1)+ndebug),stat=ierror)
!!$  include 'allocate-inc.f90'
!!$end subroutine il1_all
!!$
!!$subroutine il1_all_free(array)
!!$  use metadata_interfaces, metadata_address => getil1
!!$  implicit none
!!$  integer(kind=8), dimension(:), allocatable, intent(inout) :: array
!!$  include 'deallocate-profile-inc.f90'
!!$  include 'deallocate-inc.f90'
!!$end subroutine il1_all_free

subroutine il2_all(array,m)
  use metadata_interfaces, metadata_address => getil2
  implicit none
  type(malloc_information_all), intent(in) :: m
  integer(kind=8), dimension(:,:), allocatable, intent(inout) :: array
  include 'allocate-profile-inc.f90'
  !allocate the array
  allocate(array(m%lbounds(1):m%ubounds(1),m%lbounds(2):m%ubounds(2)+ndebug),stat=ierror)
  include 'allocate-inc.f90'
end subroutine il2_all

subroutine il2_all_free(array)
  use metadata_interfaces, metadata_address => getil2
  implicit none
  integer(kind=8), dimension(:,:), allocatable, intent(inout) :: array
  include 'deallocate-profile-inc.f90'
  include 'deallocate-inc.f90'
end subroutine il2_all_free

!character arrays
subroutine c1_all(array,m)
  use metadata_interfaces, metadata_address => getc1
  implicit none
  type(malloc_information_str_all), intent(in) :: m
  character(len=m%len), dimension(:), allocatable, intent(inout) :: array
  include 'allocate-profile-inc.f90'
  if (f_nan_pad_size > 0) then
     padding=f_nan_pad_size
     call togglepadding(product(int(m%shape(1:m%rank-1),f_long))*&
          m%len*kind(array)*(m%shape(m%rank)+padding))
  end if
  !allocate the array
  allocate(array(m%lbounds(1):m%ubounds(1)+ndebug),stat=ierror)
  include 'allocate-c-inc.f90'
  !include 'allocate-inc.f90'
end subroutine c1_all

subroutine c2_all(array,m)
  use metadata_interfaces, metadata_address => getc1
  implicit none
  type(malloc_information_str_all), intent(in) :: m
  character(len=m%len), dimension(:,:), allocatable, intent(inout) :: array
  include 'allocate-profile-inc.f90'
  if (f_nan_pad_size > 0) then
     padding=f_nan_pad_size
     call togglepadding(product(int(m%shape(1:m%rank-1),f_long))*&
          m%len*kind(array)*(m%shape(m%rank)+padding))
  end if
  !allocate the array
  allocate(array(m%lbounds(1):m%ubounds(1),m%lbounds(2):m%ubounds(2)+ndebug),stat=ierror)
  include 'allocate-c-inc.f90'
  !include 'allocate-inc.f90'
end subroutine c2_all

subroutine c1_all_free(length,array)
!subroutine f_free_str(length,array)
  use metadata_interfaces, metadata_address => getc1
  implicit none
  integer, intent(in) :: length !< need to specify length for the declaration below (sometimes fortran runtime error)
  character(len=length), dimension(:), allocatable, intent(inout) :: array
  include 'deallocate-profile-inc.f90'
  !include 'deallocate-c-inc.f90'
  include 'deallocate-inc.f90'
end subroutine c1_all_free

subroutine c2_all_free(length,array)
  use metadata_interfaces, metadata_address => getc1
  implicit none
  integer, intent(in) :: length !< need to specify length for the declaration below (sometimes fortran runtime error)
  character(len=length), dimension(:,:), allocatable, intent(inout) :: array
  include 'deallocate-profile-inc.f90'
  !include 'deallocate-c-inc.f90'
  include 'deallocate-inc.f90'
end subroutine c2_all_free

subroutine b1_all(array,m)
  use metadata_interfaces, metadata_address => getl1
  implicit none
  type(malloc_information_all), intent(in) :: m
  logical(f_byte), dimension(:), allocatable, intent(inout) :: array
  include 'allocate-profile-inc.f90'
  !allocate the array
  allocate(array(m%lbounds(1):m%ubounds(1)+ndebug),stat=ierror)
  include 'allocate-inc.f90'
end subroutine b1_all

subroutine b1_all_free(array)
  use metadata_interfaces, metadata_address => getl1
  implicit none
  logical(f_byte), dimension(:), allocatable, intent(inout) :: array
  include 'deallocate-profile-inc.f90'
  include 'deallocate-inc.f90'
end subroutine b1_all_free

subroutine b2_all(array,m)
  use metadata_interfaces, metadata_address => getl2
  implicit none
  type(malloc_information_all), intent(in) :: m
  logical(f_byte), dimension(:,:), allocatable, intent(inout) :: array
  include 'allocate-profile-inc.f90'
  !allocate the array
  allocate(array(m%lbounds(1):m%ubounds(1),m%lbounds(2):m%ubounds(2)+ndebug),stat=ierror)
  include 'allocate-inc.f90'
end subroutine b2_all

subroutine b2_all_free(array)
  use metadata_interfaces, metadata_address => getl2
  implicit none
  logical(f_byte), dimension(:,:), allocatable, intent(inout) :: array
  include 'deallocate-profile-inc.f90'
  include 'deallocate-inc.f90'
end subroutine b2_all_free

subroutine b3_all(array,m)
  use metadata_interfaces, metadata_address => getl3
  implicit none
  type(malloc_information_all), intent(in) :: m
  logical(f_byte), dimension(:,:,:), allocatable, intent(inout) :: array
  include 'allocate-profile-inc.f90'
  !allocate the array
  allocate(array(m%lbounds(1):m%ubounds(1),m%lbounds(2):m%ubounds(2),&
       m%lbounds(3):m%ubounds(3)+ndebug),stat=ierror)
  include 'allocate-inc.f90'
end subroutine b3_all

subroutine b3_all_free(array)
  use metadata_interfaces, metadata_address => getl3
  implicit none
  logical(f_byte), dimension(:,:,:), allocatable, intent(inout) :: array
  include 'deallocate-profile-inc.f90'
  include 'deallocate-inc.f90'
end subroutine b3_all_free


subroutine l1_all(array,m)
  use metadata_interfaces, metadata_address => getl1
  implicit none
  type(malloc_information_all), intent(in) :: m
  logical, dimension(:), allocatable, intent(inout) :: array
  include 'allocate-profile-inc.f90'
  !allocate the array
  allocate(array(m%lbounds(1):m%ubounds(1)+ndebug),stat=ierror)
  include 'allocate-inc.f90'
end subroutine l1_all

subroutine l1_all_free(array)
  use metadata_interfaces, metadata_address => getl1
  implicit none
  logical, dimension(:), allocatable, intent(inout) :: array
  include 'deallocate-profile-inc.f90'
  include 'deallocate-inc.f90'
end subroutine l1_all_free

subroutine l2_all(array,m)
  use metadata_interfaces, metadata_address => getl2
  implicit none
  type(malloc_information_all), intent(in) :: m
  logical, dimension(:,:), allocatable, intent(inout) :: array
  include 'allocate-profile-inc.f90'
  !allocate the array
  allocate(array(m%lbounds(1):m%ubounds(1),m%lbounds(2):m%ubounds(2)+ndebug),stat=ierror)
  include 'allocate-inc.f90'
end subroutine l2_all

subroutine l2_all_free(array)
  use metadata_interfaces, metadata_address => getl2
  implicit none
  logical, dimension(:,:), allocatable, intent(inout) :: array
  include 'deallocate-profile-inc.f90'
  include 'deallocate-inc.f90'
end subroutine l2_all_free

subroutine l3_all(array,m)
  use metadata_interfaces, metadata_address => getl3
  implicit none
  type(malloc_information_all), intent(in) :: m
  logical, dimension(:,:,:), allocatable, intent(inout) :: array
  include 'allocate-profile-inc.f90'
  !allocate the array
  allocate(array(m%lbounds(1):m%ubounds(1),m%lbounds(2):m%ubounds(2),&
       m%lbounds(3):m%ubounds(3)+ndebug),stat=ierror)
  include 'allocate-inc.f90'
end subroutine l3_all

subroutine l3_all_free(array)
  use metadata_interfaces, metadata_address => getl3
  implicit none
  logical, dimension(:,:,:), allocatable, intent(inout) :: array
  include 'deallocate-profile-inc.f90'
  include 'deallocate-inc.f90'
end subroutine l3_all_free

subroutine r1_all(array,m)
  use metadata_interfaces, metadata_address => getr1
  implicit none
  type(malloc_information_all), intent(in) :: m
  real, dimension(:), allocatable, intent(inout) :: array
  include 'allocate-profile-inc.f90'
  !allocate the array
  allocate(array(m%lbounds(1):m%ubounds(1)+ndebug),stat=ierror)
  include 'allocate-inc.f90'
end subroutine r1_all

subroutine r1_all_free(array)
  use metadata_interfaces, metadata_address => getr1
  implicit none
  real, dimension(:), allocatable, intent(inout) :: array
  include 'deallocate-profile-inc.f90'
  include 'deallocate-inc.f90'
end subroutine r1_all_free

subroutine r2_all(array,m)
  use metadata_interfaces, metadata_address => getr2
  implicit none
  type(malloc_information_all), intent(in) :: m
  real, dimension(:,:), allocatable, intent(inout) :: array
  include 'allocate-profile-inc.f90'
  allocate(array(m%lbounds(1):m%ubounds(1),m%lbounds(2):m%ubounds(2)+ndebug),stat=ierror)
  include 'allocate-inc.f90'
end subroutine r2_all

subroutine r2_all_free(array)
  use metadata_interfaces, metadata_address => getr2
  implicit none
  real, dimension(:,:), allocatable, intent(inout) :: array
  include 'deallocate-profile-inc.f90'
  include 'deallocate-inc.f90'
end subroutine r2_all_free

subroutine r3_all(array,m)
  use metadata_interfaces, metadata_address => getr3
  implicit none
  type(malloc_information_all), intent(in) :: m
  real, dimension(:,:,:), allocatable, intent(inout) :: array
  !local variables
  include 'allocate-profile-inc.f90'
  allocate(array(m%lbounds(1):m%ubounds(1),&
       m%lbounds(2):m%ubounds(2),m%lbounds(3):m%ubounds(3)+ndebug),stat=ierror)
  include 'allocate-inc.f90'
end subroutine r3_all

subroutine r3_all_free(array)
  use metadata_interfaces, metadata_address => getr3
  implicit none
  real, dimension(:,:,:), allocatable, intent(inout) :: array
  include 'deallocate-profile-inc.f90'
  include 'deallocate-inc.f90'
end subroutine r3_all_free

subroutine r4_all(array,m)
  use metadata_interfaces, metadata_address => getr4
  implicit none
  type(malloc_information_all), intent(in) :: m
  real, dimension(:,:,:,:), allocatable, intent(inout) :: array
  !local variables
  include 'allocate-profile-inc.f90'
  allocate(array(m%lbounds(1):m%ubounds(1),&
       m%lbounds(2):m%ubounds(2),m%lbounds(3):m%ubounds(3),&
       m%lbounds(4):m%ubounds(4)+ndebug),stat=ierror)
  include 'allocate-inc.f90'
end subroutine r4_all

subroutine r4_all_free(array)
  use metadata_interfaces, metadata_address => getr4
  implicit none
  real, dimension(:,:,:,:), allocatable, intent(inout) :: array
  include 'deallocate-profile-inc.f90'
  include 'deallocate-inc.f90'
end subroutine r4_all_free


subroutine d1_all(array,m)
  use metadata_interfaces, metadata_address => getdp1
  implicit none
  type(malloc_information_all), intent(in) :: m
  double precision, dimension(:), allocatable, intent(inout) :: array
  include 'allocate-profile-inc.f90'
  !allocate the array
  allocate(array(m%lbounds(1):m%ubounds(1)+ndebug),stat=ierror)
  include 'allocate-inc.f90'
end subroutine d1_all

subroutine d1_all_free(array)
  use metadata_interfaces, metadata_address => getdp1
  implicit none
  double precision, dimension(:), allocatable, intent(inout) :: array
  include 'deallocate-profile-inc.f90'
  include 'deallocate-inc.f90'
end subroutine d1_all_free


subroutine d2_all(array,m)
  use metadata_interfaces, metadata_address => getdp2
  implicit none
  type(malloc_information_all), intent(in) :: m
  double precision, dimension(:,:), allocatable, intent(inout) :: array
  include 'allocate-profile-inc.f90'
  allocate(array(m%lbounds(1):m%ubounds(1),m%lbounds(2):m%ubounds(2)+ndebug),stat=ierror)
  include 'allocate-inc.f90'
end subroutine d2_all

subroutine d2_all_free(array)
  use metadata_interfaces, metadata_address => getdp2
  implicit none
  double precision, dimension(:,:), allocatable, intent(inout) :: array
  include 'deallocate-profile-inc.f90'
  include 'deallocate-inc.f90'
end subroutine d2_all_free

subroutine d3_all(array,m)
  use metadata_interfaces, metadata_address => getdp3
  implicit none
  type(malloc_information_all), intent(in) :: m
  double precision, dimension(:,:,:), allocatable, intent(inout) :: array
  !local variables
  include 'allocate-profile-inc.f90'
  allocate(array(m%lbounds(1):m%ubounds(1),&
       m%lbounds(2):m%ubounds(2),m%lbounds(3):m%ubounds(3)+ndebug),stat=ierror)
  include 'allocate-inc.f90'
end subroutine d3_all

subroutine d3_all_free(array)
  use metadata_interfaces, metadata_address => getdp3
  implicit none
  double precision, dimension(:,:,:), allocatable, intent(inout) :: array
  include 'deallocate-profile-inc.f90'
  include 'deallocate-inc.f90'
end subroutine d3_all_free

subroutine d4_all(array,m)
  use metadata_interfaces, metadata_address => getdp4
  implicit none
  type(malloc_information_all), intent(in) :: m
  double precision, dimension(:,:,:,:), allocatable, intent(inout) :: array
  !local variables
  include 'allocate-profile-inc.f90'
  allocate(array(m%lbounds(1):m%ubounds(1),&
       m%lbounds(2):m%ubounds(2),m%lbounds(3):m%ubounds(3),&
       m%lbounds(4):m%ubounds(4)+ndebug),stat=ierror)
  include 'allocate-inc.f90'
end subroutine d4_all

subroutine d4_all_free(array)
  use metadata_interfaces, metadata_address => getdp4
  implicit none
  double precision, dimension(:,:,:,:), allocatable, intent(inout) :: array
  include 'deallocate-profile-inc.f90'
  include 'deallocate-inc.f90'
end subroutine d4_all_free

subroutine d5_all(array,m)
  use metadata_interfaces, metadata_address => getdp5
  implicit none
  type(malloc_information_all), intent(in) :: m
  double precision, dimension(:,:,:,:,:), allocatable, intent(inout) :: array
  !local variables
  include 'allocate-profile-inc.f90'
  allocate(array(m%lbounds(1):m%ubounds(1),&
       m%lbounds(2):m%ubounds(2),m%lbounds(3):m%ubounds(3),&
       m%lbounds(4):m%ubounds(4),m%lbounds(5):m%ubounds(5)+ndebug),stat=ierror)
  include 'allocate-inc.f90'
end subroutine d5_all

subroutine d5_all_free(array)
  use metadata_interfaces, metadata_address => getdp5
  implicit none
  double precision, dimension(:,:,:,:,:), allocatable, intent(inout) :: array
  include 'deallocate-profile-inc.f90'
  include 'deallocate-inc.f90'
end subroutine d5_all_free

subroutine d6_all(array,m)
  use metadata_interfaces, metadata_address => getdp6
  implicit none
  type(malloc_information_all), intent(in) :: m
  double precision, dimension(:,:,:,:,:,:), allocatable, intent(inout) :: array
  !local variables
  include 'allocate-profile-inc.f90'
  allocate(array(m%lbounds(1):m%ubounds(1),&
       m%lbounds(2):m%ubounds(2),m%lbounds(3):m%ubounds(3),&
       m%lbounds(4):m%ubounds(4),m%lbounds(5):m%ubounds(5),&
       m%lbounds(6):m%ubounds(6)+ndebug),stat=ierror)
  include 'allocate-inc.f90'
end subroutine d6_all

subroutine d6_all_free(array)
  use metadata_interfaces, metadata_address => getdp6
  implicit none
  double precision, dimension(:,:,:,:,:,:), allocatable, intent(inout) :: array
  include 'deallocate-profile-inc.f90'
  include 'deallocate-inc.f90'
end subroutine d6_all_free

subroutine d7_all(array,m)
  use metadata_interfaces, metadata_address => getdp7
  implicit none
  type(malloc_information_all), intent(in) :: m
  double precision, dimension(:,:,:,:,:,:,:), allocatable, intent(inout) :: array
  !local variables
  include 'allocate-profile-inc.f90'
  allocate(array(m%lbounds(1):m%ubounds(1),&
       m%lbounds(2):m%ubounds(2),m%lbounds(3):m%ubounds(3),&
       m%lbounds(4):m%ubounds(4),m%lbounds(5):m%ubounds(5),&
       m%lbounds(6):m%ubounds(6),m%lbounds(7):m%ubounds(7)+ndebug),stat=ierror)
  include 'allocate-inc.f90'
end subroutine d7_all

subroutine d7_all_free(array)
  use metadata_interfaces, metadata_address => getdp7
  implicit none
  double precision, dimension(:,:,:,:,:,:,:), allocatable, intent(inout) :: array
  include 'deallocate-profile-inc.f90'
  include 'deallocate-inc.f90'
end subroutine d7_all_free

subroutine li1_all(array,m)
  use metadata_interfaces, metadata_address => geti1
  implicit none
  type(malloc_information_all), intent(in) :: m
  integer(kind=8), dimension(:), allocatable, intent(inout) :: array
  include 'allocate-profile-inc.f90'
  !allocate the array
  allocate(array(m%lbounds(1):m%ubounds(1)+ndebug),stat=ierror)
  include 'allocate-inc.f90'
end subroutine li1_all

subroutine li1_all_free(array)
  use metadata_interfaces, metadata_address => geti1
  implicit none
  integer(kind=8), dimension(:), allocatable, intent(inout) :: array
  include 'deallocate-profile-inc.f90'
  include 'deallocate-inc.f90'
end subroutine li1_all_free

subroutine li2_all(array,m)
  use metadata_interfaces, metadata_address => geti2
  implicit none
  type(malloc_information_all), intent(in) :: m
  integer(kind=8), dimension(:,:), allocatable, intent(inout) :: array
  include 'allocate-profile-inc.f90'
  !allocate the array
  allocate(array(m%lbounds(1):m%ubounds(1),m%lbounds(2):m%ubounds(2)+ndebug),stat=ierror)
  include 'allocate-inc.f90'
end subroutine li2_all

subroutine li2_all_free(array)
  use metadata_interfaces, metadata_address => geti2
  implicit none
  integer(kind=8), dimension(:,:), allocatable, intent(inout) :: array
  include 'deallocate-profile-inc.f90'
  include 'deallocate-inc.f90'
end subroutine li2_all_free

subroutine li3_all(array,m)
  use metadata_interfaces, metadata_address => geti3
  implicit none
  type(malloc_information_all), intent(in) :: m
  integer(kind=8), dimension(:,:,:), allocatable, intent(inout) :: array
  include 'allocate-profile-inc.f90'
  !allocate the array
  allocate(array(m%lbounds(1):m%ubounds(1),m%lbounds(2):m%ubounds(2),m%lbounds(3):m%ubounds(3)+ndebug),stat=ierror)
  include 'allocate-inc.f90'
end subroutine li3_all

subroutine li3_all_free(array)
  use metadata_interfaces, metadata_address => geti3
  implicit none
  integer(kind=8), dimension(:,:,:), allocatable, intent(inout) :: array
  include 'deallocate-profile-inc.f90'
  include 'deallocate-inc.f90'
end subroutine li3_all_free

subroutine li4_all(array,m)
  use metadata_interfaces, metadata_address => geti4
  implicit none
  type(malloc_information_all), intent(in) :: m
  integer(kind=8), dimension(:,:,:,:), allocatable, intent(inout) :: array
  include 'allocate-profile-inc.f90'
  !allocate the array
  allocate(array(m%lbounds(1):m%ubounds(1),m%lbounds(2):m%ubounds(2),&
       m%lbounds(3):m%ubounds(3),m%lbounds(4):m%ubounds(4)+ndebug),stat=ierror)
  include 'allocate-inc.f90'
end subroutine li4_all

subroutine li4_all_free(array)
  use metadata_interfaces, metadata_address => geti4
  implicit none
  integer(kind=8), dimension(:,:,:,:), allocatable, intent(inout) :: array
  include 'deallocate-profile-inc.f90'
  include 'deallocate-inc.f90'
end subroutine li4_all_free

!test to see if this is convenient
subroutine d1_all_free_multi(arrayA,arrayB,arrayC,arrayD,arrayE,arrayF,arrayG,arrayH)
  implicit none
  double precision, dimension(:), allocatable, intent(inout) :: arrayA
  double precision, dimension(:), allocatable, intent(inout) :: arrayB
  double precision, dimension(:), allocatable, optional, intent(inout) :: arrayC
  double precision, dimension(:), allocatable, optional, intent(inout) :: arrayD
  double precision, dimension(:), allocatable, optional, intent(inout) :: arrayE
  double precision, dimension(:), allocatable, optional, intent(inout) :: arrayF
  double precision, dimension(:), allocatable, optional, intent(inout) :: arrayG
  double precision, dimension(:), allocatable, optional, intent(inout) :: arrayH

  include 'deallocate-multiple-inc.f90'

end subroutine d1_all_free_multi

!test to see if this is convenient
subroutine i1_all_free_multi(arrayA,arrayB,arrayC,arrayD,arrayE,arrayF,arrayG,arrayH)
  implicit none
  integer, dimension(:), allocatable, intent(inout) :: arrayA
  integer, dimension(:), allocatable, intent(inout) :: arrayB
  integer, dimension(:), allocatable, optional, intent(inout) :: arrayC
  integer, dimension(:), allocatable, optional, intent(inout) :: arrayD
  integer, dimension(:), allocatable, optional, intent(inout) :: arrayE
  integer, dimension(:), allocatable, optional, intent(inout) :: arrayF
  integer, dimension(:), allocatable, optional, intent(inout) :: arrayG
  integer, dimension(:), allocatable, optional, intent(inout) :: arrayH

  include 'deallocate-multiple-inc.f90'

end subroutine i1_all_free_multi

! double complex
subroutine z2_all(array,m)
  use metadata_interfaces, metadata_address => getz2
  implicit none
  type(malloc_information_all), intent(in) :: m
  double complex, dimension(:,:), allocatable, intent(inout) :: array
  include 'allocate-profile-inc.f90'
  allocate(array(m%lbounds(1):m%ubounds(1),m%lbounds(2):m%ubounds(2)+ndebug),stat=ierror)
  include 'allocate-inc.f90'
end subroutine z2_all

subroutine z2_all_free(array)
  use metadata_interfaces, metadata_address => getz2
  implicit none
  double complex, dimension(:,:), allocatable, intent(inout) :: array
  include 'deallocate-profile-inc.f90'
  include 'deallocate-inc.f90'
end subroutine z2_all_free

subroutine z3_all(array,m)
  use metadata_interfaces, metadata_address => getz3
  implicit none
  type(malloc_information_all), intent(in) :: m
  double complex, dimension(:,:,:), allocatable, intent(inout) :: array
  include 'allocate-profile-inc.f90'
  allocate(array(m%lbounds(1):m%ubounds(1),m%lbounds(2):m%ubounds(2), &
&                m%lbounds(3):m%ubounds(3)+ndebug),stat=ierror)
  include 'allocate-inc.f90'
end subroutine z3_all

subroutine z3_all_free(array)
  use metadata_interfaces, metadata_address => getz3
  implicit none
  double complex, dimension(:,:,:), allocatable, intent(inout) :: array
  include 'deallocate-profile-inc.f90'
  include 'deallocate-inc.f90'
end subroutine z3_all_free

!pointers
subroutine d1_ptr(array,m)
  use metadata_interfaces, metadata_address => getdp1ptr
  implicit none
  type(malloc_information_ptr), intent(in) :: m
  double precision, dimension(:), pointer, intent(inout) :: array
  double precision :: d
  type(c_ptr) :: p
  !local variables
  include 'allocate-profile-inc.f90'
  include 'allocate-ptr-inc.f90'
  include 'allocate-aligned-inc.f90'
  include 'allocate-simgrid-inc.f90'
  !allocate the array
  allocate(array(m%lbounds(1):m%ubounds(1)+ndebug),stat=ierror)
  include 'allocate-inc.f90'
end subroutine d1_ptr

subroutine d1_ptr_free(array)
  use metadata_interfaces, metadata_address => getdp1ptr
  implicit none
  double precision, dimension(:), pointer, intent(inout) :: array
  include 'deallocate-profile-inc.f90'
  if (.not. associated(array)) return
  include 'deallocate-ptr-inc.f90'
  nullify(array)
end subroutine d1_ptr_free

subroutine i1_ptr(array,m)
  use metadata_interfaces, metadata_address => geti1ptr
  implicit none
  type(malloc_information_ptr), intent(in) :: m
  integer, dimension(:), pointer, intent(inout) :: array
  !local variables
  include 'allocate-profile-inc.f90'
  include 'allocate-ptr-inc.f90'
  !allocate the array
  allocate(array(m%lbounds(1):m%ubounds(1)+ndebug),stat=ierror)
  include 'allocate-inc.f90'
end subroutine i1_ptr

subroutine i1_ptr_free(array)
  use metadata_interfaces, metadata_address => geti1ptr
  implicit none
  integer, dimension(:), pointer, intent(inout) :: array
  include 'deallocate-profile-inc.f90'
  if (.not. associated(array)) return
  include 'deallocate-ptr-inc.f90'
  nullify(array)
end subroutine i1_ptr_free

!test to see if this is convenient
subroutine i1_ptr_free_multi(arrayA,arrayB,arrayC,arrayD,arrayE,arrayF,arrayG,arrayH)
  implicit none
  integer, dimension(:), pointer, intent(inout) :: arrayA
  integer, dimension(:), pointer, intent(inout) :: arrayB
  integer, dimension(:), pointer, optional, intent(inout) :: arrayC
  integer, dimension(:), pointer, optional, intent(inout) :: arrayD
  integer, dimension(:), pointer, optional, intent(inout) :: arrayE
  integer, dimension(:), pointer, optional, intent(inout) :: arrayF
  integer, dimension(:), pointer, optional, intent(inout) :: arrayG
  integer, dimension(:), pointer, optional, intent(inout) :: arrayH

  include 'deallocate-multiple-ptr-inc.f90'

end subroutine i1_ptr_free_multi

subroutine d2_ptr(array,m)
  use metadata_interfaces, metadata_address => getdp2ptr
  implicit none
  type(malloc_information_ptr), intent(in) :: m
  double precision :: d
  type(c_ptr) :: p
  double precision, dimension(:,:), pointer, intent(inout) :: array
  include 'allocate-profile-inc.f90'
  include 'allocate-ptr-inc.f90'
  include 'allocate-aligned-inc.f90'
  include 'allocate-simgrid-inc.f90'
  !allocate the array
  allocate(array(m%lbounds(1):m%ubounds(1),m%lbounds(2):m%ubounds(2)+ndebug),stat=ierror)
  include 'allocate-inc.f90'
end subroutine d2_ptr

subroutine d2_ptr_free(array)
  use metadata_interfaces, metadata_address => getdp2ptr
  implicit none
  double precision, dimension(:,:), pointer, intent(inout) :: array
  include 'deallocate-profile-inc.f90'
  if (.not. associated(array)) return
  include 'deallocate-ptr-inc.f90'
  nullify(array)
end subroutine d2_ptr_free

subroutine i2_ptr(array,m)
  use metadata_interfaces, metadata_address => geti2ptr
  implicit none
  type(malloc_information_ptr), intent(in) :: m
  integer, dimension(:,:), pointer, intent(inout) :: array
  include 'allocate-profile-inc.f90'
  include 'allocate-ptr-inc.f90'
  !allocate the array
  allocate(array(m%lbounds(1):m%ubounds(1),m%lbounds(2):m%ubounds(2)+ndebug),stat=ierror)
  include 'allocate-inc.f90'
end subroutine i2_ptr

subroutine i2_ptr_free(array)
  use metadata_interfaces, metadata_address => geti2ptr
  implicit none
  integer, dimension(:,:), pointer, intent(inout) :: array
  include 'deallocate-profile-inc.f90'
  if (.not. associated(array)) return
  include 'deallocate-ptr-inc.f90'
  nullify(array)
end subroutine i2_ptr_free


subroutine d3_ptr(array,m)
  use metadata_interfaces, metadata_address => getdp3ptr
  implicit none
  type(malloc_information_ptr), intent(in) :: m
  double precision, dimension(:,:,:), pointer, intent(inout) :: array
  double precision :: d
  type(c_ptr) :: p
  include 'allocate-profile-inc.f90'
  include 'allocate-ptr-inc.f90'
  include 'allocate-aligned-inc.f90'
  include 'allocate-simgrid-inc.f90'
  !allocate the array
  allocate(array(m%lbounds(1):m%ubounds(1),m%lbounds(2):m%ubounds(2),&
       m%lbounds(3):m%ubounds(3)+ndebug),stat=ierror)
  include 'allocate-inc.f90'
end subroutine d3_ptr

subroutine d3_ptr_free(array)
  use metadata_interfaces, metadata_address => getdp3ptr
  implicit none
  double precision, dimension(:,:,:), pointer, intent(inout) :: array
  include 'deallocate-profile-inc.f90'
  if (.not. associated(array)) return
  include 'deallocate-ptr-inc.f90'
  nullify(array)
end subroutine d3_ptr_free

subroutine i3_ptr(array,m)
  use metadata_interfaces, metadata_address => geti3ptr
  implicit none
  type(malloc_information_ptr), intent(in) :: m
  integer, dimension(:,:,:), pointer, intent(inout) :: array
  include 'allocate-profile-inc.f90'
  include 'allocate-ptr-inc.f90'
  !allocate the array
  allocate(array(m%lbounds(1):m%ubounds(1),m%lbounds(2):m%ubounds(2),&
       m%lbounds(3):m%ubounds(3)+ndebug),stat=ierror)
  include 'allocate-inc.f90'
end subroutine i3_ptr

subroutine i3_ptr_free(array)
  use metadata_interfaces, metadata_address => geti3ptr
  implicit none
  integer, dimension(:,:,:), pointer, intent(inout) :: array
  include 'deallocate-profile-inc.f90'
  if (.not. associated(array)) return
  include 'deallocate-ptr-inc.f90'
  nullify(array)
end subroutine i3_ptr_free

subroutine i4_ptr(array,m)
  use metadata_interfaces, metadata_address => geti4ptr
  implicit none
  type(malloc_information_ptr), intent(in) :: m
  integer, dimension(:,:,:,:), pointer, intent(inout) :: array
  include 'allocate-profile-inc.f90'
  include 'allocate-ptr-inc.f90'
  !allocate the array
  allocate(array(m%lbounds(1):m%ubounds(1),m%lbounds(2):m%ubounds(2),&
       m%lbounds(3):m%ubounds(3),m%lbounds(4):m%ubounds(4)+ndebug),stat=ierror)
  include 'allocate-inc.f90'
end subroutine i4_ptr

subroutine i4_ptr_free(array)
  use metadata_interfaces, metadata_address => geti4ptr
  implicit none
  integer, dimension(:,:,:,:), pointer, intent(inout) :: array
  include 'deallocate-profile-inc.f90'
  if (.not. associated(array)) return
  include 'deallocate-ptr-inc.f90'
  nullify(array)
end subroutine i4_ptr_free

subroutine li1_ptr(array,m)
  use metadata_interfaces, metadata_address => getl2ptr
  implicit none
  type(malloc_information_ptr), intent(in) :: m
  integer(f_long), dimension(:), pointer, intent(inout) :: array
  include 'allocate-profile-inc.f90'
  include 'allocate-ptr-inc.f90'
  !allocate the array
  allocate(array(m%lbounds(1):m%ubounds(1)+ndebug),stat=ierror)
  include 'allocate-inc.f90'
end subroutine li1_ptr

subroutine li1_ptr_free(array)
  use metadata_interfaces, metadata_address => getl2ptr
  implicit none
  integer(f_long), dimension(:), pointer, intent(inout) :: array
  include 'deallocate-profile-inc.f90'
  if (.not. associated(array)) return
  include 'deallocate-ptr-inc.f90'
end subroutine li1_ptr_free


subroutine d4_ptr(array,m)
  use metadata_interfaces, metadata_address => getdp4ptr
  implicit none
  type(malloc_information_ptr), intent(in) :: m
  double precision, dimension(:,:,:,:), pointer, intent(inout) :: array
  double precision :: d
  type(c_ptr) :: p
  include 'allocate-profile-inc.f90'
  include 'allocate-ptr-inc.f90'
  include 'allocate-aligned-inc.f90'
  include 'allocate-simgrid-inc.f90'
  !allocate the array
  allocate(array(m%lbounds(1):m%ubounds(1),m%lbounds(2):m%ubounds(2),&
       m%lbounds(3):m%ubounds(3),m%lbounds(4):m%ubounds(4)+ndebug),stat=ierror)
  include 'allocate-inc.f90'
end subroutine d4_ptr

subroutine d4_ptr_free(array)
  use metadata_interfaces, metadata_address => getdp4ptr
  implicit none
  double precision, dimension(:,:,:,:), pointer, intent(inout) :: array
  include 'deallocate-profile-inc.f90'
  if (.not. associated(array)) return
  include 'deallocate-ptr-inc.f90'
  nullify(array)
end subroutine d4_ptr_free

subroutine d5_ptr(array,m)
  use metadata_interfaces, metadata_address => getdp5ptr
  implicit none
  type(malloc_information_ptr), intent(in) :: m
  double precision, dimension(:,:,:,:,:), pointer, intent(inout) :: array
  double precision :: d
  type(c_ptr) :: p
  include 'allocate-profile-inc.f90'
  include 'allocate-ptr-inc.f90'
  include 'allocate-aligned-inc.f90'
  include 'allocate-simgrid-inc.f90'
  !allocate the array
  allocate(array(m%lbounds(1):m%ubounds(1),m%lbounds(2):m%ubounds(2),&
       m%lbounds(3):m%ubounds(3),m%lbounds(4):m%ubounds(4),&
       m%lbounds(5):m%ubounds(5)+ndebug),stat=ierror)
  include 'allocate-inc.f90'
end subroutine d5_ptr

subroutine d5_ptr_free(array)
  use metadata_interfaces, metadata_address => getdp5ptr
  implicit none
  double precision, dimension(:,:,:,:,:), pointer, intent(inout) :: array
  include 'deallocate-profile-inc.f90'
  if (.not. associated(array)) return
  include 'deallocate-ptr-inc.f90'
  nullify(array)
end subroutine d5_ptr_free

subroutine d6_ptr(array,m)
  use metadata_interfaces, metadata_address => getdp6ptr
  implicit none
  type(malloc_information_ptr), intent(in) :: m
  double precision, dimension(:,:,:,:,:,:), pointer, intent(inout) :: array
  double precision :: d
  type(c_ptr) :: p
  include 'allocate-profile-inc.f90'
  include 'allocate-ptr-inc.f90'
  include 'allocate-aligned-inc.f90'
  include 'allocate-simgrid-inc.f90'
  !allocate the array
  allocate(array(m%lbounds(1):m%ubounds(1),m%lbounds(2):m%ubounds(2),&
       m%lbounds(3):m%ubounds(3),m%lbounds(4):m%ubounds(4),&
       m%lbounds(5):m%ubounds(5),m%lbounds(6):m%ubounds(6)+ndebug),stat=ierror)
  include 'allocate-inc.f90'
end subroutine d6_ptr

subroutine d6_ptr_free(array)
  use metadata_interfaces, metadata_address => getdp6ptr
  implicit none
  double precision, dimension(:,:,:,:,:,:), pointer, intent(inout) :: array
  include 'deallocate-profile-inc.f90'
  if (.not. associated(array)) return
  include 'deallocate-ptr-inc.f90'
  nullify(array)
end subroutine d6_ptr_free

subroutine l1_ptr(array,m)
  use metadata_interfaces, metadata_address => getl2ptr
  implicit none
  type(malloc_information_ptr), intent(in) :: m
  logical, dimension(:), pointer, intent(inout) :: array
  include 'allocate-profile-inc.f90'
  include 'allocate-ptr-inc.f90'
  !allocate the array
  allocate(array(m%lbounds(1):m%ubounds(1)+ndebug),stat=ierror)
  include 'allocate-inc.f90'
end subroutine l1_ptr

subroutine l1_ptr_free(array)
  use metadata_interfaces, metadata_address => getl2ptr
  implicit none
  logical, dimension(:), pointer, intent(inout) :: array
  include 'deallocate-profile-inc.f90'
  if (.not. associated(array)) return
  include 'deallocate-ptr-inc.f90'
end subroutine l1_ptr_free

subroutine l2_ptr(array,m)
  use metadata_interfaces, metadata_address => getl2ptr
  implicit none
  type(malloc_information_ptr), intent(in) :: m
  logical, dimension(:,:), pointer, intent(inout) :: array
  include 'allocate-profile-inc.f90'
  include 'allocate-ptr-inc.f90'
  !allocate the array
  allocate(array(m%lbounds(1):m%ubounds(1),m%lbounds(2):m%ubounds(2)+ndebug),stat=ierror)
  include 'allocate-inc.f90'
end subroutine l2_ptr

subroutine l3_ptr(array,m)
  use metadata_interfaces, metadata_address => getl3ptr
  implicit none
  type(malloc_information_ptr), intent(in) :: m
  logical, dimension(:,:,:), pointer, intent(inout) :: array
  include 'allocate-profile-inc.f90'
  include 'allocate-ptr-inc.f90'
  !allocate the array
  allocate(array(m%lbounds(1):m%ubounds(1),m%lbounds(2):m%ubounds(2),&
       m%lbounds(3):m%ubounds(3)+ndebug),stat=ierror)
  include 'allocate-inc.f90'
end subroutine l3_ptr

subroutine l2_ptr_free(array)
  use metadata_interfaces, metadata_address => getl2ptr
  implicit none
  logical, dimension(:,:), pointer, intent(inout) :: array
  include 'deallocate-profile-inc.f90'
  if (.not. associated(array)) return
  include 'deallocate-ptr-inc.f90'
end subroutine l2_ptr_free

subroutine l3_ptr_free(array)
  use metadata_interfaces, metadata_address => getl3ptr
  implicit none
  logical, dimension(:,:,:), pointer, intent(inout) :: array
  include 'deallocate-profile-inc.f90'
  if (.not. associated(array)) return
  include 'deallocate-ptr-inc.f90'
end subroutine l3_ptr_free

subroutine ll1_ptr(array,m)
  use metadata_interfaces, metadata_address => getl2ptr
  implicit none
  type(malloc_information_ptr), intent(in) :: m
  logical(f_byte), dimension(:), pointer, intent(inout) :: array
  include 'allocate-profile-inc.f90'
  include 'allocate-ptr-inc.f90' 
  !allocate the array
  allocate(array(m%lbounds(1):m%ubounds(1)+ndebug),stat=ierror)
  include 'allocate-inc.f90'
end subroutine ll1_ptr

subroutine ll1_ptr_free(array)
  use metadata_interfaces, metadata_address => getl2ptr
  implicit none
  logical(f_byte), dimension(:), pointer, intent(inout) :: array
  include 'deallocate-profile-inc.f90' 
  if (.not. associated(array)) return
  include 'deallocate-ptr-inc.f90' 
end subroutine ll1_ptr_free

subroutine ll2_ptr(array,m)
  use metadata_interfaces, metadata_address => getl2ptr
  implicit none
  type(malloc_information_ptr), intent(in) :: m
  logical(f_byte), dimension(:,:), pointer, intent(inout) :: array
  include 'allocate-profile-inc.f90'
  include 'allocate-ptr-inc.f90' 
  !allocate the array
  allocate(array(m%lbounds(1):m%ubounds(1),m%lbounds(2):m%ubounds(2)+ndebug),stat=ierror)
  include 'allocate-inc.f90'
end subroutine ll2_ptr

subroutine ll3_ptr(array,m)
  use metadata_interfaces, metadata_address => getl3ptr
  implicit none
  type(malloc_information_ptr), intent(in) :: m
  logical(f_byte), dimension(:,:,:), pointer, intent(inout) :: array
  include 'allocate-profile-inc.f90'
  include 'allocate-ptr-inc.f90' 
  !allocate the array
  allocate(array(m%lbounds(1):m%ubounds(1),m%lbounds(2):m%ubounds(2),&
       m%lbounds(3):m%ubounds(3)+ndebug),stat=ierror)
  include 'allocate-inc.f90'
end subroutine ll3_ptr

subroutine ll2_ptr_free(array)
  use metadata_interfaces, metadata_address => getl2ptr
  implicit none
  logical(f_byte), dimension(:,:), pointer, intent(inout) :: array
  include 'deallocate-profile-inc.f90' 
  if (.not. associated(array)) return
  include 'deallocate-ptr-inc.f90' 
end subroutine ll2_ptr_free

subroutine ll3_ptr_free(array)
  use metadata_interfaces, metadata_address => getl3ptr
  implicit none
  logical(f_byte), dimension(:,:,:), pointer, intent(inout) :: array
  include 'deallocate-profile-inc.f90' 
  if (.not. associated(array)) return
  include 'deallocate-ptr-inc.f90' 
end subroutine ll3_ptr_free

!character arrays
subroutine c1_ptr(array,m)
  use metadata_interfaces, metadata_address => getc1ptr
  implicit none
  type(malloc_information_str_ptr), intent(in) :: m
  character(len=m%len), dimension(:), pointer, intent(inout) :: array
  include 'allocate-profile-inc.f90'
  !include 'allocate-ptr-inc.f90'
  if (f_nan_pad_size > 0) then
     padding=f_nan_pad_size
     call togglepadding(product(int(m%shape(1:m%rank-1),f_long))*&
          m%len*kind(array)*(m%shape(m%rank)+padding))
  end if
  if (m%srcdata_add == int(-1,kind=8)) then
     call f_free_str_ptr(m%len,array) !to avoid memory leaks
     call f_timer_resume()!TCAT_ARRAY_ALLOCATIONS
     return
  end if
  !allocate the array
  allocate(array(m%lbounds(1):m%ubounds(1)+ndebug),stat=ierror)
  include 'allocate-c-inc.f90'
  !include 'allocate-inc.f90'
end subroutine c1_ptr

!subroutine c1_ptr_free(length,array)
subroutine f_free_str_ptr(length,array)
  use metadata_interfaces, metadata_address => getc1ptr
  implicit none
  integer, intent(in) :: length
  character(len=length), dimension(:), pointer, intent(inout) :: array
  include 'deallocate-profile-inc.f90'
  if (.not. associated(array)) return
  !include 'deallocate-c-inc.f90'
  include 'deallocate-ptr-inc.f90'
  nullify(array)
end subroutine f_free_str_ptr

subroutine z1_ptr(array,m)
  use metadata_interfaces, metadata_address => getz1ptr
  implicit none
  type(malloc_information_ptr), intent(in) :: m
  double complex, dimension(:), pointer, intent(inout) :: array
  !local variables
  include 'allocate-profile-inc.f90'
  include 'allocate-ptr-inc.f90'
  !allocate the array
  allocate(array(m%lbounds(1):m%ubounds(1)+ndebug),stat=ierror)
  include 'allocate-inc.f90'
end subroutine z1_ptr

subroutine z1_ptr_free(array)
  use metadata_interfaces, metadata_address => getz1ptr
  implicit none
  double complex, dimension(:), pointer, intent(inout) :: array
  include 'deallocate-profile-inc.f90'
  if (.not. associated(array)) return
  include 'deallocate-ptr-inc.f90'
  nullify(array)
end subroutine z1_ptr_free

subroutine z2_ptr(array,m)
  use metadata_interfaces, metadata_address => getz1ptr
  implicit none
  type(malloc_information_ptr), intent(in) :: m
  double complex, dimension(:,:), pointer, intent(inout) :: array
  !local variables
  include 'allocate-profile-inc.f90'
  include 'allocate-ptr-inc.f90'
  !allocate the array
  allocate(array(m%lbounds(1):m%ubounds(1),m%lbounds(2):m%ubounds(2)+ndebug),stat=ierror)
  include 'allocate-inc.f90'
end subroutine z2_ptr

subroutine z2_ptr_free(array)
  use metadata_interfaces, metadata_address => getz1ptr
  implicit none
  double complex, dimension(:,:), pointer, intent(inout) :: array
  include 'deallocate-profile-inc.f90'
  if (.not. associated(array)) return
  include 'deallocate-ptr-inc.f90'
  nullify(array)
end subroutine z2_ptr_free
