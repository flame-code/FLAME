!> @file
!! Include fortran file for allocation templates
!! file included in module dynamic_memory.f90
!! @author
!!    Copyright (C) 2013-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS

subroutine c1_all(array,m)
  use metadata_interfaces, metadata_address => getc1
  implicit none
  type(malloc_information_all), intent(in) :: m
  character(len=*), dimension(:), allocatable, intent(inout) :: array
  include 'allocate-profile-inc.f90' 
  !allocate the array
  allocate(array(m%lbounds(1):m%ubounds(1)+ndebug),stat=ierror)
  include 'allocate-inc.f90'
end subroutine c1_all

subroutine c1_all_free(array)
  use metadata_interfaces, metadata_address => getc1
  implicit none
  character(len=*), dimension(:), allocatable, intent(inout) :: array
  include 'deallocate-profile-inc.f90' 
  include 'deallocate-inc.f90' 
end subroutine c1_all_free

!!$subroutine c1_2_all(array,m)
!!$  use metadata_interfaces, metadata_address => getc1_2
!!$  implicit none
!!$  type(malloc_information_all), intent(in) :: m
!!$  character*2, dimension(:), allocatable, intent(inout) :: array
!!$  include 'allocate-profile-inc.f90' 
!!$  !allocate the array
!!$  allocate(array(m%lbounds(1):m%ubounds(1)+ndebug),stat=ierror)
!!$  include 'allocate-inc.f90'
!!$end subroutine c1_2_all
!!$
!!$subroutine c1_2_all_free(array)
!!$  use metadata_interfaces, metadata_address => getc1_2
!!$  implicit none
!!$  character*2, dimension(:), allocatable, intent(inout) :: array
!!$  include 'deallocate-profile-inc.f90' 
!!$  include 'deallocate-inc.f90' 
!!$end subroutine c1_2_all_free
!!$
!!$subroutine c1_3_all(array,m)
!!$  use metadata_interfaces, metadata_address => getc1_3
!!$  implicit none
!!$  type(malloc_information_all), intent(in) :: m
!!$  character(len=3), dimension(:), allocatable, intent(inout) :: array
!!$  include 'allocate-profile-inc.f90' 
!!$  !allocate the array
!!$  allocate(array(m%lbounds(1):m%ubounds(1)+ndebug),stat=ierror)
!!$  include 'allocate-inc.f90'
!!$end subroutine c1_3_all
!!$
!!$subroutine c1_3_all_free(array)
!!$  use metadata_interfaces, metadata_address => getc1_3
!!$  implicit none
!!$  character(len=3), dimension(:), allocatable, intent(inout) :: array
!!$  include 'deallocate-profile-inc.f90' 
!!$  include 'deallocate-inc.f90' 
!!$end subroutine c1_3_all_free
!!$
!!$subroutine c1_4_all(array,m)
!!$  use metadata_interfaces, metadata_address => getc1_4
!!$  implicit none
!!$  type(malloc_information_all), intent(in) :: m
!!$  character(len=4), dimension(:), allocatable, intent(inout) :: array
!!$  include 'allocate-profile-inc.f90' 
!!$  !allocate the array
!!$  allocate(array(m%lbounds(1):m%ubounds(1)+ndebug),stat=ierror)
!!$  include 'allocate-inc.f90'
!!$end subroutine c1_4_all
!!$
!!$subroutine c1_4_all_free(array)
!!$  use metadata_interfaces, metadata_address => getc1_4
!!$  implicit none
!!$  character(len=4), dimension(:), allocatable, intent(inout) :: array
!!$  include 'deallocate-profile-inc.f90' 
!!$  include 'deallocate-inc.f90' 
!!$end subroutine c1_4_all_free
!!$
!!$subroutine c1_5_all(array,m)
!!$  use metadata_interfaces, metadata_address => getc1_5
!!$  implicit none
!!$  type(malloc_information_all), intent(in) :: m
!!$  character(len=5), dimension(:), allocatable, intent(inout) :: array
!!$  include 'allocate-profile-inc.f90' 
!!$  !allocate the array
!!$  allocate(array(m%lbounds(1):m%ubounds(1)+ndebug),stat=ierror)
!!$  include 'allocate-inc.f90'
!!$end subroutine c1_5_all
!!$
!!$subroutine c1_5_all_free(array)
!!$  use metadata_interfaces, metadata_address => getc1_5
!!$  implicit none
!!$  character(len=5), dimension(:), allocatable, intent(inout) :: array
!!$  include 'deallocate-profile-inc.f90' 
!!$  include 'deallocate-inc.f90' 
!!$end subroutine c1_5_all_free
!!$
!!$subroutine c1_6_all(array,m)
!!$  use metadata_interfaces, metadata_address => getc1_6
!!$  implicit none
!!$  type(malloc_information_all), intent(in) :: m
!!$  character(len=6), dimension(:), allocatable, intent(inout) :: array
!!$  include 'allocate-profile-inc.f90' 
!!$  !allocate the array
!!$  allocate(array(m%lbounds(1):m%ubounds(1)+ndebug),stat=ierror)
!!$  include 'allocate-inc.f90'
!!$end subroutine c1_6_all
!!$
!!$subroutine c1_6_all_free(array)
!!$  use metadata_interfaces, metadata_address => getc1_6
!!$  implicit none
!!$  character(len=6), dimension(:), allocatable, intent(inout) :: array
!!$  include 'deallocate-profile-inc.f90' 
!!$  include 'deallocate-inc.f90' 
!!$end subroutine c1_6_all_free
!!$
!!$subroutine c1_7_all(array,m)
!!$  use metadata_interfaces, metadata_address => getc1_7
!!$  implicit none
!!$  type(malloc_information_all), intent(in) :: m
!!$  character(len=7), dimension(:), allocatable, intent(inout) :: array
!!$  include 'allocate-profile-inc.f90' 
!!$  !allocate the array
!!$  allocate(array(m%lbounds(1):m%ubounds(1)+ndebug),stat=ierror)
!!$  include 'allocate-inc.f90'
!!$end subroutine c1_7_all
!!$
!!$subroutine c1_7_all_free(array)
!!$  use metadata_interfaces, metadata_address => getc1_7
!!$  implicit none
!!$  character(len=7), dimension(:), allocatable, intent(inout) :: array
!!$  include 'deallocate-profile-inc.f90' 
!!$  include 'deallocate-inc.f90' 
!!$end subroutine c1_7_all_free
!!$
!!$subroutine c1_8_all(array,m)
!!$  use metadata_interfaces, metadata_address => getc1_8
!!$  implicit none
!!$  type(malloc_information_all), intent(in) :: m
!!$  character(len=8), dimension(:), allocatable, intent(inout) :: array
!!$  include 'allocate-profile-inc.f90' 
!!$  !allocate the array
!!$  allocate(array(m%lbounds(1):m%ubounds(1)+ndebug),stat=ierror)
!!$  include 'allocate-inc.f90'
!!$end subroutine c1_8_all
!!$
!!$subroutine c1_8_all_free(array)
!!$  use metadata_interfaces, metadata_address => getc1_8
!!$  implicit none
!!$  character(len=8), dimension(:), allocatable, intent(inout) :: array
!!$  include 'deallocate-profile-inc.f90' 
!!$  include 'deallocate-inc.f90' 
!!$end subroutine c1_8_all_free
!!$
!!$subroutine c1_9_all(array,m)
!!$  use metadata_interfaces, metadata_address => getc1_9
!!$  implicit none
!!$  type(malloc_information_all), intent(in) :: m
!!$  character(len=9), dimension(:), allocatable, intent(inout) :: array
!!$  include 'allocate-profile-inc.f90' 
!!$  !allocate the array
!!$  allocate(array(m%lbounds(1):m%ubounds(1)+ndebug),stat=ierror)
!!$  include 'allocate-inc.f90'
!!$end subroutine c1_9_all
!!$
!!$subroutine c1_9_all_free(array)
!!$  use metadata_interfaces, metadata_address => getc1_9
!!$  implicit none
!!$  character(len=9), dimension(:), allocatable, intent(inout) :: array
!!$  include 'deallocate-profile-inc.f90' 
!!$  include 'deallocate-inc.f90' 
!!$end subroutine c1_9_all_free
!!$
!!$subroutine c1_10_all(array,m)
!!$  use metadata_interfaces, metadata_address => getc1_10
!!$  implicit none
!!$  type(malloc_information_all), intent(in) :: m
!!$  character(len=10), dimension(:), allocatable, intent(inout) :: array
!!$  include 'allocate-profile-inc.f90' 
!!$  !allocate the array
!!$  allocate(array(m%lbounds(1):m%ubounds(1)+ndebug),stat=ierror)
!!$  include 'allocate-inc.f90'
!!$end subroutine c1_10_all
!!$
!!$subroutine c1_10_all_free(array)
!!$  use metadata_interfaces, metadata_address => getc1_10
!!$  implicit none
!!$  character(len=10), dimension(:), allocatable, intent(inout) :: array
!!$  include 'deallocate-profile-inc.f90' 
!!$  include 'deallocate-inc.f90' 
!!$end subroutine c1_10_all_free
!!$
!!$subroutine c1_11_all(array,m)
!!$  use metadata_interfaces, metadata_address => getc1_11
!!$  implicit none
!!$  type(malloc_information_all), intent(in) :: m
!!$  character(len=11), dimension(:), allocatable, intent(inout) :: array
!!$  include 'allocate-profile-inc.f90' 
!!$  !allocate the array
!!$  allocate(array(m%lbounds(1):m%ubounds(1)+ndebug),stat=ierror)
!!$  include 'allocate-inc.f90'
!!$end subroutine c1_11_all
!!$
!!$subroutine c1_11_all_free(array)
!!$  use metadata_interfaces, metadata_address => getc1_11
!!$  implicit none
!!$  character(len=11), dimension(:), allocatable, intent(inout) :: array
!!$  include 'deallocate-profile-inc.f90' 
!!$  include 'deallocate-inc.f90' 
!!$end subroutine c1_11_all_free
!!$
!!$subroutine c1_12_all(array,m)
!!$  use metadata_interfaces, metadata_address => getc1_12
!!$  implicit none
!!$  type(malloc_information_all), intent(in) :: m
!!$  character(len=12), dimension(:), allocatable, intent(inout) :: array
!!$  include 'allocate-profile-inc.f90' 
!!$  !allocate the array
!!$  allocate(array(m%lbounds(1):m%ubounds(1)+ndebug),stat=ierror)
!!$  include 'allocate-inc.f90'
!!$end subroutine c1_12_all
!!$
!!$subroutine c1_12_all_free(array)
!!$  use metadata_interfaces, metadata_address => getc1_12
!!$  implicit none
!!$  character(len=12), dimension(:), allocatable, intent(inout) :: array
!!$  include 'deallocate-profile-inc.f90' 
!!$  include 'deallocate-inc.f90' 
!!$end subroutine c1_12_all_free
!!$
!!$subroutine c1_13_all(array,m)
!!$  use metadata_interfaces, metadata_address => getc1_13
!!$  implicit none
!!$  type(malloc_information_all), intent(in) :: m
!!$  character(len=13), dimension(:), allocatable, intent(inout) :: array
!!$  include 'allocate-profile-inc.f90' 
!!$  !allocate the array
!!$  allocate(array(m%lbounds(1):m%ubounds(1)+ndebug),stat=ierror)
!!$  include 'allocate-inc.f90'
!!$end subroutine c1_13_all
!!$
!!$subroutine c1_13_all_free(array)
!!$  use metadata_interfaces, metadata_address => getc1_13
!!$  implicit none
!!$  character(len=13), dimension(:), allocatable, intent(inout) :: array
!!$  include 'deallocate-profile-inc.f90' 
!!$  include 'deallocate-inc.f90' 
!!$end subroutine c1_13_all_free
!!$
!!$subroutine c1_14_all(array,m)
!!$  use metadata_interfaces, metadata_address => getc1_14
!!$  implicit none
!!$  type(malloc_information_all), intent(in) :: m
!!$  character(len=14), dimension(:), allocatable, intent(inout) :: array
!!$  include 'allocate-profile-inc.f90' 
!!$  !allocate the array
!!$  allocate(array(m%lbounds(1):m%ubounds(1)+ndebug),stat=ierror)
!!$  include 'allocate-inc.f90'
!!$end subroutine c1_14_all
!!$
!!$subroutine c1_14_all_free(array)
!!$  use metadata_interfaces, metadata_address => getc1_14
!!$  implicit none
!!$  character(len=14), dimension(:), allocatable, intent(inout) :: array
!!$  include 'deallocate-profile-inc.f90' 
!!$  include 'deallocate-inc.f90' 
!!$end subroutine c1_14_all_free
!!$
!!$subroutine c1_15_all(array,m)
!!$  use metadata_interfaces, metadata_address => getc1_15
!!$  implicit none
!!$  type(malloc_information_all), intent(in) :: m
!!$  character(len=15), dimension(:), allocatable, intent(inout) :: array
!!$  include 'allocate-profile-inc.f90' 
!!$  !allocate the array
!!$  allocate(array(m%lbounds(1):m%ubounds(1)+ndebug),stat=ierror)
!!$  include 'allocate-inc.f90'
!!$end subroutine c1_15_all
!!$
!!$subroutine c1_15_all_free(array)
!!$  use metadata_interfaces, metadata_address => getc1_15
!!$  implicit none
!!$  character(len=15), dimension(:), allocatable, intent(inout) :: array
!!$  include 'deallocate-profile-inc.f90' 
!!$  include 'deallocate-inc.f90' 
!!$end subroutine c1_15_all_free
!!$
!!$subroutine c1_16_all(array,m)
!!$  use metadata_interfaces, metadata_address => getc1_16
!!$  implicit none
!!$  type(malloc_information_all), intent(in) :: m
!!$  character(len=16), dimension(:), allocatable, intent(inout) :: array
!!$  include 'allocate-profile-inc.f90' 
!!$  !allocate the array
!!$  allocate(array(m%lbounds(1):m%ubounds(1)+ndebug),stat=ierror)
!!$  include 'allocate-inc.f90'
!!$end subroutine c1_16_all
!!$
!!$subroutine c1_16_all_free(array)
!!$  use metadata_interfaces, metadata_address => getc1_16
!!$  implicit none
!!$  character(len=16), dimension(:), allocatable, intent(inout) :: array
!!$  include 'deallocate-profile-inc.f90' 
!!$  include 'deallocate-inc.f90' 
!!$end subroutine c1_16_all_free
!!$
!!$subroutine c1_17_all(array,m)
!!$  use metadata_interfaces, metadata_address => getc1_17
!!$  implicit none
!!$  type(malloc_information_all), intent(in) :: m
!!$  character(len=17), dimension(:), allocatable, intent(inout) :: array
!!$  include 'allocate-profile-inc.f90' 
!!$  !allocate the array
!!$  allocate(array(m%lbounds(1):m%ubounds(1)+ndebug),stat=ierror)
!!$  include 'allocate-inc.f90'
!!$end subroutine c1_17_all
!!$
!!$subroutine c1_17_all_free(array)
!!$  use metadata_interfaces, metadata_address => getc1_17
!!$  implicit none
!!$  character(len=17), dimension(:), allocatable, intent(inout) :: array
!!$  include 'deallocate-profile-inc.f90' 
!!$  include 'deallocate-inc.f90' 
!!$end subroutine c1_17_all_free
!!$
!!$subroutine c1_18_all(array,m)
!!$  use metadata_interfaces, metadata_address => getc1_18
!!$  implicit none
!!$  type(malloc_information_all), intent(in) :: m
!!$  character(len=18), dimension(:), allocatable, intent(inout) :: array
!!$  include 'allocate-profile-inc.f90' 
!!$  !allocate the array
!!$  allocate(array(m%lbounds(1):m%ubounds(1)+ndebug),stat=ierror)
!!$  include 'allocate-inc.f90'
!!$end subroutine c1_18_all
!!$
!!$subroutine c1_18_all_free(array)
!!$  use metadata_interfaces, metadata_address => getc1_18
!!$  implicit none
!!$  character(len=18), dimension(:), allocatable, intent(inout) :: array
!!$  include 'deallocate-profile-inc.f90' 
!!$  include 'deallocate-inc.f90' 
!!$end subroutine c1_18_all_free
!!$
!!$subroutine c1_19_all(array,m)
!!$  use metadata_interfaces, metadata_address => getc1_19
!!$  implicit none
!!$  type(malloc_information_all), intent(in) :: m
!!$  character(len=19), dimension(:), allocatable, intent(inout) :: array
!!$  include 'allocate-profile-inc.f90' 
!!$  !allocate the array
!!$  allocate(array(m%lbounds(1):m%ubounds(1)+ndebug),stat=ierror)
!!$  include 'allocate-inc.f90'
!!$end subroutine c1_19_all
!!$
!!$subroutine c1_19_all_free(array)
!!$  use metadata_interfaces, metadata_address => getc1_19
!!$  implicit none
!!$  character(len=19), dimension(:), allocatable, intent(inout) :: array
!!$  include 'deallocate-profile-inc.f90' 
!!$  include 'deallocate-inc.f90' 
!!$end subroutine c1_19_all_free
!!$
!!$subroutine c1_20_all(array,m)
!!$  use metadata_interfaces, metadata_address => getc1_20
!!$  implicit none
!!$  type(malloc_information_all), intent(in) :: m
!!$  character(len=20), dimension(:), allocatable, intent(inout) :: array
!!$  include 'allocate-profile-inc.f90' 
!!$  !allocate the array
!!$  allocate(array(m%lbounds(1):m%ubounds(1)+ndebug),stat=ierror)
!!$  include 'allocate-inc.f90'
!!$end subroutine c1_20_all
!!$
!!$subroutine c1_20_all_free(array)
!!$  use metadata_interfaces, metadata_address => getc1_20
!!$  implicit none
!!$  character(len=20), dimension(:), allocatable, intent(inout) :: array
!!$  include 'deallocate-profile-inc.f90' 
!!$  include 'deallocate-inc.f90' 
!!$end subroutine c1_20_all_free
!!$
!!$subroutine c1_21_all(array,m)
!!$  use metadata_interfaces, metadata_address => getc1_21
!!$  implicit none
!!$  type(malloc_information_all), intent(in) :: m
!!$  character(len=21), dimension(:), allocatable, intent(inout) :: array
!!$  include 'allocate-profile-inc.f90' 
!!$  !allocate the array
!!$  allocate(array(m%lbounds(1):m%ubounds(1)+ndebug),stat=ierror)
!!$  include 'allocate-inc.f90'
!!$end subroutine c1_21_all
!!$
!!$subroutine c1_21_all_free(array)
!!$  use metadata_interfaces, metadata_address => getc1_21
!!$  implicit none
!!$  character(len=21), dimension(:), allocatable, intent(inout) :: array
!!$  include 'deallocate-profile-inc.f90' 
!!$  include 'deallocate-inc.f90' 
!!$end subroutine c1_21_all_free
!!$
!!$subroutine c1_22_all(array,m)
!!$  use metadata_interfaces, metadata_address => getc1_22
!!$  implicit none
!!$  type(malloc_information_all), intent(in) :: m
!!$  character(len=22), dimension(:), allocatable, intent(inout) :: array
!!$  include 'allocate-profile-inc.f90' 
!!$  !allocate the array
!!$  allocate(array(m%lbounds(1):m%ubounds(1)+ndebug),stat=ierror)
!!$  include 'allocate-inc.f90'
!!$end subroutine c1_22_all
!!$
!!$subroutine c1_22_all_free(array)
!!$  use metadata_interfaces, metadata_address => getc1_22
!!$  implicit none
!!$  character(len=22), dimension(:), allocatable, intent(inout) :: array
!!$  include 'deallocate-profile-inc.f90' 
!!$  include 'deallocate-inc.f90' 
!!$end subroutine c1_22_all_free
!!$
!!$subroutine c1_23_all(array,m)
!!$  use metadata_interfaces, metadata_address => getc1_23
!!$  implicit none
!!$  type(malloc_information_all), intent(in) :: m
!!$  character(len=23), dimension(:), allocatable, intent(inout) :: array
!!$  include 'allocate-profile-inc.f90' 
!!$  !allocate the array
!!$  allocate(array(m%lbounds(1):m%ubounds(1)+ndebug),stat=ierror)
!!$  include 'allocate-inc.f90'
!!$end subroutine c1_23_all
!!$
!!$subroutine c1_23_all_free(array)
!!$  use metadata_interfaces, metadata_address => getc1_23
!!$  implicit none
!!$  character(len=23), dimension(:), allocatable, intent(inout) :: array
!!$  include 'deallocate-profile-inc.f90' 
!!$  include 'deallocate-inc.f90' 
!!$end subroutine c1_23_all_free
!!$
!!$subroutine c1_24_all(array,m)
!!$  use metadata_interfaces, metadata_address => getc1_24
!!$  implicit none
!!$  type(malloc_information_all), intent(in) :: m
!!$  character(len=24), dimension(:), allocatable, intent(inout) :: array
!!$  include 'allocate-profile-inc.f90' 
!!$  !allocate the array
!!$  allocate(array(m%lbounds(1):m%ubounds(1)+ndebug),stat=ierror)
!!$  include 'allocate-inc.f90'
!!$end subroutine c1_24_all
!!$
!!$subroutine c1_24_all_free(array)
!!$  use metadata_interfaces, metadata_address => getc1_24
!!$  implicit none
!!$  character(len=24), dimension(:), allocatable, intent(inout) :: array
!!$  include 'deallocate-profile-inc.f90' 
!!$  include 'deallocate-inc.f90' 
!!$end subroutine c1_24_all_free
!!$
!!$subroutine c1_25_all(array,m)
!!$  use metadata_interfaces, metadata_address => getc1_25
!!$  implicit none
!!$  type(malloc_information_all), intent(in) :: m
!!$  character(len=25), dimension(:), allocatable, intent(inout) :: array
!!$  include 'allocate-profile-inc.f90' 
!!$  !allocate the array
!!$  allocate(array(m%lbounds(1):m%ubounds(1)+ndebug),stat=ierror)
!!$  include 'allocate-inc.f90'
!!$end subroutine c1_25_all
!!$
!!$subroutine c1_25_all_free(array)
!!$  use metadata_interfaces, metadata_address => getc1_25
!!$  implicit none
!!$  character(len=25), dimension(:), allocatable, intent(inout) :: array
!!$  include 'deallocate-profile-inc.f90' 
!!$  include 'deallocate-inc.f90' 
!!$end subroutine c1_25_all_free
!!$
!!$subroutine c1_26_all(array,m)
!!$  use metadata_interfaces, metadata_address => getc1_26
!!$  implicit none
!!$  type(malloc_information_all), intent(in) :: m
!!$  character(len=26), dimension(:), allocatable, intent(inout) :: array
!!$  include 'allocate-profile-inc.f90' 
!!$  !allocate the array
!!$  allocate(array(m%lbounds(1):m%ubounds(1)+ndebug),stat=ierror)
!!$  include 'allocate-inc.f90'
!!$end subroutine c1_26_all
!!$
!!$subroutine c1_26_all_free(array)
!!$  use metadata_interfaces, metadata_address => getc1_26
!!$  implicit none
!!$  character(len=26), dimension(:), allocatable, intent(inout) :: array
!!$  include 'deallocate-profile-inc.f90' 
!!$  include 'deallocate-inc.f90' 
!!$end subroutine c1_26_all_free
!!$
!!$subroutine c1_27_all(array,m)
!!$  use metadata_interfaces, metadata_address => getc1_27
!!$  implicit none
!!$  type(malloc_information_all), intent(in) :: m
!!$  character(len=27), dimension(:), allocatable, intent(inout) :: array
!!$  include 'allocate-profile-inc.f90' 
!!$  !allocate the array
!!$  allocate(array(m%lbounds(1):m%ubounds(1)+ndebug),stat=ierror)
!!$  include 'allocate-inc.f90'
!!$end subroutine c1_27_all
!!$
!!$subroutine c1_27_all_free(array)
!!$  use metadata_interfaces, metadata_address => getc1_27
!!$  implicit none
!!$  character(len=27), dimension(:), allocatable, intent(inout) :: array
!!$  include 'deallocate-profile-inc.f90' 
!!$  include 'deallocate-inc.f90' 
!!$end subroutine c1_27_all_free
!!$
!!$subroutine c1_28_all(array,m)
!!$  use metadata_interfaces, metadata_address => getc1_28
!!$  implicit none
!!$  type(malloc_information_all), intent(in) :: m
!!$  character(len=28), dimension(:), allocatable, intent(inout) :: array
!!$  include 'allocate-profile-inc.f90' 
!!$  !allocate the array
!!$  allocate(array(m%lbounds(1):m%ubounds(1)+ndebug),stat=ierror)
!!$  include 'allocate-inc.f90'
!!$end subroutine c1_28_all
!!$
!!$subroutine c1_28_all_free(array)
!!$  use metadata_interfaces, metadata_address => getc1_28
!!$  implicit none
!!$  character(len=28), dimension(:), allocatable, intent(inout) :: array
!!$  include 'deallocate-profile-inc.f90' 
!!$  include 'deallocate-inc.f90' 
!!$end subroutine c1_28_all_free
!!$
!!$subroutine c1_29_all(array,m)
!!$  use metadata_interfaces, metadata_address => getc1_29
!!$  implicit none
!!$  type(malloc_information_all), intent(in) :: m
!!$  character(len=29), dimension(:), allocatable, intent(inout) :: array
!!$  include 'allocate-profile-inc.f90' 
!!$  !allocate the array
!!$  allocate(array(m%lbounds(1):m%ubounds(1)+ndebug),stat=ierror)
!!$  include 'allocate-inc.f90'
!!$end subroutine c1_29_all
!!$
!!$subroutine c1_29_all_free(array)
!!$  use metadata_interfaces, metadata_address => getc1_29
!!$  implicit none
!!$  character(len=29), dimension(:), allocatable, intent(inout) :: array
!!$  include 'deallocate-profile-inc.f90' 
!!$  include 'deallocate-inc.f90' 
!!$end subroutine c1_29_all_free
!!$
!!$subroutine c1_30_all(array,m)
!!$  use metadata_interfaces, metadata_address => getc1_30
!!$  implicit none
!!$  type(malloc_information_all), intent(in) :: m
!!$  character(len=30), dimension(:), allocatable, intent(inout) :: array
!!$  include 'allocate-profile-inc.f90' 
!!$  !allocate the array
!!$  allocate(array(m%lbounds(1):m%ubounds(1)+ndebug),stat=ierror)
!!$  include 'allocate-inc.f90'
!!$end subroutine c1_30_all
!!$
!!$subroutine c1_30_all_free(array)
!!$  use metadata_interfaces, metadata_address => getc1_30
!!$  implicit none
!!$  character(len=30), dimension(:), allocatable, intent(inout) :: array
!!$  include 'deallocate-profile-inc.f90' 
!!$  include 'deallocate-inc.f90' 
!!$end subroutine c1_30_all_free
!!$
!!$subroutine c1_31_all(array,m)
!!$  use metadata_interfaces, metadata_address => getc1_31
!!$  implicit none
!!$  type(malloc_information_all), intent(in) :: m
!!$  character(len=31), dimension(:), allocatable, intent(inout) :: array
!!$  include 'allocate-profile-inc.f90' 
!!$  !allocate the array
!!$  allocate(array(m%lbounds(1):m%ubounds(1)+ndebug),stat=ierror)
!!$  include 'allocate-inc.f90'
!!$end subroutine c1_31_all
!!$
!!$subroutine c1_31_all_free(array)
!!$  use metadata_interfaces, metadata_address => getc1_31
!!$  implicit none
!!$  character(len=31), dimension(:), allocatable, intent(inout) :: array
!!$  include 'deallocate-profile-inc.f90' 
!!$  include 'deallocate-inc.f90' 
!!$end subroutine c1_31_all_free
!!$
!!$subroutine c1_32_all(array,m)
!!$  use metadata_interfaces, metadata_address => getc1_32
!!$  implicit none
!!$  type(malloc_information_all), intent(in) :: m
!!$  character(len=32), dimension(:), allocatable, intent(inout) :: array
!!$  include 'allocate-profile-inc.f90' 
!!$  !allocate the array
!!$  allocate(array(m%lbounds(1):m%ubounds(1)+ndebug),stat=ierror)
!!$  include 'allocate-inc.f90'
!!$end subroutine c1_32_all
!!$
!!$subroutine c1_32_all_free(array)
!!$  use metadata_interfaces, metadata_address => getc1_32
!!$  implicit none
!!$  character(len=32), dimension(:), allocatable, intent(inout) :: array
!!$  include 'deallocate-profile-inc.f90' 
!!$  include 'deallocate-inc.f90' 
!!$end subroutine c1_32_all_free
!!$
!!$subroutine c1_33_all(array,m)
!!$  use metadata_interfaces, metadata_address => getc1_33
!!$  implicit none
!!$  type(malloc_information_all), intent(in) :: m
!!$  character(len=33), dimension(:), allocatable, intent(inout) :: array
!!$  include 'allocate-profile-inc.f90' 
!!$  !allocate the array
!!$  allocate(array(m%lbounds(1):m%ubounds(1)+ndebug),stat=ierror)
!!$  include 'allocate-inc.f90'
!!$end subroutine c1_33_all
!!$
!!$subroutine c1_33_all_free(array)
!!$  use metadata_interfaces, metadata_address => getc1_33
!!$  implicit none
!!$  character(len=33), dimension(:), allocatable, intent(inout) :: array
!!$  include 'deallocate-profile-inc.f90' 
!!$  include 'deallocate-inc.f90' 
!!$end subroutine c1_33_all_free
!!$
!!$subroutine c1_34_all(array,m)
!!$  use metadata_interfaces, metadata_address => getc1_34
!!$  implicit none
!!$  type(malloc_information_all), intent(in) :: m
!!$  character(len=34), dimension(:), allocatable, intent(inout) :: array
!!$  include 'allocate-profile-inc.f90' 
!!$  !allocate the array
!!$  allocate(array(m%lbounds(1):m%ubounds(1)+ndebug),stat=ierror)
!!$  include 'allocate-inc.f90'
!!$end subroutine c1_34_all
!!$
!!$subroutine c1_34_all_free(array)
!!$  use metadata_interfaces, metadata_address => getc1_34
!!$  implicit none
!!$  character(len=34), dimension(:), allocatable, intent(inout) :: array
!!$  include 'deallocate-profile-inc.f90' 
!!$  include 'deallocate-inc.f90' 
!!$end subroutine c1_34_all_free
!!$
!!$subroutine c1_35_all(array,m)
!!$  use metadata_interfaces, metadata_address => getc1_35
!!$  implicit none
!!$  type(malloc_information_all), intent(in) :: m
!!$  character(len=35), dimension(:), allocatable, intent(inout) :: array
!!$  include 'allocate-profile-inc.f90' 
!!$  !allocate the array
!!$  allocate(array(m%lbounds(1):m%ubounds(1)+ndebug),stat=ierror)
!!$  include 'allocate-inc.f90'
!!$end subroutine c1_35_all
!!$
!!$subroutine c1_35_all_free(array)
!!$  use metadata_interfaces, metadata_address => getc1_35
!!$  implicit none
!!$  character(len=35), dimension(:), allocatable, intent(inout) :: array
!!$  include 'deallocate-profile-inc.f90' 
!!$  include 'deallocate-inc.f90' 
!!$end subroutine c1_35_all_free
!!$
!!$subroutine c1_36_all(array,m)
!!$  use metadata_interfaces, metadata_address => getc1_36
!!$  implicit none
!!$  type(malloc_information_all), intent(in) :: m
!!$  character(len=36), dimension(:), allocatable, intent(inout) :: array
!!$  include 'allocate-profile-inc.f90' 
!!$  !allocate the array
!!$  allocate(array(m%lbounds(1):m%ubounds(1)+ndebug),stat=ierror)
!!$  include 'allocate-inc.f90'
!!$end subroutine c1_36_all
!!$
!!$subroutine c1_36_all_free(array)
!!$  use metadata_interfaces, metadata_address => getc1_36
!!$  implicit none
!!$  character(len=36), dimension(:), allocatable, intent(inout) :: array
!!$  include 'deallocate-profile-inc.f90' 
!!$  include 'deallocate-inc.f90' 
!!$end subroutine c1_36_all_free
!!$
!!$subroutine c1_37_all(array,m)
!!$  use metadata_interfaces, metadata_address => getc1_37
!!$  implicit none
!!$  type(malloc_information_all), intent(in) :: m
!!$  character(len=37), dimension(:), allocatable, intent(inout) :: array
!!$  include 'allocate-profile-inc.f90' 
!!$  !allocate the array
!!$  allocate(array(m%lbounds(1):m%ubounds(1)+ndebug),stat=ierror)
!!$  include 'allocate-inc.f90'
!!$end subroutine c1_37_all
!!$
!!$subroutine c1_37_all_free(array)
!!$  use metadata_interfaces, metadata_address => getc1_37
!!$  implicit none
!!$  character(len=37), dimension(:), allocatable, intent(inout) :: array
!!$  include 'deallocate-profile-inc.f90' 
!!$  include 'deallocate-inc.f90' 
!!$end subroutine c1_37_all_free
!!$
!!$subroutine c1_38_all(array,m)
!!$  use metadata_interfaces, metadata_address => getc1_38
!!$  implicit none
!!$  type(malloc_information_all), intent(in) :: m
!!$  character(len=38), dimension(:), allocatable, intent(inout) :: array
!!$  include 'allocate-profile-inc.f90' 
!!$  !allocate the array
!!$  allocate(array(m%lbounds(1):m%ubounds(1)+ndebug),stat=ierror)
!!$  include 'allocate-inc.f90'
!!$end subroutine c1_38_all
!!$
!!$subroutine c1_38_all_free(array)
!!$  use metadata_interfaces, metadata_address => getc1_38
!!$  implicit none
!!$  character(len=38), dimension(:), allocatable, intent(inout) :: array
!!$  include 'deallocate-profile-inc.f90' 
!!$  include 'deallocate-inc.f90' 
!!$end subroutine c1_38_all_free
!!$
!!$subroutine c1_39_all(array,m)
!!$  use metadata_interfaces, metadata_address => getc1_39
!!$  implicit none
!!$  type(malloc_information_all), intent(in) :: m
!!$  character(len=39), dimension(:), allocatable, intent(inout) :: array
!!$  include 'allocate-profile-inc.f90' 
!!$  !allocate the array
!!$  allocate(array(m%lbounds(1):m%ubounds(1)+ndebug),stat=ierror)
!!$  include 'allocate-inc.f90'
!!$end subroutine c1_39_all
!!$
!!$subroutine c1_39_all_free(array)
!!$  use metadata_interfaces, metadata_address => getc1_39
!!$  implicit none
!!$  character(len=39), dimension(:), allocatable, intent(inout) :: array
!!$  include 'deallocate-profile-inc.f90' 
!!$  include 'deallocate-inc.f90' 
!!$end subroutine c1_39_all_free
!!$
!!$subroutine c1_40_all(array,m)
!!$  use metadata_interfaces, metadata_address => getc1_40
!!$  implicit none
!!$  type(malloc_information_all), intent(in) :: m
!!$  character(len=40), dimension(:), allocatable, intent(inout) :: array
!!$  include 'allocate-profile-inc.f90' 
!!$  !allocate the array
!!$  allocate(array(m%lbounds(1):m%ubounds(1)+ndebug),stat=ierror)
!!$  include 'allocate-inc.f90'
!!$end subroutine c1_40_all
!!$
!!$subroutine c1_40_all_free(array)
!!$  use metadata_interfaces, metadata_address => getc1_40
!!$  implicit none
!!$  character(len=40), dimension(:), allocatable, intent(inout) :: array
!!$  include 'deallocate-profile-inc.f90' 
!!$  include 'deallocate-inc.f90' 
!!$end subroutine c1_40_all_free
!!$
