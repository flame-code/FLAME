!> @file
!! Define the maxdiff operation via mpi wrappers
!! @author
!!    Copyright (C) 2012-2018 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
!!$module f_maxdiff
!!$  use f_precisions
!!$  use f_allreduce
!!$  implicit none
!!$
!!$  private
!!$
!!$
!!$
!!$contains
  
  !> Detect the maximum difference between arrays all over a given communicator
  function mpimaxdiff_i0(n,array,root,source,comm,bcast) result(maxdiff)
    use dynamic_memory
    implicit none
    integer, intent(in) :: n !<number of elements to be controlled
    integer(f_integer) :: array !< starting point of the array
    integer(f_integer), dimension(:,:), allocatable :: array_glob
    integer(f_integer) :: maxdiff
    include 'maxdiff-decl-inc.f90'

    ndims = n
    maxdiff=0

    include 'maxdiff-inc.f90'
  end function mpimaxdiff_i0

  !> Detect the maximum difference between arrays all over a given communicator
  function mpimaxdiff_li0(n,array,root,source,comm,bcast) result(maxdiff)
    use dynamic_memory
    implicit none
    integer(f_long), intent(in) :: n !<number of elements to be controlled
    integer(f_long) :: array !< starting point of the array
    integer(f_long), dimension(:,:), allocatable :: array_glob
    integer(f_long) :: maxdiff
    include 'maxdiff-decl-inc.f90'
    ndims = int(n,kind=4)
    maxdiff=int(0,f_long)
    include 'maxdiff-inc.f90'
  end function mpimaxdiff_li0

  function mpimaxdiff_c1(array,root,source,comm,bcast) result(maxdiff)
    use dynamic_memory
    implicit none
    !> array to be checked
    character, dimension(:), intent(inout) :: array
    integer, dimension(:,:), allocatable :: array_glob
    integer :: maxdiff
    include 'maxdiff-decl-inc.f90'
    ndims = size(array)
    maxdiff=0
    include 'maxdiff-arr-inc.f90'
  end function mpimaxdiff_c1

  function mpimaxdiff_d0(n,array,root,source,comm,bcast) result(maxdiff)
    use dynamic_memory
    implicit none
    integer, intent(in) :: n !<number of elements to be controlled
    double precision :: array !< starting point of the array
    double precision, dimension(:,:), allocatable :: array_glob
    double precision :: maxdiff
    include 'maxdiff-decl-inc.f90'

    ndims = n
    maxdiff=0.d0

    include 'maxdiff-inc.f90'
  end function mpimaxdiff_d0

  function mpimaxdiff_d1(array,root,source,comm,bcast) result(maxdiff)
    use dynamic_memory
    implicit none
    !> array to be checked
    double precision, dimension(:), intent(inout) :: array
    double precision, dimension(:,:), allocatable :: array_glob
    double precision :: maxdiff
    include 'maxdiff-decl-inc.f90'

    ndims = size(array)

    maxdiff=0.d0

    include 'maxdiff-arr-inc.f90'
  end function mpimaxdiff_d1

  function mpimaxdiff_i1(array,root,source,comm,bcast) result(maxdiff)
    use dynamic_memory
    implicit none
    !> array to be checked
    integer(f_integer), dimension(:), intent(inout) :: array
    integer(f_integer), dimension(:,:), allocatable :: array_glob
    integer(f_integer) :: maxdiff
    include 'maxdiff-decl-inc.f90'

    ndims = size(array)

    maxdiff=0

    include 'maxdiff-arr-inc.f90'
  end function mpimaxdiff_i1

  function mpimaxdiff_li1(array,root,source,comm,bcast) result(maxdiff)
    use dynamic_memory
    implicit none
    !> array to be checked
    integer(f_long), dimension(:), intent(inout) :: array
    integer(f_long), dimension(:,:), allocatable :: array_glob
    integer(f_long) :: maxdiff
    include 'maxdiff-decl-inc.f90'

    ndims = size(array)

    maxdiff=0

    include 'maxdiff-arr-inc.f90'
  end function mpimaxdiff_li1

  function mpimaxdiff_i2(array,root,source,comm,bcast) result(maxdiff)
    use dynamic_memory
    implicit none
    !> array to be checked
    integer(f_integer), dimension(:,:), intent(inout) :: array
    integer(f_integer), dimension(:,:), allocatable :: array_glob
    integer :: maxdiff
    include 'maxdiff-decl-inc.f90'

    ndims = size(array)

    maxdiff=0

    include 'maxdiff-arr-inc.f90'
  end function mpimaxdiff_i2

  function mpimaxdiff_d2(array,root,source,comm,bcast) result(maxdiff)
    use dynamic_memory
    implicit none
    !> array to be checked
    double precision, dimension(:,:), intent(inout) :: array
    double precision, dimension(:,:), allocatable :: array_glob
    double precision :: maxdiff
    include 'maxdiff-decl-inc.f90'

    ndims = size(array)

    maxdiff=0.d0

    include 'maxdiff-arr-inc.f90'
  end function mpimaxdiff_d2

!!$end module f_maxdiff
