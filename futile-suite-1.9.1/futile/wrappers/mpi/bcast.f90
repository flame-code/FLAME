!> @file
!! Wrapper for mpi_bcast
!! Use error handling
!! @author
!!    Copyright (C) 2018-2018 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
module f_bcast
  use f_precisions
  use fmpi_types
  use yaml_strings
  use f_gather
  use f_allreduce
  use dictionaries, only: f_err_throw
  implicit none
  private

  interface fmpi_bcast
     module procedure mpibcast_i0,mpibcast_li0,mpibcast_d0,mpibcast_c0
     module procedure mpibcast_c1,mpibcast_d1,mpibcast_d2,mpibcast_d3,mpibcast_i1,mpibcast_i2, mpibcast_i3
  end interface fmpi_bcast

  interface fmpi_maxdiff
     module procedure mpimaxdiff_i0,mpimaxdiff_li0,mpimaxdiff_d0
     module procedure mpimaxdiff_i1,mpimaxdiff_li1,mpimaxdiff_d1
     module procedure mpimaxdiff_i2,mpimaxdiff_d2
  end interface fmpi_maxdiff


  public :: fmpi_bcast,fmpi_maxdiff

contains
  
  recursive subroutine mpibcast_i0(buffer,count,root,comm,check,maxdiff)
    use dynamic_memory
    use dictionaries, only: f_err_throw
    use yaml_output !for check=.true.
    use f_utils, only: f_zero
    implicit none
    integer(f_integer) ::  buffer
    integer(f_integer), intent(out), optional :: maxdiff
    integer(f_integer), dimension(:), allocatable :: array_diff
    include 'bcast-decl-inc.f90'
    include 'bcast-inc.f90'
  end subroutine mpibcast_i0

  subroutine mpibcast_li0(buffer,count,root,comm,check,maxdiff)
    use dynamic_memory
    use dictionaries, only: f_err_throw
    use yaml_output !for check=.true.
    use f_utils, only: f_zero
    implicit none
    integer(f_long) ::  buffer
    integer(f_long), intent(out), optional :: maxdiff
    integer(f_long), dimension(:), allocatable :: array_diff
    include 'bcast-decl-inc.f90'
    include 'bcast-inc.f90'
  end subroutine mpibcast_li0

  recursive subroutine mpibcast_d0(buffer,count,root,comm,check,maxdiff)
    use dynamic_memory
    use dictionaries, only: f_err_throw
    use yaml_output !for check=.true.
    use f_utils, only: f_zero
    implicit none
    real(f_double) ::  buffer
    real(f_double), intent(out), optional :: maxdiff
    real(f_double), dimension(:), allocatable :: array_diff
    include 'bcast-decl-inc.f90'
    include 'bcast-inc.f90'
  end subroutine mpibcast_d0

  subroutine mpibcast_c0(buffer,root,comm,check,maxdiff)
    use dynamic_memory
    use dictionaries, only: f_err_throw
    use yaml_output !for check=.true.
    use f_utils, only: f_zero
    implicit none
    character(len=*) ::  buffer
    integer(f_integer), intent(out), optional :: maxdiff
    integer(f_integer), dimension(:), allocatable :: array_diff !<the difference is performed with ascii value
    ! 'bcast-decl-arr-inc.f90'
    integer, intent(in), optional :: root  !< @copydoc doc::root
    integer, intent(in), optional :: comm  !< @copydoc doc::comm
    logical, intent(in), optional :: check !< performs the check of the arguments
    !local variables
    logical chk
    integer :: n,iroot,mpi_comm,ierr
    integer, dimension(3) :: iarg_check
    external :: MPI_BCAST

    chk=.false.
    n=len(buffer)
    if (present(maxdiff)) then
       call f_zero(maxdiff)
       array_diff=f_malloc(n,id='array_diff')
       call f_memcpy(src=buffer,dest=array_diff)
    end if
    ! end bcast_decl
    include 'bcast-inc.f90'
  end subroutine mpibcast_c0

  subroutine mpibcast_c1(buffer,root,comm,check,maxdiff)
    use dynamic_memory
    use dictionaries, only: f_err_throw
    use yaml_output !for check=.true.
    use f_utils, only: f_zero
    implicit none
    character, dimension(:), intent(inout) ::  buffer
    integer, intent(out), optional :: maxdiff
    integer, dimension(:), allocatable :: array_diff !<the difference is performed with ascii value
    include 'bcast-decl-arr-inc.f90'
    include 'bcast-inc.f90'
  end subroutine mpibcast_c1

  subroutine mpibcast_i1(buffer,root,comm,check,maxdiff)
    use dynamic_memory
    use dictionaries, only: f_err_throw
    use yaml_output !for check=.true.
    use f_utils, only: f_zero
    implicit none
    integer, dimension(:), intent(inout) ::  buffer
    integer, intent(out), optional :: maxdiff
    integer, dimension(:), allocatable :: array_diff
    include 'bcast-decl-arr-inc.f90'
    include 'bcast-inc.f90'
  end subroutine mpibcast_i1

  subroutine mpibcast_i2(buffer,root,comm,check,maxdiff)
    use dynamic_memory
    use dictionaries, only: f_err_throw
    use yaml_output !for check=.true.
    use f_utils, only: f_zero
    implicit none
    integer, dimension(:,:), intent(inout) ::  buffer
    integer, intent(out), optional :: maxdiff
    integer, dimension(:), allocatable :: array_diff
    include 'bcast-decl-arr-inc.f90'
    include 'bcast-inc.f90'
  end subroutine mpibcast_i2

  subroutine mpibcast_i3(buffer,root,comm,check,maxdiff)
    use dynamic_memory
    use dictionaries, only: f_err_throw
    use yaml_output !for check=.true.
    use f_utils, only: f_zero
    implicit none
    integer, dimension(:,:,:), intent(inout) ::  buffer
    integer, intent(out), optional :: maxdiff
    integer, dimension(:), allocatable :: array_diff
    include 'bcast-decl-arr-inc.f90'
    include 'bcast-inc.f90'
  end subroutine mpibcast_i3

  subroutine mpibcast_d1(buffer,root,comm,check,maxdiff)
    use dynamic_memory
    use dictionaries, only: f_err_throw
    use yaml_output !for check=.true.
    use f_utils, only: f_zero
    implicit none
    real(f_double), dimension(:), intent(inout) ::  buffer
    real(f_double), intent(out), optional :: maxdiff
    real(f_double), dimension(:), allocatable :: array_diff
    include 'bcast-decl-arr-inc.f90'
    include 'bcast-inc.f90'
  end subroutine mpibcast_d1

  subroutine mpibcast_d2(buffer,root,comm,check,maxdiff)
    use dynamic_memory
    use dictionaries, only: f_err_throw
    use yaml_output !for check=.true.
    use f_utils, only: f_zero
    implicit none
    real(f_double), dimension(:,:), intent(inout) ::  buffer
    real(f_double), intent(out), optional :: maxdiff
    real(f_double), dimension(:), allocatable :: array_diff
    include 'bcast-decl-arr-inc.f90'
    include 'bcast-inc.f90'
  end subroutine mpibcast_d2

  subroutine mpibcast_d3(buffer,root,comm,check,maxdiff)
    use dynamic_memory
    use dictionaries, only: f_err_throw
    use yaml_output !for check=.true.
    use f_utils, only: f_zero
    implicit none
    real(f_double), dimension(:,:,:), intent(inout) ::  buffer
    real(f_double), intent(out), optional :: maxdiff
    real(f_double), dimension(:), allocatable :: array_diff
    include 'bcast-decl-arr-inc.f90'
    include 'bcast-inc.f90'
  end subroutine mpibcast_d3

  include 'maxdiff-routines-inc.f90'

end module f_bcast
