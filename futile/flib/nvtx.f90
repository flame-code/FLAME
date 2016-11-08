!> @file
!! Fortran bindings for NVTX
!! @author
!!    Copyright (C) 2016-2016 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
module nvtx
  use iso_c_binding
  public :: nvtxrangepusha, nvtxrangepop
  interface

  subroutine nvtxrangepusha(cat) bind(C, name="nvtxRangePushA")
   use iso_c_binding , only : c_char
   character(kind=c_char) :: cat(*)
  end subroutine nvtxrangepusha

  subroutine nvtxrangepop() bind(C, name="nvtxRangePop")
   end subroutine

  subroutine nvtxmarka(cat) bind(C, name="nvtxMarkA")
   use iso_c_binding , only : c_char
   character(kind=c_char) :: cat(*)
  end subroutine

  end interface
end module nvtx
