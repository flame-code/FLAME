! Fortran bindings for NVTX
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
