!> fake module to substitute the nvidia profiler module
module nvtx
  implicit none
  contains 

    pure subroutine nvtxrangepusha(str)
      implicit none
      character(len=*), intent(in) :: str
    end subroutine nvtxrangepusha

    pure subroutine nvtxrangepop()
      implicit none
    end subroutine nvtxrangepop
    
end module nvtx
