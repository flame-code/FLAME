!> Module used in the builtin_rand function (see razero.f90)
module randomData
  implicit none

  integer, parameter :: ntab=32

  logical :: start = .true.
  integer :: iy = 0
  integer, dimension(NTAB) :: iv
end module randomData


!> Random Number generator from Numerical Recipes
!! To be used for reproducibility of the results
module random
  implicit none

  private

  public :: builtin_rand

  contains

    function builtin_rand(idum, reset)
      use randomData, only : ntab, iy, iv, start
    
      implicit none
    
      integer, intent(inout) :: idum
      real(kind=4) :: builtin_rand
      logical,intent(in),optional :: reset
      !local variables
      integer, parameter :: ia=16807,im=2147483647,iq=127773,ir=2836,ndiv=1+(im-1)/ntab
      real(kind=4), parameter :: am=1.e0/im,eps=1.2e-7,rnmx=1.-eps
      integer :: j,k
      logical :: reset_
    
      reset_ = .false.
      if (present(reset)) reset_ = reset
    
      if (reset_) start = .true.
    
      if (start) then
         iv(:) = 0
         start = .false.
      end if
      if (idum <= 0.or. iy == 0) then
         idum=max(-idum,1)
         do j=ntab+8,1,-1
            k=idum/iq
            idum=ia*(idum-k*iq)-ir*k
            if (idum < 0) idum=idum+im
            if (j <= ntab) iv(j)=idum
         end do
         iy=iv(1)
      endif
      k=idum/iq
      idum=ia*(idum-k*iq)-ir*k
      if (idum <= 0) idum=idum+im
      j=1+iy/ndiv
      iy=iv(j)
      iv(j)=idum
      builtin_rand=min(am*iy,rnmx)
    END FUNCTION builtin_rand

end module random
