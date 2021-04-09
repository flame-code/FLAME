!> Random Number generator from Numerical Recipes
!! To be used for reproducibility of the results
module f_random
  use f_precisions
  implicit none

  private

  integer, parameter :: ntab=32

  logical :: start = .true.
  integer :: iy = 0
  integer :: idum = 0
  integer, dimension(NTAB), save :: iv

!  public :: builtin_rand 
  interface f_random_number
     module procedure f_random_number_i0,f_random_number_d0,f_random_number_r0,f_random_number_d3
     module procedure f_random_number_i1,f_random_number_d1,f_random_number_i1_i1
     !module procedure f_random_number_d4
  end interface f_random_number

  public :: f_random_number,f_random_seed

  contains

    subroutine f_random_seed(seed)
      implicit none
      integer, intent(in) :: seed
      idum=seed
    end subroutine f_random_seed

    subroutine f_random_number_d0(harvest,seed,reset)
      implicit none
      real(f_double), intent(out) :: harvest
      integer, intent(in), optional :: seed
      logical, intent(in), optional :: reset

      if (present(seed)) idum=seed
      harvest = real(builtin_rand(idum,reset=reset),f_double)
      
    end subroutine f_random_number_d0

    subroutine f_random_number_r0(harvest,seed,reset)
      implicit none
      real(f_simple), intent(out) :: harvest
      integer, intent(in), optional :: seed
      logical, intent(in), optional :: reset

      if (present(seed)) idum=seed
      harvest = builtin_rand(idum,reset=reset)

    end subroutine f_random_number_r0

    subroutine f_random_number_d1(harvest,seed,reset)
      implicit none
      real(f_double), dimension(:), intent(out) :: harvest
      integer, intent(in), optional :: seed
      logical, intent(in), optional :: reset
      !local variables
      integer :: i

      if (present(seed)) idum=seed
      
      do i=1,size(harvest)
         harvest(i) = real(builtin_rand(idum,reset=reset),f_double)
      end do
    end subroutine f_random_number_d1

    subroutine f_random_number_d3(harvest,seed,reset)
      implicit none
      real(f_double), dimension(:,:,:), intent(out) :: harvest
      integer, intent(in), optional :: seed
      logical, intent(in), optional :: reset
      !local variables
      integer :: i1,i2,i3

      if (present(seed)) idum=seed
      
      do i3=lbound(harvest,3),ubound(harvest,3)
         do i2=lbound(harvest,2),ubound(harvest,2)
            do i1=lbound(harvest,1),ubound(harvest,1)
               harvest(i1,i2,i3) = real(builtin_rand(idum,reset=reset),f_double)
            end do
         end do
      end do
    end subroutine f_random_number_d3


    subroutine f_random_number_i0(harvest,range,seed,reset)
      implicit none
      integer, intent(in) :: range
      integer, intent(out) :: harvest
      integer, intent(in), optional :: seed
      logical, intent(in), optional :: reset
      
      if (present(seed)) idum=seed
      harvest=nint(builtin_rand(idum,reset=reset)*real(range,f_simple))

    end subroutine f_random_number_i0

    subroutine f_random_number_i1(harvest,range,seed,reset)
      implicit none
      integer, intent(in) :: range
      integer, dimension(:), intent(out) :: harvest
      integer, intent(in), optional :: seed
      logical, intent(in), optional :: reset
      !local variables
      integer :: i

      if (present(seed)) idum=seed

      do i=1,size(harvest)
         harvest(i) =nint(builtin_rand(idum,reset=reset)*real(range,f_simple))
      end do
    end subroutine f_random_number_i1

    subroutine f_random_number_i1_i1(harvest,ranges,seed,reset)
      implicit none
      integer, dimension(:), intent(in) :: ranges
      integer, dimension(:), intent(out) :: harvest
      integer, intent(in), optional :: seed
      logical, intent(in), optional :: reset
      !local variables
      integer :: i

      if (present(seed)) idum=seed

      do i=1,size(harvest)
         harvest(i) =nint(builtin_rand(idum,reset=reset)*real(ranges(i),f_simple))
      end do
    end subroutine f_random_number_i1_i1


    function builtin_rand(idum, reset)
      !use randomData, only : ntab, iy, iv, start
    
      implicit none
    
      integer, intent(inout) :: idum
      real(f_simple) :: builtin_rand
      logical,intent(in),optional :: reset
      !local variables
      integer, parameter :: ia=16807,im=2147483647,iq=127773,ir=2836,ndiv=1+(im-1)/ntab
      real(f_simple), parameter :: am=1.e0/real(im,kind=f_simple),eps=1.2e-7,rnmx=1.0-eps
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

end module f_random
