!> @file
!!  Routines to initialize to zero arrays
!! @author
!!    Copyright (C) 2009-2015 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 

!SM: overload subroutines such that they can be called with size n being a 4 or 8 byte integer.
!!   Define the interface with 
module module_razero

  private

  public :: razero, razero_complex, razero_simple, razero_integer, razero_integerlong

  interface razero
    module procedure razero_i_interface, razero_li_interface
  end interface razero

  interface razero_complex
    module procedure razero_complex_i_interface, razero_complex_li_interface
  end interface razero_complex

  interface razero_simple
    module procedure razero_simple_i_interface, razero_simple_li_interface
  end interface razero_simple

  interface razero_integer
    module procedure razero_integer_i_interface, razero_integer_li_interface
  end interface razero_integer

  interface razero_integerlong
    module procedure razero_integerlong_i_interface, razero_integerlong_li_interface
  end interface razero_integerlong


  contains

    subroutine razero_i_interface(n,x)
      implicit none
      !Arguments
      integer(kind=4), intent(in) :: n
      double precision, intent(out) :: x
      call razero_i(n,x)
    end subroutine razero_i_interface

    subroutine razero_li_interface(n,x)
      implicit none
      !Arguments
      integer(kind=8), intent(in) :: n
      double precision, intent(out) :: x
      call razero_li(n,x)
    end subroutine razero_li_interface


    subroutine razero_complex_i_interface(n,x)
      implicit none
      !Arguments
      integer(kind=4), intent(in) :: n
      double complex, intent(out) :: x
      call razero_complex_i(n,x)
    end subroutine razero_complex_i_interface

    subroutine razero_complex_li_interface(n,x)
      implicit none
      !Arguments
      integer(kind=8), intent(in) :: n
      double complex, intent(out) :: x
      call razero_complex_li(n,x)
    end subroutine razero_complex_li_interface


    subroutine razero_simple_i_interface(n,x)
      implicit none
      !Arguments
      integer(kind=4), intent(in) :: n
      real(kind=4), intent(out) :: x
      call razero_simple_i(n,x)
    end subroutine razero_simple_i_interface

    subroutine razero_simple_li_interface(n,x)
      implicit none
      !Arguments
      integer(kind=8), intent(in) :: n
      real(kind=4), intent(out) :: x
      call razero_simple_li(n,x)
    end subroutine razero_simple_li_interface


    subroutine razero_integer_i_interface(n,x)
      implicit none
      !Arguments
      integer(kind=4), intent(in) :: n
      integer(kind=4), intent(out) :: x
      call razero_integer_i(n,x)
    end subroutine razero_integer_i_interface

    subroutine razero_integer_li_interface(n,x)
      implicit none
      !Arguments
      integer(kind=8), intent(in) :: n
      integer(kind=4), intent(out) :: x
      call razero_integer_li(n,x)
    end subroutine razero_integer_li_interface


    subroutine razero_integerlong_i_interface(n,x)
      implicit none
      !Arguments
      integer(kind=4), intent(in) :: n
      integer(kind=8), intent(out) :: x
      call razero_integerlong_i(n,x)
    end subroutine razero_integerlong_i_interface

    subroutine razero_integerlong_li_interface(n,x)
      implicit none
      !Arguments
      integer(kind=8), intent(in) :: n
      integer(kind=8), intent(out) :: x
      call razero_integerlong_li(n,x)
    end subroutine razero_integerlong_li_interface

end module module_razero

!> Routine initialize double precision arrays to zero
subroutine razero_i(n,x)
  implicit none
  !Arguments
  integer(kind=4), intent(in) :: n
  double precision, dimension(n), intent(out) :: x
  !Local variables
  integer :: i
  !$ logical :: omp_in_parallel,do_omp
  !$ do_omp = n > 1024
  !$ if (do_omp) do_omp= .not. omp_in_parallel()
  !$omp parallel if (do_omp) shared(x,n) private(i)
  !$omp do
  do i=1,n
      x(i)=0.d0
  end do
  !$omp enddo
  !$omp end parallel
end subroutine razero_i

subroutine razero_li(n,x)
  implicit none
  !Arguments
  integer(kind=8), intent(in) :: n
  double precision, dimension(n), intent(out) :: x
  !Local variables
  integer(kind=8) :: i
  !$ logical :: omp_in_parallel,do_omp
  !$ do_omp = n > 1024
  !$ if (do_omp) do_omp= .not. omp_in_parallel()
  !$omp parallel if (do_omp) shared(x,n) private(i)
  !$omp do
  do i=1,n
      x(i)=0.d0
  end do
  !$omp enddo
  !$omp end parallel
end subroutine razero_li


!> Routine initialize double complex arrays to zero
subroutine razero_complex_i(n,x)
  implicit none
  !Arguments
  integer(kind=4), intent(in) :: n
  double complex, dimension(n), intent(out) :: x
  !Local variables
  integer :: i
  !$ logical :: omp_in_parallel,do_omp
  !$ do_omp = n > 1024
  !$ if (do_omp) do_omp= .not. omp_in_parallel()
  !$omp parallel if (do_omp) shared(x,n) private(i)
  !$omp do
  do i=1,n
      x(i)=(0.d0,0.d0)
  end do
  !$omp enddo
  !$omp end parallel
end subroutine razero_complex_i

subroutine razero_complex_li(n,x)
  implicit none
  !Arguments
  integer(kind=8), intent(in) :: n
  double complex, dimension(n), intent(out) :: x
  !Local variables
  integer(kind=8) :: i
  !$ logical :: omp_in_parallel,do_omp
  !$ do_omp = n > 1024
  !$ if (do_omp) do_omp= .not. omp_in_parallel()
  !$omp parallel if (do_omp) shared(x,n) private(i)
  !$omp do
  do i=1,n
      x(i)=(0.d0,0.d0)
  end do
  !$omp enddo
  !$omp end parallel
end subroutine razero_complex_li


!> Set to zero an array x(n)
subroutine razero_simple_i(n,x)
  implicit none
  !Arguments
  integer(kind=4), intent(in) :: n
  real(kind=4), intent(out) :: x(n)
  !Local variables
  integer :: i
  !integer :: m
  !$ logical :: omp_in_parallel,do_omp
  !$ do_omp = n > 1024
  !$ if (do_omp) do_omp= .not. omp_in_parallel()
  !$omp parallel if (do_omp) shared(x,n) private(i)
  !$omp do
  do i=1,n
      x(i)=0.e0
  end do
  !$omp enddo
  !$omp end parallel
END SUBROUTINE razero_simple_i


subroutine razero_simple_li(n,x)
  implicit none
  !Arguments
  integer(kind=8), intent(in) :: n
  real(kind=4), intent(out) :: x(n)
  !Local variables
  integer(kind=8) :: i
  !$ logical :: omp_in_parallel,do_omp
  !$ do_omp = n > 1024
  !$ if (do_omp) do_omp= .not. omp_in_parallel()
  !$omp parallel if (do_omp) shared(x,n) private(i)
  !$omp do
  do i=1,n
      x(i)=0.e0
  end do
  !$omp enddo
  !$omp end parallel
END SUBROUTINE razero_simple_li


!> Set to zero an array x(n)
subroutine razero_integer_i(n,x)
  implicit none
  !Arguments
  integer(kind=4), intent(in) :: n
  integer(kind=4), dimension(n), intent(out) :: x
  !Local variables
  integer :: i
  !integer :: m
  !$ logical :: omp_in_parallel,do_omp
  !$ do_omp = n > 1024
  !$ if (do_omp) do_omp= .not. omp_in_parallel()
  !$omp parallel if (do_omp) shared(x,n) private(i)
  !$omp do
  do i=1,n
      x(i)=0
  end do
  !$omp enddo
  !$omp end parallel
END SUBROUTINE razero_integer_i

subroutine razero_integer_li(n,x)
  implicit none
  !Arguments
  integer(kind=8), intent(in) :: n
  integer(kind=4), dimension(n), intent(out) :: x
  !Local variables
  integer(kind=8) :: i
  !$ logical :: omp_in_parallel,do_omp
  !$ do_omp = n > 1024
  !$ if (do_omp) do_omp= .not. omp_in_parallel()
  !$omp parallel if (do_omp) shared(x,n) private(i)
  !$omp do
  do i=1,n
      x(i)=0
  end do
  !$omp enddo
  !$omp end parallel
END SUBROUTINE razero_integer_li



!> Set to zero an array x(n)
subroutine razero_integerlong_i(n,x)
  implicit none
  !Arguments
  integer(kind=4), intent(in) :: n
  integer(kind=8), dimension(n), intent(out) :: x
  !Local variables
  integer :: i
  !integer :: m
  !$ logical :: omp_in_parallel,do_omp
  !$ do_omp = n > 1024
  !$ if (do_omp) do_omp= .not. omp_in_parallel()
  !$omp parallel if (do_omp) shared(x,n) private(i)
  !$omp do
  do i=1,n
      x(i)=int(0,kind=8)
  end do
  !$omp enddo
  !$omp end parallel
END SUBROUTINE razero_integerlong_i

subroutine razero_integerlong_li(n,x)
  implicit none
  !Arguments
  integer(kind=8), intent(in) :: n
  integer(kind=8), dimension(n), intent(out) :: x
  !Local variables
  integer(kind=8) :: i
  !$ logical :: omp_in_parallel,do_omp
  !$ do_omp = n > 1024
  !$ if (do_omp) do_omp= .not. omp_in_parallel()
  !$omp parallel if (do_omp) shared(x,n) private(i)
  !$omp do
  do i=1,n
      x(i)=int(0,kind=8)
  end do
  !$omp enddo
  !$omp end parallel
END SUBROUTINE razero_integerlong_li


!!!> Set to zero an array x(n): omp version of razero
!!subroutine omp_razero(n,x)
!!  use module_base
!!  implicit none
!!  !Arguments
!!  integer, intent(in) :: n
!!  real(kind=8), intent(out) :: x(n)
!!  !Local variables
!!  integer :: i,is
!!
!!
!!!!!$omp do
!!      do i=1,n-7,8
!!      x(i+0)=0.d0
!!      x(i+1)=0.d0
!!      x(i+2)=0.d0
!!      x(i+3)=0.d0
!!      x(i+4)=0.d0
!!      x(i+5)=0.d0
!!      x(i+6)=0.d0
!!      x(i+7)=0.d0
!!      x(i+8)=0.d0
!!      end do
!!!!!$omp enddo
!!      is=i
!!      do i=is,n
!!      x(i)=0.d0
!!      end do
!!END SUBROUTINE omp_razero


!> Set to 10^-20 an array x(n) for exchange-correlation function of ABINIT
subroutine tenminustwenty(n,x,nproc)
  implicit none
! Arguments
  integer :: n,nproc
  real(kind=8) :: x(n)
! Local variables
  integer :: i
  do i=1,n
     x(i)=1.d-20/real(nproc,kind=8)
  end do
END SUBROUTINE tenminustwenty


!> Set to 10^-10 an array x(n) for exchange-correlation function of ABINIT.
!! We use 10^-10 here since the array will be squared later and we then arrive at
!! the desired 10^-20.
subroutine tenminusten(n,x,nproc)
  implicit none
! Arguments
  integer :: n,nproc
  real(kind=8) :: x(n)
! Local variables
  integer :: i
  do i=1,n
     x(i)=1.d-10/real(nproc,kind=8)
  end do
END SUBROUTINE tenminusten


!> Routine doing daxpy with dx in real(kind=4)
subroutine dasxpdy(n,da,dx,incx,dy,incy)
  implicit none
  integer, intent(in) :: n,incx,incy
  real(kind=8), intent(in) :: da
  real(kind=4), dimension(*), intent(in) :: dx
  real(kind=8), dimension(*), intent(inout) :: dy
  !local variables
  integer :: i,ix,iy
  
  ix=1
  iy=1
  do i=1,n
     dy(iy)=dy(iy)+da*real(dx(ix),kind=8)
     ix=ix+incx
     iy=iy+incy
  end do
end subroutine dasxpdy


!> Copy from real(kind=4) into real(kind=8)
subroutine dscopy(n,dx,incx,dy,incy)
  implicit none
  integer, intent(in) :: n,incx,incy
  real(kind=8), dimension(*), intent(in) :: dx
  real(kind=4), dimension(*), intent(out) :: dy
  !local variables
  integer :: i,ix,iy
  
  ix=1
  iy=1
  do i=1,n
     dy(iy)=real(dx(ix),kind=4)
     ix=ix+incx
     iy=iy+incy
  end do

end subroutine dscopy


!> dcopy for integer arrays
subroutine icopy(n,dx,incx,dy,incy)
  implicit none
  integer, intent(in) :: n,incx,incy
  integer, dimension(*), intent(in) :: dx
  integer, dimension(*), intent(out) :: dy
  !local variables
  integer :: i,ix,iy
  
  ix=1
  iy=1
  do i=1,n
     dy(iy)=dx(ix)
     ix=ix+incx
     iy=iy+incy
  end do

end subroutine icopy

subroutine diff_i(n,a,b,diff,idiff)
  use f_precisions, only: f_long
  implicit none
  integer(f_long), intent(in) :: n
  integer, dimension(n), intent(in) :: a
  integer, dimension(n), intent(in) :: b
  integer, intent(out) :: diff
  integer(f_long), intent(out) :: idiff
  !local variables
  integer(f_long) :: i

  diff=0
  do i=1,n
     if (diff < abs(a(i)-b(i))) then
        diff=max(diff,abs(a(i)-b(i)))
        idiff=i
     end if
  end do
end subroutine diff_i


subroutine diff_li(n,a,b,diff,idiff)
  use f_precisions, only: f_long
  implicit none
  integer(f_long), intent(in) :: n
  integer(kind=8), dimension(n), intent(in) :: a
  integer(kind=8), dimension(n), intent(in) :: b
  integer(kind=8), intent(out) :: diff
  integer(f_long), intent(out) :: idiff
  !local variables
  integer(f_long) :: i

  diff=int(0,kind=8)
  do i=1,n
     if (diff < abs(a(i)-b(i))) then
        diff=max(diff,abs(a(i)-b(i)))
        idiff=i
     end if
  end do
end subroutine diff_li


subroutine diff_r(n,a,b,diff,idiff)
  use f_precisions, only: f_long
  implicit none
  integer(f_long), intent(in) :: n
  real, dimension(n), intent(in) :: a
  real, dimension(n), intent(in) :: b
  real, intent(out) :: diff
  integer(f_long), intent(out) :: idiff
  !local variables
  integer(f_long) :: i

  diff=0.0e0
  do i=1,n
     if (diff < abs(a(i)-b(i))) then
        diff=max(diff,abs(a(i)-b(i)))
        idiff=i
     end if
  end do
end subroutine diff_r


subroutine diff_d(n,a,b,diff,idiff)
  use f_precisions, only: f_long
  implicit none
  integer(f_long), intent(in) :: n
  double precision, dimension(n), intent(in) :: a
  double precision, dimension(n), intent(in) :: b
  double precision, intent(out) :: diff
  integer(f_long), intent(out) :: idiff
  !local variables
  integer(f_long) :: i

  diff=0.0d0
  do i=1,n
     if (diff < abs(a(i)-b(i))) then
        diff=max(diff,abs(a(i)-b(i)))
        idiff=i
     end if
  end do
end subroutine diff_d


subroutine diff_l(n,a,b,diff)
  use f_precisions, only: f_long
  implicit none
  integer(f_long), intent(in) :: n
  logical, dimension(n), intent(in) :: a
  logical, dimension(n), intent(in) :: b
  logical, intent(out) :: diff
  !local variables
  integer(f_long) :: i

  diff=.false.
  do i=1,n
     diff=a(i) .eqv. b(i)
     if (diff) exit
  end do
end subroutine diff_l


subroutine diff_ci(n,a,b,diff,idiff)
  use f_precisions, only: f_long
  implicit none
  integer(f_long), intent(in) :: n
  character, dimension(n), intent(in) :: a
  integer, dimension(n), intent(in) :: b
  integer, intent(out) :: diff
  integer(f_long), intent(out) :: idiff
  !local variables
  integer(f_long) :: i

  diff=0
  do i=1,n
     if (diff < abs(ichar(a(i))-b(i))) then
        diff=max(diff,abs(ichar(a(i))-b(i)))
        idiff=i
     end if
  end do
end subroutine diff_ci


subroutine f_itoa(n,src,dest)
  use f_precisions, only: f_long
  implicit none
  integer(f_long), intent(in) :: n
  integer(kind=4), dimension(n), intent(in) :: src
  character, dimension(n), intent(out) :: dest
  !local variables
  integer(f_long) :: i
  
  do i=1,n
     dest(i)=achar(src(i))
  end do

end subroutine f_itoa


subroutine f_litoa(n,src,dest)
  use f_precisions, only: f_long
  implicit none
  integer(f_long), intent(in) :: n
  integer(f_long), dimension(n), intent(in) :: src
  character, dimension(n), intent(out) :: dest
  !local variables
  integer(f_long) :: i
  
  do i=1,n
     dest(i)=achar(src(i))
  end do

end subroutine f_litoa

subroutine f_atoi(n,src,dest)
  use f_precisions, only: f_integer,f_long
  implicit none
  integer(f_long), intent(in) :: n
  character, dimension(n), intent(in) :: src
  integer(f_integer), dimension(n), intent(out) :: dest
  !local variables
  integer(f_long) :: i
  
  do i=1,n
     dest(i)=ichar(src(i))
  end do

end subroutine f_atoi

subroutine f_atoli(n,src,dest)
  use f_precisions, only: f_long
  implicit none
  integer(f_long), intent(in) :: n
  character, dimension(n), intent(in) :: src
  integer(f_long), dimension(n), intent(out) :: dest
  !local variables
  integer(f_long) :: i
  
  do i=1,n
     dest(i)=ichar(src(i))
  end do

end subroutine f_atoli

!> set nbytes term to zero
subroutine setzero(nbytes,x)
  use f_precisions, only: f_long
  implicit none
  integer(f_long), intent(in) :: nbytes
  character, dimension(nbytes), intent(inout) :: x
  !local variables
  integer :: nthreads,ithread
  !$ integer omp_get_max_threads,omp_get_thread_num
  integer(f_long) :: nbt,it,nt,nb

  if (nbytes == 0_f_long) return
  nthreads=1
  ithread=0
  !$ nthreads=omp_get_max_threads()
  nt=int(nthreads,f_long)
  nbt=(nbytes+nt-1)/nt
  !$omp parallel private(ithread,it,nb) 
  !$ ithread=omp_get_thread_num()
  !calculate the number of elements for each thread
  it=min(nbt*ithread,nbytes-1)
  nb=min(nbytes-it,nbt)
  call memsetzero(x(it+1),nb)
  !$omp end parallel
end subroutine setzero
