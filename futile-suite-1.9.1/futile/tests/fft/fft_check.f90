!> @file
!! Check the FFT routine.
!!
!! @copyright
!!    Copyright (C) Stefan Goedecker, CEA Grenoble, 2002, Basel University, 2009
!!    Copyright (C) 2009-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!>   Check the 3-dimensional complex-complex FFT routine: 
!!   When compared to the best vendor implementations on RISC architectures 
!!   it gives close to optimal performance (perhaps loosing 20 percent in speed)
!!   and it is significanly faster than many not so good vendor implementations 
!!   as well as other portable FFT's. 
!!   On all vector machines tested so far (Cray, NEC, Fujitsu) is 
!!   was significantly faster than the vendor routines
!!
!! The theoretical background is described in :
!! 1) S. Goedecker: Rotating a three-dimensional array in optimal
!! positions for vector processing: Case study for a three-dimensional Fast
!! Fourier Transform, Comp. Phys. Commun. 76, 294 (1993)
!! Citing of this reference is greatly appreciated if the routines are used 
!! for scientific work.
!!
!!   Presumably good compiler flags:
!!   IBM, serial power 2: xlf -qarch=pwr2 -O2 -qmaxmem=-1
!!   with OpenMP: IBM: xlf_r -qfree -O4 -qarch=pwr3 -qtune=pwr3 -qsmp=omp -qmaxmem=-1 ; 
!!                     a.out
!!   DEC: f90 -O3 -arch ev67 -pipeline
!!   with OpenMP: DEC: f90 -O3 -arch ev67 -pipeline -omp -lelan ; 
!!                     prun -N1 -c4 a.out
program fft_check

   use module_fft_sg
   implicit none
   real(kind=8), dimension(:,:), allocatable :: reference
   real(kind=8), dimension(:,:,:), allocatable :: zinout
   integer :: i,ndim,j,inzee,ii
   real(kind=4) :: tt
   real(kind=8) :: maxdiff
   write(*,'(a)') 'FFT test: (n1,n2,n3)'
   !do i=1,ndata
   !   call do_fft(i_data(i), 3, 3)
   !end do

   !test several dimensions of the FFT
   do i=1,ndata

      ndim=i_data(i)
      allocate(zinout(2,ndim,2),reference(2,ndim))
      do j=1,ndim
         call random_number(tt)
         reference(1,j)=real(tt,kind=8)
         reference(2,j)=0.d0!real(tt,kind=8)
      end do
      
      call dcopy(2*ndim,reference(1,1),1,zinout(1,1,1),1)

      !forward FFT
      call fft_1d_ctoc(1,1,ndim,zinout,inzee)
      !copy the results in the first part
      if (inzee == 2) then
         call dcopy(2*ndim,zinout(1,1,inzee),1,zinout(1,1,1),1)
      end if

      !backward FFT
      call fft_1d_ctoc(-1,1,ndim,zinout,inzee)

      call dscal(2*ndim,1.d0/real(ndim,kind=8),zinout(1,1,inzee),1)

      !check
      maxdiff=0.d0
      do j=1,ndim
         do ii=1,2
            maxdiff=max(maxdiff,abs(reference(ii,j)-zinout(ii,j,inzee)))
         end do
      end do
   
      if (maxdiff < 1.d-15) then
         write(*,'(a,i8,1pe25.17)') 'FFT test: dimension, maxdiff:',&
              ndim,maxdiff
      else
         write(*,'(a,i8,1pe25.17,a)') 'FFT test: dimension, maxdiff:',&
              ndim,maxdiff,' <<<<WARNING'
      end if

      
   deallocate(zinout,reference)

   end do

   call do_fft(  7, 16,128)
   call do_fft(  3, 16,128)
   call do_fft(128,128,128)
   !call do_one_fft(128,128,128)


contains

   subroutine do_fft(n1,n2,n3)

      implicit none

      ! dimension parameters
      integer, intent(in) :: n1, n2, n3   
      ! Local variables
      integer :: count1,count2,count_rate,count_max,i,inzee,i_sign
      real(kind=8) :: ttm,tta,t1,t2,tela,time,flops
      ! parameters for FFT
      integer :: nd1, nd2, nd3
      ! general array
      real(kind=8), allocatable :: zin(:,:)
      ! arrays for FFT 
      real(kind=8), allocatable :: z(:,:,:)
      character(len=10) :: message

      nd1=n1+1
      nd2=n2+1
      nd3=n3+1
      write(6,'(3(i5))',advance='no') n1,n2,n3

      ! Allocations
      allocate(zin(2,n1*n2*n3))
      allocate(z(2,nd1*nd2*nd3,2))

      do i=1,nd1*nd2*nd3
         z(1,i,1)=0.d0
         z(2,i,1)=0.d0
         z(1,i,2)=0.d0
         z(2,i,2)=0.d0
      end do

      call init(n1,n2,n3,nd1,nd2,nd3,zin,z)

      i_sign=-1
      inzee=1
      call fft(n1,n2,n3,nd1,nd2,nd3,z,i_sign,inzee)

      call cpu_time(t1)
      call system_clock(count1,count_rate,count_max)      

      i_sign=1
      call fft(n1,n2,n3,nd1,nd2,nd3,z,i_sign,inzee)

      call cpu_time(t2)
      call system_clock(count2,count_rate,count_max)      
      time=(t2-t1)
      tela=(count2-count1)/real(count_rate,kind=8)

      call vgl(n1,n2,n3,nd1,nd2,nd3,z(1,1,inzee), &
                     n1,n2,n3,zin,1.d0/real(n1*n2*n3,kind=8),tta,ttm)
      if (ttm.gt.1.d-10) then
         message = 'Failed'
      else
         message = 'Succeeded'
      end if
      flops=5*n1*n2*n3*log(1.d0*n1*n2*n3)/log(2.d0)
      !write(6,'(a,2(x,e11.4),x,i4)')  'Time (CPU,ELA) per FFT call (sec):' ,time,tela
      !write(6,*) 'Estimated floating point operations per FFT call',flops
      !write(6,*) 'CPU MFlops',1.d-6*flops/time
      write(6,'(1x,a,2(1pg9.2),1x,a)') 'Backw<>Forw:ttm=,tta=',ttm,tta,message

      ! De-allocations
      deallocate(z)
      deallocate(zin)

   end subroutine do_fft

   subroutine init(n1,n2,n3,nd1,nd2,nd3,zin,z)
      implicit none
      !Arguments
      integer, intent(in) :: n1,n2,n3,nd1,nd2,nd3
      real(kind=8) :: zin(2,n1,n2,n3),z(2,nd1,nd2,nd3)
      !Local variables
      integer :: i1,i2,i3
      do i3=1,n3
         do i2=1,n2
            do i1=1,n1
               zin(1,i1,i2,i3) = cos(1.23d0*real(i1*111 + i2*11 + i3,kind=8))
               zin(2,i1,i2,i3) = sin(3.21d0*real(i3*111 + i2*11 + i1,kind=8))
               z(1,i1,i2,i3) = zin(1,i1,i2,i3) 
               z(2,i1,i2,i3) = zin(2,i1,i2,i3) 
            end do
         end do
      end do
   end subroutine init

   subroutine vgl(n1,n2,n3,nd1,nd2,nd3,x,md1,md2,md3,y,scale,tta,ttm)
      implicit none
      !Arguments
      integer, intent(in) :: n1,n2,n3,nd1,nd2,nd3,md1,md2,md3
      real(kind=8), intent(in) :: x(2,nd1,nd2,nd3),y(2,md1,md2,md3)
      real(kind=8), intent(in) :: scale
      !Local variables
      real(kind=8) :: ttm,tta,ttr,tti
      integer :: i1,i2,i3
      ttm=0.d0
      tta=0.d0
      do i3=1,n3
         do i2=1,n2
            do i1=1,n1
               ttr=abs(x(1,i1,i2,i3)*scale-y(1,i1,i2,i3))/abs(y(1,i1,i2,i3))
               tti=abs(x(2,i1,i2,i3)*scale-y(2,i1,i2,i3))/abs(y(2,i1,i2,i3))
               ttm=max(ttr,tti,ttm)
               tta=tta+ttr+tti
            end do
         end do
      end do
      tta=tta/(n1*n2*n3)
   end subroutine vgl

   subroutine do_one_fft(n1,n2,n3)

      implicit none

      ! dimension parameters
      integer, intent(in) :: n1, n2, n3   
      ! Local variables
      integer :: count1,count2,count_rate,count_max,i,inzee,i_sign
      real(kind=8) :: ttm,tta,t1,t2,tela,time
      ! parameters for FFT
      integer :: nd1, nd2, nd3
      ! general array
      real(kind=8), allocatable :: zin(:,:),zout(:,:)
      ! arrays for FFT 
      real(kind=8), allocatable :: z(:,:,:)
      character(len=10) :: message

      nd1=n1+1
      nd2=n2+1
      nd3=n3+1
      write(6,'(3(i5))',advance='no') n1,n2,n3

      ! Allocations
      allocate(zout(2,n1*n2*n3))
      allocate(zin(2,n1*n2*n3))
      allocate(z(2,nd1*nd2*nd3,2))

      do i=1,nd1*nd2*nd3
         z(1,i,1)=0.d0
         z(2,i,1)=0.d0
         z(1,i,2)=0.d0
         z(2,i,2)=0.d0
      end do

      call init_gaussian(n1,n2,n3,nd1,nd2,nd3,zin,z,zout)

      call cpu_time(t1)
      call system_clock(count1,count_rate,count_max)

      i_sign=-1
      inzee=1
      call fft(n1,n2,n3,nd1,nd2,nd3,z,i_sign,inzee)

      call cpu_time(t2)
      call system_clock(count2,count_rate,count_max)

      time=(t2-t1)
      tela=(count2-count1)/real(count_rate,kind=8)

      call vgl_gaussian(n1,n2,n3,nd1,nd2,nd3,z(1,1,inzee), &
                     n1,n2,n3,zout,1.d0/real(n1*n2*n3,kind=8),tta,ttm)
      if (ttm.gt.1.d-10) then
         message = 'Failed'
      else
         message = 'Succeeded'
      end if
      write(6,'(1x,a,2(1pg9.2),1x,a)') 'Backw<>Forw:ttm=,tta=',ttm,tta,message

      ! De-allocations
      deallocate(z)
      deallocate(zin)
      deallocate(zout)

   end subroutine do_one_fft

   !> Stefan convention
   !! Between 1 to n3/2 (0 -> n3/2 - 1)
   !! Between n3/2 to n3 (-n3/2 + 1 -> -1) 
   subroutine init_gaussian(n1,n2,n3,nd1,nd2,nd3,zin,z,zout)
      implicit none
      !Arguments
      integer, intent(in) :: n1,n2,n3,nd1,nd2,nd3
      real(kind=8), intent(inout) :: z(2,nd1,nd2,nd3)
      real(kind=8), intent(out)   :: zin(2,n1,n2,n3),zout(2,n1,n2,n3)
      !Local variables
      real(kind=8), parameter :: a_gauss = 1.0d0,a2 = a_gauss**2
      integer :: i1,i2,i3,j1,j2,j3
      real(kind=8) :: pi,pi2,factor,factor_fft
      real(kind=8) :: r2,x1,x2,x3,g2,p1,p2,p3
      !gaussian 1D e^(- x^2/ a^2) : FFT sqrt(a pi) e^(-pi^2 a^2 k^2)
      pi = 4.d0*atan(1.d0)
      pi2 = pi*pi
      !Normalisation
      factor = 1.d0/(a_gauss*a2*pi*sqrt(pi))
      factor_fft = 1.d0
      !gaussian function
      do i3=1,n3
         x3 = real(i3-n3/2,kind=8)
         j3 = i3+(i3/(n3/2+2))*(n3+2-2*i3)
         p3 = real(j3-1,kind=8)
         do i2=1,n2
            x2 = real(i2-n2/2,kind=8)
            j2 = i2+(i2/(n2/2+2))*(n2+2-2*i2)
            p2 = real(j2-1,kind=8)
            do i1=1,n1
               x1 = real(i1-n1/2,kind=8)
               p1 = real(i1-1,kind=8)
               j1 = i1+(i1/(n1/2+2))*(n1+2-2*i1)
               p1 = real(j1-1,kind=8)
               r2 = x1*x1+x2*x2+x3*x3
               g2 = pi*(p1*p1+p2*p2+p3*p3)
               zin(1,i1,i2,i3) = factor*exp(-r2/a2)
               zin(2,i1,i2,i3) = 0.d0
               z(1,i1,i2,i3) = zin(1,i1,i2,i3)
               z(2,i1,i2,i3) = zin(2,i1,i2,i3)
               zout(1,i1,i2,i3) = factor_fft*exp(-g2*a2)
               zout(2,i1,i2,i3) = 0.d0
            end do
         end do
      end do
   end subroutine init_gaussian

   subroutine vgl_gaussian(n1,n2,n3,nd1,nd2,nd3,x,md1,md2,md3,y,scale,tta,ttm)
      implicit none
      !Arguments
      integer, intent(in) :: n1,n2,n3,nd1,nd2,nd3,md1,md2,md3
      real(kind=8), intent(in) :: x(2,nd1,nd2,nd3),y(2,md1,md2,md3)
      real(kind=8), intent(in) :: scale
      !Local variables
      real(kind=8) :: ttm,tta,ttr,tti
      integer :: i1,i2,i3
      ttm=0.d0
      tta=0.d0
      open(unit=18,form="formatted")
      do i3=1,n3
         do i2=1,n2
            do i1=1,n1
               ttr=abs(x(1,i1,i2,i3)*scale-y(1,i1,i2,i3))
               tti=abs(x(2,i1,i2,i3)*scale-y(2,i1,i2,i3))
               ttm=max(ttr,tti,ttm)
               tta=tta+ttr+tti
               if (i1 == 1 .and. i2 == 1) then
                  write(18,*) i1,i2,i3,x(1,i1,i2,i3),y(1,i1,i2,i3)
                  print *,i1,i2,i3,x(:,i1,i2,i3),y(:,i1,i2,i3),ttr,tti,ttm,tta
               end if
            end do
         end do
      end do
      tta=tta/(n1*n2*n3)
      close(unit=18)
   end subroutine vgl_gaussian

end program fft_check
