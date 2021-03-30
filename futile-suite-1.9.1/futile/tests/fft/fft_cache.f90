!> @file
!!   Calculate the optimal cache parameter for the FFT routine.
!! @copyright
!!   Copyright (C) Stefan Goedecker, CEA Grenoble, 2002, Basel University, 2009
!!   Copyright (C) 2009-2013 BigDFT group
!!   This file is distributed under the terms of the
!!   GNU General Public License, see ~/COPYING file
!!   or http://www.gnu.org/copyleft/gpl.txt .
!!   For the list of contributors, see ~/AUTHORS


!>   Check the best cache size
!!   3-dimensional complex-complex FFT routine: 
!!   When compared to the best vendor implementations on RISC architectures 
!!   it gives close to optimal performance (perhaps loosing 20 percent in speed)
!!   and it is significanly faster than many not so good vendor implementations 
!!   as well as other portable FFT's. 
!!   On all vector machines tested so far (Cray, NEC, Fujitsu) is 
!!   was significantly faster than the vendor routines
!!
!!  The theoretical background is described in :
!!  S. Goedecker: Rotating a three-dimensional array in optimal
!!  positions for vector processing: Case study for a three-dimensional Fast
!!  Fourier Transform, Comp. Phys. Commun. 76, 294 (1993)
!!  Citing of this reference is greatly appreciated if the routines are used 
!!  for scientific work.
program fft_cache

   use module_fft_sg
   implicit none
   integer, parameter :: nsize_dat = 20000000
   integer, parameter :: iunit=21
   character (len=100) :: tatonam
   integer :: ip,jp,n3,ndat,ntime
   real(kind=8) :: time,tela
   logical :: success

   ! Get arguments
   call getarg(1,tatonam)
   
   if(trim(tatonam)=='') then
      write(*,'(1x,a)')&
           'Usage: ./fft_cache <cache_size>'
      write(*,'(1x,a)')&
           'Indicate the cache size in kB for FFT test'
      stop
   else
      !Read the cache size
      read(unit=tatonam,fmt=*) ncache
      ncache = ncache*1024
      write(unit=6,fmt="(a,i0,a)",advance="no") "Cache size=",ncache,":"
   end if

   open(unit=iunit,file="fft_cache.dat",status="unknown",position="append")
   do ip=1,ndata,10
      n3 = i_data(ip)
      !n3 = 315
      ndat = nsize_dat/8/n3
      ntime=1!3
      do jp=1,ntime
         call do_fft(ndat, n3, time, tela, success)
         if (.not. success) exit
      end do
      if (.not. success) exit
      write(unit=6,fmt="(a,i0,a)",advance="no") "[",n3,"]"
      write(unit=iunit,fmt=*) ncache,n3,ndat,time,tela
   end do
   write(unit=iunit,fmt=*)
   close(unit=iunit)
   write(unit=6,fmt="(a)") 

contains

   subroutine do_fft(ndat,n3,time,tela,success)

      implicit none

      ! dimension parameters
      integer, intent(in) :: ndat, n3
      real(kind=8), intent(out) :: time,tela
      logical, intent(out) :: success
      ! Local variables
      integer :: count1,count2,count_rate,count_max,i,inzee,i_sign
      real(kind=8) :: t1,t2
      ! parameters for FFT
      integer :: nd3, nddat
      ! general array
      real(kind=8), allocatable :: zin(:,:)
      ! arrays for FFT 
      real(kind=8), allocatable :: z(:,:,:)

      nd3=n3+1
      nddat=ndat+1

      ! Allocations
      allocate(zin(2,ndat*n3))
      allocate(z(2,nddat*nd3,2))

      do i=1,nddat*nd3
         z(1,i,1)=0.d0
         z(2,i,1)=0.d0
         z(1,i,2)=0.d0
         z(2,i,2)=0.d0
      end do

      call init(ndat,n3,nddat,nd3,zin,z)

      i_sign=-1
      inzee=1
      call cpu_time(t1)
      call system_clock(count1,count_rate,count_max)      
      call fft1(ndat,n3,nddat,nd3,z,i_sign,inzee,success)
      call system_clock(count2,count_rate,count_max)      
      call cpu_time(t2)
      time=(t2-t1)
      tela=(count2-count1)/real(count_rate,kind=8)

      ! De-allocations
      deallocate(z)
      deallocate(zin)

   end subroutine do_fft

   subroutine init(ndat,n3,nddat,nd3,zin,z)
      implicit none
      !Arguments
      integer, intent(in) :: ndat,n3,nddat,nd3
      real(kind=8) :: zin(2,ndat,n3),z(2,nddat,nd3)
      !Local variables
      integer :: id,i3
      do i3=1,n3
         do id=1,ndat
            zin(1,id,i3) = cos(1.23d0*real(id*11 + i3,kind=8))
            zin(2,id,i3) = sin(3.21d0*real(i3*11 + id,kind=8))
            z(1,id,i3) = zin(1,id,i3) 
            z(2,id,i3) = zin(2,id,i3) 
         end do
      end do
   end subroutine init

   !> Do one basic FFT (transform along Z)
   subroutine fft1(ndat,n3,nddat,nd3,z,i_sign,inzee, success)

      use module_fft_sg
      implicit real(kind=8) (a-h,o-z), integer (i-n)

      !!!$      interface
      !!!!$        integer ( kind=4 ) function omp_get_num_threads ( )
      !!!!$        end function omp_get_num_threads
      !!!$      end interface
      !!!!$      interface
      !!!!$        integer ( kind=4 ) function omp_get_thread_num ( )
      !!!!$        end function omp_get_thread_num
      !!!!$      end interface

      !Arguments
      integer, intent(in) :: ndat,n3,nddat,nd3,i_sign
      integer, intent(inout) :: inzee
      logical, intent(out) :: success
      real(kind=8), intent(inout) :: z(2,nddat*nd3,2)
      !Local variables
      real(kind=8), dimension(:,:), allocatable :: trig
      integer, dimension(n_factors) :: after,now,before
      real(kind=8), allocatable, dimension(:,:,:) :: zw  

      if (n3.gt.nfft_max) then
         write(*,*) 'Dimension bigger than ', nfft_max
         stop
      end if

      success = .true.

      ntrig=n3
      allocate(trig(2,ntrig))

      ! check whether input values are reasonable
      if (inzee.le.0 .or. inzee.ge.3) stop 'wrong inzee'
      if (i_sign.ne.1 .and. i_sign.ne.-1) stop 'wrong i_sign'
      if (n3.gt.nd3) stop 'n3>nd3'
       
      ! vector computer with memory banks:
      if (ncache.eq.0) then
         call ctrig_sg(n3,ntrig,trig,after,before,now,i_sign,ic)
         nfft=ndat
         mm=nddat
         do i=1,ic-1
            call fftstp_sg(mm,nfft,nd3,mm,nd3,z(1,1,inzee),z(1,1,3-inzee), &
                           ntrig,trig,after(i),now(i),before(i),i_sign)
            inzee=3-inzee
         end do
         i=ic
         call fftrot_sg(mm,nfft,nd3,mm,nd3,z(1,1,inzee),z(1,1,3-inzee), &
              ntrig,trig,after(i),now(i),before(i),i_sign)
         inzee=3-inzee


      ! RISC machine with cache:
      else
         ! Intel IFC does not understand default(private)
         !!!!!$omp parallel  default(private) &
         !!!!$omp parallel & 
         !!!!$omp private(zw,trig,before,after,now,i,j,iam,npr,jj,ma,mb,mm,ic,n,m,jompa,jompb,lot,lotomp,inzeep,inzet,nn,nfft) &
         !!!!$omp shared(n1,n2,n3,nd1,nd2,nd3,z,i_sign,inzee,ncache) 
         npr=1
         !!!!$       npr=omp_get_num_threads()
         iam=0
         !!!!$       iam=omp_get_thread_num()
         !      write(6,*) 'npr,iam',npr,iam
         ! Critical section only necessary on Intel
         !!!!$omp critical
         allocate(zw(2,ncache/4,2))
         !!!!$omp end critical

         inzet=inzee
         ! TRANSFORM ALONG Z AXIS

         mm=nddat
         m=nd3
         lot=max(1,ncache/(4*n3))
         nn=lot
         n=n3
         if (2*n*lot*2.gt.ncache) then
            write(*,"(A)",advance="no") 'STOP ncache1'
            success = .false.
            return
         end if

         call ctrig_sg(n3,ntrig,trig,after,before,now,i_sign,ic)

         if (ic.eq.1) then
            i=ic
            lotomp=(ndat)/npr+1
            ma=iam*lotomp+1
            mb=min((iam+1)*lotomp,ndat)
            nfft=mb-ma+1
            j=ma
            jj=j*nd3-nd3+1
            call fftrot_sg(mm,nfft,m,mm,m,z(1,j,inzet),z(1,jj,3-inzet), &
                           ntrig,trig,after(i),now(i),before(i),i_sign)

         else

            lotomp=(ndat)/npr+1
            jompa=iam*lotomp+1
            jompb=min((iam+1)*lotomp,ndat)
            do j=jompa,jompb,lot
               ma=j
               mb=min(j+(lot-1),jompb)
               nfft=mb-ma+1
               jj=j*nd3-nd3+1
              
               i=1
               inzeep=2
               call fftstp_sg(mm,nfft,m,nn,n,z(1,j,inzet),zw(1,1,3-inzeep), &
                              ntrig,trig,after(i),now(i),before(i),i_sign)
               inzeep=1
              
               do i=2,ic-1
                  call fftstp_sg(nn,nfft,n,nn,n,zw(1,1,inzeep),zw(1,1,3-inzeep), &
                                 ntrig,trig,after(i),now(i),before(i),i_sign)
                  inzeep=3-inzeep
               end do
               i=ic
               call fftrot_sg(nn,nfft,n,mm,m,zw(1,1,inzeep),z(1,jj,3-inzet), &
                              ntrig,trig,after(i),now(i),before(i),i_sign)
            end do
         end if

         inzet=3-inzet

         !!!!!!!!!$omp barrier

         deallocate(zw)
         if (iam.eq.0) inzee=inzet

      end if
      deallocate(trig)

   end subroutine fft1


end program fft_cache
