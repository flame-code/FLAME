!> @file
!!  Define routne FFT 2D
!! @author
!!  Copyright (C) Stefan Goedecker, Lausanne, Switzerland, August 1, 1991
!!  Copyright (C) Stefan Goedecker, Cornell University, Ithaca, USA, 1994
!!  Copyright (C) Stefan Goedecker, MPI Stuttgart, Germany, 1999
!!  Copyright (C) Stefan Goedecker, CEA Grenoble, 2002
!!  Copyright (C) 2002-2009 BigDFT group 
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~/COPYING file
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!  For the list of contributors, see ~/AUTHORS 


!> CALCULATES THE DISCRETE FOURIER TRANSFORM F(I1,I2)=
!!    S_(j1,j2) EXP(isign*i*2*pi*(j1*i1/n1+j2*i2/n2)) R(j1,j2)
!!    INPUT:
!!        n1,n2:physical dimension of the transform. It must be a 
!!          product of the prime factors 2,3,5, but greater than 3. 
!!                 If two ni's are equal it is recommended to place them 
!!                 behind each other.
!!        nd1,nd2:memory dimension of Z. ndi must always be greater or 
!!            equal than ni. It is recomended to chose ndi=ni 
!!            if ni is odd and ndi=ni+1 if ni is even to obtain 
!!                   optimal execution speed. For small ni it can however 
!!                   sometimes be advantageous to set ndi=ni
!!           inzee=1: first part of Z is data array, second part work array
!!           inzee=2: first part of Z is work array, second part data array
!!        Z(1,i1,i2,inzee)=real(R(i1,i2))
!!        Z(2,i1,i2,inzee)=imag(R(i1,i2))
!!    OUTPUT:
!!           inzee=1: first part of Z is data array, second part work array
!!           inzee=2: first part of Z is work array, second part data array
!!        real(F(i1,i2))=Z(1,i1,i2,inzee)
!!        imag(F(i1,i2))=Z(2,i1,i2,inzee)
!!            inzee on output is in general different from inzee on input
!!    THE ARRAYELEMENTS Z( , , , ,3-inzee) ARE USED AS WORKING SPACE
!!       On a RISC machine with cache, it is very important to find the optimal 
!!       value of NCACHE. NCACHE determines the size of the work array zw, that
!!       has to fit into cache. It has therefore to be chosen to equal roughly 
!!    half the size of the physical cache in units of real*8 numbers.
!!       The optimal value of ncache can easily be determined by numerical 
!!       experimentation. A too large value of ncache leads to a dramatic 
!!       and sudden decrease of performance, a too small value to a to a 
!!       slow and less dramatic decrease of performance. If NCACHE is set 
!!       to a value so small, that not even a single one dimensional transform 
!!       can be done in the workarray zw, the program stops with an error message.
!!    On a vector machine ncache has to be put to 0 
!!
subroutine FFT2d(n1,n2,nd1,nd2,z,isign,inzee,zw,ncache)

   implicit none
   !Arguments
   integer, intent(in) :: n1,n2,nd1,nd2,isign,ncache
   integer, intent(inout) :: inzee
   integer, dimension(20) :: after,before,now
   real(kind=8), dimension(2,nd1*nd2,2) :: z
   real(kind=8), dimension(2,ncache/4,2) :: zw
   !Local variables
   integer, parameter :: nmax = 1024
   real(kind=8), dimension(2,nmax) :: trig
   integer :: ntrig,nn,nfft,n,mm,mb,ma,j,m,lot,jp,jj,jb,ja,inzeep
   integer :: ic,i

   if (max(n1,n2).gt.nmax) stop '1024'
   ntrig=nmax
   ! vector computer with memory banks:
   if (ncache.eq.0) then

      ! TRANSFORM ALONG Y AXIS
      call ctrig_sg(n2,ntrig,trig,after,before,now,isign,ic)
      nfft=n1
      mm=nd1
      do i=1,ic-1
         call fftstp_sg(mm,nfft,nd2,mm,nd2,z(1,1,inzee),z(1,1,3-inzee),ntrig,trig,after(i),now(i),before(i),isign)
         inzee=3-inzee
      end do
      i=ic
      call fftrot_sg(mm,nfft,nd2,mm,nd2,z(1,1,inzee),z(1,1,3-inzee),ntrig,trig,after(i),now(i),before(i),isign)
      inzee=3-inzee

      ! TRANSFORM ALONG X AXIS
      if (n1.ne.n2) call ctrig_sg(n1,ntrig,trig,after,before,now,isign,ic)
      nfft=n2
      mm=nd2
      do i=1,ic-1
         call fftstp_sg(mm,nfft,nd1,mm,nd1,z(1,1,inzee),z(1,1,3-inzee),ntrig,trig,after(i),now(i),before(i),isign)
         inzee=3-inzee
      end do
      i=ic
      call fftrot_sg(mm,nfft,nd1,mm,nd1,z(1,1,inzee),z(1,1,3-inzee),ntrig,trig,after(i),now(i),before(i),isign)
      inzee=3-inzee

      ! RISC machine with cache:
   else


      ! TRANSFORM ALONG Y AXIS
      mm=nd1
      m=nd2
      lot=max(1,ncache/(4*n2))
      !    print*,'lot',lot
      nn=lot
      n=n2
      if (2*n*lot*2.gt.ncache) stop 'enlarge ncache :2'
      call ctrig_sg(n2,ntrig,trig,after,before,now,isign,ic)

      ja=1
      jb=n1

      if (ic.eq.1) then
         i=ic
         jj=ja*nd2-nd2+1
         nfft=jb-ja+1
         call fftrot_sg(mm,nfft,m,mm,m,z(1,ja,inzee),z(1,jj,3-inzee),ntrig,trig,after(i),now(i),before(i),isign)

      else

         do j=ja,jb,lot
            ma=j
            mb=min(j+(lot-1),jb)
            nfft=mb-ma+1
            jj=j*nd2-nd2+1

            i=1
            inzeep=2
            call fftstp_sg(mm,nfft,m,nn,n,z(1,j,inzee),zw(1,1,3-inzeep),ntrig,trig,after(i),now(i),before(i),isign)
            inzeep=1

            do i=2,ic-1
               call fftstp_sg(nn,nfft,n,nn,n,zw(1,1,inzeep),zw(1,1,3-inzeep),ntrig,trig,after(i),now(i),before(i),isign)
               inzeep=3-inzeep
            end do

            i=ic
            call fftrot_sg(nn,nfft,n,mm,m,zw(1,1,inzeep),z(1,jj,3-inzee),ntrig,trig,after(i),now(i),before(i),isign)
         end do
      endif
      inzee=3-inzee


      ! TRANSFORM ALONG X AXIS
      mm=nd2
      m=nd1
      lot=max(1,ncache/(4*n1))
      !    print*,'lot',lot
      nn=lot
      n=n1
      if (2*n*lot*2.gt.ncache) stop 'enlarge ncache :1'
      if (n1.ne.n2) call ctrig_sg(n1,ntrig,trig,after,before,now,isign,ic)

      if (ic.eq.1) then
         i=ic
         j=1
         jp=nd2
         jj=j*nd1-nd1+1
         nfft=jp-j+1
         call fftrot_sg(mm,nfft,m,mm,m,z(1,j,inzee),z(1,jj,3-inzee),ntrig,trig,after(i),now(i),before(i),isign)

      else

         do j=1,n2,lot
            ma=j
            mb=min(j+(lot-1),n2)
            nfft=mb-ma+1
            jj=j*nd1-nd1+1

            i=1
            inzeep=2
            call fftstp_sg(mm,nfft,m,nn,n,z(1,j,inzee),zw(1,1,3-inzeep),ntrig,trig,after(i),now(i),before(i),isign)
            inzeep=1

            do i=2,ic-1
               call fftstp_sg(nn,nfft,n,nn,n,zw(1,1,inzeep),zw(1,1,3-inzeep),ntrig,trig,after(i),now(i),before(i),isign)
               inzeep=3-inzeep
            end do
            i=ic
            call fftrot_sg(nn,nfft,n,mm,m,zw(1,1,inzeep),z(1,jj,3-inzee),ntrig,trig,after(i),now(i),before(i),isign)
         end do
      endif
      inzee=3-inzee

   endif

END SUBROUTINE fft2d
