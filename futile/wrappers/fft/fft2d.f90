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


!!$!>two dimensional full convolution with a kernel given in the reciprocal space
!!$subroutine 2D_fullconv
!!$!$omp do schedule(static)
!!$do j3 = 1, maxIter
!!$   !this condition ensures that we manage only the interesting part for the FFT
!!$   !if (iproc*(nd3/nproc)+j3 <= n3/2+1) then
!!$   Jp2stb=1
!!$   J2stb=1
!!$   Jp2stf=1
!!$   J2stf=1
!!$   ! transform along x axis
!!$   lot=ncache/(4*n1)
!!$   if (lot < 1)  stop 'n1'
!!$
!!$   do j=1,n2dimp/n3pr1,lot
!!$      nfft=min(j+(lot-1), n2dimp/n3pr1) -j +1
!!$
!!$      !input: I1,J2,j3,Jp2,(jp3)
!!$      if (nproc > 1) then
!!$         call fft_parallel_block(n1,md2/(n3pr1*n3pr2),n3pr2,&
!!$              n1dim,md2/(n3pr1*n3pr2),nd3/n3pr2,&
!!$              nfft,lot,lzt/n3pr1,J2stb,Jp2stb,j3,&
!!$              ntrig,btrig1,after1,now1,before1,ic1,1,&
!!$              zmpi1,zw(1,1,1,ithread),zt(1,j,1,ithread))
!!$      else
!!$         call fft_parallel_block(n1,md2,1,&
!!$              n1dim,md2,nd3,&
!!$              nfft,lot,lzt/n3pr1,J2stb,Jp2stb,j3,&
!!$              ntrig,btrig1,after1,now1,before1,ic1,1,&
!!$              zmpi2,zw(1,1,1,ithread),zt(1,j,1,ithread))
!!$      end if
!!$      !output: I2,i1,j3,(jp3)
!!$   end do
!!$
!!$   !transform along y axis
!!$   lot=ncache/(4*n2)
!!$   if (lot < 1) stop 'n2'
!!$
!!$   do j=1,n1p/n3pr1,lot
!!$      i_one=1
!!$      nfft=min(j+(lot-1),n1p/n3pr1)-j+1
!!$      !reverse ordering
!!$      !input: I2,i1,j3,(jp3)
!!$      call fft_parallel_block(n2,n1,1,&
!!$           n2dim,n1,1,&
!!$           nfft,lot,lot,one,one,1,&
!!$           ntrig,btrig2,after2,now2,before2,ic2,1,&
!!$           zt(1,1,j,ithread),zw(1,1,1,ithread),zw(1,1,1,ithread),inzee)
!!$      we should insert inzee and eliminate zout in this case
!!$      !output: i1,i2,j3,(jp3)
!!$      !Multiply with kernel in fourier space
!!$      i3=mod(iproc,n3pr2)*(nd3/n3pr2)+j3
!!$
!!$      j1start=0
!!$      if (n3pr1>1) j1start=(n1p/n3pr1)*iproc_inplane
!!$
!!$      if (geocode == 'P') then
!!$         call P_multkernel(nd1,nd2,n1,n2,n3,lot,nfft,j+j1start,pot(1,1,j3),zw(1,1,inzee,ithread),&
!!$              i3,hx,hy,hz,offset,scal,strten_omp)
!!$      else
!!$         !write(*,*) 'pot(1,1,j3) = ', pot(1,1,j3)
!!$         call multkernel(nd1,nd2,n1,n2,lot,nfft,j+j1start,pot(1,1,j3),zw(1,1,inzee,ithread))
!!$      end if
!!$
!!$      !TRANSFORM BACK IN REAL SPACE
!!$      !transform along y axis
!!$      !input: i1,i2,j3,(jp3)
!!$      do i=1,ic2
!!$         call fftstp_sg(lot,nfft,n2,lot,n2,zw(1,1,inzee,ithread),zw(1,1,3-inzee,ithread),&
!!$              ntrig,ftrig2,after2(i),now2(i),before2(i),-1)
!!$         !zw(:,:,3-inzee)=zw(:,:,inzee)
!!$         inzee=3-inzee
!!$      end do
!!$
!!$      !reverse ordering
!!$      !input: i1,I2,j3,(jp3)
!!$      if (n3pr1 == 1) then
!!$         call G_unswitch_downcorn(nfft,n2,n2dim,lot,n1,lzt, &
!!$              zw(1,1,inzee,ithread),zt(1,1,j,ithread))
!!$      else
!!$         call G_unswitch_downcorn2(nfft,n2,n2dim,lot,n1p,lzt, &
!!$              zw(1,1,inzee,ithread),zt_t(1,1,1),n3pr1,j)
!!$      endif
!!$      !output: I2,i1,j3,(jp3)
!!$   end do
!!$   !transform along x axis
!!$   !input: I2,i1,j3,(jp3)
!!$
!!$   !LG: this MPI_ALLTOALL is inside a loop. I think that it will rapidly become unoptimal
!!$   if (n3pr1 > 1 .and. inplane_comm/=MPI_COMM_NULL) then
!!$      call f_timing(TCAT_PSOLV_COMPUT,'OF')
!!$      call f_timing(TCAT_PSOLV_COMMUN,'ON')
!!$
!!$      call MPI_ALLTOALL(zt_t,2*(n1p/n3pr1)*(lzt/n3pr1),MPI_double_precision,&
!!$           zt(1,1,1,ithread),2*(n1p/n3pr1)*(lzt/n3pr1), &
!!$           MPI_double_precision,inplane_comm,ierr)
!!$
!!$      call f_timing(TCAT_PSOLV_COMMUN,'OF')
!!$      call f_timing(TCAT_PSOLV_COMPUT,'ON')
!!$   endif
!!$
!!$   lot=ncache/(4*n1)
!!$   do j=1,n2dimp/n3pr1,lot
!!$      nfft=min(j+(lot-1),n2dimp/n3pr1)-j+1
!!$
!!$      !performing FFT
!!$      i=1
!!$      if (n3pr1 > 1) then
!!$         call fftstp_sg(lzt/n3pr1,nfft,n1,lot,n1,zt(1,j,1,ithread),zw(1,1,1,ithread),&
!!$              ntrig,ftrig1,after1(i),now1(i),before1(i),-1)
!!$      else
!!$         call fftstp_sg(lzt,nfft,n1,lot,n1,zt(1,j,1,ithread),zw(1,1,1,ithread),&
!!$              ntrig,ftrig1,after1(i),now1(i),before1(i),-1)
!!$      endif
!!$      inzee=1
!!$      do i=2,ic1
!!$         call fftstp_sg(lot,nfft,n1,lot,n1,zw(1,1,inzee,ithread),zw(1,1,3-inzee,ithread),&
!!$              ntrig,ftrig1,after1(i),now1(i),before1(i),-1)
!!$         inzee=3-inzee
!!$      enddo
!!$
!!$      !output: I2,I1,j3,(jp3)
!!$      !reverse ordering
!!$      !input: J2,Jp2,I1,j3,(jp3)
!!$      if (nproc == 1) then
!!$         call G_unmpiswitch_downcorn(j3,nfft,Jp2stf,J2stf,lot,n1,&
!!$              n1dim,md2,nd3,nproc,zw(1,1,inzee,ithread),zmpi2)
!!$      else
!!$         call G_unmpiswitch_downcorn2(j3,nfft,Jp2stf,J2stf,lot,n1,&
!!$              n1dim,md2,nd3,n3pr1,n3pr2,zw(1,1,inzee,ithread),zmpi1)  
!!$      endif
!!$      ! output: I1,J2,j3,Jp2,(jp3)
!!$   end do
!!$   !endif
!!$
!!$   !END OF TRANSFORM FOR X AND Z
!!$
!!$end do
!!$!$omp end do

subroutine fill_second_part(i1s,ld1,n1,npad,zw,zin)
  use f_precisions, only: dp => f_double
  implicit none
  integer, intent(in) :: ld1,n1,i1s,npad
  real(dp),intent(inout) ::  zw(2,ld1,*)
  real(dp), intent(in) :: zin(2,*)
  !local variables
  integer :: I

  do I=1,n1
     zw(1,i1s,I+npad)=zin(1,I)
     zw(2,i1s,I+npad)=zin(2,I)
  end do

end subroutine fill_second_part

subroutine get_first_part(i1s,ld1,n1,zw,zout)
  use f_precisions, only: dp => f_double
  implicit none
  integer, intent(in) :: ld1,n1,i1s
  real(dp),intent(in) ::  zw(2,ld1,*)
  real(dp), intent(inout) :: zout(2,*)
  !local variables
  integer :: i

  do i=1,n1
     zout(1,i)=zw(1,i1s,i)
     zout(2,i)=zw(2,i1s,i)
  end do
  
end subroutine get_first_part

subroutine transpose_output(ndims,lds,n1dim,lot,nfft,j3,&
     i2s,i2ps,zw,zout)
  use f_precisions, only: dp => f_double
  implicit none
  integer, intent(in) :: n1dim,lot,nfft,j3
  integer, dimension(3), intent(in) :: ndims,lds
  real(dp), dimension(*), intent(in) ::  zw
  real(dp), dimension(2,lds(1),lds(2),lds(3),*), intent(inout) ::  zout
  integer, intent(inout) :: i2s,i2ps
  !local variables
  integer :: j2,mfft,jp2,n1,n2,n2p

  n1=ndims(1)
  n2=ndims(2)
  n2p=ndims(3)
  
  mfft=0
  do Jp2=i2ps,n2p
     do J2=i2s,n2
        mfft=mfft+1
        if (mfft > nfft) then
           i2ps=Jp2
           i2s=J2
           return
        end if
        call get_first_part(mfft,lot,n1dim,zw,zout(1,1,J2,j3,Jp2))
     end do
     i2s=1
  end do

end subroutine transpose_output

subroutine transpose_and_pad_input(ndims,lds,n1dim,lot,nfft,j3,&
     i2s,i2ps,zin,zw)
  use f_precisions, only: dp => f_double
  implicit none
  integer, intent(in) :: n1dim,lot,nfft,j3
  integer, dimension(3), intent(in) :: ndims,lds
  real(dp), dimension(2,lds(1),lds(2),lds(3),*), intent(in) ::  zin
  integer, intent(inout) :: i2s,i2ps
  real(dp), dimension(2,lot,*), intent(inout) ::  zw
  !local variables
  integer :: j2,mfft,jp2,ish,i1,imfft,n1,n2,n2p

  n1=ndims(1)
  n2=ndims(2)
  n2p=ndims(3)

  ish=n1-n1dim
  mfft=0
  loop_Jp2: do Jp2=i2ps,n2p
     do J2=i2s,n2
        mfft=mfft+1
        if (mfft > nfft) then
           i2ps=Jp2
           i2s=J2
           exit loop_Jp2
        end if
        call fill_second_part(mfft,lot,n1dim,ish,&
             zw,zin(1,1,J2,j3,Jp2))
     end do
     i2s=1
  end do loop_Jp2

  !then zero the rest of the cache array if the padding is required
  do i1=1,ish
     do imfft=1,mfft-1
        zw(1,imfft,i1)=0.0_dp
        zw(2,imfft,i1)=0.0_dp
     end do
  end do

end subroutine transpose_and_pad_input


!!$subroutine fft_1d_cache(ndat_in,ld_in,ndat_out,ld_out,nfft,ninout,n,&
!!$     ncache,ntrig,trig,after,now,before,ic,&
!!$     i_sign,inzee,transpose,iam,nthread,z,zw)
!!$  use f_precisions
!!$  use module_fft_sg, only: n_factors
!!$  implicit none
!!$  logical, intent(in) :: transpose
!!$  integer, intent(in) :: ndat_in,ld_in,ndat_out,ld_out,nfft,ninout,n
!!$  integer, intent(in) :: ncache,ntrig,ic
!!$  integer, intent(in) :: i_sign,iam,nthread
!!$  integer, intent(inout) :: inzee
!!$  integer, dimension(n_factors) :: after,now,before
!!$  real(f_double), dimension(2,ntrig), intent(in) :: trig
!!$  real(f_double), dimension(2,ninout,2), intent(inout) :: z
!!$  real(f_double), dimension(2,ncache/4,2), intent(inout) :: zw
!!$  !local variables
!!$  logical :: real_input = .false.
!!$  integer :: i,lotomp,ma,mb,nfft_th,j,jj,jompa,jompb,inzeep,inzet,lot,nn
!!$
!!$  !set of fft to be treated by the present thread
!!$  lotomp=nfft/nthread+1
!!$  jompa=iam*lotomp+1
!!$  jompb=min((iam+1)*lotomp,nfft)
!!$
!!$  !here we might put the treatment for the input and the output
!!$  !in case they are necessary
!!$
!!$  inzet=inzee
!!$  if (ic == 1 .or. ncache == 0) then
!!$     i=ic
!!$     ma=jompa
!!$     mb=jompb
!!$     nfft_th=mb-ma+1
!!$     j=ma
!!$     if (transpose) then
!!$        jj=j*ld_out-ld_out+1
!!$        !input z(2,ndat_in,ld_in,inzet)
!!$        !output z(2,ld_out,ndat_out,3-inzet)          
!!$        call fftrot_sg(ndat_in,nfft_th,ld_in,ndat_out,ld_out,&
!!$             z(1,j,inzet),z(1,jj,3-inzet), &
!!$             ntrig,trig,after(i),now(i),before(i),i_sign)
!!$        if (real_input) then !only works when ld_out==n
!!$           inzet=3-inzet
!!$           call unpack_rfft(ndat_out,n,&
!!$                z(1,jj,inzet),z(1,jj,3-inzet))
!!$        end if
!!$     else
!!$        !input z(2,ndat_in,ld_in,inzet)
!!$        !output z(2,ndat_out,ld_out,3-inzet)
!!$        call fftstp_sg(ndat_in,nfft_th,ld_in,ndat_out,ld_out,&
!!$             z(1,j,inzet),z(1,j,3-inzet), &
!!$             ntrig,trig,after(i),now(i),before(i),i_sign)
!!$        if (real_input) then !only works when ld_out==n
!!$           inzet=3-inzet
!!$           !here maybe the ndat_out has to be rethought
!!$           call unpack_rfft_t((ndat_out-1)/2+1,ndat_out,n,&
!!$                z(1,j,inzet),z(1,j,3-inzet))
!!$        end if
!!$     end if
!!$  else
!!$     lot=max(1,ncache/(4*n))
!!$     nn=lot
!!$     do j=jompa,jompb,lot
!!$        ma=j
!!$        mb=min(j+(lot-1),jompb)
!!$        nfft_th=mb-ma+1
!!$        i=1
!!$        inzeep=2
!!$        call fftstp_sg(ndat_in,nfft_th,ld_in,nn,n,&
!!$             z(1,j,inzet),zw(1,1,3-inzeep), &
!!$             ntrig,trig,after(i),now(i),before(i),i_sign)
!!$        inzeep=1
!!$        do i=2,ic-1
!!$           call fftstp_sg(nn,nfft_th,n,nn,n,&
!!$                zw(1,1,inzeep),zw(1,1,3-inzeep), &
!!$                ntrig,trig,after(i),now(i),before(i),i_sign)
!!$           inzeep=3-inzeep
!!$        end do
!!$        i=ic
!!$        if (transpose) then
!!$           jj=j*ld_out-ld_out+1
!!$           call fftrot_sg(nn,nfft_th,n,ndat_out,ld_out,&
!!$                zw(1,1,inzeep),z(1,jj,3-inzet), &
!!$                ntrig,trig,after(i),now(i),before(i),i_sign)
!!$           if (real_input) then !only works when ld_out==n
!!$              inzet=3-inzet
!!$              call unpack_rfft(ndat_out,n,&
!!$                   z(1,jj,inzet),z(1,jj,3-inzet))
!!$           end if
!!$        else
!!$           call fftstp_sg(nn,nfft_th,n,ndat_out,ld_out,&
!!$                zw(1,1,inzeep),z(1,j,3-inzet), &
!!$                ntrig,trig,after(i),now(i),before(i),i_sign)
!!$           if (real_input) then !only works when ld_out==n
!!$              inzet=3-inzet
!!$              call unpack_rfft_t((ndat_out-1)/2+1,ndat_out,n,&
!!$                   z(1,j,inzet),z(1,j,3-inzet))
!!$           end if
!!$        end if
!!$     end do
!!$  end if
!!$  inzet=3-inzet
!!$  if (iam==0) inzee=inzet !as it is a shared variable
!!$
!!$end subroutine fft_1d_base
!!$
!!$!> base wrapper for FFT, not to be called in public routines
!!$!! @warning: for transa='t' nwork must be > 0
!!$subroutine fft_1d_base(transa,transb,ndat_in,ld_in,ndat_out,ld_out,ninout,nfft,n,&
!!$     ntrig,trig,after,now,before,ic,i_sign,&
!!$     input_provided,za,&
!!$     output_provided,zb,&
!!$     inout_provided,zab,inzee,&
!!$     nwork,zw,iam,nthread)
!!$  use f_precisions
!!$  use module_fft_sg, only: n_factors
!!$  implicit none
!!$  integer, intent(in) :: ndat_in,ld_in,ndat_out,ld_out,nfft,ninout,n
!!$  integer, intent(in) :: nwork,ntrig,ic
!!$  integer, intent(in) :: i_sign,iam,nthread
!!$  character(len=1), intent(in) :: transa,transb
!!$  integer, intent(inout) :: inzee
!!$  integer, dimension(n_factors) :: after,now,before
!!$  real(f_double), dimension(2,ntrig), intent(in) :: trig
!!$  real(f_double), dimension(2,ndat_in*ld_in), intent(in) :: za
!!$  real(f_double), dimension(2,ndata_out*ld_out), intent(inout) :: zb
!!$  real(f_double), dimension(2,ninout,2), intent(inout) :: zab
!!$  real(f_double), dimension(2,nwork/4,2), intent(inout) :: zw
!!$  !local variables
!!$  logical :: real_input = .false.,transpose_output,transpose_input,nowork
!!$  integer :: i,lotomp,ma,mb,nfft_th,j,jj,jompa,jompb,inzeep,inzet,lot,nn
!!$
!!$  !set of fft to be treated by the present thread
!!$  transpose_input=transa=='T' .or. transa=='t'
!!$  transpose_output=transb=='T' .or. transb=='t'
!!$
!!$  !check on arguments
!!$  call f_assert(input_provided .and. (output_provided .or. nwork>0),id='Error on input work')
!!$  call f_assert(transpose_input .and. input_provided .and. nwork >0,id='Transposition needs workarray')
!!$
!!$
!!$  lotomp=nfft/nthread+1
!!$  jompa=iam*lotomp+1
!!$  jompb=min((iam+1)*lotomp,nfft)
!!$  nowork= ic == 1 .or. nwork == 0 !work array cannot be used for only one step
!!$  if (nowork) then
!!$     lot=jompb-jompa 
!!$  else
!!$     lot=max(1,nwork/(4*n))
!!$  end if
!!$
!!$  !here we might put the treatment for the input and the output
!!$  !in case they are necessary
!!$
!!$  inzet=inzee
!!$  loop_on_lines: do j=jompa,jompb,lot !this loop is internal to the threads, it might be further parallelised
!!$     ma=j
!!$     mb=min(j+(lot-1),jompb)
!!$     nfft_th=mb-ma+1
!!$
!!$     !then it follows the treatments for the different cases
!!$     if (nowork) then
!!$        if (
!!$        !we need that z is large enough to contain max(ndat_in*ld_in,ndat_out*ld_out)
!!$        do i=1,ic-1
!!$           call fftstp_sg(ndat_in,nfft_th,ld_in,ndat_out,ld_out,&
!!$                zab(1,j,inzet),zab(1,j,3-inzet), &
!!$                ntrig,trig,after(i),now(i),before(i),i_sign)
!!$           inzet=3-inzet
!!$        end do
!!$        i=ic    
!!$        if (transpose_output) then
!!$           jj=j*ld_out-ld_out+1
!!$           !input z(2,ndat_in,ld_in,inzet)
!!$           !output z(2,ld_out,ndat_out,3-inzet)          
!!$           call fftrot_sg(ndat_in,nfft_th,ld_in,ndat_out,ld_out,&
!!$                zab(1,j,inzet),zab(1,jj,3-inzet), &
!!$                ntrig,trig,after(i),now(i),before(i),i_sign)
!!$           if (real_input) then !only works when ld_out==n
!!$              inzet=3-inzet
!!$              call unpack_rfft(ndat_out,n,&
!!$                   z(1,jj,inzet),z(1,jj,3-inzet))
!!$           end if
!!$        else
!!$           !input z(2,ndat_in,ld_in,inzet)
!!$           !output z(2,ndat_out,ld_out,3-inzet)
!!$           call fftstp_sg(ndat_in,nfft_th,ld_in,ndat_out,ld_out,&
!!$                zab(1,j,inzet),zab(1,j,3-inzet), &
!!$                ntrig,trig,after(i),now(i),before(i),i_sign)
!!$           if (real_input) then !only works when ld_out==n
!!$              inzet=3-inzet
!!$              !here maybe the ndat_out has to be rethought
!!$              call unpack_rfft_t((ndat_out-1)/2+1,ndat_out,n,&
!!$                   z(1,j,inzet),z(1,j,3-inzet))
!!$           end if
!!$        end if
!!$     else
!!$        nn=lot
!!$        i=1
!!$        inzeep=2
!!$        if (input_provided) then
!!$           call fftstp_sg(ndat_in,nfft_th,ld_in,nn,n,&
!!$                za,zw(1,1,3-inzeep), &
!!$                ntrig,trig,after(i),now(i),before(i),i_sign)
!!$        else if (inout_provided) then
!!$           call fftstp_sg(ndat_in,nfft_th,ld_in,nn,n,&
!!$                zab(1,j,inzet),zw(1,1,3-inzeep), &
!!$                ntrig,trig,after(i),now(i),before(i),i_sign)
!!$        end if
!!$        inzeep=1
!!$        do i=2,ic-1
!!$           call fftstp_sg(nn,nfft_th,n,nn,n,&
!!$                zw(1,1,inzeep),zw(1,1,3-inzeep), &
!!$                ntrig,trig,after(i),now(i),before(i),i_sign)
!!$           inzeep=3-inzeep
!!$        end do
!!$        i=ic
!!$        if (transpose_output) then
!!$           jj=j*ld_out-ld_out+1
!!$           call fftrot_sg(nn,nfft_th,n,ndat_out,ld_out,&
!!$                zw(1,1,inzeep),zab(1,jj,3-inzet), &
!!$                ntrig,trig,after(i),now(i),before(i),i_sign)
!!$           if (real_input) then !only works when ld_out==n
!!$              inzet=3-inzet
!!$              call unpack_rfft(ndat_out,n,&
!!$                   z(1,jj,inzet),z(1,jj,3-inzet))
!!$           end if
!!$        else
!!$           call fftstp_sg(nn,nfft_th,n,ndat_out,ld_out,&
!!$                zw(1,1,inzeep),zab(1,j,3-inzet), &
!!$                ntrig,trig,after(i),now(i),before(i),i_sign)
!!$           if (real_input) then !only works when ld_out==n
!!$              inzet=3-inzet
!!$              call unpack_rfft_t((ndat_out-1)/2+1,ndat_out,n,&
!!$                   z(1,j,inzet),z(1,j,3-inzet))
!!$           end if
!!$        end if
!!$     end if
!!$  end do loop_on_lines
!!$  inzet=3-inzet
!!$  if (iam==0) inzee=inzet !as it is a shared variable
!!$
!!$end subroutine fft_1d_base
