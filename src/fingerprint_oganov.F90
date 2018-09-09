module parameters
     implicit none
     save
     integer, allocatable:: Kinds(:)
     real(8), allocatable:: rxyz(:,:)
     real(8):: rcut(2,2)
end module parameters

subroutine fingerprint_oganov(parini)
use mod_parini, only: typ_parini
use parameters, only: kinds,rxyz
implicit none
type(typ_parini), intent(in):: parini
integer:: i,j,k,l,m,iat,jat,nat,nkinds,nec,nexp,fp_size,i_bin,near_bin,n_peaks,fp_dim
integer, allocatable:: nkinds_sum1(:),nkinds_sum2(:),kinds1(:),kinds2(:)
real(8):: latvec(3,3)
real(8), allocatable:: rxyzexp(:,:,:,:,:),transvecall(:,:,:,:),fp1(:),fp2(:)
real(8):: sigma,r_cut,d_bin,vol,r,value,val_all,rij
real(8):: phi(600),floor_x,ceil_x,int_val,x,distance
real(8), parameter :: pi=3.141592653589793238462643383279502884197d0
character(40):: filename1,filename2
integer:: n1,n2,n3
integer, allocatable:: ggt_0(:,:,:)
!Set parameters
sigma=0.02d0
r_cut=30.00d0
d_bin=0.05d0
write(*,*) "How many atom kinds?"
read(*,*)nkinds
!!if(nkinds.ne.2 .and. nkinds.ne.1) stop "Wrong number of atom kinds"

!n1=20;n2=20;n3=20
!allocate(ggt_0(n1,n2,n3))
!call create_ggtlist(n1,n2,n3,ggt_0)


!!do i=1,5
!!  do j=1,i
!!  write(*,*) i,j,(i*(i+1)/2+j-i)
!!  enddo
!!enddo
!!write(*,*)
!!do k=1,15
!!   i=ceiling(-0.5d0+0.5d0*sqrt(1.d0+8.d0*real(k,8)))
!!   j=(i*(1-i))/2+k
!!   write(*,*) i,j,k
!!enddo
!!stop



!Allocate fingerprint array
fp_size=ceiling(r_cut/d_bin)
fp_dim=nkinds*(nkinds+1)/2
write(*,*) fp_size,fp_dim

!Read input file
write(*,*) "Reading input file..."
!filename1="run002.ascii"
read(*,*) filename1
call readfile(nat,latvec,nkinds,filename1)
!!if(nkinds==1) then 
!!call readfile_single(nat,latvec,nkinds,filename1)
!!elseif(nkinds==2) then
!!call readfile_double(nat,latvec,nkinds,filename1)
!!else
!!stop "Wrong number NKINDS"
!!endif
allocate(nkinds_sum1(nkinds))
nkinds_sum1=0
do iat=1,nat
  do i=1,nkinds
  if(Kinds(iat)==i) nkinds_sum1(i)= nkinds_sum1(i)+1
  enddo
enddo

!Check if the arrays have been allocated properly
if(.not.allocated(rxyz)) stop "Initial allocation failed,RXYZ"
if(.not.allocated(kinds)) stop "Initial allocation failed,KINDS"

!Allocate the fingerprint funktion array and normalization matrix
allocate(fp1(fp_size*fp_dim),kinds1(nat))
kinds1=kinds
call get_fp_oganov(nat,rxyz,latvec,r_cut,sigma,d_bin,kinds1,nkinds,nkinds_sum1,fp_size,fp_dim,fp1)
deallocate(rxyz,kinds)


!Read input file
write(*,*) "Reading input file..."
!filename2="run003.ascii"
read(*,*) filename2
call readfile(nat,latvec,nkinds,filename2)
!!if(nkinds==1) then 
!!call readfile_single(nat,latvec,nkinds,filename2)
!!elseif(nkinds==2) then
!!call readfile_double(nat,latvec,nkinds,filename2)
!!else
!!stop "Wrong number NKINDS"
!!endif
allocate(nkinds_sum2(nkinds))
nkinds_sum2=0
do iat=1,nat
  do i=1,nkinds
  if(Kinds(iat)==i) nkinds_sum2(i)= nkinds_sum2(i)+1
  enddo
enddo

!Check if the arrays have been allocated properly
if(.not.allocated(rxyz)) stop "Initial allocation failed,RXYZ"
if(.not.allocated(kinds)) stop "Initial allocation failed,KINDS"

!Allocate the fingerprint funktion array and normalization matrix
allocate(fp2(fp_size*fp_dim),kinds2(nat))
kinds2=kinds
call get_fp_oganov(nat,rxyz,latvec,r_cut,sigma,d_bin,kinds2,nkinds,nkinds_sum2,fp_size,fp_dim,fp2)
deallocate(rxyz,kinds)


do i=1,fp_size
write(11,*) (real(i,8)-1.d0)*d_bin,fp1(i)
write(12,*) (real(i,8)-1.d0)*d_bin,fp1(fp_size+i)
write(13,*) (real(i,8)-1.d0)*d_bin,fp1(2*fp_size+i)
write(21,*) (real(i,8)-1.d0)*d_bin,fp2(i)
write(22,*) (real(i,8)-1.d0)*d_bin,fp2(fp_size+i)
write(23,*) (real(i,8)-1.d0)*d_bin,fp2(2*fp_size+i)
!!write(14,*) (real(i,8)-1.d0)*d_bin,fp1(i,2,2)
enddo
!
!
!
!!do i=1,fp_size*fp_dim
!!write(21,*) (real(i,8)-1.d0)*d_bin,fp2(i)
!!!write(22,*) (real(i,8)-1.d0)*d_bin,fp2(i,1,2)
!!!write(23,*) (real(i,8)-1.d0)*d_bin,fp2(i,2,1)
!!!write(24,*) (real(i,8)-1.d0)*d_bin,fp2(i,2,2)
!!enddo


do i=1,nkinds
if(nkinds_sum1(i).ne.nkinds_sum2(i)) write(*,*) "WARNING: The number of atoms A and B do not coincide"
enddo
call get_cosinedistance(parini,fp1,fp2,fp_size,fp_dim,nkinds,nkinds_sum1,distance)
write(*,*) distance
end subroutine

subroutine get_fp_oganov(nat,rxyz,latvec,r_cut,sigma,d_bin,kinds,nkinds,nkinds_sum,fp_size,fp_dim,fp)
!This subroutine will compute the F-fingeprint function discretized into bins of size d_bin with 
!smoothed out delta functions with sigma and generate the output fp. All conventions and 
!methods are from J.Chem.Phys, 130, 104504 (2009) and IEEE Symposium, Okt 21-23. (2008) (M.Valle and A.Oganov)
use yaml_output
implicit none
integer:: i,j,k,l,m,iat,jat,nat,nec,nec1,nec2,nec3,nexp,fp_size,i_bin,near_bin,n_peaks,fp_dim,iarr
integer:: nkinds_sum(nkinds),kinds(nat),nkinds,imin,imax
real(8):: latvec(3,3),rxyz(3,nat),fp_norm(fp_dim),fp(fp_size*fp_dim)
real(8), allocatable:: rxyzexp(:,:,:,:,:),transvecall(:,:,:,:)
real(8):: sigma,r_cut,d_bin,vol,r,value,val_all,rij
real(8):: phi(600),floor_x,ceil_x,int_val,x
real(8), parameter:: pi=acos(-1.d0)
real(8):: drx,dry,drz,t(2),lx,ly,lz,rmax,drx2,dry2,drz2,lx2,ly2,lz2,rmax2,l2,nom_f,discr,rxyzj(3),transvec(3)
integer, allocatable:: list(:,:)
integer:: ilist,nlist,t_it 
integer, allocatable:: ggt_0(:,:,:)

!Set up indexing table (cant i do it easier???)
!We could: k=i(i+1)/2+j-i
!!k=0
!!do i=1,ndkins
!!do j=1,i
!!   itable(i,j)=k
!!   itable(j,i)=itable(i,j)
!!   k=k+1
!!enddo
!!enddo



!Call gaussian integral parameters
call set_phi(phi)

!Calculate the volume of the cell
call getvol(latvec,vol)

!Create the expanded unit cells
!call n_rep(latvec,r_cut,nec)
!write(*,'(a,i1.1,a)') " Creating expansion ",nec,"..."
call n_rep_dim(latvec,r_cut+6.d0*sigma,nec1,nec2,nec3)
call yaml_mapping_open('Creating expansion for periodic',flow=.true.)
call yaml_map('method','Oganov')
call yaml_map('nec',(/nec1,nec2,nec3/))
call yaml_mapping_close()
!write(*,'(a,3(i3,1x),a)') " Creating expansion for periodic Oganov FP ",nec1,nec2,nec3,"..."
!allocate(rxyzexp(3,nat,nec1+1,nec2+1,nec3+1),transvecall(3,nec1+1,nec2+1,nec3+1))

allocate(rxyzexp(3,nat,2*nec1+1,2*nec2+1,2*nec3+1),transvecall(3,2*nec1+1,2*nec2+1,2*nec3+1))
!write(*,*) "Allocated"

call expand_dim(rxyz,rxyzexp,transvecall,latvec,nat,nec1,nec2,nec3)
!write(*,'(a,i3,a)') " Creating expansion ",nec,"..."
!write(15,*) rxyzexp
!write(16,*) transvecall 
!deallocate(rxyzexp,transvecall)


!    if(nec==1) then
!allocate(rxyzexp(3,nat,3,3,3),transvecall(3,3,3,3))
!   call expand(rxyz,rxyzexp,transvecall,latvec,nat)
!   nexp=3
!elseif(nec==2) then
!allocate(rxyzexp(3,nat,5,5,5),transvecall(3,5,5,5))
!   call expand_2(rxyz,rxyzexp,transvecall,latvec,nat)
!   nexp=5
!elseif(nec==3) then
!allocate(rxyzexp(3,nat,7,7,7),transvecall(3,7,7,7))
!   call expand_3(rxyz,rxyzexp,transvecall,latvec,nat)
!   nexp=7
!elseif(nec==4) then
!allocate(rxyzexp(3,nat,9,9,9),transvecall(3,9,9,9))
!   call expand_4(rxyz,rxyzexp,transvecall,latvec,nat)
!   nexp=9
!elseif(nec==5) then
!allocate(rxyzexp(3,nat,11,11,11),transvecall(3,11,11,11))
!   call expand_5(rxyz,rxyzexp,transvecall,latvec,nat)
!   nexp=11
!elseif(nec==6) then
!allocate(rxyzexp(3,nat,13,13,13),transvecall(3,13,13,13))
!   call expand_6(rxyz,rxyzexp,transvecall,latvec,nat)
!   nexp=13
!elseif(nec==7) then
!allocate(rxyzexp(3,nat,15,15,15),transvecall(3,15,15,15))
!   call expand_7(rxyz,rxyzexp,transvecall,latvec,nat)
!   nexp=15
!elseif(nec==8) then
!allocate(rxyzexp(3,nat,17,17,17),transvecall(3,17,17,17))
!   call expand_8(rxyz,rxyzexp,transvecall,latvec,nat)
!   nexp=17
!elseif(nec==9) then
!allocate(rxyzexp(3,nat,19,19,19),transvecall(3,19,19,19))
!   call  expand_9(rxyz,rxyzexp,transvecall,latvec,nat)
!   nexp=19
!else
!   write(*,*) "r_cut is too large or a wrong value, reduce it!"
!   stop
!endif

!write(25,*) rxyzexp(:,:,nec+1:nexp,nec+1:nexp,nec+1:nexp)
!write(26,*) transvecall(:,nec+1:nexp,nec+1:nexp,nec+1:nexp)








rmax=r_cut
rmax2=r_cut*r_cut
n_peaks=0.d0
fp=0.d0

do iat=1,nat
  do jat=1,iat
     do k=1,2*nec1+1
     do l=1,2*nec2+1
     do m=1,2*nec3+1

      if(iat==jat.and.k==nec1+1.and.l==nec2+1.and.m==nec3+1) goto 2020
!      if(iat==jat) goto 2020
        
        rxyzj(:)=rxyz(:,jat)+transvecall(:,k,l,m)
        rij=sqrt((rxyz(1,iat)-rxyzj(1))**2+(rxyz(2,iat)-rxyzj(2))**2+(rxyz(3,iat)-rxyzj(3))**2)
        if(rij==0.d0) write(*,*) k,l,m,iat,jat,nec1+1,nec2+1,nec3+1
        if(rij.lt.r_cut+sigma*6.d0) then
!        if(rij.lt.r_cut) then
        n_peaks=n_peaks+1
        if(jat.ne.iat)        n_peaks=n_peaks+1
        near_bin=int(rij/d_bin)+1        
        do i_bin=max(1,near_bin-int(6.d0*sigma/d_bin)-1),min(fp_size,near_bin+int(6.d0*sigma/d_bin)+2)
!        do i_bin=1,fp_size
        floor_x=(real(i_bin,8)-1.d0)*d_bin
        ceil_x=(real(i_bin,8))*d_bin
        call integrate_gaussian(floor_x,ceil_x,rij,sigma,phi,int_val)
!        call integrate_gaussian_exact(floor_x,ceil_x,rij,sigma,int_val)
        imax=max(Kinds(iat),Kinds(jat))
        imin=min(Kinds(iat),Kinds(jat))
        iarr=imax*(imax+1)
        iarr=i_bin+(iarr/2+imin-imax-1)*fp_size
        if(iarr.gt.fp_size*fp_dim.or.iarr.le.0) stop "Something wrong!!!"
        fp(iarr)=fp(iarr)+int_val/rij**2
        if(iat.ne.jat) fp(iarr)=fp(iarr)+int_val/rij**2*real(1/(abs(imin-imax)+1),8)
        enddo
        endif
      2020 continue
     enddo
     enddo
     enddo


  enddo
enddo




!rmax=r_cut
!rmax2=r_cut*r_cut
!n_peaks=0.d0
!fp=0.d0
!allocate(ggt_0(nec1,nec2,nec3),list(3,nec1*nec2*nec3))
!call create_ggtlist(nec1,nec2,nec3,ggt_0,list,nlist)
!
!do iat=1,nat
!  do jat=1,nat
!     drx=rxyz(1,jat)-rxyz(1,iat);dry=rxyz(2,jat)-rxyz(2,iat);drz=rxyz(3,jat)-rxyz(3,iat)
!     drx2=drx*drx;dry2=dry*dry;drz2=drz*drz
!!     do k=1,nec1+1
!!     do l=1,nec2+1
!!     do m=1,nec3+1
!     do ilist=1,nlist
!      do k=1,2
!      do l=1,2
!      do m=1,2
!
!
!      transvec(:)=list(1,ilist)*latvec(:,1)*(-1)**k+list(2,ilist)*latvec(:,2)*(-1)**l+list(3,ilist)*latvec(:,3)*(-1)**m
!      
!      lx=transvec(1);ly=transvec(2);lz=transvec(3)
!      lx2=lx*lx;ly2=ly*ly;lz2=lz*lz
!      nom_f=drx*lx+dry*ly+drz*lz
!      l2=lx2+ly2+lz2
!      discr=sqrt(l2*rmax2-dry2*lx2-drz2*lx2-drx2*ly2-drz2*ly2-drx2*lz2-dry2*lz2+drx*dry*lx*ly*2+drx*drz*lx*lz*2+dry*drz*ly*lz*2)
!     
!      T(1) = -(nom_f-discr)/l2
!      T(2) = -(nom_f+discr)/l2
!!      write(*,*) "T1",T(1)
!!      write(*,*) "T2",T(2)
!    
!      do t_it=int(T(2)),int(T(1))
!        if(t_it==0.and.iat==jat) goto 2020     
!!      if(iat==jat.and.k==int(nexp*0.5d0)+1.and.l==int(nexp*0.5d0)+1.and.l==int(nexp*0.5d0)+1) goto 2020
!        rxyzj(:)=rxyz(:,jat)+t_it*transvec(:)
!        rij=sqrt((rxyz(1,iat)-rxyzj(1))**2+(rxyz(2,iat)-rxyzj(2))**2+(rxyz(3,iat)-rxyzj(3))**2)
!!        if(rij.lt.r_cut-sigma*5.d0) then
!!        if(rij.lt.r_cut) then
!        n_peaks=n_peaks+1
!        near_bin=int(rij/d_bin)+1        
!        do i_bin=max(1,near_bin-int(6.d0*sigma/d_bin)-1),min(fp_size,near_bin+int(6.d0*sigma/d_bin)+2)
!!        do i_bin=1,fp_size
!        floor_x=(real(i_bin,8)-1.d0)*d_bin
!        ceil_x=(real(i_bin,8))*d_bin
!        call integrate_gaussian(floor_x,ceil_x,rij,sigma,phi,int_val)
!        fp(i_bin,Kinds(iat),Kinds(jat))=fp(i_bin,Kinds(iat),Kinds(jat))+int_val/rij**2
!        enddo
!!        endif
!      2020 continue
!     enddo
!     enddo
!     enddo
!
!
!!    enddo
!    enddo
!    enddo
!  enddo
!enddo




!n_peaks=0.d0
!fp=0.d0
!do iat=1,nat
!  do jat=1,nat
!     do k=1,nexp
!     do l=1,nexp
!     do m=1,nexp
!      if(iat==jat.and.k==int(nexp*0.5d0)+1.and.l==int(nexp*0.5d0)+1.and.l==int(nexp*0.5d0)+1) goto 2020
!        rij=sqrt((rxyz(1,iat)-rxyzexp(1,jat,k,l,m))**2+(rxyz(2,iat)-rxyzexp(2,jat,k,l,m))**2+(rxyz(3,iat)-rxyzexp(3,jat,k,l,m))**2)
!        if(rij.lt.r_cut-sigma*5.d0) then
!!        if(rij.lt.r_cut) then
!        n_peaks=n_peaks+1
!        near_bin=int(rij/d_bin)+1        
!        do i_bin=max(1,near_bin-int(6.d0*sigma/d_bin)-1),min(fp_size,near_bin+int(6.d0*sigma/d_bin)+2)
!!        do i_bin=1,fp_size
!        floor_x=(real(i_bin,8)-1.d0)*d_bin
!        ceil_x=(real(i_bin,8))*d_bin
!        call integrate_gaussian(floor_x,ceil_x,rij,sigma,phi,int_val)
!        fp(i_bin,Kinds(iat),Kinds(jat))=fp(i_bin,Kinds(iat),Kinds(jat))+int_val/rij**2
!        enddo
!        endif
!      2020 continue
!    enddo
!    enddo
!    enddo
!  enddo
!enddo

!Calculate the Normalization constant
do k=1,fp_dim
!Stephan's clever formula...
   i=ceiling(-0.5d0+0.5d0*sqrt(1.d0+8.d0*real(k,8)))
   j=(i*(1-i))/2+k
!write(*,*) i,j
   fp_norm(k)=4.d0*pi*nkinds_sum(i)*nkinds_sum(j)*d_bin/vol 
   fp((k-1)*fp_size+1:k*fp_size)=fp((k-1)*fp_size+1:k*fp_size)/fp_norm(k)-1.d0
!write(*,*) fp_norm(k)
enddo

!!do i=1,nkinds
!!  do j=1,nkinds
!!  fp_norm(i,j)=4.d0*pi*nkinds_sum(i)*nkinds_sum(j)*d_bin/vol
!!  fp(:,i,j)=fp(:,i,j)/fp_norm(i,j)-1.d0
!!  enddo
!!enddo
deallocate(rxyzexp,transvecall)

end subroutine


subroutine get_cosinedistance(parini,fp1,fp2,fp_size,fp_dim,nkinds,nkinds_sum,distance)
use mod_parini, only: typ_parini
!This sunroutine will compute the cosine distance between two fingerprints fp1 and fp2 and 
!store the output distance.  All conventions and 
!methods are from J.Chem.Phys, 130, 104504 (2009) and IEEE Symposium, Okt 21-23. (2008) (M.Valle and A.Oganov)
implicit none
type(typ_parini), intent(in):: parini
integer, intent(IN) :: nkinds, nkinds_sum(nkinds), fp_size,fp_dim
!real(8), intent(IN) :: fp1(fp_size,nkinds,nkinds),fp2(fp_size,nkinds,nkinds)
real(8), intent(INOUT) :: fp1(fp_size*fp_dim),fp2(fp_size*fp_dim)
real(8), intent(OUT):: distance
integer             :: i_kind,j_kind,i_fp,i,j,k,imin,imax
real(8)             :: w_ab(fp_dim),w_norm,num,denom,tmp_1,tmp_2

!!TEST
!open(unit=2,file="plot.14")
!do i_fp=1,600
!  read(2,*)fp1(i_fp,1,1)
!enddo 
!close(2)
!open(unit=2,file="plot.13")
!do i_fp=1,600
!  read(2,*)fp1(i_fp,1,2)
!  fp1(i_fp,2,1)=fp1(i_fp,1,2)
!enddo 
!close(2)
!open(unit=2,file="plot.11")
!do i_fp=1,600
!  read(2,*)fp1(i_fp,2,2)
!enddo 
!close(2)
!!TEST
!open(unit=2,file="plot.24")
!do i_fp=1,600
!  read(2,*)fp2(i_fp,1,1)
!enddo 
!close(2)
!open(unit=2,file="plot.23")
!do i_fp=1,600
!  read(2,*)fp2(i_fp,1,2)
!  fp2(:,2,1)=fp2(:,1,2)
!enddo 
!close(2)
!open(unit=2,file="plot.21")
!do i_fp=1,600
!  read(2,*)fp2(i_fp,2,2)
!enddo 
!close(2)





!Compute the weight of each fingerprint
w_ab=0.d0
w_norm=0.d0
do k=1,fp_dim
   i=ceiling(-0.5d0+0.5d0*sqrt(1.d0+8.d0*real(k,8)))
   j=(i*(1-i))/2+k
   w_ab(k)=real(nkinds_sum(i),8)*real(nkinds_sum(j),8)
   w_norm=w_norm+real(nkinds_sum(i),8)*real(nkinds_sum(j),8)
enddo
!!do i=1,nkinds
!!  do j=1,i!nkinds
!!  w_ab(i,j)=real(nkinds_sum(i),8)*real(nkinds_sum(j),8)
!!  w_norm=w_norm+real(nkinds_sum(i),8)*real(nkinds_sum(j),8)
!!  enddo
!!enddo
w_ab=w_ab/w_norm
!write(*,*) w_ab

!Compute the distance
num=0.d0;  denom=0.d0;  tmp_1=0.d0;  tmp_2=0.d0

if(trim(parini%potential_potential).ne."lenosky_tb_lj") then
!Usual case: atoms in general
  do k=1,fp_dim
       imin=(k-1)*fp_size+1
       imax=k*fp_size
       num=num+dot_product(fp1(imin:imax),fp2(imin:imax))*w_ab(k)**2
       tmp_1=tmp_1+dot_product(fp1(imin:imax),fp1(imin:imax))*w_ab(k)**2
       tmp_2=tmp_2+dot_product(fp2(imin:imax),fp2(imin:imax))*w_ab(k)**2
  enddo
else
!Only for lenosky_tb_lj
  do k=1,fp_dim
!Stephan's formula to identify the LJ indexes...
       i=ceiling(-0.5d0+0.5d0*sqrt(1.d0+8.d0*real(k,8)))
       j=(i*(1-i))/2+k
       if(.not.(parini%znucl(i).ge.200.or.parini%znucl(j).ge.200)) then 
         imin=(k-1)*fp_size+1
         imax=k*fp_size
         num=num+dot_product(fp1(imin:imax),fp2(imin:imax))*w_ab(k)**2
         tmp_1=tmp_1+dot_product(fp1(imin:imax),fp1(imin:imax))*w_ab(k)**2
         tmp_2=tmp_2+dot_product(fp2(imin:imax),fp2(imin:imax))*w_ab(k)**2
       endif
  enddo
endif


!!do i_kind=1,nkinds
!!  do j_kind=1,i_kind
!!      num=num+dot_product(fp1(:,i_kind,j_kind),fp2(:,i_kind,j_kind))*w_ab(i_kind,j_kind)**2
!!      tmp_1=tmp_1+dot_product(fp1(:,i_kind,j_kind),fp1(:,i_kind,j_kind))*w_ab(i_kind,j_kind)**2
!!      tmp_2=tmp_2+dot_product(fp2(:,i_kind,j_kind),fp2(:,i_kind,j_kind))*w_ab(i_kind,j_kind)**2
!!!     do i_fp=1,fp_size
!!!     num=num+fp1(i_fp,i_kind,j_kind)*fp2(i_fp,i_kind,j_kind)*w_ab(i_kind,j_kind) 
!!!     tmp_1=tmp_1+fp1(i_fp,i_kind,j_kind)**2*w_ab(i_kind,j_kind)
!!!     tmp_2=tmp_2+fp2(i_fp,i_kind,j_kind)**2*w_ab(i_kind,j_kind)
!!!     enddo
!!  enddo
!!enddo
denom=sqrt(tmp_1*tmp_2)
!write(*,*) num, denom
distance=0.5d0*(1.d0-num/denom)
end subroutine get_cosinedistance




subroutine readfile_single(nat,latvec,nkinds,filename)
use parameters, only: kinds,rxyz
implicit none
integer:: nat,iat,nkinds
real(8):: latvec(3,3),dproj(6)
character(40):: filename
if(nkinds.ne.1) stop "Wrong number of kinds!"
open(unit=22,file=trim(filename))
read(22,*) nat
allocate(rxyz(3,nat),kinds(nat))
read(22,*) dproj(1:3)
read(22,*) dproj(4:6)
do iat=1,nat
read(22,*) rxyz(:,iat)
enddo
close(22)
kinds=1
call dproj2latvec(dproj,latvec)
endsubroutine readfile_single


subroutine readfile_double(nat,latvec,nkinds,filename)
use parameters, only: kinds,rxyz
implicit none
integer:: nat,iat,nkinds
real(8):: latvec(3,3),dproj(6)
character(4):: atom
character(40):: filename
if(nkinds.ne.2) stop "Wrong number of kinds!"
open(unit=22,file=trim(filename))
read(22,*) nat
allocate(rxyz(3,nat),kinds(nat))
read(22,*) dproj(1:3)
read(22,*) dproj(4:6)
do iat=1,nat
  read(22,*) rxyz(1,iat), rxyz(2,iat), rxyz(3,iat), atom
  if(trim(atom)=="A") then
     Kinds(iat)=1
  elseif(trim(atom)=="B") then
     Kinds(iat)=2
  else
     stop "Wrong atom types"
  endif
enddo
close(22)
call dproj2latvec(dproj,latvec)
endsubroutine readfile_double

subroutine readfile(nat,latvec,nkinds,filename)
use parameters, only: kinds,rxyz
implicit none
integer:: nat,iat,nkinds
real(8):: latvec(3,3),dproj(6)
character(4):: atom
character(40):: filename
open(unit=22,file=trim(filename))
read(22,*) nat
allocate(rxyz(3,nat),kinds(nat))
read(22,*) dproj(1:3)
read(22,*) dproj(4:6)
if(nkinds.gt.4) stop "Not more than 4 atom types supported"
do iat=1,nat
  read(22,*) rxyz(1,iat), rxyz(2,iat), rxyz(3,iat), atom
  if(trim(atom)=="A") then
     Kinds(iat)=1
  elseif(trim(atom)=="B") then
     Kinds(iat)=2
  elseif(trim(atom)=="C") then
     Kinds(iat)=3
  elseif(trim(atom)=="D") then
     Kinds(iat)=4
  else
     stop "Wrong atom types"
  endif
enddo
close(22)
call dproj2latvec(dproj,latvec)
endsubroutine readfile

subroutine n_rep(latvec,cut,nec)
!This subroutine will return how many periodic expansions are necessary for the periodic boundary conditions
!with for the given cut. 
implicit none
real*8 :: latvec(3,3),cut,nvec(3,3),point(3),point0(3),dist(3),eps,dd
integer:: i
integer:: nec
! eps=1.d-6
nec=0
call nveclatvec(latvec,nvec)
point0=(/0.d0,0.d0,0.d0/)
do i=1,3
call dist2plane(latvec(:,mod(i+1,3)+1),nvec(:,i),point0,dist(i))
! write(*,*) "cut",i,cut, dist
enddo
dd=minval(abs(dist))
nec=int(cut/dd)+1
end subroutine


subroutine n_rep_dim(latvec,cut,nec1,nec2,nec3)
!This subroutine will return how many periodic expansions for each lattice vector direction are necessary for the periodic boundary conditions
!with for the given cut. nec1,nec2,nec3 for latvec(:,1),latvec(:,2),latvec(:,3)
implicit none
real*8 :: latvec(3,3),cut,nvec(3,3),point(3),point0(3),dist(3),eps,dd
integer:: i
integer:: nec1,nec2,nec3
! eps=1.d-6
nec1=0
nec2=0
nec3=0
call nveclatvec(latvec,nvec)
point0=(/0.d0,0.d0,0.d0/)
do i=1,3
call dist2plane(latvec(:,mod(i+1,3)+1),nvec(:,i),point0,dist(i))
! write(*,*) "cut",i,cut, dist
enddo
dist=abs(dist)
nec1=int(cut/dist(2))+1
nec2=int(cut/dist(3))+1
nec3=int(cut/dist(1))+1
end subroutine

subroutine expand_dim(rxyz,rxyzout,transvecall,latvec,nat,nec1,nec2,nec3)
!This subroutine will expand the unit cell into the necessary periodic cells and store them in rxyzout. Only the first octant is computed
implicit none
real*8, intent(in)  :: rxyz(3,nat),latvec(3,3)
integer, intent(in) :: nat,nec1,nec2,nec3
real*8, intent(out) :: rxyzout(3,nat,2*nec1+1,2*nec2+1,2*nec3+1) !only necessary periodic images in the first octant plus the main cell
integer             :: iat,m,k,l
real*8,intent(inout):: transvecall(3,2*nec1+1,2*nec2+1,2*nec3+1)
do m=-nec3,nec3
   do k=-nec2,nec2
      do l=-nec1,nec1
      transvecall(:,l+nec1+1,k+nec2+1,m+nec3+1)=real(l,8)*latvec(:,1)+real(k,8)*latvec(:,2)+real(m,8)*latvec(:,3)
      enddo
   enddo
enddo

do m=0,2*nec3
   do k=0,2*nec2
      do l=0,2*nec1
      do iat=1,nat
      rxyzout(:,iat,l+1,k+1,m+1)=rxyz(:,iat)+transvecall(:,l+1,k+1,m+1)
      enddo
      enddo
   enddo
enddo
end


!subroutine expand_dim(rxyz,rxyzout,transvecall,latvec,nat,nec1,nec2,nec3)
!!This subroutine will expand the unit cell into the necessary periodic cells and store them in rxyzout. Only the first octant is computed
!implicit none
!real*8, intent(in)  :: rxyz(3,nat),latvec(3,3)
!integer, intent(in) :: nat,nec1,nec2,nec3
!real*8, intent(out) :: rxyzout(3,nat,nec1+1,nec2+1,nec3+1) !only necessary periodic images in the first octant plus the main cell
!integer             :: iat,m,k,l
!real*8,intent(inout):: transvecall(3,nec1+1,nec2+1,nec3+1)
!do m=0,nec3
!   do k=0,nec2
!      do l=0,nec1
!      transvecall(:,l+1,k+1,m+1)=real(l,8)*latvec(:,1)+real(k,8)*latvec(:,2)+real(m,8)*latvec(:,3)
!      enddo
!   enddo
!enddo
!
!do m=0,nec3
!   do k=0,nec2
!      do l=0,nec1
!      do iat=1,nat
!      rxyzout(:,iat,l+1,k+1,m+1)=rxyz(:,iat)+transvecall(:,l+1,k+1,m+1)
!      enddo
!      enddo
!   enddo
!enddo
!end


subroutine ggt_dif(a0,b0,c)
   INTEGER:: a0,b0,c,a,b
   a=a0
   b=b0
10 IF (A .EQ. B) GO TO 20
   IF (A .GT. B) THEN
     A = A - B
   ELSE
     B = B - A
   END IF
   GO TO 10
20 c=a
end subroutine

subroutine create_ggtlist(n1,n2,n3,ggt_0,list,nlist)
!This subroutine will create a ggt_0 grid with 0, if this node is a linear combination of a node vectors with n1'<n1 and n2'<n2 and n3'<n3
!else there will be a 1 in that grid point 
implicit none
integer:: n1,n2,n3, ggt_0(n1,n2,n3),k,l,m,n_tmp,n
integer:: nlist,list(3,n1*n2*n3)
ggt_0=0
ggt_0(1:2,1:2,1:2)=1
nlist=7
list(:,1)=(/1,0,0/)
list(:,2)=(/0,1,0/)
list(:,3)=(/0,0,1/)
list(:,4)=(/1,1,0/)
list(:,5)=(/1,0,1/)
list(:,6)=(/0,1,1/)
list(:,7)=(/1,1,1/)

do k=2,n1
  do l=2,n2
    do m=2,n3
    call ggt_dif(k-1,l-1,n_tmp)
    call ggt_dif(m-1,n_tmp,n)
!    call ggt_dif(k-1,l-1,n)
    if (n==1) then
    ggt_0(k,l,m) = 1 
    nlist=nlist+1
    list(1,nlist)=k
    list(2,nlist)=l
    list(3,nlist)=m
!    if (n==1) ggt_0(k,l,2) = 1  
    endif
    enddo
  enddo
enddo

!open(unit=60, file="ggtlist")
!write(60,*) n1,n2,n3
!do k=1,n1
!  do l=1,n2
!    do m=1,n3
!    write(60,*) ggt_0(k,l,m)
!    enddo
!  enddo
!enddo


end subroutine





!subroutine get_volume(latvec,vol)
!implicit none
!real(8):: latvec(3,3),v(3,3),vol
! v=latvec
! vol= v(1,1)*v(2,2)*v(3,3)-v(1,1)*v(2,3)*v(3,2)-v(1,2)*v(2,1)*v(3,3)+&
!      v(1,2)*v(2,3)*v(3,1)+v(1,3)*v(2,1)*v(3,2)-v(1,3)*v(2,2)*v(3,1)
!end subroutine

subroutine gaussian(x,sigma,norm,value)
implicit none
real(8)::x,sigma,norm,value
value=norm*exp(-0.5d0*x*x/(sigma*sigma))
end subroutine


subroutine integrate_gaussian_exact(floor_x0,ceil_x0,x,sigma,int_val)
implicit none
real(8):: floor_x,ceil_x,floor_x0,ceil_x0,x,mu,sigma,phi(600),int_val
real(8):: val_floor, val_ceil
mu=x
val_floor=0.5d0*(1.d0+erf((floor_x0-mu)/(sigma*sqrt(2.d0))))
val_ceil =0.5d0*(1.d0+erf((ceil_x0 -mu)/(sigma*sqrt(2.d0))))
int_val=val_ceil-val_floor
end subroutine


subroutine integrate_gaussian(floor_x0,ceil_x0,x,sigma,phi,int_val)
implicit none
real(8):: floor_x,ceil_x,floor_x0,ceil_x0,x,mu,sigma,phi(600),int_val
real(8):: val_floor, val_ceil
mu=x
if(floor_x0.gt.ceil_x0) stop "Wrong floor/ceiling values in integrate_gaussian"
floor_x=(floor_x0-mu)/sigma
ceil_x=(ceil_x0-mu)/sigma

floor_x=max(floor_x,-5.99d0)
floor_x=min(floor_x,5.99d0)
ceil_x=max(ceil_x,-5.99d0)
ceil_x=min(ceil_x,5.99d0)

if(floor_x.lt.0.d0) val_floor=1.d0-phi(int(-floor_x*100.d0+1.d0))
if(floor_x.ge.0.d0) val_floor=phi(int(floor_x*100.d0+1.d0))

if(ceil_x.lt.0.d0) val_ceil=1.d0-phi(int(-ceil_x*100.d0+1.d0))
if(ceil_x.ge.0.d0) val_ceil=phi(int(ceil_x*100.d0+1.d0))

int_val=val_ceil-val_floor
end subroutine


subroutine set_phi(phi)
implicit none
real(8):: phi(600)
!Phi contains the values of Phi(z), z=0,6,0.01 of the standard normal distribution, size of array is 600
phi=(/&
.5000000000000000000d0,.5039893563146315980d0,.5079783137169019414d0,.5119664734141126106d0,.5159534368528307935d0,&
.5199388058383724864d0,.5239221826541068383d0,.5279031701805211307d0,.5318813720139873302d0,.5358563925851721477d0,&
.5398278372770289879d0,.5437953125423168332d0,.5477584260205838884d0,.5517167866545611421d0,.5556700048059064478d0,&
.5596176923702425032d0,.5635594628914328830d0,.5674949316750383943d0,.5714237159009006861d0,.5753454347347954911d0,&
.5792597094391029877d0,.5831661634824423235d0,.5870644226482145678d0,.5909541151420059091d0,.5948348716977958084d0,&
.5987063256829237012d0,.6025681132017605135d0,.6064198731980394719d0,.6102612475557972482d0,.6140918811988773651d0,&
.6179114221889526748d0,.6217195218220192832d0,.6255158347233200633d0,.6293000189406535716d0,.6330717360360280654d0,&
.6368306511756191002d0,.6405764332179912923d0,.6443087548005467236d0,.6480272924241627930d0,.6517317265359823253d0,&
.6554217416103241822d0,.6590970262276774072d0,.6627572731517504812d0,.6664021794045423830d0,.6700314463394063669d0,&
.6736447797120799219d0,.6772418897496522705d0,.6808224912174442034d0,.6843863034837773807d0,.6879330505826094511d0,&
.6914624612740131182d0,.6949742691024806129d0,.6984682124530338099d0,.7019440346051235569d0,.7054014837843018970d0,&
.7088403132116536387d0,.7122602811509729515d0,.7156611509536758842d0,.7190426911014355937d0,.7224046752465350663d0,&
.7257468822499264505d0,.7290690962169943390d0,.7323711065310170021d0,.7356527078843224654d0,.7389137003071384324d0,&
.7421538891941352745d0,.7453730853286638647d0,.7485711049046899213d0,.7517477695464294118d0,.7549029063256904593d0,&
.7580363477769269664d0,.7611479319100131757d0,.7642375022207488211d0,.7673049076991025341d0,.7703500028352093798d0,&
.7733726476231317370d0,.7763727075624005103d0,.7793500536573503279d0,.7823045624142668242d0,.7852361158363628801d0,&
.7881446014166032521d0,.7910299121283983492d0,.7938919464141869220d0,.7967306081719316424d0,.7995458067395502244d0,&
.8023374568773076199d0,.8051054787481916053d0,.8078497978963037340d0,.8105703452232878581d0,.8132670569628273061d0,&
.8159398746532404711d0,.8185887451082027866d0,.8212136203856281735d0,.8238144577547420466d0,.8263912196613754091d0,&
.8289438736915182293d0,.8314723925331621857d0,.8339767539364704163d0,.8364569406723076916d0,.8389129404891690900d0,&
.8413447460685428148d0,.8437523549787453447d0,.8461357696272651108d0,.8484949972116562211d0,.8508300496690185399d0,&
.8531409436241039757d0,.8554277003360903908d0,.8576903456440607698d0,.8599289099112308321d0,.8621434279679645041d0,&
.8643339390536173283d0,.8665004867572527747d0,.8686431189572692002d0,.8707618877599820895d0,.8728568494372017650d0,&
.8749280643628496446d0,.8769755969486565661d0,.8789995155789818160d0,.8809998925447992679d0,.8829768039768912669d0,&
.8849303297782917799d0,.8868605535560226683d0,.8887675625521653799d0,.8906514475743080306d0,.8925123029254130591d0,&
.8943502263331446489d0,.8961653188786995461d0,.8979576849251809101d0,.8997274320455579399d0,.9014746709502521327d0,&
.9031995154143896976d0,.9049020822047608714d0,.9065824910065282127d0,.9082408643497191791d0,.9098773275355475088d0,&
.9114920085625979329d0,.9130850380529149657d0,.9146565491780329626d0,.9162066775849857514d0,.9177355613223310282d0,&
.9192433407662290445d0,.9207301585466076688d0,.9221961594734535694d0,.9236414904632609391d0,.9250663004656729527d0,&
.9264707403903515992d0,.9278549630341061949d0,.9292191230083144404d0,.9305633766666683293d0,.9318878820332745505d0,&
.9331927987311419148d0,.9344782879110835605d0,.9357445121810641364d0,.9369916355360214943d0,.9382198232881880928d0,&
.9394292419979409781d0,.9406200594052069874d0,.9417924443614469343d0,.9429465667622458636d0,.9440825974805305831d0,&
.9452007083004421162d0,.9463010718518802822d0,.9473838615457479406d0,.9484492515099106624d0,.9494974165258962540d0,&
.9505285319663518973d0,.9515427737332772251d0,.9525403181970526489d0,.9535213421362800368d0,.9544860226784501744d0,&
.9554345372414569937d0,.9563670634759681155d0,.9572837792086710262d0,.9581848623864051007d0,.9590704910211926837d0,&
.9599408431361828864d0,.9607960967125173113d0,.9616364296371288090d0,.9624620196514832582d0,.9632730443012738064d0,&
.9640696808870742318d0,.9648521064159611971d0,.9656204975541100577d0,.9663750305803716634d0,.9671158813408361477d0,&
.9678432252043862594d0,.9685572370192473413d0,.9692580910705340669d0,.9699459610388002639d0,.9706210199595906030d0,&
.9712834401839982590d0,.9719333933402274361d0,.9725710502961631976d0,.9731965811229450480d0,.9738101550595472666d0,&
.9744119404783613270d0,.9750021048517796274d0,.9755808147197775337d0,.9761482356584915143d0,.9767045322497882598d0,&
.9772498680518207914d0,.9777844055705685600d0,.9783083062323532086d0,.9788217303573276684d0,.9793248371339300373d0,&
.9798177845942955821d0,.9803007295906231988d0,.9807738277724826759d0,.9812372335650623167d0,.9816911001483410448d0,&
.9821355794371834369d0,.9825708220623429190d0,.9829969773523672405d0,.9834141933163950133d0,.9838226166278338791d0,&
.9842223926089095354d0,.9846136652160746294d0,.9849965770262678610d0,.9853712692240107485d0,.9857378815893311774d0,&
.9860965524865014098d0,.9864474188535800048d0,.9867906161927437747d0,.9871262785613980073d0,.9874545385640534079d0,&
.9877755273449553286d0,.9880893745814529616d0,.9883962084780963941d0,.9886961557614472040d0,.9889893416755886069d0,&
.9892758899783242743d0,.9895559229380488375d0,.9898295613312804164d0,.9900969244408357461d0,.9903581300546416832d0,&
.9906132944651614425d0,.9908625324694273484d0,.9911059573696632263d0,.9913436809744834433d0,.9915758136006542767d0,&
.9918024640754039556d0,.9920237397392662748d0,.9922397464494463470d0,.9924505885836907293d0,.9926563690446517096d0,&
.9928571892647286568d0,.9930531492113756631d0,.9932443473928593836d0,.9934308808644531918d0,.9936128452350567741d0,&
.9937903346742238408d0,.9939634419195872983d0,.9941322582846674472d0,.9942968736670493302d0,.9944573765569173496d0,&
.9946138540459332766d0,.9947663918364442193d0,.9949150742510088907d0,.9950599842422294117d0,.9952012034028738796d0,&
.9953388119762812680d0,.9954728888670326681d0,.9956035116518786587d0,.9957307565909105929d0,.9958546986389640310d0,&
.9959754114572416661d0,.9960929674251473021d0,.9962074376523144537d0,.9963188919908250174d0,.9964273990476002485d0,&
.9965330261969593817d0,.9966358395933307968d0,.9967359041841086231d0,.9968332837226422383d0,.9969280407813494449d0,&
.9970202367649454445d0,.9971099319237738401d0,.9971971853672350061d0,.9972820550772987236d0,.9973645979220950863d0,&
.9974448696695721317d0,.9975229250012140891d0,.9975988175258108104d0,.9976725997932684997d0,.9977443233084577479d0,&
.9978140385450867678d0,.9978817949595953918d0,.9979476410050602819d0,.9980116241451056913d0,.9980737908678121162d0,&
.9981341866996159551d0,.9981928562191935139d0,.9982498430713239168d0,.9983051899807227070d0,.9983589387658430292d0,&
.9984111303526350678d0,.9984618047882619640d0,.9985110012547625535d0,.9985587580826600362d0,.9986051127645076964d0,&
.9986501019683698965d0,.9986937615512305744d0,.9987361265723276871d0,.9987772313064077201d0,.9988171092568955967d0,&
.9988557931689773239d0,.9988933150425907126d0,.9989297061453210613d0,.9989649970251971434d0,.9989992175233859406d0,&
.9990323967867815735d0,.9990645632804859844d0,.9990957448001775987d0,.9991259684843684097d0,.9991552608265413804d0,&
.9991836476871713835d0,.9992111543056243494d0,.9992378053119328474d0,.9992636247384460990d0,.9992886360313546490d0,&
.9993128620620841396d0,.9993363251385600776d0,.9993590470163399297d0,.9993810489096131011d0,.9994023515020655779d0,&
.9994229749576092336d0,.9994429389309753553d0,.9994622625781702796d0,.9994809645667930287d0,.9994990630862142789d0,&
.9995165758576162185d0,.9995335201438924067d0,.9995499127594078548d0,.9995657700796183320d0,.9995811080505496715d0,&
.9995959421981359672d0,.9996102876374179935d0,.9996241590815999611d0,.9996375708509669389d0,.9996505368816620551d0,&
.9996630707343231448d0,.9996751856025811733d0,.9996868943214187730d0,.9996982093753914445d0,.9997091429067093138d0,&
.9997197067231837764d0,.9997299123060365833d0,.9997397708175725928d0,.9997492931087195167d0,.9997584897264322201d0,&
.9997673709209644599d0,.9997759466530089512d0,.9997842266007053169d0,.9997922201665193631d0,.9997999364839926795d0,&
.9998073844243643427d0,.9998145726030667202d0,.9998215093860951530d0,.9998282028962540702d0,.9998346610192798689d0,&
.9998408914098424471d0,.9998469014974262770d0,.9998526984920925731d0,.9998582893901242219d0,.9998636809795542479d0,&
.9998688798455794835d0,.9998738923758614394d0,.9998787247657145993d0,.9998833830231845798d0,.9998878729740177107d0,&
.9998922002665225905d0,.9998963703763259492d0,.9999003886110240380d0,.9999042601147311027d0,.9999079898725258264d0,&
.9999115827147991853d0,.9999150433215020506d0,.9999183762262973119d0,.9999215858206164098d0,.9999246763576212782d0,&
.9999276519560749144d0,.9999305166041201343d0,.9999332741629702870d0,.9999359283705111512d0,.9999384828448167895d0,&
.9999409410875810256d0,.9999433064874657662d0,.9999455823233662777d0,.9999477717675981925d0,.9999498778890038020d0,&
.9999519036559824103d0,.9999538519394437497d0,.9999557255156878988d0,.9999575270692112605d0,.9999592591954413745d0,&
.9999609244034022293d0,.9999625251183088537d0,.9999640636840971819d0,.9999655423658849740d0,.9999669633523706747d0,&
.9999683287581668800d0,.9999696406260734083d0,.9999709009292880868d0,.9999721115735593635d0,.9999732743992805206d0,&
.9999743911835259347d0,.9999754636420336018d0,.9999764934311314857d0,.9999774821496114630d0,.9999784313405517544d0,&
.9999793424930873975d0,.9999802170441317584d0,.9999810563800495267d0,.9999818618382818602d0,.9999826347089264544d0,&
.9999833762362704270d0,.9999840876202809037d0,.9999847700180519716d0,.9999854245452091117d0,.9999860522772731075d0,&
.9999866542509840972d0,.9999872314655862127d0,.9999877848840748040d0,.9999883154344053615d0,.9999888240106677983d0,&
.9999893114742250955d0,.9999897786548159750d0,.9999902263516271539d0,.9999906553343298476d0,.9999910663440871872d0,&
.9999914600945289944d0,.9999918372726972482d0,.9999921985399619073d0,.9999925445329086449d0,.9999928758641984938d0,&
.9999931931234007365d0,.9999934968777990374d0,.9999937876731730402d0,.9999940660345543186d0,.9999943324669582356d0,&
.9999945874560922654d0,.9999948314690427775d0,.9999950649549373960d0,.9999952883455880404d0,.9999955020561114294d0,&
.9999957064855300448d0,.9999959020173534441d0,.9999960890201397001d0,.9999962678480394107d0,.9999964388413203897d0,&
.9999966023268752613d0,.9999967586187126223d0,.9999969080184309966d0,.9999970508156771354d0,.9999971872885882185d0,&
.9999973177042202899d0,.9999974423189605943d0,.9999975613789262585d0,.9999976751203500935d0,.9999977837699518535d0,&
.9999978875452975036d0,.9999979866551451657d0,.9999980812997799617d0,.9999981716713364222d0,.9999982579541096817d0,&
.9999983403248555724d0,.9999984189530810585d0,.9999984940013224577d0,.9999985656254155586d0,.9999986339747554132d0,&
.9999986991925460256d0,.9999987614160426030d0,.9999988207767834814d0,.9999988774008146120d0,.9999989314089055004d0,&
.9999989829167574840d0,.9999990320352039053d0,.9999990788704038458d0,.9999991235240270893d0,.9999991660934340887d0,&
.9999992066718480510d0,.9999992453485209154d0,.9999992822088930033d0,.9999993173347474507d0,.9999993508043572010d0,&
.9999993826926280027d0,.9999994130712355211d0,.9999994420087567892d0,.9999994695707969949d0,.9999994958201117168d0,&
.9999995208167233862d0,.9999995446180351966d0,.9999995672789381285d0,.9999995888519161991d0,.9999996093871457159d0,&
.9999996289325920884d0,.9999996475341017543d0,.9999996652354916638d0,.9999996820786338780d0,.9999996981035375043d0,&
.9999997133484281875d0,.9999997278498227171d0,.9999997416426023022d0,.9999997547600819603d0,.9999997672340770194d0,&
.9999997790949677334d0,.9999997903717610104d0,.9999998010921489211d0,.9999998112825658758d0,.9999998209682428030d0,&
.9999998301732593298d0,.9999998389205938532d0,.9999998472321717236d0,.9999998551289106530d0,.9999998626307655680d0,&
.9999998697567704653d0,.9999998765250788235d0,.9999998829530025723d0,.9999998890570498400d0,.9999998948529597032d0,&
.9999999003557368260d0,.9999999055796842118d0,.9999999105384346221d0,.9999999152449801088d0,.9999999197117013239d0,&
.9999999239503948312d0,.9999999279722991963d0,.9999999317881205219d0,.9999999354080567615d0,.9999999388418200352d0,&
.9999999420986596110d0,.9999999451873824441d0,.9999999481163739379d0,.9999999508936165959d0,.9999999535267092288d0,&
.9999999560228840512d0,.9999999583890240018d0,.9999999606316789524d0,.9999999627570805849d0,.9999999647711581563d0,&
.9999999666795514885d0,.9999999684876255124d0,.9999999702004822577d0,.9999999718229739543d0,.9999999733597144669d0,&
.9999999748150899537d0,.9999999761932708564d0,.9999999774982211154d0,.9999999787337081614d0,.9999999799033127967d0,&
.9999999810104375220d0,.9999999820583150845d0,.9999999830500168052d0,.9999999839884605723d0,.9999999848764176136d0,&
.9999999857165200456d0,.9999999865112674247d0,.9999999872630331854d0,.9999999879740707476d0,.9999999886465193999d0,&
.9999999892824097403d0,.9999999898836692269d0,.9999999904521270633d0,.9999999909895189720d0,.9999999914974917470d0,&
.9999999919776081381d0,.9999999924313502930d0,.9999999928601239763d0,.9999999932652628987d0,.9999999936480313822d0,&
.9999999940096285789d0,.9999999943511912459d0,.9999999946737974099d0,.9999999949784682540d0,.9999999952661724478d0,&
.9999999955378275907d0,.9999999957943032092d0,.9999999960364234219d0,.9999999962649687157d0,.9999999964806787212d0,&
.9999999966842540999d0,.9999999968763578773d0,.9999999970576185504d0,.9999999972286311989d0,.9999999973899589278d0,&
.9999999975421349774d0,.9999999976856641659d0,.9999999978210243334d0,.9999999979486675628d0,.9999999980690219559d0,&
.9999999981824920781d0,.9999999982894614003d0,.9999999983902918554d0,.9999999984853266133d0,.9999999985748896369d0,&
.9999999986592875700d0,.9999999987388104028d0,.9999999988137320273d0,.9999999988843119025d0,.9999999989507948328d0/)
end subroutine

