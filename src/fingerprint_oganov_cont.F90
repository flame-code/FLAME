subroutine get_fp_coganov_atomic(nat,rxyz,latvec,r_cut,sigma,rcov,kinds,nkinds,nkinds_sum,fp_size,fp_dim,fp_coganov_atomic)
implicit none
integer:: i,j,k,l,m,iat,jat,nat,nec,nec1,nec2,nec3,nexp,fp_size,i_bin,near_bin,n_peaks,fp_dim,iarr,nn
integer:: nkinds_sum(nkinds),kinds(nat),nkinds,imin,imax,np1,np2,npeaks(fp_dim)
real(8):: latvec(3,3),rxyz(3,nat),fp_norm(fp_dim),fp_coganov_atomic(3,fp_size,nkinds,nat),fp_integrated(3*fp_size*fp_dim),fp_coganov(3*fp_size*fp_dim)
real(8), allocatable:: rxyzexp(:,:,:,:,:),transvecall(:,:,:,:)
real(8):: sigma,r_cut,d_bin,vol,r,value,val_all,rij,val,val2
real(8), parameter:: pi=acos(-1.d0)
real(8):: drx,dry,drz,t(2),lx,ly,lz,rmax,drx2,dry2,drz2,lx2,ly2,lz2,rmax2,l2,nom_f,discr,rxyzj(3),transvec(3),prefac
integer, allocatable:: list(:,:)
integer:: ilist,nlist,t_it 
integer, allocatable:: ggt_0(:,:,:)
character(3):: fn

real(8):: rcov(nkinds),sigma_cov


!Calculate the volume of the cell
call getvol(latvec,vol)

!Create the expanded unit cells
!call n_rep(latvec,r_cut,nec)
!write(*,'(a,i1.1,a)') " Creating expansion ",nec,"..."
call n_rep_dim(latvec,r_cut+6.d0*sigma,nec1,nec2,nec3)
write(*,'(a,3(i3,1x),a)') " Creating expansion for periodic Oganov FP ",nec1,nec2,nec3,"..."
!allocate(rxyzexp(3,nat,nec1+1,nec2+1,nec3+1),transvecall(3,nec1+1,nec2+1,nec3+1))

allocate(rxyzexp(3,nat,2*nec1+1,2*nec2+1,2*nec3+1),transvecall(3,2*nec1+1,2*nec2+1,2*nec3+1))
!write(*,*) "Allocated"

call expand_dim(rxyz,rxyzexp,transvecall,latvec,nat,nec1,nec2,nec3)

rmax=r_cut
rmax2=r_cut*r_cut
n_peaks=0.d0
fp_coganov_atomic=0.d0
fp_integrated=0.d0
npeaks=0
!Doing per atom fingerprint
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
        if(rij.lt.r_cut) then
!       write(*,*) rij
        imax=max(Kinds(iat),Kinds(jat))
        imin=min(Kinds(iat),Kinds(jat))
        np1=fp_coganov_atomic(1,fp_size,Kinds(jat),iat)
        np2=fp_coganov_atomic(1,fp_size,Kinds(iat),jat)
        np1=np1+1 
        np2=np2+1 
        if(np1.ge.fp_size) stop "Fingerprint size exceeded np1"
        if(np2.ge.fp_size) stop "Fingerprint size exceeded np2"
        fp_coganov_atomic(1,fp_size,Kinds(jat),iat)=np1
        fp_coganov_atomic(1,fp_size,Kinds(iat),jat)=np2
        sigma_cov=0.01d0*(rcov(Kinds(iat))+rcov(Kinds(jat)))
        fp_coganov_atomic(1,np1,Kinds(jat),iat)=vol/(4.d0*pi*rij**2*nkinds_sum(Kinds(jat)))!*real(1/(abs(imin-imax)+1),8)
        fp_coganov_atomic(2,np1,Kinds(jat),iat)=rij
        fp_coganov_atomic(3,np1,Kinds(jat),iat)=sigma!_cov
        fp_coganov_atomic(1,np2,Kinds(iat),jat)=vol/(4.d0*pi*rij**2*nkinds_sum(Kinds(iat)))!*real(1/(abs(imin-imax)+1),8)
        fp_coganov_atomic(2,np2,Kinds(iat),jat)=rij
        fp_coganov_atomic(3,np2,Kinds(iat),jat)=sigma!_cov

        iarr=imax*(imax+1)
        iarr=(iarr/2+imin-imax)
        npeaks(iarr)=npeaks(iarr)+1
        nn=npeaks(iarr)
        iarr=(iarr-1)*3*fp_size+(nn-1)*3
        prefac=1.d0
        if(iat.ne.jat.and.imax==imin) prefac=2.d0
        fp_integrated(iarr+1)=vol/(4.d0*pi*rij**2*nkinds_sum(Kinds(jat))*nkinds_sum(Kinds(iat)))*prefac
        fp_integrated(iarr+2)=rij
        fp_integrated(iarr+3)=sigma!_cov
!!        if(iat.ne.jat) then
!!          iarr=imax*(imax+1)
!!          iarr=(iarr/2+imin-imax)
!!          npeaks(iarr)=npeaks(iarr)+1
!!          nn=npeaks(iarr)
!!          iarr=(iarr-1)*3*fp_size+(nn-1)*3
!!          write(*,*) real(1/(abs(imin-imax)+1),8)
!!          fp_integrated(iarr+1)=vol/(4.d0*pi*rij**2*nkinds_sum(Kinds(jat))*nkinds_sum(Kinds(iat)))*real(1/(abs(imin-imax)+1),8)
!!          fp_integrated(iarr+2)=rij
!!          fp_integrated(iarr+3)=sigma!_cov
!!        endif
        endif
      2020 continue
     enddo
     enddo
     enddo


  enddo
enddo


!!!!!!!!Calculate the Normalization constant
!!!!!!!fp_norm=1.d0
!!!!!!!!!do k=1,fp_dim
!!!!!!!!!!Stephan's clever formula...
!!!!!!!!!   i=ceiling(-0.5d0+0.5d0*sqrt(1.d0+8.d0*real(k,8)))
!!!!!!!!!   j=(i*(1-i))/2+k
!!!!!!!!!!write(*,*) i,j
!!!!!!!!!   fp_norm(k)=nkinds_sum(i)*nkinds_sum(j) 
!!!!!!!!!!   fp_integrated((k-1)*3*fp_size+1:k*3*fp_size)=fp_integrated((k-1)*3*fp_size+1:k*3*fp_size)/fp_norm(k)
!!!!!!!!!!   fp((k-1)*fp_size+1:k*fp_size)=fp((k-1)*fp_size+1:k*fp_size)/fp_norm(k)-1.d0
!!!!!!!!!!write(*,*) fp_norm(k)
!!!!!!!!!enddo
!!!!!!!
!!!!!!!!Plot the FP
!!!!!!!write(*,*) fp_size,fp_dim
!!!!!!!write(*,*) size(fp_integrated)
!!!!!!!call fp_atomic2coganov(fp_coganov_atomic,fp_coganov,nkinds_sum,nat,nkinds,kinds,fp_size,fp_dim)
!!!!!!!do k=1,fp_dim
!!!!!!!do i=1,fp_size
!!!!!!!write(22,*) fp_integrated((k-1)*3*fp_size+(i-1)*3+2),fp_integrated((k-1)*3*fp_size+(i-1)*3+3)
!!!!!!!enddo
!!!!!!!enddo
!!!!!!!do k=1,fp_dim
!!!!!!! write(fn,'(i3.3)') k
!!!!!!! open(unit=54,file="ContFP"//fn)
!!!!!!! open(unit=55,file="ContFP_sum"//fn)
!!!!!!! do r=0,r_cut,0.01d0
!!!!!!!   val=0.d0
!!!!!!!   val2=0.d0
!!!!!!!   do i=1,fp_size
!!!!!!!     val=val+1.d0/sigma/sqrt(2.d0*pi)*exp(-0.5d0*((r-fp_integrated((k-1)*3*fp_size+(i-1)*3+2))/fp_integrated((k-1)*3*fp_size+(i-1)*3+3))**2)/fp_norm(k)*fp_integrated((k-1)*3*fp_size+(i-1)*3+1)
!!!!!!!     val2=val2+1.d0/sigma/sqrt(2.d0*pi)*exp(-0.5d0*((r-fp_coganov((k-1)*3*fp_size+(i-1)*3+2))/fp_coganov((k-1)*3*fp_size+(i-1)*3+3))**2)/fp_norm(k)*fp_coganov((k-1)*3*fp_size+(i-1)*3+1)
!!!!!!!   enddo  
!!!!!!! write(54,*) r,val-1.d0
!!!!!!! write(55,*) r,val2-1.d0
!!!!!!! enddo
!!!!!!! close(54)
!!!!!!! close(55)
!!!!!!!enddo
!!!!!!!stop

deallocate(rxyzexp,transvecall)

end subroutine

subroutine get_fp_coganov(nat,rxyz,latvec,r_cut,sigma,rcov,kinds,nkinds,nkinds_sum,fp_size,fp_dim,fp_coganov)
implicit none
integer:: i,j,k,l,m,iat,jat,nat,nec,nec1,nec2,nec3,nexp,fp_size,i_bin,near_bin,n_peaks,fp_dim,iarr,nn
integer:: nkinds_sum(nkinds),kinds(nat),nkinds,imin,imax,np1,np2,npeaks(fp_dim)
real(8):: latvec(3,3),rxyz(3,nat),fp_norm(fp_dim),fp_coganov(3*fp_size*fp_dim)
real(8), allocatable:: rxyzexp(:,:,:,:,:),transvecall(:,:,:,:)
real(8):: sigma,r_cut,d_bin,vol,r,value,val_all,rij,val,val2
real(8), parameter:: pi=acos(-1.d0)
real(8):: drx,dry,drz,t(2),lx,ly,lz,rmax,drx2,dry2,drz2,lx2,ly2,lz2,rmax2,l2,nom_f,discr,rxyzj(3),transvec(3),prefac
integer, allocatable:: list(:,:)
integer:: ilist,nlist,t_it 
integer, allocatable:: ggt_0(:,:,:)
character(3):: fn

real(8):: rcov(nkinds),sigma_cov


!Calculate the volume of the cell
call getvol(latvec,vol)
call estimate_nmax_per_atom(vol,nat,nkinds,r_cut,pi,i)

!Create the expanded unit cells
!call n_rep(latvec,r_cut,nec)
!write(*,'(a,i1.1,a)') " Creating expansion ",nec,"..."
call n_rep_dim(latvec,r_cut+6.d0*sigma,nec1,nec2,nec3)
write(*,'(a,3(i3,1x),a)') " Creating expansion for periodic Oganov FP ",nec1,nec2,nec3,"..."
!allocate(rxyzexp(3,nat,nec1+1,nec2+1,nec3+1),transvecall(3,nec1+1,nec2+1,nec3+1))

allocate(rxyzexp(3,nat,2*nec1+1,2*nec2+1,2*nec3+1),transvecall(3,2*nec1+1,2*nec2+1,2*nec3+1))
!write(*,*) "Allocated"

call expand_dim(rxyz,rxyzexp,transvecall,latvec,nat,nec1,nec2,nec3)

rmax=r_cut
rmax2=r_cut*r_cut
n_peaks=0.d0
fp_coganov=0.d0
npeaks=0
!Doing per atom fingerprint
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
        if(rij.lt.r_cut) then
!       write(*,*) rij
        imax=max(Kinds(iat),Kinds(jat))
        imin=min(Kinds(iat),Kinds(jat))
        sigma_cov=0.01d0*(rcov(Kinds(iat))+rcov(Kinds(jat)))
        iarr=imax*(imax+1)
        iarr=(iarr/2+imin-imax)
        npeaks(iarr)=npeaks(iarr)+1
        nn=npeaks(iarr)
        fp_coganov(iarr*3*fp_size-2)=nn
        if(nn.ge.fp_size) stop "Increase the fingerprint size!"
        iarr=(iarr-1)*3*fp_size+(nn-1)*3
        prefac=1.d0
        if(iat.ne.jat.and.imax==imin) prefac=2.d0
        fp_coganov(iarr+1)=vol/(4.d0*pi*rij**2*nkinds_sum(Kinds(jat))*nkinds_sum(Kinds(iat)))*prefac
        fp_coganov(iarr+2)=rij
        fp_coganov(iarr+3)=sigma!_cov
        endif
      2020 continue
     enddo
     enddo
     enddo


  enddo
enddo

deallocate(rxyzexp,transvecall)

end subroutine

subroutine fp_atomic2coganov(fp_atomic,fp_coganov,nkinds_sum,nat,nkinds,kinds,fp_size,fp_dim)
implicit none
integer:: fp_size,fp_dim,nat,nkinds,kinds(nat),nkinds_sum(nkinds)
integer:: iat,ikind,n_gauss,igauss,npeaks(fp_dim),imax,imin,nn,iarr
real(8):: fp_atomic(3,fp_size,nkinds,nat),fp_coganov(3*fp_size*fp_dim),prefac
npeaks=0
do iat=1,nat
   do ikind=1,nkinds
      n_gauss=fp_atomic(1,fp_size,ikind,iat)
      write(*,*) n_gauss
      do igauss=1,n_gauss
        prefac=0.5d0
        imax=max(Kinds(iat),ikind)
        imin=min(Kinds(iat),ikind)
        if(imax==imin) prefac=1.d0
        iarr=imax*(imax+1)
        iarr=(iarr/2+imin-imax)
        npeaks(iarr)=npeaks(iarr)+1
        nn=npeaks(iarr)
        if(nn.gt.fp_size-1) then 
           write(*,*) nn,fp_size,iarr
           stop "Increase the fp size!"
        endif
        fp_coganov(iarr*3*fp_size-2)=nn
        iarr=(iarr-1)*3*fp_size+(nn-1)*3
        fp_coganov(iarr+1)=fp_atomic(1,igauss,ikind,iat)/nkinds_sum(kinds(iat))*prefac
        fp_coganov(iarr+2)=fp_atomic(2,igauss,ikind,iat)
        fp_coganov(iarr+3)=fp_atomic(3,igauss,ikind,iat)
      enddo
   enddo
enddo
end subroutine


subroutine get_cosinedistance_coganov_atomic(fp1,fp2,nat,fp_size,fp_dim,kinds,nkinds,nkinds_sum,r_cut,pi,distance)
!This sunroutine will compute the cosine distance between two fingerprints fp1 and fp2 and 
!store the output distance.  All conventions and 
!methods are from J.Chem.Phys, 130, 104504 (2009) and IEEE Symposium, Okt 21-23. (2008) (M.Valle and A.Oganov)
implicit none
integer, intent(IN) :: nkinds, nkinds_sum(nkinds), fp_size,fp_dim,kinds(nat),nat
real(8), intent(INOUT) :: fp1(3,fp_size,nkinds,nat),fp2(3,fp_size,nkinds,nat)
real(8), intent(OUT):: distance
integer             :: i_kind,j_kind,i_fp,i,j,k,imin,imax,n_gauss1,n_gauss2,n1t,n1a,n1e,n2t,n2a,n2e
real(8)             :: w_ab(nkinds,nkinds),w_norm,num,denom,tmp_1,tmp_2,denom_tmp1,denom_tmp2,num_tmp(4),r_cut,pi,cost
real(8),allocatable :: A(:,:)
INTEGER:: iat,jat, f(nat)


!Compute the weight of each fingerprint
w_ab=0.d0
w_norm=0.d0
do i=1,nkinds
  do j=1,nkinds
   w_ab(i,j)=real(nkinds_sum(i),8)*real(nkinds_sum(j),8)
   if(i.lt.j) w_norm=w_norm+real(nkinds_sum(i),8)*real(nkinds_sum(j),8)
  enddo
enddo
w_ab=w_ab/w_norm

      allocate(A(nat,nat))
      do iat=1,nat
      do jat=1,nat
         if(kinds(iat) .ne. kinds(jat))  then
             A(iat,jat)=1.d7   ! do not permute!
         else
         tmp_1=0.d0 
         tmp_2=0.d0 
         num=0.d0
             do k=1,nkinds
               n_gauss1=fp1(1,fp_size,k,iat)
               n_gauss2=fp2(1,fp_size,k,jat)
               num_tmp=0.d0
               denom_tmp1=0.d0
               denom_tmp2=0.d0
                   call gauss_int_at(fp1(:,1:n_gauss1,k,iat),fp2(:,1:n_gauss2,k,jat),n_gauss1,n_gauss2,pi,num_tmp(1))
                   num_tmp(2)=-1.d0*sum(fp1(1,1:n_gauss1,k,iat))
                   num_tmp(3)=-1.d0*sum(fp2(1,1:n_gauss2,k,jat))
                   num_tmp(4)=r_cut
                   num=num+sum(num_tmp(:))*w_ab(kinds(iat),k)**2
             !      call gauss_int_at(fp1(:,1:n_gauss1,k),fp1(:,1:n_gauss1,k),n_gauss1,n_gauss1,pi,denom_tmp1)
             !      call gauss_int_at(fp2(:,1:n_gauss2,k),fp2(:,1:n_gauss2,k),n_gauss2,n_gauss2,pi,denom_tmp2)
                   call gauss_int_at_self(fp1(:,1:n_gauss1,k,iat),n_gauss1,pi,denom_tmp1)
                   call gauss_int_at_self(fp2(:,1:n_gauss2,k,jat),n_gauss2,pi,denom_tmp2)
                   tmp_1=tmp_1+(denom_tmp1+2.d0*num_tmp(2)+num_tmp(4))*w_ab(kinds(iat),k)**2
                   tmp_2=tmp_2+(denom_tmp2+2.d0*num_tmp(3)+num_tmp(4))*w_ab(kinds(jat),k)**2
             enddo
             denom=sqrt(tmp_1*tmp_2)
             distance=abs(0.5d0*(1.d0-num/denom))
             A(iat,jat)=distance
        endif
      enddo
      enddo
      CALL APC(Nat,A,F,cost)
distance=cost/real(nat,8)
end subroutine



subroutine get_cosinedistance_coganov(fp1,fp2,fp_size,fp_dim,nkinds,nkinds_sum,r_cut,pi,distance)
!This sunroutine will compute the cosine distance between two fingerprints fp1 and fp2 and 
!store the output distance.  All conventions and 
!methods are from J.Chem.Phys, 130, 104504 (2009) and IEEE Symposium, Okt 21-23. (2008) (M.Valle and A.Oganov)
implicit none
integer, intent(IN) :: nkinds, nkinds_sum(nkinds), fp_size,fp_dim
real(8), intent(INOUT) :: fp1(3,fp_size,fp_dim),fp2(3,fp_size,fp_dim)
real(8), intent(OUT):: distance
integer             :: i_kind,j_kind,i_fp,i,j,k,imin,imax,n_gauss1,n_gauss2,n1t,n1a,n1e,n2t,n2a,n2e
real(8)             :: w_ab(fp_dim),w_norm,num,denom,tmp_1,tmp_2,denom_tmp1,denom_tmp2,num_tmp(4),r_cut,pi

!Compute the weight of each fingerprint
w_ab=0.d0
w_norm=0.d0
do k=1,fp_dim
   i=ceiling(-0.5d0+0.5d0*sqrt(1.d0+8.d0*real(k,8)))
   j=(i*(1-i))/2+k
   w_ab(k)=real(nkinds_sum(i),8)*real(nkinds_sum(j),8)
   w_norm=w_norm+real(nkinds_sum(i),8)*real(nkinds_sum(j),8)
enddo
w_ab=w_ab/w_norm

!Compute the distance
num=0.d0;  denom=0.d0;  tmp_1=0.d0;  tmp_2=0.d0

!if(trim(code).ne."lenosky_tb_lj") then
!Usual case: atoms in general

do k=1,fp_dim
  n_gauss1=fp1(1,fp_size,k)
  n_gauss2=fp2(1,fp_size,k)
  num_tmp=0.d0
  denom_tmp1=0.d0
  denom_tmp2=0.d0
  n1t=fp_size*fp_dim
  n2t=n1t
  n1a=(k-1)*fp_size
  n2a=n1a
  n1e=n1a+n_gauss1
  n1a=n1a+1
  n2e=n2a+n_gauss2
  n2a=n2a+1
      call gauss_int_at2(fp1,fp2,n1t,n1a,n1e,n2t,n2a,n2e,pi,num_tmp(1))
!      call gauss_int_at(fp1(:,1:n_gauss1,k),fp2(:,1:n_gauss2,k),n_gauss1,n_gauss2,pi,num_tmp(1))
      num_tmp(2)=-1.d0*sum(fp1(1,1:n_gauss1,k))
      num_tmp(3)=-1.d0*sum(fp2(1,1:n_gauss2,k))
      num_tmp(4)=r_cut
      num=num+sum(num_tmp(:))*w_ab(k)**2
!      call gauss_int_at(fp1(:,1:n_gauss1,k),fp1(:,1:n_gauss1,k),n_gauss1,n_gauss1,pi,denom_tmp1)
!      call gauss_int_at(fp2(:,1:n_gauss2,k),fp2(:,1:n_gauss2,k),n_gauss2,n_gauss2,pi,denom_tmp2)
      call gauss_int_at_self(fp1(:,1:n_gauss1,k),n_gauss1,pi,denom_tmp1)
      call gauss_int_at_self(fp2(:,1:n_gauss2,k),n_gauss2,pi,denom_tmp2)
      tmp_1=tmp_1+(denom_tmp1+2.d0*num_tmp(2)+num_tmp(4))*w_ab(k)**2
      tmp_2=tmp_2+(denom_tmp2+2.d0*num_tmp(3)+num_tmp(4))*w_ab(k)**2
enddo

denom=sqrt(tmp_1*tmp_2)
distance=abs(0.5d0*(1.d0-num/denom))
end subroutine get_cosinedistance_coganov

subroutine fp_coganov_local_order(fp_atomic,local_order,nkinds_sum,nat,nkinds,kinds,fp_size,r_cut,vol,pi)
!We have a sum of gaussians (with some forefactors) sum(gi)
!The atomic fingerprint is then f_at=sum(gi)-1
!The local order is a sum of |f_at|^2, which is <f_at|f_at>=<sum(gi)-1|sum(gi)-1>=<sum(gi)|sum(gi)>-2<1|sum(gi)>+<1|1>
!The second term are 2xsum of the prefactors of the gaussians
!The last term cannot be integrated from -inf to inf, so we only integrate from 0 to r_cut
implicit none
integer:: fp_size,nat,nkinds,kinds(nat),nkinds_sum(nkinds)
integer:: iat,ikind,n_gauss,igauss,imax,imin,nn,iarr
real(8):: fp_atomic(3,fp_size,nkinds,nat),local_order(nat),pi,r_cut,locord_at(3),vol,locord_sum
local_order=0.d0
do iat=1,nat
   do ikind=1,nkinds
      n_gauss=fp_atomic(1,fp_size,ikind,iat)
      call gauss_int_at(fp_atomic(:,1:n_gauss,ikind,iat),fp_atomic(:,1:n_gauss,ikind,iat),n_gauss,n_gauss,pi,locord_at(1))
      locord_at(2)=-2.d0*sum(fp_atomic(1,1:n_gauss,ikind,iat))
      locord_at(3)=r_cut
      locord_sum=locord_at(1)+locord_at(2)+locord_at(3)
      local_order(iat)=local_order(iat)+nkinds_sum(ikind)/(nat*(vol/nat)**(1.d0/3.d0))*locord_sum 
   enddo
   local_order(iat)=sqrt(local_order(iat))
enddo
end subroutine

subroutine gauss_int_penalty(fp1_all,fp2_all,nat,nmax,pi,penalty)
!This routine computes the penalty function of two fingerprint arrays
!fp_all has 3 first entries for the a, mu and sigma, then the list with nmax-1 neighbor interactions, then
!the dimenson of the number of atoms for each individual atom. The number of interaction per
!atom is given in fp_all(1,nmax,iat)
implicit none
integer:: n1,n2,nat,nmax,iat,jat
real(8):: fp1_all(3,nmax,nat)
real(8):: fp2_all(3,nmax,nat)
real(8):: pi,penalty,int_at_1(nat),int_at_ij,int_at_2(nat)
penalty=0.d0
do iat=1,nat
   n1=int(fp1_all(1,nmax,iat))
   call gauss_int_at(fp1_all(:,1:n1,iat),fp1_all(:,1:n1,iat),n1,n1,pi,int_at_1(iat))
   n2=int(fp2_all(1,nmax,iat))
   call gauss_int_at(fp2_all(:,1:n2,iat),fp2_all(:,1:n2,iat),n2,n2,pi,int_at_2(iat))
   do jat=1,iat-1
      n2=int(fp2_all(1,nmax,jat))
      call gauss_int_at(fp1_all(:,1:n1,iat),fp2_all(:,1:n2,jat),n1,n2,pi,int_at_ij)
      penalty=penalty+int_at_1(iat)-2.d0*int_at_ij+int_at_2(jat)
   enddo
enddo 
end subroutine

subroutine gauss_int_at2(fp1_at,fp2_at,n1t,n1a,n1e,n2t,n2a,n2e,pi,int_at)
!This routine computes the integral of two arrays containing sums
!of gaussians, usually used for per-atoms array of gaussians
implicit none
integer:: n1t,n1a,n1e,n2t,n2a,n2e,nat,nmax,i,j
real(8):: fp1_at(3*n1t)
real(8):: fp2_at(3*n2t)
real(8):: pi,integral,int_at
int_at=0.d0
do i=3*(n1a-1),3*(n1e-1),3
  do j=3*(n2a-1),3*(n2e-1),3
  call gauss_int(fp1_at(i+1),fp1_at(i+2),fp1_at(i+3),fp2_at(j+1),fp2_at(j+2),fp2_at(j+3),pi,integral)
  int_at=int_at+integral
  enddo
enddo
end subroutine

subroutine gauss_int_at(fp1_at,fp2_at,n1,n2,pi,int_at)
!This routine computes the integral of two arrays containing sums
!of gaussians, usually used for per-atoms array of gaussians
implicit none
integer:: n1,n2,nat,nmax,i,j
real(8):: fp1_at(3,n1)
real(8):: fp2_at(3,n2)
real(8):: pi,integral,int_at
int_at=0.d0
do i=1,n1
  do j=1,n2
  call gauss_int(fp1_at(1,i),fp1_at(2,i),fp1_at(3,i),fp2_at(1,j),fp2_at(2,j),fp2_at(3,j),pi,integral)
  int_at=int_at+integral
  enddo
enddo
end subroutine

subroutine gauss_int(a1,mu1,sigma1,a2,mu2,sigma2,pi,integral)
!This subroutine will compute the integral from -inf to inf from the product
!of two normalized gaussian curves with prefactors a1 and a2
!G1(a1,mu1,sigma1)=a1*exp(-(1/2)*((r-mu1)/sigma1)^2)/(sqrt(2*pi)*sigma1)
implicit none
real(8):: a1,mu1,sigma1,a2,mu2,sigma2,pi,integral
real(8):: mm1,mm2,ssigma1,ssigma2,ssum,e,sdenom,denom
mm1 = mu1 ** 2
mm2 = mu2 ** 2
ssigma1 = sigma1 ** 2
ssigma2 = sigma2 ** 2
ssum = ssigma1+ssigma2
e = exp(-0.5d0*((-2.d0 * mu1 * mu2 + mm1 + mm2) / ssum) )
sdenom = ssum/(ssigma1*ssigma2)
denom = sqrt(2.d0*pi*sdenom)*sigma1*sigma2
integral = a1*a2*e/denom
end subroutine

subroutine gauss_int_at_self(fp1_at,n1,pi,int_at)
!This routine computes the integral of two arrays containing sums
!of gaussians, usually used for per-atoms array of gaussians
implicit none
integer:: n1,n2,nat,nmax,i,j
real(8):: fp1_at(3,n1)
real(8):: pi,integral,int_at
int_at=0.d0
do i=1,n1
  do j=1,i
  call gauss_int(fp1_at(1,i),fp1_at(2,i),fp1_at(3,i),fp1_at(1,j),fp1_at(2,j),fp1_at(3,j),pi,integral)
!  call gauss_int_self(fp1_at(1,i),fp1_at(2,i),fp1_at(3,i),pi,integral)
  int_at=int_at+integral*(2.d0-real(1/(abs(i-j)+1),8))
  enddo
enddo
end subroutine


subroutine gauss_int_self(a1,mu1,sigma1,pi,integral)
!This subroutine will compute the integral from -inf to inf from the product
!of two normalized gaussian curves with prefactors a1 and a2
!G1(a1,mu1,sigma1)=a1*exp(-(1/2)*((r-mu1)/sigma1)^2)/(sqrt(2*pi)*sigma1)
implicit none
real(8):: a1,mu1,sigma1,a2,mu2,sigma2,pi,integral
real(8):: mm1,mm2,ssigma1,ssigma2,ssum,e,sdenom,denom
denom = sqrt(pi)*sigma1
integral = 0.5d0*a1**2/denom
end subroutine



subroutine estimate_nmax_per_atom(vol,nat,ntypat,r_cut,pi,nmax)
implicit none
integer:: nat,nmax,ntypat
real(8):: vol,r_cut,r_vol,pi,den
den=real(nat,8)/vol
r_vol=4.d0/3.d0*pi*r_cut**3
r_vol=r_vol*1.1d0
nmax=int(den*r_vol)*nat/ntypat
end subroutine
