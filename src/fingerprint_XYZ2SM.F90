!!subroutine diffsm(nat,cell1,rxyz1,xat1,cell2,rxyz2,xat2,diff)
!!! calculates th difference between two structures based on xyz2sm fingerprint 
!!implicit none
!!real*8  cell1(3,3),cell2(3,3)
!!integer nat
!!real*8  rxyz1(3,nat), rxyz2(3,nat)
!!character*2  xat1(nat),xat2(nat)   ! chemical symbol (used for covalent radius)
!!real*8  diff
!!
!!
!!real*8  w1 , w2  ! relative  widths of the GTO's on each atom (e.g. use w1=1, w2=1.5)
!!integer m  !  m = maximum power
!!parameter (w1=1.d0, w2=1.5d0, m=3)
!!
!!integer, parameter :: lfp=3*m
!!real*8  fp1(lfp, nat), fp2(lfp,nat)  
!!real*8  covrad, rcov1(nat), rcov2(nat)
!!integer iat
!!
!!
!!
!!do iat=1,nat
!!    rcov1(iat)=covrad(xat1(iat))
!!    rcov2(iat)=covrad(xat2(iat))
!!enddo
!!
!!call  xyz2sm(nat, cell1,rxyz1,rcov1,w1,w2,m,fp1)
!!print '(<lfp>f11.5)',fp1
!!print*
!!call  xyz2sm(nat, cell2,rxyz2,rcov2,w1,w2,m,fp2)
!!print '(<lfp>f11.5)',fp2
!!
!!call hungarian(nat,rxyz1,rxyz2,lfp,fp1,fp2,diff,xat1,xat2)
!!
!!diff=sqrt(sum((fp2-fp1)**2)/dble(lfp*nat))
!!
!!
!!endsubroutine diffsm

!!subroutine diffsm(nat,cell1,rxyz1,xat1,cell2,rxyz2,xat2,diff)
!!! calculates th difference between two structures based on xyz2sm fingerprint 
!!implicit none
!!real*8  cell1(3,3),cell2(3,3)
!!integer nat
!!real*8  rxyz1(3,nat), rxyz2(3,nat)
!!character*2  xat1(nat),xat2(nat)   ! chemical symbol (used for covalent radius)
!!real*8  diff
!!
!!
!!real*8  w1 , w2  ! relative  widths of the GTO's on each atom (e.g. use w1=1, w2=1.5)
!!integer m  !  m = maximum power
!!parameter (w1=1.d0, w2=1.5d0, m=3)
!!
!!integer, parameter :: lfp=3*m
!!real*8  fp1(lfp, nat), fp2(lfp,nat)  , ddot
!!real*8  covrad, rcov1(nat), rcov2(nat)
!!integer iat
!!
!!
!!
!!do iat=1,nat
!!    rcov1(iat)=covrad(xat1(iat))
!!    rcov2(iat)=covrad(xat2(iat))
!!enddo
!!
!!call  xyz2sm(nat, cell1,rxyz1,rcov1,w1,w2,m,fp1)
!!print '(<lfp>f11.5)',fp1
!!print*
!!call  xyz2sm(nat, cell2,rxyz2,rcov2,w1,w2,m,fp2)
!!print '(<lfp>f11.5)',fp2
!!
!!call hungarian(nat,rxyz1,rxyz2,lfp,fp1,fp2,diff,xat1,xat2)
!!
!!endsubroutine diffsm
! =============================================
subroutine xyz2sm(nat, cell,rxyz,rcov,w1,w2,m,fngprt)
implicit none
! in/out
real*8  cell(3,3)
integer nat
real*8  rxyz(3,nat), rcov(nat)
real*8  w1 , w2  ! relative  widths of the GTO's on each atom (e.g. use w1=1, w2=1.5)
integer m  !  m = maximum power
real*8  fngprt(3*m, nat)  ! Fingerprint array:  2 s-GTOs per atom ==> a symmetric 2x2 matrix per atom which has 3 independent elements 

! Local variables:
integer norb 
integer isph, jsph, nsph  ! # counter, and # of atoms in the sphere
integer, parameter :: nsphx=10000   ! max of nsphr; later on should be estimated on-the-fly
real*8 rtmp(3), rsph(3,nsphx)  ! rxyz in the sphere
integer  isph2iat(nsphx) , iat2sph(nat)
real*8,allocatable ::  sijsph(:,:), sijat(:,:), spow(:,:), stmp(:,:)
real*8,allocatable :: eval(:)

integer iat,jat,is,js,ip,jp, kxyz, kat, iorb, jorb, lat
integer iatcell
logical  periodic
real*8  ai,aj, xi,yi,zi, xj,yj,zj, xij, yij, zij, r2 ,sij, dtau 
!real*8  ovrlp(norb,norb), evec(norb,norb), domdr(norb,norb,3)  , tmpvec(norb)
real*8, parameter :: pi=3.14159265358979323846264338d0
real*8  tt, t1, t2, t3, t4, t5, t
real*8  dij
integer i1,i2,i3, n1,n2,n3
real*8  r0xyz(3,nat),tau(3),r2cutoff,rcm(3),tmpxyz(3)
real*8  factor,evalx
integer i,j
real*8  ddot
real*8  d2,wi,wj, cell2(3,3),celleval(3)
integer im, ifngprt, jfngprt, ii,jj
integer,parameter ::  ns=2 !# of s-type GTO's (length of array rwidth),  
real*8 rwidth(ns)

rwidth(1)=w1
rwidth(2)=w2
 
if(sum(cell(:,:)**2)==0.d0) then  ! non-periodic case
  periodic=.false.
else
  periodic=.true.
endif

! exp(-.5 (r_cutoff/rcov)**2)<10^-20  ==> r_cutoff**2 > 40*log(10)*rcov_max**2
!r2cutoff=40.d0*log(10.d0)* (maxval(rcov(:))*max(w1,w2))**2
r2cutoff=20.d0*log(10.d0)* (maxval(rcov(:))*max(w1,w2))**2

if (periodic)  then
  n1=nint(sqrt(r2cutoff/dot_product(cell(:,1),cell(:,1))))
  rtmp(1:3)=cell(:,2) - cell(:,1)*dot_product(cell(:,2),cell(:,1))/dot_product(cell(:,1),cell(:,1) ) !normalize to cell(:,1)
  n2=nint(sqrt(r2cutoff/dot_product(rtmp,rtmp))) 
  rtmp(1:3)=cell(:,3) - rtmp(1:3)*dot_product(cell(:,3),rtmp(1:3))/dot_product(rtmp(1:3),rtmp(1:3) ) &  !normalize to  cell(:,2) already normalized to cell(:1)
                 &    - cell(:,1)*dot_product(cell(:,3),cell(:,1))/dot_product(cell(:,1),cell(:,1) )    !and to cell(:,1)
  n3=nint(sqrt(r2cutoff/dot_product(rtmp,rtmp)))  
call n_rep_dim(cell,sqrt(r2cutoff),n1,n2,n3)
else
  n1=0; n2=0; n3=0
endif

!print*, sqrt(r2cutoff) ,n1,n2,n3 

!call DGEMM ('T','N',3,3,3, 1.d0,cell,3,cell,3,0.d0,cell2,3)
!call diagonalizematrix (3,cell2,celleval) 
 
!if (periodic)
!  n1=nint(sqrt(r2cutoff/celleval(1)))  ! + 1 
!  n2=nint(sqrt(r2cutoff/celleval(2)))  ! + 1 
!  n3=nint(sqrt(r2cutoff/celleval(3)))  ! + 1 
!else
!  n1=0; n2=0; n3=0
!endif
!
!print*, sqrt(r2cutoff) ,n1,n2,n3 ; stop


print*, " making super cluster; numbers of images: ",n1,n2,n3

nsph=0
do i1= -n1,n1
do i2= -n2,n2
do i3= -n3,n3

   
   tau(:)=cell(:,1)*dble(i1) + cell(:,2)*dble(i2) + cell(:,3)*dble(i3)
   !rtmp= [i1,i2,i3] ; call dgemv ('N', 3, 3, 1.d0, cell, 3, rtmp, 1, 0.d0, tau,1)
   !print*, tau; stop 
   
   ! check if the atom should be considered
   do iat=1,nat
      do jat=1,nat
          rtmp(1:3) = rxyz(1:3,iat) + tau (1:3) - rxyz(1:3,jat)
         if( (rtmp(1)**2+rtmp(2)**2+rtmp(3)**2) < r2cutoff) goto 1000
      enddo
      cycle ! the atom is not within any of the spheres centered at atoms
   
   1000 continue  ! consider this atom, becasue it is in the sphere
      
      nsph = nsph + 1
      isph = nsph
   !print*,isph
      if(nsph>nsphx) stop 'nsph>nsphx' 
      rsph (1:3,isph) = rxyz(1:3,iat) + tau (1:3)
      isph2iat(isph)=iat
      if(i1==0 .and. i2==0 .and. i3==0) iat2sph(iat)=isph
   enddo
   
    
enddo  ! i3 
enddo  ! i2
enddo  ! i1

!do iat=1,nat
!  do isph=1,nsph
!     print*, sqrt(sum((rxyz(:,iat)-rsph(:,isph))**2))
!  enddo
!enddo
!return
!print*, rsph(:,1:nsph)
!stop
!
!print*,  nsph
!print*, rsph(:,1:nsph) * 0.52917720859d0
!stop



 print*, " Constructing the overlap matrix for ",nsph, " atoms within the sphere ... "
! 2- setup the overlap matrix for the atoms within the sphere

norb= nsph*ns

allocate (stmp(norb,norb))

!  <s|s>
do jat=1,nsph
xj=rsph(1,jat) ; yj=rsph(2,jat); zj=rsph(3,jat)
do js =1,ns
      jorb=jat+(js-1)*nsph
      aj= rcov(isph2iat(jat))*rwidth(js)

      do iat=1  ,nsph
      xi=rsph(1,iat) ; yi=rsph(2,iat); zi=rsph(3,iat)
      xij=xi-xj; yij=yi-yj; zij=zi-zj
      r2=xij**2 + yij**2 + zij**2 
      do is =1,ns
        iorb=iat+(is-1)*nsph
        ai= rcov(isph2iat(iat))*rwidth(is)
        ! normalized GTOs:
        t1=.5d0/(ai**2 + aj**2)
        stmp(jorb,iorb)=sqrt(4.d0*t1*ai*aj)**3 * exp(-r2*t1)

     enddo
     enddo
enddo  
enddo  



!allocate (eval(norb))
!call diagonalizematrix(norb,stmp,eval)
!print '(es)', eval
!deallocate(eval)
!return
!!! stop

allocate (sijsph(norb,norb))


if(ns .ne. 2) stop 'ns /= 2'
! orthogonalizing


t=sqrt(2.d0*rwidth(1)*rwidth(2)/(rwidth(1)**2 + rwidth(2)**2 ) )**3
tt=sqrt(1.d0-t*t)
if(t==1.d0) stop 'cannot orthogonalize '
do jat=1,nsph
      do iat=1,nsph
        sij = stmp(iat,jat) 

        sijsph(iat     ,jat     ) =  sij
        sijsph(iat+nsph,jat     ) = (stmp(iat+nsph,jat     )- sij*t)/tt 
        sijsph(iat     ,jat+nsph) = (stmp(iat     ,jat+nsph)- sij*t)/tt !
        sijsph(iat+nsph,jat+nsph) = (stmp(iat+nsph,jat+nsph)+ sij*t*t -t*(stmp(iat+nsph,jat) + stmp(iat,jat+nsph)))/tt**2 
     enddo
enddo  

deallocate(stmp)


!do iorb=1,norb
!print '(100f12.7)' , sijsph(:,iorb)
!enddo
!stop

allocate (sijat(norb,norb))
allocate (spow(norb,norb), stmp(norb,norb))

print '(" calcaulating the  powers ... " )' 

do iatcell = 1, nat

  if (periodic) then 
      do jat=1,nsph
      rtmp = rxyz(:,iatcell) - rsph(:,jat)
      ! d2=dot_product(rtmp,rtmp)
      d2=rtmp(1)**2+rtmp(2)**2+rtmp(3)**2 

!      wj=exp(-0.5d0* d2/r2cutoff)

      d2=d2/r2cutoff
      if(d2<1.d0) then 
          wj=(1.d0-d2/2)*(1.d0-d2/2)
      else
          wj=0.d0
      endif

      do js=1,ns
            jorb=jat+(js-1)*nsph
            do iat=1,nsph
            do is=1,ns
              iorb=iat+(is-1)*nsph
              rtmp = rxyz(:,iatcell) - rsph(:,iat)
              d2=rtmp(1)**2+rtmp(2)**2+rtmp(3)**2 

              !wi=exp(-0.5d0*d2/r2cutoff)

              d2=d2/r2cutoff
              if(d2<1.d0) then 
                  wi=(1.d0-d2/2)*(1.d0-d2/2)
              else
                  wi=0.d0
              endif
              sijat(iorb,jorb)=sijsph(iorb,jorb)   * wi*wj
           enddo
           enddo
      enddo  
      enddo  
  else  !wi=1;wj=1  ! for clusters
      sijat(:,:)=sijsph(:,:)
 endif


!print*, " ================  iat=",iatcell

! spow(:,:)=0.d0
! do iorb=1,norb
!   spow(iorb,iorb)=1.d0
! enddo

allocate (eval(norb))
stmp=sijat
call diagonalizematrix(norb,stmp,eval)
!print '("lambda_max= ", es)', eval(norb)
evalx=eval(norb)
deallocate (eval)

spow=sijat
do im=1,m
print '(i3)', im
   call DGEMM ('n', 'n', norb, norb , norb, 1.d0, spow, norb, sijat, norb, 0.d0, stmp, norb )

   spow=stmp


!print*,1+(im-1)*3,iatcell,iat2sph(iatcell)      ,iat2sph(iatcell)
!print*,2+(im-1)*3,iatcell,iat2sph(iatcell)+nsph ,iat2sph(iatcell)
!print*,3+(im-1)*3,iatcell,iat2sph(iatcell)+nsph ,iat2sph(iatcell)+nsph
factor=1.d0/(evalx**im * dble(im))
fngprt(1+(im-1)*3,iatcell)=spow(iat2sph(iatcell)      ,iat2sph(iatcell))       * factor    
fngprt(2+(im-1)*3,iatcell)=spow(iat2sph(iatcell)+nsph ,iat2sph(iatcell))       * factor
fngprt(3+(im-1)*3,iatcell)=spow(iat2sph(iatcell)+nsph ,iat2sph(iatcell)+nsph)  * factor

!i=(im-1)*3; print*, fngprt(1+i:3+i,iatcell),i
enddo   ! im
!stop ! for cluster
enddo ! iatcell
print '(" done ")' 

end subroutine xyz2sm

!!! =============================================
!!subroutine xyz2sm(nat, cell,rxyz,rcov,w1,w2,m,fngprt)
!!implicit none
!!! in/out
!!real*8  cell(3,3)
!!integer nat
!!real*8  rxyz(3,nat), rcov(nat)
!!real*8  w1 , w2  ! relative  widths of the GTO's on each atom (e.g. use w1=1, w2=1.5)
!!integer m  !  m = maximum power
!!real*8  fngprt(3*m, nat)  ! Fingerprint array:  2 s-GTOs per atom ==> a symmetric 2x2 matrix per atom which has 3 independent elements 
!!
!!! Local variables:
!!integer norb 
!!integer isph, jsph, nsph  ! # counter, and # of atoms in the sphere
!!integer, parameter :: nsphx=10000   ! max of nsphr; later on should be estimated on-the-fly
!!real*8 rtmp(3), rsph(3,nsphx)  ! rxyz in the sphere
!!integer  isph2iat(nsphx) , iat2sph(nat)
!!real*8,allocatable ::  sijsph(:,:), sijat(:,:), spow(:,:), stmp(:,:)
!!real*8,allocatable :: eval(:)
!!
!!integer iat,jat,is,js,ip,jp, kxyz, kat, iorb, jorb, lat
!!integer iatcell
!!logical  periodic
!!real*8  ai,aj, xi,yi,zi, xj,yj,zj, xij, yij, zij, r2 ,sij, dtau 
!!!real*8  ovrlp(norb,norb), evec(norb,norb), domdr(norb,norb,3)  , tmpvec(norb)
!!real*8, parameter :: pi=3.14159265358979323846264338d0
!!real*8  tt, t1, t2, t3, t4, t5, t
!!real*8  dij
!!integer i1,i2,i3, n1,n2,n3
!!real*8  r0xyz(3,nat),tau(3),r2cutoff,rcm(3),tmpxyz(3)
!!integer i,j
!!real*8  ddot
!!real*8  d2,wi,wj, cell2(3,3),celleval(3)
!!integer im, ifngprt, jfngprt, ii,jj
!!integer,parameter ::  ns=2 !# of s-type GTO's (length of array rwidth),  
!!real*8 rwidth(ns)
!!
!!rwidth(1)=w1
!!rwidth(2)=w2
!! 
!!if(sum(cell(:,:)**2)==0.d0) then  ! non-periodic case
!!  periodic=.false.
!!else
!!  periodic=.true.
!!endif
!!
!!! exp(-.5 (r_cutoff/rcov)**2)<10^-20  ==> r_cutoff**2 > 40*log(10)*rcov_max**2
!!!r2cutoff=40.d0*log(10.d0)* (maxval(rcov(:))*max(w1,w2))**2
!!r2cutoff=20.d0*log(10.d0)* (maxval(rcov(:))*max(w1,w2))**2
!!
!!if (periodic)  then
!!  call n_rep_dim(cell,sqrt(r2cutoff),n1,n2,n3)
!!!!  n1=nint(sqrt(r2cutoff/dot_product(cell(:,1),cell(:,1))))
!!!!  rtmp(1:3)=cell(:,2) - cell(:,1)*dot_product(cell(:,2),cell(:,1))/dot_product(cell(:,1),cell(:,1) ) !normalize to cell(:,1)
!!!!  n2=nint(sqrt(r2cutoff/dot_product(rtmp,rtmp))) 
!!!!  rtmp(1:3)=cell(:,3) - rtmp(1:3)*dot_product(cell(:,3),rtmp(1:3))/dot_product(rtmp(1:3),rtmp(1:3) ) &  !normalize to  cell(:,2) already normalized to cell(:1)
!!!!                 &    - cell(:,1)*dot_product(cell(:,3),cell(:,1))/dot_product(cell(:,1),cell(:,1) )    !and to cell(:,1)
!!!!  n3=nint(sqrt(r2cutoff/dot_product(rtmp,rtmp)))  
!!else
!!  n1=0; n2=0; n3=0
!!endif
!!
!!!print*, sqrt(r2cutoff) ,n1,n2,n3 
!!
!!!call DGEMM ('T','N',3,3,3, 1.d0,cell,3,cell,3,0.d0,cell2,3)
!!!call diagonalizematrix (3,cell2,celleval) 
!! 
!!!if (periodic)
!!!  n1=nint(sqrt(r2cutoff/celleval(1)))  ! + 1 
!!!  n2=nint(sqrt(r2cutoff/celleval(2)))  ! + 1 
!!!  n3=nint(sqrt(r2cutoff/celleval(3)))  ! + 1 
!!!else
!!!  n1=0; n2=0; n3=0
!!!endif
!!!
!!!print*, sqrt(r2cutoff) ,n1,n2,n3 ; stop
!!
!!
!!print*, " making super cluster; numbers of images: ",n1,n2,n3
!!
!!nsph=0
!!do i1= -n1,n1
!!do i2= -n2,n2
!!do i3= -n3,n3
!!
!!   
!!   tau(:)=cell(:,1)*dble(i1) + cell(:,2)*dble(i2) + cell(:,3)*dble(i3)
!!   !rtmp= [i1,i2,i3] ; call dgemv ('N', 3, 3, 1.d0, cell, 3, rtmp, 1, 0.d0, tau,1)
!!   !print*, tau; stop 
!!   
!!   ! check if the atom should be considered
!!   do iat=1,nat
!!      do jat=1,nat
!!          rtmp(1:3) = rxyz(1:3,iat) + tau (1:3) - rxyz(1:3,jat)
!!         if( (rtmp(1)**2+rtmp(2)**2+rtmp(3)**2) < r2cutoff) goto 1000
!!      enddo
!!      cycle ! the atom is not within any of the spheres centered at atoms
!!   
!!   1000 continue  ! consider this atom, becasue it is in the sphere
!!      
!!      nsph = nsph + 1
!!      isph = nsph
!!   !print*,isph
!!      if(nsph>nsphx) stop 'nsph>nsphx' 
!!      rsph (1:3,isph) = rxyz(1:3,iat) + tau (1:3)
!!      isph2iat(isph)=iat
!!      if(i1==0 .and. i2==0 .and. i3==0) iat2sph(iat)=isph
!!   enddo
!!   
!!    
!!enddo  ! i3 
!!enddo  ! i2
!!enddo  ! i1
!!
!!!do iat=1,nat
!!!  do isph=1,nsph
!!!     print*, sqrt(sum((rxyz(:,iat)-rsph(:,isph))**2))
!!!  enddo
!!!enddo
!!!return
!!!print*, rsph(:,1:nsph)
!!!stop
!!!
!!!print*,  nsph
!!!print*, rsph(:,1:nsph) * 0.52917720859d0
!!!stop
!!
!!
!!
!! print*, " Constructing the overlap matrix for ",nsph, " atoms within the sphere ... "
!!! 2- setup the overlap matrix for the atoms within the sphere
!!
!!norb= nsph*ns
!!
!!allocate (stmp(norb,norb))
!!
!!!  <s|s>
!!do jat=1,nsph
!!xj=rsph(1,jat) ; yj=rsph(2,jat); zj=rsph(3,jat)
!!do js =1,ns
!!      jorb=jat+(js-1)*nsph
!!      aj= rcov(isph2iat(jat))*rwidth(js)
!!
!!      do iat=1  ,nsph
!!      xi=rsph(1,iat) ; yi=rsph(2,iat); zi=rsph(3,iat)
!!      xij=xi-xj; yij=yi-yj; zij=zi-zj
!!      r2=xij**2 + yij**2 + zij**2 
!!      do is =1,ns
!!        iorb=iat+(is-1)*nsph
!!        ai= rcov(isph2iat(iat))*rwidth(is)
!!        ! normalized GTOs:
!!        t1=.5d0/(ai**2 + aj**2)
!!        stmp(jorb,iorb)=sqrt(4.d0*t1*ai*aj)**3 * exp(-r2*t1)
!!
!!     enddo
!!     enddo
!!enddo  
!!enddo  
!!
!!
!!
!!!allocate (eval(norb))
!!!call diagonalizematrix(norb,stmp,eval)
!!!print '(es)', eval
!!!deallocate(eval)
!!!return
!!!!! stop
!!
!!allocate (sijsph(norb,norb))
!!
!!
!!if(ns .ne. 2) stop 'ns /= 2'
!!! orthogonalizing
!!
!!
!!t=sqrt(2.d0*rwidth(1)*rwidth(2)/(rwidth(1)**2 + rwidth(2)**2 ) )**3
!!tt=sqrt(1.d0-t*t)
!!if(t==1.d0) stop 'cannot orthogonalize '
!!do jat=1,nsph
!!      do iat=1,nsph
!!        sij = stmp(iat,jat) 
!!
!!        sijsph(iat     ,jat     ) =  sij
!!        sijsph(iat+nsph,jat     ) = (stmp(iat+nsph,jat     )- sij*t)/tt 
!!        sijsph(iat     ,jat+nsph) = (stmp(iat     ,jat+nsph)- sij*t)/tt !
!!        sijsph(iat+nsph,jat+nsph) = (stmp(iat+nsph,jat+nsph)+ sij*t*t -t*(stmp(iat+nsph,jat) + stmp(iat,jat+nsph)))/tt**2 
!!     enddo
!!enddo  
!!
!!deallocate(stmp)
!!
!!
!!!do iorb=1,norb
!!!print '(100f12.7)' , sijsph(:,iorb)
!!!enddo
!!!stop
!!
!!!allocate (eval(norb))
!!!call diagonalizematrix(norb,sijsph,eval)
!!!print '(es)', eval; stop
!!
!!allocate (sijat(norb,norb))
!!allocate (spow(norb,norb), stmp(norb,norb))
!!
!!print '(" calcaulating the  powers ... " \)' 
!!
!!do iatcell = 1, nat
!!
!!  if (periodic) then 
!!      do jat=1,nsph
!!      rtmp = rxyz(:,iatcell) - rsph(:,jat)
!!      ! d2=dot_product(rtmp,rtmp)
!!      d2=rtmp(1)**2+rtmp(2)**2+rtmp(3)**2 
!!
!!!      wj=exp(-0.5d0* d2/r2cutoff)
!!
!!      d2=d2/r2cutoff
!!      if(d2<1.d0) then 
!!          wj=(1.d0-d2/2)*(1.d0-d2/2)
!!      else
!!          wj=0.d0
!!      endif
!!
!!      do js=1,ns
!!            jorb=jat+(js-1)*nsph
!!            do iat=1,nsph
!!            do is=1,ns
!!              iorb=iat+(is-1)*nsph
!!              rtmp = rxyz(:,iatcell) - rsph(:,iat)
!!              d2=rtmp(1)**2+rtmp(2)**2+rtmp(3)**2 
!!
!!              !wi=exp(-0.5d0*d2/r2cutoff)
!!
!!              d2=d2/r2cutoff
!!              if(d2<1.d0) then 
!!                  wi=(1.d0-d2/2)*(1.d0-d2/2)
!!              else
!!                  wi=0.d0
!!              endif
!!              sijat(iorb,jorb)=sijsph(iorb,jorb)   * wi*wj
!!           enddo
!!           enddo
!!      enddo  
!!      enddo  
!!  else  !wi=1;wj=1  ! for clusters
!!      sijat(:,:)=sijsph(:,:)
!! endif
!!
!!
!!!print*, " ================  iat=",iatcell
!!
!!! spow(:,:)=0.d0
!!! do iorb=1,norb
!!!   spow(iorb,iorb)=1.d0
!!! enddo
!!
!!spow=sijat
!!do im=1,m
!!print '(i3\)', im
!!   call DGEMM ('n', 'n', norb, norb , norb, 1.d0, spow, norb, sijat, norb, 0.d0, stmp, norb )
!!
!!   spow=stmp
!!
!!!print*,1+(im-1)*3,iatcell,iat2sph(iatcell)      ,iat2sph(iatcell)
!!!print*,2+(im-1)*3,iatcell,iat2sph(iatcell)+nsph ,iat2sph(iatcell)
!!!print*,3+(im-1)*3,iatcell,iat2sph(iatcell)+nsph ,iat2sph(iatcell)+nsph
!!fngprt(1+(im-1)*3,iatcell)=spow(iat2sph(iatcell)      ,iat2sph(iatcell))
!!fngprt(2+(im-1)*3,iatcell)=spow(iat2sph(iatcell)+nsph ,iat2sph(iatcell))
!!fngprt(3+(im-1)*3,iatcell)=spow(iat2sph(iatcell)+nsph ,iat2sph(iatcell)+nsph)
!!
!!!i=(im-1)*3; print*, fngprt(1+i:3+i,iatcell),i
!!enddo   ! im
!!!stop ! for cluster
!!enddo ! iatcell
!!print '(" done ")' 
!!
!!end subroutine xyz2sm


subroutine get_distance_xyz2sm(nat,typat,lfp,fp1,fp2,diff)
  implicit none
  integer  :: nat,iat,jat,lfp,typat(nat)
  real(8)  :: cost,r2,fp1(lfp,nat), fp2(lfp,nat),fp2copy(lfp,nat),diff
  INTEGER f(nat)
  REAL*8,allocatable :: A(:,:)

      allocate(A(nat,nat))
      do jat=1,nat
      do iat=1,nat
         if(typat(iat) .ne. typat(jat))  then
             A(iat,jat)=1.d7   ! do not permute!
         else
!         A(iat,jat)=sum(rat1(:,iat)-rat2(:,jat))**2
         A(iat,jat)=sum((fp1(1:lfp,iat)-fp2(1:lfp,jat))**2)/dble(nat*lfp)
!         A(iat,jat)=sum(((fp1(1:lfp,iat)-fp2(1:lfp,jat))/(fp1(1:lfp,iat)+fp2(1:lfp,jat)))**2)
!         A(iat,jat)=1.d0 - (ddot(lfp, fp1(:,iat),1, fp2(:,jat),1))**2/ (ddot(lfp, fp1(:,iat),1, fp1(:,iat),1)*ddot(lfp, fp2(:,jat),1, fp2(:,jat),1))
        endif
      enddo
      enddo

      CALL APC(Nat,A,F,cost)
!m=nat;       call assndx(1,a, nat, m, F, cost)   ! this is slower by a factor of ~3

!!      xat2copy =xat2
!!      rat2copy=rat2
!!      fp2copy =fp2
!!
!!      do iat=1,nat
!!        xat2(iat)  =xat2copy(F(iat))
!!        rat2(1:3,iat)  =rat2copy(1:3,F(iat))
!!        fp2 (1:lfp,iat)=fp2copy (1:lfp,F(iat))
!!      enddo

do iat=1,nat
  write(24,*) iat,F(iat)
enddo 

!!!      allocate(A(nat,nat))
!!!      do jat=1,nat
!!!      do iat=1,nat
!!!         if(typat(iat) .ne. typat(jat)) then
!!!             A(iat,jat)=1.d7
!!!         else
!!!             !A(iat,jat)=sum(((fp1(1:lfp,iat)-fp2(1:lfp,jat))/(fp1(1:lfp,iat)+fp2(1:lfp,jat)))**2)
!!!             !A(iat,jat)=1.d0 - (ddot(lfp, fp1(:,iat),1, fp2(:,jat),1))**2/ (ddot(lfp, fp1(:,iat),1, fp1(:,iat),1)*ddot(lfp, fp2(:,jat),1, fp2(:,jat),1))
!!!             A(iat,jat)=1.d0 - (dot_product(fp1(:,iat),fp2(:,jat)))**2/ (dot_product(fp1(:,iat),fp1(:,iat))*dot_product(fp2(:,jat),fp2(:,jat)))
!!!         endif
!!!      enddo
!!!      enddo
!!!
!!!
!!!      CALL APC(Nat,A,F,cost)
!!!!m=nat;       call assndx(1,a, nat, m, F, cost)   ! this is slower by a factor of ~3
!!!
!!!!      fp2copy =fp2 
!!!
!!!!!write(*,*) "fp1"
!!!!!write(*,*) fp1
!!!!!write(*,*) "fp2"
!!!!!write(*,*) fp2
!!!!!write(*,*) "fp2copy"
!!!!!write(*,*) fp2copy
!!!!      do iat=1,nat
!!!!        fp2copy (1:lfp,iat)=fp2 (1:lfp,F(iat))
!!!!      enddo
!!!!diff=sqrt(sum((fp2copy-fp1)**2)/dble(lfp*nat))
!!!!diff=sqrt(sum(((fp2copy-fp1)/(fp2copy+fp1))**2)/dble(lfp*nat))
diff=cost
!!write(*,*) "diff",diff

end subroutine
!!
!!! =============================================
!!subroutine xyz2sm(nat, cell,rxyz,rcov,w1,w2,m,fngprt)
!!implicit none
!!! in/out
!!real*8  cell(3,3)
!!integer nat
!!real*8  rxyz(3,nat), rcov(nat)
!!real*8  w1 , w2  ! relative  widths of the GTO's on each atom (e.g. use w1=1, w2=1.5)
!!integer m  !  m = maximum power
!!real*8  fngprt(3*m, nat)  ! Fingerprint array:  2 s-GTOs per atom ==> a symmetric 2x2 matrix per atom which has 3 independent elements 
!!
!!! Local variables:
!!integer norb 
!!integer isph, jsph, nsph  ! # counter, and # of atoms in the sphere
!!integer, parameter :: nsphx=10000   ! max of nsphr; later on should be estimated on-the-fly
!!real*8 rtmp(3), rsph(3,nsphx)  ! rxyz in the sphere
!!integer  isph2iat(nsphx) , iat2sph(nat)
!!real*8,allocatable ::  sijsph(:,:), sijat(:,:), spow(:,:), stmp(:,:)
!!real*8,allocatable :: eval(:)
!!
!!integer iat,jat,is,js,ip,jp, kxyz, kat, iorb, jorb, lat
!!integer iatcell
!!logical  periodic
!!real*8  ai,aj, xi,yi,zi, xj,yj,zj, xij, yij, zij, r2 ,sij, dtau 
!!!real*8  ovrlp(norb,norb), evec(norb,norb), domdr(norb,norb,3)  , tmpvec(norb)
!!real*8, parameter :: pi=3.14159265358979323846264338d0
!!real*8  tt, t1, t2, t3, t4, t5, t
!!real*8  dij
!!integer i1,i2,i3, n1,n2,n3
!!real*8  r0xyz(3,nat),tau(3),r2cutoff,rcm(3),tmpxyz(3)
!!integer i,j
!!real*8  ddot
!!real*8  d2,wi,wj, cell2(3,3),celleval(3)
!!integer im, ifngprt, jfngprt, ii,jj
!!integer,parameter ::  ns=2 !# of s-type GTO's (length of array rwidth),  
!!real*8 rwidth(ns)
!!
!!rwidth(1)=w1
!!rwidth(2)=w2
!! 
!!if(sum(cell(:,:)**2)==0.d0) then  ! non-periodic case
!!  periodic=.false.
!!else
!!  periodic=.true.
!!endif
!!
!!! exp(-.5 (r_cutoff/rcov)**2)<10^-20  ==> r_cutoff**2 > 40*log(10)*rcov_max**2
!!!r2cutoff=40.d0*log(10.d0)* (maxval(rcov(:))*max(w1,w2))**2
!!r2cutoff=20.d0*log(10.d0)* (maxval(rcov(:))*max(w1,w2))**2
!!
!!if (periodic)  then
!!  n1=nint(sqrt(r2cutoff/sum(cell(:,1)**2)))  ! + 1 
!!  n2=nint(sqrt(r2cutoff/sum(cell(:,2)**2)))  ! + 1 
!!  n3=nint(sqrt(r2cutoff/sum(cell(:,3)**2)))  ! + 1 
!!  call n_rep_dim(cell,sqrt(r2cutoff),n1,n2,n3)
!!  n1=n1+2
!!  n2=n2+2
!!  n3=n3+2
!!else
!!  n1=0; n2=0; n3=0
!!endif
!!
!!!print*, sqrt(r2cutoff) ,n1,n2,n3 
!!
!!!call DGEMM ('T','N',3,3,3, 1.d0,cell,3,cell,3,0.d0,cell2,3)
!!!call diagonalizematrix (3,cell2,celleval) 
!! 
!!!if (periodic)
!!!  n1=nint(sqrt(r2cutoff/celleval(1)))  ! + 1 
!!!  n2=nint(sqrt(r2cutoff/celleval(2)))  ! + 1 
!!!  n3=nint(sqrt(r2cutoff/celleval(3)))  ! + 1 
!!!else
!!!  n1=0; n2=0; n3=0
!!!endif
!!!
!!!print*, sqrt(r2cutoff) ,n1,n2,n3 ; stop
!!
!!write(*,*) cell(:,1)
!!write(*,*) cell(:,2)
!!write(*,*) cell(:,3)
!!do iat=1,nat
!!write(*,*) rxyz(1:3,iat)
!!enddo
!!return
!!print*, " making super cluster; numbers of images: ",n1,n2,n3
!!
!!nsph=0
!!do i1= -n1,n1
!!do i2= -n2,n2
!!do i3= -n3,n3
!!
!!   
!!   tau(:)=cell(:,1)*dble(i1) + cell(:,2)*dble(i2) + cell(:,3)*dble(i3)
!!   !rtmp= [i1,i2,i3] ; call dgemv ('N', 3, 3, 1.d0, cell, 3, rtmp, 1, 0.d0, tau,1)
!!   !print*, tau; stop 
!!   
!!   ! check if the atom should be considered
!!   do iat=1,nat
!!      do jat=1,nat
!!          rtmp(1:3) = rxyz(1:3,iat) + tau (1:3) - rxyz(1:3,jat)
!!         if( (rtmp(1)**2+rtmp(2)**2+rtmp(3)**2) < r2cutoff) goto 1000
!!      enddo
!!      cycle ! the atom is not within any of the spheres centered at atoms
!!   
!!   1000 continue  ! consider this atom, becasue it is in the sphere
!!      
!!      nsph = nsph + 1
!!      isph = nsph
!!   !print*,isph
!!      if(nsph>nsphx) stop 'nsph>nsphx' 
!!      rsph (1:3,isph) = rxyz(1:3,iat) + tau (1:3)
!!      isph2iat(isph)=iat
!!      if(i1==0 .and. i2==0 .and. i3==0) iat2sph(iat)=isph
!!   enddo
!!   
!!    
!!enddo  ! i3 
!!enddo  ! i2
!!enddo  ! i1
!!
!!!do iat=1,nat
!!!  do isph=1,nsph
!!!     print*, sqrt(sum((rxyz(:,iat)-rsph(:,isph))**2))
!!!  enddo
!!!enddo
!!!return
!!!print*, rsph(:,1:nsph)
!!!stop
!!!
!!!print*,  nsph
!!!print*, rsph(:,1:nsph)
!!!stop
!!!
!!
!!
!! print*, " Constructing the overlap matrix for ",nsph, " atoms within the sphere ... "
!!! 2- setup the overlap matrix for the atoms within the sphere
!!
!!norb= nsph*ns
!!
!!allocate (stmp(norb,norb))
!!
!!!  <s|s>
!!do jat=1,nsph
!!xj=rsph(1,jat) ; yj=rsph(2,jat); zj=rsph(3,jat)
!!do js =1,ns
!!      jorb=jat+(js-1)*nsph
!!      aj= rcov(isph2iat(jat))*rwidth(js)
!!
!!      do iat=1  ,nsph
!!      xi=rsph(1,iat) ; yi=rsph(2,iat); zi=rsph(3,iat)
!!      xij=xi-xj; yij=yi-yj; zij=zi-zj
!!      r2=xij**2 + yij**2 + zij**2 
!!      do is =1,ns
!!        iorb=iat+(is-1)*nsph
!!        ai= rcov(isph2iat(iat))*rwidth(is)
!!        ! normalized GTOs:
!!        t1=.5d0/(ai**2 + aj**2)
!!        stmp(jorb,iorb)=sqrt(4.d0*t1*ai*aj)**3 * exp(-r2*t1)
!!
!!     enddo
!!     enddo
!!enddo  
!!enddo  
!!
!!
!!
!!!allocate (eval(norb))
!!!call diagonalizematrix(norb,stmp,eval)
!!!print '(es)', eval
!!!deallocate(eval)
!!!return
!!!!! stop
!!
!!allocate (sijsph(norb,norb))
!!
!!
!!if(ns .ne. 2) stop 'ns /= 2'
!!! orthogonalizing
!!
!!
!!t=sqrt(2.d0*rwidth(1)*rwidth(2)/(rwidth(1)**2 + rwidth(2)**2 ) )**3
!!tt=sqrt(1.d0-t*t)
!!if(t==1.d0) stop 'cannot orthogonalize '
!!do jat=1,nsph
!!      do iat=1,nsph
!!        sij = stmp(iat,jat) 
!!
!!        sijsph(iat     ,jat     ) =  sij
!!        sijsph(iat+nsph,jat     ) = (stmp(iat+nsph,jat     )- sij*t)/tt 
!!        sijsph(iat     ,jat+nsph) = (stmp(iat     ,jat+nsph)- sij*t)/tt !
!!        sijsph(iat+nsph,jat+nsph) = (stmp(iat+nsph,jat+nsph)+ sij*t*t -t*(stmp(iat+nsph,jat) + stmp(iat,jat+nsph)))/tt**2 
!!     enddo
!!enddo  
!!
!!deallocate(stmp)
!!
!!
!!!do iorb=1,norb
!!!print '(100f12.7)' , sijsph(:,iorb)
!!!enddo
!!!stop
!!
!!!allocate (eval(norb))
!!!call diagonalizematrix(norb,sijsph,eval)
!!!print '(es)', eval; stop
!!
!!allocate (sijat(norb,norb))
!!allocate (spow(norb,norb), stmp(norb,norb))
!!
!!print '(" calcaulating the  powers ... " \)' 
!!
!!do iatcell = 1, nat
!!
!!  if (periodic) then 
!!      do jat=1,nsph
!!      rtmp = rxyz(:,iatcell) - rsph(:,jat)
!!      ! d2=dot_product(rtmp,rtmp)
!!      d2=rtmp(1)**2+rtmp(2)**2+rtmp(3)**2 
!!      wj=exp(-0.5d0* d2/r2cutoff)
!!      do js=1,ns
!!            jorb=jat+(js-1)*nsph
!!            do iat=1,nsph
!!            do is=1,ns
!!              iorb=iat+(is-1)*nsph
!!              rtmp = rxyz(:,iatcell) - rsph(:,iat)
!!              d2=rtmp(1)**2+rtmp(2)**2+rtmp(3)**2 
!!              wi=exp(-0.5d0*d2/r2cutoff)
!!              sijat(iorb,jorb)=sijsph(iorb,jorb)   * wi*wj
!!           enddo
!!           enddo
!!      enddo  
!!      enddo  
!!  else  !wi=1;wj=1  ! for clusters
!!      sijat(:,:)=sijsph(:,:)
!! endif
!!
!!
!!!print*, " ================  iat=",iatcell
!!
!!! spow(:,:)=0.d0
!!! do iorb=1,norb
!!!   spow(iorb,iorb)=1.d0
!!! enddo
!!
!!spow=sijat
!!do im=1,m
!!print '(i3\)', im
!!   call DGEMM ('n', 'n', norb, norb , norb, 1.d0, spow, norb, sijat, norb, 0.d0, stmp, norb )
!!
!!   spow=stmp
!!
!!!print*,1+(im-1)*3,iatcell,iat2sph(iatcell)      ,iat2sph(iatcell)
!!!print*,2+(im-1)*3,iatcell,iat2sph(iatcell)+nsph ,iat2sph(iatcell)
!!!print*,3+(im-1)*3,iatcell,iat2sph(iatcell)+nsph ,iat2sph(iatcell)+nsph
!!fngprt(1+(im-1)*3,iatcell)=spow(iat2sph(iatcell)      ,iat2sph(iatcell))
!!fngprt(2+(im-1)*3,iatcell)=spow(iat2sph(iatcell)+nsph ,iat2sph(iatcell))
!!fngprt(3+(im-1)*3,iatcell)=spow(iat2sph(iatcell)+nsph ,iat2sph(iatcell)+nsph)
!!
!!!i=(im-1)*3; print*, fngprt(1+i:3+i,iatcell),i
!!enddo   ! im
!!!stop ! for cluster
!!enddo ! iatcell
!!print '(" done ")' 
!!
!!end subroutine xyz2sm









subroutine  diagonalizeMatrix(n, mat, eval)
  implicit none
  
  ! Calling arguments
  integer,intent(in):: n
  real(8),dimension(n,n),intent(inout):: mat
  real(8),dimension(n),intent(out):: eval
  
  ! Local variables
  integer:: lwork, info
  real(8),dimension(:),allocatable:: work
  
  lwork=100*n
  allocate(work(lwork))
  call dsyev('n','l', n, mat, n, eval, work, lwork, info)
  if(info/=0) stop ' ERROR in dsyev'
  deallocate(work)

end subroutine diagonalizeMatrix

!=====================================================================================================

subroutine hungarian(nat,rat1,rat2,lfp,fp1,fp2,cost,xat1,xat2)
implicit none
integer  :: nat,iat,jat,lfp 
  real(8),dimension(3,nat)         :: rat1,rat2,rat2copy       ! atomic position
  real(8)                          :: cost,r2,fp1(lfp,nat), fp2(lfp,nat),fp2copy(lfp,nat)   
  INTEGER f(nat)
  REAL*8,allocatable :: A(:,:)
  real*8  :: ddot
character(2) xat1(nat),xat2(nat),xat2copy(nat)


!print*,xat1; stop

      allocate(A(nat,nat))
      do jat=1,nat
      do iat=1,nat
         if(xat1(iat) .ne. xat2(jat))  then 
             A(iat,jat)=1.d7   ! do not permute!
         else
!         A(iat,jat)=sum(rat1(:,iat)-rat2(:,jat))**2
         A(iat,jat)=sum((fp1(1:lfp,iat)-fp2(1:lfp,jat))**2)/dble(nat*lfp)
!         A(iat,jat)=sum(((fp1(1:lfp,iat)-fp2(1:lfp,jat))/(fp1(1:lfp,iat)+fp2(1:lfp,jat)))**2)
!         A(iat,jat)=1.d0 - (ddot(lfp, fp1(:,iat),1, fp2(:,jat),1))**2/ (ddot(lfp, fp1(:,iat),1, fp1(:,iat),1)*ddot(lfp, fp2(:,jat),1, fp2(:,jat),1))
        endif
      enddo
      enddo

      CALL APC(Nat,A,F,cost)
!m=nat;       call assndx(1,a, nat, m, F, cost)   ! this is slower by a factor of ~3

      xat2copy =xat2 
      rat2copy=rat2
      fp2copy =fp2 

      do iat=1,nat
        xat2(iat)  =xat2copy(F(iat))
        rat2(1:3,iat)  =rat2copy(1:3,F(iat))
        fp2 (1:lfp,iat)=fp2copy (1:lfp,F(iat))
      enddo

!print '(i)',f(1:nat)
end subroutine hungarian



SUBROUTINE APC(N,A,F,Z)
implicit none
! Modified by Ali Sadeghi to get real*8 matrix A(N,N) and converted to F90 
!
! SOLUTION OF THE LINEAR MIN-SUM ASSIGNMENT PROBLEM.
! HUNGARIAN METHOD. COMPLEXITY O(N**3).
!
! MEANING OF THE INPUT PARAMETERS:
! N      = NUMBER OF ROWS AND COLUMNS OF THE COST MATRIX.
! A(I,J) = COST OF THE ASSIGNMENT OF ROW  I  TO COLUMN  J .
! ON RETURN, THE INPUT PARAMETERS ARE UNCHANGED.
!
! MEANING OF THE OUTPUT PARAMETERS:
! F(I) = COLUMN ASSIGNED TO ROW  I .
! Z    = COST OF THE OPTIMAL ASSIGNMENT =
!      = A(1,F(1)) + A(2,F(2)) + ... + A(N,F(N)) .
!
!
! THE CODE IS BASED ON THE HUNGARIAN METHOD AS DESCRIBED BY
! LAWLER (COMBINATORIAL OPTIMIZATION : NETWORKS AND
! MATROIDS, HOLT, RINEHART AND WINSTON, NEW YORK, 1976).
! THE ALGORITHMIC PASCAL-LIKE DESCRIPTION OF THE CODE IS
! GIVEN IN G.CARPANETO, S.MARTELLO AND P.TOTH, ALGORITHMS AND
! CODES FOR THE ASSIGNMENT PROBLEM, ANNALS OF OPERATIONS
! RESEARCH 7, 1988.
!
! SUBROUTINE APC DETERMINES THE INITIAL DUAL AND PARTIAL
! PRIMAL SOLUTIONS AND THEN SEARCHES FOR AUGMENTING PATHS
! UNTIL ALL ROWS AND COLUMNS ARE ASSIGNED.
!
! MEANING OF THE MAIN INTERNAL VARIABLES:
! FB(J) = ROW ASSIGNED TO COLUMN  J .
! M     = NUMBER OF INITIAL ASSIGNMENTS.
! U(I)  = DUAL VARIABLE ASSOCIATED WITH ROW  I .
! V(J)  = DUAL VARIABLE ASSOCIATED WITH COLUMN  J .
!
! APC NEEDS THE FOLLOWING SUBROUTINES: INCR
!                                      INIT
!                                      PATH
!
! THIS WORK WAS SUPPORTED BY  C.N.R. , ITALY.
      INTEGER N
      REAL*8  A(n,n),Z,U(n),V(n)
      integer F(N),FB(n), RC(n)
      INTEGER M,I,J,K
! SEARCH FOR THE INITIAL DUAL AND PARTIAL PRIMAL SOLUTIONS.
      CALL INIT(N,A,F,M,U,V,FB,RC)
      IF ( M .NE. N ) then 
! SOLUTION OF THE REDUCED PROBLEM.
      DO  I=1,N
        IF ( F(I) == 0 ) THEN 
! DETERMINATION OF AN AUGMENTING PATH STARTING FROM ROW  I .
        CALL PATH(N,A,I,F,J,U,V,FB,RC)
! ASSIGNMENT OF ROW  I  AND COLUMN  J .
        CALL INCR(n,F,J,FB,RC)
        ENDIF
      ENDDO    
      ENDIF

! COMPUTATION OF THE SOLUTION COST  Z .
      Z = sum(u(1:N)) + sum(V(1:N))
      END


      SUBROUTINE INCR(n,F,J,FB,RC)
!
! ASSIGNMENT OF COLUMN  J .
!
      INTEGER n,I,J,JJ,  F(n),FB(n),RC(n)
   10 I = RC(J)
      FB(J) = I
      JJ = F(I)
      F(I) = J
      J = JJ
      IF ( J > 0 ) GO TO 10
      RETURN
      END


      SUBROUTINE INIT(N,A,F,M,U,V,FB,P)
!
! SEARCH FOR THE INITIAL DUAL AND PARTIAL PRIMAL SOLUTIONS.
!
! P(I) = FIRST UNSCANNED COLUMN OF ROW  I .
       IMPLICIT NONE
!
      INTEGER n,m, F(n),FB(n),P(n)
      real*8 A(n,n) , U(n),V(n)
      REAL*8, parameter :: INF = 1.d9
      real*8 min, IA
      integer i,j, k,R, JMIN, KK
! PHASE 1 .
      M = 0
      F(1:N)=0
      FB(1:N)=0
! SCANNING OF THE COLUMNS ( INITIALIZATION OF  V(J) ).
      DO 40 J=1,N
        MIN = INF
        DO 30 I=1,N
          IA = A(I,J)
          IF ( IA .GT. MIN ) GO TO 30
          IF ( IA .LT. MIN ) GO TO 20
          IF ( F(I) .NE. 0 ) GO TO 30
   20     MIN = IA
          R = I
   30   CONTINUE
        V(J) = MIN
        IF ( F(R) .NE. 0 ) GO TO 40
! ASSIGNMENT OF COLUMN  J  TO ROW  R .
        M = M + 1
        FB(J) = R
        F(R) = J
        U(R) = 0.d0
        P(R) = J + 1
   40 CONTINUE
! PHASE 2 .
! SCANNING OF THE UNASSIGNED ROWS ( UPDATING OF  U(I) ).
      DO 110 I=1,N
        IF ( F(I) .NE. 0 ) GO TO 110
        MIN = INF
        DO 60 K=1,N
          IA = A(I,K) - V(K)
          IF ( IA .GT. MIN )  GO TO 60
          IF ( IA .LT. MIN )  GO TO 50
          IF ( FB(K) .NE. 0 ) GO TO 60
          IF ( FB(J) .EQ. 0 ) GO TO 60
   50     MIN = IA
          J = K
   60   CONTINUE
        U(I) = MIN
        JMIN = J
        IF ( FB(J) .EQ. 0 ) GO TO 100
        DO 80 J=JMIN,N
          IF ( A(I,J) - V(J) .GT. MIN ) GO TO 80
          R = FB(J)
          KK = P(R)
          IF ( KK .GT. N ) GO TO 80
          DO 70 K=KK,N
            IF ( FB(K) .GT. 0 ) GO TO 70
            IF ( A(R,K) - U(R) - V(K) .EQ. 0.d0 ) GO TO 90
   70     CONTINUE
          P(R) = N + 1
   80   CONTINUE
        GO TO 110
! REASSIGNMENT OF ROW  R  AND COLUMN  K .
   90   F(R) = K
        FB(K) = R
        P(R) = K + 1
! ASSIGNMENT OF COLUMN  J  TO ROW  I .
  100   M = M + 1
        F(I) = J
        FB(J)= I
        P(I) = J + 1
  110 CONTINUE
      RETURN
      END
      SUBROUTINE PATH(N,A,II,F,JJ,U,V,FB,RC)
!
! DETERMINATION OF AN AUGMENTING PATH STARTING FROM
! UNASSIGNED ROW  II  AND TERMINATING AT UNASSIGNED COLUMN
! JJ , WITH UPDATING OF DUAL VARIABLES  U(I)  AND  V(J) .
!
! MEANING OF THE MAIN INTERNAL VARIABLES:
! LR(L) = L-TH LABELLED ROW ( L=1,NLR ).
! PI(J) = MIN ( A(I,J) - U(I) - V(J) , SUCH THAT ROW  I  IS
!         LABELLED AND NOT EQUAL TO  FB(J) ).
! RC(J) = ROW PRECEDING COLUMN  J  IN THE CURRENT
!         ALTERNATING PATH.
! UC(L) = L-TH UNLABELLED COLUMN ( L=1,NUC ).
!
      implicit none
      INTEGER N 
      real*8  A(n,n),Z, U(n),V(N),PI(n), IA, MIN
      INTEGER F(N),LR(n),UC(n)
      INTEGER FB(n),RC(n)
      REAL*8, parameter :: INF = 1.d9
      integer  i,j,k,L, ii,jj, NUC , NLR, R
! INITIALIZATION.
      LR(1) = II
      DO 10 K=1,N
        PI(K) = A(II,K) - U(II) - V(K)
        RC(K) = II
        UC(K) = K
   10 CONTINUE
      NUC = N
      NLR = 1
      GO TO 40
! SCANNING OF THE LABELLED ROWS.
   20 R = LR(NLR)
      DO 30 L=1,NUC
        J = UC(L)
        IA = A(R,J) - U(R) - V(J)
        IF ( IA .GE. PI(J) ) GO TO 30
        PI(J) = IA
        RC(J) = R
   30 CONTINUE
! SEARCH FOR A ZERO ELEMENT IN AN UNLABELLED COLUMN.
   40 DO 50 L=1,NUC
        J = UC(L)
        IF ( PI(J) .EQ. 0.d0 ) GO TO 100
   50 CONTINUE
! UPDATING OF THE DUAL VARIABLES  U(I)  AND  V(J) .
      MIN = INF
      DO 60 L=1,NUC
        J = UC(L)
        IF ( MIN .GT. PI(J) ) MIN = PI(J)
   60 CONTINUE
      DO 70 L=1,NLR
        R = LR(L)
        U(R) = U(R) + MIN
   70 CONTINUE
      DO 90 J=1,N
        IF ( PI(J) .EQ. 0.d0 ) GO TO 80
        PI(J) = PI(J) - MIN
        GO TO 90
   80   V(J) = V(J) - MIN
   90 CONTINUE
      GO TO 40
  100 IF ( FB(J) .EQ. 0 ) GO TO 110
! LABELLING OF ROW  FB(J)  AND REMOVAL OF THE LABEL  OF COLUMN  J .
      NLR = NLR + 1
      LR(NLR) = FB(J)
      UC(L) = UC(NUC)
      NUC = NUC - 1
      GO TO 20
! DETERMINATION OF THE UNASSIGNED COLUMN  J .
  110 JJ = J
      RETURN
      END
 
!COSINE
!!!!!!MODIFIED METHOD TO INCLUDE COSINE DISTANCE
!!!!!!=====================================================================================================
!!!!!
!!!!!subroutine hungarian(nat,rat1,rat2,lfp,fp1,fp2,cost,xat1,xat2)
!!!!!implicit none
!!!!!integer  :: nat,iat,jat,lfp 
!!!!!  real(8),dimension(3,nat)         :: rat1,rat2,rat2copy       ! atomic position
!!!!!  real(8)                          :: cost,r2,fp1(lfp,nat), fp2(lfp,nat),fp2copy(lfp,nat)   
!!!!!  INTEGER f(nat)
!!!!!  REAL*8,allocatable :: A(:,:)
!!!!!  real*8  :: ddot
!!!!!character(2) xat1(nat),xat2(nat),xat2copy(nat)
!!!!!
!!!!!
!!!!!!print*,xat1; stop
!!!!!
!!!!!      allocate(A(nat,nat))
!!!!!      do jat=1,nat
!!!!!      do iat=1,nat
!!!!!         if(xat1(iat) .ne. xat2(jat))  then 
!!!!!             A(iat,jat)=1.d7   ! do not permute!
!!!!!         else
!!!!!!         A(iat,jat)=sum(rat1(:,iat)-rat2(:,jat))**2
!!!!!!         A(iat,jat)=sum((fp1(1:lfp,iat)-fp2(1:lfp,jat))**2)/dble(nat*lfp)
!!!!!!         A(iat,jat)=sum(((fp1(1:lfp,iat)-fp2(1:lfp,jat))/(fp1(1:lfp,iat)+fp2(1:lfp,jat)))**2)
!!!!!!         A(iat,jat)=1.d0 - (dot_product(fp1(:,iat),fp2(:,jat)))**2/ (dot_product(fp1(:,iat),fp1(:,iat)) * dot_product(fp2(:,jat),fp2(:,jat)))
!!!!!         A(iat,jat)=1.d0 - (ddot(lfp, fp1(:,iat),1, fp2(:,jat),1))**2/ (ddot(lfp, fp1(:,iat),1, fp1(:,iat),1)*ddot(lfp, fp2(:,jat),1, fp2(:,jat),1))
!!!!!        endif
!!!!!      enddo
!!!!!      enddo
!!!!!
!!!!!      CALL APC(Nat,A,F,cost)
!!!!!!m=nat;       call assndx(1,a, nat, m, F, cost)   ! this is slower by a factor of ~3
!!!!!
!!!!!      xat2copy =xat2 
!!!!!      rat2copy=rat2
!!!!!      fp2copy =fp2 
!!!!!
!!!!!      do iat=1,nat
!!!!!        xat2(iat)  =xat2copy(F(iat))
!!!!!        rat2(1:3,iat)  =rat2copy(1:3,F(iat))
!!!!!        fp2 (1:lfp,iat)=fp2copy (1:lfp,F(iat))
!!!!!      enddo
!!!!!
!!!!!!print '(i)',f(1:nat)
!!!!!end subroutine hungarian
!!!!!
!!!!!
!!!!!
!!!!!SUBROUTINE APC(N,A,F,Z)
!!!!!implicit none
!!!!!! Modified by Ali Sadeghi to get real*8 matrix A(N,N) and converted to F90 
!!!!!!
!!!!!! SOLUTION OF THE LINEAR MIN-SUM ASSIGNMENT PROBLEM.
!!!!!! HUNGARIAN METHOD. COMPLEXITY O(N**3).
!!!!!!
!!!!!! MEANING OF THE INPUT PARAMETERS:
!!!!!! N      = NUMBER OF ROWS AND COLUMNS OF THE COST MATRIX.
!!!!!! A(I,J) = COST OF THE ASSIGNMENT OF ROW  I  TO COLUMN  J .
!!!!!! ON RETURN, THE INPUT PARAMETERS ARE UNCHANGED.
!!!!!!
!!!!!! MEANING OF THE OUTPUT PARAMETERS:
!!!!!! F(I) = COLUMN ASSIGNED TO ROW  I .
!!!!!! Z    = COST OF THE OPTIMAL ASSIGNMENT =
!!!!!!      = A(1,F(1)) + A(2,F(2)) + ... + A(N,F(N)) .
!!!!!!
!!!!!!
!!!!!! THE CODE IS BASED ON THE HUNGARIAN METHOD AS DESCRIBED BY
!!!!!! LAWLER (COMBINATORIAL OPTIMIZATION : NETWORKS AND
!!!!!! MATROIDS, HOLT, RINEHART AND WINSTON, NEW YORK, 1976).
!!!!!! THE ALGORITHMIC PASCAL-LIKE DESCRIPTION OF THE CODE IS
!!!!!! GIVEN IN G.CARPANETO, S.MARTELLO AND P.TOTH, ALGORITHMS AND
!!!!!! CODES FOR THE ASSIGNMENT PROBLEM, ANNALS OF OPERATIONS
!!!!!! RESEARCH 7, 1988.
!!!!!!
!!!!!! SUBROUTINE APC DETERMINES THE INITIAL DUAL AND PARTIAL
!!!!!! PRIMAL SOLUTIONS AND THEN SEARCHES FOR AUGMENTING PATHS
!!!!!! UNTIL ALL ROWS AND COLUMNS ARE ASSIGNED.
!!!!!!
!!!!!! MEANING OF THE MAIN INTERNAL VARIABLES:
!!!!!! FB(J) = ROW ASSIGNED TO COLUMN  J .
!!!!!! M     = NUMBER OF INITIAL ASSIGNMENTS.
!!!!!! U(I)  = DUAL VARIABLE ASSOCIATED WITH ROW  I .
!!!!!! V(J)  = DUAL VARIABLE ASSOCIATED WITH COLUMN  J .
!!!!!!
!!!!!! APC NEEDS THE FOLLOWING SUBROUTINES: INCR
!!!!!!                                      INIT
!!!!!!                                      PATH
!!!!!!
!!!!!! THIS WORK WAS SUPPORTED BY  C.N.R. , ITALY.
!!!!!      INTEGER N
!!!!!      REAL*8  A(n,n),Z,U(n),V(n)
!!!!!      integer F(N),FB(n), RC(n)
!!!!!      INTEGER M,I,J,K
!!!!!! SEARCH FOR THE INITIAL DUAL AND PARTIAL PRIMAL SOLUTIONS.
!!!!!      CALL INIT(N,A,F,M,U,V,FB,RC)
!!!!!      IF ( M .NE. N ) then 
!!!!!! SOLUTION OF THE REDUCED PROBLEM.
!!!!!      DO  I=1,N
!!!!!        IF ( F(I) == 0 ) THEN 
!!!!!! DETERMINATION OF AN AUGMENTING PATH STARTING FROM ROW  I .
!!!!!        CALL PATH(N,A,I,F,J,U,V,FB,RC)
!!!!!! ASSIGNMENT OF ROW  I  AND COLUMN  J .
!!!!!        CALL INCR(n,F,J,FB,RC)
!!!!!        ENDIF
!!!!!      ENDDO    
!!!!!      ENDIF
!!!!!
!!!!!! COMPUTATION OF THE SOLUTION COST  Z .
!!!!!      Z = sum(u(1:N)) + sum(V(1:N))
!!!!!      END
!!!!!
!!!!!
!!!!!      SUBROUTINE INCR(n,F,J,FB,RC)
!!!!!!
!!!!!! ASSIGNMENT OF COLUMN  J .
!!!!!!
!!!!!      INTEGER n,I,J,JJ,  F(n),FB(n),RC(n)
!!!!!   10 I = RC(J)
!!!!!      FB(J) = I
!!!!!      JJ = F(I)
!!!!!      F(I) = J
!!!!!      J = JJ
!!!!!      IF ( J > 0 ) GO TO 10
!!!!!      RETURN
!!!!!      END
!!!!!
!!!!!
!!!!!      SUBROUTINE INIT(N,A,F,M,U,V,FB,P)
!!!!!!
!!!!!! SEARCH FOR THE INITIAL DUAL AND PARTIAL PRIMAL SOLUTIONS.
!!!!!!
!!!!!! P(I) = FIRST UNSCANNED COLUMN OF ROW  I .
!!!!!       IMPLICIT NONE
!!!!!!
!!!!!      INTEGER n,m, F(n),FB(n),P(n)
!!!!!      real*8 A(n,n) , U(n),V(n)
!!!!!      REAL*8, parameter :: INF = 1.d9
!!!!!      real*8 min, IA
!!!!!      integer i,j, k,R, JMIN, KK
!!!!!! PHASE 1 .
!!!!!      M = 0
!!!!!      F(1:N)=0
!!!!!      FB(1:N)=0
!!!!!! SCANNING OF THE COLUMNS ( INITIALIZATION OF  V(J) ).
!!!!!      DO 40 J=1,N
!!!!!        MIN = INF
!!!!!        DO 30 I=1,N
!!!!!          IA = A(I,J)
!!!!!          IF ( IA .GT. MIN ) GO TO 30
!!!!!          IF ( IA .LT. MIN ) GO TO 20
!!!!!          IF ( F(I) .NE. 0 ) GO TO 30
!!!!!   20     MIN = IA
!!!!!          R = I
!!!!!   30   CONTINUE
!!!!!        V(J) = MIN
!!!!!        IF ( F(R) .NE. 0 ) GO TO 40
!!!!!! ASSIGNMENT OF COLUMN  J  TO ROW  R .
!!!!!        M = M + 1
!!!!!        FB(J) = R
!!!!!        F(R) = J
!!!!!        U(R) = 0.d0
!!!!!        P(R) = J + 1
!!!!!   40 CONTINUE
!!!!!! PHASE 2 .
!!!!!! SCANNING OF THE UNASSIGNED ROWS ( UPDATING OF  U(I) ).
!!!!!      DO 110 I=1,N
!!!!!        IF ( F(I) .NE. 0 ) GO TO 110
!!!!!        MIN = INF
!!!!!        DO 60 K=1,N
!!!!!          IA = A(I,K) - V(K)
!!!!!          IF ( IA .GT. MIN )  GO TO 60
!!!!!          IF ( IA .LT. MIN )  GO TO 50
!!!!!          IF ( FB(K) .NE. 0 ) GO TO 60
!!!!!          IF ( FB(J) .EQ. 0 ) GO TO 60
!!!!!   50     MIN = IA
!!!!!          J = K
!!!!!   60   CONTINUE
!!!!!        U(I) = MIN
!!!!!        JMIN = J
!!!!!        IF ( FB(J) .EQ. 0 ) GO TO 100
!!!!!        DO 80 J=JMIN,N
!!!!!          IF ( A(I,J) - V(J) .GT. MIN ) GO TO 80
!!!!!          R = FB(J)
!!!!!          KK = P(R)
!!!!!          IF ( KK .GT. N ) GO TO 80
!!!!!          DO 70 K=KK,N
!!!!!            IF ( FB(K) .GT. 0 ) GO TO 70
!!!!!            IF ( A(R,K) - U(R) - V(K) .EQ. 0.d0 ) GO TO 90
!!!!!   70     CONTINUE
!!!!!          P(R) = N + 1
!!!!!   80   CONTINUE
!!!!!        GO TO 110
!!!!!! REASSIGNMENT OF ROW  R  AND COLUMN  K .
!!!!!   90   F(R) = K
!!!!!        FB(K) = R
!!!!!        P(R) = K + 1
!!!!!! ASSIGNMENT OF COLUMN  J  TO ROW  I .
!!!!!  100   M = M + 1
!!!!!        F(I) = J
!!!!!        FB(J)= I
!!!!!        P(I) = J + 1
!!!!!  110 CONTINUE
!!!!!      RETURN
!!!!!      END
!!!!!      SUBROUTINE PATH(N,A,II,F,JJ,U,V,FB,RC)
!!!!!!
!!!!!! DETERMINATION OF AN AUGMENTING PATH STARTING FROM
!!!!!! UNASSIGNED ROW  II  AND TERMINATING AT UNASSIGNED COLUMN
!!!!!! JJ , WITH UPDATING OF DUAL VARIABLES  U(I)  AND  V(J) .
!!!!!!
!!!!!! MEANING OF THE MAIN INTERNAL VARIABLES:
!!!!!! LR(L) = L-TH LABELLED ROW ( L=1,NLR ).
!!!!!! PI(J) = MIN ( A(I,J) - U(I) - V(J) , SUCH THAT ROW  I  IS
!!!!!!         LABELLED AND NOT EQUAL TO  FB(J) ).
!!!!!! RC(J) = ROW PRECEDING COLUMN  J  IN THE CURRENT
!!!!!!         ALTERNATING PATH.
!!!!!! UC(L) = L-TH UNLABELLED COLUMN ( L=1,NUC ).
!!!!!!
!!!!!      implicit none
!!!!!      INTEGER N 
!!!!!      real*8  A(n,n),Z, U(n),V(N),PI(n), IA, MIN
!!!!!      INTEGER F(N),LR(n),UC(n)
!!!!!      INTEGER FB(n),RC(n)
!!!!!      REAL*8, parameter :: INF = 1.d9
!!!!!      integer  i,j,k,L, ii,jj, NUC , NLR, R
!!!!!! INITIALIZATION.
!!!!!      LR(1) = II
!!!!!      DO 10 K=1,N
!!!!!        PI(K) = A(II,K) - U(II) - V(K)
!!!!!        RC(K) = II
!!!!!        UC(K) = K
!!!!!   10 CONTINUE
!!!!!      NUC = N
!!!!!      NLR = 1
!!!!!      GO TO 40
!!!!!! SCANNING OF THE LABELLED ROWS.
!!!!!   20 R = LR(NLR)
!!!!!      DO 30 L=1,NUC
!!!!!        J = UC(L)
!!!!!        IA = A(R,J) - U(R) - V(J)
!!!!!        IF ( IA .GE. PI(J) ) GO TO 30
!!!!!        PI(J) = IA
!!!!!        RC(J) = R
!!!!!   30 CONTINUE
!!!!!! SEARCH FOR A ZERO ELEMENT IN AN UNLABELLED COLUMN.
!!!!!   40 DO 50 L=1,NUC
!!!!!        J = UC(L)
!!!!!        IF ( PI(J) .EQ. 0.d0 ) GO TO 100
!!!!!   50 CONTINUE
!!!!!! UPDATING OF THE DUAL VARIABLES  U(I)  AND  V(J) .
!!!!!      MIN = INF
!!!!!      DO 60 L=1,NUC
!!!!!        J = UC(L)
!!!!!        IF ( MIN .GT. PI(J) ) MIN = PI(J)
!!!!!   60 CONTINUE
!!!!!      DO 70 L=1,NLR
!!!!!        R = LR(L)
!!!!!        U(R) = U(R) + MIN
!!!!!   70 CONTINUE
!!!!!      DO 90 J=1,N
!!!!!        IF ( PI(J) .EQ. 0.d0 ) GO TO 80
!!!!!        PI(J) = PI(J) - MIN
!!!!!        GO TO 90
!!!!!   80   V(J) = V(J) - MIN
!!!!!   90 CONTINUE
!!!!!      GO TO 40
!!!!!  100 IF ( FB(J) .EQ. 0 ) GO TO 110
!!!!!! LABELLING OF ROW  FB(J)  AND REMOVAL OF THE LABEL  OF COLUMN  J .
!!!!!      NLR = NLR + 1
!!!!!      LR(NLR) = FB(J)
!!!!!      UC(L) = UC(NUC)
!!!!!      NUC = NUC - 1
!!!!!      GO TO 20
!!!!!! DETERMINATION OF THE UNASSIGNED COLUMN  J .
!!!!!  110 JJ = J
!!!!!      RETURN
!!!!!      END
!!!!! 
!!!!!
!!!!!
!!!!!
!ORIGINAL WITH EUCLIDIAN DISTANCE
!!!!!=====================================================================================================
!!!!
!!!!subroutine hungarian(nat,rat1,rat2,lfp,fp1,fp2,cost,xat1,xat2)
!!!!implicit none
!!!!integer  :: nat,iat,jat,lfp 
!!!!  real(8),dimension(3,nat)         :: rat1,rat2,rat2copy       ! atomic position
!!!!  real(8)                          :: cost,r2,fp1(lfp,nat), fp2(lfp,nat),fp2copy(lfp,nat)   
!!!!  INTEGER f(nat)
!!!!  REAL*8,allocatable :: A(:,:)
!!!!character(2) xat1(nat),xat2(nat),xat2copy(nat)
!!!!
!!!!
!!!!!print*,xat1; stop
!!!!
!!!!      allocate(A(nat,nat))
!!!!      do jat=1,nat
!!!!      do iat=1,nat
!!!!!         A(iat,jat)=sum(rat1(:,iat)-rat2(:,jat))**2
!!!!!         A(iat,jat)=sum((fp1(1:lfp,iat)-fp2(1:lfp,jat))**2)/dble(nat*lfp)
!!!!         A(iat,jat)=sum(((fp1(1:lfp,iat)-fp2(1:lfp,jat))/(fp1(1:lfp,iat)+fp2(1:lfp,jat)))**2)
!!!!         if(xat1(iat) .ne. xat2(jat)) A(iat,jat)=A(iat,jat)+1.d7
!!!!      enddo
!!!!      enddo
!!!!
!!!!
!!!!      CALL APC(Nat,A,F,cost)
!!!!!m=nat;       call assndx(1,a, nat, m, F, cost)   ! this is slower by a factor of ~3
!!!!
!!!!      xat2copy =xat2 
!!!!      rat2copy=rat2
!!!!      fp2copy =fp2 
!!!!
!!!!      do iat=1,nat
!!!!        xat2(iat)  =xat2copy(F(iat))
!!!!        rat2(1:3,iat)  =rat2copy(1:3,F(iat))
!!!!        fp2 (1:lfp,iat)=fp2copy (1:lfp,F(iat))
!!!!      enddo
!!!!
!!!!!print '(i)',f(1:nat)
!!!!end subroutine hungarian
!!!!
!!!!
!!!!
!!!!SUBROUTINE APC(N,A,F,Z)
!!!!implicit none
!!!!! Modified by Ali Sadeghi to get real*8 matrix A(N,N) and converted to F90 
!!!!!
!!!!! SOLUTION OF THE LINEAR MIN-SUM ASSIGNMENT PROBLEM.
!!!!! HUNGARIAN METHOD. COMPLEXITY O(N**3).
!!!!!
!!!!! MEANING OF THE INPUT PARAMETERS:
!!!!! N      = NUMBER OF ROWS AND COLUMNS OF THE COST MATRIX.
!!!!! A(I,J) = COST OF THE ASSIGNMENT OF ROW  I  TO COLUMN  J .
!!!!! ON RETURN, THE INPUT PARAMETERS ARE UNCHANGED.
!!!!!
!!!!! MEANING OF THE OUTPUT PARAMETERS:
!!!!! F(I) = COLUMN ASSIGNED TO ROW  I .
!!!!! Z    = COST OF THE OPTIMAL ASSIGNMENT =
!!!!!      = A(1,F(1)) + A(2,F(2)) + ... + A(N,F(N)) .
!!!!!
!!!!!
!!!!! THE CODE IS BASED ON THE HUNGARIAN METHOD AS DESCRIBED BY
!!!!! LAWLER (COMBINATORIAL OPTIMIZATION : NETWORKS AND
!!!!! MATROIDS, HOLT, RINEHART AND WINSTON, NEW YORK, 1976).
!!!!! THE ALGORITHMIC PASCAL-LIKE DESCRIPTION OF THE CODE IS
!!!!! GIVEN IN G.CARPANETO, S.MARTELLO AND P.TOTH, ALGORITHMS AND
!!!!! CODES FOR THE ASSIGNMENT PROBLEM, ANNALS OF OPERATIONS
!!!!! RESEARCH 7, 1988.
!!!!!
!!!!! SUBROUTINE APC DETERMINES THE INITIAL DUAL AND PARTIAL
!!!!! PRIMAL SOLUTIONS AND THEN SEARCHES FOR AUGMENTING PATHS
!!!!! UNTIL ALL ROWS AND COLUMNS ARE ASSIGNED.
!!!!!
!!!!! MEANING OF THE MAIN INTERNAL VARIABLES:
!!!!! FB(J) = ROW ASSIGNED TO COLUMN  J .
!!!!! M     = NUMBER OF INITIAL ASSIGNMENTS.
!!!!! U(I)  = DUAL VARIABLE ASSOCIATED WITH ROW  I .
!!!!! V(J)  = DUAL VARIABLE ASSOCIATED WITH COLUMN  J .
!!!!!
!!!!! APC NEEDS THE FOLLOWING SUBROUTINES: INCR
!!!!!                                      INIT
!!!!!                                      PATH
!!!!!
!!!!! THIS WORK WAS SUPPORTED BY  C.N.R. , ITALY.
!!!!      INTEGER N
!!!!      REAL*8  A(n,n),Z,U(n),V(n)
!!!!      integer F(N),FB(n), RC(n)
!!!!      INTEGER M,I,J,K
!!!!! SEARCH FOR THE INITIAL DUAL AND PARTIAL PRIMAL SOLUTIONS.
!!!!      CALL INIT(N,A,F,M,U,V,FB,RC)
!!!!      IF ( M .NE. N ) then 
!!!!! SOLUTION OF THE REDUCED PROBLEM.
!!!!      DO  I=1,N
!!!!        IF ( F(I) == 0 ) THEN 
!!!!! DETERMINATION OF AN AUGMENTING PATH STARTING FROM ROW  I .
!!!!        CALL PATH(N,A,I,F,J,U,V,FB,RC)
!!!!! ASSIGNMENT OF ROW  I  AND COLUMN  J .
!!!!        CALL INCR(n,F,J,FB,RC)
!!!!        ENDIF
!!!!      ENDDO    
!!!!      ENDIF
!!!!
!!!!! COMPUTATION OF THE SOLUTION COST  Z .
!!!!      Z = sum(u(1:N)) + sum(V(1:N))
!!!!      END
!!!!
!!!!
!!!!      SUBROUTINE INCR(n,F,J,FB,RC)
!!!!!
!!!!! ASSIGNMENT OF COLUMN  J .
!!!!!
!!!!      INTEGER n,I,J,JJ,  F(n),FB(n),RC(n)
!!!!   10 I = RC(J)
!!!!      FB(J) = I
!!!!      JJ = F(I)
!!!!      F(I) = J
!!!!      J = JJ
!!!!      IF ( J > 0 ) GO TO 10
!!!!      RETURN
!!!!      END
!!!!
!!!!
!!!!      SUBROUTINE INIT(N,A,F,M,U,V,FB,P)
!!!!!
!!!!! SEARCH FOR THE INITIAL DUAL AND PARTIAL PRIMAL SOLUTIONS.
!!!!!
!!!!! P(I) = FIRST UNSCANNED COLUMN OF ROW  I .
!!!!       IMPLICIT NONE
!!!!!
!!!!      INTEGER n,m, F(n),FB(n),P(n)
!!!!      real*8 A(n,n) , U(n),V(n)
!!!!      REAL*8, parameter :: INF = 1.d9
!!!!      real*8 min, IA
!!!!      integer i,j, k,R, JMIN, KK
!!!!! PHASE 1 .
!!!!      M = 0
!!!!      F(1:N)=0
!!!!      FB(1:N)=0
!!!!! SCANNING OF THE COLUMNS ( INITIALIZATION OF  V(J) ).
!!!!      DO 40 J=1,N
!!!!        MIN = INF
!!!!        DO 30 I=1,N
!!!!          IA = A(I,J)
!!!!          IF ( IA .GT. MIN ) GO TO 30
!!!!          IF ( IA .LT. MIN ) GO TO 20
!!!!          IF ( F(I) .NE. 0 ) GO TO 30
!!!!   20     MIN = IA
!!!!          R = I
!!!!   30   CONTINUE
!!!!        V(J) = MIN
!!!!        IF ( F(R) .NE. 0 ) GO TO 40
!!!!! ASSIGNMENT OF COLUMN  J  TO ROW  R .
!!!!        M = M + 1
!!!!        FB(J) = R
!!!!        F(R) = J
!!!!        U(R) = 0.d0
!!!!        P(R) = J + 1
!!!!   40 CONTINUE
!!!!! PHASE 2 .
!!!!! SCANNING OF THE UNASSIGNED ROWS ( UPDATING OF  U(I) ).
!!!!      DO 110 I=1,N
!!!!        IF ( F(I) .NE. 0 ) GO TO 110
!!!!        MIN = INF
!!!!        DO 60 K=1,N
!!!!          IA = A(I,K) - V(K)
!!!!          IF ( IA .GT. MIN )  GO TO 60
!!!!          IF ( IA .LT. MIN )  GO TO 50
!!!!          IF ( FB(K) .NE. 0 ) GO TO 60
!!!!          IF ( FB(J) .EQ. 0 ) GO TO 60
!!!!   50     MIN = IA
!!!!          J = K
!!!!   60   CONTINUE
!!!!        U(I) = MIN
!!!!        JMIN = J
!!!!        IF ( FB(J) .EQ. 0 ) GO TO 100
!!!!        DO 80 J=JMIN,N
!!!!          IF ( A(I,J) - V(J) .GT. MIN ) GO TO 80
!!!!          R = FB(J)
!!!!          KK = P(R)
!!!!          IF ( KK .GT. N ) GO TO 80
!!!!          DO 70 K=KK,N
!!!!            IF ( FB(K) .GT. 0 ) GO TO 70
!!!!            IF ( A(R,K) - U(R) - V(K) .EQ. 0.d0 ) GO TO 90
!!!!   70     CONTINUE
!!!!          P(R) = N + 1
!!!!   80   CONTINUE
!!!!        GO TO 110
!!!!! REASSIGNMENT OF ROW  R  AND COLUMN  K .
!!!!   90   F(R) = K
!!!!        FB(K) = R
!!!!        P(R) = K + 1
!!!!! ASSIGNMENT OF COLUMN  J  TO ROW  I .
!!!!  100   M = M + 1
!!!!        F(I) = J
!!!!        FB(J)= I
!!!!        P(I) = J + 1
!!!!  110 CONTINUE
!!!!      RETURN
!!!!      END
!!!!      SUBROUTINE PATH(N,A,II,F,JJ,U,V,FB,RC)
!!!!!
!!!!! DETERMINATION OF AN AUGMENTING PATH STARTING FROM
!!!!! UNASSIGNED ROW  II  AND TERMINATING AT UNASSIGNED COLUMN
!!!!! JJ , WITH UPDATING OF DUAL VARIABLES  U(I)  AND  V(J) .
!!!!!
!!!!! MEANING OF THE MAIN INTERNAL VARIABLES:
!!!!! LR(L) = L-TH LABELLED ROW ( L=1,NLR ).
!!!!! PI(J) = MIN ( A(I,J) - U(I) - V(J) , SUCH THAT ROW  I  IS
!!!!!         LABELLED AND NOT EQUAL TO  FB(J) ).
!!!!! RC(J) = ROW PRECEDING COLUMN  J  IN THE CURRENT
!!!!!         ALTERNATING PATH.
!!!!! UC(L) = L-TH UNLABELLED COLUMN ( L=1,NUC ).
!!!!!
!!!!      implicit none
!!!!      INTEGER N 
!!!!      real*8  A(n,n),Z, U(n),V(N),PI(n), IA, MIN
!!!!      INTEGER F(N),LR(n),UC(n)
!!!!      INTEGER FB(n),RC(n)
!!!!      REAL*8, parameter :: INF = 1.d9
!!!!      integer  i,j,k,L, ii,jj, NUC , NLR, R
!!!!! INITIALIZATION.
!!!!      LR(1) = II
!!!!      DO 10 K=1,N
!!!!        PI(K) = A(II,K) - U(II) - V(K)
!!!!        RC(K) = II
!!!!        UC(K) = K
!!!!   10 CONTINUE
!!!!      NUC = N
!!!!      NLR = 1
!!!!      GO TO 40
!!!!! SCANNING OF THE LABELLED ROWS.
!!!!   20 R = LR(NLR)
!!!!      DO 30 L=1,NUC
!!!!        J = UC(L)
!!!!        IA = A(R,J) - U(R) - V(J)
!!!!        IF ( IA .GE. PI(J) ) GO TO 30
!!!!        PI(J) = IA
!!!!        RC(J) = R
!!!!   30 CONTINUE
!!!!! SEARCH FOR A ZERO ELEMENT IN AN UNLABELLED COLUMN.
!!!!   40 DO 50 L=1,NUC
!!!!        J = UC(L)
!!!!        IF ( PI(J) .EQ. 0.d0 ) GO TO 100
!!!!   50 CONTINUE
!!!!! UPDATING OF THE DUAL VARIABLES  U(I)  AND  V(J) .
!!!!      MIN = INF
!!!!      DO 60 L=1,NUC
!!!!        J = UC(L)
!!!!        IF ( MIN .GT. PI(J) ) MIN = PI(J)
!!!!   60 CONTINUE
!!!!      DO 70 L=1,NLR
!!!!        R = LR(L)
!!!!        U(R) = U(R) + MIN
!!!!   70 CONTINUE
!!!!      DO 90 J=1,N
!!!!        IF ( PI(J) .EQ. 0.d0 ) GO TO 80
!!!!        PI(J) = PI(J) - MIN
!!!!        GO TO 90
!!!!   80   V(J) = V(J) - MIN
!!!!   90 CONTINUE
!!!!      GO TO 40
!!!!  100 IF ( FB(J) .EQ. 0 ) GO TO 110
!!!!! LABELLING OF ROW  FB(J)  AND REMOVAL OF THE LABEL  OF COLUMN  J .
!!!!      NLR = NLR + 1
!!!!      LR(NLR) = FB(J)
!!!!      UC(L) = UC(NUC)
!!!!      NUC = NUC - 1
!!!!      GO TO 20
!!!!! DETERMINATION OF THE UNASSIGNED COLUMN  J .
!!!!  110 JJ = J
!!!!      RETURN
!!!!      END
!!!! 
