subroutine get_fp_calypso(nat,rxyz,latvec,r_cut,kinds,nkinds,fp_dim,nl,fp)
!This routine will compute the Bond Characterization Matrix as proposed by 
!Yanchao Wang in computer physics communications 183, 2063
!The dimension of fp is a matrix with dimension of AB interactions and nl spherical harmonics (be careful: if nl=10 it means lmax=(nl-1)*2
!fp_dim=ntypat*(ntypat+1)/2
implicit none
integer:: nl !Number of l components, here only even ones 
integer:: fp_dim !Number of AB interactions, doublecounting eliminated
integer:: nat,nkinds
integer:: nbond(fp_dim),kinds(nat)
integer:: lmax,yll,ymm,yl,ym
integer:: imax,imin
integer:: nec1,nec2,nec3,iat,jat,k,l,m,iarr,iat_kind,jat_kind,i,j
integer:: llist(nl)
real(8):: latvec(3,3),rxyz(3,nat),fp(fp_dim,nl)
real(8), allocatable:: rxyzexp(:,:,:,:,:),transvecall(:,:,:,:)
real(8):: sigma,r_cut(fp_dim) !Cutoff for each AB interaction
real(8), parameter :: pi=3.141592653589793238462643383279502884197d0
real(8):: min_bond(fp_dim),alpha(fp_dim),ylm_r,ylm_i
real(8):: rij, rij_vec(3),exp_fac,rxyzj(3),theta,phi,eps,tr,ti
real(8),allocatable:: qlm(:,:,:,:)
real(8):: double_counting
!Initiallize fp
fp=0.d0

!Maximum l-component
if(nl.gt.10) stop "Not more than 10 nl implemented"

!Radial decay function rate
eps=1.d-4 

!We can set a list of l components, with nl elements
!!!All l till nl
!!do i=1,nl
!!  llist(i)=i-1
!!enddo
!Only even l
do i=1,nl
  llist(i)=(i-1)*2
enddo

!Allocate qml
lmax=maxval(llist)
allocate(qlm(2,nl,2*lmax+1,fp_dim))
qlm=0.d0

!Create the expanded unit cells
call n_rep_dim(latvec,maxval(r_cut),nec1,nec2,nec3)
write(*,'(a,3(i3,1x),a)') " Creating expansion for periodic BCM FP ",nec1,nec2,nec3,"..."
allocate(rxyzexp(3,nat,2*nec1+1,2*nec2+1,2*nec3+1),transvecall(3,2*nec1+1,2*nec2+1,2*nec3+1))
call expand_dim(rxyz,rxyzexp,transvecall,latvec,nat,nec1,nec2,nec3)


!First find min_bond, which is b_AB in the paper
min_bond=1.d10
do iat=1,nat
  do jat=1,iat
     do k=1,2*nec1+1
     do l=1,2*nec2+1
     do m=1,2*nec3+1
      if(iat==jat.and.k==nec1+1.and.l==nec2+1.and.m==nec3+1) goto 2010
        rxyzj(:)=rxyz(:,jat)+transvecall(:,k,l,m)
        rij=(rxyz(1,iat)-rxyzj(1))**2+(rxyz(2,iat)-rxyzj(2))**2+(rxyz(3,iat)-rxyzj(3))**2
        imax=max(Kinds(iat),Kinds(jat))
        imin=min(Kinds(iat),Kinds(jat))
        iarr=imax*(imax+1)
        iarr=(iarr/2+imin-imax)
        if(min_bond(iarr).gt.rij) min_bond(iarr)=rij
      2010 continue
     enddo
     enddo
     enddo
  enddo
enddo
min_bond=sqrt(min_bond)

!Determine alpha
do i=1,fp_dim
  alpha(i)=-log(eps)/(r_cut(i)-min_bond(i))
enddo
!!alpha=0.8d0



!Get the Ylm and Q_lm^dab
nbond=0
fp=0.d0
do iat=1,nat
  do jat=1,nat
     do k=1,2*nec1+1
     do l=1,2*nec2+1
     do m=1,2*nec3+1
      if(iat==jat.and.k==nec1+1.and.l==nec2+1.and.m==nec3+1) goto 2020
        rxyzj(:)=rxyz(:,jat)+transvecall(:,k,l,m)
        rij_vec=rxyzj(:)-rxyz(:,iat)
        call cart2polar(rij_vec,rij,theta,phi)
        imax=max(Kinds(iat),Kinds(jat))
        imin=min(Kinds(iat),Kinds(jat))
        iarr=imax*(imax+1)
        iarr=(iarr/2+imin-imax)
        if(rij.lt.r_cut(iarr)) then
!!           double_counting=1.d0/min(real(imax-imin,8)+1.d0,2.d0)
           double_counting=1.d0
!!           write(*,'(a,i5,i5,f5.2)') "DOUBLE",imax,imin,double_counting
           nbond(iarr)=nbond(iarr)+1
!!                exp_fac=exp(-alpha(iarr)*(rij-min_bond(iarr)))*double_counting
                exp_fac=1.d0*double_counting
!!                exp_fac=1.d0/(rij-min_bond(iarr)+1.d0)**4*double_counting
!!                  call poly1(rij**2,min_bond(iarr),r_cut(iarr),exp_fac)
                do yll=1,nl
                   yl=llist(yll)
                   ym=-yl 
                   do ymm=1,2*yl+1
!!                      call ylm_ri (yl, ym, theta, phi, ylm_r, ylm_i)
                      call ylm_mathematica (yl, ym, theta, phi, ylm_r, ylm_i)
!!!                      write(*,'(a,3f25.16)') "xcart",rij_vec
!!!                      write(*,'(a,3f25.16)') "r,theta, phi",rij,theta,phi
!!!                      write(*,*) "l,m",yl,ym
!!!                      write(*,*) "yml",ylm_r,ylm_i
!                      if(modulo(yl,2)==0) then
                        qlm(1,yll,ymm,iarr)=qlm(1,yll,ymm,iarr)+ylm_r*exp_fac
                        qlm(2,yll,ymm,iarr)=qlm(2,yll,ymm,iarr)+ylm_i*exp_fac
!                      else
!                        qlm(1,yll,ymm,iarr)=qlm(1,yll,ymm,iarr)+abs(ylm_r)*exp_fac
!                        qlm(2,yll,ymm,iarr)=qlm(2,yll,ymm,iarr)+abs(ylm_i)*exp_fac
!                      endif
                      ym=ym+1
                   enddo
                enddo
        endif
      2020 continue
     enddo
     enddo
     enddo
  enddo
enddo


!Normalize
do i=1,fp_dim
 if(nbond(i).gt.0) qlm(:,:,:,i)=qlm(:,:,:,i)/real(nbond(i),8)
enddo
!!write(*,*) "qlm"
!!iarr=1
!!do l=1,nl
!!write(*,'(10(es15.7))') (qlm(:,l,m,iarr), m=1,2*lmax+1)
!!enddo

!Sum over all m components and get fp
do i=1,fp_dim
  do yll=1,nl
     yl=llist(yll)
     do ymm=1,2*yl+1
        fp(i,yll)=fp(i,yll)+qlm(1,yll,ymm,iarr)**2+qlm(2,yll,ymm,iarr)**2
     enddo
     fp(i,yll)=fp(i,yll)*4.d0*pi/(2.d0*real(yl,8)+1.d0)
     fp(i,yll)=sqrt(fp(i,yll))
  enddo
enddo
!!!write(*,*) "fp",fp
!!!write(*,*) "nbond",nbond
!!!write(*,*) "minbond",min_bond
!Set the zero moment to zero
!fp(:,1)=0.d0


deallocate(qlm,rxyzexp,transvecall)
end subroutine

subroutine get_distance_calypso(parini,fp1,fp2,fp_dim,nl,dist)
use mod_parini, only: typ_parini
implicit none
type(typ_parini), intent(in):: parini
integer:: fp_dim,nl,yll,i,j
real(8):: fp1(fp_dim*nl),fp2(fp_dim*nl),dist
!This sunroutine will compute the cosine distance between two fingerprints fp1 and fp2 and 
!store the output distance.  All conventions and 
!methods are from J.Chem.Phys, 130, 104504 (2009) and IEEE Symposium, Okt 21-23. (2008) (M.Valle and A.Oganov)
integer :: nkinds_sum(parini%ntypat_global), fp_size
real(8) :: fp1_1(nl*fp_dim),fp2_1(nl*fp_dim)
real(8) :: distance
integer             :: i_kind,j_kind,i_fp,k,imin,imax,iat
real(8)             :: w_ab(fp_dim),w_norm,num,denom,tmp_1,tmp_2

call res(fp1,fp1_1,fp_dim,nl)
call res(fp2,fp2_1,fp_dim,nl)

fp_size=nl

dist=0.d0
!Euclidic distance
do i=1,fp_dim*nl
    dist=dist+(fp1(i)-fp2(i))**2
enddo
dist=sqrt(dist) 
!!!Cosine distance
!!dist=abs(0.5d0*(1.d0-dot_product(fp1(:),fp2(:))/(sqrt(dot_product(fp1(:),fp1(:)))*sqrt(dot_product(fp2(:),fp2(:))))))
!!!!!!end subroutine
!subroutine get_cosinedistance(fp1,fp2,fp_size,fp_dim,nkinds,nkinds_sum,distance)
!This sunroutine will compute the cosine distance between two fingerprints fp1 and fp2 and 
!store the output distance.  All conventions and 
!methods are from J.Chem.Phys, 130, 104504 (2009) and IEEE Symposium, Okt 21-23. (2008) (M.Valle and A.Oganov)


!!!nkinds_sum=0
!!!do iat=1,nat
!!!  do i=1,ntypat
!!!  if(typat(iat)==i) nkinds_sum(i)= nkinds_sum(i)+1
!!!  enddo
!!!enddo
!!!!Compute the weight of each fingerprint
!!!w_ab=0.d0
!!!w_norm=0.d0
!!!do k=1,fp_dim
!!!   i=ceiling(-0.5d0+0.5d0*sqrt(1.d0+8.d0*real(k,8)))
!!!   j=(i*(1-i))/2+k
!!!   w_ab(k)=real(nkinds_sum(i),8)*real(nkinds_sum(j),8)
!!!   w_norm=w_norm+real(nkinds_sum(i),8)*real(nkinds_sum(j),8)
!!!enddo
!!!call scaleab(fp1_1,nl,fp_dim,w_ab)
!!!call scaleab(fp2_1,nl,fp_dim,w_ab)
!!!!!do i=1,nkinds
!!!!!  do j=1,i!nkinds
!!!!!  w_ab(i,j)=real(nkinds_sum(i),8)*real(nkinds_sum(j),8)
!!!!!  w_norm=w_norm+real(nkinds_sum(i),8)*real(nkinds_sum(j),8)
!!!!!  enddo
!!!!!enddo
!!!w_ab=w_ab/w_norm
!!!!write(*,*) w_ab
!!!
!!!!Compute the distance
!!!num=0.d0;  denom=0.d0;  tmp_1=0.d0;  tmp_2=0.d0
!!!
!!!if(trim(code).ne."lenosky_tb_lj") then
!!!!Usual case: atoms in general
!!!  do k=1,fp_dim
!!!       imin=(k-1)*fp_size+1
!!!       imax=k*fp_size
!!!       num=num+dot_product(fp1_1(imin:imax),fp2_1(imin:imax))*w_ab(k)**2
!!!       tmp_1=tmp_1+dot_product(fp1_1(imin:imax),fp1_1(imin:imax))*w_ab(k)**2
!!!       tmp_2=tmp_2+dot_product(fp2_1(imin:imax),fp2_1(imin:imax))*w_ab(k)**2
!!!  enddo
!!!else
!!!!Only for lenosky_tb_lj
!!!  do k=1,fp_dim
!!!!Stephan's formula to identify the LJ indexes...
!!!       i=ceiling(-0.5d0+0.5d0*sqrt(1.d0+8.d0*real(k,8)))
!!!       j=(i*(1-i))/2+k
!!!       if(.not.(znucl(i).ge.200.or.znucl(j).ge.200)) then 
!!!         imin=(k-1)*fp_size+1
!!!         imax=k*fp_size
!!!         num=num+dot_product(fp1_1(imin:imax),fp2_1(imin:imax))*w_ab(k)**2
!!!         tmp_1=tmp_1+dot_product(fp1_1(imin:imax),fp1_1(imin:imax))*w_ab(k)**2
!!!         tmp_2=tmp_2+dot_product(fp2_1(imin:imax),fp2_1(imin:imax))*w_ab(k)**2
!!!       endif
!!!  enddo
!!!endif
!!!
!!!
!!!denom=sqrt(tmp_1*tmp_2)
!!!dist=0.5d0*(1.d0-num/denom)
contains

subroutine res(a,b,n,m)
implicit none
real(8):: a(n,m),b(m,n)
integer:: n,m,i,j
do i=1,n
  do j=1,m
  b(j,i)=a(i,j)
  enddo
enddo
end subroutine
subroutine scaleab(fp,n,m,w_ab)
implicit none
integer:: n,m,i
real(8):: fp(n,m),w_ab(m)
do i=1,m
fp(:,i)=fp(:,i)/w_ab(i)
enddo
end subroutine

end subroutine




!!subroutine n_rep(latvec,cut,nec)
!!!This subroutine will return how many periodic expansions are necessary for the periodic boundary conditions
!!!with for the given cut. 
!!implicit none
!!real*8 :: latvec(3,3),cut,nvec(3,3),point(3),point0(3),dist(3),eps,dd
!!integer:: i
!!integer:: nec
!!! eps=1.d-6
!!nec=0
!!call nveclatvec(latvec,nvec)
!!point0=(/0.d0,0.d0,0.d0/)
!!do i=1,3
!!call dist2plane(latvec(:,mod(i+1,3)+1),nvec(:,i),point0,dist(i))
!!! write(*,*) "cut",i,cut, dist
!!enddo
!!dd=minval(abs(dist))
!!nec=int(cut/dd)+1
!!end subroutine
!!
!!
!!subroutine n_rep_dim(latvec,cut,nec1,nec2,nec3)
!!!This subroutine will return how many periodic expansions for each lattice vector direction are necessary for the periodic boundary conditions
!!!with for the given cut. nec1,nec2,nec3 for latvec(:,1),latvec(:,2),latvec(:,3)
!!implicit none
!!real*8 :: latvec(3,3),cut,nvec(3,3),point(3),point0(3),dist(3),eps,dd
!!integer:: i
!!integer:: nec1,nec2,nec3
!!! eps=1.d-6
!!nec1=0
!!nec2=0
!!nec3=0
!!call nveclatvec(latvec,nvec)
!!point0=(/0.d0,0.d0,0.d0/)
!!do i=1,3
!!call dist2plane(latvec(:,mod(i+1,3)+1),nvec(:,i),point0,dist(i))
!!! write(*,*) "cut",i,cut, dist
!!enddo
!!dist=abs(dist)
!!nec1=int(cut/dist(2))+1
!!nec2=int(cut/dist(3))+1
!!nec3=int(cut/dist(1))+1
!!end subroutine
!!
!!subroutine expand_dim(rxyz,rxyzout,transvecall,latvec,nat,nec1,nec2,nec3)
!!!This subroutine will expand the unit cell into the necessary periodic cells and store them in rxyzout. Only the first octant is computed
!!implicit none
!!real*8, intent(in)  :: rxyz(3,nat),latvec(3,3)
!!integer, intent(in) :: nat,nec1,nec2,nec3
!!real*8, intent(out) :: rxyzout(3,nat,2*nec1+1,2*nec2+1,2*nec3+1) !only necessary periodic images in the first octant plus the main cell
!!integer             :: iat,m,k,l
!!real*8,intent(inout):: transvecall(3,2*nec1+1,2*nec2+1,2*nec3+1)
!!do m=-nec3,nec3
!!   do k=-nec2,nec2
!!      do l=-nec1,nec1
!!      transvecall(:,l+nec1+1,k+nec2+1,m+nec3+1)=real(l,8)*latvec(:,1)+real(k,8)*latvec(:,2)+real(m,8)*latvec(:,3)
!!      enddo
!!   enddo
!!enddo
!!
!!do m=0,2*nec3
!!   do k=0,2*nec2
!!      do l=0,2*nec1
!!      do iat=1,nat
!!      rxyzout(:,iat,l+1,k+1,m+1)=rxyz(:,iat)+transvecall(:,l+1,k+1,m+1)
!!      enddo
!!      enddo
!!   enddo
!!enddo
!!end

subroutine cart2polar(rxyz,r,theta,phi)
implicit none
real(8):: rxyz(3),t,theta,phi,r,d2,phi2
real(8), parameter :: pi=3.141592653589793238462643383279502884197d0

r=sqrt(rxyz(1)**2+rxyz(2)**2+rxyz(3)**2)
d2=rxyz(1)**2+rxyz(2)**2
theta=acos(rxyz(3)/r)
if(abs(d2).lt.1.d-15) then
phi=0.d0
else
phi=atan(abs(rxyz(2)/rxyz(1)))
if(rxyz(1).lt.0.d0) then
  if(rxyz(2).lt.0.d0) then
    phi=phi+pi
  else
    phi=pi-phi
  endif
else
  if(rxyz(2).lt.0.d0) then
    phi=2.d0*pi-phi
  endif  
endif
endif
phi2=acos(rxyz(1)/sqrt(d2))
!!if(rxyz(2).lt.0.d0) phi2=2.d0*pi-phi2
!!if(abs(phi-phi2).gt.1.d-10) then
!!write(*,*)abs(phi-phi2)
!!write(*,*) "phi",phi,phi2
!!write(*,*) rxyz(1),rxyz(2)
!!stop
!!endif
end subroutine


!!!
!!!      IMPLICIT NONE
!!!
!!!      DOUBLE PRECISION P(3), THETA, PHI
!!!
!!!      DOUBLE PRECISION X, Y, Z, D2
!!!
!!!* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!!
!!!      X = P(1)
!!!      Y = P(2)
!!!      Z = P(3)
!!!      D2 = X*X + Y*Y
!!!
!!!      IF ( D2 .EQ. 0D0 ) THEN
!!!         THETA = 0D0
!!!      ELSE
!!!         THETA = ATAN2(Y,X)
!!!      END IF
!!!
!!!      IF ( Z .EQ. 0D0 ) THEN
!!!         PHI = 0D0
!!!      ELSE
!!!         PHI = ATAN2(Z,SQRT(D2))
!!!      END IF
!!!



  subroutine ylm_ri (l, m, thrad, phirad, ylm_r, ylm_i)
  implicit none
!--------------------------------------------------------------------
! Computes the spherical harmonic Y_lm (theta,phi) using the
! reduced rotation matrix d^l_{m 0} (theta) and using the
! external function fac10(n) = factorial(n)/10**n
!--------------------------------------------------------------------
! input: angular momentum quantum numbers l, m (integers)
!        angles theta and phi (radian)
! -------------------------------------------------------------------
! Reference: D.M. Brink and G.R. Satchler, Angular Momentum,
!            second edition, Oxford University Press, p.22 and p. 145
! -------------------------------------------------------------------
  integer, parameter  :: wp = kind(1.0d0)  ! working precision = double (portable)
!--------------------------------------------------------------------
!   local constants
!--------------------------------------------------------------------
  real(wp), parameter :: pi = 3.14159265358979323846_wp    ! Schaum's Math handbook
!  complex(wp), parameter :: eye = (0.0_wp,1.0_wp)
!--------------------------------------------------------------------
!   formal arguments
!--------------------------------------------------------------------
  integer, intent(in)  :: l, m
  real(wp), intent(in) :: thrad, phirad
!--------------------------------------------------------------------
!   local variables
!--------------------------------------------------------------------
  integer :: itmin1, itmin2, itmin, itmax1, itmax2, itmax, it, iphase, &
             ia, ib, ic
  real(wp) :: sqrt_fac, sumt, denom, term, dlm0, const, cosb2, sinb2
  real(wp) :: ylm_r, ylm_i, exphi_r, exphi_i
!--------------------------------------------------------------------
!   external function
!--------------------------------------------------------------------
  real(wp), external :: fac10
!--------------------------------------------------------------------
!  program starts here
!  first calculate d^l_{m 0} (theta)
!--------------------------------------------------------------------
  cosb2 = cos(thrad/2.0_wp)
  sinb2 = sin(thrad/2.0_wp)
!--------------------------------------------------------------------
! determine lower and upper limits for summation index it; these
! are derived from the requirement that all factorials n! in the
! denominator are restricted to values with n >=0.
!--------------------------------------------------------------------
  itmin1 = 0
  itmin2 = m
  itmin = max(itmin1,itmin2)
  itmax1 = l+m
  itmax2 = l
  itmax = min(itmax1,itmax2)
!  write (6,'(10X,A,2I6)') ' itmin, itmax = ', itmin, itmax
  sqrt_fac = sqrt( fac10(l+m) * fac10(l-m) * fac10(l) * fac10(l) )
!
  sumt = 0.0_wp
  do it = itmin, itmax
     iphase = (-1)**it
     ia = l + m - it
     ib = l - it
     ic = it - m
!     write (6,'(10X,A,5I6)') ' it, iphase, ia, ib, ic  = ', it, iphase, ia, ib, ic
     denom = fac10(ia) * fac10(ib) * fac10(it) * fac10(ic)
     term = iphase * cosb2**(ia+ib) * sinb2**(it+ic) / denom
     sumt = sumt + term
  end do
  dlm0 = sqrt_fac * sumt
!--------------------------------------------------------------------
!  now compute Y_{l m} (theta,phi) from d^l_{m 0} (theta)
!--------------------------------------------------------------------
  const = sqrt( (2.0_wp *l + 1.0_wp) / (4.0_wp * pi) )
!  exphi = exp( eye * m * phirad )
  exphi_r = cos( m * phirad )
  exphi_i = sin( m * phirad )
  ylm_r = const * exphi_r * dlm0
  ylm_i = const * exphi_i * dlm0
!
  return
  end subroutine

  function fac10 (n)
  implicit none
! -----------------------------------------------
! function fac10(n) calculates factorial(n)/10**n
! -----------------------------------------------
! input: integer n >= 0 (you may want to check this
!        in the program calling this function)
! -----------------------------------------------
  integer, parameter :: wp = kind(1.0d0)  ! working precision = double (portable)
!------------------------------------------------
!      formal arguments
!------------------------------------------------
  integer, intent(in) :: n
!------------------------------------------------
!      local variables
!------------------------------------------------
  integer :: i
  real(wp) :: fac10, q
! -----------------------------------------------
  if (n == 0) then
     fac10 = 1.0_wp
  else
     fac10 = 1.0_wp
     q = 1.0_wp
     do i = 1, n
        fac10 = fac10 * q / 10.0_wp
        q = q + 1.0_wp
     end do
  endif
!
  return
  end function fac10
!  program main_ylm_test
!  implicit none
!!--------------------------------------------------------------------
!! test of function ylm (l, m, thrad, phirad)
!! last edited by V. Oberacker, October 13, 2010
!!--------------------------------------------------------------------
!  integer, parameter :: wp = kind(1.0d0)  ! working precision = double (portable)
!!--------------------------------------------------------------------
!!   local constants
!!--------------------------------------------------------------------
!  real(wp), parameter :: pi = 3.14159265358979323846_wp    ! Schaum's Math handbook
!  complex(wp), parameter :: eye = (0.0_wp,1.0_wp)
!!--------------------------------------------------------------------
!!   L o c a l   V a r i a b l e s
!!--------------------------------------------------------------------
!  integer :: l, m, itheta, iphi
!  real(wp) :: th_rad, phi_rad, sinth_rad, costh_rad, sphr, sphi
!  complex(wp) :: ylm_analyt, ylm_numeric, error, ylm_max
!!----- --------------------------------------------------------------
!!   E x t e r n a l   F u n c t i o n s
!!--------------------------------------------------------------------
!  complex(wp), external :: ylm
!! -------------------------------------------------------------------
!  open(unit=6, file='ylm_test.out', form='formatted', status='old', access='sequential')
!  write (6, '(/70(''*'')/A/70(''*''))') ' YLM_TEST Fortran 95 code, V.E. Oberacker'
!! -------------------------------------------------------------------
!! test Y_lm routine for various (l,m) combinations
!! For analytical values see Jackson, Classical Electrodynamics, p.109
!! -------------------------------------------------------------------
!  l = 0 ; m = 0
!  write (6, '(5/,A,2I5,2/)') ' spherical harmonic Y_lm for l,m = ', l, m
!  do itheta = 0, 180, 30
!     write (6, '(A,I5)') ' theta (deg) = ', itheta
!     th_rad = itheta * pi / 180.0_wp
!     write (6, '(A)') ' phi (deg), ylm_analyt, ylm_numeric, error'
!     do iphi = 0, 360, 60
!        phi_rad = iphi * pi / 180.0_wp
!        ylm_analyt = cmplx( 1.0_wp / sqrt(4.0_wp*pi) )
!        ylm_numeric = ylm (l, m, th_rad, phi_rad)
!        error = ylm_analyt - ylm_numeric
!        write (6, '(I5,3(5X,2ES14.6))') iphi, ylm_analyt, ylm_numeric, error
!        call ylm_ri (l, m, th_rad, phi_rad, sphr, sphi)
!        ylm_max=CMPLX(sphr, sphi)
!        error = ylm_numeric - ylm_max
!        write (6, '(I5,3(5X,2ES14.6))') iphi, ylm_analyt, ylm_max, error
!     end do
!  end do
!! -------------------------------------------------------------------
!  l = 1 ; m = 0
!  write (6, '(5/,A,2I5,2/)') ' spherical harmonic Y_lm for l,m = ', l, m
!  do itheta = 0, 180, 30
!     write (6, '(A,I5)') ' theta (deg) = ', itheta
!     th_rad = itheta * pi / 180.0_wp
!     costh_rad = cos(th_rad)
!     write (6, '(A)') ' phi (deg), ylm_analyt, ylm_numeric, error'
!     do iphi = 0, 360, 60
!        phi_rad = iphi * pi / 180.0_wp
!        ylm_analyt = cmplx( sqrt(3.0_wp/(4.0_wp*pi)) * costh_rad )
!        ylm_numeric = ylm (l, m, th_rad, phi_rad)
!        error = ylm_analyt - ylm_numeric
!        write (6, '(I5,3(5X,2ES14.6))') iphi, ylm_analyt, ylm_numeric, error
!        call ylm_ri (l, m, th_rad, phi_rad, sphr, sphi)
!        ylm_max=CMPLX(sphr, sphi)
!        error = ylm_numeric - ylm_max
!        write (6, '(I5,3(5X,2ES14.6))') iphi, ylm_analyt, ylm_max, error
!     end do
!  end do
!! -------------------------------------------------------------------
!  l = 2 ; m = 2
!  write (6, '(5/,A,2I5,2/)') ' spherical harmonic Y_lm for l,m = ', l, m
!  do itheta = 0, 180, 30
!     write (6, '(A,I5)') ' theta (deg) = ', itheta
!     th_rad = itheta * pi / 180.0_wp
!     sinth_rad = sin(th_rad)
!     write (6, '(A)') ' phi (deg), ylm_analyt, ylm_numeric, error'
!     do iphi = 0, 360, 60
!        phi_rad = iphi * pi / 180.0_wp
!        ylm_analyt = 0.25_wp * sqrt(15.0_wp/(2.0_wp*pi)) * sinth_rad * sinth_rad &
!                     * exp( eye * 2 * phi_rad )
!        ylm_numeric = ylm (l, m, th_rad, phi_rad)
!        error = ylm_analyt - ylm_numeric
!        write (6, '(I5,3(5X,2ES14.6))') iphi, ylm_analyt, ylm_numeric, error
!        call ylm_ri (l, m, th_rad, phi_rad, sphr, sphi)
!        ylm_max=CMPLX(sphr, sphi)
!        error = ylm_numeric - ylm_max
!        write (6, '(I5,3(5X,2ES14.6))') iphi, ylm_analyt, ylm_max, error
!     end do
!  end do
!! ---------------------------------
!  l = 3 ; m = 1
!  l = 3 ; m = -3
!  write (6, '(5/,A,2I5,2/)') ' spherical harmonic Y_lm for l,m = ', l, m
!  do itheta = 0, 180, 30
!     write (6, '(A,I5)') ' theta (deg) = ', itheta
!     th_rad = itheta * pi / 180.0_wp
!     sinth_rad = sin(th_rad)
!     costh_rad = cos(th_rad)
!     write (6, '(A)') ' phi (deg), ylm_analyt, ylm_numeric, error'
!     do iphi = 0, 360, 60
!        phi_rad = iphi * pi / 180.0_wp
!        ylm_analyt = -0.25_wp * sqrt(21.0_wp/(4.0_wp*pi)) * sinth_rad &
!                     * (5.0_wp*costh_rad*costh_rad -1.0_wp) * exp(eye*phi_rad)
!        ylm_numeric = ylm (l, m, th_rad, phi_rad)
!        error = ylm_analyt - ylm_numeric
!        write (6, '(I5,3(5X,2ES14.6))') iphi, ylm_analyt, ylm_numeric, error
!        call ylm_ri (l, m, th_rad, phi_rad, sphr, sphi)
!        ylm_max=CMPLX(sphr, sphi)
!        error = ylm_numeric - ylm_max
!        write (6, '(I5,3(5X,2ES14.6))') iphi, ylm_analyt, ylm_max, error
!     end do
!  end do
!!
!! ---------------------------------
!  l = 30 ; m = -25
!  write (6, '(5/,A,2I5,2/)') ' spherical harmonic Y_lm for l,m = ', l, m
!  do itheta = 0, 180, 30
!     write (6, '(A,I5)') ' theta (deg) = ', itheta
!     th_rad = itheta * pi / 180.0_wp
!     sinth_rad = sin(th_rad)
!     costh_rad = cos(th_rad)
!     write (6, '(A)') ' phi (deg), ylm_analyt, ylm_numeric, error'
!     do iphi = 0, 360, 60
!        phi_rad = iphi * pi / 180.0_wp
!        ylm_analyt = -0.25_wp * sqrt(21.0_wp/(4.0_wp*pi)) * sinth_rad &
!                     * (5.0_wp*costh_rad*costh_rad -1.0_wp) * exp(eye*phi_rad)
!        ylm_numeric = ylm (l, m, th_rad, phi_rad)
!        error = ylm_analyt - ylm_numeric
!        write (6, '(I5,3(5X,2ES14.6))') iphi, ylm_analyt, ylm_numeric, error
!        call ylm_ri (l, m, th_rad, phi_rad, sphr, sphi)
!        ylm_max=CMPLX(sphr, sphi)
!        error = ylm_numeric - ylm_max
!        write (6, '(I5,3(5X,2ES14.6))') iphi, ylm_analyt, ylm_max, error
!     end do
!  end do
!  stop
!  end program main_ylm_test

subroutine poly1(r2,rc1,rc2,val)
!This subroutine computes the parameters for a 3rd order polynomial for which the value is
!1 at rc1 and 0 at rc2, and the first derivatives vanish at rc1 and rc2.
!r2 is the value r^2
implicit none
real(8):: r2,a0,a1,a2,a3,rc1,rc2,r12_3,r,val
if(r2.le.rc1*rc1) then
val=1.d0
elseif(r2.ge.rc2*rc2) then
val=0.d0
else
r12_3=1.d0/(rc1-rc2)**3
a0=(3.d0*rc1*rc2**2-rc2**3)*r12_3
a1=-(6.d0*rc1*rc2)*r12_3
a2=3.d0*(rc1+rc2)*r12_3
a3=-2.d0*r12_3
r=sqrt(r2)
val=a0+a1*r+a2*r2+a3*r2*r
endif
end subroutine

subroutine poly2(r2,rc1,rc2,val,scaling)
!This subroutine computes the matrix element but based on a simple 1/r function
!scaled by r1 and cutoff at r2
implicit none
real(8):: r2,a0,a1,a2,a3,rc1,rc2,r12_3,r,val,scaling
if(r2.le.0.1d0*0.1d0*rc1*rc1) then
r=sqrt(r2)
val=scaling/(0.1d0*rc1)
elseif(r2.ge.rc2*rc2) then
val=0.d0
else
r=sqrt(r2)
val=scaling/r
endif
end subroutine
