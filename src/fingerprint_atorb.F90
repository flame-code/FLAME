subroutine get_fp_malypso(nat,rxyz,rcov,latvec,r_cut_in,kinds,nkinds,fp_dim,nl,fp)
!This routine will compute the Bond Characterization Matrix as proposed by 
!Yanchao Wang in computer physics communications 183, 2063
!The dimension of fp is a matrix with dimension of AB interactions and nl spherical harmonics (be careful: if nl=10 it means lmax=(nl-1)*2
!fp_dim=ntypat*(ntypat+1)/2
implicit none
integer:: nl !Number of l components, here only even ones 
integer:: fp_dim !Number of AB interactions, doublecounting eliminated
integer:: nat,nkinds
integer:: nbond(fp_dim,nat),kinds(nat)
integer:: lmax,yll,ymm,yl,ym
integer:: imax,imin
integer:: nec1,nec2,nec3,iat,jat,k,l,m,iat_kind,jat_kind,i,j,iarr
integer:: llist(nl)
real(8):: latvec(3,3),rxyz(3,nat),fp(nl,fp_dim,nat),fp_ri(2,nl,nat),rcov(nkinds)
real(8), allocatable:: rxyzexp(:,:,:,:,:),transvecall(:,:,:,:)
real(8):: sigma,r_cut_in(fp_dim) !Cutoff for each AB interaction
real(8):: r_cut(fp_dim,fp_dim)
real(8), parameter :: pi=3.141592653589793238462643383279502884197d0
real(8):: min_bond(fp_dim,fp_dim),ylm_r,ylm_i
real(8):: rij, rij_vec(3),exp_fac,rxyzj(3),theta,phi,eps,tr,ti
real(8),allocatable:: qlm(:,:,:,:,:)
real(8):: double_counting,rcmax,minmax,tmp
!Initiallize fp
fp=0.d0
fp_ri=0.d0

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
allocate(qlm(2,nl,2*lmax+1,fp_dim,nat))
qlm=0.d0


!Create the expanded unit cells
rcmax=maxval(r_cut_in)
call n_rep_dim(latvec,rcmax,nec1,nec2,nec3)
write(*,'(a,3(i3,1x),a)') " Creating expansion for periodic ATORB FP ",nec1,nec2,nec3,"..."
allocate(rxyzexp(3,nat,2*nec1+1,2*nec2+1,2*nec3+1),transvecall(3,2*nec1+1,2*nec2+1,2*nec3+1))
call expand_dim(rxyz,rxyzexp,transvecall,latvec,nat,nec1,nec2,nec3)


!Min Bond feed with rcov
do i=1,fp_dim
do j=1,fp_dim
!!        imax=max(i,j)
!!        imin=min(i,j)
!!        iarr=imax*(imax+1)
!!        iarr=(iarr/2+imin-imax)
!!        min_bond(iarr)=rcov(imin)+rcov(imax)
        min_bond(i,j)=rcov(i)+rcov(j)
enddo
enddo

!Define the rcut based on rcov...
minmax=maxval(min_bond(:,:))
tmp=rcmax/minmax
do i=1,fp_dim
do j=1,fp_dim
   r_cut(i,j)=min_bond(i,j)*tmp
enddo
enddo

!Get the Ylm and Q_lm^dab
nbond=0
do iat=1,nat
  do jat=1,nat
     do k=1,2*nec1+1
     do l=1,2*nec2+1
     do m=1,2*nec3+1
      if(iat==jat.and.k==nec1+1.and.l==nec2+1.and.m==nec3+1) goto 2020
        rxyzj(:)=rxyz(:,jat)+transvecall(:,k,l,m)
        rij_vec=rxyzj(:)-rxyz(:,iat)
        call cart2polar(rij_vec,rij,theta,phi)
!!        imax=max(Kinds(iat),Kinds(jat))
!!        imin=min(Kinds(iat),Kinds(jat))
!!        iarr=imax*(imax+1)
!!        iarr=(iarr/2+imin-imax)
!!        if(rij.lt.r_cut(iarr)) then
        if(rij.lt.r_cut(kinds(iat),kinds(jat))) then
!!           nbond(iarr,iat)=nbond(iarr,iat)+1
           nbond(kinds(jat),iat)=nbond(kinds(jat),iat)+1
!!           double_counting=1.d0/min(real(imax-imin,8)+1.d0,2.d0)
           double_counting=1.d0
!!                exp_fac=exp(-alpha(iarr)*(rij-min_bond(iarr)))*double_counting
!!                exp_fac=1.d0*double_counting
                call poly1(rij**2,min_bond(kinds(iat),kinds(jat)),r_cut(kinds(iat),kinds(jat)),exp_fac)
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
                        qlm(1,yll,ymm,kinds(jat),iat)=qlm(1,yll,ymm,kinds(jat),iat)+ylm_r*exp_fac
                        qlm(2,yll,ymm,kinds(jat),iat)=qlm(2,yll,ymm,kinds(jat),iat)+ylm_i*exp_fac
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
do iat=1,nat
do i=1,fp_dim
  qlm(:,:,:,i,iat)=qlm(:,:,:,i,iat)/real(nbond(i,iat),8)
enddo
enddo

!Sum over all m components and get fp
do iat=1,nat
  do i=1,fp_dim
  do yll=1,nl
     yl=llist(yll)
     do ymm=1,2*yl+1
        fp(yll,i,iat)=fp(yll,i,iat)+qlm(1,yll,ymm,i,iat)**2+qlm(2,yll,ymm,i,iat)**2
     enddo
  fp(yll,i,iat)=fp(yll,i,iat)*4.d0*pi/(2.d0*real(yl,8)+1.d0)
  fp(yll,i,iat)=sqrt(fp(yll,i,iat))
  enddo
  enddo
enddo

!!write(*,*) "fp",fp
!!write(*,*) "nbond",nbond
!!write(*,*) "minbond",min_bond
deallocate(qlm,rxyzexp,transvecall)
end subroutine

subroutine get_distance_malypso(fp1,fp2,fp_dim,nat,kinds,nl,dist)
implicit none
integer:: fp_dim,nl,yll,i,j,mode,nat,iarr,i_dim,ii,jj,kinds(nat),k(nat),iii,jjj,nmat,imax,imin
real(8):: fp1(nl,fp_dim,nat),fp2(nl,fp_dim,nat),dist,a(nat,nat),summ,vec(nl),vec1(nl),vec2(nl),norm1,norm2
real(8),allocatable:: aa(:,:,:)
integer,allocatable:: aa_dim(:),aa_list(:,:)
integer:: nn_i,nn_j


nmat=fp_dim*(fp_dim+1)/2
allocate(aa_dim(nmat),aa(nat,nat,nmat),aa_list(nat,nmat))

!!write(*,*) "fp1",fp1
!!write(*,*) "fp2",fp2
dist=0.d0
aa_dim=0
aa=0.d0
aa_list=0

!Set up matrices
do i=1,nat
   do j=1,fp_dim
        imax=max(Kinds(i),j)
        imin=min(Kinds(i),j)
        iarr=imax*(imax+1)
        iarr=(iarr/2+imin-imax)
        aa_dim(iarr)=aa_dim(iarr)+1
        aa_list(aa_dim(iarr),iarr)=i
   enddo
enddo

!Euclidic distance
do iarr=1,nmat !Loop for all interactions
aa=0.d0
!Stephan's clever formula...
ii=ceiling(-0.5d0+0.5d0*sqrt(1.d0+8.d0*real(iarr,8)))
jj=(ii*(1-ii))/2+iarr
do i=1,aa_dim(iarr) !Loop over all atoms in this matrix
do j=1,aa_dim(iarr) !Loop over all atoms in this matrix
   iii=ii
   if(kinds(aa_list(i,iarr)).gt.jj) iii=jj 
   jjj=ii
   if(kinds(aa_list(j,iarr)).gt.jj) jjj=jj 
!Euclidic
   vec=fp1(:,iii,aa_list(i,iarr))-fp2(:,jjj,aa_list(j,iarr))
   aa(i,j,iarr)=dot_product(vec,vec) !Euclydic
!Cosine
!!   vec1=fp1(:,iii,aa_list(i,iarr))
!!   vec2=fp2(:,jjj,aa_list(j,iarr))
!!   norm1=dot_product(vec1,vec1) 
!!   norm2=dot_product(vec2,vec2) 
!!   aa(i,j,iarr)=abs(0.5d0*(1.d0-dot_product(vec1,vec2)/sqrt(norm1*norm2))/real(aa_dim(iarr),8))
enddo
enddo
mode=1
!write(*,*)  aa(1:aa_dim(iarr),1:aa_dim(iarr),iarr)
call assndx(mode, aa(1:aa_dim(iarr),1:aa_dim(iarr),iarr), aa_dim(iarr), aa_dim(iarr), k(1:aa_dim(iarr)), summ)
!dist=dist+sqrt(real(summ,8)*1.d-6)/real(fp_dim,8)
dist=dist+summ!/real(nmat,8)
!dist=dist+sqrt(real(summ,8))/real(fp_dim,8)
enddo
!Euclidic
dist=dist/real(nat,8)
dist=sqrt(abs(dist))
!Cosine
!dist=dist/real(nmat,8)
contains
   





!!!ii=ceiling(-0.5d0+0.5d0*sqrt(1.d0+8.d0*real(i_dim,8)))
!!!jj=(ii*(1-ii))/2+i_dim
!!!a=0.d0
!!!nn_i=0
!!!nn_j=0
!!!   do i=1,nat
!!!   do j=1,nat
!!!!    if(min(typat(i),typat(j))==jj.and.max(typat(i),typat(j))==ii) a(i,j)=dot_product(fp1(:,i_dim,i)-fp2(:,i_dim,j),fp1(:,i_dim,i)-fp2(:,i_dim,j))
!!!    vec=fp1(:,i_dim,i)-fp2(:,i_dim,j)
!!!    a(i,j)=dot_product(vec,vec)
!!!!    write(*,*) i,j,vec
!!!!    a(j,i)=a(i,j)
!!!!    write(*,*) ii,jj,typat(i),typat(j)
!!!   enddo 
!!!   enddo 
!!!!!write(*,'(a,4(es15.7))') "A",A(:,1)
!!!!!write(*,'(a,4(es15.7))') "A",A(:,2)
!!!!!write(*,'(a,4(es15.7))') "A",A(:,3)
!!!!!write(*,'(a,4(es15.7))') "A",A(:,4)
!!!mode=1
!!!call assndx(mode, a, nat, nat, k, summ)
!!!!!aa=int(1.d6*A)
!!!!!call APC(nat,aa,K,summ)
!!!write(*,*)"SUMM",summ
!!!!dist=dist+sqrt(real(summ,8)*1.d-6)/real(fp_dim,8)
!!!dist=dist+summ/real(fp_dim,8)
!!!!dist=dist+sqrt(real(summ,8))/real(fp_dim,8)
!!!enddo
!!!dist=sqrt(dist)
!!!contains
SUBROUTINE assndx(mode, a, n, m, k, sum)
! N.B. Arguments IDA, IW & IDW have been removed.

! If MODE = 1, then it finds k(1), k(2), ..., k(n) to minimize
!        S = Sum(i=1, .., n) a(i, k(i))
! If MODE = 2,  then it finds k(1), k(2), ..., k(m) to minimize
!        S = Sum(j=1, .., m) a(k(j), j)
! given the array a(n,m).

! References:
! Munkres, J. (1957) `Algorithms for the assignment and transportation problems',
!                    J. SIAM, vol.5, 32-38.
! Silver, R. (1960) `An algorithm for the assignment problem', Comm. ACM, vol.3,
!                   605-606.   The algorithm (CACM 27) is in Algol.

IMPLICIT NONE

INTEGER, INTENT(IN)   :: mode
REAL(8), INTENT(IN OUT)  :: a(:,:)
INTEGER, INTENT(IN)   :: n
INTEGER, INTENT(IN)   :: m
INTEGER, INTENT(OUT)  :: k(:)
REAL(8), INTENT(OUT)  :: sum

LOGICAL  :: lsw
INTEGER  :: i, icbl, icl, icl0, iflag, imax, imin, ipp, irl, irs, &
            j, j1, jsv, new
REAL(8)  :: rmin
INTEGER, ALLOCATABLE  :: iw(:,:)

IF (n < 1 .OR. m < 1) THEN
  WRITE(*, '(a, 2i8)') ' ** Error in call to ASSNDX; m, n = ', m, n
  RETURN
END IF

imax = MAX(n,m)
imin = MIN(n,m)
ALLOCATE( iw(imax,6) )
sum = 0.0
IF (n <= m) THEN
  DO  i = 1, n
    rmin = a(i,1)
    DO  j = 1, m
      rmin = MIN(rmin, a(i,j))
    END DO
    sum = sum + rmin
    a(i,1:m) = a(i,1:m) - rmin
  END DO
END IF
IF (n >= m) THEN
  DO  j = 1, m
    rmin = a(1,j)
    DO  i = 1, n
      rmin = MIN(rmin,a(i,j))
    END DO
    sum = sum + rmin
    a(1:n,j) = a(1:n,j) - rmin
  END DO
END IF

DO  i = 1, imax
  k(i) = 0
  iw(i,1) = 0
END DO

loop90:  DO  i = 1, n
  DO  j = 1, m
    IF (a(i,j)+iw(j,1) == 0) THEN
      k(i) = j
      iw(j,1) = i
      CYCLE loop90
    END IF
  END DO
END DO loop90

100 iflag = n
irl = 0
icl = 0
irs = 1

DO  i = 1, n
  iw(i,5) = 0
  IF (k(i) == 0) THEN
    irl = irl + 1
    iw(irl,6) = i
    iw(i,5) = -1
    iflag = iflag - 1
  END IF
END DO
IF (iflag == imin) THEN
  IF (mode == 2) k(1:imax) = iw(1:imax,1)
  RETURN
END IF

iw(1:m,4) = 0

140 i = iw(irs,6)
irs = irs + 1
DO  j = 1, m
  IF (a(i,j)+iw(j,4) == 0) THEN
    iw(j,4) = i
    icl = icl + 1
    iw(icl,2) = j
    NEW = iw(j,1)
    IF (NEW == 0) THEN
      j1 = j
      DO
        iw(j1,1) = iw(j1,4)
        i = iw(j1,4)
        IF (k(i) == 0) THEN
          k(i) = j1
          GO TO 100
        END IF
        jsv = j1
        j1 = k(i)
        k(i) = jsv
      END DO
    END IF
    irl = irl + 1
    iw(irl,6) = NEW
    iw(NEW,5) = i
  END IF
END DO
IF (irs <= irl) GO TO 140

lsw = .true.
icl0 = icl
icbl = 0
DO  j = 1, m
  IF (iw(j,4) == 0) THEN
    icbl = icbl + 1
    iw(icbl,3) = j
  END IF
END DO
rmin = a(iw(1,6),iw(1,3))
DO  i = 1, irl
  DO  j = 1, icbl
    rmin = MIN(rmin, a(iw(i,6), iw(j,3)))
  END DO
END DO
sum = sum + rmin * (irl+icbl-imax)

DO  i = 1, n
  IF (iw(i,5) == 0) THEN
    DO  ipp = 1, icl0
      a(i,iw(ipp,2)) = a(i,iw(ipp,2)) + rmin
    END DO
    CYCLE
  END IF
  DO  ipp = 1, icbl
    NEW = iw(ipp,3)
    a(i,NEW) = a(i,NEW) - rmin
    IF (lsw.AND.a(i,NEW)+iw(NEW,4) == 0) THEN
      iw(NEW,4) = i
      IF (iw(NEW,1) == 0) THEN
        j1 = NEW
        lsw = .false.
      ELSE
        icl = icl + 1
        iw(icl,2) = NEW
        irl = irl + 1
        iw(irl,6) = iw(NEW,1)
      END IF
    END IF
  END DO
END DO

IF (lsw) THEN
  DO  i = icl0 + 1, icl
    iw(iw(iw(i,2),1),5) = iw(i,2)
  END DO
  GO TO 140
ELSE
  DO
    iw(j1,1) = iw(j1,4)
    i = iw(j1,4)
    IF (k(i) == 0) THEN
      k(i) = j1
      GO TO 100
    END IF
    jsv = j1
    j1 = k(i)
    k(i) = jsv
  END DO
END IF

RETURN
END SUBROUTINE assndx
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

!!!!subroutine cart2polar(rxyz,r,theta,phi)
!!!!implicit none
!!!!real(8):: rxyz(3),t,theta,phi,r,d2,phi2
!!!!real(8), parameter :: pi=3.141592653589793238462643383279502884197d0
!!!!
!!!!r=sqrt(rxyz(1)**2+rxyz(2)**2+rxyz(3)**2)
!!!!d2=rxyz(1)**2+rxyz(2)**2
!!!!theta=acos(rxyz(3)/r)
!!!!if(abs(d2).lt.1.d-15) then
!!!!phi=0.d0
!!!!else
!!!!phi=atan(abs(rxyz(2)/rxyz(1)))
!!!!if(rxyz(1).lt.0.d0) then
!!!!  if(rxyz(2).lt.0.d0) then
!!!!    phi=phi+pi
!!!!  else
!!!!    phi=pi-phi
!!!!  endif
!!!!else
!!!!  if(rxyz(2).lt.0.d0) then
!!!!    phi=2.d0*pi-phi
!!!!  endif  
!!!!endif
!!!!endif
!!!!phi2=acos(rxyz(1)/sqrt(d2))
!!!!!!if(rxyz(2).lt.0.d0) phi2=2.d0*pi-phi2
!!!!!!if(abs(phi-phi2).gt.1.d-10) then
!!!!!!write(*,*)abs(phi-phi2)
!!!!!!write(*,*) "phi",phi,phi2
!!!!!!write(*,*) rxyz(1),rxyz(2)
!!!!!!stop
!!!!!!endif
!!!!end subroutine
!!!!

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



!!!  subroutine ylm_ri (l, m, thrad, phirad, ylm_r, ylm_i)
!!!  implicit none
!!!!--------------------------------------------------------------------
!!!! Computes the spherical harmonic Y_lm (theta,phi) using the
!!!! reduced rotation matrix d^l_{m 0} (theta) and using the
!!!! external function fac10(n) = factorial(n)/10**n
!!!!--------------------------------------------------------------------
!!!! input: angular momentum quantum numbers l, m (integers)
!!!!        angles theta and phi (radian)
!!!! -------------------------------------------------------------------
!!!! Reference: D.M. Brink and G.R. Satchler, Angular Momentum,
!!!!            second edition, Oxford University Press, p.22 and p. 145
!!!! -------------------------------------------------------------------
!!!  integer, parameter  :: wp = kind(1.0d0)  ! working precision = double (portable)
!!!!--------------------------------------------------------------------
!!!!   local constants
!!!!--------------------------------------------------------------------
!!!  real(wp), parameter :: pi = 3.14159265358979323846_wp    ! Schaum's Math handbook
!!!!  complex(wp), parameter :: eye = (0.0_wp,1.0_wp)
!!!!--------------------------------------------------------------------
!!!!   formal arguments
!!!!--------------------------------------------------------------------
!!!  integer, intent(in)  :: l, m
!!!  real(wp), intent(in) :: thrad, phirad
!!!!--------------------------------------------------------------------
!!!!   local variables
!!!!--------------------------------------------------------------------
!!!  integer :: itmin1, itmin2, itmin, itmax1, itmax2, itmax, it, iphase, &
!!!             ia, ib, ic
!!!  real(wp) :: sqrt_fac, sumt, denom, term, dlm0, const, cosb2, sinb2
!!!  real(wp) :: ylm_r, ylm_i, exphi_r, exphi_i
!!!!--------------------------------------------------------------------
!!!!   external function
!!!!--------------------------------------------------------------------
!!!  real(wp), external :: fac10
!!!!--------------------------------------------------------------------
!!!!  program starts here
!!!!  first calculate d^l_{m 0} (theta)
!!!!--------------------------------------------------------------------
!!!  cosb2 = cos(thrad/2.0_wp)
!!!  sinb2 = sin(thrad/2.0_wp)
!!!!--------------------------------------------------------------------
!!!! determine lower and upper limits for summation index it; these
!!!! are derived from the requirement that all factorials n! in the
!!!! denominator are restricted to values with n >=0.
!!!!--------------------------------------------------------------------
!!!  itmin1 = 0
!!!  itmin2 = m
!!!  itmin = max(itmin1,itmin2)
!!!  itmax1 = l+m
!!!  itmax2 = l
!!!  itmax = min(itmax1,itmax2)
!!!!  write (6,'(10X,A,2I6)') ' itmin, itmax = ', itmin, itmax
!!!  sqrt_fac = sqrt( fac10(l+m) * fac10(l-m) * fac10(l) * fac10(l) )
!!!!
!!!  sumt = 0.0_wp
!!!  do it = itmin, itmax
!!!     iphase = (-1)**it
!!!     ia = l + m - it
!!!     ib = l - it
!!!     ic = it - m
!!!!     write (6,'(10X,A,5I6)') ' it, iphase, ia, ib, ic  = ', it, iphase, ia, ib, ic
!!!     denom = fac10(ia) * fac10(ib) * fac10(it) * fac10(ic)
!!!     term = iphase * cosb2**(ia+ib) * sinb2**(it+ic) / denom
!!!     sumt = sumt + term
!!!  end do
!!!  dlm0 = sqrt_fac * sumt
!!!!--------------------------------------------------------------------
!!!!  now compute Y_{l m} (theta,phi) from d^l_{m 0} (theta)
!!!!--------------------------------------------------------------------
!!!  const = sqrt( (2.0_wp *l + 1.0_wp) / (4.0_wp * pi) )
!!!!  exphi = exp( eye * m * phirad )
!!!  exphi_r = cos( m * phirad )
!!!  exphi_i = sin( m * phirad )
!!!  ylm_r = const * exphi_r * dlm0
!!!  ylm_i = const * exphi_i * dlm0
!!!!
!!!  return
!!!  end subroutine
!!!
!!!  function fac10 (n)
!!!  implicit none
!!!! -----------------------------------------------
!!!! function fac10(n) calculates factorial(n)/10**n
!!!! -----------------------------------------------
!!!! input: integer n >= 0 (you may want to check this
!!!!        in the program calling this function)
!!!! -----------------------------------------------
!!!  integer, parameter :: wp = kind(1.0d0)  ! working precision = double (portable)
!!!!------------------------------------------------
!!!!      formal arguments
!!!!------------------------------------------------
!!!  integer, intent(in) :: n
!!!!------------------------------------------------
!!!!      local variables
!!!!------------------------------------------------
!!!  integer :: i
!!!  real(wp) :: fac10, q
!!!! -----------------------------------------------
!!!  if (n == 0) then
!!!     fac10 = 1.0_wp
!!!  else
!!!     fac10 = 1.0_wp
!!!     q = 1.0_wp
!!!     do i = 1, n
!!!        fac10 = fac10 * q / 10.0_wp
!!!        q = q + 1.0_wp
!!!     end do
!!!  endif
!!!!
!!!  return
!!!  end function fac10
!!!!  program main_ylm_test
!!!!  implicit none
!!!!!--------------------------------------------------------------------
!!!!! test of function ylm (l, m, thrad, phirad)
!!!!! last edited by V. Oberacker, October 13, 2010
!!!!!--------------------------------------------------------------------
!!!!  integer, parameter :: wp = kind(1.0d0)  ! working precision = double (portable)
!!!!!--------------------------------------------------------------------
!!!!!   local constants
!!!!!--------------------------------------------------------------------
!!!!  real(wp), parameter :: pi = 3.14159265358979323846_wp    ! Schaum's Math handbook
!!!!  complex(wp), parameter :: eye = (0.0_wp,1.0_wp)
!!!!!--------------------------------------------------------------------
!!!!!   L o c a l   V a r i a b l e s
!!!!!--------------------------------------------------------------------
!!!!  integer :: l, m, itheta, iphi
!!!!  real(wp) :: th_rad, phi_rad, sinth_rad, costh_rad, sphr, sphi
!!!!  complex(wp) :: ylm_analyt, ylm_numeric, error, ylm_max
!!!!!----- --------------------------------------------------------------
!!!!!   E x t e r n a l   F u n c t i o n s
!!!!!--------------------------------------------------------------------
!!!!  complex(wp), external :: ylm
!!!!! -------------------------------------------------------------------
!!!!  open(unit=6, file='ylm_test.out', form='formatted', status='old', access='sequential')
!!!!  write (6, '(/70(''*'')/A/70(''*''))') ' YLM_TEST Fortran 95 code, V.E. Oberacker'
!!!!! -------------------------------------------------------------------
!!!!! test Y_lm routine for various (l,m) combinations
!!!!! For analytical values see Jackson, Classical Electrodynamics, p.109
!!!!! -------------------------------------------------------------------
!!!!  l = 0 ; m = 0
!!!!  write (6, '(5/,A,2I5,2/)') ' spherical harmonic Y_lm for l,m = ', l, m
!!!!  do itheta = 0, 180, 30
!!!!     write (6, '(A,I5)') ' theta (deg) = ', itheta
!!!!     th_rad = itheta * pi / 180.0_wp
!!!!     write (6, '(A)') ' phi (deg), ylm_analyt, ylm_numeric, error'
!!!!     do iphi = 0, 360, 60
!!!!        phi_rad = iphi * pi / 180.0_wp
!!!!        ylm_analyt = cmplx( 1.0_wp / sqrt(4.0_wp*pi) )
!!!!        ylm_numeric = ylm (l, m, th_rad, phi_rad)
!!!!        error = ylm_analyt - ylm_numeric
!!!!        write (6, '(I5,3(5X,2ES14.6))') iphi, ylm_analyt, ylm_numeric, error
!!!!        call ylm_ri (l, m, th_rad, phi_rad, sphr, sphi)
!!!!        ylm_max=CMPLX(sphr, sphi)
!!!!        error = ylm_numeric - ylm_max
!!!!        write (6, '(I5,3(5X,2ES14.6))') iphi, ylm_analyt, ylm_max, error
!!!!     end do
!!!!  end do
!!!!! -------------------------------------------------------------------
!!!!  l = 1 ; m = 0
!!!!  write (6, '(5/,A,2I5,2/)') ' spherical harmonic Y_lm for l,m = ', l, m
!!!!  do itheta = 0, 180, 30
!!!!     write (6, '(A,I5)') ' theta (deg) = ', itheta
!!!!     th_rad = itheta * pi / 180.0_wp
!!!!     costh_rad = cos(th_rad)
!!!!     write (6, '(A)') ' phi (deg), ylm_analyt, ylm_numeric, error'
!!!!     do iphi = 0, 360, 60
!!!!        phi_rad = iphi * pi / 180.0_wp
!!!!        ylm_analyt = cmplx( sqrt(3.0_wp/(4.0_wp*pi)) * costh_rad )
!!!!        ylm_numeric = ylm (l, m, th_rad, phi_rad)
!!!!        error = ylm_analyt - ylm_numeric
!!!!        write (6, '(I5,3(5X,2ES14.6))') iphi, ylm_analyt, ylm_numeric, error
!!!!        call ylm_ri (l, m, th_rad, phi_rad, sphr, sphi)
!!!!        ylm_max=CMPLX(sphr, sphi)
!!!!        error = ylm_numeric - ylm_max
!!!!        write (6, '(I5,3(5X,2ES14.6))') iphi, ylm_analyt, ylm_max, error
!!!!     end do
!!!!  end do
!!!!! -------------------------------------------------------------------
!!!!  l = 2 ; m = 2
!!!!  write (6, '(5/,A,2I5,2/)') ' spherical harmonic Y_lm for l,m = ', l, m
!!!!  do itheta = 0, 180, 30
!!!!     write (6, '(A,I5)') ' theta (deg) = ', itheta
!!!!     th_rad = itheta * pi / 180.0_wp
!!!!     sinth_rad = sin(th_rad)
!!!!     write (6, '(A)') ' phi (deg), ylm_analyt, ylm_numeric, error'
!!!!     do iphi = 0, 360, 60
!!!!        phi_rad = iphi * pi / 180.0_wp
!!!!        ylm_analyt = 0.25_wp * sqrt(15.0_wp/(2.0_wp*pi)) * sinth_rad * sinth_rad &
!!!!                     * exp( eye * 2 * phi_rad )
!!!!        ylm_numeric = ylm (l, m, th_rad, phi_rad)
!!!!        error = ylm_analyt - ylm_numeric
!!!!        write (6, '(I5,3(5X,2ES14.6))') iphi, ylm_analyt, ylm_numeric, error
!!!!        call ylm_ri (l, m, th_rad, phi_rad, sphr, sphi)
!!!!        ylm_max=CMPLX(sphr, sphi)
!!!!        error = ylm_numeric - ylm_max
!!!!        write (6, '(I5,3(5X,2ES14.6))') iphi, ylm_analyt, ylm_max, error
!!!!     end do
!!!!  end do
!!!!! ---------------------------------
!!!!  l = 3 ; m = 1
!!!!  l = 3 ; m = -3
!!!!  write (6, '(5/,A,2I5,2/)') ' spherical harmonic Y_lm for l,m = ', l, m
!!!!  do itheta = 0, 180, 30
!!!!     write (6, '(A,I5)') ' theta (deg) = ', itheta
!!!!     th_rad = itheta * pi / 180.0_wp
!!!!     sinth_rad = sin(th_rad)
!!!!     costh_rad = cos(th_rad)
!!!!     write (6, '(A)') ' phi (deg), ylm_analyt, ylm_numeric, error'
!!!!     do iphi = 0, 360, 60
!!!!        phi_rad = iphi * pi / 180.0_wp
!!!!        ylm_analyt = -0.25_wp * sqrt(21.0_wp/(4.0_wp*pi)) * sinth_rad &
!!!!                     * (5.0_wp*costh_rad*costh_rad -1.0_wp) * exp(eye*phi_rad)
!!!!        ylm_numeric = ylm (l, m, th_rad, phi_rad)
!!!!        error = ylm_analyt - ylm_numeric
!!!!        write (6, '(I5,3(5X,2ES14.6))') iphi, ylm_analyt, ylm_numeric, error
!!!!        call ylm_ri (l, m, th_rad, phi_rad, sphr, sphi)
!!!!        ylm_max=CMPLX(sphr, sphi)
!!!!        error = ylm_numeric - ylm_max
!!!!        write (6, '(I5,3(5X,2ES14.6))') iphi, ylm_analyt, ylm_max, error
!!!!     end do
!!!!  end do
!!!!!
!!!!! ---------------------------------
!!!!  l = 30 ; m = -25
!!!!  write (6, '(5/,A,2I5,2/)') ' spherical harmonic Y_lm for l,m = ', l, m
!!!!  do itheta = 0, 180, 30
!!!!     write (6, '(A,I5)') ' theta (deg) = ', itheta
!!!!     th_rad = itheta * pi / 180.0_wp
!!!!     sinth_rad = sin(th_rad)
!!!!     costh_rad = cos(th_rad)
!!!!     write (6, '(A)') ' phi (deg), ylm_analyt, ylm_numeric, error'
!!!!     do iphi = 0, 360, 60
!!!!        phi_rad = iphi * pi / 180.0_wp
!!!!        ylm_analyt = -0.25_wp * sqrt(21.0_wp/(4.0_wp*pi)) * sinth_rad &
!!!!                     * (5.0_wp*costh_rad*costh_rad -1.0_wp) * exp(eye*phi_rad)
!!!!        ylm_numeric = ylm (l, m, th_rad, phi_rad)
!!!!        error = ylm_analyt - ylm_numeric
!!!!        write (6, '(I5,3(5X,2ES14.6))') iphi, ylm_analyt, ylm_numeric, error
!!!!        call ylm_ri (l, m, th_rad, phi_rad, sphr, sphi)
!!!!        ylm_max=CMPLX(sphr, sphi)
!!!!        error = ylm_numeric - ylm_max
!!!!        write (6, '(I5,3(5X,2ES14.6))') iphi, ylm_analyt, ylm_max, error
!!!!     end do
!!!!  end do
!!!!  stop
!!!!  end program main_ylm_test

!!      SUBROUTINE APC(N,A,F,Z)
!!!
!!! SOLUTION OF THE LINEAR MIN-SUM ASSIGNMENT PROBLEM.
!!!
!!! HUNGARIAN METHOD. COMPLEXITY O(N**3).
!!!
!!!
!!! MEANING OF THE INPUT PARAMETERS:
!!! N      = NUMBER OF ROWS AND COLUMNS OF THE COST MATRIX.
!!! A(I,J) = COST OF THE ASSIGNMENT OF ROW  I  TO COLUMN  J .
!!! ON RETURN, THE INPUT PARAMETERS ARE UNCHANGED.
!!!
!!! MEANING OF THE OUTPUT PARAMETERS:
!!! F(I) = COLUMN ASSIGNED TO ROW  I .
!!! Z    = COST OF THE OPTIMAL ASSIGNMENT =
!!!      = A(1,F(1)) + A(2,F(2)) + ... + A(N,F(N)) .
!!!
!!! ALL THE PARAMETERS ARE INTEGERS.
!!! VECTOR  F  MUST BE DIMENSIONED AT LEAST AT  N , MATRIX  A
!!! AT LEAST AT  (N,N) . AS CURRENTLY DIMENSIONED, THE SIZE
!!! LIMITATION IS  N .LE. 260 . IN ALL THE SUBROUTINES, THE
!!! INTERNAL VARIABLES WHICH ARE PRESENTLY DIMENSIONED AT
!!! 260 MUST BE DIMENSIONED AT LEAST AT  N .
!!!
!!! THE ONLY MACHINE-DEPENDENT CONSTANT USED IS  INF (DEFINED
!!! BY THE FIRST EXECUTABLE STATEMENT OF THIS SUBROUTINE). INF
!!! MUST BE SET TO A VERY LARGE INTEGER VALUE.
!!!
!!! THE CODE IS BASED ON THE HUNGARIAN METHOD AS DESCRIBED BY
!!! LAWLER (COMBINATORIAL OPTIMIZATION : NETWORKS AND
!!! MATROIDS, HOLT, RINEHART AND WINSTON, NEW YORK, 1976).
!!! THE ALGORITHMIC PASCAL-LIKE DESCRIPTION OF THE CODE IS
!!! GIVEN IN G.CARPANETO, S.MARTELLO AND P.TOTH, ALGORITHMS AND
!!! CODES FOR THE ASSIGNMENT PROBLEM, ANNALS OF OPERATIONS
!!! RESEARCH 7, 1988.
!!!
!!! SUBROUTINE APC DETERMINES THE INITIAL DUAL AND PARTIAL
!!! PRIMAL SOLUTIONS AND THEN SEARCHES FOR AUGMENTING PATHS
!!! UNTIL ALL ROWS AND COLUMNS ARE ASSIGNED.
!!!
!!! MEANING OF THE MAIN INTERNAL VARIABLES:
!!! FB(J) = ROW ASSIGNED TO COLUMN  J .
!!! M     = NUMBER OF INITIAL ASSIGNMENTS.
!!! U(I)  = DUAL VARIABLE ASSOCIATED WITH ROW  I .
!!! V(J)  = DUAL VARIABLE ASSOCIATED WITH COLUMN  J .
!!!
!!! APC NEEDS THE FOLLOWING SUBROUTINES: INCR
!!!                                      INIT
!!!                                      PATH
!!!
!!! ALL THE SUBROUTINES ARE WRITTEN IN AMERICAN NATIONAL
!!! STANDARD FORTRAN AND ARE ACCEPTED BY THE PFORT VERIFIER.
!!!
!!!
!!! THIS WORK WAS SUPPORTED BY  C.N.R. , ITALY.
!!!
!!      INTEGER A(260,260),F(260),Z
!!      INTEGER U,V,FB
!!      COMMON U(260),V(260),FB(260)
!!      INF = 10**9
!!! SEARCH FOR THE INITIAL DUAL AND PARTIAL PRIMAL SOLUTIONS.
!!      CALL INIT(N,A,F,M,INF)
!!      IF ( M .EQ. N ) GO TO 20
!!! SOLUTION OF THE REDUCED PROBLEM.
!!      DO 10 I=1,N
!!        IF ( F(I) .GT. 0 ) GO TO 10
!!! DETERMINATION OF AN AUGMENTING PATH STARTING FROM ROW  I .
!!        CALL PATH(N,A,I,F,INF,J)
!!! ASSIGNMENT OF ROW  I  AND COLUMN  J .
!!        CALL INCR(F,J)
!!   10 CONTINUE
!!! COMPUTATION OF THE SOLUTION COST  Z .
!!   20 Z = 0
!!      DO 30 K=1,N
!!        Z = Z + U(K) + V(K)
!!   30 CONTINUE
!!      RETURN
!!      END
!!      SUBROUTINE INCR(F,J)
!!!
!!! ASSIGNMENT OF COLUMN  J .
!!!
!!      INTEGER F(260)
!!      INTEGER U,V,FB,RC
!!      COMMON U(260),V(260),FB(260),RC(260)
!!   10 I = RC(J)
!!      FB(J) = I
!!      JJ = F(I)
!!      F(I) = J
!!      J = JJ
!!      IF ( J .GT. 0 ) GO TO 10
!!      RETURN
!!      END
!!      SUBROUTINE INIT(N,A,F,M,INF)
!!!
!!! SEARCH FOR THE INITIAL DUAL AND PARTIAL PRIMAL SOLUTIONS.
!!!
!!! P(I) = FIRST UNSCANNED COLUMN OF ROW  I .
!!!
!!      INTEGER A(260,260),F(260)
!!      INTEGER U,V,FB,P,R
!!      COMMON U(260),V(260),FB(260),P(260)
!!! PHASE 1 .
!!      M = 0
!!      DO 10 K=1,N
!!        F(K) = 0
!!        FB(K) = 0
!!   10 CONTINUE
!!! SCANNING OF THE COLUMNS ( INITIALIZATION OF  V(J) ).
!!      DO 40 J=1,N
!!        MIN = INF
!!        DO 30 I=1,N
!!          IA = A(I,J)
!!          IF ( IA .GT. MIN ) GO TO 30
!!          IF ( IA .LT. MIN ) GO TO 20
!!          IF ( F(I) .NE. 0 ) GO TO 30
!!   20     MIN = IA
!!          R = I
!!   30   CONTINUE
!!        V(J) = MIN
!!        IF ( F(R) .NE. 0 ) GO TO 40
!!! ASSIGNMENT OF COLUMN  J  TO ROW  R .
!!        M = M + 1
!!        FB(J) = R
!!        F(R) = J
!!        U(R) = 0
!!        P(R) = J + 1
!!   40 CONTINUE
!!! PHASE 2 .
!!! SCANNING OF THE UNASSIGNED ROWS ( UPDATING OF  U(I) ).
!!      DO 110 I=1,N
!!        IF ( F(I) .NE. 0 ) GO TO 110
!!        MIN = INF
!!        DO 60 K=1,N
!!          IA = A(I,K) - V(K)
!!          IF ( IA .GT. MIN ) GO TO 60
!!          IF ( IA .LT. MIN ) GO TO 50
!!          IF ( FB(K) .NE. 0 ) GO TO 60
!!          IF ( FB(J) .EQ. 0 ) GO TO 60
!!   50     MIN = IA
!!          J = K
!!   60   CONTINUE
!!        U(I) = MIN
!!        JMIN = J
!!        IF ( FB(J) .EQ. 0 ) GO TO 100
!!        DO 80 J=JMIN,N
!!          IF ( A(I,J) - V(J) .GT. MIN ) GO TO 80
!!          R = FB(J)
!!          KK = P(R)
!!          IF ( KK .GT. N ) GO TO 80
!!          DO 70 K=KK,N
!!            IF ( FB(K) .GT. 0 ) GO TO 70
!!            IF ( A(R,K) - U(R) - V(K) .EQ. 0 ) GO TO 90
!!   70     CONTINUE
!!          P(R) = N + 1
!!   80   CONTINUE
!!        GO TO 110
!!! REASSIGNMENT OF ROW  R  AND COLUMN  K .
!!   90   F(R) = K
!!        FB(K) = R
!!        P(R) = K + 1
!!! ASSIGNMENT OF COLUMN  J  TO ROW  I .
!!  100   M = M + 1
!!        F(I) = J
!!        FB(J)= I
!!        P(I) = J + 1
!!  110 CONTINUE
!!      RETURN
!!      END
!!      SUBROUTINE PATH(N,A,II,F,INF,JJ)
!!!
!!! DETERMINATION OF AN AUGMENTING PATH STARTING FROM
!!! UNASSIGNED ROW  II  AND TERMINATING AT UNASSIGNED COLUMN
!!! JJ , WITH UPDATING OF DUAL VARIABLES  U(I)  AND  V(J) .
!!!
!!! MEANING OF THE MAIN INTERNAL VARIABLES:
!!! LR(L) = L-TH LABELLED ROW ( L=1,NLR ).
!!! PI(J) = MIN ( A(I,J) - U(I) - V(J) , SUCH THAT ROW  I  IS
!!!         LABELLED AND NOT EQUAL TO  FB(J) ).
!!! RC(J) = ROW PRECEDING COLUMN  J  IN THE CURRENT
!!!         ALTERNATING PATH.
!!! UC(L) = L-TH UNLABELLED COLUMN ( L=1,NUC ).
!!!
!!      INTEGER A(260,260),F(260),Z
!!      INTEGER PI(260),LR(260),UC(260)
!!      INTEGER U,V,FB,RC,R
!!      COMMON U(260),V(260),FB(260),RC(260)
!!! INITIALIZATION.
!!      LR(1) = II
!!      DO 10 K=1,N
!!        PI(K) = A(II,K) - U(II) - V(K)
!!        RC(K) = II
!!        UC(K) = K
!!   10 CONTINUE
!!      NUC = N
!!      NLR = 1
!!      GO TO 40
!!! SCANNING OF THE LABELLED ROWS.
!!   20 R = LR(NLR)
!!      DO 30 L=1,NUC
!!        J = UC(L)
!!        IA = A(R,J) - U(R) - V(J)
!!        IF ( IA .GE. PI(J) ) GO TO 30
!!        PI(J) = IA
!!        RC(J) = R
!!   30 CONTINUE
!!! SEARCH FOR A ZERO ELEMENT IN AN UNLABELLED COLUMN.
!!   40 DO 50 L=1,NUC
!!        J = UC(L)
!!        IF ( PI(J) .EQ. 0 ) GO TO 100
!!   50 CONTINUE
!!! UPDATING OF THE DUAL VARIABLES  U(I)  AND  V(J) .
!!      MIN = INF
!!      DO 60 L=1,NUC
!!        J = UC(L)
!!        IF ( MIN .GT. PI(J) ) MIN = PI(J)
!!   60 CONTINUE
!!      DO 70 L=1,NLR
!!        R = LR(L)
!!        U(R) = U(R) + MIN
!!   70 CONTINUE
!!      DO 90 J=1,N
!!        IF ( PI(J) .EQ. 0 ) GO TO 80
!!        PI(J) = PI(J) - MIN
!!        GO TO 90
!!   80   V(J) = V(J) - MIN
!!   90 CONTINUE
!!      GO TO 40
!!  100 IF ( FB(J) .EQ. 0 ) GO TO 110
!!! LABELLING OF ROW  FB(J)  AND REMOVAL OF THE LABEL  OF
!!! COLUMN  J .
!!      NLR = NLR + 1
!!      LR(NLR) = FB(J)
!!      UC(L) = UC(NUC)
!!      NUC = NUC - 1
!!      GO TO 20
!!! DETERMINATION OF THE UNASSIGNED COLUMN  J .
!!  110 JJ = J
!!      RETURN
!!      END
!! 
!!! Code converted using TO_F90 by Alan Miller
!!! Date: 2002-03-06  Time: 08:36:31
!!
!!! Converted with permission, from the F77 code in the CERN MATHLIB library.
!!! $Id: assndx.F,v 1.1.1.1 1996/04/01 15:02:49 mclareni Exp $
!!
!!! $Log: assndx.F,v $
!!! Revision 1.1.1.1  1996/04/01 15:02:49  mclareni
!!! Mathlib gen/H (H301)
!!! Author: F. Bourgeois, 15 February 1994
!!
